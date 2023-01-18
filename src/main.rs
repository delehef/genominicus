use anyhow::*;
use anyhow::{Context, Result};
use clap::*;
use colored::Colorize;
use log::*;
use rusqlite::*;
use utils::*;

mod align;
mod render;
mod utils;

#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    #[clap(flatten)]
    verbose: clap_verbosity_flag::Verbosity,

    /// The half-length of the syntenic context
    #[clap(short, long, default_value_t = 15)]
    window: usize,

    #[clap(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Create a syntenic database of the provided genomes
    BuildDatabase {
        /// the files and/or directories containing the gene families to process
        #[clap(long, required = true)]
        families: Vec<String>,

        /// where to write the database
        #[clap(long = "out", short = 'o')]
        outfile: String,

        /// the directory containing the GFF3 files to process; those can be gzipped
        #[clap(long)]
        gffs: Vec<String>,

        /// the features to extract from GFF files
        #[clap(long, default_value_t = String::from("gene"))]
        id_type: String,

        /// regex to extract feature name from GFF ID field; must contain a named capture group `id`
        #[clap(long, default_value_t = String::from("gene:(?P<id>.*)"))]
        id_pattern: String,

        /// regex to extract species name from GFF file name; must contain a named capture group `species`
        #[clap(long, default_value_t = String::from("(?P<species>.*).gff3"))]
        species_pattern: String,
    },
    /// Render one or more gene trees, with their syntenic environment stored in the provided database
    Plot {
        /// The gene trees to render
        #[arg(required = true)]
        files: Vec<String>,

        #[arg(short = 'o')]
        out: Option<String>,

        /// The database containing the syntenic environment of each gene, as built with `build-database`
        #[arg(short = 'D', long = "database")]
        database: String,

        /// The species tree
        #[arg(short = 'S', required_if_eq("graph_type", "barcode"))]
        species_tree: Option<String>,

        #[arg(short = 'T', long = "type", default_value = "flat", value_parser=["flat", "html", "barcode", "skeleton"])]
        graph_type: String,

        #[arg(
            short = 'I',
            long = "id",
            help = "the column name mapping to the IDs in the gene trees",
            default_value = "id"
        )]
        id_column: String,

        /// If set, introduce a new set of color gradients at each duplication node
        #[arg()]
        colorize_per_duplication: bool,

        /// Ensure that all genes are colorized
        #[arg()]
        colorize_all: bool,

        /// Filter out species present in the species tree but not in the gene tree
        #[arg(long = "filter-species")]
        filter_species_tree: bool,

        /// Additional annotations to the plot
        #[arg(long="annotations", value_delimiter = ',', value_parser=["links", "inner-nodes", "cs", "elc", "ellc"])]
        annotations: Vec<String>,
    },
}

fn main() -> Result<()> {
    let args = Args::parse();
    buche::new()
        .timestamp(buche::Timestamp::Off)
        .verbosity(args.verbose.log_level_filter())
        .init()
        .unwrap();

    match args.command {
        Commands::BuildDatabase {
            families,
            outfile,
            gffs,
            id_type,
            id_pattern,
            species_pattern,
        } => syntesuite::dbmaker::db_from_gffs(
            &families,
            &gffs,
            &outfile,
            &species_pattern,
            &id_type,
            &id_pattern,
            args.window as isize,
        ),
        Commands::Plot {
            files,
            out,
            database,
            species_tree,
            graph_type,
            id_column,
            colorize_per_duplication,
            colorize_all,
            filter_species_tree,
            annotations,
        } => {
            utils::set_reference(&id_column);

            let mut render_settings = RenderSettings::default();
            for annotation in annotations {
                match annotation.as_str() {
                    "links" => render_settings.links = true,
                    "cs" => render_settings.cs = true,
                    "elc" => render_settings.elc = true,
                    "ellc" => render_settings.ellc = true,
                    "inner-nodes" => render_settings.inner_nodes = true,
                    _ => unreachable!(),
                }
            }

            for filename in files.iter() {
                info!(
                    "Rendering {} as {}",
                    filename.bold().bright_white(),
                    graph_type.bold().yellow()
                );
                let mut out_filename =
                    std::path::PathBuf::from(out.clone().unwrap_or_else(|| filename.to_string()));
                out_filename.set_file_name(
                    out_filename
                        .file_stem()
                        .with_context(|| {
                            anyhow!(
                                "invalid file name: {}",
                                out_filename.to_str().unwrap().bold().yellow()
                            )
                        })?
                        .to_owned(),
                );
                let out_filename = out_filename.to_str().unwrap();
                let t = newick::one_from_filename(filename)
                    .context(format!("failed to read `{}`", filename))?;
                match graph_type.as_str() {
                    "flat" => {
                        let mut db = Connection::open_with_flags(
                            &database,
                            OpenFlags::SQLITE_OPEN_READ_ONLY,
                        )
                        .unwrap();
                        let genes = make_genes_cache(&t, &mut db);
                        let colormap = if colorize_per_duplication {
                            make_colormap_per_duplication(&t, &genes, colorize_all)
                        } else {
                            make_colormap(&t, &genes)
                        };
                        render::flat::render(
                            &t,
                            &genes,
                            &colormap,
                            &format!("{}-flat.svg", out_filename),
                            &render_settings,
                        );
                    }
                    "html" => {
                        let mut db = Connection::open_with_flags(
                            &database,
                            OpenFlags::SQLITE_OPEN_READ_ONLY,
                        )
                        .unwrap();
                        let genes = make_genes_cache(&t, &mut db);
                        let colormap = if colorize_per_duplication {
                            make_colormap_per_duplication(&t, &genes, colorize_all)
                        } else {
                            make_colormap(&t, &genes)
                        };
                        render::html::render(
                            &t,
                            &genes,
                            &colormap,
                            &format!("{}.html", out_filename),
                        )
                    }
                    "barcode" => {
                        render::barcode::render(
                            &t,
                            species_tree.as_ref().unwrap(),
                            &format!("{}-barcode.svg", out_filename),
                            filter_species_tree,
                            &render_settings,
                        );
                    }
                    "skeleton" => {
                        render::skeleton::render(
                            &t,
                            &format!("{}-skeleton.svg", out_filename),
                            &render_settings,
                        );
                    }
                    _ => unimplemented!(),
                };
            }
            Ok(())
        }
    }
}
