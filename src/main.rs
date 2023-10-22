use anyhow::*;
use anyhow::{Context, Result};
use clap::*;
use colored::Colorize;
use log::*;
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

        /// the directory containing the genomes to process; those can be either in the GFF3 or BED format, and may be gzipped
        #[clap(long)]
        genome_files: Vec<String>,

        /// the features to extract from GFF files
        #[clap(long, default_value_t = String::from("gene"))]
        id_type: String,

        /// regex to extract feature name from either a GFF ID attribute or a BED name field; must contain a named capture group `id`
        #[clap(long, default_value_t = String::from("gene:(?P<id>.*)"))]
        id_pattern: String,

        /// regex to extract species name from genome files; must contain a named capture group `species`
        #[clap(long, default_value_t = String::from("(?P<species>.*)\\.(gff3|bed)"))]
        species_pattern: String,
    },
    /// Render one or more gene trees, with their syntenic environment stored in the provided database
    Plot {
        /// The gene trees to render
        #[arg(required = true)]
        files: Vec<String>,

        /// Explicitely set an output file name
        #[arg(short, long)]
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
        #[arg(long)]
        colorize_per_duplication: bool,

        /// Ensure that all genes are colorized
        #[arg(long)]
        colorize_all: bool,

        /// Filter out species present in the species tree but not in the gene tree
        #[arg(long = "filter-species")]
        filter_species_tree: bool,

        /// Additional annotations to the plot
        #[arg(long="annotations", value_delimiter = ',', value_parser=["links", "inner-nodes", "cs", "elc", "ellc", "dids"])]
        annotations: Vec<String>,

        /// Display the plot after creation. If a program name is passed, use it to open the plot; otherwise use
        /// the system default
        #[arg(short = 'O', long)]
        open: Option<Option<String>>,
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
            genome_files: gffs,
            id_type,
            id_pattern,
            species_pattern,
        } => syntesuite::dbmaker::db_from_files(
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
            open,
        } => {
            let mut render_settings = RenderSettings::default();
            for annotation in annotations {
                match annotation.as_str() {
                    "links" => render_settings.links = true,
                    "cs" => render_settings.cs = true,
                    "elc" => render_settings.elc = true,
                    "ellc" => render_settings.ellc = true,
                    "inner-nodes" => render_settings.inner_nodes = true,
                    "dids" => render_settings.duplication_ids = true,
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
                let out = match graph_type.as_str() {
                    "flat" => {
                        let genes = make_genes_cache(&t, &database, &id_column)?;
                        let colormap = if colorize_per_duplication {
                            make_colormap_per_duplication(&t, &genes, colorize_all)
                        } else {
                            make_colormap(&t, &genes)
                        };
                        let petmap = make_petnamemap(&t, &genes);
                        let out = format!("{}-flat.svg", out_filename);
                        render::flat::render(
                            &t,
                            &genes,
                            &colormap,
                            &petmap,
                            &out,
                            &render_settings,
                        );
                        out
                    }
                    "html" => {
                        let genes = make_genes_cache(&t, &database, &id_column)?;
                        let colormap = if colorize_per_duplication {
                            make_colormap_per_duplication(&t, &genes, colorize_all)
                        } else {
                            make_colormap(&t, &genes)
                        };
                        let out = format!("{}.html", out_filename);
                        render::html::render(&t, &genes, &colormap, &out);
                        out
                    }
                    "barcode" => {
                        let out = format!("{}-barcode.svg", out_filename);
                        render::barcode::render(
                            &t,
                            species_tree.as_ref().unwrap(),
                            &out,
                            filter_species_tree,
                            &render_settings,
                        );
                        out
                    }
                    "skeleton" => {
                        let out = format!("{}-skeleton.svg", out_filename);
                        render::skeleton::render(&t, &out, &render_settings);
                        out
                    }
                    _ => unimplemented!(),
                };
                if let Some(open_with) = open.as_ref() {
                    if let Some(program) = open_with.as_ref() {
                        open::with(&out, program)?;
                    } else {
                        open::that(&out)?;
                    };
                };
            }
            Ok(())
        }
    }
}
