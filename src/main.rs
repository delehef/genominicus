use anyhow::{Context, Result};
use clap::*;
use colored::Colorize;
use rusqlite::*;
use utils::*;

mod align;
mod render;
mod utils;

#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    #[arg(required = true)]
    files: Vec<String>,

    #[arg(short = 'o')]
    out: Option<String>,

    #[arg(short = 'T', long = "type", default_value = "flat", value_parser=["flat", "html", "barcode", "skeleton"])]
    graph_type: String,

    #[arg(
        short = 'S',
        long = "species-tree",
        required_if_eq("graph_type", "barcode")
    )]
    species_tree: Option<String>,

    #[arg(
        short='D',
        long="database",
        required_if_eq_any([("graph_type", "barcode"), ("graph_type", "html"), ("graph_type", "flat")])
    )]
    database: Option<String>,

    #[arg(
        long = "colorize-per-duplication",
        help = "Introduce a new set of color gradients at each duplicated node"
    )]
    colorize_per_duplication: bool,

    #[arg(long = "colorize-all", help = "Ensure that all genes are colorized")]
    colorize_all: bool,

    #[arg(
        long = "filter-species",
        help = "Filter out species present in the species tree but not in the gene tree"
    )]
    filter_species_tree: bool,

    #[arg(long="annotations", value_delimiter = ',', value_parser=["links", "inner-nodes", "cs", "elc", "ellc"])]
    annotations: Vec<String>,

    #[arg(
        short = 'R',
        long = "reference",
        help = "whether the forest trees are annotated with genes or protein",
        default_value = "gene"
    )]
    reference: String,
}

fn main() -> Result<()> {
    let args = Args::parse();
    utils::set_reference(&args.reference);

    let mut render_settings = RenderSettings::default();
    for annotation in args.annotations {
        match annotation.as_str() {
            "links" => render_settings.links = true,
            "cs" => render_settings.cs = true,
            "elc" => render_settings.elc = true,
            "ellc" => render_settings.ellc = true,
            "inner-nodes" => render_settings.inner_nodes = true,
            _ => unreachable!(),
        }
    }

    for filename in args.files.iter() {
        println!(
            "Rendering {} as {}",
            filename.bold().bright_white(),
            args.graph_type.bold().yellow()
        );
        let mut out_filename =
            std::path::PathBuf::from(args.out.clone().unwrap_or_else(|| filename.to_string()));
        out_filename.set_file_name(out_filename.file_stem().unwrap().to_owned());
        let out_filename = out_filename.to_str().unwrap();
        let t = newick::one_from_filename(filename)
            .context(format!("failed to read `{}`", filename))?;
        match args.graph_type.as_str() {
            "flat" => {
                let mut db = Connection::open_with_flags(
                    args.database.as_ref().unwrap(),
                    OpenFlags::SQLITE_OPEN_READ_ONLY,
                )
                .unwrap();
                let genes = make_genes_cache(&t, &mut db);
                let colormap = if args.colorize_per_duplication {
                    make_colormap_per_duplication(&t, &genes, args.colorize_all)
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
                    &args.database.as_ref().unwrap(),
                    OpenFlags::SQLITE_OPEN_READ_ONLY,
                )
                .unwrap();
                let genes = make_genes_cache(&t, &mut db);
                let colormap = if args.colorize_per_duplication {
                    make_colormap_per_duplication(&t, &genes, args.colorize_all)
                } else {
                    make_colormap(&t, &genes)
                };
                render::html::render(&t, &genes, &colormap, &format!("{}.html", out_filename))
            }
            "barcode" => {
                render::barcode::render(
                    &t,
                    args.species_tree.as_ref().unwrap(),
                    &format!("{}-barcode.svg", out_filename),
                    args.filter_species_tree,
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
