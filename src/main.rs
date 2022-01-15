use clap::*;
use newick::*;
use rusqlite::*;
use utils::*;

mod align;
mod render;
mod utils;

fn main() {
    let args = App::new("Genominicus")
        .version(clap::crate_version!())
        .author(clap::crate_authors!())
        .arg(
            Arg::with_name("type")
                .help("The graph type to use")
                .short("T")
                .long("type")
                .possible_values(&["flat", "html", "barcode"])
                .default_value("flat")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("species-tree")
                .short("S")
                .long("species-tree")
                .help("The species tree to plot against")
                .required_if("type", "barcode")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("database")
                .short("D")
                .long("database")
                .help("The database to use")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("FILE")
                .help("Sets the input file to use")
                .required(true)
                .multiple(true),
        )
        .arg(
            Arg::with_name("OUT")
                .short("o")
                .long("out")
                .help("Output filename")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("colorize_per_duplication")
                .help("Introduce a new set of color gradients at each duplicated node")
                .long("colorize-duplications"),
        )
        .arg(
            Arg::with_name("colorize_all")
                .help("Ensure that all genes are colorized")
                .long("colorize-all"),
        )
        .arg(
            Arg::with_name("filter_species_tree")
                .help("Filter out species present in the species tree but not in the gene tree")
                .long("filter-species"),
        )
        .get_matches();

    let db_filename = value_t!(args, "database", String).unwrap();
    let graph_type = value_t!(args, "type", String).unwrap();
    let colorize_per_duplication = args.is_present("colorize_per_duplication");
    let colorize_all = args.is_present("colorize_all");

    let mut db =
        Connection::open_with_flags(db_filename, OpenFlags::SQLITE_OPEN_READ_ONLY).unwrap();

    for filename in values_t!(args, "FILE", String).unwrap().iter() {
        println!("Processing {}", filename);
        let mut out_filename =
            std::path::PathBuf::from(value_t!(args, "OUT", String).unwrap_or_else(|_| filename.to_string()));
        out_filename.set_file_name(out_filename.file_stem().unwrap().to_owned());
        let out_filename = out_filename.to_str().unwrap();

        let t = Tree::from_filename(filename).unwrap();
        let genes = make_genes_cache(&t, &mut db);
        let colormap = if colorize_per_duplication {
            make_colormap_per_duplication(&t, &genes, colorize_all)
        } else {
            make_colormap(&t, &genes)
        };
        match graph_type.as_str() {
            "flat" => {
                render::flat::render(&t, &genes, &colormap, &format!("{}-flat.svg", out_filename));
            }
            "html" => {
                render::html::render(&t, &genes, &colormap, &format!("{}.html", out_filename))
            }
            "barcode" => {
                render::barcode::render(
                    &t,
                    &value_t!(args, "species-tree", String).unwrap(),
                    &format!("{}-barcode.svg", out_filename),
                    args.is_present("filter_species_tree"),
                );
            }
            _ => unimplemented!(),
        };
    }
}
