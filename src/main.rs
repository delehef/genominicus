use clap::*;
use colored::Colorize;
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
            Arg::new("type")
                .help("The graph type to use")
                .short('T')
                .long("type")
                .possible_values(&["flat", "html", "barcode", "skeleton"])
                .default_value("flat")
                .takes_value(true),
        )
        .arg(
            Arg::new("species-tree")
                .short('S')
                .long("species-tree")
                .help("The species tree to plot against")
                .required_if_eq("type", "barcode")
                .takes_value(true),
        )
        .arg(
            Arg::new("database")
                .short('D')
                .long("database")
                .help("The database to use")
                .required_if_eq_any(&[("type", "barcode"), ("type", "html")])
                .takes_value(true),
        )
        .arg(
            Arg::new("FILE")
                .help("Sets the input file to use")
                .required(true)
                .multiple_values(true),
        )
        .arg(
            Arg::new("OUT")
                .short('o')
                .long("out")
                .help("Output filename")
                .takes_value(true),
        )
        .arg(
            Arg::new("colorize_per_duplication")
                .help("Introduce a new set of color gradients at each duplicated node")
                .long("colorize-duplications"),
        )
        .arg(
            Arg::new("colorize_all")
                .help("Ensure that all genes are colorized")
                .long("colorize-all"),
        )
        .arg(
            Arg::new("filter_species_tree")
                .help("Filter out species present in the species tree but not in the gene tree")
                .long("filter-species"),
        )
        .arg(
            Arg::new("annotations")
                .long("annotations")
                .value_delimiter(',')
                .multiple_values(true)
                .possible_values(["links", "inner-nodes", "cs", "elc", "ellc"]),
        )
        .get_matches();

    let graph_type = args.value_of("type").unwrap();
    let colorize_per_duplication = args.is_present("colorize_per_duplication");
    let colorize_all = args.is_present("colorize_all");

    let mut render_settings = RenderSettings::default();
    if let Some(annotations) = args.values_of("annotations") {
        for annotation in annotations {
            match annotation {
                "links" => render_settings.links = true,
                "cs" => render_settings.cs = true,
                "elc" => render_settings.elc = true,
                "ellc" => render_settings.ellc = true,
                "inner-nodes" => render_settings.inner_nodes = true,
                _ => unreachable!(),
            }
        }
    };

    for filename in args.values_of_t::<String>("FILE").unwrap().iter() {
        println!("Rendering {} as {}", filename.bold().magenta(), graph_type.bold().yellow());
        let mut out_filename = std::path::PathBuf::from(
            args.value_of_t("OUT")
                .unwrap_or_else(|_| filename.to_string()),
        );
        out_filename.set_file_name(out_filename.file_stem().unwrap().to_owned());
        let out_filename = out_filename.to_str().unwrap();
        let t = newick::from_filename(filename).unwrap();
        match graph_type {
            "flat" => {
                let db_filename = args.value_of("database").unwrap();
                let mut db =
                    Connection::open_with_flags(db_filename, OpenFlags::SQLITE_OPEN_READ_ONLY)
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
                let db_filename = args.value_of("database").unwrap();
                let mut db =
                    Connection::open_with_flags(db_filename, OpenFlags::SQLITE_OPEN_READ_ONLY)
                        .unwrap();
                let genes = make_genes_cache(&t, &mut db);
                let colormap = if colorize_per_duplication {
                    make_colormap_per_duplication(&t, &genes, colorize_all)
                } else {
                    make_colormap(&t, &genes)
                };
                render::html::render(&t, &genes, &colormap, &format!("{}.html", out_filename))
            }
            "barcode" => {
                render::barcode::render(
                    &t,
                    args.value_of("species-tree").unwrap(),
                    &format!("{}-barcode.svg", out_filename),
                    args.is_present("filter_species_tree"),
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
}
