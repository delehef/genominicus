use clap::*;
use newick::*;
use rusqlite::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use svarog::*;
use utils::*;

mod nw;
mod poa;
mod utils;

pub type SeqID = usize;
pub type Nucleotide = String;
pub type Sequence = Vec<Nucleotide>;
pub type Sequences = HashMap<SeqID, Sequence>;

const DROP_EMPTY: bool = true;

#[derive(Serialize, Deserialize)]
struct PolyGene {
    genes: Vec<(Gene, f32)>,
}
#[derive(Serialize, Deserialize)]
struct Gene {
    color: String,
    name: String,
}
#[derive(Serialize, Deserialize)]
struct Landscape {
    lefts: Vec<Gene>,
    me: Gene,
    rights: Vec<Gene>,
}
#[derive(Serialize, Deserialize)]
struct HtmlNode {
    species: String,
    chr: String,
    gene: String,
    ancestral: String,
    t_len: i32,
    color: String,
    children: Vec<HtmlNode>,

    isDuplication: bool,
    confidence: f32,
    repr: Landscape,
    clustered: Option<Vec<PolyGene>>,
}

fn nuc_to_str(nuc: &str) -> String {
    nuc.to_string()
}

fn draw_background(
    svg: &mut SvgDrawing,
    genes: &GeneCache,
    depth: f32,
    tree: &Tree,
    node: &Node,
    xoffset: f32,
    yoffset: f32,
    xlabels: f32,
    width: f32,
) -> f32 {
    let mut y = yoffset;

    let mut children = node
        .children
        .as_ref()
        .map(|children| children.iter().map(|i| &tree[*i]).collect::<Vec<_>>())
        .unwrap_or_default();
    children.sort_by_cached_key(|x| {
        if let Some(DbGene { species, .. }) = x
            .name
            .as_ref()
            .and_then(|name| name.split('#').next())
            .and_then(|gene_name| genes.get(&gene_name.to_string()))
        {
            species.as_str()
        } else {
            "zzz"
        }
    });

    for &child in children.iter() {
        let new_y = if child.is_leaf() {
            y + 20.
        } else {
            draw_background(
                svg,
                genes,
                depth,
                tree,
                child,
                xoffset + BRANCH_WIDTH,
                y,
                xlabels,
                width,
            )
        };

        if node.is_duplication() {
            let d = xoffset / depth;
            svg.polygon()
                .from_pos_dims(
                    xoffset + BRANCH_WIDTH / 2.,
                    y - 6.,
                    width - xoffset - d * BRANCH_WIDTH,
                    new_y - y - 6.,
                )
                .style(|s| {
                    s.fill_color(StyleColor::Percent(0.5, 0.5, 1.))
                        .fill_opacity(0.1 + 0.9 * d)
                });
        }
        y = new_y;
    }
    y
}

fn draw_gene<'a>(
    svg: &'a mut SvgDrawing,
    x: f32,
    y: f32,
    right: bool,
    color: &StyleColor,
    name: &str,
) -> &'a mut Polygon {
    if right {
        svg.polygon()
            .add_point(x, y - 5.)
            .add_point(x + GENE_WIDTH - 3., y - 5.)
            .add_point(x + GENE_WIDTH, y)
            .add_point(x + GENE_WIDTH - 3., y + 5.)
            .add_point(x, y + 5.)
            .set_hover(name)
            .style(|s| {
                s.fill_color(color.clone())
                    .stroke_width(0.5)
                    .stroke_color(StyleColor::Percent(0.2, 0.2, 0.2))
            })
    } else {
        svg.polygon()
            .add_point(x, y)
            .add_point(x + 3., y - 5.)
            .add_point(x + GENE_WIDTH, y - 5.)
            .add_point(x + GENE_WIDTH, y + 5.)
            .add_point(x + 3., y + 5.)
            .set_hover(name)
            .style(|s| {
                s.fill_color(color.clone())
                    .stroke_width(0.5)
                    .stroke_color(StyleColor::Percent(0.2, 0.2, 0.2))
            })
    }
}

fn draw_tree(
    svg: &mut SvgDrawing,
    genes: &GeneCache,
    colormap: &ColorMap,
    depth: f32,
    tree: &Tree,
    node: &Node,
    xoffset: f32,
    yoffset: f32,
    xlabels: f32,
    width: f32,
    links: &mut Vec<(Vec<String>, String, Vec<String>)>,
) -> f32 {
    let mut y = yoffset;
    let mut old_y = 0.;
    let children = node
        .children
        .as_ref()
        .map(|children| children.iter().map(|i| &tree[*i]).collect::<Vec<_>>())
        .unwrap_or_default();
    for (i, child) in children.iter().enumerate() {
        if i > 0 {
            svg.line()
                .from_coords(xoffset, old_y, xoffset, y)
                .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5));
        }
        old_y = y;

        if child.is_leaf() {
            // Leaf branch
            svg.line()
                .from_coords(xoffset, y, depth, y)
                .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5));

            // Landscape support line
            svg.line()
                .from_points([
                    (xlabels - 5., y),
                    (
                        xlabels + (GENE_WIDTH + GENE_SPACING) * (2. * WINDOW as f32 + 1.)
                            - GENE_SPACING
                            + 5.,
                        y,
                    ),
                ])
                .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5));

            if let Some(name) = &child.name {
                let gene_name = name.split('#').next().unwrap();
                if let Some(DbGene {
                    ancestral,
                    species,
                    chr,
                    left_tail,
                    right_tail,
                    ..
                }) = genes.get(gene_name)
                {
                    // Gene/protein name
                    svg.text()
                        .pos(depth, y + 5.)
                        .text(format!("{} {}/{}", gene_name, species, chr))
                        .style(|s| s.fill_color(name2color(species)));

                    // Left tail
                    let xbase = xlabels + (WINDOW as f32 - 1.) * (GENE_WIDTH + GENE_SPACING);
                    for (k, g) in left_tail.iter().enumerate() {
                        let xstart = xbase - (k as f32) * (GENE_WIDTH + GENE_SPACING);
                        let drawn = draw_gene(
                            svg,
                            xstart,
                            y,
                            true,
                            colormap
                                .get(&g.clone())
                                .unwrap_or(&StyleColor::String("#aaa".to_string())),
                            g,
                        );
                        if g == ancestral {
                            drawn.style(|s| {
                                s.stroke_width(2.)
                                    .stroke_color(StyleColor::Percent(0.1, 0.1, 0.1))
                            });
                        }
                    }

                    // The Gene
                    draw_gene(
                        svg,
                        xlabels + WINDOW as f32 * (GENE_WIDTH + GENE_SPACING),
                        y,
                        true,
                        &gene2color(&ancestral),
                        &ancestral,
                    )
                    .style(|s| {
                        s.stroke_width(2.)
                            .stroke_color(StyleColor::Percent(0.1, 0.1, 0.1))
                    });

                    // Right tail
                    let xbase = xlabels + (WINDOW as f32 + 1.) * (GENE_WIDTH + GENE_SPACING);
                    for (k, g) in right_tail.iter().enumerate() {
                        let xstart = xbase + (k as f32) * (GENE_WIDTH + GENE_SPACING);
                        let drawn = draw_gene(
                            svg,
                            xstart,
                            y,
                            true, // g.1.to_string() == direction,
                            colormap
                                .get(&g.clone())
                                .unwrap_or(&StyleColor::String("#aaa".to_string())),
                            g,
                        );
                        if g == ancestral {
                            drawn.style(|s| {
                                s.stroke_width(2.)
                                    .stroke_color(StyleColor::Percent(0.1, 0.1, 0.1))
                            });
                        }
                    }
                    links.push((left_tail.to_vec(), ancestral.into(), right_tail.to_vec()));
                } else {
                    // The node was not found in the database
                    eprintln!("{} -- {} not found", name, gene_name);
                    links.push((Vec::new(), name.into(), Vec::new()));
                }
                y += 20.;
            }
        } else {
            svg.line()
                .from_coords(xoffset, y, xoffset + BRANCH_WIDTH, y)
                .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5));
            y = draw_tree(
                svg,
                genes,
                colormap,
                depth,
                tree,
                child,
                xoffset + BRANCH_WIDTH,
                y,
                xlabels,
                width,
                links,
            );
        }
    }

    if node.is_duplication() {
        let dcs = node.data.get("DCS").and_then(|s| s.parse::<f32>().ok());

        svg.text()
            .pos(xoffset - FONT_SIZE, yoffset - FONT_SIZE)
            .text(
                dcs.map(|s| format!("DCS = {:.2}", s))
                    .unwrap_or("?".to_string()),

            );
        let dcs = dcs.unwrap_or(0.0);
        svg.polygon()
            .from_pos_dims(xoffset - 3., yoffset - 3., 6., 6.)
            .style(|s| s.fill_color(StyleColor::Percent(1.0 - dcs, dcs, 0.)));
    }
    y
}

fn draw_links(
    svg: &mut SvgDrawing,
    links: &[(Vec<String>, String, Vec<String>)],
    yoffset: f32,
    xlabels: f32,
) {
    let mut y = yoffset;
    for w in links.windows(2) {
        let xbase = xlabels + (WINDOW as f32 - 1.) * (GENE_WIDTH + GENE_SPACING);
        for (i, ancestral) in w[0].0.iter().enumerate() {
            let x1 = xbase - i as f32 * (GENE_WIDTH + GENE_SPACING) + GENE_WIDTH / 2.;
            for j in
                w[1].0
                    .iter()
                    .enumerate()
                    .filter_map(|(j, name)| if name == ancestral { Some(j) } else { None })
            {
                let x2 = xbase - j as f32 * (GENE_WIDTH + GENE_SPACING) + GENE_WIDTH / 2.;
                svg.line()
                    .from_points([(x1, y + 5.), (x2, y + 20. - 5.)])
                    .style(|s| {
                        s.stroke_color(StyleColor::String("#000".into()))
                            .stroke_width(1.0)
                            .dashed(&[2, 2])
                    });
            }
        }

        let xbase = xlabels + (WINDOW as f32 + 1.) * (GENE_WIDTH + GENE_SPACING);
        for (i, ancestral) in w[0].2.iter().enumerate() {
            let x1 = xbase + i as f32 * (GENE_WIDTH + GENE_SPACING) + GENE_WIDTH / 2.;
            for j in
                w[1].2
                    .iter()
                    .enumerate()
                    .filter_map(|(j, name)| if name == ancestral { Some(j) } else { None })
            {
                let x2 = xbase + j as f32 * (GENE_WIDTH + GENE_SPACING) + GENE_WIDTH / 2.;
                svg.line()
                    .from_points([(x1, y + 5.), (x2, y + 20. - 5.)])
                    .style(|s| {
                        s.stroke_color(StyleColor::String("#000".into()))
                            .stroke_width(1.0)
                            .dashed(&[2, 2])
                    });
            }
        }

        y += 20.;
    }
}

fn draw_html(tree: &Tree, genes: &GeneCache, colormap: &ColorMap) -> HtmlNode {
    fn process(tree: &Tree, node: usize, genes: &GeneCache, colormap: &ColorMap) -> HtmlNode {
        let descendants = tree.descendants(node);
        let mut common_ancestral = String::new();
        let clustered = {
            // Node ID to tail mapping
            let tails = descendants
                .iter()
                .filter_map(|&d| {
                    if let Some(name) = &tree[d].name {
                        let gene_name = name.split('#').next().unwrap();
                        if let Some(DbGene {
                            ancestral,
                            left_tail,
                            right_tail,
                            ..
                        }) = genes.get(gene_name)
                        {
                            common_ancestral = ancestral.to_owned();
                            Some((
                                d,
                                left_tail
                                    .iter()
                                    .rev() // XXX Pour que les POA partent bien du bout
                                    .chain([MARKER.to_owned()].iter())
                                    .chain(right_tail.iter())
                                    .cloned()
                                    .collect::<Vec<String>>(),
                            ))
                        } else {
                            eprintln!("{} -- {} not found", name, gene_name);
                            None
                        }
                    } else {
                        None
                    }
                })
                .collect::<HashMap<_, _>>();

            if !tails.is_empty() {
                let (g, heads) = nw::align(&tails);
                let mut alignment = nw::poa_to_strings(&g, &heads)
                    .values()
                    .cloned()
                    .collect::<Vec<_>>();
                if false {
                    // Differentiate between tail of shorter alignments and actual indels
                    // Left tail
                    alignment.iter_mut().for_each(|a| {
                        for pos in a.iter_mut() {
                            if pos == INDEL {
                                *pos = "".to_string();
                            } else {
                                break;
                            }
                        }
                    });
                    // Right tail
                    alignment.iter_mut().for_each(|a| {
                        for pos in a.iter_mut().rev() {
                            if pos == INDEL {
                                *pos = "".to_string();
                            } else {
                                break;
                            }
                        }
                    });
                }

                Some(
                    (0..alignment[0].len())
                        .map(|i| {
                            // for kk in alignment.iter().map(|a| &a[i]) {
                            //     let kkk = kk.chars().take(18).collect::<String>();
                            //     print!(" {:<18} ", kkk);
                            // }
                            let count = if false {
                                alignment.iter().filter(|a| a[i] != EMPTY).count() as f32
                            } else {
                                alignment.len() as f32
                            };

                            let mut counts: HashMap<String, i32> = HashMap::new();
                            alignment
                                .iter()
                                .map(|a| &a[i])
                                .for_each(|g| *counts.entry(g.clone()).or_insert(0) += 1);
                            PolyGene {
                                genes: counts
                                    .into_iter()
                                    .filter_map(|(name, v)| {
                                        // if DROP_EMPTY && (name == EMPTY || name == INDEL) {
                                        //     None
                                        // } else {
                                        let name = if name != MARKER {
                                            &name
                                        } else {
                                            &common_ancestral
                                        };
                                        Some((
                                            Gene {
                                                name: name.clone(),
                                                color: if name == EMPTY || name == INDEL {
                                                    "#ccc".to_string()
                                                } else {
                                                    colormap
                                                        .get(&name.clone())
                                                        .map(|c| c.to_hex_string())
                                                        .unwrap_or_else(|| "#aaa".to_string())
                                                },
                                            },
                                            v as f32 / count,
                                        ))
                                        // }
                                    })
                                    .collect::<Vec<(Gene, f32)>>(),
                            }
                        })
                        .collect::<Vec<PolyGene>>(),
                )
            } else {
                None
            }
        };

        let ((species, chr, gene, ancestral, t_len), (lefts, rights)) =
            if let Some(name) = &tree[node].name {
                let gene_name = name.split('#').next().unwrap();
                if let Some(DbGene {
                    ancestral,
                    species,
                    chr,
                    t_len,
                    left_tail,
                    right_tail,
                    ..
                }) = genes.get(gene_name)
                {
                    common_ancestral = ancestral.to_owned();
                    let (proto_lefts, proto_rights) = (left_tail, right_tail);
                    let (lefts, rights) = if true {
                        (proto_rights, proto_lefts)
                    } else {
                        (proto_lefts, proto_rights)
                    };
                    (
                        (
                            species.to_owned(),
                            chr.to_owned(),
                            gene_name.to_owned(),
                            ancestral.to_owned(),
                            *t_len,
                        ),
                        (
                            lefts
                                .iter()
                                .rev()
                                .map(|g| Gene {
                                    name: g.clone(),
                                    color: if g == EMPTY {
                                        "#111".into()
                                    } else {
                                        colormap
                                            .get(&g.clone())
                                            .map(|c| c.to_hex_string())
                                            .unwrap_or_else(|| "#aaa".to_string())
                                    },
                                })
                                .collect::<Vec<_>>(),
                            rights
                                .iter()
                                .map(|g| Gene {
                                    name: g.clone(),
                                    color: if g == EMPTY {
                                        "#111".into()
                                    } else {
                                        colormap
                                            .get(&g.clone())
                                            .map(|c| c.to_hex_string())
                                            .unwrap_or_else(|| "#aaa".to_string())
                                    },
                                })
                                .collect::<Vec<_>>(),
                        ),
                    )
                } else {
                    eprintln!("{} -- {} not found", name, gene_name);
                    (
                        (
                            String::new(),
                            String::new(),
                            String::new(),
                            String::new(),
                            0,
                        ),
                        (vec![], vec![]),
                    )
                }
            } else {
                (
                    (
                        String::new(),
                        String::new(),
                        String::new(),
                        String::new(),
                        0,
                    ),
                    (vec![], vec![]),
                )
            };

        let color = name2color(&species).to_hex_string();
        HtmlNode {
            species,
            chr,
            gene,
            ancestral,
            t_len,
            color,
            children: tree[node]
                .children
                .as_ref()
                .map(|children| {
                    children
                        .iter()
                        .map(|n| process(tree, *n, genes, colormap))
                        .collect()
                })
                .unwrap_or_default(),
            isDuplication: tree[node].is_duplication(),
            confidence: tree[node]
                .data
                .get("DCS")
                .and_then(|x| x.parse::<f32>().ok())
                .unwrap_or(0.0),
            repr: Landscape {
                lefts,
                rights,
                me: Gene {
                    color: gene2color(&common_ancestral).to_hex_string(),
                    name: common_ancestral,
                },
            },
            clustered,
        }
    }

    process(tree, 0, genes, colormap)
}

struct Settings {
    colorize_per_duplication: bool,
    colorize_all: bool,
}

fn process_file(filename: &str, db_filename: &str, graph_type: &str, settings: Settings) {
    println!("Processing {}", filename);
    let t = Tree::from_filename(filename).unwrap();
    let mut db =
        Connection::open_with_flags(db_filename, OpenFlags::SQLITE_OPEN_READ_ONLY).unwrap();
    let genes = make_genes_cache(&t, &mut db);
    let colormap = if settings.colorize_per_duplication {
        make_colormap_per_duplication(&t, &genes, settings.colorize_all)
    } else {
        make_colormap(&t, &genes)
    };

    let depth = BRANCH_WIDTH * (t.topological_depth().1 + 1.);
    let longest_name = t
        .leaf_names()
        .iter()
        .map(|(_, name)| name.map_or(0, |s| s.len()))
        .max()
        .unwrap() as f32
        * FONT_SIZE;
    let xlabels = 0.85 * (10. + depth + longest_name + 50.);
    let width = xlabels + (2. * WINDOW as f32 + 1.) * (GENE_WIDTH + GENE_SPACING) + 60.;
    let mut svg = SvgDrawing::new();
    svg.text()
        .pos(FONT_SIZE, FONT_SIZE)
        .text(Path::new(filename).file_stem().unwrap().to_str().unwrap());

    let t_lens = t
        .leaves()
        .filter_map(|n| {
            t[n].name
                .as_ref()
                .and_then(|name| name.split('#').next())
                .and_then(|gene_name| genes.get(gene_name))
        })
        .map(|x| x.t_len)
        .filter(|&x| x >= 0)
        .collect::<Vec<_>>();

    let min_t_len = t_lens.iter().min().unwrap_or(&0);
    let max_t_len = t_lens.iter().max().unwrap_or(&0);

    match graph_type {
        "flat" => {
            draw_background(
                &mut svg, &genes, depth, &t, &t[0], 10.0, 50.0, xlabels, width,
            );
            let mut links = Vec::new();
            draw_tree(
                &mut svg, &genes, &colormap, depth, &t, &t[0], 10.0, 50.0, xlabels, width,
                &mut links,
            );
            draw_links(&mut svg, &links, 50.0, xlabels);
            svg.auto_fit();
            let mut out = File::create(&format!("{}.svg", filename)).unwrap();
            out.write_all(svg.render_svg().as_bytes()).unwrap();
        }
        "html" => {
            use tera::{Context, Tera};

            let mut tera = match Tera::new("templates/*.{html,css,js}") {
                Ok(t) => t,
                Err(e) => {
                    println!("Parsing error(s): {}", e);
                    ::std::process::exit(1);
                }
            };
            let mut out = File::create(&format!("{}.html", &filename)).unwrap();
            let mut context = Context::new();
            tera.autoescape_on(vec![]);
            context.insert(
                "title",
                &Path::new(filename).file_name().unwrap().to_str().unwrap(),
            );
            context.insert(
                "comment",
                &format!("Transcript length: {} to {}", min_t_len, max_t_len),
            );
            context.insert("min_t_len", min_t_len);
            context.insert("max_t_len", max_t_len);
            context.insert(
                "data",
                &serde_json::to_string_pretty(&draw_html(&t, &genes, &colormap)).unwrap(),
            );
            tera.render_to("genominicus.html", &context, &mut out)
                .unwrap();
        }
        _ => unimplemented!(),
    };
}

fn main() {
    let args = App::new("Genominicus")
        .version(clap::crate_version!())
        .author(clap::crate_authors!())
        .arg(
            Arg::with_name("type")
                .help("The graph type to use")
                .required(true)
                .possible_values(&["flat", "condensed", "html"])
                .index(1),
        )
        .arg(
            Arg::with_name("database")
                .help("The database to use")
                .required(true)
                .index(2),
        )
        .arg(
            Arg::with_name("FILE")
                .help("Sets the input file to use")
                .required(true)
                .multiple(true),
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
        .get_matches();

    for filename in values_t!(args, "FILE", String).unwrap().iter() {
        process_file(
            filename,
            &value_t!(args, "database", String).unwrap(),
            &value_t!(args, "type", String).unwrap(),
            Settings {
                colorize_per_duplication: args.is_present("colorize_per_duplication"),
                colorize_all: args.is_present("colorize_all"),
            },
        );
    }
}
