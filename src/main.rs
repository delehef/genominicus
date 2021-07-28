use crate::nhx::*;
use bimap::BiMap;
use clap::*;
use colorsys::{Hsl, Rgb};
use petgraph::Direction;
use petgraph::{
    algo::dijkstra,
    graph::{DiGraph, NodeIndex, UnGraph},
};
use petgraph::{
    dot::{Config, Dot},
    visit::EdgeRef,
};
use poa::POA;
use postgres::{Client, NoTls};
use std::collections::HashSet;
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use std::{collections::HashMap, usize};
use svarog::*;
use Direction::{Incoming, Outgoing};
mod nhx;
mod nw;
mod poa;

pub type SeqID = usize;
pub type Nucleotide = String;
pub type Sequence = Vec<Nucleotide>;
pub type Sequences = HashMap<SeqID, Sequence>;

fn nuc_to_str(nuc: &Nucleotide) -> String {
    nuc.to_string()
}

const WINDOW: i64 = 15;
const GENE_WIDTH: f32 = 15.;
const GENE_SPACING: f32 = 5.;
const BRANCH_WIDTH: f32 = 20.;
const FONT_SIZE: f32 = 10.;

const ANCESTRAL_QUERY: &str = concat!(
    "select name, ancestral, annotations.species, direction, chr, start from annotations ",
    "left join mapping on name=modern where name=",
    "(select parent from annotations where name=(select parent from annotations where name=$1 limit 1))"
);
const LEFTS_QUERY: &str = concat!(
    "select coalesce(ancestral, concat($1, ':', $2, ':', $3)) as ancestral, direction from ",
    "(select * from annotations left join mapping on name=modern where ",
    "annotations.species=$1 and chr=$2 and type='gene' and start<$3 order by start desc) as x limit $4"
);
const RIGHTS_QUERY: &str = concat!(
    "select coalesce(ancestral, concat($1, ':', $2, ':', $3)) as ancestral, direction from ",
    "(select * from annotations left join mapping on name=modern where ",
    "annotations.species=$1 and chr=$2 and type='gene' and start>$3 order by start asc) as x limit $4"
);

fn left_tail(
    db: &mut Client,
    species: &str,
    chr: &str,
    pos: i32,
    span: i64,
) -> Vec<(String, char)> {
    db.query(LEFTS_QUERY, &[&species, &chr, &pos, &WINDOW])
        .unwrap()
        .into_iter()
        .map(|row| {
            let ancestral: String = row.get("ancestral");
            let direction: char = row.get::<_, String>("direction").chars().next().unwrap();
            (ancestral, direction)
        })
        .collect::<Vec<_>>()
}

fn right_tail(
    db: &mut Client,
    species: &str,
    chr: &str,
    pos: i32,
    span: i64,
) -> Vec<(String, char)> {
    db.query(RIGHTS_QUERY, &[&species, &chr, &pos, &WINDOW])
        .unwrap()
        .into_iter()
        .map(|row| {
            let ancestral: String = row.get("ancestral");
            let direction: char = row.get::<_, String>("direction").chars().next().unwrap();
            (ancestral, direction)
        })
        .collect::<Vec<_>>()
}

fn tails(
    db: &mut Client,
    species: &str,
    chr: &str,
    pos: i32,
    span: i64,
) -> (Vec<(String, char)>, Vec<(String, char)>) {
    (
        left_tail(db, species, chr, pos, span),
        right_tail(db, species, chr, pos, span),
    )
}

fn draw_background(
    svg: &mut SvgDrawing,
    depth: f32,
    tree: &Tree,
    n_: usize,
    xoffset: f32,
    yoffset: f32,
    xlabels: f32,
    width: f32,
) -> f32 {
    let mut y = yoffset;

    for &n in tree[n_].children().iter() {
        let child = &tree[n];
        let new_y = if child.is_leaf() {
            y + 20.
        } else {
            draw_background(
                svg,
                depth,
                tree,
                n,
                xoffset + BRANCH_WIDTH,
                y,
                xlabels,
                width,
            )
        };

        if tree[n_].is_duplication() {
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

fn name2color<S: AsRef<str>>(name: S) -> StyleColor {
    let bytes: [u8; 16] = md5::compute(name.as_ref().as_bytes()).into();
    let rgb = Rgb::from((bytes[0] as f32, bytes[1] as f32, bytes[2] as f32));

    let mut hsl: Hsl = rgb.into();
    hsl.set_lightness(hsl.lightness().clamp(30., 40.));

    let rgb: Rgb = hsl.into();
    StyleColor::Percent(
        rgb.red() as f32 / 255.,
        rgb.green() as f32 / 255.,
        rgb.blue() as f32 / 255.,
    )
}
fn gene2color<S: AsRef<str>>(name: S) -> StyleColor {
    let bytes: [u8; 16] = md5::compute(name.as_ref().as_bytes()).into();
    let r = (bytes[0] as f32 / 255.).clamp(0.1, 0.9);
    let g = (bytes[1] as f32 / 255.).clamp(0.1, 0.9);
    let b = (bytes[2] as f32 / 255.).clamp(0.1, 0.9);
    StyleColor::Percent(r, g, b)
}

fn draw_gene<'a>(
    svg: &'a mut SvgDrawing,
    x: f32,
    y: f32,
    right: bool,
    color: StyleColor,
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
                s.fill_color(color)
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
                s.fill_color(color)
                    .stroke_width(0.5)
                    .stroke_color(StyleColor::Percent(0.2, 0.2, 0.2))
            })
    }
}

fn draw_tree(
    svg: &mut SvgDrawing,
    depth: f32,
    tree: &Tree,
    node: &Node,
    xoffset: f32,
    yoffset: f32,
    xlabels: f32,
    width: f32,
    links: &mut Vec<(Vec<String>, String, Vec<String>)>,
) -> f32 {
    let mut pg = Client::connect(
        "host=localhost user=franklin dbname=duplications password='notanumber'",
        NoTls,
    )
    .unwrap();
    let mut y = yoffset;
    let mut old_y = 0.;
    for (i, child) in node.children().iter().map(|i| &tree[*i]).enumerate() {
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
                .start(xlabels - 5., y)
                .end(
                    xlabels + (GENE_WIDTH + GENE_SPACING) * (2. * WINDOW as f32 + 1.)
                        - GENE_SPACING
                        + 5.,
                    y,
                )
                .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5));

            if let Some(name) = &child.name {
                let protein_name = name.split('_').next().unwrap();
                if let Ok(r) = pg.query_one(ANCESTRAL_QUERY, &[&protein_name]) {
                    let gene_name: &str = r.get("name");
                    let ancestral_name: &str = r.get("ancestral");
                    let species: &str = r.get("species");
                    let chr: &str = r.get("chr");
                    let pos: i32 = r.get("start");
                    let direction: &str = r.get("direction");

                    let proto_lefts = pg
                        .query(LEFTS_QUERY, &[&species, &chr, &pos, &WINDOW])
                        .unwrap()
                        .into_iter()
                        .map(|row| {
                            let ancestral: String = row.get("ancestral");
                            let direction: String = row.get("direction");
                            (ancestral, direction)
                        })
                        .collect::<Vec<_>>();
                    let proto_rights = pg
                        .query(RIGHTS_QUERY, &[&species, &chr, &pos, &WINDOW])
                        .unwrap()
                        .into_iter()
                        .map(|row| {
                            let ancestral: String = row.get("ancestral");
                            let direction: String = row.get("direction");
                            (ancestral, direction)
                        })
                        .collect::<Vec<_>>();
                    let (lefts, rights) = if true {
                        (proto_lefts, proto_rights)
                    } else {
                        (proto_rights, proto_lefts)
                    };

                    // Gene/protein name
                    svg.text()
                        .pos(depth, y + 5.)
                        .text(format!("{} {} {}/{}", child.id, protein_name, species, chr))
                        .style(|s| s.fill_color(name2color(species)));

                    // Left tail
                    let xbase = xlabels + (WINDOW as f32 - 1.) * (GENE_WIDTH + GENE_SPACING);
                    for (k, g) in lefts.iter().enumerate() {
                        let xstart = xbase - (k as f32) * (GENE_WIDTH + GENE_SPACING);
                        draw_gene(svg, xstart, y, g.1 == "+", gene2color(&g.0), &g.0);
                    }

                    // The Gene
                    draw_gene(
                        svg,
                        xlabels + WINDOW as f32 * (GENE_WIDTH + GENE_SPACING),
                        y,
                        true,
                        gene2color(&ancestral_name),
                        &ancestral_name,
                    )
                    .style(|s| {
                        s.stroke_width(2.)
                            .stroke_color(StyleColor::Percent(0.1, 0.1, 0.1))
                    });

                    // Right tail
                    let xbase = xlabels + (WINDOW as f32 + 1.) * (GENE_WIDTH + GENE_SPACING);
                    for (k, g) in rights.iter().enumerate() {
                        let xstart = xbase + (k as f32) * (GENE_WIDTH + GENE_SPACING);
                        draw_gene(svg, xstart, y, g.1 == "+", gene2color(&g.0), &g.0);
                    }
                    links.push((
                        lefts.iter().map(|x| x.0.clone()).collect(),
                        ancestral_name.into(),
                        rights.iter().map(|x| x.0.clone()).collect(),
                    ));
                } else {
                    // The node was not found in the database
                    eprintln!("{} not found", name);
                    links.push((Vec::new(), name.into(), Vec::new()));
                }
            } else { // The node does not have a name
            }
            y += 20.;
        } else {
            svg.line()
                .from_coords(xoffset, y, xoffset + BRANCH_WIDTH, y)
                .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5));
            y = draw_tree(
                svg,
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
        svg.polygon()
            .from_pos_dims(xoffset - 3., yoffset - 3., 6., 6.)
            .style(|s| s.fill_color(StyleColor::Percent(0.8, 0., 0.)));
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
                    .start(x1, y + 5.)
                    .end(x2, y + 20. - 5.)
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
                    .start(x1, y + 5.)
                    .end(x2, y + 20. - 5.)
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

fn draw_clustered(
    svg: &mut SvgDrawing,
    depth: f32,
    tree: &Tree,
    node: &Node,
    xoffset: f32,
    yoffset: f32,
    xlabels: f32,
    width: f32,
) -> f32 {
    let mut y = yoffset;
    if node.is_duplication() {
        assert!(node.children().len() == 2);

        svg.line()
            .from_coords(xoffset, y, xoffset + 2. * BRANCH_WIDTH, y)
            .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5));

        // Left arm
        let old_y = y;
        y = draw_clustered(
            svg,
            depth,
            tree,
            &tree[node.children()[0]],
            xoffset + 2. * BRANCH_WIDTH,
            y,
            xlabels,
            width,
        );

        // Vertical rake
        svg.line()
            .from_coords(xoffset + BRANCH_WIDTH, old_y, xoffset + BRANCH_WIDTH, y)
            .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5));

        // Right arm
        svg.line()
            .from_coords(xoffset + BRANCH_WIDTH, y, xoffset + 2. * BRANCH_WIDTH, y)
            .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5));
        y = draw_clustered(
            svg,
            depth,
            tree,
            &tree[node.children()[1]],
            xoffset + 2. * BRANCH_WIDTH,
            y,
            xlabels,
            width,
        );
        // The little red square
        svg.polygon()
            .from_pos_dims(BRANCH_WIDTH + xoffset - 3., yoffset - 3., 6., 6.)
            .style(|s| s.fill_color(StyleColor::Percent(0.8, 0., 0.)));
    } else {
        let leaves = tree.non_d_descendants(node.id);
        let dups = tree.d_descendants(node.id);
        let mut pg = Client::connect(
            "host=localhost user=franklin dbname=duplications password='notanumber'",
            NoTls,
        )
        .unwrap();

        let old_y = y;
        for &child in dups.iter() {
            svg.line()
                .from_coords(xoffset, y, xoffset + BRANCH_WIDTH, y)
                .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5));
            y = draw_clustered(
                svg,
                depth,
                tree,
                &tree[child],
                xoffset + BRANCH_WIDTH,
                y,
                xlabels,
                width,
            );
        }

        let mut tails = leaves
            .iter()
            .filter_map(|&l| {
                if let Some(name) = &tree[l].name {
                    let protein_name = name.split('_').next().unwrap();
                    if let Ok(r) = pg.query_one(ANCESTRAL_QUERY, &[&protein_name]) {
                        let gene_name: &str = r.get("name");
                        let ancestral_name: &str = r.get("ancestral");
                        let species: &str = r.get("species");
                        let chr: &str = r.get("chr");
                        let pos: i32 = r.get("start");
                        let direction: char = r.get::<_, &str>("direction").chars().next().unwrap();
                        let (left, right) = tails(&mut pg, species, chr, pos, WINDOW);
                        Some(
                            left.iter()
                                .map(|g| &g.0)
                                .rev()
                                .chain(["=======FINAL=======".to_owned()].iter())
                                .chain(right.iter().map(|g| &g.0))
                                .cloned()
                                .collect::<Vec<String>>(),
                        )

                        // Some((left, right))
                    } else {
                        // The node was not found in the database
                        eprintln!("{} not found", name);
                        None
                    }
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        if tails.is_empty() {
            eprintln!("EMPTY");
        } else {
            let h_tails = tails
                .into_iter()
                .enumerate()
                .collect::<HashMap<usize, Vec<String>>>();
            let (g, starts) = nw::align(&h_tails);
            let labels: HashMap<usize, String> = leaves
                .iter()
                .filter_map(|&l| tree[l].name.clone())
                .enumerate()
                .collect();
            let xbase = xlabels;
            if node.id == 5 {
                let _ = File::create(&format!("poa-{}.dot", node.id))
                    .unwrap()
                    .write_all(format!("{:?}", Dot::new(&g)).as_bytes());

                let mut ranks = HashMap::<NodeIndex, (usize, usize)>::new();
                let sources: Vec<NodeIndex> = g
                    .node_indices()
                    .filter(|n| g.neighbors_directed(*n, Incoming).next().is_none())
                    .collect();
                let sinks: Vec<NodeIndex> = g
                    .node_indices()
                    .filter(|n| g.neighbors_directed(*n, Outgoing).next().is_none())
                    .collect();
                let furthest_sink: HashMap<NodeIndex, i32> = g
                    .node_indices()
                    .map(|s| {
                        (
                            s,
                            *dijkstra(&g, s, None, |_| 1)
                                .iter()
                                .filter(|(n, _d)| sinks.contains(n))
                                .max_by(|x, y| x.1.cmp(&y.1))
                                .unwrap()
                                .1,
                        )
                    })
                    .collect();
                fn rank(
                    n: NodeIndex,
                    x: usize,
                    y: usize,
                    _max_y: usize,

                    g: &poa::POAGraph,
                    furthest_sink: &HashMap<NodeIndex, i32>,
                    ranks: &mut HashMap<NodeIndex, (usize, usize)>,
                ) -> usize {
                    let mut max_y = _max_y;
                    if !ranks.contains_key(&n) {
                        ranks.insert(n, (x, y));
                        let mut children: Vec<_> = g.neighbors_directed(n, Outgoing).collect();
                        children.sort_by(|x, y| furthest_sink[y].cmp(&furthest_sink[x]));
                        let mut children = children.iter();
                        if let Some(first_child) = children.next() {
                            rank(*first_child, x + 1, y, max_y, &g, furthest_sink, ranks);
                        }

                        let mut offset = 0;
                        for &c in children {
                            rank(c, x + 1, y + offset, max_y, &g, furthest_sink, ranks);
                            offset += cw;
                        }
                        current_width + offset
                    } else {
                        current_width
                    }
                    max_y
                }
                let mut offset = 0;
                for source in sources {
                    offset = rank(source, 0, offset, 1, &g, &furthest_sink, &mut ranks);
                }
                dbg!(&ranks);
                let xbase = xlabels;
                for e in g.edge_references() {
                    let (_x1, _y1) = *ranks.get(&e.source()).unwrap();
                    let (_x2, _y2) = *ranks.get(&e.target()).unwrap();
                    let x1 = std::cmp::min(_x1, _x2) as f32;
                    let x2 = std::cmp::max(_x1, _x2) as f32;
                    let y1 = std::cmp::min(_y1, _y2) as f32;
                    let y2 = std::cmp::max(_y1, _y2) as f32;

                    // if x1 - x2 == 1. {
                    svg.line()
                        .from_coords(
                            xbase + x1 as f32 * 20. + 5.,
                            y + y1 as f32 * 20. + 5.,
                            xbase + x2 as f32 * 20. + 5.,
                            y + y2 as f32 * 20. + 5.,
                        )
                        .style(|s| {
                            s.stroke_color(StyleColor::RGB(0, 0, 0))
                                .stroke_width(e.weight().len() as f32)
                        });
                    // } else {
                    //     let dx = x1 as f32 - x2 as f32;
                    //     let mid_y = y + if y2 != y1 {
                    //         (y2 as f32 - y1 as f32)/2.
                    //     } else {
                    //         y2 as f32 + 10.
                    //     };
                    //     svg.line()
                    //         .from_coords(
                    //             xbase + x1 as f32 * 20. + 5.,
                    //             y + y1 as f32 * 20. + 5.,
                    //             xbase + x2 as f32 * 20. + 5.,
                    //             y + y2 as f32 * 20. + 5.,
                    //         )
                    //         .style(|s| {
                    //             s.stroke_color(StyleColor::RGB(0, 0, 0))
                    //                 .stroke_width(e.weight().len() as f32)
                    //         });
                    //     // svg.polygon()
                    //     //     .from_coords([
                    //     //         (
                    //     //             xbase + x1 as f32 * 20. + 5.,
                    //     //             y + y1 as f32 * 20. + 5.,
                    //     //         ),

                    //     //         // (
                    //     //         //     xbase + x1 as f32 * 20. + 5.,
                    //     //         //     mid_y
                    //     //         // ),
                    //     //         // (
                    //     //         //     xbase + x2 as f32 * 20. + 5.,
                    //     //         //     mid_y
                    //     //         // ),

                    //     //         (
                    //     //             xbase + x2 as f32 * 20. + 5.,
                    //     //             y + y2 as f32 * 20. + 5.,
                    //     //         ),
                    //     //     ])
                    //     //     .style(|s| {
                    //     //         s.stroke_color(StyleColor::RGB(0, 0, 0))
                    //     //             .stroke_width(e.weight().len() as f32)
                    //     //             .fill_opacity(0.)
                    //     //     });
                    // }
                }
                for (n, (xi, yi)) in ranks.iter() {
                    svg.polygon()
                        .from_pos_dims(xbase + *xi as f32 * 20., y + *yi as f32 * 20., 10., 10.)
                        .set_hover(g[*n].nucs.values().next().unwrap())
                        .style(|s| s.fill_color(gene2color(g[*n].nucs.values().next().unwrap())));
                }
            }

            y += 20.;

            svg.line()
                .from_coords(xoffset, old_y, xoffset, y - 20.)
                .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5));
        }
    };
    y
}

fn process_file(filename: &str) {
    println!("Processing {}", filename);
    let t = Tree::from_filename(filename).unwrap();

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

    // draw_background(&mut svg, depth, &t, 0, 10.0, 50.0, xlabels, width);
    // let mut links = Vec::new();
    // draw_tree(
    //     &mut svg, depth, &t, &t[0], 10.0, 50.0, xlabels, width, &mut links,
    // );
    // draw_links(&mut svg, &links, 50.0, xlabels);

    draw_clustered(&mut svg, depth, &t, &t[0], 10.0, 50.0, xlabels, width);

    svg.auto_fit();

    let mut out = File::create(&format!("{}.svg", filename)).unwrap();
    out.write_all(svg.render_svg().as_bytes()).unwrap();
}

fn main() {
    let args = App::new("Genominicus")
        .version(clap::crate_version!())
        .author(clap::crate_authors!())
        .arg(
            Arg::with_name("FILE")
                .help("Sets the input file to use")
                .required(true)
                .multiple(true),
        )
        .get_matches();

    for filename in values_t!(args, "FILE", String).unwrap().iter() {
        process_file(filename);
    }
}
