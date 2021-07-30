use crate::nhx::*;
use clap::*;
use colorsys::{Hsl, Rgb};
use petgraph::{
    algo::dijkstra,
    graph::{DiGraph, NodeIndex, UnGraph},
};
use petgraph::{algo::toposort, Direction};
use petgraph::{
    dot::{Config, Dot},
    visit::EdgeRef,
};
use poa::POA;
use rusqlite::*;
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

const ANCESTRAL_QUERY: &str =
    concat!("select gene, ancestral, species, chr, start, direction from genomes where protein=?",);
const LEFTS_QUERY: &str = "select ancestral, direction from genomes where species=? and chr=? and start<? order by start desc limit ?";
const RIGHTS_QUERY: &str = "select ancestral, direction from genomes where species=? and chr=? and start>? order by start asc limit ?";

fn left_tail(
    db: &mut Connection,
    species: &str,
    chr: &str,
    pos: i32,
    span: i64,
) -> Vec<(String, char)> {
    let mut r = Vec::new();
    for g in db
        .prepare(LEFTS_QUERY)
        .unwrap()
        .query_map(params![&species, &chr, pos, WINDOW], |row| {
            let ancestral: String = row.get("ancestral").unwrap();
            let direction: char = row
                .get::<_, String>("direction")
                .unwrap()
                .chars()
                .next()
                .unwrap();
            Ok((ancestral, direction))
        })
        .unwrap()
    {
        r.push(g.unwrap());
    }
    r
}

fn right_tail(
    db: &mut Connection,
    species: &str,
    chr: &str,
    pos: i32,
    span: i64,
) -> Vec<(String, char)> {
    let mut r = Vec::new();
    for g in db
        .prepare(RIGHTS_QUERY)
        .unwrap()
        .query_map(params![&species, &chr, pos, WINDOW], |row| {
            let ancestral: String = row.get("ancestral").unwrap();
            let direction: char = row
                .get::<_, String>("direction")
                .unwrap()
                .chars()
                .next()
                .unwrap();
            Ok((ancestral, direction))
        })
        .unwrap()
    {
        r.push(g.unwrap());
    }
    r
}

fn tails(
    db: &mut Connection,
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
    let mut db = Connection::open_with_flags(
        "/home/franklin/work/duplications/data/db.sqlite",
        OpenFlags::SQLITE_OPEN_READ_ONLY,
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
                if let Ok((gene, ancestral, species, chr, pos, direction)) =
                    db.query_row(ANCESTRAL_QUERY, &[&protein_name], |r| {
                        let gene: String = r.get("gene").unwrap();
                        let ancestral: String = r.get("ancestral").unwrap();
                        let species: String = r.get("species").unwrap();
                        let chr: String = r.get("chr").unwrap();
                        let pos: i32 = r.get("start").unwrap();
                        let direction: String = r.get("direction").unwrap();
                        Ok((gene, ancestral, species, chr, pos, direction))
                    })
                {
                    let proto_lefts = left_tail(&mut db, &species, &chr, pos, WINDOW);
                    let proto_rights = right_tail(&mut db, &species, &chr, pos, WINDOW);
                    let (lefts, rights) = if true {
                        (proto_lefts, proto_rights)
                    } else {
                        (proto_rights, proto_lefts)
                    };

                    // Gene/protein name
                    svg.text()
                        .pos(depth, y + 5.)
                        .text(format!("{} {}/{}", gene, species, chr))
                        .style(|s| s.fill_color(name2color(species)));

                    // Left tail
                    let xbase = xlabels + (WINDOW as f32 - 1.) * (GENE_WIDTH + GENE_SPACING);
                    for (k, g) in lefts.iter().enumerate() {
                        let xstart = xbase - (k as f32) * (GENE_WIDTH + GENE_SPACING);
                        draw_gene(svg, xstart, y, g.1 == '+', gene2color(&g.0), &g.0);
                    }

                    // The Gene
                    draw_gene(
                        svg,
                        xlabels + WINDOW as f32 * (GENE_WIDTH + GENE_SPACING),
                        y,
                        true,
                        gene2color(&ancestral),
                        &ancestral,
                    )
                    .style(|s| {
                        s.stroke_width(2.)
                            .stroke_color(StyleColor::Percent(0.1, 0.1, 0.1))
                    });

                    // Right tail
                    let xbase = xlabels + (WINDOW as f32 + 1.) * (GENE_WIDTH + GENE_SPACING);
                    for (k, g) in rights.iter().enumerate() {
                        let xstart = xbase + (k as f32) * (GENE_WIDTH + GENE_SPACING);
                        draw_gene(svg, xstart, y, g.1 == '+', gene2color(&g.0), &g.0);
                    }
                    links.push((
                        lefts.iter().map(|x| x.0.clone()).collect(),
                        ancestral.into(),
                        rights.iter().map(|x| x.0.clone()).collect(),
                    ));
                } else {
                    // The node was not found in the database
                    eprintln!("{} -- {} not found", name, protein_name);
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
        let mut db = Connection::open_with_flags(
            "/home/franklin/work/duplications/data/db.sqlite",
            OpenFlags::SQLITE_OPEN_READ_ONLY,
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
                    if let Ok((gene, ancestral, species, chr, pos, direction)) =
                        db.query_row(ANCESTRAL_QUERY, &[&protein_name], |r| {
                            let gene: String = r.get("gene").unwrap();
                            let ancestral: String = r.get("ancestral").unwrap();
                            let species: String = r.get("species").unwrap();
                            let chr: String = r.get("chr").unwrap();
                            let pos: i32 = r.get("start").unwrap();
                            let direction: char = r
                                .get::<_, String>("direction")
                                .unwrap()
                                .chars()
                                .next()
                                .unwrap();
                            Ok((gene, ancestral, species, chr, pos, direction))
                        })
                    {
                        let (left, right) = tails(&mut db, &species, &chr, pos, WINDOW);
                        Some(
                            left.iter()
                                .map(|g| &g.0)
                                .rev()
                                .chain(["=======FINAL=======".to_owned()].iter())
                                .chain(right.iter().map(|g| &g.0))
                                .cloned()
                                .collect::<Vec<String>>(),
                        )
                    } else {
                        eprintln!("{} -- {} not found", name, protein_name);
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
            if node.id != 57907890890 {
                let _ = File::create(&format!("poa-{}.dot", node.id))
                    .unwrap()
                    .write_all(format!("{:?}", Dot::new(&g)).as_bytes());

                let centromere = g
                    .node_indices()
                    .filter(|n| g[*n].nucs.values().any(|v| v == "=======FINAL======="))
                    .next()
                    .unwrap();
                dbg!(centromere);
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
                let mut sources: Vec<NodeIndex> = g
                    .node_indices()
                    .filter(|n| g.neighbors_directed(*n, Incoming).next().is_none())
                    .collect();
                sources.sort_by(|x, y| furthest_sink[y].cmp(&furthest_sink[x]));

                fn rank(
                    n: NodeIndex,
                    x: i32,
                    y: i32,

                    g: &poa::POAGraph,
                    furthest_sink: &HashMap<NodeIndex, i32>,
                    ranks: &mut HashMap<NodeIndex, (i32, i32)>,
                    direction: bool,
                ) {
                    if !ranks.contains_key(&n) {
                        ranks.insert(n, (x, y));
                        let mut children: Vec<_> = g.neighbors_directed(n, if direction { Outgoing } else { Incoming }).collect();
                        children.sort_by(|y, x| furthest_sink[y].cmp(&furthest_sink[x]));

                        let mut top;
                        for (k, &c) in children.iter().enumerate() {
                            top = ranks.values().map(|v| v.1).max().unwrap_or(0) as i32;
                            let next_x = if direction { x+1 } else { x - 1 };
                            if k > 0 {
                                rank(c, next_x, top + 1, &g, furthest_sink, ranks, direction);
                            } else {
                                rank(c, next_x, y, &g, furthest_sink, ranks, direction);
                            }
                        }
                    }
                }

                let mut ranks_left = HashMap::<NodeIndex, (i32, i32)>::new();
                let mut ranks_right = HashMap::<NodeIndex, (i32, i32)>::new();
                // for source in sources {
                    // let top = ranks.values().map(|v| v.1).max().unwrap_or(1);
                rank(centromere, 20, 0, &g, &furthest_sink, &mut ranks_left, false);
                ranks_left.remove(&centromere);
                rank(centromere, 20, 0, &g, &furthest_sink, &mut ranks_right, true);
                // }
                // for (k, n) in toposort(&g, None).unwrap().iter().enumerate() {
                //     ranks.insert(*n, (k, 0));
                // }
                let ranks: HashMap::<NodeIndex, (i32, i32)> = ranks_left.into_iter().chain(ranks_right).collect();

                let xbase = xlabels;
                for e in g.edge_references() {
                    let (_x1, _y1) = *ranks.get(&e.source()).unwrap_or(&(0, 0));
                    let (_x2, _y2) = *ranks.get(&e.target()).unwrap_or(&(0, 0));
                    let ((x1, y1), (x2, y2)) = if _x1 < _x2 {
                        ((_x1 as f32, _y1 as f32), (_x2 as f32, _y2 as f32))
                    } else {
                        ((_x2 as f32, _y2 as f32), (_x1 as f32, _y1 as f32))
                    };

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

                    // if (x1 - x2).abs() > 1. {
                    //     let dy = 15. + 20. * (x1 - x2).abs().log10();
                    //     svg.polygon()
                    //         .from_coords([
                    //             (xbase + x1 as f32 * 20. + 5., y + y1 as f32 * 20. + 5.),
                    //             (xbase + x1 as f32 * 20. + 5., y + y1 as f32 * 20. + 5. + dy),
                    //             (xbase + x2 as f32 * 20. + 5., y + y2 as f32 * 20. + 5. + dy),
                    //             (xbase + x2 as f32 * 20. + 5., y + y2 as f32 * 20. + 5.),
                    //         ])
                    //         .style(|s| {
                    //             s.stroke_color(StyleColor::RGB(0, 0, 0))
                    //                 .stroke_width(e.weight().len() as f32)
                    //                 .fill_opacity(0.)
                    //         });
                    // }
                }
                for (n, (xi, yi)) in ranks.iter() {
                    svg.polygon()
                        .from_pos_dims(xbase + *xi as f32 * 20., y + *yi as f32 * 20., 10., 10.)
                        .set_hover(g[*n].nucs.values().next().unwrap())
                        .style(|s| s.fill_color(gene2color(g[*n].nucs.values().next().unwrap())));
                }
            }

            y += 300.;

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
