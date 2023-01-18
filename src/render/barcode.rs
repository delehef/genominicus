use crate::utils::*;
use newick::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::prelude::*;
use svarog::*;

const K: f32 = FONT_SIZE;

fn draw_stripes(svg: &mut SvgDrawing, n: usize, width: f32) {
    for i in 0..n {
        if i % 2 != 0 {
            svg.polygon()
                .from_corners((0., K * i as f32), (width, K * i as f32 + FONT_SIZE))
                .style(|s| {
                    s.fill_color(if i % 4 == 1 {
                        Some(StyleColor::String("#fcf7d9".into()))
                    } else {
                        Some(StyleColor::String("#e5f8f9".into()))
                    })
                });
        }
    }
}

fn draw_nodes_in_tree(
    svg: &mut SvgDrawing,
    nodes: &HashMap<String, Vec<f32>>,
    species_map: &HashMap<String, (f32, f32)>,
) {
    for mrca in nodes.keys() {
        let dups = &nodes[mrca];
        let opacity = 1. / dups.len() as f32;
        let (mut x, mut y) = species_map.get(mrca).unwrap().clone();
        for dcs in dups {
            let c = StyleColor::Percent(1. - dcs, *dcs, 0.);
            svg.polygon()
                .from_pos_dims(x - 3., y - 3. + K / 2., 6., 6.)
                .style(|s| s.fill_color(Some(c)).fill_opacity(opacity));
            x += 1.;
            y += 1.;
        }
    }
}

// Returns (SvgGroup, map speciesname -> (coords))
fn draw_species_tree(
    species_tree: &NewickTree,
    species_to_render: &[&str],
    present_species: &[&str],
) -> (Group, HashMap<String, (f32, f32)>) {
    fn render_node(
        svg: &mut Group,
        x: f32,
        xlabels: f32,
        y: f32,
        t: &NewickTree,
        n: usize,
        species_to_render: &[&str],
        present_species: &[&str],
        species_map: &mut HashMap<String, (f32, f32)>,
    ) -> f32 {
        let mut y = y;
        if t[n].is_leaf() {
            t[n].data.name.as_ref().map(|name| {
                svg.line()
                    .from_coords(x, y + K, xlabels, y + K)
                    .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5))
                    .shift(0., -K / 2.);
                svg.text()
                    .pos(xlabels + K, y + FONT_SIZE)
                    .text(name)
                    .style(|s| {
                        s.fill_color(Some(name2color(name))).fill_opacity(
                            if present_species.contains(&name.as_str()) {
                                1.0
                            } else {
                                0.3
                            },
                        )
                    });
                species_map.insert(name.to_owned(), (xlabels, y))
            });
            y += K;
        } else {
            let old_y = y;

            if let Some(name) = &t[n].data.name {
                species_map.insert(name.to_owned(), (x, y));
            }
            for (i, c) in t.children(n).iter().enumerate() {
                if t.leaves_of(*c).iter().any(|l| {
                    t[*l]
                        .data
                        .name
                        .as_ref()
                        .map(|name| species_to_render.contains(&name.as_str()))
                        .unwrap_or(false)
                }) {
                    if i == 0 {
                        svg.line()
                            .from_coords(x, y + K, x + K, y + K)
                            .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5))
                            .shift(0., -K / 2.);
                    } else {
                        svg.line()
                            .from_coords(x, old_y + K, x, y + K)
                            .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5))
                            .shift(0., -K / 2.);
                        svg.line()
                            .from_coords(x, y + K, x + K, y + K)
                            .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5))
                            .shift(0., -K / 2.);
                    }
                    y = render_node(
                        svg,
                        x + K,
                        xlabels,
                        y,
                        t,
                        *c,
                        species_to_render,
                        present_species,
                        species_map,
                    );
                }
            }
        }
        y
    }

    let mut species_map = HashMap::<String, (f32, f32)>::new();
    let mut out = Group::new();
    render_node(
        &mut out,
        0.,
        species_tree.topological_depth().1 as f32 * K,
        0.,
        species_tree,
        species_tree.root(),
        species_to_render,
        present_species,
        &mut species_map,
    );
    (out, species_map)
}

pub fn draw_duplications_blocks(
    t: &NewickTree,
    species_tree: &NewickTree,
    species_map: &mut HashMap<String, (f32, f32)>,
    render: &RenderSettings,
) -> (Group, HashMap<String, Vec<f32>>) {
    fn species_name(t: &NewickTree, n: usize) -> String {
        t[n].data.attrs["S"].to_string()
    }

    let mut out = Group::new();
    let mut xoffset = 0.;
    let mut duplication_sets: Vec<(Vec<(HashSet<String>, i32, i32)>, f32, usize)> = t
        .inners()
        .filter(|&n| t.is_duplication(n))
        .map(|n| {
            // assert!(children.len() == 2, "has {} children", children.len());
            let arms = t[n]
                .children()
                .iter()
                .map(|&c| {
                    (
                        t.leaves_of(c)
                            .iter()
                            .map(|&n| species_name(t, n))
                            .collect::<HashSet<_>>(),
                        0,
                        0,
                    )
                })
                .collect::<Vec<_>>();
            let all_species = t
                .leaves_of(n)
                .iter()
                .map(|o| {
                    species_tree
                        .find_leaf(|n| {
                            n.name
                                .as_ref()
                                .map(|n| *n == species_name(t, *o))
                                .unwrap_or(false)
                        })
                        .expect(&format!(
                            "{} not found in species tree",
                            species_name(t, *o)
                        ))
                })
                .collect::<Vec<_>>();

            // let dcs = n.data.attrs["DCS"].parse::<f32>().unwrap();
            // let elc_all = n.data.attrs["ELC"].parse::<i32>().unwrap();
            // let elc_large = n.data.attrs["ELLC"].parse::<i32>().unwrap();

            let mrca = species_tree.mrca(&all_species).unwrap();
            (
                arms,
                t[n].data.attrs["DCS"].parse::<f32>().unwrap_or_default(),
                mrca,
            )
        })
        .collect::<Vec<_>>();
    duplication_sets.sort_by_cached_key(|a| {
        (
            -(species_tree.leaves_of(a.2).len() as i64),
            -(a.0.iter().map(|x| x.0.len()).sum::<usize>() as i64),
        )
    });

    let mut dup_nodes: HashMap<String, Vec<f32>> = HashMap::new();
    for d in duplication_sets.iter() {
        let dcs = d.1;
        let c = StyleColor::Percent(1. - dcs, dcs, 0.);
        // let elc_all = d.3;
        // let elc_large = d.4;
        let mrca = d.2;
        let mrca_name = species_tree[mrca].data.name.as_ref().unwrap();
        dup_nodes
            .entry(mrca_name.to_owned())
            .or_insert(vec![])
            .push(dcs);

        let y_min =
            d.0.iter()
                .flat_map(|a| a.0.iter().map(|s| species_map.get(s).unwrap().1))
                .fold(f32::INFINITY, f32::min);
        let y_max =
            d.0.iter()
                .flat_map(|a| a.0.iter().map(|s| species_map.get(s).unwrap().1))
                .fold(f32::NEG_INFINITY, f32::max);

        out.polygon()
            .from_corners(
                (xoffset, y_min),
                (xoffset + d.0.len() as f32 * K, y_max + K),
            )
            .style(|s| {
                s.fill_color(Some(c.clone()));
                s.fill_opacity(0.3)
            });

        for (shift, arm) in d.0.iter().enumerate() {
            for species in arm.0.iter() {
                let y = species_map.get(species).unwrap().1;
                out.polygon()
                    .from_pos_dims(xoffset + shift as f32 * K, y, K, K)
                    .style(|s| s.fill_color(Some(c.clone())));
            }
        }

        let label_offset = 0.;
        if render.cs {
            out.text()
                .pos(xoffset + 1.1 * label_offset * K, y_min)
                .text(format!("DCS:{:.1}%", 100. * dcs))
                .transform(|t| t.rotate_from(-45., xoffset, y_min));
        }
        // if render.elc {
        //     out.text()
        //         .pos(xoffset + 1.1 * label_offset * K, y_min)
        //         .text(format!("ELC:{}", elc_all))
        //         .transform(|t| t.rotate_from(-45., xoffset + 1.2 * K, y_min));
        //     label_offset += 1.;
        // }
        // if render.ellc {
        //     out.text()
        //         .pos(xoffset + 1.1 * label_offset * K, y_min)
        //         .text(format!("ELLC:{}", elc_large))
        //         .transform(|t| t.rotate_from(-45., xoffset + 2.4 * K, y_min));
        //     label_offset += 1.;
        // }

        xoffset += d.0.len() as f32 * K + 10.;
    }
    (out, dup_nodes)
}

pub fn render(
    t: &NewickTree,
    species_tree_filename: &str,
    out_filename: &str,
    filter_species_tree: bool,
    render: &RenderSettings,
) {
    let mut svg = SvgDrawing::new();
    let species_tree = newick::one_from_filename(species_tree_filename).unwrap();
    let species_in_tree = t
        .leaves()
        .map(|s| t[s].data.attrs["S"].as_str())
        .collect::<HashSet<&str>>();
    let species_to_render = species_tree
        .leaf_names()
        .filter(|s| !filter_species_tree || species_in_tree.contains(s))
        .collect::<Vec<_>>();
    let present_species = species_tree
        .leaf_names()
        .filter(|s| species_in_tree.contains(s))
        .collect::<Vec<_>>();

    let (tree_group, mut present_species_map) =
        draw_species_tree(&species_tree, &species_to_render, &present_species);
    let (mut dups_group, dups_nodes) =
        draw_duplications_blocks(t, &species_tree, &mut present_species_map, render);
    dups_group.shift(tree_group.bbox().x2, 0.);
    draw_stripes(&mut svg, species_to_render.len(), dups_group.bbox().x2);
    svg.push(Box::new(tree_group));
    svg.push(Box::new(dups_group));
    draw_nodes_in_tree(&mut svg, &dups_nodes, &present_species_map);

    svg.auto_fit();
    let mut out = File::create(out_filename).unwrap();
    out.write_all(svg.render_svg().as_bytes()).unwrap();
}
