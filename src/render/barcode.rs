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
                        StyleColor::String("#fcf7d9".into())
                    } else {
                        StyleColor::String("#e5f8f9".into())
                    })
                });
        }
    }
}

fn draw_species_tree(
    species_tree: &Tree,
    present_species: &[String],
) -> (Group, HashMap<String, f32>) {
    fn render_node(
        svg: &mut Group,
        x: f32,
        xlabels: f32,
        y: f32,
        t: &Tree,
        n: usize,
        present_species: &[String],
        species_map: &mut HashMap<String, f32>,
    ) -> f32 {
        let mut y = y;
        if t[n].is_leaf() {
            t[n].name.as_ref().map(|name| {
                svg.line()
                    .from_coords(x, y + K, xlabels, y + K)
                    .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5))
                    .shift(0., -K / 2.);
                svg.text()
                    .pos(xlabels + K, y + FONT_SIZE)
                    .text(name)
                    .style(|s| s.fill_color(name2color(name)));
                species_map.insert(name.to_owned(), y)
            });
            y += K;
        } else {
            let old_y = y;
            for (i, c) in t.children(n).iter().enumerate() {
                if t.leaves_of(*c).iter().any(|l| {
                    t[*l]
                        .name
                        .as_ref()
                        .map(|name| present_species.contains(name))
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
                    y = render_node(svg, x + K, xlabels, y, t, *c, present_species, species_map);
                }
            }
        }
        y
    }

    let mut species_map = HashMap::<String, f32>::new();
    let mut out = Group::new();
    render_node(
        &mut out,
        0.,
        species_tree.topological_depth().1 * K,
        0.,
        species_tree,
        0,
        present_species,
        &mut species_map,
    );
    (
        out,
        present_species
            .iter()
            .enumerate()
            .map(|(i, s)| (s.to_owned(), i as f32 * K))
            .collect::<HashMap<_, _>>(),
    )
}

pub fn draw_duplications_blocks(t: &Tree, species: &HashMap<String, f32>) -> Group {
    let mut out = Group::new();
    let mut xoffset = 0.;
    let mut duplication_sets = t
        .inners()
        .filter(|n| t[*n].is_duplication())
        .map(|n| {
            fn species_name(t: &Tree, n: usize) -> String {
                t[n].name
                    .as_ref()
                    .unwrap()
                    .split('#')
                    .nth(1)
                    .unwrap()
                    .to_owned()
            }

            let n = &t[n];
            let children = n.children.as_ref().unwrap();
            assert!(children.len() == 2);
            let lefts = t
                .leaves_of(children[0])
                .iter()
                .map(|n| species_name(t, *n))
                .collect::<HashSet<_>>();
            let rights = t
                .leaves_of(children[1])
                .iter()
                .map(|n| species_name(t, *n))
                .collect::<HashSet<_>>();
            let dcs = n.data["DCS"].parse::<f32>().unwrap();
            let elc_all = n.data["ELC"].parse::<i32>().unwrap();
            let elc_large = n.data["ELLC"].parse::<i32>().unwrap();

            (lefts, rights, dcs, elc_all, elc_large)
        })
        .collect::<Vec<_>>();
    duplication_sets.sort_by(|a, b| {
        let a_species_ys = a.0.iter().chain(a.1.iter()).map(|s| species[s] as i32);
        let b_species_ys = b.0.iter().chain(b.1.iter()).map(|s| species[s] as i32);
        let a_span = a_species_ys.clone().max().unwrap() - a_species_ys.clone().min().unwrap();
        let b_span = b_species_ys.clone().max().unwrap() - b_species_ys.clone().min().unwrap();
        match b_span.cmp(&a_span) {
            std::cmp::Ordering::Equal => a_species_ys.min().unwrap().cmp(&b_species_ys.min().unwrap()),
            x => x
        }
    });

    for d in duplication_sets {
        let lefts = &d.0;
        let rights = &d.1;
        let dcs = d.2;
        let elc_all = d.3;
        let elc_large = d.4;
        let c = StyleColor::Percent(1. - dcs, dcs, 0.);

        let left_min = lefts
            .iter()
            .map(|s| *species.get(s).unwrap())
            .fold(f32::INFINITY, f32::min);
        let left_max = lefts
            .iter()
            .map(|s| *species.get(s).unwrap())
            .fold(f32::NEG_INFINITY, f32::max);
        let right_min = rights
            .iter()
            .map(|s| *species.get(s).unwrap())
            .fold(f32::INFINITY, f32::min);
        let right_max = rights
            .iter()
            .map(|s| *species.get(s).unwrap())
            .fold(f32::NEG_INFINITY, f32::max);
        let y_min = left_min.min(right_min);
        let y_max = left_max.max(right_max);
        out.polygon()
            .from_corners((xoffset, y_min), (xoffset + 2. * K, y_max + K))
            .style(|s| {
                s.fill_color(c.clone());
                s.fill_opacity(0.3)
            });

        for s in lefts.iter() {
            let y = *species.get(s).unwrap();
            out.polygon()
                .from_pos_dims(xoffset, y, K, K)
                .style(|s| s.fill_color(c.clone()));
        }

        for s in rights.iter() {
            let y = *species.get(s).unwrap();
            out.polygon()
                .from_pos_dims(xoffset + K, y, K, K)
                .style(|s| s.fill_color(c.clone()));
        }

        out.text()
            .pos(xoffset, y_min)
            .text(format!("DCS:{:.1}%", 100. * dcs))
            .transform(|t| t.rotate_from(-45., xoffset, y_min));
        out.text()
            .pos(xoffset + 1.2 * K, y_min)
            .text(format!("ELC:{}", elc_all))
            .transform(|t| t.rotate_from(-45., xoffset + 1.2 * K, y_min));
        out.text()
            .pos(xoffset + 2.4 * K, y_min)
            .text(format!("ELLC:{}", elc_large))
            .transform(|t| t.rotate_from(-45., xoffset + 2.4 * K, y_min));

        xoffset += 2. * K + 10.;
    }
    out
}

pub fn render(
    t: &Tree,
    species_tree_filename: &str,
    out_filename: &str,
    filter_species_tree: bool,
) {
    let mut svg = SvgDrawing::new();
    let species_tree = Tree::from_filename(species_tree_filename).unwrap();
    let species_in_tree = t
        .leaf_names()
        .iter()
        .map(|(_, s)| s.as_ref().unwrap().split('#').nth(1).unwrap())
        .collect::<HashSet<&str>>();
    let present_species = species_tree
        .leaf_names()
        .iter()
        .map(|(_, s)| s.cloned().unwrap())
        .filter(|s| !filter_species_tree || species_in_tree.contains(s.as_str()))
        .collect::<Vec<_>>();

    let (tree_group, present_species_map) = draw_species_tree(&species_tree, &present_species);
    let mut dups_group = draw_duplications_blocks(t, &present_species_map);
    dups_group.shift(tree_group.bbox().x2, 0.);
    draw_stripes(&mut svg, present_species.len(), dups_group.bbox().x2);
    svg.push(Box::new(tree_group));
    svg.push(Box::new(dups_group));

    svg.auto_fit();
    let mut out = File::create(out_filename).unwrap();
    out.write_all(svg.render_svg().as_bytes()).unwrap();
}
