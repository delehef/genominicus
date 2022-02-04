use crate::utils::*;
use newick::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::prelude::*;
use svarog::*;

pub fn draw_species_blocks(
    svg: &mut SvgDrawing,
    t: &Tree,
    species_tree: &Tree,
    filter_species_tree: bool,
) {
    const K: f32 = FONT_SIZE;
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
    duplication_sets.sort_by(|a, b| (b.0.len() + b.1.len()).cmp(&(a.0.len() + a.1.len())));
    let species_in_tree = t
        .leaf_names()
        .iter()
        .map(|(_, s)| s.as_ref().unwrap().split('#').nth(1).unwrap())
        .collect::<HashSet<&str>>();

    let species = species_tree
        .leaf_names()
        .iter()
        .map(|(_, s)| s.as_ref().unwrap())
        .filter(|s| !filter_species_tree || species_in_tree.contains(s.as_str()))
        .enumerate()
        .map(|(i, s)| (s.to_owned(), i as f32 * K))
        .collect::<HashMap<_, _>>();
    let longest_species =
        species.iter().map(|(name, _)| name.len()).max().unwrap() as f32 * FONT_SIZE;

    let mut xoffset = longest_species + 50.;
    let x_max = xoffset + duplication_sets.len() as f32 * 3. * K;
    for i in (0..species.len()).step_by(2) {
        svg.polygon()
            .from_corners((0., K*i as f32), (x_max, K*i as f32 + FONT_SIZE))
            .style(|s| {
                s.fill_color(StyleColor::String("#fcf7d9".into()))
            });
    }
    for (species, y) in species.iter() {
        svg.text()
            .pos(0., *y + FONT_SIZE)
            .text(species)
            .style(|s| s.fill_color(name2color(species)));
    }

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
        svg.polygon()
            .from_corners((xoffset, y_min), (xoffset + 2. * K, y_max + K))
            .style(|s| {
                s.fill_color(c.clone());
                s.fill_opacity(0.3)
            });

        for s in lefts.iter() {
            let y = *species.get(s).unwrap();
            svg.polygon()
                .from_pos_dims(xoffset, y, K, K)
                .style(|s| s.fill_color(c.clone()));
        }

        for s in rights.iter() {
            let y = *species.get(s).unwrap();
            svg.polygon()
                .from_pos_dims(xoffset + K, y, K, K)
                .style(|s| s.fill_color(c.clone()));
        }

        svg.text()
            .pos(xoffset + 2. * K, y_min)
            .text(format!("DCS:{:.1}%", 100. * dcs));
        svg.text()
            .pos(xoffset + 2. * K, y_min + K)
            .text(format!("ELC:{}", elc_all));
        svg.text()
            .pos(xoffset + 2. * K, y_min + 2. * K)
            .text(format!("ELLC:{}", elc_large));

        xoffset += 2. * K + 10.;
    }
}

pub fn render(
    t: &Tree,
    species_tree_filename: &str,
    out_filename: &str,
    filter_species_tree: bool,
) {
    let species_tree = Tree::from_filename(species_tree_filename).unwrap();
    let mut svg = SvgDrawing::new();
    draw_species_blocks(&mut svg, t, &species_tree, filter_species_tree);
    svg.auto_fit();
    let mut out = File::create(out_filename).unwrap();
    out.write_all(svg.render_svg().as_bytes()).unwrap();
}
