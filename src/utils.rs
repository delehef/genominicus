#![allow(dead_code)]
use colorsys::{Hsl, Rgb};
use newick::*;
use palette::*;
use rusqlite::*;
use std::collections::{HashMap, HashSet};
use std::iter::FromIterator;
use svarog::*;

const ANCESTRAL_QUERY: &str = concat!(
    "select gene, ancestral, species, chr, start, t_len, direction, left_tail_names, right_tail_names from genomes where gene=?",
);
const LEFTS_QUERY: &str = "select ancestral, direction from genomes where species=? and chr=? and start<? order by start desc limit ?";
const RIGHTS_QUERY: &str = "select ancestral, direction from genomes where species=? and chr=? and start>? order by start asc limit ?";

pub const WINDOW: i64 = 15;
pub const GENE_WIDTH: f32 = 15.;
pub const GENE_SPACING: f32 = 5.;
pub const BRANCH_WIDTH: f32 = 20.;
pub const FONT_SIZE: f32 = 10.;

pub const MARKER: &str = "=======FINAL=======";
pub const INDEL: &str = "-";
pub const EMPTY: &str = "";

pub struct DbGene {
    pub ancestral: String,
    pub species: String,
    pub chr: String,
    pub pos: i32,
    pub t_len: i32,
    pub left_tail: Vec<String>,
    pub right_tail: Vec<String>,
}

pub type GeneCache = HashMap<String, DbGene>;
pub type ColorMap = HashMap<String, StyleColor>;

// Creates a color for a string while trying to ensure it remains readable
pub fn name2color<S: AsRef<str>>(name: S) -> StyleColor {
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

pub fn gene2color<S: AsRef<str>>(name: S) -> StyleColor {
    let bytes: [u8; 16] = md5::compute(name.as_ref().as_bytes()).into();
    let r = (bytes[0] as f32 / 255.).clamp(0.1, 0.9);
    let g = (bytes[1] as f32 / 255.).clamp(0.1, 0.9);
    let b = (bytes[2] as f32 / 255.).clamp(0.1, 0.9);
    StyleColor::Percent(r, g, b)
}

fn get_gene(db: &mut Connection, name: &str) -> std::result::Result<DbGene, rusqlite::Error> {
    db.query_row(ANCESTRAL_QUERY, &[&name], |r| {
        let ancestral: String = r.get("ancestral").unwrap();
        let species: String = r.get("species").unwrap();
        let chr: String = r.get("chr").unwrap();
        let pos: i32 = r.get("start").unwrap();
        let t_len: i32 = r.get("t_len").unwrap();

        let left_tail: String = r.get("left_tail_names").unwrap();
        let left_tail = left_tail
            .split(".")
            .collect::<Vec<_>>()
            .iter()
            .rev()
            .take(WINDOW as usize)
            .map(|x| x.to_string())
            .collect();

        let right_tail: String = r.get("right_tail_names").unwrap();
        let right_tail = right_tail
            .split(".")
            .take(WINDOW as usize)
            .map(|x| x.to_string())
            .collect();

        Ok(DbGene {
            ancestral,
            species,
            chr,
            pos,
            t_len,
            left_tail,
            right_tail,
        })
    })
}

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
        .query_map(params![&species, &chr, pos, span], |row| {
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
        .query_map(params![&species, &chr, pos, span], |row| {
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

pub fn make_colormap(tree: &Tree, genes: &GeneCache) -> ColorMap {
    let mut colormap = ColorMap::new();
    for l in tree.leaves() {
        tree[l]
            .name
            .as_ref()
            .and_then(|name| name.split('#').next())
            .and_then(|name| genes.get(name))
            .map(|g| {
                g.left_tail.iter().chain(g.right_tail.iter()).for_each(|g| {
                    colormap
                        .entry(g.to_string())
                        .or_insert_with(|| gene2color(&g));
                })
            });
    }
    colormap
}

pub fn make_colormap_per_duplication(
    tree: &Tree,
    genes: &GeneCache,
    colorize_all: bool,
) -> ColorMap {
    fn fill_colors(t: &Tree, leave_nodes: &[usize], genes: &GeneCache, colormap: &mut ColorMap) {
        // No need to create a custom gradient for single-gene dups
        if leave_nodes.len() < 2 {
            return;
        }
        // Chose a random leaf from the leaves featuring the longest tails as a reference
        let mut leaves = leave_nodes
            .iter()
            .filter_map(|l| t[*l].name.as_ref().and_then(|name| name.split('#').next()))
            .filter_map(|name| genes.get(name))
            .map(|g| (g.species.clone(), g.left_tail.clone(), g.right_tail.clone()))
            .collect::<Vec<_>>();
        leaves.sort_by(|a, b| {
            let len_a = a.1.len() + a.2.len();
            let len_b = b.1.len() + b.2.len();
            // Select by species name to discriminate between equal syntenic landscapes lenghts
            if len_a == len_b {
                b.0.cmp(&a.0)
            } else {
                len_b.cmp(&len_a)
            }
        });
        let ref_left_tail = &leaves[0].1;
        let ref_right_tail = &leaves[0].2;

        // Use its syntenic landscape to fill the color map
        // and the reference left and right tails
        let start_color = name2color(&ref_left_tail[0]).to_percent();
        let end_color = name2color(&ref_right_tail[ref_right_tail.len() - 1]).to_percent();
        let gradient = Gradient::new(vec![
            LinSrgb::new(start_color.0, start_color.1, start_color.2),
            LinSrgb::new(end_color.0, end_color.1, end_color.2),
        ])
        .take(ref_left_tail.len() + ref_right_tail.len())
        .collect::<Vec<_>>();
        for (i, gene_name) in ref_left_tail
            .iter()
            .chain(ref_right_tail.iter())
            .enumerate()
        {
            let color = &gradient[i];
            colormap
                .entry(gene_name.to_string())
                .or_insert(StyleColor::Percent(color.red, color.green, color.blue));
        }
    }

    fn rec_fill_colormap(tree: &Tree, node: usize, genes: &GeneCache, colormap: &mut ColorMap) {
        if node == 0 || tree[node].is_duplication() {
            tree[node].children.as_ref().map(|children| {
                let members = children
                    .iter()
                    .filter(|c| tree[**c].is_leaf())
                    .cloned()
                    .collect::<Vec<_>>();
                fill_colors(tree, &members, genes, colormap);

                children
                    .iter()
                    .filter(|c| !tree[**c].is_leaf())
                    .for_each(|&c| fill_colors(tree, &tree.leaves_of(c), genes, colormap));
            });
        }

        tree[node].children.as_ref().map(|children| {
            children
                .iter()
                .for_each(|&c| rec_fill_colormap(tree, c, genes, colormap))
        });
    }

    let mut colormap = ColorMap::new();
    rec_fill_colormap(tree, 0, genes, &mut colormap);
    if colorize_all {
        for l in tree.leaves() {
            tree[l]
                .name
                .as_ref()
                .and_then(|name| name.split('#').next())
                .and_then(|name| genes.get(name))
                .map(|g| {
                    g.left_tail.iter().chain(g.right_tail.iter()).for_each(|g| {
                        colormap
                            .entry(g.to_string())
                            .or_insert_with(|| gene2color(&g));
                    })
                });
        }
    }
    colormap
}

pub fn make_genes_cache(t: &Tree, db: &mut Connection) -> HashMap<String, DbGene> {
    fn reorder_tails(tree: &Tree, node: usize, genes: &mut GeneCache) {
        fn reorder_leaves(t: &Tree, leave_nodes: &[usize], genes: &mut GeneCache) {
            fn jaccard<T: std::hash::Hash + Eq>(x: &HashSet<T>, y: &HashSet<T>) -> f32 {
                x.intersection(&y).collect::<Vec<_>>().len() as f32
                    / x.union(&y).collect::<Vec<_>>().len() as f32
            }

            if leave_nodes.len() < 2 {
                return;
            }
            // Chose a random leaf from the leaves featuring the longest tails as a reference
            let mut leaves = leave_nodes
                .iter()
                .filter_map(|l| t[*l].name.as_ref().and_then(|name| name.split('#').next()))
                .filter_map(|name| genes.get(name))
                .map(|g| (g.species.clone(), g.left_tail.clone(), g.right_tail.clone()))
                .collect::<Vec<_>>();
            leaves.sort_by(|a, b| {
                let len_a = a.1.len() + a.2.len();
                let len_b = b.1.len() + b.2.len();
                // Select by species name to discriminate between equal syntenic landscapes lenghts
                if len_a == len_b {
                    b.0.cmp(&a.0)
                } else {
                    len_b.cmp(&len_a)
                }
            });

            let ref_left_tail: HashSet<&String> = HashSet::from_iter(&leaves[0].1);
            let ref_right_tail: HashSet<&String> = HashSet::from_iter(&leaves[0].2);

            leave_nodes
                .iter()
                .filter_map(|l| t[*l].name.as_ref().and_then(|name| name.split('#').next()))
                .for_each(|l_name| {
                    genes.get_mut(l_name).map(|gene| {
                        let left_tail: HashSet<&String> = HashSet::from_iter(&gene.left_tail);
                        let right_tail: HashSet<&String> = HashSet::from_iter(&gene.right_tail);

                        let direct_score = jaccard(&left_tail, &ref_left_tail)
                            + jaccard(&right_tail, &ref_left_tail);
                        let reverse_score = jaccard(&left_tail, &ref_right_tail)
                            + jaccard(&right_tail, &ref_left_tail);

                        if reverse_score > direct_score {
                            std::mem::swap(&mut gene.left_tail, &mut gene.right_tail);
                        }
                    });
                })
        }

        if node == 0 || tree[node].is_duplication() {
            tree[node].children.as_ref().map(|children| {
                let members = children
                    .iter()
                    .filter(|c| tree[**c].is_leaf())
                    .cloned()
                    .collect::<Vec<_>>();
                reorder_leaves(tree, &members, genes);

                children
                    .iter()
                    .filter(|c| !tree[**c].is_leaf())
                    .for_each(|&c| reorder_leaves(tree, &tree.leaves_of(c), genes));
            });
        }

        tree[node]
            .children
            .as_ref()
            .map(|children| children.iter().for_each(|&c| reorder_tails(tree, c, genes)));
    }

    let mut r = t
        .leaves()
        .filter_map(|n| {
            t[n].name
                .as_ref()
                .map(|name| name.split('#').next().unwrap())
        })
        .map(|name| {
            (
                name.to_owned(),
                get_gene(db, &name).expect(&format!("{} not found", &name)),
            )
        })
        .collect();

    reorder_tails(t, 0, &mut r);
    r
}
