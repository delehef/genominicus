#![allow(dead_code)]
use colorsys::{Hsl, Rgb};
use newick::*;
use once_cell::sync::OnceCell;
use palette::*;
use rand::prelude::*;
use rusqlite::*;
use std::collections::{HashMap, HashSet};
use std::iter::FromIterator;
use svarog::*;

static ANCESTRAL_QUERY: OnceCell<String> = OnceCell::new();
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

#[derive(Debug)]
pub struct RenderSettings {
    pub inner_nodes: bool,
    pub cs: bool,
    pub elc: bool,
    pub ellc: bool,
    pub links: bool,
}
impl Default for RenderSettings {
    fn default() -> Self {
        RenderSettings {
            inner_nodes: false,
            cs: false,
            elc: false,
            ellc: false,
            links: false,
        }
    }
}

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

pub fn nuc_to_str(nuc: &str) -> String {
    nuc.to_string()
}

fn jaccard<T: std::hash::Hash + Eq>(x: &HashSet<T>, y: &HashSet<T>) -> f32 {
    x.intersection(y).count() as f32 / x.union(y).count() as f32
}

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

pub fn set_reference(reference: &str) {
    ANCESTRAL_QUERY.set(format!(
        "select ancestral, species, chr, start, direction, left_tail_names, right_tail_names from genomes where {}=?", reference)).unwrap();
}

fn get_gene(db: &mut Connection, name: &str) -> std::result::Result<DbGene, rusqlite::Error> {
    db.query_row(ANCESTRAL_QUERY.get().unwrap(), &[&name], |r| {
        let ancestral: String = r.get("ancestral").unwrap();
        let species: String = r.get("species").unwrap();
        let chr: String = r.get("chr").unwrap();
        let pos: i32 = r.get("start").unwrap();
        let t_len: i32 = r.get("t_len").unwrap_or_default();

        let left_tail: String = r.get("left_tail_names").unwrap();
        let left_tail = left_tail
            .split('.')
            .collect::<Vec<_>>()
            .iter()
            .rev()
            .take(WINDOW as usize)
            .map(|x| x.to_string())
            .filter(|x| !x.is_empty())
            .collect();

        let right_tail: String = r.get("right_tail_names").unwrap();
        let right_tail = right_tail
            .split('.')
            .take(WINDOW as usize)
            .map(|x| x.to_string())
            .filter(|x| !x.is_empty())
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

pub fn make_colormap(tree: &NewickTree, genes: &GeneCache) -> ColorMap {
    let mut colormap = ColorMap::new();
    for l in tree.leaves() {
        if let Some(g) = tree[l].data.name.as_ref().and_then(|name| genes.get(name)) {
            g.left_tail.iter().chain(g.right_tail.iter()).for_each(|g| {
                colormap
                    .entry(g.to_string())
                    .or_insert_with(|| gene2color(&g));
            })
        }
    }
    colormap
}

pub fn make_colormap_per_duplication(
    tree: &NewickTree,
    genes: &GeneCache,
    colorize_all: bool,
) -> ColorMap {
    fn create_gradient(
        t: &NewickTree,
        leave_nodes: &[usize],
        genes: &GeneCache,
        colormap: &mut ColorMap,
    ) {
        if leave_nodes.len() < 2 {
            return;
        }
        let tails = leave_nodes
            .iter()
            .filter_map(|l| t[*l].data.name.as_ref())
            .filter_map(|name| genes.get(name))
            .map(|g| {
                g.left_tail
                    .iter()
                    .chain(g.right_tail.iter())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let tailsets = tails
            .iter()
            .map(|t| HashSet::<_>::from_iter(t.iter().map(md5::compute)))
            .collect::<Vec<_>>();
        let scores = tails
            .iter()
            .enumerate()
            .map(|(i, _)| {
                let mut s = 0;
                for (j, _) in tails.iter().enumerate() {
                    if i != j {
                        s += (&tailsets[i] & &tailsets[j]).len();
                    }
                }
                s
            })
            .collect::<Vec<_>>();
        let ref_tail = &tails[scores
            .iter()
            .enumerate()
            .max_by_key(|(_, s)| *s)
            .map(|(i, _)| i)
            .unwrap_or(0)];

        let mut rng = rand::thread_rng();
        let start = Hsv::new(
            360. * rng.gen::<f64>(),
            0.5 + rng.gen::<f64>() / 2.,
            0.5 + rng.gen::<f64>() / 2.,
        );
        let end = Hsv::new(
            360. * rng.gen::<f64>(),
            0.5 + rng.gen::<f64>() / 2.,
            0.5 + rng.gen::<f64>() / 2.,
        );

        let gradient = Gradient::new(vec![start, end])
            .take(ref_tail.len())
            .collect::<Vec<_>>();
        for (i, gene_name) in ref_tail.iter().enumerate() {
            let color = palette::rgb::Rgb::from_color(gradient[i]);
            colormap
                .entry(gene_name.to_string())
                .or_insert(StyleColor::Percent(
                    color.red as f32,
                    color.green as f32,
                    color.blue as f32,
                ));
        }
    }

    fn rec_fill_colormap(
        tree: &NewickTree,
        node: usize,
        genes: &GeneCache,
        colormap: &mut ColorMap,
    ) {
        if node == 0 || tree.is_duplication(node) {
            let children = tree[node].children();
            let members = children
                .iter()
                .filter(|c| tree[**c].is_leaf())
                .cloned()
                .collect::<Vec<_>>();
            create_gradient(tree, &members, genes, colormap);

            children
                .iter()
                .filter(|c| !tree[**c].is_leaf())
                .for_each(|&c| create_gradient(tree, &tree.leaves_of(c), genes, colormap));
        }

        tree[node]
            .children()
            .iter()
            .for_each(|&c| rec_fill_colormap(tree, c, genes, colormap))
    }

    let mut colormap = ColorMap::new();
    rec_fill_colormap(tree, 0, genes, &mut colormap);
    if colorize_all {
        for l in tree.leaves() {
            if let Some(g) = tree[l].data.name.as_ref().and_then(|name| genes.get(name)) {
                g.left_tail.iter().chain(g.right_tail.iter()).for_each(|g| {
                    colormap
                        .entry(g.to_string())
                        .or_insert_with(|| gene2color(&g));
                })
            }
        }
    }
    colormap
}

pub fn make_genes_cache(t: &NewickTree, db: &mut Connection) -> HashMap<String, DbGene> {
    fn reorder_tails(tree: &NewickTree, node: usize, genes: &mut GeneCache) {
        fn reorder_leaves(t: &NewickTree, leave_nodes: &[usize], genes: &mut GeneCache) {
            if leave_nodes.len() < 2 {
                return;
            }
            // Chose a random leaf from the leaves featuring the longest tails as a reference
            let tails = leave_nodes
                .iter()
                .filter_map(|l| t[*l].data.name.as_ref())
                .filter_map(|name| genes.get(name))
                .map(|g| (g.species.clone(), g.left_tail.clone(), g.right_tail.clone()))
                .collect::<Vec<_>>();
            let tailsets = tails
                .iter()
                .map(|t| HashSet::<_>::from_iter(t.1.iter().chain(t.2.iter()).map(md5::compute)))
                .collect::<Vec<_>>();
            let scores = tails
                .iter()
                .enumerate()
                .map(|(i, _)| {
                    let mut s = 0;
                    for (j, _) in tails.iter().enumerate() {
                        if i != j {
                            s += (&tailsets[i] & &tailsets[j]).len();
                        }
                    }
                    s
                })
                .collect::<Vec<_>>();

            let ref_id = scores
                .iter()
                .enumerate()
                .max_by_key(|(_, s)| *s)
                .map(|(i, _)| i)
                .unwrap_or(0);
            let ref_left_tail: HashSet<&String> = HashSet::from_iter(&tails[ref_id].1);
            let ref_right_tail: HashSet<&String> = HashSet::from_iter(&tails[ref_id].2);

            leave_nodes
                .iter()
                .filter_map(|l| t[*l].data.name.as_ref())
                .for_each(|l_name| {
                    if let Some(gene) = genes.get_mut(l_name) {
                        let left_tail: HashSet<&String> = HashSet::from_iter(&gene.left_tail);
                        let right_tail: HashSet<&String> = HashSet::from_iter(&gene.right_tail);

                        let direct_score = jaccard(&left_tail, &ref_left_tail)
                            + jaccard(&right_tail, &ref_left_tail);
                        let reverse_score = jaccard(&left_tail, &ref_right_tail)
                            + jaccard(&right_tail, &ref_left_tail);

                        if reverse_score > direct_score {
                            std::mem::swap(&mut gene.left_tail, &mut gene.right_tail);
                        }
                    }
                })
        }

        if tree.is_root(node) || tree.is_duplication(node) {
            let children = tree[node].children();
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
        }

        tree[node]
            .children()
            .iter()
            .for_each(|&c| reorder_tails(tree, c, genes))
    }

    let mut r = t
        .leaves()
        .filter_map(|n| t[n].data.name.as_ref())
        .map(|name| {
            (
                name.to_owned(),
                get_gene(db, name).expect(&format!("{} not found", &name)),
            )
        })
        .collect();

    reorder_tails(t, t.root(), &mut r);
    r
}
