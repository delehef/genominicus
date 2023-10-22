#![allow(dead_code)]
use anyhow::*;
use colorsys::{Hsl, Rgb};
use newick::*;
use once_cell::sync::OnceCell;
use palette::*;
use rand::prelude::*;
use std::collections::{HashMap, HashSet};
use std::iter::FromIterator;
use svarog::*;
use syntesuite::genebook::{FamilyID, Gene, GeneBook};

static ANCESTRAL_QUERY: OnceCell<String> = OnceCell::new();
const LEFTS_QUERY: &str = "select ancestral, direction from genomes where species=? and chr=? and start<? order by start desc limit ?";
const RIGHTS_QUERY: &str = "select ancestral, direction from genomes where species=? and chr=? and start>? order by start asc limit ?";

pub const WINDOW: usize = 15;
pub const GENE_WIDTH: f32 = 15.;
pub const GENE_SPACING: f32 = 5.;
pub const BRANCH_WIDTH: f32 = 20.;
pub const FONT_SIZE: f32 = 10.;

#[derive(Copy, Clone, Hash)]
pub enum PoaElt {
    Gene(FamilyID),
    Marker,
    Indel,
    Empty,
}
impl std::fmt::Display for PoaElt {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                PoaElt::Gene(id) => id.to_string(),
                PoaElt::Marker => "===MARKER===".to_string(),
                PoaElt::Indel => "-".to_string(),
                PoaElt::Empty => "".to_string(),
            }
        )
    }
}
impl std::cmp::PartialEq for PoaElt {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (PoaElt::Gene(x), PoaElt::Gene(y)) => x.eq(y),
            (PoaElt::Marker, PoaElt::Marker)
            | (PoaElt::Indel, PoaElt::Indel)
            | (PoaElt::Empty, PoaElt::Empty) => true,
            _ => false,
        }
    }
}
impl Eq for PoaElt {}

#[derive(Debug, Default)]
pub struct RenderSettings {
    pub inner_nodes: bool,
    pub cs: bool,
    pub elc: bool,
    pub ellc: bool,
    pub links: bool,
    pub duplication_ids: bool,
}

pub type GeneCache = HashMap<String, Gene>;
pub type ColorMap = HashMap<usize, StyleColor>;
pub type PetnameMap = HashMap<usize, String>;

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

pub fn gene2color(id: &[u8]) -> StyleColor {
    let bytes: [u8; 16] = md5::compute(id).into();
    let r = (bytes[0] as f32 / 255.).clamp(0.1, 0.9);
    let g = (bytes[1] as f32 / 255.).clamp(0.1, 0.9);
    let b = (bytes[2] as f32 / 255.).clamp(0.1, 0.9);
    StyleColor::Percent(r, g, b)
}

pub fn set_reference(reference: &str) {
    ANCESTRAL_QUERY.set(format!(
        "select ancestral, species, chr, start, direction, left_tail_names, right_tail_names from genomes where {}=?", reference)).unwrap();
}

pub fn make_petnamemap(tree: &NewickTree, genes: &GeneCache) -> PetnameMap {
    let mut petmap = PetnameMap::new();
    for l in tree.leaves() {
        if let Some(g) = tree
            .name(l)
            .as_ref()
            .and_then(|name| genes.get(name.as_str()))
        {
            for f in g
                .left_landscape
                .iter()
                .map(|tg| tg.family)
                .chain(std::iter::once(g.family))
                .chain(g.right_landscape.iter().map(|tg| tg.family))
            {
                petmap
                    .entry(f)
                    .or_insert_with(|| petname::Petnames::default().generate_one(2, "-"));
            }
        }
    }
    petmap
}

pub fn make_colormap(tree: &NewickTree, genes: &GeneCache) -> ColorMap {
    let mut colormap = ColorMap::new();
    for l in tree.leaves() {
        if let Some(g) = tree
            .name(l)
            .as_ref()
            .and_then(|name| genes.get(name.as_str()))
        {
            for tg in g.left_landscape.iter().chain(g.right_landscape.iter()) {
                colormap
                    .entry(tg.family)
                    .or_insert_with(|| gene2color(&tg.family.to_ne_bytes()));
            }
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
            .filter_map(|l| t.name(*l))
            .filter_map(|name| genes.get(name.as_str()))
            .map(|g| {
                g.left_landscape
                    .iter()
                    .chain(g.right_landscape.iter())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let tailsets = tails
            .iter()
            .map(|t| {
                HashSet::<_>::from_iter(t.iter().map(|tg| md5::compute(tg.family.to_ne_bytes())))
            })
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
        for (i, tail_gene) in ref_tail.iter().enumerate() {
            let color = palette::rgb::Rgb::from_color(gradient[i]);
            colormap
                .entry(tail_gene.family)
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

            for c in children.iter().filter(|&c| !tree[*c].is_leaf()) {
                create_gradient(tree, &tree.leaves_of(*c), genes, colormap)
            }
        }

        for c in tree[node].children().iter() {
            rec_fill_colormap(tree, *c, genes, colormap)
        }
    }

    let mut colormap = ColorMap::new();
    rec_fill_colormap(tree, 0, genes, &mut colormap);
    if colorize_all {
        for l in tree.leaves() {
            if let Some(g) = tree.name(l).and_then(|name| genes.get(name.as_str())) {
                for tg in g.left_landscape.iter().chain(g.right_landscape.iter()) {
                    colormap
                        .entry(tg.family)
                        .or_insert_with(|| gene2color(&tg.family.to_ne_bytes()));
                }
            }
        }
    }
    colormap
}

pub fn make_genes_cache(
    t: &NewickTree,
    db_file: &str,
    id_column: &str,
) -> Result<HashMap<String, Gene>> {
    fn reorder_tails(tree: &NewickTree, node: usize, genes: &mut GeneBook) {
        fn reorder_leaves(t: &NewickTree, leave_nodes: &[usize], genes: &mut GeneBook) {
            if leave_nodes.len() < 2 {
                return;
            }
            // Chose a random leaf from the leaves featuring the longest tails as a reference
            let tails = leave_nodes
                .iter()
                .filter_map(|l| t.name(*l))
                .filter_map(|name| genes.get(name).ok())
                .map(|g| {
                    (
                        g.strand,
                        g.left_landscape
                            .iter()
                            .map(|g| g.family)
                            .collect::<Vec<_>>(),
                        g.right_landscape
                            .iter()
                            .map(|g| g.family)
                            .collect::<Vec<_>>(),
                    )
                })
                .collect::<Vec<_>>();
            let tailsets = tails
                .iter()
                .map(|t| {
                    HashSet::<_>::from_iter(
                        t.1.iter()
                            .chain(t.2.iter())
                            .map(|tg| md5::compute(tg.to_ne_bytes())),
                    )
                })
                .collect::<Vec<_>>();
            // scores: Vec<(id: usize, scores: usize)>
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
            let ref_left_tail: HashSet<_> = HashSet::from_iter(tails[ref_id].1.iter().cloned());
            let ref_right_tail: HashSet<_> = HashSet::from_iter(tails[ref_id].2.iter().cloned());

            for l_name in leave_nodes.iter().filter_map(|l| t.name(*l)) {
                if let Result::Ok(gene) = genes.get_mut(l_name) {
                    let left_tail: HashSet<FamilyID> =
                        HashSet::from_iter(gene.left_landscape.iter().map(|tg| tg.family));
                    let right_tail: HashSet<FamilyID> =
                        HashSet::from_iter(gene.right_landscape.iter().map(|tg| tg.family));

                    let direct_score =
                        jaccard(&left_tail, &ref_left_tail) + jaccard(&right_tail, &ref_right_tail);
                    let reverse_score =
                        jaccard(&left_tail, &ref_right_tail) + jaccard(&right_tail, &ref_left_tail);

                    if reverse_score > direct_score {
                        std::mem::swap(&mut gene.left_landscape, &mut gene.right_landscape);
                        gene.left_landscape.reverse();
                        gene.right_landscape.reverse();
                        gene.strand.reverse();
                    }
                }
            }
        }

        if tree.is_root(node) || tree.is_duplication(node) {
            let children = tree[node].children();
            let members = children
                .iter()
                .filter(|c| tree[**c].is_leaf())
                .cloned()
                .collect::<Vec<_>>();
            reorder_leaves(tree, &members, genes);

            for c in children.iter().filter(|c| !tree[**c].is_leaf()) {
                reorder_leaves(tree, &tree.leaves_of(*c), genes);
            }
        }

        for c in tree[node].children().iter() {
            reorder_tails(tree, *c, genes);
        }
    }

    let leaves = t.leaves().filter_map(|n| t.name(n)).collect::<Vec<_>>();
    let mut gene_book =
        GeneBook::cached(db_file, WINDOW, id_column, &leaves).map_err(|e| anyhow!(e))?;
    reorder_tails(t, t.root(), &mut gene_book);
    let r = leaves
        .into_iter()
        .map(|g| {
            gene_book.get(g).map(|mut gene| {
                // XXX: left tails are drawn right to left, so they must be reversed
                gene.left_landscape.reverse();
                (g.to_owned(), gene)
            })
        })
        .collect::<Result<HashMap<_, _>>>()?;
    Ok(r)
}
