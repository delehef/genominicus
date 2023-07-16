use crate::align;
use crate::utils::*;
use askama::Template;
use newick::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use syntesuite::genebook::Gene;

#[derive(Serialize, Deserialize)]
struct PolyGene {
    genes: Vec<(HtmlGene, f32)>,
}
#[derive(Serialize, Deserialize)]
struct HtmlGene {
    color: String,
    name: String,
}
#[derive(Serialize, Deserialize)]
struct Landscape {
    lefts: Vec<HtmlGene>,
    me: HtmlGene,
    rights: Vec<HtmlGene>,
}
#[derive(Serialize, Deserialize)]
struct HtmlNode {
    species: String,
    chr: String,
    gene: String,
    ancestral: String,
    color: String,
    children: Vec<HtmlNode>,

    #[serde(rename = "isDuplication")]
    is_duplication: bool,
    confidence: f32,
    repr: Landscape,
    clustered: Option<Vec<PolyGene>>,
}

fn draw_html(tree: &NewickTree, genes: &GeneCache, colormap: &ColorMap) -> HtmlNode {
    fn process(tree: &NewickTree, node: usize, genes: &GeneCache, colormap: &ColorMap) -> HtmlNode {
        let descendants = tree.descendants(node);
        let mut common_ancestral = 0;
        let clustered = {
            // Node ID to tail mapping
            let tails = descendants
                .iter()
                .filter_map(|&d| {
                    if let Some(gene_name) = tree.name(d) {
                        if let Some(Gene {
                            family,
                            left_landscape,
                            right_landscape,
                            ..
                        }) = genes.get(gene_name)
                        {
                            common_ancestral = *family;
                            Some((
                                d,
                                left_landscape
                                    .iter()
                                    .map(|tg| PoaElt::Gene(tg.family))
                                    .rev() // XXX Pour que les POA partent bien du bout
                                    .chain(std::iter::once(PoaElt::Marker))
                                    .chain(right_landscape.iter().map(|tg| PoaElt::Gene(tg.family)))
                                    .collect::<Vec<_>>(),
                            ))
                        } else {
                            None
                        }
                    } else {
                        None
                    }
                })
                .collect::<HashMap<_, _>>();

            if !tails.is_empty() {
                let (g, heads) = align::align(&tails);
                let mut alignment = align::poa_to_strings(&g, &heads)
                    .values()
                    .cloned()
                    .collect::<Vec<_>>();
                if false {
                    // Differentiate between tail of shorter alignments and actual indels
                    alignment.iter_mut().for_each(|a| {
                        for pos in a.iter_mut() {
                            if *pos == PoaElt::Indel {
                                *pos = PoaElt::Empty;
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
                                alignment.iter().filter(|a| a[i] != PoaElt::Empty).count() as f32
                            } else {
                                alignment.len() as f32
                            };

                            let mut counts: HashMap<PoaElt, i32> = HashMap::new();
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
                                        let name = if name != PoaElt::Marker {
                                            name
                                        } else {
                                            PoaElt::Gene(common_ancestral)
                                        };
                                        Some((
                                            HtmlGene {
                                                name: name.to_string(),
                                                color: match name {
                                                    PoaElt::Gene(family) => colormap
                                                        .get(&family)
                                                        .map(|c| c.to_hex_string())
                                                        .unwrap_or_else(|| "#aaa".to_string()),
                                                    PoaElt::Marker => todo!(),
                                                    PoaElt::Indel | PoaElt::Empty => {
                                                        "#ccc".to_string()
                                                    }
                                                },
                                            },
                                            v as f32 / count,
                                        ))
                                        // }
                                    })
                                    .collect::<Vec<(HtmlGene, f32)>>(),
                            }
                        })
                        .collect::<Vec<PolyGene>>(),
                )
            } else {
                None
            }
        };

        let ((species, chr, gene, ancestral), (lefts, rights)) =
            if let Some(gene_name) = &tree.name(node) {
                if let Some(Gene {
                    family,
                    species,
                    chr,
                    left_landscape,
                    right_landscape,
                    ..
                }) = genes.get(gene_name.as_str())
                {
                    common_ancestral = *family;
                    let (proto_lefts, proto_rights) = (left_landscape, right_landscape);
                    let (lefts, rights) = if true {
                        (proto_rights, proto_lefts)
                    } else {
                        (proto_lefts, proto_rights)
                    };
                    (
                        (
                            species.to_owned(),
                            chr.to_owned(),
                            gene_name.to_string(),
                            *family,
                        ),
                        (
                            lefts
                                .iter()
                                .rev()
                                .map(|g| HtmlGene {
                                    name: g.family.to_string(),
                                    color: colormap
                                        .get(&g.family)
                                        .map(|c| c.to_hex_string())
                                        .unwrap_or_else(|| "#aaa".to_string()),
                                })
                                .collect::<Vec<_>>(),
                            rights
                                .iter()
                                .map(|g| HtmlGene {
                                    name: g.family.to_string(),
                                    color: colormap
                                        .get(&g.family)
                                        .map(|c| c.to_hex_string())
                                        .unwrap_or_else(|| "#aaa".to_string()),
                                })
                                .collect::<Vec<_>>(),
                        ),
                    )
                } else {
                    (
                        (String::new(), String::new(), String::new(), 0),
                        (vec![], vec![]),
                    )
                }
            } else {
                (
                    (String::new(), String::new(), String::new(), 0),
                    (vec![], vec![]),
                )
            };

        let color = name2color(&species).to_hex_string();
        HtmlNode {
            species,
            chr,
            gene,
            ancestral: ancestral.to_string(), // FIXME: random name
            color,
            children: tree[node]
                .children()
                .as_ref()
                .iter()
                .map(|n| process(tree, *n, genes, colormap))
                .collect(),
            is_duplication: tree.is_duplication(node),
            confidence: tree
                .attrs(node)
                .get("DCS")
                .and_then(|x| x.parse::<f32>().ok())
                .unwrap_or(0.0),
            repr: Landscape {
                lefts,
                rights,
                me: HtmlGene {
                    color: gene2color(&common_ancestral.to_ne_bytes()).to_hex_string(),
                    name: common_ancestral.to_string(),
                },
            },
            clustered,
        }
    }

    process(tree, tree.root(), genes, colormap)
}

pub fn render(t: &NewickTree, genes: &GeneCache, colormap: &ColorMap, out_filename: &str) {
    #[derive(Template)]
    #[template(path = "genominicus.html", escape = "none")]
    struct GenominicusTemplate<'a> {
        css: &'a str,
        js_svg: &'a str,
        js_genominicus: &'a str,

        title: &'a str,
        comment: &'a str,

        data: &'a str,
    }

    let html = GenominicusTemplate {
        css: include_str!("../../templates/genominicus.css"),
        js_genominicus: include_str!("../../templates/genominicus.js"),
        js_svg: include_str!("../../templates/svg.min.js"),
        title: out_filename,
        comment: "".into(),
        data: &serde_json::to_string_pretty(&draw_html(t, genes, colormap)).unwrap(),
    };
    let mut out = File::create(out_filename).unwrap();
    let _ = out.write(html.render().unwrap().as_bytes()).unwrap();
}
