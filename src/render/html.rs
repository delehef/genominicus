use crate::align;
use crate::utils::*;
use askama::Template;
use newick::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;

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

    #[serde(rename = "isDuplication")]
    is_duplication: bool,
    confidence: f32,
    repr: Landscape,
    clustered: Option<Vec<PolyGene>>,
}

fn draw_html(tree: &NewickTree, genes: &GeneCache, colormap: &ColorMap) -> HtmlNode {
    fn process(tree: &NewickTree, node: usize, genes: &GeneCache, colormap: &ColorMap) -> HtmlNode {
        let descendants = tree.descendants(node);
        let mut common_ancestral = String::new();
        let clustered = {
            // Node ID to tail mapping
            let tails = descendants
                .iter()
                .filter_map(|&d| {
                    if let Some(gene_name) = tree.name(d) {
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
            if let Some(gene_name) = &tree.name(node) {
                if let Some(DbGene {
                    ancestral,
                    species,
                    chr,
                    t_len,
                    left_tail,
                    right_tail,
                    ..
                }) = genes.get(gene_name.as_str())
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
                            gene_name.to_string(),
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
                me: Gene {
                    color: gene2color(&common_ancestral).to_hex_string(),
                    name: common_ancestral,
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

        min_t_len: &'a i32,
        max_t_len: &'a i32,

        data: &'a str,
    }

    let t_lens = t
        .leaves()
        .filter_map(|n| t.name(n).and_then(|gene_name| genes.get(gene_name)))
        .map(|x| x.t_len)
        .filter(|&x| x >= 0)
        .collect::<Vec<_>>();

    let min_t_len = t_lens.iter().min().unwrap_or(&0);
    let max_t_len = t_lens.iter().max().unwrap_or(&0);

    let html = GenominicusTemplate {
        css: include_str!("../../templates/genominicus.css"),
        js_genominicus: include_str!("../../templates/genominicus.js"),
        js_svg: include_str!("../../templates/svg.min.js"),
        title: out_filename,
        comment: &format!("Transcript length: {} to {}", min_t_len, max_t_len),
        min_t_len,
        max_t_len,
        data: &serde_json::to_string_pretty(&draw_html(t, genes, colormap)).unwrap(),
    };
    let mut out = File::create(out_filename).unwrap();
    let _ = out.write(html.render().unwrap().as_bytes()).unwrap();
}
