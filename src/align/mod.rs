#![allow(dead_code)]
#![allow(non_snake_case)] // I like my matrices to be in capitals
use maplit::hashmap;
use std::cmp::max;
use std::collections::HashMap;

use petgraph::prelude::*;

use super::*;
use crate::align::poa::*;

mod poa;

const SIMDW: usize = 8;
const NEG_INF: i32 = -3_000_000;

pub type SeqID = usize;
pub type Nucleotide = PoaElt;
pub type Sequence = Vec<Nucleotide>;
pub type Sequences = HashMap<SeqID, Sequence>;

type Alignment = (Vec<Option<NodeIndex>>, Vec<Option<SeqID>>);
struct AffineNWSettings {
    matches: i32,
    mismatches: i32,
    open_gap: i32,
    extend_gap: i32,
}

fn insert_hanging_seq(
    g: &mut POAGraph,
    seq: &[Nucleotide],
    seq_id: SeqID,
) -> Option<(NodeIndex, NodeIndex)> {
    let mut first_node: Option<NodeIndex> = None;
    let mut last_node: Option<NodeIndex> = None;
    if seq.is_empty() {
        return None;
    }

    for n in seq.iter() {
        let node = g.add_node(POANode {
            nucs: hashmap! {seq_id => *n},
        });
        if first_node.is_none() {
            first_node = Some(node)
        }
        if last_node.is_some() {
            poa::update_edge(g, last_node, Some(node), seq_id);
        }
        last_node = Some(node);
    }

    Some((first_node.unwrap(), last_node.unwrap()))
}

fn add_alignment(
    g: &mut POAGraph,
    alignment: &Alignment,
    seq: &Sequence,
    seq_id: SeqID,
) -> Option<NodeIndex> {
    let graph_ids = &alignment.0;
    let seq_idxs = &alignment.1;

    let valid_seq_idxs = seq_idxs
        .iter()
        .filter_map(|x| x.as_ref())
        .collect::<Vec<_>>();
    if alignment.0.is_empty() {
        return None;
    }
    let start_seq_idx = valid_seq_idxs[0];
    let end_seq_idx = valid_seq_idxs.last().unwrap();

    let (mut first_id, mut head_id) = if *start_seq_idx > 0 {
        insert_hanging_seq(g, &seq[0..*start_seq_idx], seq_id)
            .map_or((None, None), |(first_id, head_id)| {
                (Some(first_id), Some(head_id))
            })
    } else {
        (None, None)
    };

    let tail_id = if **end_seq_idx < seq.len() {
        insert_hanging_seq(g, &seq[*end_seq_idx + 1..], seq_id).map(|(tail_id, _)| tail_id)
    } else {
        None
    };

    for (seq_idx, graph_id) in seq_idxs.iter().zip(graph_ids.iter()) {
        // TODO with a match
        if seq_idx.is_none() {
            continue;
        }

        let n = &seq[seq_idx.unwrap()];

        let new_node = if graph_id.is_none() {
            // Not aligned, create a new node
            g.add_node(POANode {
                nucs: hashmap! {seq_id => *n},
            })
        } else {
            // Aligned, add our nuc if necessary, and graft us on the edge
            g[graph_id.unwrap()].nucs.insert(seq_id, *n);
            // Add edge
            graph_id.unwrap()
        };

        update_edge(g, head_id, Some(new_node), seq_id);
        head_id = Some(new_node);
        if first_id.is_none() {
            first_id = head_id;
        }
    }

    update_edge(g, head_id, tail_id, seq_id);
    Some(first_id.unwrap())
}

fn build_matrix(
    g: &POAGraph,
    seq: &Sequence,
    settings: &AffineNWSettings,
    ranks_to_nodes: &[NodeIndex],
    nodes_to_ranks: &[usize],
    nucs: &[Vec<Nucleotide>],
    H: &mut [i32],
    F: &mut [i32],
    E: &mut [i32],
) {
    let m = settings.matches;
    let n = settings.mismatches;
    let _g = settings.open_gap;
    let e = settings.extend_gap;
    let m_width = seq.len() + 1;
    let m_height = ranks_to_nodes.len() + 1;

    for j in 1..m_width {
        F[j] = NEG_INF;
        E[j] = _g + (j as i32 - 1) * e;
    }
    for i in 1..m_height {
        let edges: Vec<_> = g
            .edges_directed(NodeIndex::new(i - 1), Direction::Incoming)
            .collect();
        let mut penalty = if edges.is_empty() { _g - e } else { NEG_INF };
        for edge in edges {
            let pred_i = nodes_to_ranks[edge.source().index()] + 1;
            penalty = max(penalty, F[pred_i * m_width]);
        }
        F[i * m_width] = penalty + e;
        E[i * m_width] = NEG_INF;
    }

    H[1..m_width].clone_from_slice(&E[1..m_width]);
    for i in 1..m_height {
        H[i * m_width] = F[i * m_width];
    }

    for node_id in ranks_to_nodes.iter() {
        // Process the guaranteed first predecessor (outside of edge conditions)
        let i = nodes_to_ranks[node_id.index()] + 1;
        let row = i * m_width;
        let preds = g
            .edges_directed(*node_id, Direction::Incoming)
            .collect::<Vec<_>>();
        let pred_i = if preds.is_empty() {
            0
        } else {
            nodes_to_ranks[preds[0].source().index()] + 1
        };
        let pred_row = pred_i * m_width;

        for j in 1..m_width {
            H[row + j] = H[pred_row + j - 1]
                + if nucs[node_id.index()].contains(&seq[j - 1]) {
                    if seq[j - 1] == PoaElt::Marker {
                        100
                    } else {
                        m
                    }
                } else {
                    n
                };
            F[row + j] = max(H[pred_row + j] + _g, F[pred_row + j] + e);
        }

        // Then the other putative predecessors
        for p in preds.iter().skip(1) {
            let pred_i = nodes_to_ranks[p.source().index()] + 1;
            let pred_row = pred_i * m_width;

            for j in 1..m_width {
                H[row + j] = max(
                    H[row + j],
                    H[pred_row + j - 1]
                        + if nucs[node_id.index()].contains(&seq[j - 1]) {
                            if seq[j - 1] == PoaElt::Marker {
                                100
                            } else {
                                m
                            }
                        } else {
                            n
                        },
                );
                F[row + j] = max(F[row + j], max(H[pred_row + j] + _g, F[pred_row + j] + e));
            }
        }

        for j in 1..m_width {
            E[row + j] = max(H[row + j - 1] + _g, E[row + j - 1] + e);
            H[row + j] = max(H[row + j], max(F[row + j], E[row + j]));
        }
    }
}

fn print_matrix(
    M: &[i32],
    m_width: usize,
    m_height: usize,
    seq: &Sequence,
    ranks_to_nodes: &[NodeIndex],
    g: &POAGraph,
) {
    print!("                    ");
    for j in 1..m_width {
        print!("{:>6}", j);
    }
    println!();

    print!("                    ");
    for j in 1..m_width {
        print!("{:>20}", seq[j - 1]);
    }
    println!();

    for i in 0..m_height {
        if i > 0 {
            print!(
                "{:6}  {:20}",
                i,
                g[ranks_to_nodes[i - 1]]
                    .nucs
                    .values()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(" ")
            );
        } else {
            print!("                          ")
        }
        for j in 0..m_width {
            print!("{:6}", M[i * m_width + j]);
        }
        println!();
    }
    println!();
}

fn affine_sw(g: &POAGraph, seq: &Sequence, settings: &AffineNWSettings) -> (i32, Alignment) {
    let m = settings.matches;
    let n = settings.mismatches;
    let _g = settings.open_gap;
    let e = settings.extend_gap;
    // Affine gap penalty alignment
    let ranks_to_nodes = petgraph::algo::toposort(g, None).ok().unwrap();
    let mut nodes_to_ranks = vec![0; ranks_to_nodes.len()];
    for i in 0..ranks_to_nodes.len() {
        nodes_to_ranks[ranks_to_nodes[i].index()] = i;
    }

    let m_width = seq.len() + 1;
    let m_height = ranks_to_nodes.len() + 1;
    #[allow(non_snake_case)]
    let mut H: Vec<i32> = vec![0; m_width * m_height];
    #[allow(non_snake_case)]
    let mut F: Vec<i32> = vec![0; m_width * m_height];
    #[allow(non_snake_case)]
    let mut E: Vec<i32> = vec![0; m_width * m_height];

    let mut nucs = vec![vec![]; ranks_to_nodes.len()];
    for node in ranks_to_nodes.iter() {
        // dbg!(g[*node].nucs.values());
        nucs[node.index()] = g[*node].nucs.values().cloned().collect::<Vec<_>>();
    }

    build_matrix(
        g,
        seq,
        settings,
        &ranks_to_nodes,
        &nodes_to_ranks,
        &nucs,
        &mut H,
        &mut F,
        &mut E,
    );

    // Backtrack
    let max_j = m_width - 1;
    let (max_score, max_i) = ranks_to_nodes.iter().fold((NEG_INF, 0), |ax, &node_id| {
        let (max_score, _) = ax;
        let i = nodes_to_ranks[node_id.index()] + 1;
        let score = H[i * m_width + m_width - 1];
        let is_a_tail = g
            .edges_directed(node_id, Direction::Outgoing)
            .next()
            .is_none();

        if is_a_tail && score > max_score {
            (score, i)
        } else {
            ax
        }
    });

    let mut seq_idxs = Vec::with_capacity(seq.len());
    let mut graph_idxs = Vec::with_capacity(nodes_to_ranks.len());

    let mut i = max_i;
    let mut j = max_j;
    let mut prev_j = 0;
    let mut prev_i = 0;
    while i != 0 && j != 0 {
        // print_matrix(&H, m_width, m_height, seq, &ranks_to_nodes, g);
        // dbg!(i, j);
        let hij = H[i * m_width + j];
        let mut pred_found = false;
        let mut extend_left = false;
        let mut extend_up = false;

        // Looking for a direct match...
        if i != 0 && j != 0 {
            let node_id = ranks_to_nodes[i - 1];
            let preds = g
                .edges_directed(node_id, Direction::Incoming)
                .collect::<Vec<_>>();

            // ...first in the directly preceding node...
            let match_cost = if nucs[node_id.index()].contains(&seq[j - 1]) {
                if seq[j - 1] == PoaElt::Marker {
                    100
                } else {
                    m
                }
            } else {
                n
            };
            let pred_i = if preds.is_empty() {
                0
            } else {
                nodes_to_ranks[preds[0].source().index()] + 1
            };
            if hij == H[pred_i * m_width + j - 1] + match_cost {
                prev_i = pred_i;
                prev_j = j - 1;
                pred_found = true;
            } else {
                // ...or in the other preceding nodes.
                for p in preds.iter().skip(1) {
                    let pred_i = nodes_to_ranks[p.source().index()] + 1;

                    if hij == H[pred_i * m_width + j - 1] + match_cost {
                        prev_i = pred_i;
                        prev_j = j - 1;
                        pred_found = true;
                        break;
                    }
                }
            }
        }

        // No match found. Try to add an insertion on the sequence...
        if !pred_found && i != 0 {
            let node_id = ranks_to_nodes[i - 1];
            let preds = g
                .edges_directed(node_id, Direction::Incoming)
                .collect::<Vec<_>>();
            let pred_i = if preds.is_empty() {
                0
            } else {
                nodes_to_ranks[preds[0].source().index()] + 1
            };
            // ...first from the directly preceding node...
            if hij == F[pred_i * m_width + j] + e {
                extend_up = true;
                prev_i = pred_i;
                prev_j = j;
                pred_found = true;
            } else if hij == H[pred_i * m_width + j] + _g {
                prev_i = pred_i;
                prev_j = j;
                pred_found = true;
            } else {
                // ...or from the other preceding ones.
                for p in preds.iter().skip(1) {
                    let pred_i = nodes_to_ranks[p.source().index()] + 1;
                    if hij == F[pred_i * m_width + j] + e {
                        extend_up = true;
                        prev_i = pred_i;
                        prev_j = j;
                        pred_found = true;
                        break;
                    } else if hij == H[pred_i * m_width + j] + _g {
                        prev_i = pred_i;
                        prev_j = j;
                        pred_found = true;
                        break;
                    }
                }
            }
        }

        // No insertion either; add a deletion on the sequence.
        if !pred_found && j != 0 {
            if hij == E[i * m_width + j - 1] + e {
                extend_left = true;
                prev_i = i;
                prev_j = j - 1;
                // pred_found = true;
            } else if hij == H[i * m_width + j - 1] + _g {
                prev_i = i;
                prev_j = j - 1;
                // pred_found = true;
            }
        }

        graph_idxs.push(if i == prev_i {
            None
        } else {
            Some(ranks_to_nodes[i - 1])
        });
        seq_idxs.push(if j == prev_j { None } else { Some(j - 1) });
        i = prev_i;
        j = prev_j;

        // Extend the gap if we opened one
        if extend_left {
            loop {
                graph_idxs.push(None);
                seq_idxs.push(Some(j - 1));
                j -= 1;
                if E[i * m_width + j] + e != E[i * m_width + j + 1] {
                    break;
                }
            }
        } else if extend_up {
            loop {
                let mut stop = false;
                prev_i = 0;
                for p in g.edges_directed(ranks_to_nodes[i - 1], Direction::Incoming) {
                    let pred_i = nodes_to_ranks[p.source().index()] + 1;
                    if F[i * m_width + j] == H[pred_i * m_width + j] + _g {
                        stop = true;
                        prev_i = pred_i;
                        break;
                    } else if F[i * m_width + j] == F[pred_i * m_width + j] + e {
                        prev_i = pred_i;
                        break;
                    }
                }

                graph_idxs.push(Some(ranks_to_nodes[i - 1]));
                seq_idxs.push(None);
                i = prev_i;
                if stop || i == 0 {
                    break;
                }
            }
        }
    }

    graph_idxs.reverse();
    seq_idxs.reverse();

    (max_score, (graph_idxs, seq_idxs))
}

pub fn align(seqs: &Sequences) -> (POAGraph, HashMap<SeqID, NodeIndex>) {
    let settings = AffineNWSettings {
        matches: 10,
        mismatches: -0,
        open_gap: -1,
        extend_gap: -1,
    };
    let mut g = POAGraph::new();
    // We sort the input sequences for two reasons:
    // 1. align the longer ones first, so that the resulting NW MSA is more resilient to
    //    small sequences insertion;
    // 2. Saving memory is always welcome.
    let mut sorted_seqs = seqs.iter().collect::<Vec<_>>();
    sorted_seqs.sort_by(|(_, seq1), (_, seq2)| seq2.len().cmp(&seq1.len()));

    let starts = sorted_seqs
        .iter()
        .enumerate()
        .filter_map(|(i, (&id, seq))| {
            if i == 0 {
                insert_hanging_seq(&mut g, seq, id).map(|new| (id, new.0))
            } else {
                let rev_seq = seq.iter().cloned().rev().collect();
                let (direct_score, direct_alignment) = affine_sw(&g, seq, &settings);
                let (reverse_score, reverse_alignment) = affine_sw(&g, &rev_seq, &settings);
                if direct_score >= reverse_score {
                    add_alignment(&mut g, &direct_alignment, seq, id).map(|new| (id, new))
                } else {
                    add_alignment(&mut g, &reverse_alignment, &rev_seq, id).map(|new| (id, new))
                }
            }
        })
        .collect::<HashMap<_, _>>();

    (g, starts)
}

pub fn poa_to_strings(g: &POAGraph, starts: &Heads) -> HashMap<usize, Vec<PoaElt>> {
    let nodes = petgraph::algo::toposort(g, None).ok().unwrap();
    let rank_to_column = nodes
        .iter()
        .enumerate()
        .map(|(i, n)| (n, i))
        .collect::<HashMap<_, _>>();

    starts
        .iter()
        .map(|(seq_id, start)| {
            let mut seq_out = vec![PoaElt::Indel; nodes.len()];
            let mut node = *start;

            loop {
                let nucs = &g[node].nucs;
                seq_out[rank_to_column[&node]] = nucs[seq_id];

                if let Some(next) = g
                    .edges_directed(node, Direction::Outgoing)
                    .find(|e| e.weight().contains(seq_id))
                    .map(|e| e.target())
                // Check there is only one?
                {
                    // i += 1;
                    node = next;
                } else {
                    break;
                }
            }
            (seq_id, seq_out)
        })
        .fold(HashMap::new(), |mut ax, (&seq_id, seq_out)| {
            ax.insert(seq_id, seq_out);
            ax
        })
}
