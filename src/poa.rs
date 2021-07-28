use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;

use petgraph::dot::*;
use petgraph::prelude::*;

use crate::*;
pub struct POANode {
    pub nucs: HashMap<SeqID, Nucleotide>,
}
impl std::fmt::Display for POANode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let x = &self
            .nucs
            .values()
            .map(|x| nuc_to_str(x))
            .collect::<Vec<_>>()
            .join(",");
        write!(f, "{}", x)
    }
}
impl std::fmt::Debug for POANode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let x = &self
            .nucs
            .values()
            .map(|x| nuc_to_str(x))
            .collect::<Vec<_>>()
            .join(",");
        write!(f, "{}", x)
    }
}

pub type POAGraph = Graph<POANode, Vec<SeqID>, petgraph::Directed>;
pub type Heads = HashMap<SeqID, NodeIndex>;
pub type Tails = HashMap<SeqID, NodeIndex>;

pub fn update_edge(
    g: &mut POAGraph,
    start: Option<NodeIndex>,
    end: Option<NodeIndex>,
    label: SeqID,
) {
    if !(start.is_some() && end.is_some()) {
        return;
    }
    let start = start.unwrap();
    let end = end.unwrap();

    if let Some(e) = g.find_edge(start, end) {
        g[e].push(label);
    } else {
        g.update_edge(start, end, vec![label]);
    }
}

fn find_head(g: &POAGraph, tail: NodeIndex, seq_id: SeqID) -> NodeIndex {
    let mut node_idx = tail;

    loop {
        if let Some(next_node) = g
            .edges_directed(node_idx, Direction::Outgoing)
            .filter(|e| e.weight().contains(&seq_id))
            .nth(0)
            .map(|e| e.target())
        // Check there is only one?
        {
            node_idx = next_node
        } else {
            break;
        }
    }

    node_idx
}

pub struct POA {
    pub g: POAGraph,
    pub heads: Heads,
    pub tails: Tails,
    pub seq_ids: Vec<SeqID>,
}

impl POA {
    // Create a new POAGraph from a slice of MSABlocks
    pub fn from_nw(a: (POAGraph, HashMap<SeqID, NodeIndex>), seqs: &Sequences) -> Self {
        let mut g = Graph::new();
        let first_node = g.add_node(POANode {
            nucs: HashMap::new(),
        });
        let heads = seqs
            .iter()
            .enumerate()
            .map(|(i, _)| (i, first_node))
            .collect::<HashMap<SeqID, NodeIndex>>();
        let tails = heads.clone();
        let seq_ids = seqs.keys().copied().collect::<Vec<_>>();
        let mut poa = POA {
            g,
            tails,
            heads,
            seq_ids,
        };

        let (h, h_heads) = a;
        poa.extend_with_nw(&h, &h_heads);

        dbg!(&poa.tails);
        for (seq_id, old_tail) in poa.tails.iter_mut() {
            if let Some(tail) = poa
                .g
                .edges_directed(*old_tail, Direction::Outgoing)
                .filter(|e| e.weight().contains(&seq_id))
                .nth(0)
                .map(|e| e.target())
            {
                *old_tail = tail;
                dbg!(&tail);
                dbg!(&poa.g[tail]);
            } else {
                panic!();
            }
        }
        dbg!(&poa.tails);
        for tail in poa.tails.iter() {
            dbg!(tail.1);
            dbg!(&poa.g[*tail.1]);
        }

        // poa.g.remove_node(first_node);
        for tail in poa.tails.iter() {
            dbg!(tail.1);
            dbg!(&poa.g[*tail.1]);
        }
        poa
    }

    pub fn extend_with_nw(&mut self, other: &POAGraph, other_tails: &Tails) {
        let old_nodes = petgraph::algo::toposort(other, None).ok().unwrap();
        let mut old_to_new = HashMap::new();

        // Insert all the other nodes...
        for old_node in old_nodes {
            old_to_new.insert(
                old_node,
                self.g.add_node(POANode {
                    nucs: other
                        .node_weight(old_node)
                        .unwrap()
                        .nucs
                        .iter()
                        .map(|(s, n)| ((*s, n.clone())))
                        .collect(),
                }),
            );
        }
        // ...then all the other edges
        for e in other.raw_edges() {
            self.g.add_edge(
                *old_to_new.get(&e.source()).unwrap(),
                *old_to_new.get(&e.target()).unwrap(),
                e.weight.clone(),
            );
        }

        for (seq_id, old_head) in self.heads.iter_mut() {
            if other_tails.contains_key(&seq_id) {
                update_edge(
                    &mut self.g,
                    Some(*old_head),
                    Some(*old_to_new.get(other_tails.get(&seq_id).unwrap()).unwrap()),
                    *seq_id,
                );
                let new_head = find_head(&self.g, *old_head, *seq_id);
                *old_head = new_head;
            }
        }
    }

    pub fn to_strings(&self) -> HashMap<SeqID, String> {
        let nodes = petgraph::algo::toposort(&self.g, None).ok().unwrap();
        let rank_to_column = nodes
            .iter()
            .enumerate()
            .map(|(i, n)| (n, i))
            .collect::<HashMap<_, _>>();

        self.tails
            .iter()
            .map(|(seq_id, start)| {
                let mut seq_out = vec!["                  ".to_owned(); nodes.len()];
                let mut node = *start;

                loop {
                    let nucs = &self.g[node].nucs;
                    seq_out[rank_to_column[&node]] = if nucs.len() > 1 {
                        nuc_to_str(&nucs[seq_id])
                    } else {
                        nuc_to_str(&nucs[seq_id])
                    };

                    if let Some(next) = self
                        .g
                        .edges_directed(node, Direction::Outgoing)
                        .filter(|e| e.weight().contains(&seq_id))
                        .nth(0)
                        .map(|e| e.target())
                    // Check there is only one?
                    {
                        node = next;
                    } else {
                        break;
                    }
                }

                (seq_id, seq_out)
            })
            .fold(HashMap::new(), |mut ax, (&seq_id, seq_out)| {
                // ax.insert(seq_id, seq_out.iter().collect());
                ax.insert(seq_id, seq_out.join(""));
                ax
            })
    }

    pub fn to_dot(&self, filename: &str) {
        let mut file = File::create(filename).unwrap();
        let _ = file.write_all(format!("{:?}", Dot::new(&self.g)).as_bytes());
    }
}
