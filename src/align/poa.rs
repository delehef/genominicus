use std::collections::HashMap;

use petgraph::prelude::*;

use super::*;

pub struct POANode {
    pub nucs: HashMap<SeqID, Nucleotide>,
}
impl std::fmt::Display for POANode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let x = &self
            .nucs
            .values()
            .map(|x| x.to_string())
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
            .map(|x| x.to_string())
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

    while let Some(next_node) = g
        .edges_directed(node_idx, Direction::Outgoing)
        .find(|e| e.weight().contains(&seq_id))
        .map(|e| e.target())
    // Check there is only one?
    {
        node_idx = next_node
    }

    node_idx
}
