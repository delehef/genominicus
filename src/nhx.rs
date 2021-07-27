use pest::Parser;
use std::{collections::HashMap, usize};

use pest_derive::Parser;
#[derive(Parser)]
#[grammar = "nhx.pest"]
pub struct NhxParser;

#[derive(Debug)]
pub struct Node {
    pub id: usize,
    pub name: Option<String>,
    parent: usize,
    pub children: Vec<usize>,
    length: Option<f32>,
    data: HashMap<String, String>,
}
impl Node {
    pub fn new(parent: usize, id: usize) -> Self {
        Node {
            id,
            name: None,
            parent,
            children: Vec::new(),
            length: None,
            data: HashMap::new(),
        }
    }

    pub fn is_duplication(&self) -> bool {
        self.data.get("D").map_or(false, |d| d == "Y")
    }

    pub fn children(&self) -> &[usize] {
        &self.children
    }
    pub fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }
}

pub struct Tree {
    nodes: Vec<Node>,
}
impl Tree {
    pub fn print(&self) {
        fn print_node(nodes: &[Node], n: usize, o: usize) {
            println!(
                "{}{}:{:?} - {:?}",
                str::repeat(" ", o),
                &nodes[n].name.as_ref().unwrap_or(&String::new()),
                &nodes[n].length.unwrap_or(-1.),
                &nodes[n].data
            );
            for &c in &nodes[n].children {
                print_node(nodes, c, o + 2)
            }
        }
        print_node(&self.nodes, 0, 0);
    }

    pub fn d_descendants(&self, n: usize) -> Vec<usize> {
        fn find_dups(t: &Tree, n: usize, ax: &mut Vec<usize>) {
            if t[n].is_duplication() {
                ax.push(t[n].id);
            } else {
                for &c in t[n].children.iter() {
                    find_dups(t, c, ax);
                }
            }
        }

        let mut r = vec![];
        for &c in self[n].children.iter() {
            find_dups(self, c, &mut r);
        }
        r
    }

    pub fn non_d_descendants(&self, n: usize) -> Vec<usize> {
        fn find_leaves(t: &Tree, n: usize, ax: &mut Vec<usize>) {
            if !t[n].is_duplication() {
                if t[n].is_leaf() {
                    ax.push(t[n].id);
                }
                for &c in t[n].children.iter() {
                    find_leaves(t, c, ax);
                }
            }
        }

        let mut r = vec![];
        find_leaves(self, n, &mut r);
        r
    }

    pub fn from_string(content: &str) -> Result<Self, pest::error::Error<Rule>> {
        use pest::iterators::Pair;

        fn parse_attrs(pair: Pair<Rule>, me: &mut Node) {
            match pair.as_rule() {
                Rule::float => me.length = Some(pair.as_str().parse::<f32>().unwrap()),
                Rule::NhxEntry => {
                    let mut kv = pair.into_inner();
                    let k = kv.next().unwrap().as_str().to_owned();
                    let v = kv.next().unwrap().as_str().to_owned();
                    me.data.insert(k, v);
                }
                _ => {
                    unimplemented!();
                }
            }
        }

        fn parse_node(pair: Pair<Rule>, parent: usize, storage: &mut Vec<Node>) -> usize {
            let my_id = storage.len();
            storage.push(Node::new(parent, my_id));

            pair.into_inner().for_each(|inner| match inner.as_rule() {
                Rule::Leaf | Rule::Clade => {
                    let child = parse_node(inner, my_id, storage);
                    storage[my_id].children.push(child);
                }
                Rule::name => {
                    storage[my_id].name = Some(inner.as_str().to_owned());
                }
                Rule::Attributes => {
                    for attr in inner.into_inner() {
                        parse_attrs(attr, &mut storage[my_id])
                    }
                }
                _ => unimplemented!(),
            });

            my_id
        }

        let root = NhxParser::parse(Rule::Tree, &content)?.next().unwrap();

        let mut r = Vec::new();
        let _ = parse_node(root, 0, &mut r);
        Ok(Tree { nodes: r })
    }

    pub fn from_filename(filename: &str) -> Result<Self, pest::error::Error<Rule>> {
        let content = std::fs::read_to_string(filename).expect("cannot read file");
        Self::from_string(&content)
    }

    pub fn depth(&self) -> f32 {
        todo!()
    }

    pub fn node_depth(&self, n: usize) -> f32 {
        let mut depth = self.nodes[n].length.unwrap();
        let mut parent = self.nodes[n].parent;
        while parent != 0 {
            depth += self.nodes[parent].length.unwrap();
            parent = self.nodes[parent].parent;
        }
        depth
    }

    pub fn node_topological_depth(&self, n: usize) -> f32 {
        let mut depth = 0.;
        let mut parent = self.nodes[n].parent;
        while parent != 0 {
            depth += 1.;
            parent = self.nodes[parent].parent;
        }
        depth
    }

    pub fn topological_depth(&self) -> (usize, f32) {
        self.leaves()
            .map(|n| (n, self.node_topological_depth(n)))
            .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
            .unwrap()
    }

    pub fn leaves<'a>(&'a self) -> impl Iterator<Item = usize> + 'a {
        (0..self.nodes.len()).filter(move |n| self.nodes[*n].children.is_empty())
    }

    pub fn leaf_names(&self) -> Vec<(usize, Option<&String>)> {
        self.leaves()
            .map(|n| (n, self.nodes[n].name.as_ref()))
            .collect()
    }
}
impl std::ops::Index<usize> for Tree {
    type Output = Node;
    fn index(&self, i: usize) -> &Self::Output {
        &self.nodes[i]
    }
}
