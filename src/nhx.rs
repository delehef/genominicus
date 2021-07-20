use std::{collections::HashMap, usize};
use pest::Parser;

use pest_derive::Parser;
#[derive(Parser)]
#[grammar = "nhx.pest"]
pub struct NhxParser;

#[derive(Debug)]
pub struct Node {
    name: Option<String>,
    parent: usize,
    children: Vec<usize>,
    length: Option<f32>,
    data: HashMap<String, String>,
}
impl Node {
    pub fn new(parent: usize) -> Self {
        Node {
            name: None,
            parent,
            children: Vec::new(),
            length: None,
            data: HashMap::new(),
        }
    }
}

pub struct Tree {
    nodes: Vec<Node>,
}

impl Tree {
    pub fn print(&self) {
        fn print_node(nodes: &Vec<Node>, n: usize, o: usize) {
            println!("{}{}:{:?} - {:?}",
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

    pub fn from_string(content: &str) -> Result<Self, pest::error::Error<Rule>> {
        use pest::iterators::Pair;

        fn parse_attrs(pair: Pair<Rule>, me: &mut Node) {
            match pair.as_rule() {
                Rule::float => {
                    me.length = Some(pair.as_str().parse::<f32>().unwrap())
                }
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
            storage.push(Node::new(parent));

            pair.into_inner().for_each(|inner|  {
                match inner.as_rule() {
                    Rule::Leaf | Rule::Clade => {
                        let child = parse_node(inner, 0, storage);
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
                    _ => unimplemented!()
                }
            });

            my_id
        }

        let root = NhxParser::parse(Rule::Tree, &content)?.next().unwrap();

        let mut r = Vec::new();
        let _ = parse_node(root, 0, &mut r);
        Ok(Tree{ nodes: r })
    }

    pub fn from_filename(filename: &str) -> Result<Self, pest::error::Error<Rule>> {
        let content = std::fs::read_to_string(filename).expect("cannot read file");
        Self::from_string(&content)
    }
}
