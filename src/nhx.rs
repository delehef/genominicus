use std::collections::HashMap;
use pest::Parser;
#[derive(Parser)]
#[grammar = "nhx.pest"]
pub struct NhxParser;

pub struct Node {
    name: String,
    children: Vec<usize>,
    length: f32,
    data: HashMap<String, String>
}

pub fn parse_nhx_string(content: &str) -> Result<TreeValue, pest::error::Error<Rule>> {
    use pest::iterators::Pair;
    let mut tree = Vec<Node>::new();

    fn parse_value(pair: Pair<Rule>) -> TreeValue {
        match pair.as_rule() {
            Rule::Leaf => {
                let name = pair.into_inner().next().unwrap().as_str();
                TreeValue::Node {
                    name: Some(name.into()),
                    children: None,
                }
            }
            Rule::Internal => {
                let mut inner_rules = pair.into_inner();
                let children = Some(
                    inner_rules
                        .next()
                        .unwrap()
                        .into_inner()
                        .map(parse_value)
                        .collect(),
                );
                let name = if let Some(clade) = inner_rules.next() {
                    Some(clade.as_str().into())
                } else {
                    None
                };
                TreeValue::Node { children, name }
            }

            Rule::Branch => {
                fn get_weight(mut inner: pest::iterators::Pairs<Rule>) -> f32 {
                    if let Some(weight) = inner.next() {
                        weight.as_str().parse::<f32>().unwrap()
                    } else {
                        f32::NAN
                    }
                }

                let mut inner = pair.into_inner();
                let (node, weight) = if let Some(next) = inner.next() {
                    match next.as_rule() {
                        Rule::SubTree => (parse_value(next), get_weight(inner)),
                        _ => (
                            TreeValue::Node {
                                name: None,
                                children: None,
                            },
                            next.as_str().parse::<f32>().unwrap(),
                        ),
                    }
                } else {
                    (
                        TreeValue::Node {
                            name: None,
                            children: None,
                        },
                        get_weight(inner),
                    )
                };

                TreeValue::Link {
                    weight,
                    node: Box::new(node),
                }
            }

            Rule::SubTree => parse_value(pair.into_inner().next().unwrap()),
            Rule::EOI
                | Rule::WHITESPACE
                | Rule::Tree
                | Rule::Length
                | Rule::BranchSet
                | Rule::float
                | Rule::safe
                | Rule::name => panic!("WEIRD"),
        }
    }

    let root = NewickParser::parse(Rule::Tree, &content)?.next().unwrap();

    Ok(parse_value(root))
}

pub fn parse_nhx_file(file: &str) -> Result<TreeValue, pest::error::Error<Rule>> {
    let content = std::fs::read_to_string(file).expect("cannot read file");
    parse_nhx_string(&content)
}
