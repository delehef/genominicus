use anyhow::*;
use either::Either;

use super::DispGene;

/// A Combinator operates on boolean expressions
#[derive(Clone, Debug)]
pub enum Combinator {
    And,
    Or,
    Not,
}
impl Combinator {
    fn apply(&self, args: &[bool]) -> bool {
        let a1 = args[0];
        match self {
            Combinator::And => {
                let a2 = args[1];
                a1 && a2
            }
            Combinator::Or => {
                let a2 = args[1];
                a1 || a2
            }
            Combinator::Not => !a1,
        }
    }
}
impl From<&str> for Combinator {
    fn from(s: &str) -> Self {
        match s {
            "&" => Combinator::And,
            "|" => Combinator::Or,
            "!" => Combinator::Not,
            _ => panic!("not a Combinator"),
        }
    }
}
impl std::fmt::Display for Combinator {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Combinator::And => write!(f, "&"),
            Combinator::Or => write!(f, "|"),
            Combinator::Not => write!(f, "!"),
        }
    }
}

/// A Relation operates on String values and convert them to boolean, which can be used with Combinators
#[derive(Clone, Debug)]
pub enum Relation {
    StartsWith,
    Equal,
    EndsWith,
    Contains,
}
impl Relation {
    fn apply(&self, args: &[String]) -> bool {
        match self {
            Relation::StartsWith => args[0].starts_with(&args[1]),
            Relation::Equal => args[0] == args[1],
            Relation::EndsWith => args[0].ends_with(&args[1]),
            Relation::Contains => args[0].contains(&args[1]),
        }
    }
}
impl TryFrom<&str> for Relation {
    type Error = String;
    fn try_from(s: &str) -> Result<Self, Self::Error> {
        match s {
            "=" => Result::Ok(Relation::Equal),
            "<" => Result::Ok(Relation::StartsWith),
            ">" => Result::Ok(Relation::EndsWith),
            "%" => Result::Ok(Relation::Contains),
            _ => Result::Err("not a function".to_string()),
        }
    }
}
impl std::fmt::Display for Relation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Relation::StartsWith => write!(f, "<"),
            Relation::Equal => write!(f, "="),
            Relation::EndsWith => write!(f, ">"),
            Relation::Contains => write!(f, "%"),
        }
    }
}

/// Nodes represent the parsed AST, sequentially built from a stack of Tokens.
#[derive(Clone, Debug)]
pub enum Node {
    Combinator(Combinator, Vec<Node>),
    Relation(Relation, Vec<Node>),
    Projector(String),
    Const(String),
}
impl Node {
    fn is_bool(&self) -> bool {
        matches!(self, Node::Combinator(..) | Node::Relation(..))
    }
    fn is_value(&self) -> bool {
        !self.is_bool()
    }
}
impl Node {
    /// Evaluates an AST at a given position i and returns, if any, the computed
    /// value.
    ///
    /// The computed value may be either Fr or boolean; depending on whether
    /// they stem from a column or a function call, or from a condition or a
    /// combinator. An Either monad encodes this dichotomy.
    pub fn eval(&self, gene: &DispGene) -> Option<Either<String, bool>> {
        match self {
            Node::Combinator(c, args) => {
                let args = args
                    .iter()
                    .map(|a| a.eval(gene).map(|x| x.right().unwrap()))
                    .collect::<Option<Vec<_>>>();
                args.map(|args| Either::Right(c.apply(&args)))
            }
            Node::Relation(f, args) => {
                let args = args
                    .iter()
                    .map(|a| a.eval(gene).map(|x| x.left().unwrap()))
                    .collect::<Option<Vec<_>>>();
                args.map(|args| Either::Right(f.apply(&args)))
            }
            // Node::Column(_, column) => project(i, column).map(Either::Left),
            Node::Projector(field) => match field.as_str() {
                "species" => Some(Either::Left(gene.species.clone())),
                "id" => Some(Either::Left(gene.name.clone())),
                _ => unreachable!(),
            },
            Node::Const(x) => Some(Either::Left(x.clone())),
        }
    }
}
impl std::fmt::Display for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Node::Combinator(c, args) => match c {
                Combinator::Not => write!(f, "({} {})", c, args[0]),
                _ => {
                    write!(f, "({} {} {})", c, args[0], args[1])
                }
            },
            // Node::Comparison(r, args) => write!(f, "({} {} {})", r, args[0], args[1]),
            Node::Relation(ff, args) => write!(f, "({} {} {})", ff, args[0], args[1]),
            Node::Projector(field) => write!(f, "gene.{}", field),
            Node::Const(x) => write!(f, "{}", x),
        }
    }
}

/// Tokens are used to parse a sequence of Forth tokens from a string
enum Token {
    Combinator(Combinator),
    Relation(Relation),
    Projector(String),
    Const(String),
}
fn parse_token(s: &str) -> Result<Token> {
    match s {
        "&" | "|" | "!" => Ok(Token::Combinator(s.into())),
        "<" | "=" | ">" | "%" => Ok(Token::Relation(s.try_into().unwrap())),
        _ => {
            if let Some(field) = s.strip_prefix('.') {
                match field {
                    "species" | "id" => Ok(Token::Projector(field.to_owned())),
                    _ => bail!("unknown field: {}", field),
                }
            } else {
                Ok(Token::Const(s.to_owned()))
            }
        }
    }
}

/// Pops & returns an argument of a stack, returns an error is none are available
fn take_one(stack: &mut Vec<Node>, fname: &str) -> Result<Node> {
    let r1 = stack
        .pop()
        .ok_or_else(|| anyhow!("{} expects an argument", fname))?;
    Ok(r1)
}

/// Pops & returns two arguments of a stack, returns an error is two are not available
fn take_two(stack: &mut Vec<Node>, fname: &str) -> Result<Vec<Node>> {
    let r2 = stack
        .pop()
        .ok_or_else(|| anyhow!("{} expects two arguments", fname))?;
    let r1 = stack
        .pop()
        .ok_or_else(|| anyhow!("{} expects two arguments", fname))?;
    Ok(vec![r1, r2])
}

/// Returns a Node representing the root of the AST parsed from the string representation of a Forth program
pub fn parse(s: &str) -> Result<Node> {
    let tokens = s.split_whitespace();
    let mut stack = Vec::new();

    for token in tokens.map(parse_token) {
        match token? {
            Token::Combinator(c) => match c {
                Combinator::And | Combinator::Or => {
                    let args = take_two(&mut stack, &c.to_string())?;
                    if !args.iter().all(|n| n.is_bool()) {
                        bail!("{} expects conditions", c.to_string());
                    }
                    stack.push(Node::Combinator(c, args))
                }
                Combinator::Not => {
                    let arg = take_one(&mut stack, &c.to_string())?;
                    if !arg.is_bool() {
                        bail!("{} expects conditions", c.to_string());
                    }
                    stack.push(Node::Combinator(c, vec![arg]))
                }
            },
            Token::Relation(f) => {
                let args = take_two(&mut stack, &f.to_string())?;
                if !args.iter().all(|n| n.is_value()) {
                    bail!("{} expects values", f.to_string());
                }
                stack.push(Node::Relation(f, args));
            }
            Token::Const(x) => stack.push(Node::Const(x)),
            Token::Projector(field) => stack.push(Node::Projector(field)),
        }
    }

    if stack.is_empty() || stack.len() > 1 {
        bail!("expected a single final expression")
    }
    if stack[0].is_value() {
        bail!("expected a comparison")
    }

    Ok(stack[0].to_owned())
}
