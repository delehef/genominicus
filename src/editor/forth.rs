use either::Either;
use itertools::Itertools;
use thiserror::Error;

use super::widgets::treeview::DispGene;

#[derive(Error, Debug)]
pub enum Error {
    #[error("field {0} unknown")]
    FieldUnknown(String),

    #[error("{0}: {1} expects {2} arguments")]
    TooFewElements(String, String, usize),

    #[error("{0}: {1} expects a condition")]
    ConditionExpected(String, String),

    #[error("{0}: {1} expects a value")]
    ValueExpected(String, String),

    #[error("{0}: expected a single final expression, found {1}")]
    SingleValueExpected(String, usize),
}

/// A Combinator operates on boolean expressions
#[derive(Clone, Debug)]
pub enum Combinator {
    And,
    Or,
    Not,
}
impl Combinator {
    fn is_combinator(s: &str) -> bool {
        <&str as TryInto<Combinator>>::try_into(s).is_ok()
    }

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
impl TryFrom<&str> for Combinator {
    type Error = String;
    fn try_from(s: &str) -> Result<Self, Self::Error> {
        match s {
            "&" => Result::Ok(Combinator::And),
            "|" => Result::Ok(Combinator::Or),
            "!" => Result::Ok(Combinator::Not),
            _ => Result::Err("not a combinator".into()),
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
    fn is_relation(s: &str) -> bool {
        <&str as TryInto<Relation>>::try_into(s).is_ok()
    }

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
            "^" => Result::Ok(Relation::StartsWith),
            "$" => Result::Ok(Relation::EndsWith),
            "%" => Result::Ok(Relation::Contains),
            _ => Result::Err("not a function".to_string()),
        }
    }
}
impl ToString for Relation {
    fn to_string(&self) -> String {
        match self {
            Relation::StartsWith => "starts-with",
            Relation::Equal => "equals",
            Relation::EndsWith => "ends-with",
            Relation::Contains => "contains",
        }
        .into()
    }
}

/// Nodes represent the parsed AST, sequentially built from a stack of Tokens.
#[derive(Clone, Debug)]
pub enum ForthExpr {
    Combinator(Combinator, Vec<ForthExpr>),
    Relation(Relation, Vec<ForthExpr>),
    Projector(String),
    Const(String),
}
impl ForthExpr {
    fn is_bool(&self) -> bool {
        matches!(self, ForthExpr::Combinator(..) | ForthExpr::Relation(..))
    }
    fn is_value(&self) -> bool {
        !self.is_bool()
    }
}
impl ForthExpr {
    /// Evaluates an AST at a given position i and returns, if any, the computed
    /// value.
    ///
    /// The computed value may be either Fr or boolean; depending on whether
    /// they stem from a column or a function call, or from a condition or a
    /// combinator. An Either monad encodes this dichotomy.
    pub fn eval(&self, gene: &DispGene) -> Option<Either<String, bool>> {
        match self {
            ForthExpr::Combinator(c, args) => {
                let args = args
                    .iter()
                    .map(|a| a.eval(gene).map(|x| x.right().unwrap()))
                    .collect::<Option<Vec<_>>>();
                args.map(|args| Either::Right(c.apply(&args)))
            }
            ForthExpr::Relation(f, args) => {
                let args = args
                    .iter()
                    .map(|a| a.eval(gene).map(|x| x.left().unwrap()))
                    .collect::<Option<Vec<_>>>();
                args.map(|args| Either::Right(f.apply(&args)))
            }
            // Node::Column(_, column) => project(i, column).map(Either::Left),
            ForthExpr::Projector(field) => match field.as_str() {
                "species" => Some(Either::Left(gene.species.clone())),
                "id" => Some(Either::Left(gene.name.clone())),
                _ => unreachable!(),
            },
            ForthExpr::Const(x) => Some(Either::Left(x.clone())),
        }
    }
}
impl std::fmt::Display for ForthExpr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ForthExpr::Combinator(c, args) => match c {
                Combinator::Not => write!(f, "({} {})", c.to_string(), args[0]),
                _ => {
                    write!(f, "({} {} {})", args[0], c.to_string(), args[1])
                }
            },
            ForthExpr::Relation(ff, args) => {
                write!(f, "({} {} {})", args[0], ff.to_string(), args[1])
            }
            ForthExpr::Projector(field) => write!(f, "gene.{}", field),
            ForthExpr::Const(x) => write!(f, "\"{}\"", x.clone()),
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
fn parse_token(s: &str) -> Result<Token, Error> {
    match s {
        _ if Combinator::is_combinator(s) => Ok(Token::Combinator(s.try_into().unwrap())),
        _ if Relation::is_relation(s) => Ok(Token::Relation(s.try_into().unwrap())),
        _ => {
            if let Some(field) = s.strip_prefix('.') {
                match field {
                    "species" | "id" => Ok(Token::Projector(field.to_owned())),
                    _ => Result::Err(Error::FieldUnknown(field.to_owned())),
                }
            } else {
                Ok(Token::Const(s.to_owned()))
            }
        }
    }
}

fn pretty_stack(stack: &[ForthExpr]) -> String {
    stack.iter().map(|x| x.to_string()).join(" ")
}

/// Pops & returns an argument of a stack, returns an error is none are available
fn take_one(stack: &mut Vec<ForthExpr>, fname: &str) -> Result<ForthExpr, Error> {
    let r1 = stack
        .pop()
        .ok_or_else(|| Error::TooFewElements(pretty_stack(stack), fname.to_owned(), 1))?;
    Ok(r1)
}

/// Pops & returns two arguments of a stack, returns an error is two are not available
fn take_two(stack: &mut Vec<ForthExpr>, fname: &str) -> Result<Vec<ForthExpr>, Error> {
    let r2 = stack
        .pop()
        .ok_or_else(|| Error::TooFewElements(pretty_stack(stack), fname.to_owned(), 2))?;
    let r1 = stack
        .pop()
        .ok_or_else(|| Error::TooFewElements(pretty_stack(stack), fname.to_owned(), 2))?;
    Ok(vec![r1, r2])
}

/// Returns a Node representing the root of the AST parsed from the string representation of a Forth program
pub fn parse(s: &str) -> Result<ForthExpr, Error> {
    let tokens = s.split_whitespace();
    let mut stack = Vec::new();

    for token in tokens.map(parse_token) {
        match token? {
            Token::Combinator(c) => match c {
                Combinator::And | Combinator::Or => {
                    let args = take_two(&mut stack, &c.to_string())?;
                    if !args.iter().all(|n| n.is_bool()) {
                        return Err(Error::ConditionExpected(
                            pretty_stack(&stack),
                            c.to_string(),
                        ));
                    }
                    stack.push(ForthExpr::Combinator(c, args))
                }
                Combinator::Not => {
                    let arg = take_one(&mut stack, &c.to_string())?;
                    if !arg.is_bool() {
                        return Err(Error::ConditionExpected(
                            pretty_stack(&stack),
                            c.to_string(),
                        ));
                    }
                    stack.push(ForthExpr::Combinator(c, vec![arg]))
                }
            },
            Token::Relation(f) => {
                let args = take_two(&mut stack, &f.to_string())?;
                if !args.iter().all(|n| n.is_value()) {
                    return Err(Error::ValueExpected(pretty_stack(&stack), f.to_string()));
                }
                stack.push(ForthExpr::Relation(f, args));
            }
            Token::Const(x) => stack.push(ForthExpr::Const(x)),
            Token::Projector(field) => stack.push(ForthExpr::Projector(field)),
        }
    }

    if stack.is_empty() || stack.len() > 1 {
        return Err(Error::SingleValueExpected(
            pretty_stack(&stack),
            stack.len(),
        ));
    }
    if stack[0].is_value() {
        return Err(Error::ConditionExpected(
            pretty_stack(&stack),
            String::new(),
        ));
    }

    Ok(stack[0].to_owned())
}
