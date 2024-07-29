use newick::{Newick, NewickTree, NodeID};
use ratatui::{
    layout::{Constraint, Margin, Rect},
    style::{Color, Style, Stylize},
    text::{Line, Span},
    widgets::{Cell, Row, Scrollbar, ScrollbarOrientation, ScrollbarState, Table, TableState},
    Frame,
};
use std::{collections::HashMap, ops::Range, rc::Rc, sync::OnceLock};
use syntesuite::genebook::Gene;

use crate::{
    editor::forth::ForthExpr, name2color, shiftreg::ShiftRegister, ColorMap, GeneCache, WINDOW,
};

const BLOCKS: &[Range<u32>] = &[
    // Greek letters
    0x391..0x39f,
    0x3b0..0x3ff,
    // Cyrillic
    0x400..0x44f,
    // Armenian
    0x531..0x54f,
    // // Shapes -- tend to be too wide in most fonts
    // 0x25a0..0x25ff,
];
static GENABET: OnceLock<Vec<char>> = OnceLock::new();

#[derive(Clone, Copy)]
pub struct TreeViewSettings {
    pub use_symbols: bool,
}

fn family_to_char(id: usize) -> char {
    let chars_len = GENABET.get().unwrap().len();

    GENABET.get().unwrap()[id % chars_len]
}

fn gene_to_char(family: usize, strand: syntesuite::Strand, symbol: bool) -> char {
    if symbol {
        family_to_char(family)
    } else {
        match strand {
            syntesuite::Strand::Direct => '▶',
            syntesuite::Strand::Reverse => '◀',
            syntesuite::Strand::Unknown => '■',
        }
    }
}

#[derive(Default, Clone)]
struct FoldingPoint {
    clade: Vec<Vec<NodeID>>,
    point: usize,
}
impl FoldingPoint {
    fn fold(&mut self) -> Option<&[NodeID]> {
        if self.point == self.clade.len() {
            return None;
        } else {
            self.point += 1;
            return Some(&self.clade[self.point - 1]);
        }
    }

    fn unfold(&mut self) -> Option<&[NodeID]> {
        if self.point == 0 {
            return None;
        } else {
            self.point -= 1;
            return Some(&self.clade[self.point]);
        }
    }
}

#[derive(Debug, Clone)]
pub struct DispGene {
    pub name: String,
    pub species: String,
}

pub struct LandscapeData {
    pub book: GeneCache,
    pub colors: ColorMap,
}

#[derive(Debug)]
enum Clade {
    Taxon {
        graph_line: usize,
        dup_nesting: Vec<f32>,
        gene: DispGene,
        id: NodeID,
    },
    SubClade {
        subclades: Vec<usize>,
        folded: bool,
        id: NodeID,
    },
}

#[derive(Debug)]
struct CladeHierarchy {
    clades: Vec<Clade>,
}
impl CladeHierarchy {
    fn new() -> Self {
        Self {
            clades: vec![Clade::SubClade {
                subclades: vec![],
                folded: false,
                id: 1,
            }],
        }
    }
    fn append_in(&mut self, clade: Clade, parent: usize) -> usize {
        let new = self.clades.len();
        self.clades.push(clade);
        if let Clade::SubClade {
            ref mut subclades, ..
        } = &mut self.clades[parent]
        {
            subclades.push(new);
        } else {
            unreachable!()
        }
        new
    }

    pub fn is_folded(&self, i: usize) -> bool {
        match &self.clades[i] {
            Clade::Taxon { .. } => false,
            Clade::SubClade { folded, .. } => *folded,
        }
    }

    fn get(&self, i: usize) -> &Clade {
        &self.clades[i]
    }

    fn get_mut(&mut self, i: usize) -> &mut Clade {
        &mut self.clades[i]
    }

    fn find_first_taxon(&self, i: usize) -> &Clade {
        match &self.clades[i] {
            taxon @ Clade::Taxon { .. } => taxon,
            Clade::SubClade { subclades, .. } => self.find_first_taxon(subclades[0]),
        }
    }
}

const DEPTH_FACTOR: usize = 2;

#[derive(PartialEq, Eq)]
enum Position {
    First,
    Last,
    Other,
}
struct NodeContext {
    id: NodeID,
    depth: i64,
    position: Position,
}

#[derive(Default)]
struct DuplicationsCache {
    nestings: HashMap<NodeID, Vec<DupNesting>>,
    max_nesting: usize,
}

#[derive(Default)]
struct FoldCache {
    fold_level: HashMap<NodeID, usize>,
    folding_points: Vec<FoldingPoint>,
}

struct Caches {
    genes: HashMap<NodeID, DispGene>,
    lineages: HashMap<NodeID, Vec<NodeContext>>,
    tree: HashMap<NodeID, String>,
    duplications: DuplicationsCache,
    folding: FoldCache,
}

struct States {
    gene_table: TableState,
    scrollbar: ScrollbarState,
}
impl States {
    fn new(size: usize) -> Self {
        States {
            gene_table: TableState::new().with_selected(0),
            scrollbar: ScrollbarState::new(size - 1),
        }
    }
}

enum DupNesting {
    Head(f32),
    Body(f32),
    Tail(f32),
}
impl DupNesting {
    fn score(&self) -> f32 {
        match self {
            DupNesting::Head(x) | DupNesting::Body(x) | DupNesting::Tail(x) => *x,
        }
    }

    fn to_span(&self) -> Span {
        let score = self.score();
        Span::from(match self {
            DupNesting::Head(_) => "┬",
            DupNesting::Body(_) => "│",
            DupNesting::Tail(_) => "┴",
        })
        .fg(if score < 0. {
            Color::LightBlue
        } else {
            Color::Rgb(
                ((1. - score) * 255.0).ceil() as u8,
                (score * 255.).ceil() as u8,
                0,
            )
        })
    }
}

pub struct TreeView {
    cache: Caches,
    tree: Rc<NewickTree>,
    pub settings: TreeViewSettings,
    landscape_data: Option<LandscapeData>,
    current_len: usize,
    // screen coordinate -> inner nodes IDs
    screen_to_nodes: HashMap<usize, Vec<usize>>,
    pub highlighters: Vec<ForthExpr>,
    states: States,
}
impl TreeView {
    pub fn from_newick(
        tree: Rc<NewickTree>,
        settings: TreeViewSettings,
        landscape_data: Option<LandscapeData>,
    ) -> Self {
        let _ = GENABET.get_or_init(|| {
            BLOCKS
                .iter()
                .flat_map(|b| b.clone())
                .map(|c| char::from_u32(c).unwrap())
                .collect()
        });

        let leave_count = tree.leaves().count();

        let genes = tree
            .leaves()
            .map(|n| {
                (
                    n,
                    DispGene {
                        name: tree.name(n).cloned().unwrap_or("UNKNWN".into()),
                        species: tree.attrs(n).get("S").cloned().unwrap_or("UNKNWN".into()),
                    },
                )
            })
            .collect();

        let lineages = tree
            .leaves()
            .map(|n| {
                (
                    n,
                    tree.ascendance(n)
                        .into_iter()
                        .map(|n| NodeContext {
                            id: n,
                            depth: tree.node_topological_depth(n).unwrap(),
                            position: {
                                if let Some(parent) = tree.parent(n) {
                                    if tree.children(parent).unwrap()[0] == n {
                                        Position::First
                                    } else {
                                        Position::Last
                                    }
                                } else {
                                    Position::Last
                                }
                            },
                        })
                        .collect(),
                )
            })
            .collect();

        let mut r = Self {
            cache: Caches {
                genes,
                lineages,
                tree: Default::default(),
                duplications: Default::default(),
                folding: Default::default(),
            },
            tree,
            landscape_data,
            settings,
            current_len: leave_count,
            screen_to_nodes: Default::default(),
            highlighters: Vec::new(),
            states: States::new(leave_count),
        };
        r.cache_tree_graph();
        r.cache_dup_nesting();
        r.cache_folding();
        r
    }

    pub fn len(&self) -> usize {
        self.current_len
    }

    fn cache_tree_graph(&mut self) {
        self.cache.tree = self
            .tree
            .leaves()
            .map(|n| (n, self.make_tree_line(n)))
            .collect();
    }

    fn cache_folding(&mut self) {
        self.cache.folding.folding_points = Vec::with_capacity(self.tree.len());
        for (y, n) in self.tree.leaves().enumerate() {
            self.cache.folding.folding_points.push(FoldingPoint {
                clade: self.cache.lineages[&n]
                    .iter()
                    .map(|a| self.tree.leaves_of(a.id))
                    .collect(),
                point: 0,
            });
            self.cache.folding.fold_level.insert(n, 0);
        }
    }

    fn cache_dup_nesting(&mut self) {
        self.cache.duplications.nestings.clear();
        for n in self.tree.leaves() {
            // let mut pure_head = ShiftRegister::new(1);
            let mut pure_head_broken = false;
            let mut pure_tail_broken = false;
            let mut pure_head = ShiftRegister::new(3, false);
            let mut pure_tail = ShiftRegister::new(3, false);
            let dup_nesting = self.cache.lineages[&n]
                .iter()
                .filter_map(|n| {
                    if n.position != Position::First {
                        pure_head_broken = true;
                    }
                    if n.position != Position::Last {
                        pure_tail_broken = true;
                    }
                    pure_tail.write(n.position == Position::Last && !pure_tail_broken);
                    pure_head.write(n.position == Position::First && !pure_head_broken);

                    if self.tree.is_duplication(n.id) {
                        let dcs = self
                            .tree
                            .attrs(n.id)
                            .get("DCS")
                            .map(|x| x.parse::<f32>().unwrap())
                            .unwrap_or_default();

                        Some(if *pure_head.read() {
                            DupNesting::Head(dcs)
                        } else if *pure_tail.read() {
                            DupNesting::Tail(dcs)
                        } else {
                            DupNesting::Body(dcs)
                        })
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>();

            self.cache.duplications.max_nesting =
                self.cache.duplications.max_nesting.max(dup_nesting.len());
            self.cache.duplications.nestings.insert(n, dup_nesting);
        }
    }

    fn make_tree_line(&self, n: NodeID) -> String {
        let lineage = &self.cache.lineages[&n];
        let last_branch_length =
            self.tree.topological_depth().0 - self.tree.node_topological_depth(n).unwrap() as usize;
        let mut r = "─".repeat(last_branch_length);

        let mut on_my_line = true;
        for n in lineage {
            let is_duplication = self
                .tree
                .parent(n.id)
                .map(|x| self.tree.is_duplication(x))
                .unwrap_or(false);

            if on_my_line {
                if n.position != Position::First {
                    on_my_line = false;
                }

                if on_my_line {
                    if is_duplication {
                        r.push_str("─D");
                    } else {
                        r.push_str("─┬");
                    }
                } else {
                    r.push_str(if is_duplication { "─╙" } else { "─└" })
                }
            } else {
                if n.position == Position::First {
                    if is_duplication {
                        r.push_str(" ║");
                    } else {
                        r.push_str(" │");
                    }
                } else {
                    r.push_str("  ");
                }
            }
        }

        r.chars().rev().collect()
    }

    fn gene_to_row<'a>(
        graph_line: &'a str,
        landscape_data: Option<&'a LandscapeData>,
        gene: DispGene,
        dups_nesting: &'a [DupNesting],
        with_fold_indicator: bool,
        use_symbols: bool,
        highlighters: &[ForthExpr],
    ) -> Row<'a> {
        const HL_COLORS: [Color; 7] = [
            Color::LightBlue,
            Color::LightRed,
            Color::LightCyan,
            Color::LightGreen,
            Color::LightYellow,
            Color::LightMagenta,
            Color::Gray,
        ];
        let landscape = if let Some(Gene {
            strand,
            left_landscape,
            right_landscape,
            family,
            ..
        }) =
            landscape_data.and_then(|landscape_data| landscape_data.book.get(&gene.name))
        {
            Line::from_iter(
                std::iter::once(Span::from("- ".repeat(WINDOW - left_landscape.len())).dark_gray())
                    .chain(left_landscape.iter().map(|g| {
                        Span::from(format!(
                            "{} ",
                            gene_to_char(g.family, g.strand, use_symbols)
                        ))
                        .fg({
                            let color = landscape_data
                                .unwrap()
                                .colors
                                .get(&g.family)
                                .unwrap()
                                .to_percent();

                            Color::Rgb(
                                (color.0 * 255.0).floor() as u8,
                                (color.1 * 255.0).floor() as u8,
                                (color.2 * 255.0).floor() as u8,
                            )
                        })
                    }))
                    .chain(
                        std::iter::once(
                            format!(" {} ", gene_to_char(*family, *strand, use_symbols)).fg({
                                let color = landscape_data
                                    .unwrap()
                                    .colors
                                    .get(family)
                                    .unwrap()
                                    .to_percent();
                                Color::Rgb(
                                    (color.0 * 255.0).floor() as u8,
                                    (color.1 * 255.0).floor() as u8,
                                    (color.2 * 255.0).floor() as u8,
                                )
                            }),
                        )
                        .chain(right_landscape.iter().map(|g| {
                            Span::from(format!(
                                " {}",
                                gene_to_char(g.family, g.strand, use_symbols)
                            ))
                            .fg({
                                let color = landscape_data
                                    .unwrap()
                                    .colors
                                    .get(&g.family)
                                    .unwrap()
                                    .to_percent();
                                Color::Rgb(
                                    (color.0 * 255.0).floor() as u8,
                                    (color.1 * 255.0).floor() as u8,
                                    (color.2 * 255.0).floor() as u8,
                                )
                            })
                        }))
                        .chain(std::iter::once(
                            Span::from(" -".repeat(WINDOW - right_landscape.len())).dark_gray(),
                        )),
                    ),
            )
        } else {
            Line::from("")
        };

        let highlighted = highlighters
            .iter()
            .enumerate()
            .filter_map(|(i, h)| {
                if h.eval(&gene).unwrap().right().unwrap() {
                    Some(i)
                } else {
                    None
                }
            })
            .next();
        let species_color = name2color(&gene.species).to_percent();
        Row::new(vec![
            if with_fold_indicator {
                Cell::from("⋮".to_string()).bold()
            } else {
                "".into()
            },
            graph_line.into(),
            Cell::from(Line::from(
                dups_nesting
                    .iter()
                    .rev()
                    .map(|x| x.to_span())
                    .collect::<Vec<_>>(),
            )),
            gene.species
                .fg(Color::Rgb(
                    (species_color.0 * 255.0).floor() as u8,
                    (species_color.1 * 255.0).floor() as u8,
                    (species_color.2 * 255.0).floor() as u8,
                ))
                .into(),
            if let Some(i) = highlighted {
                gene.name
                    .clone()
                    .bold()
                    .fg(HL_COLORS[i % HL_COLORS.len()])
                    .reversed()
                    .into()
            } else {
                gene.name.clone().into()
            },
            landscape.into(),
        ])
    }

    fn to_rows(&mut self, f: &mut Frame, t: Rect) {
        self.screen_to_nodes.clear();

        let mut rows = Vec::new();
        let mut y = 0;
        for n in self.tree.leaves() {
            let fold_level = *self.cache.folding.fold_level.get(&n).unwrap_or(&0);
            let folded = fold_level > 0;
            let lineage_len = self.cache.lineages[&n].len();
            let first_in_fold = folded
                && self.cache.lineages[&n][lineage_len - fold_level].position == Position::First;
            if !folded || first_in_fold {
                let ancestors = self.cache.lineages[&n]
                    .iter()
                    .map(|n| n.id)
                    .collect::<Vec<_>>();
                self.screen_to_nodes.insert(y, ancestors);
                let row = Self::gene_to_row(
                    &self.cache.tree[&n],
                    self.landscape_data.as_ref(),
                    self.cache.genes.get(&n).unwrap().clone(),
                    &self.cache.duplications.nestings[&n],
                    false,
                    self.settings.use_symbols,
                    &self.highlighters,
                );
                rows.push(row);
                y += 1;
            }
        }
        self.current_len = rows.len();

        let tree_depth = self.tree.topological_depth().1;
        let widths = [
            Constraint::Length(1),
            Constraint::Length((DEPTH_FACTOR * tree_depth) as u16),
            Constraint::Length((self.cache.duplications.max_nesting).try_into().unwrap()),
            Constraint::Fill(1),
            Constraint::Fill(1),
            Constraint::Fill(3),
        ];

        let table = Table::new(rows, widths)
            .column_spacing(1)
            .header(
                Row::new(vec!["", "", "Dup.", "Species", "Gene", "Synteny"])
                    .style(Style::new().bold())
                    .bottom_margin(1),
            )
            .highlight_symbol(">>")
            .highlight_style(Style::new().underlined());
        f.render_stateful_widget(table, t, &mut self.states.gene_table);
    }

    pub fn toggle_current(&mut self) {
        // let screen_y = self.states.gene_table.selected().unwrap();

        // let target_state = !self.screen_to_clade[&screen_y]
        //     .iter()
        //     .any(|c| self.clades.is_folded(*c));

        // for clade in self.screen_to_clade.get(&screen_y).unwrap().iter().rev() {
        //     if let Clade::SubClade { ref mut folded, .. } = self.clades.get_mut(*clade) {
        //         *folded = target_state;
        //     } else {
        //         unreachable!()
        //     }
        // }
    }

    pub fn fold_current(&mut self) {
        // let screen_y = self.states.gene_table.selected().unwrap();
        // if let Some(leaves) = self.cache.folding.folding_points[screen_y].fold() {
        //     for l in leaves {
        //         self.cache
        //             .folding
        //             .fold_level
        //             .entry(*l)
        //             .and_modify(|x| *x += 1);
        //     }
        // }

        // for clade in self.screen_to_clade.get(&screen_y).unwrap().iter().rev() {
        //     if let Clade::SubClade { ref mut folded, .. } = self.clades.get_mut(*clade) {
        //         if !*folded {
        //             *folded = true;
        //             return;
        //         }
        //     } else {
        //         unreachable!()
        //     }
        // }
    }

    pub fn unfold_current(&mut self) {
        let screen_y = self.states.gene_table.selected().unwrap();
        if let Some(leaves) = self.cache.folding.folding_points[screen_y].unfold() {
            for l in leaves {
                self.cache
                    .folding
                    .fold_level
                    .entry(*l)
                    .and_modify(|x| *x -= 1);
            }
        }
        // for clade in self.screen_to_clade.get(&screen_y).unwrap().iter() {
        //     if let Clade::SubClade { ref mut folded, .. } = self.clades.get_mut(*clade) {
        //         if *folded {
        //             *folded = false;
        //             return;
        //         }
        //     } else {
        //         unreachable!()
        //     }
        // }
    }

    // fn max_dup_nesting(&self) -> u16 {
    //     // self.dup_level.iter().map(Vec::len).max().unwrap_or(0) as u16
    // }

    pub fn move_to(&mut self, i: usize) {
        self.states.gene_table.select(Some(i));
        self.states.scrollbar = self.states.scrollbar.position(i);
    }

    pub fn prev(&mut self, count: usize) {
        let i = self
            .states
            .gene_table
            .selected()
            .map(|i| i.saturating_sub(count))
            .unwrap_or_default();
        self.move_to(i);
    }

    pub fn next(&mut self, count: usize) {
        let i = self
            .states
            .gene_table
            .selected()
            .map(|i| (i + count).clamp(0, self.len() - 1))
            .unwrap_or_default();
        self.move_to(i);
    }

    pub fn top(&mut self) {
        self.move_to(0);
    }

    pub fn bottom(&mut self) {
        self.move_to(self.len() - 1);
    }

    pub fn render(&mut self, f: &mut Frame, t: Rect) {
        self.to_rows(f, t);

        f.render_stateful_widget(
            Scrollbar::default()
                .orientation(ScrollbarOrientation::VerticalRight)
                .begin_symbol(Some("^"))
                .end_symbol(Some("v")),
            t.inner(Margin {
                vertical: 1,
                horizontal: 1,
            }),
            &mut self.states.scrollbar,
        );
    }
}
