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
    editor::forth::ForthExpr,
    shiftreg::ShiftRegister,
    utils::{name2color, ColorMap, GeneCache, WINDOW},
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
/// An “alphabet” of symbols to represent syntenic families.
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

#[derive(Debug, Clone)]
pub struct DispGene {
    pub name: String,
    pub species: String,
}

pub struct LandscapeData {
    pub book: GeneCache,
    pub colors: ColorMap,
}

const DEPTH_FACTOR: usize = 2;

#[derive(PartialEq, Eq, Debug)]
enum Position {
    UpperBranch,
    MiddleBranch,
    LowerBranch,
}
struct NodeContext {
    id: NodeID,
    position: Position,
}

#[derive(Default)]
struct DuplicationsCache {
    nestings: HashMap<NodeID, Vec<DupNesting>>,
    max_nesting: usize,
}

#[derive(Default)]
struct Caches {
    narrowed_tree: Rc<NewickTree>,
    genes: HashMap<NodeID, DispGene>,
    lineages: HashMap<NodeID, Vec<NodeContext>>,
    tree_chars: HashMap<NodeID, String>,
    duplications: DuplicationsCache,
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

    fn to_span(&'_ self) -> Span<'_> {
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
    tree: NewickTree,
    pub settings: TreeViewSettings,
    landscape_data: Option<LandscapeData>,
    /// screen coordinate -> inner nodes IDs
    screen_to_nodes: HashMap<usize, Vec<usize>>,
    /// A list of selectors to highlight the matching genes
    pub highlighters: Vec<ForthExpr>,
    /// A list of filters to focus on selected clades/genes
    pub narrowing: Option<ForthExpr>,
    /// UI state
    states: States,
}
impl TreeView {
    pub fn from_newick(
        tree: NewickTree,
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

        let leaves_count = tree.len();
        let mut r = Self {
            cache: Caches::default(),
            tree,
            landscape_data,
            settings,
            screen_to_nodes: Default::default(),
            highlighters: Default::default(),
            narrowing: Default::default(),
            states: States::new(leaves_count),
        };
        r.update_caches();
        r
    }

    fn update_caches(&mut self) {
        let tree = self.tree.clone();

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
                            position: {
                                if let Some(parent) = tree.parent(n) {
                                    if tree
                                        .children(parent)
                                        .unwrap()
                                        .first()
                                        .map(|first| *first == n)
                                        .unwrap_or(false)
                                    {
                                        Position::UpperBranch
                                    } else if tree
                                        .children(parent)
                                        .unwrap()
                                        .last()
                                        .map(|last| *last == n)
                                        .unwrap_or(false)
                                    {
                                        Position::LowerBranch
                                    } else {
                                        Position::MiddleBranch
                                    }
                                } else {
                                    Position::LowerBranch
                                }
                            },
                        })
                        .collect(),
                )
            })
            .collect();

        self.cache = Caches {
            narrowed_tree: Rc::new(tree),
            genes,
            lineages,
            tree_chars: Default::default(),
            duplications: Default::default(),
        };

        self.cache_tree_graph();
        self.cache_dup_nesting();
    }

    pub(crate) fn set_narrowing(&mut self, narrowing: ForthExpr) {
        self.narrowing = Some(narrowing);
        self.update_caches();
    }

    fn tree(&self) -> Rc<NewickTree> {
        self.cache.narrowed_tree.clone()
    }

    pub fn len(&self) -> usize {
        self.tree().leaves().count()
    }

    fn cache_tree_graph(&mut self) {
        self.cache.tree_chars = self
            .tree
            .leaves()
            .map(|n| (n, self.make_tree_line(n)))
            .collect();
    }

    fn cache_dup_nesting(&mut self) {
        self.cache.duplications.nestings.clear();
        for n in self.tree().leaves() {
            let mut pure_head_broken = false;
            let mut pure_tail_broken = false;
            let mut pure_head = ShiftRegister::new(3, false);
            let mut pure_tail = ShiftRegister::new(3, false);
            let dup_nesting = self.cache.lineages[&n]
                .iter()
                .filter_map(|n| {
                    if n.position != Position::UpperBranch {
                        pure_head_broken = true;
                    }
                    if n.position != Position::LowerBranch {
                        pure_tail_broken = true;
                    }
                    pure_tail.write(n.position == Position::LowerBranch && !pure_tail_broken);
                    pure_head.write(n.position == Position::UpperBranch && !pure_head_broken);

                    if self.tree().is_duplication(n.id) {
                        let dcs = self
                            .tree()
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

        let leaf_length = self.tree().topological_depth().1
            - self.tree().node_topological_depth(n).unwrap() as usize;
        let mut r = "─".repeat(leaf_length * DEPTH_FACTOR);

        let mut on_my_line = true;
        for n in lineage {
            let is_duplication = self
                .tree()
                .parent(n.id)
                .map(|x| self.tree().is_duplication(x))
                .unwrap_or(false);

            if on_my_line {
                match n.position {
                    Position::UpperBranch => {
                        if is_duplication {
                            r.push_str("─D");
                        } else {
                            r.push_str("─┬");
                        }
                    }
                    Position::MiddleBranch => {
                        on_my_line = false;
                        r.push_str(if is_duplication { "─╟" } else { "─├" })
                    }
                    Position::LowerBranch => {
                        on_my_line = false;
                        r.push_str(if is_duplication { "─╙" } else { "─└" })
                    }
                }
            } else {
                match n.position {
                    Position::UpperBranch | Position::MiddleBranch => {
                        if is_duplication {
                            r.push_str(" ║");
                        } else {
                            r.push_str(" │");
                        }
                    }
                    Position::LowerBranch => {
                        if is_duplication {
                            r.push_str("  ");
                        } else {
                            r.push_str("  ");
                        }
                    }
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
        fn percent_to_rgb((r, g, b): (f32, f32, f32)) -> Color {
            Color::Rgb(
                (r * 255.0).floor() as u8,
                (g * 255.0).floor() as u8,
                (b * 255.0).floor() as u8,
            )
        }

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
                            percent_to_rgb(
                                landscape_data
                                    .unwrap()
                                    .colors
                                    .get(&g.family)
                                    .unwrap()
                                    .to_percent(),
                            )
                        })
                    }))
                    .chain(
                        std::iter::once(
                            format!(" {} ", gene_to_char(*family, *strand, use_symbols)).fg({
                                percent_to_rgb(
                                    landscape_data
                                        .unwrap()
                                        .colors
                                        .get(family)
                                        .unwrap()
                                        .to_percent(),
                                )
                            }),
                        )
                        .chain(right_landscape.iter().map(|g| {
                            Span::from(format!(
                                " {}",
                                gene_to_char(g.family, g.strand, use_symbols)
                            ))
                            .fg({
                                percent_to_rgb(
                                    landscape_data
                                        .unwrap()
                                        .colors
                                        .get(&g.family)
                                        .unwrap()
                                        .to_percent(),
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
        for (y, n) in self.tree().leaves().enumerate() {
            let ancestors = self.cache.lineages[&n]
                .iter()
                .map(|n| n.id)
                .collect::<Vec<_>>();
            self.screen_to_nodes.insert(y, ancestors);
            let row = Self::gene_to_row(
                &self.cache.tree_chars[&n],
                self.landscape_data.as_ref(),
                self.cache.genes.get(&n).unwrap().clone(),
                &self.cache.duplications.nestings[&n],
                false,
                self.settings.use_symbols,
                &self.highlighters,
            );
            rows.push(row);
        }

        let tree_depth = self.tree().topological_depth().1;
        let widths = [
            Constraint::Length(1),
            Constraint::Length((DEPTH_FACTOR * (tree_depth + 1)) as u16),
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
        f.render_stateful_widget(&table, t, &mut self.states.gene_table);
    }

    /// Move the cursor to the given row in the table.
    pub fn move_to(&mut self, i: usize) {
        self.states.gene_table.select(Some(i));
        self.states.scrollbar = self.states.scrollbar.position(i);
    }

    /// Move the cursor one line up in the table.
    pub fn prev(&mut self, count: usize) {
        let i = self
            .states
            .gene_table
            .selected()
            .map(|i| i.saturating_sub(count))
            .unwrap_or_default();
        self.move_to(i);
    }

    /// Move the cursor one line down in the table.
    pub fn next(&mut self, count: usize) {
        let i = self
            .states
            .gene_table
            .selected()
            .map(|i| (i + count).clamp(0, self.len() - 1))
            .unwrap_or_default();
        self.move_to(i);
    }

    /// Move the cursor to the beginning of the table.
    pub fn top(&mut self) {
        self.move_to(0);
    }

    /// Move the cursor to the last line of the table.
    pub fn bottom(&mut self) {
        self.move_to(self.len() - 1);
    }

    /// Render the widget in the provided [`Rect`] within the [`Frame`].
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
