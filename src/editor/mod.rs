use anyhow::Context;
use newick::{Newick, NewickTree, NodeID};
use ratatui::{
    backend::Backend,
    crossterm::{
        event::{self, Event, KeyCode, KeyEventKind},
        execute,
        terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
    },
    layout::{Constraint, Direction, Layout},
    prelude::*,
    style::Style,
    widgets::{
        Block, Borders, Cell, Paragraph, Row, Scrollbar, ScrollbarOrientation, ScrollbarState,
        Table, TableState,
    },
    Frame, Terminal, TerminalOptions, Viewport,
};
use std::{collections::HashMap, io, iter::FromIterator, ops::Range, sync::OnceLock};
use syntesuite::genebook::Gene;

use crate::{
    name2color,
    utils::{ColorMap, GeneCache, WINDOW},
};

use self::{canvas::Canvas, forth::Node};

mod canvas;
mod forth;
mod utils;
mod widgets;

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
pub struct TreeSettings {
    pub use_symbols: bool,
}

#[derive(Clone)]
pub struct Settings {
    pub tree: TreeSettings,
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

enum Screen {
    TreeView,
}

#[derive(Debug, Clone)]
struct DispGene {
    name: String,
    species: String,
}

struct LandscapeData {
    book: GeneCache,
    colors: ColorMap,
}

struct States {
    gene_table: TableState,
    scrollbar: ScrollbarState,
    highlighter: String,
}
impl States {
    fn new(size: usize) -> Self {
        States {
            gene_table: TableState::new().with_selected(0),
            scrollbar: ScrollbarState::new(size - 1),
            highlighter: String::new(),
        }
    }
}

#[derive(Debug)]
enum Clade {
    Taxon {
        graph_line: usize,
        dup_nesting: Vec<f32>,
        gene: DispGene,
    },
    SubClade {
        subclades: Vec<usize>,
        folded: bool,
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

struct CanvassedTree {
    settings: TreeSettings,
    canvas: Canvas,
    dup_level: Vec<Vec<f32>>,
    current_len: usize,
    // screen coordinate -> clade ID
    screen_to_clade: HashMap<usize, Vec<usize>>,
    clades: CladeHierarchy,
    highlight: Option<Node>,
}
impl CanvassedTree {
    fn rec_make_tree(
        t: &NewickTree,
        n: NodeID,
        current_y: usize,
        current_depth: usize,
        graph: &mut Canvas,
        mut dups_nesting: Vec<f32>,
        dups: &mut Vec<Vec<f32>>,
        clades: &mut CladeHierarchy,
        current_clade: usize,
    ) -> usize {
        if t[n].is_leaf() {
            dups[current_y] = dups_nesting.clone();
            let my_clade = Clade::Taxon {
                graph_line: current_y,
                dup_nesting: dups_nesting,
                gene: DispGene {
                    name: t.name(n).cloned().unwrap_or("UNKNWN".into()),
                    species: t.attrs(n).get("S").cloned().unwrap_or("UNKNWN".into()),
                },
            };
            clades.append_in(my_clade, current_clade);
            current_y + 1
        } else {
            let my_clade = clades.append_in(
                Clade::SubClade {
                    subclades: vec![],
                    folded: false,
                },
                current_clade,
            );
            if t.is_duplication(n) {
                graph.write_str(current_y, current_depth, "─D");
                dups_nesting.push(
                    t.attrs(n)
                        .get("DCS")
                        .map(|x| x.parse::<f32>().unwrap())
                        .unwrap_or(-1.),
                );
            } else {
                graph.write_str(current_y, current_depth, "─┬");
            };

            let old_y = current_y;
            let mut current_y = current_y;
            let bar = if t.is_duplication(n) { '║' } else { '│' };
            let l = if t.is_duplication(n) { '╙' } else { '└' };

            for (i, child) in t[n].children().iter().enumerate() {
                if t[*child].is_leaf() {
                    for i in current_depth + 2..graph.width() {
                        graph.write(current_y, i, '─');
                    }
                }
                for y in old_y + 1..current_y {
                    graph.write(y, current_depth + 1, bar);
                }
                if i > 0 {
                    graph.write(current_y, current_depth + 1, l);
                }

                current_y = Self::rec_make_tree(
                    t,
                    *child,
                    current_y,
                    current_depth + DEPTH_FACTOR,
                    graph,
                    dups_nesting.clone(),
                    dups,
                    clades,
                    my_clade,
                );
            }
            current_y
        }
    }

    fn from_newick(t: &NewickTree, settings: TreeSettings) -> Self {
        let leave_count = t.leaves().count();
        let mut canvas = Canvas::new(leave_count, DEPTH_FACTOR * t.topological_depth().1);
        let mut dups = vec![vec![]; leave_count];
        let mut clades: CladeHierarchy = CladeHierarchy::new();

        Self::rec_make_tree(
            t,
            t.root(),
            0,
            0,
            &mut canvas,
            Vec::new(),
            &mut dups,
            &mut clades,
            0,
        );

        Self {
            settings,
            canvas,
            dup_level: dups,
            current_len: leave_count,
            screen_to_clade: Default::default(),
            clades,
            highlight: None,
        }
    }

    fn len(&self) -> usize {
        self.current_len
    }

    fn gene_to_row<'a>(
        graph_line: String,
        landscape_data: Option<&'a LandscapeData>,
        gene: DispGene,
        dups_nesting: Vec<f32>,
        with_fold_indicator: bool,
        use_symbols: bool,
        highlighter: Option<&'a Node>,
    ) -> Row<'a> {
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
        let highlighted = highlighter
            .map(|n| n.eval(&gene).unwrap().right().unwrap())
            .unwrap_or(false);
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
                    .map(|x| {
                        Span::from("│").fg(if *x < 0. {
                            Color::LightBlue
                        } else {
                            Color::Rgb(((1. - x) * 255.0).ceil() as u8, (x * 255.).ceil() as u8, 0)
                        })
                    })
                    .collect::<Vec<_>>(),
            )),
            gene.species
                .fg(Color::Rgb(
                    (species_color.0 * 255.0).floor() as u8,
                    (species_color.1 * 255.0).floor() as u8,
                    (species_color.2 * 255.0).floor() as u8,
                ))
                .into(),
            if highlighted {
                gene.name.clone().bold().blue().reversed().into()
            } else {
                gene.name.clone().into()
            },
            landscape.into(),
        ])
    }

    fn to_rows<'a>(&'a mut self, landscape_data: Option<&'a LandscapeData>) -> Vec<Row> {
        #[derive(Debug)]
        struct RowContext {
            current: usize,
            in_folded: bool,
        }

        let mut rows = Vec::new();
        self.screen_to_clade.clear();

        let mut todos = Vec::new();
        todos.push(RowContext {
            current: 0,
            in_folded: false,
        });

        while let Some(context) = todos.pop() {
            let current = self.clades.get(context.current);
            match current {
                Clade::Taxon {
                    dup_nesting,
                    gene,
                    graph_line,
                } => {
                    self.screen_to_clade.entry(rows.len()).or_default();
                    let row = Self::gene_to_row(
                        self.canvas.line(*graph_line),
                        landscape_data,
                        gene.to_owned(),
                        dup_nesting.clone(),
                        false,
                        self.settings.use_symbols,
                        self.highlight.as_ref(),
                    );
                    rows.push(row);
                }
                Clade::SubClade { subclades, folded } => {
                    self.screen_to_clade
                        .entry(rows.len())
                        .or_default()
                        .push(context.current);
                    if *folded {
                        if let Clade::Taxon {
                            graph_line,
                            dup_nesting,
                            gene,
                        } = self.clades.find_first_taxon(context.current)
                        {
                            let row = Self::gene_to_row(
                                self.canvas.line(*graph_line),
                                landscape_data,
                                gene.to_owned(),
                                dup_nesting.clone(),
                                true,
                                self.settings.use_symbols,
                                self.highlight.as_ref(),
                            );
                            rows.push(row);
                        }
                    } else {
                        for subclade in subclades.iter().rev() {
                            todos.push(RowContext {
                                current: *subclade,
                                in_folded: context.in_folded,
                            });
                        }
                    }
                }
            }
        }
        self.current_len = rows.len();
        rows
    }

    fn toggle_all_at(&mut self, screen_y: usize) {
        let target_state = !self.screen_to_clade[&screen_y]
            .iter()
            .any(|c| self.clades.is_folded(*c));

        for clade in self.screen_to_clade.get(&screen_y).unwrap().iter().rev() {
            if let Clade::SubClade { ref mut folded, .. } = self.clades.get_mut(*clade) {
                *folded = target_state;
            } else {
                unreachable!()
            }
        }
    }

    fn fold_at(&mut self, screen_y: usize) {
        for clade in self.screen_to_clade.get(&screen_y).unwrap().iter().rev() {
            if let Clade::SubClade { ref mut folded, .. } = self.clades.get_mut(*clade) {
                if !*folded {
                    *folded = true;
                    return;
                }
            } else {
                unreachable!()
            }
        }
    }

    fn unfold_at(&mut self, screen_y: usize) {
        for clade in self.screen_to_clade.get(&screen_y).unwrap().iter() {
            if let Clade::SubClade { ref mut folded, .. } = self.clades.get_mut(*clade) {
                if *folded {
                    *folded = false;
                    return;
                }
            } else {
                unreachable!()
            }
        }
    }

    fn max_dup_nesting(&self) -> u16 {
        self.dup_level.iter().map(Vec::len).max().unwrap_or(0) as u16
    }
}

struct Editor {
    name: String,
    tree: NewickTree,
    plot: CanvassedTree,
    screen: Screen,
    states: States,
    landscape_data: Option<LandscapeData>,
    minibuffer: Rect,
}
impl Editor {
    pub fn new(
        name: String,
        tree: NewickTree,
        synteny: Option<(GeneCache, ColorMap)>,
        settings: Settings,
    ) -> Self {
        let plot = CanvassedTree::from_newick(&tree, settings.tree);
        let states = States::new(plot.len());
        let landscape_data = if let Some((book, colors)) = synteny {
            Some(LandscapeData { book, colors })
        } else {
            None
        };

        Self {
            name,
            tree,
            plot,
            screen: Screen::TreeView,
            states,
            landscape_data,
            minibuffer: Default::default(),
        }
    }

    fn move_to(&mut self, i: usize) {
        self.states.gene_table.select(Some(i));
        self.states.scrollbar = self.states.scrollbar.position(i);
    }

    fn prev(&mut self, count: usize) {
        let i = self
            .states
            .gene_table
            .selected()
            .map(|i| i.saturating_sub(count))
            .unwrap_or_default();
        self.move_to(i);
    }

    fn next(&mut self, count: usize) {
        let i = self
            .states
            .gene_table
            .selected()
            .map(|i| (i + count).clamp(0, self.plot.len() - 1))
            .unwrap_or_default();
        self.move_to(i);
    }

    fn top(&mut self) {
        self.move_to(0);
    }

    fn bottom(&mut self) {
        self.move_to(self.plot.len() - 1);
    }

    fn render(&mut self, f: &mut Frame) {
        const INDENT_FACTOR: usize = 2;
        let chunks = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(3),
                Constraint::Min(1),
                Constraint::Length(3),
            ])
            .split(f.size());

        let title = Paragraph::new(Line::from(vec![
            "h".yellow().bold(),
            "ighlight".into(),
            " :: ".bold().white(),
            "toggle ".into(),
            "S".yellow().bold(),
            "ymbols".into(),
            " :: ".bold().white(),
            "[TAB]".yellow().bold(),
            " cycle fold  ".into(),
            "←".yellow().bold(),
            " fold 1×  ".into(),
            "→".yellow().bold(),
            " unfold 1×  ".into(),
        ]))
        .block(
            Block::default()
                .borders(Borders::BOTTOM)
                .style(Style::default())
                .title(Line::from(vec![
                    self.name.as_str().bold(),
                    ": ".into(),
                    self.tree.leaves().count().to_string().into(),
                    " genes".into(),
                ])),
        );
        f.render_widget(title, chunks[0]);

        self.minibuffer = chunks[2];

        let tree_depth = self.tree.topological_depth().1;
        let dups_width = self.plot.max_dup_nesting();
        let rows = self.plot.to_rows(self.landscape_data.as_ref());
        let widths = [
            Constraint::Length(1),
            Constraint::Length((INDENT_FACTOR * tree_depth) as u16),
            Constraint::Length(dups_width),
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

        f.render_stateful_widget(table, chunks[1], &mut self.states.gene_table);

        f.render_stateful_widget(
            Scrollbar::default()
                .orientation(ScrollbarOrientation::VerticalRight)
                .begin_symbol(Some("^"))
                .end_symbol(Some("v")),
            chunks[1].inner(Margin {
                vertical: 1,
                horizontal: 1,
            }),
            &mut self.states.scrollbar,
        );
    }

    fn run<B: Backend>(&mut self, terminal: &mut Terminal<B>) -> anyhow::Result<()> {
        GENABET
            .set(
                BLOCKS
                    .iter()
                    .flat_map(|b| b.clone())
                    .map(|c| char::from_u32(c).unwrap())
                    .collect(),
            )
            .unwrap();

        loop {
            terminal.draw(|term| self.render(term))?;
            if let Event::Key(key) = event::read()? {
                if key.kind == KeyEventKind::Press {
                    match key.code {
                        KeyCode::Char('q') => return Ok(()),
                        KeyCode::Char('S') => {
                            self.plot.settings.use_symbols = !self.plot.settings.use_symbols;
                        }
                        KeyCode::Char('h') => {
                            let mut t = Terminal::with_options(
                                CrosstermBackend::new(std::io::stdout()),
                                TerminalOptions {
                                    viewport: Viewport::Fixed(self.minibuffer),
                                },
                            )
                            .unwrap();
                            let expr = widgets::scan::ScanInput::new(&self.states.highlighter)
                                .run(&mut t, self.minibuffer);
                            if let Some((source, expr)) = expr {
                                self.states.highlighter = source;
                                self.plot.highlight = Some(expr);
                            }
                        }
                        KeyCode::Up => self.prev(1),
                        KeyCode::Down => self.next(1),
                        KeyCode::PageUp => self.prev(10),
                        KeyCode::PageDown => self.next(10),
                        KeyCode::Home => self.top(),
                        KeyCode::End => self.bottom(),
                        KeyCode::Left => {
                            self.plot
                                .fold_at(self.states.gene_table.selected().unwrap());
                        }
                        KeyCode::Right => {
                            self.plot
                                .unfold_at(self.states.gene_table.selected().unwrap());
                        }
                        KeyCode::Tab => self
                            .plot
                            .toggle_all_at(self.states.gene_table.selected().unwrap()),
                        _ => {}
                    }
                    terminal.draw(|term| self.render(term))?;
                }
            }
        }
    }
}

pub fn run(
    name: String,
    t: NewickTree,
    synteny: Option<(GeneCache, ColorMap)>,
    settings: Settings,
) -> anyhow::Result<()> {
    enable_raw_mode()?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen).context("failed to initialize terminal")?;
    let backend = ratatui::backend::CrosstermBackend::new(stdout);
    let mut terminal = ratatui::Terminal::new(backend)?;

    let mut editor = Editor::new(name, t, synteny, settings);
    editor.run(&mut terminal)?;

    disable_raw_mode()?;
    execute!(terminal.backend_mut(), LeaveAlternateScreen,)?;
    terminal.show_cursor()?;

    Ok(())
}
