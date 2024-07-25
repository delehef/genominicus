use anyhow::Context;
use newick::NewickTree;
use ratatui::{
    backend::Backend,
    crossterm::{
        event::{self, Event, KeyCode, KeyEvent, KeyEventKind},
        execute,
        terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
    },
    layout::{Constraint, Direction, Layout},
    prelude::*,
    style::Style,
    widgets::{Block, Borders, Paragraph},
    Frame, Terminal, TerminalOptions, Viewport,
};
use std::{io, rc::Rc};

use crate::utils::{ColorMap, GeneCache};

use self::widgets::treeview::{LandscapeData, TreeView, TreeViewSettings};

mod canvas;
mod forth;
mod utils;
pub(super) mod widgets;

enum Screen {
    TreeView,
}

struct States {
    highlighter: String,
}
impl States {
    fn new() -> Self {
        States {
            highlighter: String::new(),
        }
    }
}

#[derive(Copy, Clone)]
enum Mode {
    Root,
    Highlighter,
}
impl Mode {
    fn help(&self) -> Line {
        match self {
            Mode::Root => Line::from(vec![
                "[h]".yellow().bold(),
                "ighlight".into(),
                " :: ".bold().white(),
                "toggle ".into(),
                "[S]".yellow().bold(),
                "ymbols".into(),
                " :: ".bold().white(),
                "[TAB]".yellow().bold(),
                " cycle fold  ".into(),
                "←".yellow().bold(),
                " fold 1×  ".into(),
                "→".yellow().bold(),
                " unfold 1×  ".into(),
                "[q]".red().bold(),
                "uit ".into(),
            ]),
            Mode::Highlighter => Line::from(vec![
                "Highlights :: ".bold().white(),
                "[a]".yellow().bold(),
                "ppend ".into(),
                "[c]".yellow().bold(),
                "lear ".into(),
                "[p]".yellow().bold(),
                "op ".into(),
                "[e]".yellow().bold(),
                "dit last ".into(),
                "[q]".red().bold(),
                " back".into(),
            ]),
        }
    }
}

#[derive(Clone)]
pub struct Settings {
    pub tree: TreeViewSettings,
}

struct Editor {
    mode: Mode,
    name: String,
    tree: Rc<NewickTree>,
    plot: TreeView,
    screen: Screen,
    states: States,
    minibuffer: Rect,
}
impl Editor {
    pub fn new(
        name: String,
        tree: NewickTree,
        synteny: Option<(GeneCache, ColorMap)>,
        settings: Settings,
    ) -> Self {
        let landscape_data = if let Some((book, colors)) = synteny {
            Some(LandscapeData { book, colors })
        } else {
            None
        };
        let tree = Rc::new(tree);
        let plot = TreeView::from_newick(tree.clone(), settings.tree, landscape_data);

        Self {
            mode: Mode::Root,
            name,
            tree,
            plot,
            screen: Screen::TreeView,
            states: States::new(),
            minibuffer: Default::default(),
        }
    }

    fn render(&mut self, f: &mut Frame) {
        let chunks = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(3),
                Constraint::Min(1),
                Constraint::Length(3),
            ])
            .split(f.size());

        let title = Paragraph::new(self.mode.help()).block(
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
        self.plot.render(f, chunks[1]);
        self.minibuffer = chunks[2];
    }

    fn process_input(&mut self, key: KeyEvent) {
        match self.mode {
            Mode::Root => match key.code {
                KeyCode::Char('S') => {
                    self.plot.settings.use_symbols = !self.plot.settings.use_symbols;
                }
                KeyCode::Char('h') => self.mode = Mode::Highlighter,
                KeyCode::Up => self.plot.prev(1),
                KeyCode::Down => self.plot.next(1),
                KeyCode::PageUp => self.plot.prev(10),
                KeyCode::PageDown => self.plot.next(10),
                KeyCode::Home => self.plot.top(),
                KeyCode::End => self.plot.bottom(),
                KeyCode::Left => self.plot.fold_current(),
                KeyCode::Right => self.plot.unfold_current(),
                KeyCode::Tab => self.plot.toggle_current(),
                _ => {}
            },
            Mode::Highlighter => {
                match key.code {
                    KeyCode::Char('a') => {
                        let mut t = Terminal::with_options(
                            CrosstermBackend::new(std::io::stdout()),
                            TerminalOptions {
                                viewport: Viewport::Fixed(self.minibuffer),
                            },
                        )
                        .unwrap();
                        let expr = widgets::scan::ScanInput::new(String::new())
                            .run(&mut t, self.minibuffer);
                        if let Some((source, expr)) = expr {
                            self.states.highlighter = source;
                            self.plot.highlighters.push(expr);
                        }
                    }
                    KeyCode::Char('c') => self.plot.highlighters.clear(),
                    KeyCode::Char('p') => {
                        self.plot.highlighters.pop();
                    }
                    KeyCode::Char('e') => {
                        if self.plot.highlighters.pop().is_some() {
                            let mut t = Terminal::with_options(
                                CrosstermBackend::new(std::io::stdout()),
                                TerminalOptions {
                                    viewport: Viewport::Fixed(self.minibuffer),
                                },
                            )
                            .unwrap();

                            if let Some((source, expr)) =
                                widgets::scan::ScanInput::new(self.states.highlighter.clone())
                                    .run(&mut t, self.minibuffer)
                            {
                                self.states.highlighter = source;
                                self.plot.highlighters.push(expr);
                            }
                        }
                    }
                    _ => {}
                };
                self.mode = Mode::Root
            }
        }
    }

    fn run<B: Backend>(&mut self, terminal: &mut Terminal<B>) -> anyhow::Result<()> {
        loop {
            terminal.draw(|term| self.render(term))?;
            if let Event::Key(key) = event::read()? {
                if key.kind == KeyEventKind::Press {
                    if let KeyCode::Char('q') = key.code {
                        match self.mode {
                            Mode::Root => return Ok(()),
                            _ => self.mode = Mode::Root,
                        }
                    }
                    self.process_input(key);
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
