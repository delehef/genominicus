use crate::editor::forth::{self, Node};
use ratatui::{
    backend::Backend,
    crossterm,
    prelude::Rect,
    style::{Color, Style},
    widgets::{Block, Borders},
    Terminal,
};
use tui_textarea::{CursorMove, Input, Key, TextArea};

pub struct ScanInput<'a> {
    input: TextArea<'a>,
}
impl<'a> ScanInput<'a> {
    pub fn new(content: String) -> Self {
        let mut r = ScanInput {
            input: TextArea::from([content]),
        };
        r.input.move_cursor(CursorMove::End);
        r
    }

    fn validate(&mut self) -> anyhow::Result<Node> {
        let r = forth::parse(&self.input.lines()[0]);

        match &r {
            Err(err) => {
                self.input.set_style(Style::default().fg(Color::LightRed));
                self.input.set_block(
                    Block::default()
                        .borders(Borders::ALL)
                        .title(format!("ERROR: {}", err)),
                );
            }
            Ok(node) => {
                self.input.set_style(Style::default().fg(Color::LightGreen));
                self.input.set_block(
                    Block::default()
                        .borders(Borders::ALL)
                        .title(node.to_string()),
                );
            }
        };
        r
    }

    pub fn run<B: Backend>(
        mut self,
        term: &mut Terminal<B>,
        target: Rect,
    ) -> Option<(String, Node)> {
        self.input.set_cursor_line_style(Style::default());
        loop {
            let _ = self.validate();
            let _ = term.draw(|f| {
                f.render_widget(self.input.widget(), target);
            });

            match crossterm::event::read().unwrap().into() {
                Input {
                    key: Key::Enter, ..
                } => {
                    let _ = term.clear();
                    return self
                        .validate()
                        .ok()
                        .map(|i| (self.input.into_lines()[0].to_owned(), i));
                }
                Input { key: Key::Esc, .. } => {
                    let _ = term.clear();
                    return None;
                }
                input => {
                    self.input.input(input);
                }
            }
        }
    }
}
