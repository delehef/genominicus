pub struct Canvas {
    rows: usize,
    columns: usize,
    frame: Vec<char>,
}
impl Canvas {
    pub fn new(rows: usize, columns: usize) -> Self {
        Canvas {
            rows,
            columns,
            frame: vec![' '; rows * columns],
        }
    }

    fn index(&self, row: usize, col: usize) -> usize {
        assert!(row < self.rows, "{} > {}", row, self.rows);
        assert!(col < self.columns, "{} >= {}", col, self.columns);

        row * self.columns + col
    }

    pub fn width(&self) -> usize {
        self.columns
    }

    pub fn height(&self) -> usize {
        self.rows
    }

    pub fn write(&mut self, row: usize, col: usize, c: char) {
        let idx = self.index(row, col);
        self.frame[idx] = c;
    }

    pub fn write_str(&mut self, row: usize, col: usize, s: &str) {
        for (i, c) in s.chars().enumerate() {
            self.write(row, col + i, c);
        }
    }

    pub fn line(&self, row: usize) -> String {
        assert!(row < self.rows);
        self.frame[row * self.columns..(row + 1) * self.columns]
            .iter()
            .collect()
    }
}
