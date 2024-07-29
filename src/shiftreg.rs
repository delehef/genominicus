pub struct ShiftRegister<T> {
    reg: Vec<T>,
    i: usize,
}

impl<T: Default> ShiftRegister<T> {
    pub fn default(shift: usize) -> Self {
        assert!(shift > 0);
        Self {
            reg: std::iter::repeat_with(|| T::default())
                .take(shift)
                .collect(),
            i: 0,
        }
    }
}
impl<T: Clone> ShiftRegister<T> {
    pub fn new(shift: usize, init: T) -> Self {
        assert!(shift > 0);
        Self {
            reg: vec![init; shift],
            i: 0,
        }
    }
}
impl<T> ShiftRegister<T> {
    pub fn read(&self) -> &T {
        &self.reg[self.i]
    }

    pub fn write(&mut self, t: T) {
        self.i += 1;
        self.i %= self.reg.len();
        let new_idx = (self.i + (self.reg.len() - 1)) % self.reg.len();
        self.reg[new_idx] = t;
    }
}
