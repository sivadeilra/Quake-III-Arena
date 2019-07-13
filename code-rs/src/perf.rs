pub struct StaticCounter {
    name: &'static str,
}

impl StaticCounter {
    pub const fn new(name: &'static str) -> StaticCounter {
        StaticCounter { name }
    }

    pub fn add_one(&self) {}
    pub fn add_usize(&self, _n: usize) {}
    pub fn add_u64(&self, _n: u64) {}
}
