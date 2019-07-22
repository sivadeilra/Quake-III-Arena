pub struct cvar_t {
    pub name: &'static str,
}

impl cvar_t {
    pub const fn new(name: &'static str) -> Self {
        Self { name }
    }

    pub fn as_i32(&self) -> i32 {
        0
    }

    pub fn as_bool(&self) -> bool {
        false
    }
}
