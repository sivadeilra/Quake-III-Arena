use log::{Log, Metadata, Record};

pub struct DbgLogger;


impl  Log for DbgLogger {
    fn enabled(&self, _metadata: &Metadata) -> bool { true }
    fn log(&self, record: &Record) {
        let s = format!("{:8} {}] {}\r\n\0", record.level(),  record.module_path().unwrap_or("???"), record.args());
        unsafe {
            OutputDebugStringA(s.as_bytes().as_ptr());
        }
    }
    fn flush(&self) {

    }
}

pub static DBG_LOGGER:  DbgLogger = DbgLogger;


#[link(name = "kernel32")]
extern "stdcall" {
    fn OutputDebugStringA(s: *const  u8);
}
