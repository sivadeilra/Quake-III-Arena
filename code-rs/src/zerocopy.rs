pub unsafe trait FromBytes {}

pub unsafe trait ToBytes {}

macro_rules! to_from_bytes {
    ($($name:ident),*) => {
        $(
        unsafe impl ToBytes for $name {}
        unsafe impl FromBytes for $name {}
        )*
    }
}

to_from_bytes! {
    u8, u16, u32, u64, usize,
    i8, i16, i32, i64, isize
}

pub fn from_bytes<T: FromBytes>(bytes: &[u8]) -> &T {
    assert!(core::mem::size_of::<T>() > 0);
    unsafe {
        assert!(bytes.len() >= core::mem::size_of::<T>());
        let ptr = bytes.as_ptr();
        assert!(ptr as usize % core::mem::align_of::<T>() == 0);
        &*(ptr as *const T)
    }
}

pub fn slice_from_bytes<T: FromBytes>(bytes: &[u8]) -> &[T] {
    assert!(core::mem::size_of::<T>() > 0);
    unsafe {
        let slice_len = bytes.len() / core::mem::size_of::<T>();
        if slice_len != 0 {
            let ptr = bytes.as_ptr();
            core::slice::from_raw_parts(ptr as *const T, slice_len)
        } else {
            &[]
        }
    }
}

pub fn to_bytes<T: ToBytes>(t: &T) -> &[u8] {
    unsafe { to_bytes_unsafe(t) }
}

pub unsafe fn to_bytes_unsafe<T>(t: &T) -> &[u8] {
    core::slice::from_raw_parts::<u8>(t as *const T as *const u8, core::mem::size_of::<T>())
}

pub fn slice_to_bytes<T: ToBytes>(s: &[T]) -> &[u8] {
    unsafe { slice_to_bytes_unsafe(s) }
}

pub unsafe fn slice_to_bytes_unsafe<T>(s: &[T]) -> &[u8] {
    core::slice::from_raw_parts::<u8>(s.as_ptr() as *const u8, s.len() * core::mem::size_of::<T>())
}
