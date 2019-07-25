#pragma once
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

// implemented in Rust
void __cdecl trace_str_raw(const char* s, size_t len);
void __cdecl trace_string_copy_raw(const char* s, size_t len);
void __cdecl trace_i32(int32_t i);
void __cdecl trace_f32(float f);

// For &'static str only, NOT for dynamic strings
inline void trace_str(const char* msg) {
    size_t len = strlen(msg);
    trace_str_raw(msg, len);
}

inline void trace_string(const char* msg) {
    size_t len = strlen(msg);
    trace_string_copy_raw(msg, len);
}

inline void trace_string_f(const char* fmt, ...) {
    char buffer[0x200];
    va_list va;

    va_start(va, fmt);
    vsprintf_s(buffer, sizeof(buffer), fmt, va);
    trace_string(buffer);
    va_end(va);
}
