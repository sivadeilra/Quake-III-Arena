#pragma once
#include <stdint.h>

// implemented in Rust
void __cdecl trace_string_raw(const char* s, size_t len);
void __cdecl trace_i32(int32_t i);
