#include "types.hpp"

#include <cstdio>
#include <cstdarg>

#if __cplusplus < 201703L

#include <memory>

#endif

#ifndef WITH_LOW_PRECISION

// printing functions for __int128 and unsigned __int128
std::ostream& operator<<(std::ostream& dest, __int128_t value) {
    typedef std::ostream::char_type char_type;

    std::ostream::sentry s(dest);

    if(s) {
        auto tmp = __uint128_t(value < 0 ? -value : value);
        char_type buffer[128];
        char_type* d = std::end(buffer);

        do {
            --d;
            *d = "0123456789"[tmp % 10];
            tmp /= 10;
        } while(tmp != 0);

        if(value < 0) {
            --d;
            *d = '-';
        }

        size_t len = std::end(buffer) - d;

        if(dest.rdbuf()->sputn(d, len) != len) {
            dest.setstate(std::ios_base::badbit);
        }
    }

    return dest;
}

std::ostream& operator<<(std::ostream& dest, __uint128_t value) {
    typedef std::ostream::char_type char_type;

    std::ostream::sentry s(dest);

    if(s) {
        char_type buffer[128];
        char_type* d = std::end(buffer);

        do {
            --d;
            *d = "0123456789"[value % 10];
            value /= 10;
        } while(value != 0);

        size_t len = std::end(buffer) - d;

        if(dest.rdbuf()->sputn(d, len) != len) {
            dest.setstate(std::ios_base::badbit);
        }
    }

    return dest;
}

#endif

// Implementation (https://codereview.stackexchange.com/questions/187183/create-a-c-string-using-printf-style-formatting)
std::string format(const char* fmt, ...) {
    char buf[256];

    va_list args;
    va_start(args, fmt);
    const auto r = std::vsnprintf(buf, sizeof buf, fmt, args);
    va_end(args);

    if(r < 0)
        // conversion failed
        return {};

    auto len = size_t(r);

    if(len < sizeof buf)
        // we fit in the buffer
        return {buf, len};

#if __cplusplus >= 201703L
    // C++17: Create a string and write to its underlying array
    std::string s(len, '\0');

    va_start(args, fmt);

    std::vsnprintf(s.data(), len + 1, fmt, args);

    va_end(args);

    return s;

#else
    // C++11 or C++14: We need to allocate scratch memory
    auto vbuf = std::unique_ptr<char[]>(new char[len + 1]);

    va_start(args, fmt);

    std::vsnprintf(vbuf.get(), len + 1, fmt, args);

    va_end(args);

    return {vbuf.get(), len};

#endif
}
