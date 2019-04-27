#pragma once

#include <iostream>
#include <chrono>
#include <string>

typedef struct {
    long double unitsPerSecond;
    const char* unitString;
} TimeUnit;

constexpr TimeUnit HOURS{1.0 / 3600.0, "h"};
constexpr TimeUnit MINUTES{1.0 / 60.0, "m"};
constexpr TimeUnit SECONDS{1, "s"};
constexpr TimeUnit MILLISECONDS{1000, "ms"};
constexpr TimeUnit MICROSECONDS{1000000, "μs"};
constexpr TimeUnit NANOSECONDS{1000000000L, "ns"}; // likely inaccurate?

class HRTimer {
        typedef long double real_type;
        typedef std::chrono::duration<real_type> duration_type;

    public:
        explicit HRTimer(TimeUnit t = MILLISECONDS, const std::string& name = "");

        explicit HRTimer(const std::string& name, TimeUnit t = MILLISECONDS);

        void start();

        real_type pause();

        real_type stop();

        void reset();

        void report(std::ostream& os = std::cout);

    private:
        TimeUnit m_Unit;
        std::string m_Name;
        std::chrono::high_resolution_clock::time_point m_Start;
        //std::chrono::high_resolution_clock::time_point m_Pause;
        //std::chrono::high_resolution_clock::time_point m_Stop;
        duration_type m_Elapsed;
        bool m_Started;
        bool m_Paused;
};

