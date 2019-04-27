#include "timer.hpp"

HRTimer::HRTimer(TimeUnit t, const std::string& name) : m_Unit{t}, m_Name{name}, m_Start{}, m_Elapsed{},
                                                        m_Started{false}, m_Paused{false} {
}

HRTimer::HRTimer(const std::string& name, TimeUnit t) : m_Unit{t}, m_Name{name}, m_Start{}, m_Elapsed{},
                                                        m_Started{false}, m_Paused{false} {
}

void HRTimer::start() {
    if(!m_Started || m_Paused) {
        m_Elapsed = duration_type::zero();
        m_Started = true;
        m_Paused = false;
        m_Start = std::chrono::high_resolution_clock::now();
    }
}

HRTimer::real_type HRTimer::pause() {
    if(!m_Paused && m_Started) {
        std::chrono::high_resolution_clock::time_point pauseTime = std::chrono::high_resolution_clock::now();
        m_Elapsed += std::chrono::duration_cast<duration_type>(pauseTime - m_Start);
        m_Paused = true;
    }

    return m_Elapsed.count() * m_Unit.unitsPerSecond;
}

HRTimer::real_type HRTimer::stop() {
    if(m_Started) {
        if(!m_Paused) {
            std::chrono::high_resolution_clock::time_point stopTime = std::chrono::high_resolution_clock::now();
            m_Elapsed += std::chrono::duration_cast<duration_type>(stopTime - m_Start);
        }

        m_Started = false;
    }

    return m_Elapsed.count() * m_Unit.unitsPerSecond;
}

void HRTimer::report(std::ostream& os) {
    pause();

    if(m_Name.empty())
        os << "[" << m_Elapsed.count() * m_Unit.unitsPerSecond << " " << m_Unit.unitString << "]" << std::endl;
    else
        os << m_Name << " [" << m_Elapsed.count() * m_Unit.unitsPerSecond << " " << m_Unit.unitString << "]"
           << std::endl;

    start();
}

void HRTimer::reset() {
    m_Elapsed = duration_type::zero();
    m_Started = false;
    m_Paused = false;
}
