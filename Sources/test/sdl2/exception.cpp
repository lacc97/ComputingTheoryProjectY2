#include "exception.hpp"

#include <sstream>

#include "sdl2_internal.hpp"

SDL2::Exception::Exception(const std::string& msg, bool getSDlError) noexcept : m_What{} {
    std::stringstream ss;

    if(msg.empty()) {
        if(!getSDlError)
            ss << "SDL error";
        else
            ss << "SDL error: " << SDL_GetError();
    } else {
        ss << msg;

        if(getSDlError)
            ss << "(SDL error: " << SDL_GetError() << ")";
    }

    m_What = ss.str();
}

const char* SDL2::Exception::what() const noexcept {
    return m_What.c_str();
}

SDL2::VKException::VKException(const std::string& msg, bool getSDlError) noexcept {
    std::stringstream ss;

    ss << "Vulkan error";

    if(!msg.empty())
        ss << ": " << msg;

    if(getSDlError)
        ss << " (SDL error: " << SDL_GetError() << ")";

    m_What = ss.str();
}
