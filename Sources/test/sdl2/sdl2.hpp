#ifndef SDL2_HPP
#define SDL2_HPP

#include <cstdint>

#include "sdl2_internal.hpp"
#include "exception.hpp"

namespace SDL2 {
    inline void delay(uint32_t msec) {
        SDL_Delay(msec);
    }

    inline void init(uint32_t flags = SDL_INIT_EVERYTHING) {
        if(SDL_Init(flags) != 0)
            throw Exception("Failed to initialise SDL2");
    }

    inline void quit() {
        SDL_Quit();
    }

#ifdef ENABLE_SDL2_TTF
    inline void initTTF() {
        if(TTF_Init() != 0)
            throw Exception("Failed to initialise SDL2_ttf");
    }
    inline void quitTTF() {
        TTF_Quit();
    }
#endif

    inline const uint8_t* getKeyboardState() {
        return SDL_GetKeyboardState(nullptr);
    }

    inline uint32_t getTicks() {
        return SDL_GetTicks();
    }

    inline bool poll(SDL2::Event& ev) {
        return SDL_PollEvent(&ev);
    }

#ifdef ENABLE_SDL2_VULKAN

    PFN_vkGetInstanceProcAddr getVkGetInstanceProcAddr() {
        return reinterpret_cast<PFN_vkGetInstanceProcAddr>(SDL_Vulkan_GetVkGetInstanceProcAddr());
    }

#endif
}

#ifdef ENABLE_SDL2_IMAGE
#   include "sdl2_image.hpp"
#endif

#ifdef ENABLE_SDL2_TTF
#   include "font.hpp"
#endif

#include "renderer.hpp"
#include "surface.hpp"
#include "texture.hpp"
#include "window.hpp"

#endif
