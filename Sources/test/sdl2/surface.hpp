#ifndef SDL2_SURFACE_HPP
#define SDL2_SURFACE_HPP

#include <memory>

#include "sdl2_internal.hpp"

namespace std {
    template <> struct default_delete<SDL_Surface> {
        inline void operator()(SDL_Surface* ptr) {
            SDL_FreeSurface(ptr);
        }
    };
}

namespace SDL2 {
    class Font;
    class Renderer;

    class Surface : public Object<Surface> {
            friend class Font;
            friend class Renderer;
#ifdef ENABLE_SDL2_IMAGE
            friend Surface loadImage(const std::string&);
#endif

        public:
            Surface() : Surface(nullptr) {}

        private:
            Surface(SDL_Surface* ptr) : mp_Ptr(ptr) {}

            std::unique_ptr<SDL_Surface> mp_Ptr;
    };
}

#endif
