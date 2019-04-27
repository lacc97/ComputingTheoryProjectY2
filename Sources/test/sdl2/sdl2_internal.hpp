#ifndef SDL2_INTERNAL_HPP
#define SDL2_INTERNAL_HPP

#include <type_traits>

#include <SDL2/SDL.h>

#ifdef ENABLE_SDL2_GFX

#   include <SDL2/SDL2_gfxPrimitives.h>

#endif
#ifdef ENABLE_SDL2_IMAGE
#   include <SDL2/SDL_image.h>
#endif
#ifdef ENABLE_SDL2_TTF
#   include <SDL2/SDL_ttf.h>
#endif
#ifdef ENABLE_SDL2_VULKAN

#   include <vulkan/vulkan.hpp>
#   include <SDL2/SDL_vulkan.h>

#endif

#include "exception.hpp"

namespace SDL2 {
    typedef SDL_Colour Colour;
    typedef SDL_Event Event;
    typedef SDL_Point Point;
    typedef SDL_Rect Rect;

    template<class E>
    class Object {
        public:
            inline void free() {
                static_cast<E*>(this)->mp_Ptr.reset();
            }

            inline bool operator!() const {
                return pointer() == nullptr;
            }

            inline operator bool() const {
                return pointer() != nullptr;
            }

//             inline operator internal_ptr_t() {
//                 return pointer();
//             }
//             inline operator const internal_ptr_t() const {
//                 return pointer();
//             }

        protected:
            inline auto pointer() const {
                return static_cast<const E*>(this)->mp_Ptr.get();
            }
    };
}

#endif
