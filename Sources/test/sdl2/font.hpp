#ifndef SDL2_FONT_HPP
#define SDL2_FONT_HPP

#include <memory>
#include <string>

#include "sdl2_internal.hpp"

#include "surface.hpp"

namespace std {
    template <> struct default_delete<TTF_Font> {
        inline void operator()(TTF_Font* ptr) {
            TTF_CloseFont(ptr);
        }
    };
}

namespace SDL2 {
    class Font : public Object<Font> {
            friend class Object<Font>;
        public:
            Font() : Font(nullptr) {}
            Font(const std::string& path, int pointSize = 12) : mp_Ptr{TTF_OpenFont(path.c_str(), pointSize)} {
                if(!mp_Ptr)
                    throw Exception("Failed to open font");
            }

            inline Surface renderText(const std::string& text, Colour c) {
                SDL_Surface* sPtr = TTF_RenderText_Blended(mp_Ptr.get(), text.c_str(), c);

                if(!sPtr)
                    throw Exception("Failed to render text");

                return Surface{sPtr};
            }

        private:
            Font(TTF_Font* fPtr) : mp_Ptr{fPtr} {}

            std::unique_ptr<TTF_Font> mp_Ptr;
    };
}

#endif
