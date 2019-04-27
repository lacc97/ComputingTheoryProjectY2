#ifndef SDL2_TEXTURE_HPP
#define SDL2_TEXTURE_HPP

#include <memory>

#include "sdl2_internal.hpp"

namespace std {
    template<>
    struct default_delete<SDL_Texture> {
        inline void operator()(SDL_Texture* ptr) {
            SDL_DestroyTexture(ptr);
        }
    };
}

namespace SDL2 {
    class Renderer;

    class LockTexture;

    class Texture : public Object<Texture> {
            friend class Object<Texture>;

            friend class Renderer;

            friend class LockTexture;

        public:
            Texture() : Texture(nullptr) {}

            inline Rect size() const {
                Rect r{0, 0, 0, 0};

                if(SDL_QueryTexture(pointer(), nullptr, nullptr, &r.w, &r.h) != 0)
                    throw Exception("Failed to query texture");

                return r;
            }

            inline SDL_BlendMode blendMode() const {
                SDL_BlendMode bm;

                if(SDL_GetTextureBlendMode(pointer(), &bm) != 0)
                    throw Exception("Failed to get renderer draw blend mode");

                return bm;
            }

            inline void setBlendMode(SDL_BlendMode bm) const {
                if(SDL_SetTextureBlendMode(pointer(), bm) != 0)
                    throw Exception("Failed to set renderer draw blend mode");
            }

        private:
            explicit Texture(SDL_Texture* tPtr) : mp_Ptr{tPtr} {}

            std::unique_ptr<SDL_Texture> mp_Ptr;
    };

    class LockTexture {

        public:
            explicit LockTexture(Texture& tRef) : mr_Texture{tRef}, mp_Pixels{nullptr}, m_Pitch{0} {
                if(SDL_LockTexture(mr_Texture.pointer(), nullptr, &mp_Pixels, &m_Pitch) != 0)
                    throw Exception("Failed to lock texture");
            }

            LockTexture(Texture& tRef, const SDL2::Rect& rect) : mr_Texture{tRef}, mp_Pixels{nullptr}, m_Pitch{0} {
                if(SDL_LockTexture(mr_Texture.pointer(), &rect, &mp_Pixels, &m_Pitch) != 0)
                    throw Exception("Failed to lock texture");
            }

            ~LockTexture() {
                SDL_UnlockTexture(mr_Texture.pointer());
            }

            inline void* pixels() const {
                return mp_Pixels;
            }

            inline int pitch() const {
                return m_Pitch;
            }

        private:
            Texture& mr_Texture;
            void* mp_Pixels;
            int m_Pitch;

    };
}

#endif
