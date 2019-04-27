#ifndef SDL2_WINDOW_HPP
#define SDL2_WINDOW_HPP

#include <string>
#include <memory>

#ifdef ENABLE_SDL2_VULKAN

#   include <vulkan/vulkan.hpp>

#endif

#include "sdl2_internal.hpp"

#include "renderer.hpp"

namespace std {
    template<>
    struct default_delete<SDL_Window> {
        inline void operator()(SDL_Window* ptr) {
            SDL_DestroyWindow(ptr);
        }
    };
}

namespace SDL2 {
    class Window : public Object<Window> {
            friend class Object<Window>;

        public:
            Window() : mp_Ptr{nullptr}, m_Flags{0} {}

            Window(const std::string& title, int w, int h, uint32_t flags = SDL_WINDOW_SHOWN | SDL_WINDOW_ALLOW_HIGHDPI)
                    : Window(title, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, w, h, flags) {}

            Window(const std::string& title, int x, int y, int w, int h,
                   uint32_t flags = SDL_WINDOW_SHOWN | SDL_WINDOW_ALLOW_HIGHDPI) : mp_Ptr{
                    SDL_CreateWindow(title.c_str(), x, y, w, h, flags)}, m_Flags{flags} {
                if(!mp_Ptr)
                    throw Exception("Failed to open window");
            }

            inline Renderer createRenderer(uint32_t flags = SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC) {
                SDL_Renderer* rend = SDL_CreateRenderer(mp_Ptr.get(), -1, flags);

                if(!rend)
                    throw Exception("Failed to create renderer");

                return Renderer{rend};
            }

#ifdef ENABLE_SDL2_VULKAN

            inline vk::SurfaceKHR createSurface(const vk::Instance& inst) const {
                VkSurfaceKHR surf;

                if(SDL_Vulkan_CreateSurface(mp_Ptr.get(), inst, &surf) != SDL_TRUE)
                    throw VKException("Failed to create surface");

                return surf;
            }

#endif
#ifdef ENABLE_SDL2_OPENGL
            inline void swap() {
                SDL_GL_SwapWindow(mp_Ptr.get());
            }
#endif

            inline int drawableHeight() const {
                return drawableSize().h;
            }

            inline int drawableWidth() const {
                return drawableSize().w;
            }

            inline Rect drawableSize() const {
                Rect dim{0, 0, 0, 0};

                if((m_Flags & SDL_WINDOW_OPENGL) != 0) {
                    SDL_GL_GetDrawableSize(mp_Ptr.get(), &dim.w, &dim.h);
                }

#ifdef ENABLE_SDL2_VULKAN
                else if((m_Flags & SDL_WINDOW_VULKAN) != 0) {
                    SDL_Vulkan_GetDrawableSize(mp_Ptr.get(), &dim.w, &dim.h);
                }

#endif
                else {
                    SDL_Renderer* ren = SDL_GetRenderer(mp_Ptr.get());

                    if(ren == nullptr)
                        throw Exception("Failed to get rendering context");

                    SDL_GetRendererOutputSize(ren, &dim.w, &dim.h);
                }

                return dim;
            }

            inline int height() const {
                return size().h;
            }

#ifdef ENABLE_SDL2_VULKAN

            inline std::vector<const char*> instanceExtensions() const {
                unsigned int count;

                if(SDL_Vulkan_GetInstanceExtensions(mp_Ptr.get(), &count, nullptr) != SDL_TRUE)
                    throw VKException("Failed to get instance extensions");

                std::vector<const char*> exts;
                exts.resize(count);

                if(SDL_Vulkan_GetInstanceExtensions(mp_Ptr.get(), &count, exts.data()) != SDL_TRUE)
                    throw VKException("Failed to get instance extensions");

                return exts;
            }

#endif

            inline int width() const {
                return size().w;
            }

            inline Rect size() const {
                Rect dim{0, 0, 0, 0};
                SDL_GetWindowSize(mp_Ptr.get(), &dim.w, &dim.h);
                return dim;
            }

        private:
            std::unique_ptr<SDL_Window> mp_Ptr;
            uint32_t m_Flags;
    };
}

#endif
