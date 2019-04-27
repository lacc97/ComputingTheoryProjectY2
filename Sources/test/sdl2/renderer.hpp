#ifndef SDL2_RENDERER_HPP
#define SDL2_RENDERER_HPP

#include <memory>
#include <vector>

#include "sdl2_internal.hpp"
#include "surface.hpp"
#include "texture.hpp"

namespace std {
    template<>
    struct default_delete<SDL_Renderer> {
        inline void operator()(SDL_Renderer* ptr) {
            SDL_DestroyRenderer(ptr);
        }
    };
}

namespace SDL2 {
    class Window;

    class Renderer : public Object<Renderer> {
            friend class Object<Renderer>;

            friend class Window;

        public:
            Renderer() : Renderer(nullptr) {}

            inline void clear() {
                SDL_RenderClear(mp_Ptr.get());
            }

            inline void copy(const Texture& t) {
                if(SDL_RenderCopy(mp_Ptr.get(), t.mp_Ptr.get(), nullptr, nullptr) != 0)
                    throw Exception("Failed to copy texture to renderer");
            }

            inline void copy(const Texture& t, int x, int y) {
                Rect tSize = t.size();
                tSize.x = x;
                tSize.y = y;

                copy(t, tSize);
            }

            inline void copy(const Texture& t, const Rect& dst) {
                if(SDL_RenderCopy(mp_Ptr.get(), t.mp_Ptr.get(), nullptr, &dst) != 0)
                    throw Exception("Failed to copy texture to renderer");
            }

            inline void copy(const Texture& t, const Rect& src, int x, int y) {
                Rect dst = src;
                dst.x = x;
                dst.y = y;

                copy(t, src, dst);
            }

            inline void copy(const Texture& t, const Rect& src, const Rect& dst) {
                if(SDL_RenderCopy(mp_Ptr.get(), t.mp_Ptr.get(), &src, &dst) != 0)
                    throw Exception("Failed to copy texture to renderer");
            }

            inline void copy(const Texture& t, const Rect& dst, double angle) {
                if(SDL_RenderCopyEx(mp_Ptr.get(), t.mp_Ptr.get(), nullptr, &dst, angle, nullptr, SDL_FLIP_NONE) != 0)
                    throw Exception("Failed to copy texture to renderer");
            }

            inline void copy(const Texture& t, const Rect& src, const Rect& dst, double angle) {
                if(SDL_RenderCopyEx(mp_Ptr.get(), t.mp_Ptr.get(), &src, &dst, angle, nullptr, SDL_FLIP_NONE) != 0)
                    throw Exception("Failed to copy texture to renderer");
            }

            inline void copy(const Texture& t, const Rect& dst, double angle, const Point& center,
                             SDL_RendererFlip flip = SDL_FLIP_NONE) {
                if(SDL_RenderCopyEx(mp_Ptr.get(), t.mp_Ptr.get(), nullptr, &dst, angle, &center, flip) != 0)
                    throw Exception("Failed to copy texture to renderer");
            }

            inline void copy(const Texture& t, const Rect& src, const Rect& dst, double angle, const Point& center,
                             SDL_RendererFlip flip = SDL_FLIP_NONE) {
                if(SDL_RenderCopyEx(mp_Ptr.get(), t.mp_Ptr.get(), &src, &dst, angle, &center, flip) != 0)
                    throw Exception("Failed to copy texture to renderer");
            }

            inline Texture createTexture(const Surface& surf) {
                SDL_Texture* tPtr = SDL_CreateTextureFromSurface(mp_Ptr.get(), surf.mp_Ptr.get());

                if(!tPtr)
                    throw Exception("Failed to create texture from surface");

                return Texture{tPtr};
            }

            inline Texture createTexture(int w, int h, uint32_t format, int access = SDL_TEXTUREACCESS_STATIC) {
                SDL_Texture* tPtr = SDL_CreateTexture(mp_Ptr.get(), format, access, w, h);

                if(!tPtr)
                    throw Exception("Failed to create texture");

                return Texture{tPtr};
            }

#ifdef ENABLE_SDL2_GFX
            inline void drawQuadraticBezier(const Point& p1, const Point& anchor, const Point& p2) {
                drawQuadraticBezier(p1.x, p1.y, anchor.x, anchor.y, p2.x, p2.y, drawColour());
            }
            inline void drawQuadraticBezier(const Point& p1, const Point& anchor, const Point& p2, Colour c) {
                drawQuadraticBezier(p1.x, p1.y, anchor.x, anchor.y, p2.x, p2.y, c);
            }
            inline void drawQuadraticBezier(int16_t x1, int16_t y1, int16_t ax, int16_t ay, int16_t x2, int16_t y2) {
                drawQuadraticBezier(x1, y1, ax, ay, x2, y2, drawColour());
            }
            inline void drawQuadraticBezier(int16_t x1, int16_t y1, int16_t ax, int16_t ay, int16_t x2, int16_t y2, Colour c) {
                int16_t x[] = {x1, ax, x2};
                int16_t y[] = {y1, ay, y2};

                if(bezierRGBA(mp_Ptr.get(), x, y, 3, 2, c.r, c.g, c.b, c.a) != 0)
                    throw Exception("Failed to draw quadratic bezier curve");
            }
#endif

            inline void drawLine(const Point& p1, const Point& p2) {
                drawLine(p1.x, p1.y, p2.x, p2.y, drawColour());
            }

            inline void drawLine(const Point& p1, const Point& p2, Colour c) {
                drawLine(p1.x, p1.y, p2.x, p2.y, c);
            }

            inline void drawLine(int16_t x1, int16_t y1, int16_t x2, int16_t y2) {
                drawLine(x1, y1, x2, y2, drawColour());
            }

#ifdef ENABLE_SDL2_GFX
            inline void drawLine(int16_t x1, int16_t y1, int16_t x2, int16_t y2, Colour c) {
                if(aalineRGBA(mp_Ptr.get(), x1, y1, x2, y2, c.r, c.g, c.b, c.a) != 0)
                    throw Exception("Failed to draw line");
            }
#else

            inline void drawLine(int16_t x1, int16_t y1, int16_t x2, int16_t y2, Colour c) {
                Colour prevC = drawColour();
                setDrawColour(c);

                if(SDL_RenderDrawLine(mp_Ptr.get(), x1, y1, x2, y2) != 0)
                    throw Exception("Failed to draw line");

                setDrawColour(prevC);
            }

#endif

            inline void drawLines(const std::vector<Point>& points) {
                drawLines(points, drawColour());
            }

#ifdef ENABLE_SDL2_GFX
            inline void drawLines(const std::vector<Point>& points, Colour c) {
                for(size_t ii = 1; ii < points.size(); ii++) {
                    if(aalineRGBA(mp_Ptr.get(), points[ii - 1].x, points[ii - 1].y, points[ii].x, points[ii].y, c.r, c.g, c.b, c.a) != 0)
                        throw Exception("Failed to draw lines");
                }
            }
#else

            inline void drawLines(const std::vector<Point>& points, Colour c) {
                Colour prevC = drawColour();
                setDrawColour(c);

                if(SDL_RenderDrawLines(mp_Ptr.get(), points.data(), points.size()) != 0)
                    throw Exception("Failed to draw line");

                setDrawColour(prevC);
            }

#endif

            inline void drawPoint(const Point& p) {
                drawPoint(p.x, p.y);
            }

            inline void drawPoint(int x, int y) {
                if(SDL_RenderDrawPoint(mp_Ptr.get(), x, y) != 0)
                    throw Exception("Failed to draw point");
            }

            inline void drawPoints(const std::vector<Point>& points) {
                if(SDL_RenderDrawPoints(mp_Ptr.get(), points.data(), points.size()) != 0)
                    throw Exception("Failed to draw points");
            }

            inline void drawRectangle(const Rect& r) {
                drawRectangle(r, drawColour());
            }

#ifdef ENABLE_SDL2_GFX
            inline void drawRectangle(const Rect& r, Colour c) {
                if(rectangleRGBA(mp_Ptr.get(), r.x, r.y, r.x + r.w, r.y + r.h, c.r, c.g, c.b, c.a) != 0)
                    throw Exception("Failed to draw rectangle");
            }
#else

            inline void drawRectangle(const Rect& r, Colour c) {
                Colour prevC = drawColour();
                setDrawColour(c);

                if(SDL_RenderDrawRect(mp_Ptr.get(), &r) != 0)
                    throw Exception("Failed to draw rectangle");

                setDrawColour(prevC);
            }

#endif
#ifdef ENABLE_SDL2_GFX
            inline void drawThickLine(uint8_t w, const Point& p1, const Point& p2) {
                drawThickLine(w, p1.x, p1.y, p2.x, p2.y, drawColour());
            }
            inline void drawThickLine(uint8_t w, const Point& p1, const Point& p2, Colour c) {
                drawThickLine(w, p1.x, p1.y, p2.x, p2.y, c);
            }
            inline void drawThickLine(uint8_t w, int16_t x1, int16_t y1, int16_t x2, int16_t y2) {
                drawThickLine(w, x1, y1, x2, y2, drawColour());
            }
            inline void drawThickLine(uint8_t w, int16_t x1, int16_t y1, int16_t x2, int16_t y2, Colour c) {
                if(thickLineRGBA(mp_Ptr.get(), x1, y1, x2, y2, w, c.r, c.g, c.b, c.a) != 0)
                    throw Exception("Failed to draw thick line");
            }
#endif
#ifdef ENABLE_SDL2_GFX
            inline void fillCircle(int16_t cx, int16_t cy, int16_t r) {
                fillCircle(cx, cy, r, drawColour());
            }
            inline void fillCircle(int16_t cx, int16_t cy, int16_t r, Colour c) {
                if(filledCircleRGBA(mp_Ptr.get(), cx, cy, r, c.r, c.g, c.b, c.a) != 0)
                    throw Exception("Failed to fill circle");
            }
#endif

            inline void fill() {
                if(SDL_RenderFillRect(mp_Ptr.get(), nullptr) != 0)
                    throw Exception("Failed to fill rect");
            }

            inline void fillRect(int x, int y, int w, int h) {
                fillRect({x, y, w, h});
            }

            inline void fillRect(const Rect& rect) {
                if(SDL_RenderFillRect(mp_Ptr.get(), &rect) != 0)
                    throw Exception("Failed to fill rect");
            }

            inline void present() {
                SDL_RenderPresent(mp_Ptr.get());
            }

            inline SDL_BlendMode blendMode() const {
                SDL_BlendMode bm;

                if(SDL_GetRenderDrawBlendMode(mp_Ptr.get(), &bm) != 0)
                    throw Exception("Failed to get renderer draw blend mode");

                return bm;
            }

            inline void setBlendMode(SDL_BlendMode bm) {
                if(SDL_SetRenderDrawBlendMode(mp_Ptr.get(), bm) != 0)
                    throw Exception("Failed to set renderer draw blend mode");
            }

            inline Colour drawColour() const {
                Colour rgba;

                if(SDL_GetRenderDrawColor(mp_Ptr.get(), &rgba.r, &rgba.g, &rgba.b, &rgba.a) != 0)
                    throw Exception("Failed to get renderer draw colour");

                return rgba;
            }

            inline void setDrawColour(uint8_t r, uint8_t g, uint8_t b, uint8_t a = SDL_ALPHA_OPAQUE) {
                if(SDL_SetRenderDrawColor(mp_Ptr.get(), r, g, b, a) != 0)
                    throw Exception("Failed to set renderer draw colour");
            }

            inline void setDrawColour(Colour rgba) {
                setDrawColour(rgba.r, rgba.g, rgba.b, rgba.a);
            }

            inline void setTargetTexture(Texture& t) {
                if(SDL_SetRenderTarget(mp_Ptr.get(), t.mp_Ptr.get()) != 0)
                    throw Exception("Failed to set renderer target texture");
            }

            inline void unsetTargetTexture() {
                if(SDL_SetRenderTarget(mp_Ptr.get(), nullptr) != 0)
                    throw Exception("Failed to set renderer target texture");
            }

        private:
            explicit Renderer(SDL_Renderer* rPtr) : mp_Ptr{rPtr} {}

            std::unique_ptr<SDL_Renderer> mp_Ptr;
    };
}

#endif
