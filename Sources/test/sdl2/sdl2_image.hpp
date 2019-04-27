#ifndef SDL2_IMAGE_HPP
#define SDL2_IMAGE_HPP

#include "sdl2_internal.hpp"

#include "surface.hpp"

namespace SDL2 {
    inline void initIMG(uint32_t flags = (IMG_INIT_JPG | IMG_INIT_PNG | IMG_INIT_TIF | IMG_INIT_WEBP)) {
        if(IMG_Init(flags) != flags)
            throw Exception("Failed to initialise SDL2_image");
    }
    
    inline Surface loadImage(const std::string& filename) {
        SDL_Surface* imgSurf = IMG_Load(filename.c_str());
        
        if(!imgSurf)
            throw Exception("Failed to load image '" + filename + "'");
        
        return Surface{imgSurf};
    }
    
    inline void quitIMG() {
        IMG_Quit();
    }
}

#endif
