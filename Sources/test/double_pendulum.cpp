#include "sdl2/sdl2.hpp"

#include <vector>

#include "timer.hpp"
#include "types.hpp"

#include "../rungekutta/rk4.hpp"

/*
 * Returns the differential equation for a chaotic double pendulum.
 */
auto chaoticDoubleOscillator(real_t omega2) {
    return [=](const rvector_t<4>& y, real_t t) -> rvector_t<4> {
        real_t theta_1 = y[0];
        real_t theta_2 = y[2];
        real_t thetap_1 = y[1];
        real_t thetap_2 = y[3];
        real_t sin_1 = std::sin(theta_1);
        real_t sin_2 = std::sin(theta_2);
        real_t sin_d = std::sin(theta_1 - theta_2);
        real_t cos_d = std::cos(theta_1 - theta_2);

        return {thetap_1, -((3 * (-6 * omega2 * sin_1 + 3 * omega2 * cos_d * sin_2 -
                                  3 * cos_d * sin_d * thetap_1 * thetap_1 - 2 * sin_d * thetap_2 * thetap_2)) /
                            (-16 + 9 * cos_d * cos_d)), thetap_2, (3 *
                                                                   (-9 * omega2 * cos_d * sin_1 + 8 * omega2 * sin_2 -
                                                                    8 * sin_d * thetap_1 * thetap_1 -
                                                                    3 * cos_d * sin_d * thetap_2 * thetap_2)) /
                                                                  (-16 + 9 * cos_d * cos_d)};
    };
}

int main() {
    constexpr int WW = 1024, WH = 768;
    constexpr real_t T0 = 0, T1 = 100;
    constexpr uint_t NUM_SAMPLES = 2000000;
    constexpr uint_t SCALE = 175;

    auto fn_transform_x = [](real_t x) -> real_t {
        return SCALE * x + real_t(WW) / 2;
    };
    auto fn_transform_y = [](real_t y) -> real_t {
        return real_t(WH) / 2 - SCALE * y;
    };

    rdata_t x11, y11, x21, y21, x31, y31, x12, y12, x22, y22, x32, y32;

    {
        auto cdo = chaoticDoubleOscillator(9.81);

        {
            auto data = RK4::perform<4>({{M_PI_2l, 0, M_PI_2l - 0.001, 0}}, T0, T1, NUM_SAMPLES, cdo);

            x11 = (data.second[0].sin()).unaryExpr(fn_transform_x);
            y11 = (-data.second[0].cos()).unaryExpr(fn_transform_y);
            x12 = (data.second[0].sin() + data.second[2].sin()).unaryExpr(fn_transform_x);
            y12 = (-data.second[0].cos() - data.second[2].cos()).unaryExpr(fn_transform_y);
        }
        {
            auto data = RK4::perform<4>({{M_PI_2l, 0, M_PI_2l, 0}}, T0, T1, NUM_SAMPLES, cdo);

            x21 = (data.second[0].sin()).unaryExpr(fn_transform_x);
            y21 = (-data.second[0].cos()).unaryExpr(fn_transform_y);
            x22 = (data.second[0].sin() + data.second[2].sin()).unaryExpr(fn_transform_x);
            y22 = (-data.second[0].cos() - data.second[2].cos()).unaryExpr(fn_transform_y);
        }
        {
            auto data = RK4::perform<4>({{M_PI_2l, 0, M_PI_2l + 0.001, 0}}, T0, T1, NUM_SAMPLES, cdo);

            x31 = (data.second[0].sin()).unaryExpr(fn_transform_x);
            y31 = (-data.second[0].cos()).unaryExpr(fn_transform_y);
            x32 = (data.second[0].sin() + data.second[2].sin()).unaryExpr(fn_transform_x);
            y32 = (-data.second[0].cos() - data.second[2].cos()).unaryExpr(fn_transform_y);
        }
    }

    SDL2::init();

    SDL2::Window mainWindow("Double pendulum", WW, WH);
    SDL2::Renderer mainRend = mainWindow.createRenderer();
    mainRend.setBlendMode(SDL_BLENDMODE_BLEND);

    HRTimer timer("Frame time", {NUM_SAMPLES / (T1 - T0), ""});
    std::vector<SDL2::Point> v1, v2, v3;
    bool done = false;
    int ii = 0;

    v1.reserve(x11.size());
    v2.reserve(x21.size());
    v3.reserve(x31.size());

    while(!done) {
        timer.reset();
        timer.start();

        v1.push_back({static_cast<int>(x12[ii]), static_cast<int>(y12[ii])});
        v2.push_back({static_cast<int>(x22[ii]), static_cast<int>(y22[ii])});
        v3.push_back({static_cast<int>(x32[ii]), static_cast<int>(y32[ii])});


        mainRend.setDrawColour(255, 255, 255);
        mainRend.clear();

        mainRend.setDrawColour(0, 0, 0);
        mainRend.drawPoint(fn_transform_x(0), fn_transform_y(0));

        mainRend.setDrawColour(255, 0, 0);
        mainRend.drawLine(fn_transform_x(0), fn_transform_y(0), x11[ii], y11[ii]);
        mainRend.drawLine(x11[ii], y11[ii], x12[ii], y12[ii]);
        mainRend.setDrawColour(255, 0, 0, 100);
        mainRend.drawLines(v1);

        mainRend.setDrawColour(0, 0, 255);
        mainRend.drawLine(fn_transform_x(0), fn_transform_y(0), x31[ii], y31[ii]);
        mainRend.drawLine(x31[ii], y31[ii], x32[ii], y32[ii]);
        mainRend.setDrawColour(0, 0, 255, 100);
        mainRend.drawLines(v3);

        mainRend.setDrawColour(0, 255, 0);
        mainRend.drawLine(fn_transform_x(0), fn_transform_y(0), x21[ii], y21[ii]);
        mainRend.drawLine(x21[ii], y21[ii], x22[ii], y22[ii]);
        mainRend.setDrawColour(0, 255, 0, 100);
        mainRend.drawLines(v2);

        mainRend.present();


        SDL_Event ev;

        while(SDL_PollEvent(&ev)) {
            if(ev.type == SDL_QUIT) {
                done = SDL_TRUE;
            }
        }

        ii += timer.stop();

        if(ii >= NUM_SAMPLES)
            ii = NUM_SAMPLES - 1;
    }

    return 0;
}
