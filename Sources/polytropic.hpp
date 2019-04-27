#ifndef POLYTROPIC_HPP
#define POLYTROPIC_HPP

#include "constants.hpp"
#include "rungekutta/rk4.hpp"
#include "types.hpp"


namespace NeutronStar {
    namespace Polytropic {
        namespace {
            inline real_t I(real_t x) {
                return (3 / (8 * x * x * x)) *
                       (x * (1 + 2 * x * x) * std::sqrt(1 + x * x) - std::log(x + std::sqrt(1 + x * x)));
            }

            inline real_t Ip(real_t x) {
                return (3 * (-3 * x - std::pow(x, 3) + 2 * std::pow(x, 5) +
                             3 * std::sqrt(1 + x * x) * std::log(x + std::sqrt(1 + x * x)))) /
                       (8 * std::pow(x, 4) * std::sqrt(1 + x * x));
            }
        };

//         template <typename NonPStrucFunc>
//         auto structure(NonPStrucFunc fn_structure) {
//             constexpr real_t A_NR = 9.47818e-10L, A_R = 3.00054L;
//
//             return [=](const rvector_t<2>& y, real_t r) -> rvector_t<2> {
//                 real_t epsilon = A_NR*std::pow
//             };
//         }

        auto grStructure() {
            return [](const rvector_t<2>& y, real_t r) -> rvector_t<2> {
                constexpr real_t A_NR = 2.61290005369166L, A_R = 2.66988055730438L;     // from 0.001 to  1 in steps of 0.001
//                 constexpr real_t A_NR = 2.120866954664444L, A_R = 2.947374853289425L;   // from 0.001 to  6 in steps of 0.001
//                 constexpr real_t A_NR = 1.7271494343102791L, A_R = 2.983893795320325L;  // from 0.001 to 12 in steps of 0.001
//                 constexpr real_t A_NR = 0.94735296709173L, A_R = 2.999273659636282L;    // from 0.001 to 60 in steps of 0.001

                real_t r_2 = r * r;
                real_t epsilon = A_NR * std::pow(y[0], 3 / 5.0) + A_R * y[0];

                real_t gr1 = 1 + y[0] / epsilon;
                real_t gr2 = 1 + (r_2 * r * y[0]) / (y[1]);
                real_t gr3 = 1 / (1 - (2 * y[1] / r));

                //         std::cout << "gr1: " << gr1 << "; gr2: " << gr2 << "; gr3: " << gr3 << ";" << std::endl;

                return {-y[1] * epsilon * gr1 * gr2 * gr3 / (r * r), r * r * epsilon};
            };
        }

        auto newtonianStructure() {
            return [](const rvector_t<2>& y, real_t r) -> rvector_t<2> {
                constexpr real_t A_NR = 2.61290005369166L, A_R = 2.66988055730438L;     // from 0.001 to  1 in steps of 0.001
//                 constexpr real_t A_NR = 2.120866954664444L, A_R = 2.947374853289425L;   // from 0.001 to  6 in steps of 0.001
//                 constexpr real_t A_NR = 1.7271494343102791L, A_R = 2.983893795320325L;  // from 0.001 to 12 in steps of 0.001
//                 constexpr real_t A_NR = 0.94735296709173L, A_R = 2.999273659636282L;    // from 0.001 to 60 in steps of 0.001

                real_t epsilon = A_NR * std::pow(y[0], 3 / 5.0) + A_R * y[0];

                return {-y[1] * epsilon / (r * r), r * r * epsilon};
            };
        }

        template<typename StructureFunc>
        auto integrate(real_t startingDensity, real_t endR, uint_t steps, StructureFunc&& fn_dy) {
            using namespace SI;

            constexpr long double m3c3_over_hb3 =
                (C_C * C_NEUTRON_MASS / C_HBAR) * (C_C * C_NEUTRON_MASS / C_HBAR) * (C_C * C_NEUTRON_MASS / C_HBAR);
            constexpr long double n0 = m3c3_over_hb3 / (3 * C_PI * C_PI);
            constexpr long double rho0 = C_NEUTRON_MASS * n0;
            constexpr long double epsilon0 = C_C * C_C * rho0;

            constexpr long double pi4_e0_G = 4 * C_PI * epsilon0 * C_BIG_G;

            real_t R_0 = C_C * C_C / std::sqrt(pi4_e0_G);
            real_t M_0 = R_0 * C_C * C_C / C_BIG_G;

//             std::cout << "n0 = " << n0 << ";" << std::endl;
//             std::cout << "R0 = " << R_0 << " m; epsilon0 = " << epsilon0 << " kg/m^3; M0 = " << M_0 << " kg;" << std::endl;

            const real_t initialX = std::pow(startingDensity / rho0, 1 / 3.0);

            const real_t initialR = 1e-15;
            const real_t initialP = std::pow(initialX, 4) * Ip(initialX) / 3.0;
            const real_t initialM = 4 * C_PI * initialR * initialR * initialR * (startingDensity / rho0) / 3;

//             std::cout << "initial dimensionless pressure: " << initialP  << "; initial dimensionless mass: " << initialM <<  "; initial dimensionless radius: " << initialR << ";" << std::endl;

            uint_t lastRealIndex{0};
            real_t newRLimit{0.0};
            {
                auto data = RK4::perform<2>({{initialP, initialM}}, initialR, endR / R_0, steps, fn_dy);

                for(;
                    lastRealIndex < data.second[0].size() && isRegular(data.second[0][lastRealIndex]); lastRealIndex++);

                lastRealIndex--;

                newRLimit = 1.1 * data.first[lastRealIndex];
            }

            auto data = RK4::perform<2>({{initialP, initialM}}, initialR, newRLimit, steps, fn_dy);

            for(lastRealIndex = 0;
                lastRealIndex < data.second[0].size() && isRegular(data.second[0][lastRealIndex]); lastRealIndex++);

            lastRealIndex--;

            data.first = data.first.head(lastRealIndex + 1).eval();
            data.second[0] = data.second[0].head(lastRealIndex + 1).eval();
            data.second[1] = data.second[1].head(lastRealIndex + 1).eval();

            data.first *= R_0;
            data.second[0] *= epsilon0;
            data.second[1] *= M_0;

            return data;
        }
    }
}

#endif
