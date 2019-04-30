#ifndef NONPOLYTROPIC_HPP
#define NONPOLYTROPIC_HPP

#include "constants.hpp"
#include "rungekutta/rk4.hpp"
#include "types.hpp"

namespace NeutronStar {
    namespace NonPolytropic {
        namespace detail {
            /**
             * Return the cube root of a positive real number.
             *
             * This function exists because math::cbrt accepts negative arguments, and
             * thus the integration does not terminate. We return NaN to force the
             * integration to terminate.
             *
             * Note that in any structure equation function we only need to use this
             * function once, since NaNs are propagated through any other operation
             * on a floating point number.
             **/
            template<typename FloatT>
            inline FloatT strict_cbrt(FloatT f) {
                return (f >= 0 ? math::cbrt(f) : std::numeric_limits<FloatT>::quiet_NaN());
            }
        }

        auto newtonianStructure() {
            return [](const rvector_t<2>& y, real_t r) -> rvector_t<2> {
                real_t r_2 = r * r;

                real_t x = detail::strict_cbrt(y[0]);
                real_t x_2 = x * x;
                real_t sqrt_1_x_2 = std::sqrt(1 + x_2);

                return {-y[1] * x * sqrt_1_x_2 / r_2, r_2 * y[0]};
            };
        }

        auto grStructure() {
            return [](const rvector_t<2>& y, real_t r) -> rvector_t<2> {
                real_t r_2 = r * r;

                real_t x = detail::strict_cbrt(y[0]);
                real_t x_2 = x * x;
                real_t x_3 = y[0];
                real_t sqrt_1_x_2 = std::sqrt(1 + x_2);
                real_t log_x_sqrt = std::log(x + sqrt_1_x_2);
                real_t I = (3 / (8 * x_3)) * (x * (1 + 2 * x_2) * sqrt_1_x_2 - log_x_sqrt);
                real_t Ip =
                    (3 * (-3 * x - x_3 + 2 * x_2 * x_3 + 3 * sqrt_1_x_2 * log_x_sqrt)) / (8 * x_2 * x_2 * sqrt_1_x_2);

                real_t gr1 = 1 + (1 / 3.0) * x * Ip / I;
                real_t gr2 = 1 + (r_2 * r * x_2 * x_2 * Ip) / (3 * y[1]);
                real_t gr3 = 1 / (1 - (2 * y[1] / (3*r)));

//                std::cout << "x: " << x << "; rho: " << y[0] << ";" << std::endl;
//                std::cout << "gr1: " << gr1 << "; gr2: " << gr2 << "; gr3: " << gr3 << ";" << std::endl;

                return {-y[1] * x * sqrt_1_x_2 * gr1 * gr2 * gr3 / r_2, r_2 * y[0]};
            };
        }

        template<typename Func, typename DFunc, typename DDFunc>
        auto eNuc(real_t K0, real_t C1, const Func& tF, const DFunc& tdF, const DDFunc& td2F) {
            /// I forgot to multiply K0 by 9 while calculating the A, B, sigma parameters, so we rescale it here.
            K0 /= 9;

            constexpr long double m3c3_over_hb3 =
                (SI::C_C * SI::C_NEUTRON_MASS / SI::C_HBAR) * (SI::C_C * SI::C_NEUTRON_MASS / SI::C_HBAR) *
                (SI::C_C * SI::C_NEUTRON_MASS / SI::C_HBAR);
            constexpr long double n0 = m3c3_over_hb3 / (3 * C_PI * C_PI);
            constexpr long double rho0 = SI::C_NEUTRON_MASS * n0;

            /// Calculate the normalisation constants for r and M.
            const real_t R_0 = std::sqrt(SI::C_C * SI::C_C / (12 * C_PI * rho0 * SI::C_BIG_G));
            const real_t M_0 = (SI::C_C * SI::C_C / (3 * SI::C_BIG_G)) * R_0;

            constexpr long double N0 = 0.16e45L;           // from pdf
            constexpr long double E_B = -16e6 * SI::C_EV;  // from pdf
            constexpr long double S0 = 30e6 * SI::C_EV;   // from pdf

            constexpr long double eta = N0 / n0;
            const long double cbrt_eta = math::cbrt(eta);
            const long double e_N = 3*cbrt_eta*cbrt_eta/10;
            const long double E_N = e_N * SI::C_NUCLEON_MASS * SI::C_C * SI::C_C;

            Func F{tF};
            DFunc dF{tdF};
            DDFunc d2F{td2F};

            /// Values of F and derivatives at u = 1.
            const long double F1 = F(1);
            const long double dF1 = dF(1);
            const long double d2F1 = d2F(1);

            real_t sigma__{0};
            {
                const real_t tNum11 = math::cbrt(std::pow(2, 5)) + 9 * (math::cbrt(4) - 1) * d2F1;
                const real_t tNum12 = C1 * (-math::cbrt(std::pow(2, 8)) + 18 * (math::cbrt(4) - 1) * F1);
                const real_t tNum13 = C1 * (-18 * (math::cbrt(4) - 1) * dF1 - 9 * d2F1 + 9 * math::cbrt(4) * d2F1);
                const real_t tNum1 = E_N * (tNum11 + tNum12 + tNum13) + 9 * (K0 + C1 * (2 * E_B + K0));
                const real_t tNum2 = -9 * S0 * (d2F1 + C1 * (2 * F1 - 2 * dF1 + d2F1));
                const real_t tDen1 = 3 * (C1 - 1);
                const real_t tDen21 =
                    E_N * (-math::cbrt(4) + 3 * (math::cbrt(4) - 1) * F1 - 3 * (math::cbrt(4) - 1) * dF1);
                const real_t tDen22 = 3 * (E_B + S0 * (dF1 - F1));

                sigma__ = (tNum1 + tNum2) / (tDen1 * (tDen21 + tDen22));
            }
            const real_t sigma = sigma__;

            real_t B__{0};
            {
                const real_t tNum1 = std::pow(1 + C1, 2) * (1 + sigma);
                const real_t tNum21 = 3 * (math::cbrt(4) - 1) * F1 - 3 * (math::cbrt(4) - 1) * dF1 - math::cbrt(4);
                const real_t tNum22 = 3 * (E_B + S0 * (-F1 + dF1)) / E_N;
                const real_t tNum2 = tNum21 + tNum22;
                const real_t tDen = 3 * (sigma - 1);

                B__ = -tNum1 * tNum2 / tDen;
            }
            const real_t B = B__;

            real_t A__{0};
            {
                const real_t t1 = -2 * B * (C1 + sigma) / ((1 + C1) * (1 + C1) * (1 + sigma));
                const real_t t2 = -2 * S0 * dF1 / E_N;
                const real_t t3 = 2 * (-math::cbrt(std::pow(2, 5)) + 3 * (math::cbrt(4) - 1) * dF1) / 3;

                A__ = t1 + t2 + t3;
            }
            const real_t A = A__;

//            std::cout << "A: " << A << "; B: " << B << "; sigma: " << sigma
//                      << ";" << std::endl;

            const real_t C_eta = (10 * std::cbrt(4 / (eta * eta)) / 9);

            const real_t cbrt_2_2 = std::cbrt(4);

            return [=](real_t u) -> real_t {

//                real_t sqrt_1_x_2 = std::sqrt(1 + x_2);
//                real_t log_x_sqrt = std::log(x + sqrt_1_x_2);
//                real_t I = (3 / (8 * x_3)) * (x * (1 + 2 * x_2) * sqrt_1_x_2 - log_x_sqrt);
//                real_t Ip =
//                    (3 * (-3 * x - x_3 + 2 * x_2 * x_3 + 3 * sqrt_1_x_2 * log_x_sqrt)) / (8 * x_2 * x_2 * sqrt_1_x_2);

                const real_t cbrt_u = detail::strict_cbrt(u);
                const real_t cbrt_u_2 = cbrt_u * cbrt_u;
                const real_t u_2 = u*u;
                const real_t cbrt_u_4 = cbrt_u_2 * cbrt_u_2;
                const real_t cbrt_1_u = 1/cbrt_u;
                const real_t cbrt_1_u_4 = 1/cbrt_u_4;
                const real_t u_sigma = std::pow(u, sigma);

                const real_t x = std::cbrt(u*eta);
//                const real_t x_2 = x * x;
                const  real_t x_3 = u*eta;

                const real_t Fu = F(u);
                const real_t dFu = dF(u);
                const real_t d2Fu = d2F(u);

                const real_t u_plus_C1_u_sigma = u + C1*u_sigma;
                const real_t u_plus_C1_u_sigma__2 = u_plus_C1_u_sigma*u_plus_C1_u_sigma;
                const real_t V = u*(A/2+B*u_sigma/((1+sigma)*u_plus_C1_u_sigma));
                const real_t Vp = A/2+(B*u_sigma*(C1*u_sigma + u*sigma))/((1+sigma)*u_plus_C1_u_sigma__2);
                const real_t Vpp = -(B*u_sigma*(sigma-1)*(C1*u_sigma*(sigma-2)-u*sigma))/((sigma + 1)*u_plus_C1_u_sigma__2*u_plus_C1_u_sigma);

                const real_t tJ1 = 2*cbrt_1_u/3 + Vp + ((cbrt_2_2 - 1)*(2*cbrt_1_u/3 - dFu)+S0*dFu/E_N);
                const real_t J = u_2*tJ1;
                const real_t Jp = 2*u*tJ1 + u_2*(-2*cbrt_1_u_4/9 + Vpp + ((cbrt_2_2 - 1)*(-2*cbrt_1_u_4/9 - d2Fu)+S0*d2Fu/E_N));

                const real_t epsilon = x_3*(1+cbrt_u_2+V+((cbrt_2_2 - 1)*(cbrt_u_2 - Fu)+S0*Fu/E_N));

                return SI::C_NUCLEON_MASS*SI::C_C*SI::C_C*(epsilon/x_3 - 1)/(SI::C_EV);
            };
        }

        template<typename Func, typename DFunc, typename DDFunc>
        auto nnInteractionStructure(real_t K0, real_t C1, const Func& tF, const DFunc& tdF, const DDFunc& td2F) {
            /// I forgot to multiply K0 by 9 while calculating the A, B, sigma parameters, so we rescale it here.
            K0 /= 9;

            constexpr long double m3c3_over_hb3 =
                (SI::C_C * SI::C_NEUTRON_MASS / SI::C_HBAR) * (SI::C_C * SI::C_NEUTRON_MASS / SI::C_HBAR) *
                (SI::C_C * SI::C_NEUTRON_MASS / SI::C_HBAR);
            constexpr long double n0 = m3c3_over_hb3 / (3 * C_PI * C_PI);
            constexpr long double rho0 = SI::C_NEUTRON_MASS * n0;

            /// Calculate the normalisation constants for r and M.
            const real_t R_0 = std::sqrt(SI::C_C * SI::C_C / (12 * C_PI * rho0 * SI::C_BIG_G));
            const real_t M_0 = (SI::C_C * SI::C_C / (3 * SI::C_BIG_G)) * R_0;

            constexpr long double N0 = 0.16e45L;           // from pdf
            constexpr long double E_B = -16e6 * SI::C_EV;  // from pdf
            constexpr long double S0 = 30e6 * SI::C_EV;   // from pdf

            constexpr long double eta = N0 / n0;
            const long double cbrt_eta = math::cbrt(eta);
            const long double e_N = 3*cbrt_eta*cbrt_eta/10;
            const long double E_N = e_N * SI::C_NUCLEON_MASS * SI::C_C * SI::C_C;

            Func F{tF};
            DFunc dF{tdF};
            DDFunc d2F{td2F};

            /// Values of F and derivatives at u = 1.
            const long double F1 = F(1);
            const long double dF1 = dF(1);
            const long double d2F1 = d2F(1);

            real_t sigma__{0};
            {
                const real_t tNum1 = 9*(K0+C1*(2*E_B+K0));
                const real_t tNum2 = (2-4*C1)*E_N;
                const real_t tDen = 3 * (C1 - 1)*(3*E_B-E_N);

                sigma__ = (tNum1 + tNum2) / tDen;
            }
            const real_t sigma = sigma__;

            real_t B__{0};
            {
                const real_t tNum = std::pow(1 + C1, 2) * (1 + sigma)*(3*E_B/E_N-1);
                const real_t tDen = 3 * (sigma - 1);

                B__ = -tNum / tDen;
            }
            const real_t B = B__;

            const real_t A = -2 * B * (C1 + sigma) / ((1 + C1) * (1 + C1) * (1 + sigma)) - 4/3;

//            std::cout << "A: " << A << "; B: " << B << "; sigma: " << sigma
//                      << ";" << std::endl;

            const real_t C_eta = (10 * std::cbrt(4 / (eta * eta)) / 9);

            const real_t cbrt_2_2 = std::cbrt(4);

            return [=](const rvector_t<2>& y, real_t r) -> rvector_t<2> {
                const real_t realR = r*R_0;

                const real_t rho = y[0];
                const real_t M = y[1];

                const real_t r_2 = r * r;


                const real_t x = detail::strict_cbrt(rho);
                const  real_t x_3 = rho;

                const real_t cbrt_u = x / cbrt_eta;
                const real_t cbrt_u_2 = cbrt_u * cbrt_u;
                const real_t u = x_3 / eta;
                const real_t u_2 = u*u;
                const real_t cbrt_u_4 = cbrt_u_2 * cbrt_u_2;
                const real_t cbrt_1_u = 1/cbrt_u;
                const real_t cbrt_1_u_4 = 1/cbrt_u_4;
                const real_t u_sigma = std::pow(u, sigma);

                const real_t Fu = F(u);
                const real_t dFu = dF(u);
                const real_t d2Fu = d2F(u);

                const real_t u_plus_C1_u_sigma = u + C1*u_sigma;
                const real_t u_plus_C1_u_sigma__2 = u_plus_C1_u_sigma*u_plus_C1_u_sigma;
                const real_t V = u*(A/2+B*u_sigma/((1+sigma)*u_plus_C1_u_sigma));
                const real_t Vp = A/2+(B*u_sigma*(C1*u_sigma + u*sigma))/((1+sigma)*u_plus_C1_u_sigma__2);
                const real_t Vpp = -(B*u_sigma*(sigma-1)*(C1*u_sigma*(sigma-2)-u*sigma))/((sigma + 1)*u_plus_C1_u_sigma__2*u_plus_C1_u_sigma);
                
                const real_t tJ1 = 2*cbrt_1_u/3 + Vp + ((cbrt_2_2 - 1)*(2*cbrt_1_u/3 - dFu)+S0*dFu/E_N);
                const real_t J = u_2*tJ1;
                const real_t Jp = 2*u*tJ1 + u_2*(-2*cbrt_1_u_4/9 + Vpp + ((cbrt_2_2 - 1)*(-2*cbrt_1_u_4/9 - d2Fu)+S0*d2Fu/E_N));

                const real_t epsilon = x_3*(1+cbrt_u_2+V+((cbrt_2_2 - 1)*(cbrt_u_2 - Fu)+S0*Fu/E_N));
                const real_t P = 3*std::cbrt(std::pow(eta, 5)/4)*J/10;

                const real_t gr1 = 1 + P/epsilon;
                const real_t gr2 = 1 + (r_2 * r * P) / M;
                const real_t gr3 = 1 / (1 - (2 * M / (3*r)));

                const real_t drho = - C_eta * (M * rho / Jp) * gr1 * gr2 * gr3 / r_2;
                const real_t dM = r_2 * rho;

                return {drho, dM};
            };
        }


        /**
         * Performs a RK4 integration with starting rho = startingDensity, up to a maximum r = endR and
         * using number of steps = steps.
         *
         * The integration is performed with the given structure equation.
         **/
        template<typename StructureFunc>
        auto integrate(real_t startingDensity, real_t endR, uint_t steps, StructureFunc&& fn_dy) {
            using namespace SI;

            /// Calculate some physical constants.
            constexpr long double m3c3_over_hb3 =
                (C_C * C_NEUTRON_MASS / C_HBAR) * (C_C * C_NEUTRON_MASS / C_HBAR) * (C_C * C_NEUTRON_MASS / C_HBAR);
            constexpr long double n0 = m3c3_over_hb3 / (3 * C_PI * C_PI);
            constexpr long double rho0 = C_NEUTRON_MASS * n0;

            /// Calculate the normalisation constants for r and M.
            real_t R_0 = std::sqrt(C_C * C_C / (12 * C_PI * rho0 * C_BIG_G));
            real_t M_0 = (C_C * C_C / (3 * C_BIG_G)) * R_0;

//                 std::cout << "n0 = " << n0 << ";" << std::endl;
            //     std::cout << "R0 = " << R_0 << " m; rho0 = " << rho0 << " kg/m^3; M0 = " << M_0 << " kg;" << std::endl;

            /// Calculate initial values of the relevant quantities (radius r, density rho and mass M).
            const real_t initialR = 1e-15;
            const real_t initialRho = startingDensity / rho0;
            const real_t initialM = 4 * C_PI * initialR * initialR * initialR * initialRho / 3;

//             std::cout << "initial dimensionless density: " << initialRho  << "; initial dimensionless mass: " << initialM <<  "; initial dimensionless radius: " << initialR << "; initial x: " << std::pow(initialRho, 1 / 3.0) << ";" << std::endl;

            uint_t lastIndex{0};
            real_t newRLimit{0.0};
            {
                /// We perform an initial integration between r = initialR and r = endR/R_0.
                /// The integration terminates prematurely when a NaN is encountered.
                auto data = RK4::perform_with_test<2>({{initialRho, initialM}}, initialR, endR / R_0, steps, fn_dy,
                                                      [](const rvector_t<2>& y) {
                                                          return isRegular(y[0]);
                                                      });

                /// We find the last index before a NaN is encountered.
                for(; lastIndex < data.second[0].size() && isRegular(data.second[0][lastIndex]); lastIndex++);
                lastIndex--;

                /// We calculate a new limit for r (we multiply by 1.1 to have some margin of error)
                newRLimit = 1.1 * data.first[lastIndex];
            }

            /// We perform the integration between r = initialR and r = newRLimit, the new limit for r
            /// we just calculated. Again, the integration terminates upon encountering a NaN.
            auto data = RK4::perform_with_test<2>({{initialRho, initialM}}, initialR, newRLimit, steps, fn_dy,
                                                  [](const rvector_t<2>& y) {
                                                      return isRegular(y[0]);
                                                  });

            /// We find the last index before a NaN is encountered.
            for(lastIndex = 0; lastIndex < data.second[0].size() && isRegular(data.second[0][lastIndex]); lastIndex++);
            lastIndex--;

            /// We erase the rubbish data after lastIndex.
            data.first = data.first.head(lastIndex + 1).eval();
            data.second[0] = data.second[0].head(lastIndex + 1).eval();
            data.second[1] = data.second[1].head(lastIndex + 1).eval();

            /// Convert to SI units.
            data.first *= R_0;
            data.second[0] *= rho0;
            data.second[1] *= M_0;

            return data;
        }
    }
}

#endif
