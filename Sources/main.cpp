#include "rungekutta/rk4.hpp"

#include <cmath>

#include <fstream>

#include "matplotlib.hpp"

#include "constants.hpp"
#include "timer.hpp"
#include "types.hpp"
#include "utils.hpp"

#include "nonpolytropic.hpp"
#include "polytropic.hpp"

namespace plt = matplotlibcpp;
namespace ns = NeutronStar::NonPolytropic;
// namespace ns = NeutronStar::Polytropic;

template<typename T>
inline const std::string to_string(T val) {
    std::ostringstream ss;
    ss << val;
    return ss.str();
}

int main() {
    using namespace SI;

    plt::init();

//    auto F = [](real_t u) -> real_t {
//        return u;
//    };
//    auto dF = [](real_t u) -> real_t {
//        return 1;
//    };
//    auto d2F = [](real_t u) -> real_t {
//        return 0;
//    };

    auto F = [](real_t u) -> real_t {
        return std::sqrt(u);
    };
    auto dF = [](real_t u) -> real_t {
        return 1 / (2 * std::sqrt(u));
    };
    auto d2F = [](real_t u) -> real_t {
        return -1 / (4 * std::sqrt(u * u * u));
    };

//    auto F = [](real_t u) -> real_t {
//        return u*u;
//    };
//    auto dF = [](real_t u) -> real_t {
//        return 2*u;
//    };
//    auto d2F = [](real_t u) -> real_t {
//        return 2;
//    };

    if(/* DISABLES CODE */ (false)) {
        plt::title("Neutron star");

        constexpr uint_t N_STEPS = 200000;
        const real_t rhoI = 5e21;
//        const real_t rhoI = 1e21;

        HRTimer timer("Neutron star");
        timer.start();

        auto data = ns::integrate(rhoI, 1000000, N_STEPS, ns::nnInteractionStructure(400e6 * C_EV, 0.1, F, dF, d2F));
//        auto data = ns::integrate(rhoI, 1000000, N_STEPS, ns::grStructure());
        uint_t dataSetSize = data.first.size();

        timer.stop();
        timer.report();

        std::cout << "Density at center: " << rhoI << " kg/m^3" << std::endl;
        std::cout << "Radius: " << data.first[dataSetSize - 1] << " m ("
                  << (data.first[dataSetSize - 1] / C_SOLAR_RADIUS) << " solar radii)" << std::endl;
        std::cout << "Mass: " << data.second[1][dataSetSize - 1] << " kg ("
                  << (data.second[1][dataSetSize - 1] / C_SOLAR_MASS) << " solar masses)" << std::endl;

        plt::plot(data.first, data.second[0]);

//         plt::legend();
//         plt::save("wd.pdf");
        plt::show();
    }

    plt::clf();

    if(/* DISABLES CODE */ (false)) {
        plt::title("Neutron star");

        constexpr uint_t N_STEPS = 200000;
//         const real_t rhoI = 5e16;
        const real_t rhoI = 1e21;

        HRTimer timer("Neutron star");
        timer.start();

//        auto data = ns::integrate(rhoI, 100000, N_STEPS, ns::nnInteractionStructure(400*1e6*C_EV, 0.2));
//        uint_t dataSetSize = data.first.size();
//
//        timer.stop();
//        timer.report();
//
//        plt::plot(data.first, data.second[0]);

        for(int ii = 200; ii <= 400; ii += 25) {
            auto data = Utils::table(0, 20, 100000, ns::eNuc(ii * 1e6 * C_EV, 0.2, F, dF, d2F));
            timer.pause();

            plt::named_plot(format("K0 = %d MeV", ii), data.first, data.second);
        }
        timer.stop();
        timer.report();

        plt::legend();
        plt::save("E_nuc.pdf");
        plt::show();
    }

    plt::clf();

    if(/* DISABLES CODE */ (false)) {
        plt::title("Neutron star");

        constexpr uint_t N_STEPS = 200000;
        constexpr uint_t N_STARS = 6;
        HRTimer timer("Neutron star");

        rmdata_t tab(N_STARS, 3);

        for(uint_t ii = 0; ii < N_STARS; ii++) {
            real_t rhoI = std::pow(10, 17 + 0.5 * ii);

            timer.start();
            auto data = ns::integrate(rhoI, 100000, N_STEPS, ns::grStructure());
            timer.pause();

            uint_t dataSetSize = data.first.size();
            tab(ii, 0) = rhoI;
            tab(ii, 1) = (data.second[1][dataSetSize - 1] / C_SOLAR_MASS);
            tab(ii, 2) = data.first[dataSetSize - 1] / 1000;

//             plt::named_plot("central rho = " + to_string(rhoI) + " kg/m^3", (data.first / 1000).eval(), (data.second[1] / C_SOLAR_MASS).eval());
//             plt::named_plot("central rho = " + to_string(rhoI) + " kg/m^3", (data.first / 1000).eval(), (data.second[0]/data.second[0][0]).eval());
        }

        timer.stop();
        timer.report();

        std::cout << tab << std::endl;

        plt::legend();
        plt::save("wd.pdf");
        plt::show();
    }

    plt::clf();

    if(/* DISABLES CODE */ (false)) {
//        const std::string pref = "nw";
        const std::string pref = "gr";
//        const std::string pref = "nn";

        constexpr uint_t N_STEPS = 200000;
        constexpr uint_t N_STARS = 200;
        constexpr real_t MIN_RHO = 1e17;
        constexpr real_t MAX_RHO = 5e20;

        plt::title("Mass vs radius of neutron stars (Fermi gas model)");

        HRTimer timer("Neutron star (Fermi gas)");


        rmdata_t tab(N_STARS, 3);

        timer.start();
#pragma omp parallel for

        for(uint32_t ii = 0; ii < N_STARS; ii++) {
            real_t rhoI = MIN_RHO * std::pow(10, ii * (std::log10(MAX_RHO / MIN_RHO)) / N_STARS);

            auto data = ns::integrate(rhoI, 1000000, N_STEPS, ns::grStructure());

            uint_t dataSetSize = data.first.size();
            tab(ii, 0) = rhoI;
            tab(ii, 1) = (data.second[1][dataSetSize - 1] / C_SOLAR_MASS);
            tab(ii, 2) = data.first[dataSetSize - 1] / 1000;
        }


        timer.stop();
        timer.report();

//         std::cout << tab << std::endl;

        plt::plot(tab.col(2), tab.col(1));

        plt::xlabel("Radius (km)");
        plt::ylabel("Mass (solar masses)");
        plt::save(pref + "_mass_v_radius.pdf");
    }

    plt::clf();

    if(/* DISABLES CODE */ (true)) {
        const std::string pref = "nw";
//        const std::string pref = "gr";
//        const std::string pref = "nn";

        constexpr uint_t N_STEPS = 200000;
        constexpr uint_t N_STARS = 200;
        constexpr real_t MIN_RHO = 1e17;
        constexpr real_t MAX_RHO = 5e20;

        plt::title("Mass vs radius of neutron stars (Newtonian structure)");

        HRTimer timer("Neutron star (Newtonian structure)");


        rmdata_t tab(N_STARS, 3);

        timer.start();
#pragma omp parallel for

        for(uint32_t ii = 0; ii < N_STARS; ii++) {
            real_t rhoI = MIN_RHO * std::pow(10, ii * (std::log10(MAX_RHO / MIN_RHO)) / N_STARS);

            auto data = ns::integrate(rhoI, 1000000, N_STEPS, ns::newtonianStructure());

            uint_t dataSetSize = data.first.size();
            tab(ii, 0) = rhoI;
            tab(ii, 1) = (data.second[1][dataSetSize - 1] / C_SOLAR_MASS);
            tab(ii, 2) = data.first[dataSetSize - 1] / 1000;
        }

        timer.pause();

        plt::plot(tab.col(2), tab.col(1));

        std::ofstream(pref + "_mri.dat") << tab;

        timer.stop();
        timer.report();

//         std::cout << tab << std::endl;


        plt::xlabel("Radius (km)");
        plt::ylabel("Mass (solar masses)");
        plt::save(pref + "_mass_v_radius.pdf");

        plt::clf();

//        plt::plot((tab.col(0).log10()).eval(), tab.col(2));
//        plt::save(pref + "_radius_v_rhoi.pdf");
//
//        plt::clf();
//
//        plt::plot((tab.col(0).log10()).eval(), tab.col(1));
//        plt::save(pref + "_mass_v_rhoi.pdf");
    }

    plt::clf();

    if(/* DISABLES CODE */ (true)) {
//        const std::string pref = "nw";
        const std::string pref = "gr";
//        const std::string pref = "nn";

        constexpr uint_t N_STEPS = 200000;
        constexpr uint_t N_STARS = 200;
        constexpr real_t MIN_RHO = 1e16;
        constexpr real_t MAX_RHO = 5e20;

        plt::title("Mass vs radius of neutron stars (GR structure)");

        HRTimer timer("Neutron star (GR structure)");


        rmdata_t tab(N_STARS, 3);

        timer.start();
#pragma omp parallel for

        for(uint32_t ii = 0; ii < N_STARS; ii++) {
            real_t rhoI = MIN_RHO * std::pow(10, ii * (std::log10(MAX_RHO / MIN_RHO)) / N_STARS);

            auto data = ns::integrate(rhoI, 1000000, N_STEPS, ns::grStructure());

            uint_t dataSetSize = data.first.size();
            tab(ii, 0) = rhoI;
            tab(ii, 1) = (data.second[1][dataSetSize - 1] / C_SOLAR_MASS);
            tab(ii, 2) = data.first[dataSetSize - 1] / 1000;
        }

        timer.pause();

        plt::plot(tab.col(2), tab.col(1));

        std::ofstream(pref + "_mri.dat") << tab;

        timer.stop();
        timer.report();

//         std::cout << tab << std::endl;


        plt::xlabel("Radius (km)");
        plt::ylabel("Mass (solar masses)");
        plt::save(pref + "_mass_v_radius.pdf");

        plt::clf();

//        plt::plot((tab.col(0).log10()).eval(), tab.col(2));
//        plt::save(pref + "_radius_v_rhoi.pdf");
//
//        plt::clf();
//
//        plt::plot((tab.col(0).log10()).eval(), tab.col(1));
//        plt::save(pref + "_mass_v_rhoi.pdf");
    }

    plt::clf();

    if(/* DISABLES CODE */ (true)) {
//        const std::string pref = "nw";
//        const std::string pref = "gr";
        const std::string pref = "nn";

        constexpr uint_t N_STEPS = 200000;
        constexpr uint_t N_STARS = 200;
        constexpr real_t MIN_RHO = 1e17;
        constexpr real_t MAX_RHO = 5e20;

        plt::title("Mass vs radius of neutron stars (varying K0)");

        HRTimer timer("Neutron star (varying K0)");


        for(int jj = 100; jj <= 400; jj += 75) {
            rmdata_t tab(N_STARS, 3);

            timer.start();
#pragma omp parallel for

            for(uint32_t ii = 0; ii < N_STARS; ii++) {
                real_t rhoI = MIN_RHO * std::pow(10, ii * (std::log10(MAX_RHO / MIN_RHO)) / N_STARS);

                auto data = ns::integrate(rhoI, 1000000, N_STEPS,
                                          ns::nnInteractionStructure(jj * 1e6 * C_EV, 0.2, F, dF, d2F));

                uint_t dataSetSize = data.first.size();
                tab(ii, 0) = rhoI;
                tab(ii, 1) = (data.second[1][dataSetSize - 1] / C_SOLAR_MASS);
                tab(ii, 2) = data.first[dataSetSize - 1] / 1000;
            }

            timer.pause();

            plt::named_plot(format("K0 = %d MeV", jj), tab.col(2), tab.col(1));

            std::ofstream(format("ns_mri_(C1_0.2)_k0_%d.dat", jj)) << tab;
        }

        timer.stop();
        timer.report();

        plt::xlabel("Radius (km)");
        plt::ylabel("Mass (solar masses)");
        plt::legend();
        plt::save(pref + "_mass_v_radius_k0.pdf");
    }

    plt::clf();

    if(/* DISABLES CODE */ (true)) {
//        const std::string pref = "nw";
//        const std::string pref = "gr";
        const std::string pref = "nn";

        constexpr uint_t N_STEPS = 200000;
        constexpr uint_t N_STARS = 200;
        constexpr real_t MIN_RHO = 1e17;
        constexpr real_t MAX_RHO = 5e20;

        plt::title("Mass vs radius of neutron stars (varying C1)");

        HRTimer timer("Neutron star (varying C1)");


        for(real_t jj = 0; jj <= 0.5; jj += 0.1) {
            rmdata_t tab(N_STARS, 3);

            timer.start();
#pragma omp parallel for

            for(uint32_t ii = 0; ii < N_STARS; ii++) {
                real_t rhoI = MIN_RHO * std::pow(10, ii * (std::log10(MAX_RHO / MIN_RHO)) / N_STARS);

                auto data = ns::integrate(rhoI, 1000000, N_STEPS,
                                          ns::nnInteractionStructure(220 * 1e6 * C_EV, jj, F, dF, d2F));

                uint_t dataSetSize = data.first.size();
                tab(ii, 0) = rhoI;
                tab(ii, 1) = (data.second[1][dataSetSize - 1] / C_SOLAR_MASS);
                tab(ii, 2) = data.first[dataSetSize - 1] / 1000;
            }

            timer.pause();

            plt::named_plot("C1 = " + to_string(jj), tab.col(2), tab.col(1));

            std::ofstream("ns_mri_(K0_220_MeV)_C1_" + to_string(jj) + ".dat") << tab;
        }

        timer.stop();
        timer.report();

//         std::cout << tab << std::endl;


        plt::xlabel("Radius (km)");
        plt::ylabel("Mass (solar masses)");
        plt::legend();
        plt::save(pref + "_mass_v_radius_c1.pdf");

        plt::clf();

//        plt::plot((tab.col(0).log10()).eval(), tab.col(2));
//        plt::save(pref + "_radius_v_rhoi.pdf");
//
//        plt::clf();
//
//        plt::plot((tab.col(0).log10()).eval(), tab.col(1));
//        plt::save(pref + "_mass_v_rhoi.pdf");
    }

    plt::clf();

    if(/* DISABLES CODE */ (false)) {
        //         const std::string pref = "nw";
        const std::string pref = "gr";

        constexpr uint_t N_STEPS = 200000;
        constexpr uint_t N_STARS = 200;
        constexpr real_t MIN_RHO = 1e17;
        constexpr real_t MAX_RHO = 3e24;
        HRTimer timer("Neutron star");

        rmdata_t tabNonP(N_STARS, 3);
        rmdata_t tabPoly(N_STARS, 3);

        timer.start();

#pragma omp parallel for

        for(uint32_t ii = 0; ii < N_STARS; ii++) {
            real_t rhoI = MIN_RHO * std::pow(10, ii * (std::log10(MAX_RHO / MIN_RHO)) / N_STARS);

            {
                namespace nsns = NeutronStar::NonPolytropic;

                auto data = nsns::integrate(rhoI, 100000, N_STEPS,
                                            (pref == "nw" ? nsns::newtonianStructure() : nsns::grStructure()));

                uint_t dataSetSize = data.first.size();
                tabNonP(ii, 0) = rhoI;
                tabNonP(ii, 1) = (data.second[1][dataSetSize - 1] / C_SOLAR_MASS);
                tabNonP(ii, 2) = data.first[dataSetSize - 1] / 1000;
            }
            {
                namespace nsns = NeutronStar::Polytropic;

                auto data = nsns::integrate(rhoI, 100000, N_STEPS,
                                            (pref == "nw" ? nsns::newtonianStructure() : nsns::grStructure()));

                uint_t dataSetSize = data.first.size();
                tabPoly(ii, 0) = rhoI;
                tabPoly(ii, 1) = (data.second[1][dataSetSize - 1] / C_SOLAR_MASS);
                tabPoly(ii, 2) = data.first[dataSetSize - 1] / 1000;
            }
        }

        timer.stop();
        timer.report();

        plt::title("Mass vs radius of neutron stars");
        plt::named_plot("Non polytropic model", tabNonP.col(2), tabNonP.col(1));
        plt::named_plot("Polytropic model", tabPoly.col(2), tabPoly.col(1));
        plt::xlabel("Radius (km)");
        plt::ylabel("Mass (solar masses)");
        plt::legend();
        plt::save(pref + "_mass_v_radius.pdf");

        plt::clf();

        plt::named_plot("Non polytropic model", (tabNonP.col(0).log10()).eval(), tabNonP.col(2));
        plt::named_plot("Polytropic model", (tabPoly.col(0).log10()).eval(), tabPoly.col(2));
        plt::legend();
        plt::save(pref + "_radius_v_rhoi.pdf");

        plt::clf();

        plt::named_plot("Non polytropic model", (tabNonP.col(0).log10()).eval(), tabNonP.col(1));
        plt::named_plot("Polytropic model", (tabPoly.col(0).log10()).eval(), tabPoly.col(1));
        plt::legend();
        plt::save(pref + "_mass_v_rhoi.pdf");
    }

//    if(/* DISABLES_CODE */ (true)) {
//        constexpr real_t K0 = 100, C1 = 0.1, E_B = 1, F1 = 1, dF1 = 1, d2F1 = 1, E_N = 1, S0 = 1;
//
//        const real_t tNum11 = math::cbrt(std::pow(2, 5)) + 9 * (math::cbrt(4) - 1) * d2F1;
//        const real_t tNum12 = C1 * (-math::cbrt(std::pow(2, 8)) + 18 * (math::cbrt(4) - 1) * F1);
//        const real_t tNum13 = C1 * (-18 * (math::cbrt(4) - 1) * dF1 - 9 * d2F1 + 9 * math::cbrt(4) * d2F1);
//        const real_t tNum1 = E_N * (tNum11 + tNum12 + tNum13) + 9 * (K0 + C1 * (2 * E_B + K0));
//        const real_t tNum2 = -9 * S0 * (d2F1 + C1 * (2 * F1 - 2 * dF1 + d2F1));
//        const real_t tDen1 = 3 * (C1 - 1);
//        const real_t tDen21 = E_N * (-math::cbrt(4) + 3 * (math::cbrt(4) - 1) * F1 - 3 * (math::cbrt(4) - 1) * dF1);
//        const real_t tDen22 = 3 * (E_B + S0 * (dF1 - F1));
//
//        const real_t sigma__ = (tNum1 + tNum2) / (tDen1 * (tDen21 + tDen22));
//
//        std::cout << sigma__ << std::endl;
//    }

//    if(/* DISABLES_CODE */ (true)) {
//        constexpr real_t K0 = 100, C1 = 0.1, E_B = 1, F1 = 1, dF1 = 1, d2F1 = 0, E_N = 1, S0 = 1;
//        constexpr real_t sigma = 1.2;
//
//        const real_t tNum1 = std::pow(1 + C1, 2) * (1 + sigma);
//        const real_t tNum21 = E_N * (3 * (math::cbrt(4) - 1) * F1 - 3 * (math::cbrt(4) - 1) * dF1 - math::cbrt(4));
//        const real_t tNum22 = 3*(E_B+S0*(-F1+dF1));
//        const real_t tNum2 = tNum21+tNum22;
//        const real_t tDen = 3 * (sigma - 1);
//
//        const real_t B__ = -tNum1*tNum2/tDen;
//
//        std::cout << B__ << std::endl;
//    }

//    if(/* DISABLES_CODE */ (true)) {
//        constexpr real_t K0 = 100, C1 = 0.1, E_B = 1, F1 = 1, dF1 = 1, d2F1 = 0, E_N = 1, S0 = 1;
//        constexpr real_t sigma = 1.2, B = 1.4;
//
//        const real_t t1 = -2*B*(C1+sigma)/((1+C1)*(1+C1)*(1+sigma));
//        const real_t t2 = -2*S0*dF1/E_N;
//        const real_t t3 = 2*(-math::cbrt(std::pow(2, 5)) + 3*(math::cbrt(4) - 1)*dF1)/3;
//
//        const real_t A__ = t1 +t2 + t3;
//
//        std::cout << A__ << std::endl;
//    }
}
