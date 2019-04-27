#ifndef RK4_HPP
#define RK4_HPP

#include <array>
#include <exception>

#include "types.hpp"

namespace RK4 {
    /*
     * Performs RK4 in an interval [t0, t1) with nSteps number of steps. Function f of type Func returns the derivative of
     * the y vector, i.e. dy/dt = f(y).
     *
     */
//     template <int N, typename Numeric, typename Func>
//     std::pair<data_t<Numeric>, std::array<data_t<Numeric>, N>> perform(vector_t<Numeric, N> y0, Numeric t0, const Numeric t1, const size_t nSteps, Func&& f) {
//         // Check that we are not going over the memory limit of 1 GB.
//         if((N + 1)*nSteps * sizeof(Numeric) > 1000000000)
//             throw std::runtime_error(format("too many steps (would use at least %.3f GB)", ((N + 1)*nSteps * sizeof(Numeric)) / 10e9));
//
//         // Performs a single RK4 step.
//         auto fn_step = [](vector_t<Numeric, N>& y, Numeric & t, const Numeric h, Func && dy) {
//             vector_t<Numeric, N> k1 = dy(y, t);
//             vector_t<Numeric, N> k2 = dy(y + (h / 2) * k1, t + h / 2);
//             vector_t<Numeric, N> k3 = dy(y + (h / 2) * k2, t + h / 2);
//             vector_t<Numeric, N> k4 = dy(y + h * k3, t + h);
//             y += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
//             t += h;
//         };
//
//         // Initialize output arrays.
//         std::pair<data_t<Numeric>, std::array<data_t<Numeric>, N>> output{data_t<Numeric>{nSteps}, {}};
//         output.second.fill(data_t<Numeric> {nSteps});
//
//         const Numeric h = (t1 - t0) / nSteps;
//
//         for(size_t ii = 0; ii < nSteps; ii++) {
//             output.first[ii] = t0;
//
//             for(size_t jj = 0; jj < N; jj++)
//                 output.second[jj][ii] = y0[jj];
//
//             fn_step(y0, t0, h, f);
//         }
//
//         return output;
//     }
//
//     template <int N, typename Numeric, typename Func>
//     inline std::pair<data_t<Numeric>, std::array<data_t<Numeric>, N>> perform(std::array<Numeric, N> y0, Numeric t0, const Numeric t1, const size_t nSteps, Func&& f) {
//         return perform(vector_t<Numeric, N>(y0.data()), t0, t1, nSteps, f);
//     }

    template<int N, class Func>
    std::pair<rdata_t, std::array<rdata_t, N>>
    perform(rvector_t<N> y0, real_t t0, const real_t t1, const size_t nSteps, Func&& f) {
        // Check that we are not going over the memory limit of 1 GB.
        if((N + 1) * nSteps * sizeof(real_t) > 1000000000)
            throw std::runtime_error(
                format("too many steps (would use at least %.3f GB)", ((N + 1) * nSteps * sizeof(real_t)) / 10e9));

        // Performs a single RK4 step.
        auto fn_step = [](rvector_t<N>& y, real_t& t, const real_t h, Func&& dy) {
            rvector_t<N> k1 = dy(y, t);
            rvector_t<N> k2 = dy(y + (h / 2) * k1, t + h / 2);
            rvector_t<N> k3 = dy(y + (h / 2) * k2, t + h / 2);
            rvector_t<N> k4 = dy(y + h * k3, t + h);
            y += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
            t += h;
        };

        // Initialize output arrays.
        std::pair<rdata_t, std::array<rdata_t, N>> output{rdata_t{nSteps}, {}};
        output.second.fill(rdata_t{nSteps});

        const real_t h = (t1 - t0) / nSteps;

        for(size_t ii = 0; ii < nSteps; ii++) {
            output.first[ii] = t0;

            for(size_t jj = 0; jj < N; jj++)
                output.second[jj][ii] = y0[jj];

            fn_step(y0, t0, h, f);
        }

        return output;
    }

    template<int N, class Func>
    inline std::pair<rdata_t, std::array<rdata_t, N>>
    perform(std::array<real_t, N> y0, real_t t0, const real_t t1, const size_t nSteps, Func&& f) {
        return perform(rvector_t<N>(y0.data()), t0, t1, nSteps, f);
    }

    template<int N, typename Func, class TestFunc>
    std::pair<rdata_t, std::array<rdata_t, N>>
    perform_with_test(rvector_t<N> y0, real_t t0, const real_t t1, const size_t nSteps, Func&& f, TestFunc&& tf) {
        // Check that we are not going over the memory limit of 1 GB.
        if((N + 1) * nSteps * sizeof(real_t) > 1000000000)
            throw std::runtime_error(
                format("too many steps (would use at least %.3f GB)", ((N + 1) * nSteps * sizeof(real_t)) / 10e9));

        // Performs a single RK4 step.
        auto fn_step = [](rvector_t<N>& y, real_t& t, const real_t h, Func&& dy) {
            rvector_t<N> k1 = dy(y, t);
            rvector_t<N> k2 = dy(y + (h / 2) * k1, t + h / 2);
            rvector_t<N> k3 = dy(y + (h / 2) * k2, t + h / 2);
            rvector_t<N> k4 = dy(y + h * k3, t + h);
            y += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
            t += h;
        };

        // Initialize output arrays.
        std::pair<rdata_t, std::array<rdata_t, N>> output{rdata_t{nSteps}, {}};
        output.second.fill(rdata_t{nSteps});

        const real_t h = (t1 - t0) / nSteps;

        for(size_t ii = 0; ii < nSteps; ii++) {
            output.first[ii] = t0;

            for(size_t jj = 0; jj < N; jj++)
                output.second[jj][ii] = y0[jj];

            if(!tf(y0))
                return output;

            fn_step(y0, t0, h, f);
        }

        return output;
    }

    template<int N, class Func, class TestFunc>
    inline std::pair<rdata_t, std::array<rdata_t, N>>
    perform_with_test(std::array<real_t, N> y0, real_t t0, const real_t t1, const size_t nSteps, Func&& f,
                      TestFunc&& tf) {
        return perform_with_test(rvector_t<N>(y0.data()), t0, t1, nSteps, f, tf);
    }
}

#endif
