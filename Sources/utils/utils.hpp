#ifndef COMPUTINGPROJECT_UTILS_HPP
#define COMPUTINGPROJECT_UTILS_HPP

#include <utility>

#include "types.hpp"

namespace Utils {
//    inline constexpr uint_t factorial(uint_t n) {
//        return n == 0 ? 1 : n * factorial(n-1);
//    }
//
//    inline constexpr real_t pow(real_t x, uint_t y) {
//        return y == 0 ? 1.0 : x * pow(x, y-1);
//    }
//
//    inline constexpr real_t pow(real_t x, real_t y) {
//        return 1.0 + ;
//    }

    template<typename Func>
    inline std::pair <rdata_t, rdata_t> table(real_t t0, real_t t1, uint_t numSamples, Func&& f) {
        if(2 * numSamples * sizeof(real_t) > 100000000)
            throw std::runtime_error("too many samples");

        real_t step = (t1 - t0) / numSamples;
        std::pair <rdata_t, rdata_t> data = std::make_pair(rdata_t(numSamples), rdata_t(numSamples));

        for(uint_t ii = 0; ii < numSamples; ii++) {
            real_t t = t0 + step * ii;
            data.first[ii] = t;
            data.second[ii] = f(t);
        }

        return data;
    }
}

namespace math {
    inline real_t cbrt(std::conditional_t<std::is_fundamental<real_t>::value, real_t, const real_t&> f) {
        return std::cbrt(f);
    }
}

#endif //COMPUTINGPROJECT_UTILS_HPP
