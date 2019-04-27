#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <cmath>

#ifdef M_PIl
#   define C_PI M_PIl
#else
#   define C_PI 3.141592653589793238462643383279502884L
#endif

namespace SI {
    constexpr long double C_EV              = 1.60217662e-19L;

    constexpr long double C_C               = 299792458;
    
    constexpr long double C_PLANCK          = 6.62607015e-34L;
    constexpr long double C_HBAR            = C_PLANCK/(2*C_PI);
    
    constexpr long double C_BIG_G           = 6.674e-11L;
    
    constexpr long double C_ELECTRON_MASS   = 5.485799090e-4L*1.6605390e-27L;
    constexpr long double C_PROTON_MASS     = 1.007276466L*1.6605390e-27L;
    constexpr long double C_NEUTRON_MASS    = 1.008664915L*1.6605390e-27L;
    constexpr long double C_NUCLEON_MASS    = (C_PROTON_MASS+C_NEUTRON_MASS)/2;
    
    constexpr long double C_SOLAR_MASS      = 1.9884e30L;
    constexpr long double C_SOLAR_RADIUS    = 6.95700e8L;
    
//     namespace astro {
//         constexpr long double C_SOLAR_MASS = 
//     };
}

#endif
