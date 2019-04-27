#ifndef TYPES_HPP
#define TYPES_HPP

#include <cstdint>

#include <iostream>

#include <Eigen/Dense>

// // #ifndef WITH_LOW_PRECISION
// // #ifdef VERY_HIGH_PRECISION
// // #include <boost/multiprecision/gmp.hpp>
// 
// //very high precision
// typedef boost::multiprecision::mpz_int sint_t;
// typedef boost::multiprecision::mpz_int uint_t;
// typedef boost::multiprecision::mpf_float_100 real_t;
// // #else
//// high precision   (for actual results)
//typedef __int128            sint_t;
//typedef unsigned __int128   uint_t;
//typedef long double         real_t;
////
//std::ostream& operator<<(std::ostream& dest, __int128);
//std::ostream& operator<<(std::ostream& dest, unsigned __int128);
// // #endif
// // #else
//// low precision   (for quick prototyping and testing)
typedef int64_t sint_t;
typedef uint64_t uint_t;
typedef double real_t;
// // #endif

typedef std::complex<real_t> cplx_t;

/// Dynamic length multidimensional arrays
template <typename Numeric> using mdata_t = Eigen::Array<Numeric, Eigen::Dynamic, Eigen::Dynamic>;
typedef mdata_t<real_t> rmdata_t;
typedef mdata_t<cplx_t> cmdata_t;

/// Dynamic length arrays
template <typename Numeric> using data_t = Eigen::Array<Numeric, Eigen::Dynamic, 1>;
typedef data_t<real_t> rdata_t;
typedef data_t<cplx_t> cdata_t;

/// Matrices
template <typename Numeric, int N, int M> using matrix_t = Eigen::Matrix<Numeric, N, M>;
template <int N, int M> using rmatrix_t = matrix_t<real_t, N, M>;
typedef rmatrix_t<Eigen::Dynamic, Eigen::Dynamic> rmat_t;
template <int N, int M> using cmatrix_t = matrix_t<cplx_t, N, M>;
typedef cmatrix_t<Eigen::Dynamic, Eigen::Dynamic> cmat_t;

/// Vectors
template <typename Numeric, int N> using vector_t = Eigen::Matrix<Numeric, N, 1>;
template <int N> using rvector_t = vector_t<real_t, N>;
typedef rvector_t<Eigen::Dynamic> rvec_t;
template <int N> using cvector_t = vector_t<cplx_t, N>;
typedef cvector_t<Eigen::Dynamic> cvec_t;

template <typename Float>
inline bool isRegular(Float f) {
    return std::isnormal(f) || (f == 0.0f);
}

/// String formatting helper function
std::string format(const char* fmt, ...)
#ifdef __GNUC__
__attribute__((format(printf, 1, 2)))
#endif
;
#endif
