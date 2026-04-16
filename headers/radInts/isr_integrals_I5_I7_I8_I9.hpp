#pragma once
// isr_integrals_I5_I7_I8_I9.hpp
//
// High-performance (HPC-friendly) core for four ISR basis integrals: I5, I7, I8, I9
// (excluding I6), matching the compact closed forms derived previously.
//
// Inputs:
//   - b in (0,1)
//   - beta > 0
//   - delta finite
//
// Special functions:
//   - Li2(x) (Spence dilog) for real x: prefer GSL gsl_sf_dilog(x) if available
//   - Incomplete beta: Mathematica Beta[z,a,b] = unregularized incomplete beta B_z(a,b)
//       B_z(a,b) = B(a,b) * I_z(a,b)  (I_z regularized incomplete beta)
//
// Backends:
//   - If GSL headers are present, use GSL for Li2 and (regularized) incomplete beta.
//   - If GSL beta header is missing on your system, we fall back to Boost's ibeta for I_z.
//     (Complete beta is computed via lgamma in both cases.)
//
// Notes:
//   - atanh(u) uses std::atanh (libm); this is typically as fast as any alternative.
//   - For best throughput, use compute_all() to get all four integrals in one pass.

#include <cmath>
#include <cfloat>
#include <stdexcept>

namespace isr_int {

// ---------- Optional GSL / Boost configuration ----------
#if __has_include(<gsl/gsl_sf_dilog.h>)
  #include <gsl/gsl_sf_dilog.h>
  #define ISR_HAVE_GSL_DILOG 1
#else
  #define ISR_HAVE_GSL_DILOG 0
#endif

// Some GSL installs expose gsl_sf_beta_inc in gsl_sf_beta.h, others via gsl_sf_gamma.h.
#if __has_include(<gsl/gsl_sf_beta.h>)
  #include <gsl/gsl_sf_beta.h>
  #define ISR_HAVE_GSL_BETA 1
#elif __has_include(<gsl/gsl_sf_gamma.h>)
  #include <gsl/gsl_sf_gamma.h>
  #define ISR_HAVE_GSL_BETA 1
#else
  #define ISR_HAVE_GSL_BETA 0
#endif

#if !ISR_HAVE_GSL_BETA
  // Fall back to Boost for regularized incomplete beta if GSL beta header is missing.
  #include <boost/math/special_functions/ibeta.hpp>
#endif

#if !ISR_HAVE_GSL_DILOG
  // Fall back to Boost dilog if GSL dilog header is missing.
  #if __has_include(<boost/math/special_functions/dilog.hpp>)
    #include <boost/math/special_functions/dilog.hpp>
    #define ISR_HAVE_BOOST_DILOG 1
  #else
    #define ISR_HAVE_BOOST_DILOG 0
  #endif
#endif

// ---------- Result container ----------
struct I5I7I8I9 {
  double I5 = 0.0;
  double I7 = 0.0;
  double I8 = 0.0;
  double I9 = 0.0;
};

// ---------- Helpers ----------
inline double beta_complete(double a, double b) {
  return std::exp(std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b));
}

inline double log1m(double b) {
  return std::log1p(-b);
}

// regularized incomplete beta I_x(a,b)
inline double reg_incomplete_beta(double x, double a, double b) {
#if ISR_HAVE_GSL_BETA
  return gsl_sf_beta_inc(a, b, x);
#else
  return boost::math::ibeta(a, b, x);
#endif
}

// Mathematica Beta[x,a,b] = unregularized incomplete beta B_x(a,b)
inline double beta_incomplete_unreg(double x, double a, double b) {
  const double Ix = reg_incomplete_beta(x, a, b);
  return beta_complete(a, b) * Ix;
}

// Li2(x) for real x
inline double Li2_real(double x) {
#if ISR_HAVE_GSL_DILOG
  return gsl_sf_dilog(x);
#else
  #if ISR_HAVE_BOOST_DILOG
    return boost::math::dilog(x);
  #else
    throw std::runtime_error("Li2_real: neither GSL dilog nor Boost dilog is available.");
  #endif
#endif
}

// ---------- Core: compute all four integrals in one go ----------
inline I5I7I8I9 compute_all(double b, double beta, double delta) {
  if (!(b > 0.0 && b < 1.0)) throw std::invalid_argument("compute_all: require 0 < b < 1.");
  if (!(beta > 0.0))        throw std::invalid_argument("compute_all: require beta > 0.");

  constexpr double PILOCAL  = 3.141592653589793238462643383279502884;
  const double PI2 = PILOCAL * PILOCAL;

  const double u   = std::sqrt(1.0 - b);      // u = sqrt(1-b)
  const double u3  = u * u * u;               // (1-b)^(3/2)
  const double L1m = log1m(b);                // log(1-b)
  const double Lb  = std::log(b);             // log(b)
  const double L4  = std::log(4.0);           // log(4)
  const double artU = std::atanh(u);          // atanh(u)
  const double bp  = std::exp(beta * Lb);     // b^beta
  const double cdelta = 1.0 + delta;          // 1+delta

  // li2pm = Li2(u) - Li2(-u)
  const double li2pm = Li2_real(u) - Li2_real(-u);

  // incBeta = Beta[b, beta, 1/2] in Mathematica (unregularized incomplete beta)
  const double incBeta = beta_incomplete_unreg(b, beta, 0.5);

  // Guard divisions (b extremely close to 1). Physical expressions can be singular; this avoids FP traps.
  const double u_eps = 1e-300;
  const double invu  = 1.0 / std::max(u,  u_eps);
  const double invu3 = 1.0 / std::max(u3, u_eps);

  I5I7I8I9 out;

  // I7
  {
    const double term =
      48.0*(u - 1.0)
      - u*b*(12.0 + beta)
      + beta*(-4.0 + 4.0*u + 9.0*PI2)
      + 36.0*cdelta*incBeta
      + 12.0*beta*artU*(3.0*L1m + 8.0)
      - 48.0*beta*L4
      + 3.0*beta*u*(b - 4.0)*(3.0*L1m - 4.0*Lb)
      - 36.0*beta*li2pm;

    out.I7 = (beta/36.0) * term;
  }

  // I8
  {
    const double bracket =
      b*(-4.0 + beta)
      + 8.0*cdelta*bp
      + beta * (
          8.0*(u - 1.0)
          + u*PI2
          + (-4.0 + 3.0*b)*L1m
          + 4.0*u*artU*L1m
          - 4.0*b*Lb
          - 4.0*u*li2pm
        );

    out.I8 =
      beta*(1.0 - 2.0*beta)*cdelta*incBeta
      + (beta * 0.25) * invu * bracket;
  }

  // I9
  {
    const double pref0 = beta*cdelta / (1.0 + 2.0*beta);
    const double lead  = pref0 * (2.0*u*bp + incBeta);

    const double poly_u3 = u3*(480.0 - 180.0*b + (232.0 - 27.0*b)*beta);

    const double inner_beta =
      -1072.0 + 840.0*u + 225.0*PI2
      + 960.0*artU
      - 480.0*L4
      + 900.0*(artU - u)*L1m
      + 15.0*u3*(-8.0 + 3.0*b)*(3.0*L1m - 4.0*Lb)
      - 900.0*li2pm;

    const double tail = (beta/900.0) * ( -480.0 + poly_u3 + beta*inner_beta );

    out.I9 = lead + tail;
    // I9 no longer needed
    // out.I9 = 0.0;
  }

  // I5
  {
    const double R0 = -48.0 - 52.0*beta + 9.0*b*(4.0 + 3.0*beta);
    const double R1 = -8.0 - 12.0*L1m + (4.0 - 3.0*b)*(16.0 + 9.0*L1m - 12.0*Lb);

    const double partA =
      beta*cdelta * (
        2.0*bp*(1.0/3.0)*invu3
        + (-3.0 + 2.0*beta) * (
            -2.0*bp*(1.0/3.0)*invu
            + ((-1.0 + 2.0*beta)/3.0) * incBeta
          )
      );

    const double partB =
      (beta/36.0) * (
        48.0
        + (R0 + beta*R1)*invu3
        - beta*(24.0 + 36.0*L1m)*invu
        + beta * (
            20.0 + 9.0*PI2 + 48.0*L4
            + 36.0*artU*(L1m - 8.0/3.0)
            - 36.0*li2pm
          )
      );

    out.I5 = partA + partB;
  }

  return out;
}

// ---------- Convenience wrappers (each calls compute_all) ----------
inline double I5(double b, double beta, double delta) { return compute_all(b,beta,delta).I5; }
inline double I7(double b, double beta, double delta) { return compute_all(b,beta,delta).I7; }
inline double I8(double b, double beta, double delta) { return compute_all(b,beta,delta).I8; }
inline double I9(double b, double beta, double delta) { return compute_all(b,beta,delta).I9; }

} // namespace isr_int