#pragma once
// isr_integrals_I0_I11.hpp
//
// Unified HPC-oriented C++ core implementing the ISR integral basis I0..I11
// (with I6 kept as the propagator integral and I10/I11 the newly added
// half-/negative-integer extensions).
//
// Integrals included:
//   I0  = ∫_0^b dx F(x,W)/(1-x)
//   I1  = ∫_0^b dx F(x,W)/(A+Bx)
//   I2  = ∫_0^b dx F(x,W)
//   I3  = ∫_0^b dx F(x,W)/(1-x)^3
//   I4  = ∫_0^b dx F(x,W)/(1-x)^2
//   I5  = ∫_0^b dx F(x,W)/(1-x)^(5/2)
//   I6  = ∫_0^b dx F(x,W) sqrt(1-x)/(A+Bx)
//   I7  = ∫_0^b dx F(x,W)/(1-x)^(1/2)
//   I8  = ∫_0^b dx F(x,W)/(1-x)^(3/2)
//   I9  = ∫_0^b dx F(x,W)/(1-x)^(-1/2) = ∫_0^b dx F(x,W) sqrt(1-x)
//   I10 = ∫_0^b dx F(x,W)/(1-x)^4
//   I11 = ∫_0^b dx F(x,W)/(1-x)^(7/2)
//
// Design goals:
//   - One-pass evaluation of all integrals that reuse the same scalar blocks.
//   - HPC-friendly: minimize repeated logs/dilogs/incomplete-beta evaluations.
//   - Keep the same numerical strategy already used in your existing cores:
//       * B(b;beta,0/-1/-2) by power series in b
//       * B(b;beta,1/2) by unregularized incomplete beta
//       * I1 via a connection formula for 2F1 at large argument
//       * I6 via the existing IA + IB implementation
//   - The new I10 and I11 are implemented in refined forms that avoid the raw
//     Mathematica LerchPhi / digamma / branch-heavy representation.
//
// Build/link:
//   - Requires libm
//   - For real dilog: GSL preferred; Boost dilog fallback supported
//   - For complex dilog (needed by I1/I6): GSL is required
//
// Optional convenience overload:
//   If variables.h is available and defines alpha and me, an overload using
//   (W, Mv, Gv, Wm) is enabled automatically.

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>

#if __has_include(<gsl/gsl_sf_dilog.h>)
  #include <gsl/gsl_sf_dilog.h>
  #include <gsl/gsl_sf_result.h>
  #define ISR011_HAVE_GSL_DILOG 1
#else
  #define ISR011_HAVE_GSL_DILOG 0
#endif

#if __has_include(<gsl/gsl_sf_beta.h>)
  #include <gsl/gsl_sf_beta.h>
  #define ISR011_HAVE_GSL_BETA 1
#elif __has_include(<gsl/gsl_sf_gamma.h>)
  #include <gsl/gsl_sf_gamma.h>
  #define ISR011_HAVE_GSL_BETA 1
#else
  #define ISR011_HAVE_GSL_BETA 0
#endif


#if __has_include(<variables.h>)
  #include <variables.h>
  #define ISR011_HAVE_VARIABLES 1
#elif __has_include("variables.h")
  #include "variables.h"
  #define ISR011_HAVE_VARIABLES 1
#else
  #define ISR011_HAVE_VARIABLES 0
#endif

namespace isr_i0i11 {

using cd = std::complex<double>;
#ifndef PI
  static constexpr double PI = 3.141592653589793238462643383279502884;
#endif
struct Options {
  // B(b;beta,0/-1/-2) series truncation
  int    max_n_beta   = 200000;
  double rel_tol_beta = 1e-15;

  // I1: 2F1 series truncation for 2F1(1,1-beta;2-beta;z)
  int    max_n_2f1    = 200000;
  double rel_tol_2f1  = 1e-15;

  // I6: I_romI outer/inner series truncation
  int    max_q_romI   = 200000;
  int    max_k_romI   = 200000;
  double rel_tol_q_romI = 1e-14;
  double rel_tol_k_romI = 1e-15;

  // I6: merged II+III q-sum truncation
  int    max_q_23     = 200000;
  double rel_tol_q_23 = 1e-14;

  // I6 safety: derivation is optimized for |y| > 1, y = -bB/A
  double y_min_abs    = 1.0;

  // If true, avoid throwing and return NaNs on failure.
  bool   fast_noexcept = false;
};

struct Result {
  // Real integrals
  double I0  = 0.0;
  double I2  = 0.0;
  double I3  = 0.0;
  double I4  = 0.0;
  double I5  = 0.0;
  double I7  = 0.0;
  double I8  = 0.0;
  double I9  = 0.0;
  double I10 = 0.0;
  double I11 = 0.0;

  // Complex integrals
  cd I1  = cd(0.0, 0.0);
  cd I6  = cd(0.0, 0.0);
  cd IA  = cd(0.0, 0.0);
  cd IB  = cd(0.0, 0.0);
  cd IromI = cd(0.0, 0.0);
  cd I23   = cd(0.0, 0.0);

  // Diagnostics
  int n_used_beta = 0;
  int n_used_2f1  = 0;
  int q_used_romI = 0;
  int k_used_romI_avg = 0;
  int q_used_23 = 0;

  bool conv_beta  = true;
  bool conv_2f1   = true;
  bool converged_romI = true;
  bool converged_23   = true;
};

struct BetaSeriesOut {
  double B0  = 0.0; // B(b;beta,0)
  double Bm1 = 0.0; // B(b;beta,-1)
  double Bm2 = 0.0; // B(b;beta,-2)
  double Bhalf = 0.0; // B(b;beta,1/2)
  int n_used = 0;
  bool conv = true;
};

struct SharedScalars {
  double b = 0.0;
  double beta = 0.0;
  double delta = 0.0;
  double cdelta = 0.0;

  double om  = 0.0; // 1-b
  double om2 = 0.0;
  double om3 = 0.0;
  double u   = 0.0; // sqrt(1-b)
  double u3  = 0.0;
  double u5  = 0.0;

  double inv1mb  = 0.0;
  double inv1mb2 = 0.0;
  double inv1mb3 = 0.0;
  double invu  = 0.0;
  double invu3 = 0.0;
  double invu5 = 0.0;

  double Lb   = 0.0;
  double L1m  = 0.0;
  double L4   = 0.0;
  double artU = 0.0;

  double bp = 0.0;       // b^beta
  double Li2b = 0.0;     // Li2(b)
  double Li2_1m = 0.0;   // Li2(1-b)
  double li2pm = 0.0;    // Li2(u)-Li2(-u)
  double PI2 = 0.0;
  double Bhalf = 0.0;     // B(b;beta,1/2)

};

// -----------------------------------------------------------------------------
// Special function helpers
// -----------------------------------------------------------------------------
inline double beta_complete(double a, double b) {
  return std::exp(std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b));
}

inline double reg_incomplete_beta(double x, double a, double b) {
#if ISR011_HAVE_GSL_BETA
  return gsl_sf_beta_inc(a, b, x);
#else
  (void)x; (void)a; (void)b;
  throw std::runtime_error("reg_incomplete_beta: GSL incomplete beta backend unavailable.");
#endif
}

inline double beta_incomplete_unreg(double x, double a, double b) {
#if ISR011_HAVE_GSL_BETA
  const double Ix = reg_incomplete_beta(x, a, b);
  return beta_complete(a, b) * Ix;
#else
  (void)x; (void)a; (void)b;
  throw std::runtime_error("beta_incomplete_unreg: GSL incomplete beta backend unavailable.");
#endif
}

inline double Li2_real_series_core(double x) {
  double sum = 0.0;
  double xk = x;
  for (int k = 1; k <= 200000; ++k) {
    const double term = xk / (double(k) * double(k));
    sum += term;
    if (k > 10 && std::abs(term) <= 1e-16 * std::max(1.0, std::abs(sum))) break;
    xk *= x;
  }
  return sum;
}

inline double Li2_real(double x) {
#if ISR011_HAVE_GSL_DILOG
  return gsl_sf_dilog(x);
#else
  // Portable fallback specialized for the real arguments used here: x in [-1,1].
  if (x == 1.0)  return PI * PI / 6.0;
  if (x == -1.0) return -PI * PI / 12.0;
  if (x == 0.0)  return 0.0;
  if (x < -0.5) {
    const double y = x / (x - 1.0);
    const double l = std::log1p(-x);
    return -Li2_real(y) - 0.5 * l * l;
  }
  if (x > 0.5) {
    return PI * PI / 6.0 - std::log(x) * std::log1p(-x) - Li2_real(1.0 - x);
  }
  return Li2_real_series_core(x);
#endif
}

inline cd Li2_complex(cd z) {

#if ISR011_HAVE_GSL_DILOG
  const double r = std::abs(z);
  const double theta = std::arg(z);
  gsl_sf_result re, im;
  gsl_sf_complex_dilog_e(r, theta, &re, &im);
  return cd(re.val, im.val);
#else
  (void)z;
  throw std::runtime_error("Li2_complex: GSL complex dilog is required for I1/I6.");
#endif
}

// -----------------------------------------------------------------------------
// B(b;beta,0/-1/-2) by series
// -----------------------------------------------------------------------------
inline BetaSeriesOut beta_series_0m1m2(double b, double beta, const Options& opt) {
  BetaSeriesOut out;

  const double Lb = std::log(b);
  double pow_bn = std::exp(beta * Lb); // b^(beta+n), starts at n=0

  double S0 = 0.0, Sm1 = 0.0, Sm2 = 0.0, Sh = 0.0;
  double ch = 1.0; // (1/2)_n / n!, starts at n=0
  bool ok = true;

  for (int n = 0; n < opt.max_n_beta; ++n) {
    const double dn = static_cast<double>(n);
    const double denom = beta + dn;

    const double t0  = pow_bn / denom;
    const double tm1 = (dn + 1.0) * pow_bn / denom;
    const double tm2 = 0.5 * (dn + 1.0) * (dn + 2.0) * pow_bn / denom;
    const double th  = ch * pow_bn / denom;

    S0  += t0;
    Sm1 += tm1;
    Sm2 += tm2;
    Sh  += th;

    out.n_used = n + 1;

    const double scale = std::max({1.0, std::abs(S0), std::abs(Sm1), std::abs(Sm2), std::abs(Sh)});
    const double tmax  = std::max({std::abs(t0), std::abs(tm1), std::abs(tm2), std::abs(th)});
    if (n > 10 && tmax <= opt.rel_tol_beta * scale) break;

    pow_bn *= b;
    ch *= (2.0 * dn + 1.0) / (2.0 * dn + 2.0); // c_{n+1} = c_n * (2n+1)/(2n+2)

    if (n + 1 >= opt.max_n_beta) { ok = false; break; }
  }

  out.B0 = S0;
  out.Bm1 = Sm1;
  out.Bm2 = Sm2;
#if ISR011_HAVE_GSL_BETA
  out.Bhalf = beta_incomplete_unreg(b, beta, 0.5);
#else
  out.Bhalf = Sh;
#endif
  out.conv = ok;
  return out;
}

// -----------------------------------------------------------------------------
// Shared scalar block
// -----------------------------------------------------------------------------
inline SharedScalars make_shared_scalars(double b, double beta, double delta, double Bhalf) {
  SharedScalars s;
  s.b = b;
  s.beta = beta;
  s.delta = delta;
  s.cdelta = 1.0 + delta;

  s.om  = 1.0 - b;
  s.om2 = s.om * s.om;
  s.om3 = s.om2 * s.om;
  s.u   = std::sqrt(s.om);
  s.u3  = s.u * s.om;
  s.u5  = s.u3 * s.om;

  s.Lb   = std::log(b);
  s.L1m  = std::log1p(-b);
  s.L4   = std::log(4.0);
  s.artU = std::atanh(s.u);
  s.bp   = std::exp(beta * s.Lb);
  s.PI2  = PI * PI;

  s.Li2b   = Li2_real(b);
  s.Li2_1m = Li2_real(1.0 - b);
  s.li2pm  = Li2_real(s.u) - Li2_real(-s.u);
  s.Bhalf  = Bhalf;

  const double eps = 1e-300;
  s.inv1mb  = 1.0 / std::max(s.om, eps);
  s.inv1mb2 = s.inv1mb * s.inv1mb;
  s.inv1mb3 = s.inv1mb2 * s.inv1mb;
  s.invu  = 1.0 / std::max(s.u,  eps);
  s.invu3 = 1.0 / std::max(s.u3, eps);
  s.invu5 = 1.0 / std::max(s.u5, eps);

  return s;
}

// -----------------------------------------------------------------------------
// I0, I2, I3, I4
// -----------------------------------------------------------------------------
inline void compute_I0_I2_I3_I4(const SharedScalars& s,
                                const BetaSeriesOut& Bs,
                                Result& out)
{
  const double beta  = s.beta;
  const double beta2 = beta * beta;
  const double b     = s.b;
  const double cdelta = s.cdelta;

  // I0
  out.I0 =
      beta * cdelta * Bs.B0
    - 0.5 * b * beta
    + s.L1m * (0.25 * beta2 + 0.5 * beta + 0.375 * beta2 * b)
    + (beta2 / 16.0) * s.L1m * s.L1m
    - 0.5 * beta2 * b * s.Lb
    - 0.5 * beta2 * (s.Li2_1m - s.PI2 / 6.0 - s.Li2b);

  // I2
  out.I2 =
      cdelta * s.bp
    + 0.5 * beta2 * s.Li2b
    + b * b * (beta / 4.0 + beta2 / 32.0)
    + b * (-beta - 5.0 * beta2 / 16.0)
    + s.L1m * (-9.0 * beta2 / 16.0 + 0.75 * beta2 * b - 3.0 * beta2 * b * b / 16.0)
    + s.Lb * (beta2 * (b * b / 4.0 - b));

  // I3
  const double tA = (beta / 2.0 + beta2 / 8.0) * (b * b / 2.0) * s.inv1mb2;
  const double tB = -0.5 * (beta + 0.75 * beta2) * (s.inv1mb2 - 1.0);
  const double tC = -(beta2 / 32.0) * ((b * (2.0 - b) + 2.0 * s.L1m) * s.inv1mb2);
  const double tD = -0.75 * beta2 * s.L1m;
  const double tE = (beta2 / 4.0) * (b * ((1.0 - b) + (b - 2.0) * s.Lb) * s.inv1mb2);
  const double tF = -(beta2 / 2.0) * (b * s.Lb * s.inv1mb);
  const double tG = 0.5 * beta2 * s.Li2b + 0.25 * beta2 * s.L1m * s.L1m;
  const double tH = -(beta2 / 8.0) * ((b + s.L1m) * s.inv1mb);
  out.I3 = beta * cdelta * Bs.Bm2 + tA + tB + tC + tD + tE + tF + tG + tH;

  // I4
  out.I4 =
      beta * cdelta * Bs.Bm1
    + 0.5 * beta2 * (s.Li2b - s.Li2_1m + s.PI2 / 6.0)
    + (b * s.inv1mb) * (-0.75 * beta2 - 0.5 * beta)
    + s.L1m * (0.5 * beta - 0.375 * beta2)
    - (b * s.Lb * beta2) * (0.5 * s.inv1mb)
    + (beta2 / 16.0) * s.L1m * s.L1m
    - (beta2 / 8.0) * (s.L1m * s.inv1mb);
}

// -----------------------------------------------------------------------------
// I5, I7, I8, I9
// -----------------------------------------------------------------------------
inline void compute_I5_I7_I8_I9(const SharedScalars& s, Result& out) {
  const double b = s.b;
  const double beta = s.beta;
  const double cdelta = s.cdelta;

  // I7
  {
    const double term =
        48.0 * (s.u - 1.0)
      - s.u * b * (12.0 + beta)
      + beta * (-4.0 + 4.0 * s.u + 9.0 * s.PI2)
      + 36.0 * cdelta * s.Bhalf
      + 12.0 * beta * s.artU * (3.0 * s.L1m + 8.0)
      - 48.0 * beta * s.L4
      + 3.0 * beta * s.u * (b - 4.0) * (3.0 * s.L1m - 4.0 * s.Lb)
      - 36.0 * beta * s.li2pm;
    out.I7 = (beta / 36.0) * term;
  }

  // I8
  {
    const double bracket =
        b * (-4.0 + beta)
      + 8.0 * cdelta * s.bp
      + beta * (
            8.0 * (s.u - 1.0)
          + s.u * s.PI2
          + (-4.0 + 3.0 * b) * s.L1m
          + 4.0 * s.u * s.artU * s.L1m
          - 4.0 * b * s.Lb
          - 4.0 * s.u * s.li2pm
        );
    out.I8 = beta * (1.0 - 2.0 * beta) * cdelta * s.Bhalf
           + (beta * 0.25) * s.invu * bracket;
  }

  // I9
  {
    const double pref0 = beta * cdelta / (1.0 + 2.0 * beta);
    const double lead  = pref0 * (2.0 * s.u * s.bp + s.Bhalf);
    const double poly_u3 = s.u3 * (480.0 - 180.0 * b + (232.0 - 27.0 * b) * beta);
    const double inner_beta =
        -1072.0 + 840.0 * s.u + 225.0 * s.PI2
      + 960.0 * s.artU
      - 480.0 * s.L4
      + 900.0 * (s.artU - s.u) * s.L1m
      + 15.0 * s.u3 * (-8.0 + 3.0 * b) * (3.0 * s.L1m - 4.0 * s.Lb)
      - 900.0 * s.li2pm;
    const double tail = (beta / 900.0) * (-480.0 + poly_u3 + beta * inner_beta);
    out.I9 = lead + tail;
  }

  // I5
  {
    const double R0 = -48.0 - 52.0 * beta + 9.0 * b * (4.0 + 3.0 * beta);
    const double R1 = -8.0 - 12.0 * s.L1m + (4.0 - 3.0 * b) * (16.0 + 9.0 * s.L1m - 12.0 * s.Lb);

    const double partA =
        beta * cdelta * (
            (2.0 / 3.0) * s.bp * s.invu3
          + (-3.0 + 2.0 * beta) * (
                -(2.0 / 3.0) * s.bp * s.invu
              + ((-1.0 + 2.0 * beta) / 3.0) * s.Bhalf
            )
        );

    const double partB =
        (beta / 36.0) * (
            48.0
          + (R0 + beta * R1) * s.invu3
          - beta * (24.0 + 36.0 * s.L1m) * s.invu
          + beta * (
                20.0 + 9.0 * s.PI2 + 48.0 * s.L4
              + 36.0 * s.artU * (s.L1m - 8.0 / 3.0)
              - 36.0 * s.li2pm
            )
        );

    out.I5 = partA + partB;
  }
}

// -----------------------------------------------------------------------------
// I10, I11 (new refined forms)
// -----------------------------------------------------------------------------
inline void compute_I10_I11(const SharedScalars& s,
                            const BetaSeriesOut& Bs,
                            Result& out)
{
  const double beta = s.beta;
  const double cdelta = s.cdelta;

  // I10 refined from the raw Mathematica form
  {
    const double partA = beta * cdelta * (
        s.bp / (3.0 * s.om3)
        + (-3.0 + beta) * (
            -s.bp / (6.0 * s.om2)
            + (-2.0 + beta) * (
                s.bp / (6.0 * s.om)
                - ((-1.0 + beta) / 6.0) * Bs.B0
            )
        )
    );

    const double tail_inner =
        -64.0 - 3.0 * s.om - 24.0 * s.om2 + 91.0 * s.om3
        + 6.0 * (
            (20.0 * s.om3 - 12.0 * s.om - 8.0) * s.Lb
            - (2.0 + 3.0 * s.om + 24.0 * s.om2 + 20.0 * s.om3) * s.L1m
            + 12.0 * s.om3 * s.L1m * s.L1m
          )
        + 144.0 * s.om3 * s.Li2b;

    const double partB = (beta / (288.0 * s.om3)) * (
        -48.0 - 72.0 * s.om + 120.0 * s.om3 + beta * tail_inner
    );

    out.I10 = partA + partB;
  }

  // I11 refined from the raw Mathematica form
  {
    const double partA = beta * cdelta * (
        (2.0 / 5.0) * s.bp * s.invu5
        + (-5.0 + 2.0 * beta) * (
            -(2.0 / 15.0) * s.bp * s.invu3
            + (-3.0 + 2.0 * beta) * (
                (2.0 / 15.0) * s.bp * s.invu
                - ((-1.0 + 2.0 * beta) / 15.0) * s.Bhalf
            )
        )
    );

    const double inner =
        1088.0 + 225.0 * s.PI2 + 480.0 * s.L4
        - 960.0 * s.artU
        + 900.0 * s.artU * s.L1m
        - 900.0 * s.li2pm
        - (840.0 + 900.0 * s.L1m) * s.invu
        - (5.0 + 75.0 * s.L1m + 300.0 * s.Lb) * s.invu3
        - (243.0 + 45.0 * s.L1m + 180.0 * s.Lb) * s.invu5;

    const double partB = (beta / 900.0) * (
        480.0 - 300.0 * s.invu3 - 180.0 * s.invu5 + beta * inner
    );

    out.I11 = partA + partB;
  }
}

// -----------------------------------------------------------------------------
// I1 support: 2F1(1,1-beta;2-beta;z) series and full I1
// -----------------------------------------------------------------------------
inline cd hyp2f1_1_1mb_2mb(cd z, double beta, const Options& opt,
                          int* n_used = nullptr, bool* conv = nullptr)
{
  const double bpar = 1.0 - beta;
  const double cpar = 2.0 - beta;

  cd sum(1.0, 0.0);
  cd term(1.0, 0.0);
  bool ok = true;
  int ncount = 0;

  for (int n = 0; n < opt.max_n_2f1; ++n) {
    const double num = bpar + static_cast<double>(n);
    const double den = cpar + static_cast<double>(n);
    term *= (num / den) * z;
    sum += term;
    ncount = n + 1;

    if (n > 10 && std::abs(term) <= opt.rel_tol_2f1 * std::max(1.0, std::abs(sum))) break;
    if (ncount >= opt.max_n_2f1) { ok = false; break; }
  }

  if (n_used) *n_used = ncount;
  if (conv) *conv = ok;
  return sum;
}

inline cd compute_I1(cd A, cd B, const SharedScalars& s,
                     const Options& opt,
                     int* n_used_2f1 = nullptr, bool* conv_2f1 = nullptr)
{
  const double beta = s.beta;
  const double beta2 = beta * beta;

  const cd AB  = A + B;
  const cd ABb = A + B * cd(s.b, 0.0);

  const cd L = std::log(ABb / A);
  const cd logB_AB = std::log(B / AB);
  const cd logB_A  = std::log(B / A);

  const cd K = (Li2_complex(ABb / AB) - Li2_complex(A / AB)) + logB_AB * L;
  const cd J = cd(-s.PI2 / 6.0, 0.0) + 0.5 * L * L + Li2_complex(A / ABb) - logB_A * L;

  const cd t = cd(s.b, 0.0) * B / A;
  const cd z_small = -1.0 / t;

  bool ok = true;
  int n2 = 0;
  const cd hyp_small = hyp2f1_1_1mb_2mb(z_small, beta, opt, &n2, &ok);
  if (n_used_2f1) *n_used_2f1 = n2;
  if (conv_2f1) *conv_2f1 = ok;

  const cd t_mbeta = std::pow(t, -beta);
  const double sPI = std::sin(PI * beta);
  const cd hyp =
      (beta / (beta - 1.0)) * (1.0 / t) * hyp_small
    + (beta * PI / sPI) * t_mbeta;

  const cd term_hyp = (s.cdelta / A) * cd(s.bp, 0.0) * hyp;

  const cd invB  = 1.0 / B;
  const cd invB2 = invB * invB;
  const cd invA  = 1.0 / A;

  const cd term1 = (beta2 / 8.0 + beta / 2.0) * (cd(s.b, 0.0) * invB - A * invB2 * L);
  const cd term2 = -(beta + 0.75 * beta2) * invB * L;
  const cd term3 = 0.75 * beta2 * (-invB) * K;
  const cd term4 = -(3.0 / 8.0) * beta2 * invB * (cd(-s.b, 0.0) + cd(s.b - 1.0, 0.0) * cd(s.L1m, 0.0) + (A * invB) * K);
  const cd term5 = -beta2 * invB * J;
  const cd term6 = 0.5 * beta2 * invB * (cd(s.b * (s.Lb - 1.0), 0.0) - (A * invB) * J);
  const cd term7 = -0.5 * beta2 * invA * (cd(-s.Li2b, 0.0) + K);

  return term_hyp + term1 + term2 + term3 + term4 + term5 + term6 + term7;
}

// -----------------------------------------------------------------------------
// I6 support (ported from the existing dedicated core)
// -----------------------------------------------------------------------------
inline cd compute_IA(cd A, cd B, double b, double beta) {
  const double u  = std::sqrt(1.0 - b);
  const double Du = 1.0 - u;
  const double Lb1 = std::log1p(-b);
  const double Lu1 = std::log1p(u);

  const cd D = A + B;
  const cd S = std::sqrt(D);
  const cd sB = std::sqrt(B);
  const cd B52 = B * B * sB;

  const cd r0 = sB / S;
  const cd r1 = (sB * u) / S;
  const cd At0 = std::atanh(r0);
  const cd At1 = std::atanh(r1);

  const cd z0 = B / D;
  const cd z1 = (B * (1.0 - b)) / D;

  const cd P = 2.0 * B * (4.0 + 3.0 * beta) + A * (4.0 + beta);
  const cd Q = 3.0 * A * A + 6.0 * A * B + 4.0 * B * B;

  const double Li_1mu = Li2_real(1.0 - u);
  const double Li_minu = Li2_real(-u);
  const cd Li_z0 = Li2_complex(z0);
  const cd Li_r0 = Li2_complex(r0);
  const cd Li_r1 = Li2_complex(r1);
  const cd Li_z1 = Li2_complex(z1);

  const cd CA =
      3.0 * A * A * sB * Du * (4.0 + 7.0 * beta)
    + A * B * sB * (u * b * (4.0 + 3.0 * beta) + Du * (28.0 + 57.0 * beta))
    + B52 * beta * (PI * PI)
    - 3.0 * (
        A * S * P * At0
      + beta * (A * sB * u * ((b - 7.0) * B - 3.0 * A) - 2.0 * B52 * Lu1) * Lb1
      + S * At1 * (-A * P + Q * beta * Lb1)
    );

  const cd CB =
      8.0 * B52 * cd(Li_1mu + Li_minu, 0.0)
    + S * Q * (Li_z0 - 4.0 * Li_r0 + 4.0 * Li_r1 - Li_z1);

  return (beta / (24.0 * A * B52)) * (2.0 * CA + 3.0 * beta * CB);
}

inline cd compute_IromI(cd A, cd y, double b, double beta,
                        const Options& opt,
                        int* q_used = nullptr, int* k_used_avg = nullptr, bool* converged = nullptr)
{
  const double Lb = std::log(b);
  const cd bpow = std::exp(cd(beta * Lb, 0.0));
  const cd invy = 1.0 / y;

  const double sPI = std::sin(PI * beta);
  const cd pow_neg_y = std::pow(-y, -beta);
  const cd sqrt_term = std::sqrt(1.0 - (b * invy));
  const cd term0 = (PI / sPI) * pow_neg_y * sqrt_term;

  cd sumq(0.0, 0.0);
  cd invyq(1.0, 0.0);

  int q_count = 0;
  long long k_total = 0;
  bool ok = true;

  for (int q = 0; q < opt.max_q_romI; ++q) {
    const double denom0 = beta - 1.0 - static_cast<double>(q);
    double coeffb = 1.0;
    double sumk = 0.0;
    double termk = coeffb / denom0;
    sumk += termk;

    int k_used = 1;
    for (int k = 0; k < opt.max_k_romI - 1; ++k) {
      const double kk = static_cast<double>(k);
      coeffb *= b * ((kk - 0.5) / (kk + 1.0));
      const double denom = denom0 + static_cast<double>(k + 1);
      termk = coeffb / denom;
      sumk += termk;
      ++k_used;

      if (std::abs(termk) <= opt.rel_tol_k_romI * std::max(1.0, std::abs(sumk))) break;
      if (k_used >= opt.max_k_romI) { ok = false; break; }
    }

    k_total += k_used;
    ++q_count;

    const cd qterm = invyq * cd(sumk, 0.0);
    sumq += qterm;
    if (q > 10 && std::abs(qterm) <= opt.rel_tol_q_romI * std::max(1.0, std::abs(sumq))) break;

    invyq *= invy;
    if (q_count >= opt.max_q_romI) { ok = false; break; }
  }

  if (q_used) *q_used = q_count;
  if (k_used_avg) *k_used_avg = (q_count > 0) ? int(double(k_total) / double(q_count) + 0.5) : 0;
  if (converged) *converged = ok;

  return (bpow / A) * (term0 - invy * sumq);
}

inline cd compute_I23(cd A, cd y, double b, double beta,
                      const Options& opt,
                      int* q_used = nullptr, bool* converged = nullptr)
{
  const cd invy = 1.0 / y;
  const cd by = cd(b, 0.0) * invy;
  const double logb = std::log(b);

  cd Hq = -std::log(1.0 - y);
  cd Sq = -Li2_complex(y);
  cd wq(1.0, 0.0);
  cd powy(1.0, 0.0);
  cd sum(0.0, 0.0);

  bool ok = true;
  int q_count = 0;

  for (int q = 0; q < opt.max_q_23; ++q) {
    const int qp1 = q + 1;
    const cd powy_next = powy * y;
    const cd Hnext = Hq - powy_next / double(qp1);
    const cd Snext = Sq + powy_next / double(qp1 * qp1);

    const cd bracket =
        -(cd(logb, 0.0) * Hq + Sq)
      + (cd(b, 0.0) * invy / 4.0) * (cd(logb + 0.5, 0.0) * Hnext + 2.0 * Snext);

    const cd qterm = wq * bracket;
    sum += qterm;
    ++q_count;

    if (q > 10 && std::abs(qterm) <= opt.rel_tol_q_23 * std::max(1.0, std::abs(sum))) break;

    Hq = Hnext;
    Sq = Snext;
    powy = powy_next;

    const double fac = (double(q) - 0.5) / double(qp1);
    wq *= cd(fac, 0.0) * by;

    if (q_count >= opt.max_q_23) { ok = false; break; }
  }

  if (q_used) *q_used = q_count;
  if (converged) *converged = ok;

  return (beta * beta) * cd(b, 0.0) * (1.0 / (A * y)) * sum;
}

inline void compute_I6(cd A, cd B, double b, double beta, double delta,
                       const Options& opt, Result& out)
{
  const cd y = -cd(b, 0.0) * B / A;
  if (std::abs(y) <= opt.y_min_abs && !opt.fast_noexcept) {
    // Soft warning only: derived representation is optimized for |y| > 1.
  }

  out.IA = compute_IA(A, B, b, beta);
  out.IromI = compute_IromI(A, y, b, beta, opt, &out.q_used_romI, &out.k_used_romI_avg, &out.converged_romI);
  out.I23   = compute_I23(A, y, b, beta, opt, &out.q_used_23, &out.converged_23);
  out.IB = beta * (1.0 + delta) * out.IromI + out.I23;
  out.I6 = out.IA + out.IB;
}

// -----------------------------------------------------------------------------
// Public API
// -----------------------------------------------------------------------------
inline Result compute_all(double b, double beta, double delta, cd A, cd B,
                          const Options& opt = Options{})
{
  Result out;

  auto fail = [&](const char* msg) -> Result {
    if (!opt.fast_noexcept) throw std::runtime_error(msg);
    out.I0 = out.I2 = out.I3 = out.I4 = out.I5 = out.I7 = out.I8 = out.I9 = out.I10 = out.I11 =
      std::numeric_limits<double>::quiet_NaN();
    out.I1 = out.I6 = out.IA = out.IB = out.IromI = out.I23 = cd(NAN, NAN);
    out.conv_beta = out.conv_2f1 = out.converged_romI = out.converged_23 = false;
    return out;
  };

  if (!(b > 0.0 && b < 1.0)) return fail("compute_all(I0..I11): require 0<b<1.");
  if (!(beta > 0.0))         return fail("compute_all(I0..I11): require beta>0.");

  const auto Bs = beta_series_0m1m2(b, beta, opt);
  out.n_used_beta = Bs.n_used;
  out.conv_beta   = Bs.conv;

  const SharedScalars s = make_shared_scalars(b, beta, delta, Bs.Bhalf);

  compute_I0_I2_I3_I4(s, Bs, out);
  compute_I5_I7_I8_I9(s, out);
  compute_I10_I11(s, Bs, out);

  int n2 = 0;
  bool ok2 = true;
  try {
    out.I1 = compute_I1(A, B, s, opt, &n2, &ok2);
    out.n_used_2f1 = n2;
    out.conv_2f1 = ok2;
  } catch (...) {
    if (!opt.fast_noexcept) throw;
    out.I1 = cd(NAN, NAN);
    out.n_used_2f1 = 0;
    out.conv_2f1 = false;
  }

  try {
    compute_I6(A, B, b, beta, delta, opt, out);
  } catch (...) {
    if (!opt.fast_noexcept) throw;
    out.I6 = out.IA = out.IB = out.IromI = out.I23 = cd(NAN, NAN);
    out.q_used_romI = out.k_used_romI_avg = out.q_used_23 = 0;
    out.converged_romI = out.converged_23 = false;
  }
  return out;
}

#if ISR011_HAVE_VARIABLES
inline void make_params(double W, double Mv, double Gv, double Wm,
                        double& b, double& beta, double& delta, cd& A, cd& B)
{
  const double x = (W / Mv) * (W / Mv);
  A = cd(Gv / Mv, x - 1.0);
  B = cd(0.0, -x);
  b = 1.0 - (Wm / W) * (Wm / W);
  beta = 2.0 * alpha / PI * (2.0 * std::log(W / me) - 1.0);
  delta = 0.75 * beta + alpha / PI * (PI * PI / 3.0 - 0.5)
        + beta * beta * (9.0 / 32.0 - PI * PI / 12.0);
}

inline Result compute_all(double W, double Mv, double Gv, double Wm,
                          const Options& opt = Options{})
{
  cd A, B;
  double b, beta, delta;
  make_params(W, Mv, Gv, Wm, b, beta, delta, A, B);
  return compute_all(b, beta, delta, A, B, opt);
}
#endif

// Convenience wrappers
inline double I0(double b, double beta, double delta, cd A, cd B, const Options& opt = Options{}) { return compute_all(b, beta, delta, A, B, opt).I0; }
inline cd     I1(double b, double beta, double delta, cd A, cd B, const Options& opt = Options{}) { return compute_all(b, beta, delta, A, B, opt).I1; }
inline double I2(double b, double beta, double delta, cd A, cd B, const Options& opt = Options{}) { return compute_all(b, beta, delta, A, B, opt).I2; }
inline double I3(double b, double beta, double delta, cd A, cd B, const Options& opt = Options{}) { return compute_all(b, beta, delta, A, B, opt).I3; }
inline double I4(double b, double beta, double delta, cd A, cd B, const Options& opt = Options{}) { return compute_all(b, beta, delta, A, B, opt).I4; }
inline double I5(double b, double beta, double delta, cd A, cd B, const Options& opt = Options{}) { return compute_all(b, beta, delta, A, B, opt).I5; }
inline cd     I6(double b, double beta, double delta, cd A, cd B, const Options& opt = Options{}) { return compute_all(b, beta, delta, A, B, opt).I6; }
inline double I7(double b, double beta, double delta, cd A, cd B, const Options& opt = Options{}) { return compute_all(b, beta, delta, A, B, opt).I7; }
inline double I8(double b, double beta, double delta, cd A, cd B, const Options& opt = Options{}) { return compute_all(b, beta, delta, A, B, opt).I8; }
inline double I9(double b, double beta, double delta, cd A, cd B, const Options& opt = Options{}) { return compute_all(b, beta, delta, A, B, opt).I9; }
inline double I10(double b, double beta, double delta, cd A, cd B, const Options& opt = Options{}) { return compute_all(b, beta, delta, A, B, opt).I10; }
inline double I11(double b, double beta, double delta, cd A, cd B, const Options& opt = Options{}) { return compute_all(b, beta, delta, A, B, opt).I11; }

} // namespace isr_i0i11

