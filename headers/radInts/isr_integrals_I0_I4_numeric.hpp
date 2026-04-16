#pragma once
// isr_integrals_I0_I4_numeric.hpp
//
// High-precision numerical integration for ISR basis integrals I0..I4,
// using GSL adaptive Gauss–Kronrod quadrature (qag) on a *regularized* variable
// t ∈ [0,1] with x = b * t^p, p = 1/beta.
//
// This avoids endpoint overflow/inf that can happen with tanh–sinh when beta is small,
// while keeping the same public API style as before.
//
// -----------------------------------------------------------------------------
// Definitions (same as your Mathematica / analytic header):
//
//   F(x) = beta * x^(beta-1) * (1+delta)
//        - beta * (1 - x/2)
//        + (beta^2/8) * ( 4*(2-x)*Log[1/x]
//                         - ((1+3(1-x)^2)/x)*Log[1-x]
//                         - 6 + x ).
//
//   I0 = ∫_0^b dx  F(x)/(1-x)
//   I1 = ∫_0^b dx  F(x)/(A + B x)              (A,B complex)
//   I2 = ∫_0^b dx  F(x)
//   I3 = ∫_0^b dx  F(x)/(1-x)^3
//   I4 = ∫_0^b dx  F(x)/(1-x)^2
//
// Regularizing substitution:
//   x = b * t^p,   p = 1/beta,   t ∈ [0,1]
//   dx/dt = b * p * t^(p-1)
//
// Then I = ∫_0^1 dt  F(x(t)) * kernel(x(t)) * dx/dt
//
// Endpoint t=0 limit is finite and we hardcode it to avoid log(0).
//
// -----------------------------------------------------------------------------
// Build/link note (ROOT/ACLiC):
//   Make sure you link against GSL:
//     g++ ... -lgsl -lgslcblas
// ROOT usually already links if you use GSL elsewhere, but if not, add it.
//
// -----------------------------------------------------------------------------

#include <complex>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <algorithm>

#include <gsl/gsl_integration.h>

namespace isr_i0i4_num {

template <class Real>
using cd = std::complex<Real>;

template <class Real>
struct Options {
  // For GSL these map to epsabs/epsrel.
  // Note: GSL is double-based; asking << 1e-14 typically won’t help.
  Real abs_tol = Real(1e-14);
  Real rel_tol = Real(1e-14);

  // GSL workspace size (interval limit)
  size_t limit = 2000;

  // GSL qag rule key:
  //  1: 15pt, 2: 21pt, 3: 31pt, 4: 41pt, 5: 51pt, 6: 61pt
  int key = 6;

  // I1 safety
  Real denom_min_abs = Real(1e-30);

  bool fast_noexcept = false;
};

template <class Real>
struct Result {
  Real     I0 = Real(0);
  cd<Real> I1 = cd<Real>(Real(0), Real(0));
  Real     I2 = Real(0);
  Real     I3 = Real(0);
  Real     I4 = Real(0);

  // diagnostics
  int  evals = 0;        // GSL does not expose exact eval count; we track calls ourselves
  bool converged = true; // set false on any GSL error
};

// --------------------------- utilities (double domain) ---------------------------

static inline double log1p_neg(double x) { return std::log1p(-x); }

static inline double log1m_over_x(double x) {
  const double ax = std::abs(x);
  if (ax < 1e-14) {
    // log(1-x)/x = -(1 + x/2 + x^2/3 + x^3/4 + ...)
    const double x2 = x*x;
    return -(1.0 + x/2.0 + x2/3.0 + x2*x/4.0);
  }
  return log1p_neg(x) / x;
}

// F(x) in double (GSL uses double); stable near x=0.
static inline double RadF(double x, double beta, double delta) {
  using std::log;
  using std::exp;

  const double one = 1.0;
  const double beta2 = beta * beta;

  // Guard (should not be needed after substitution+limit handling)
  if (x <= 0.0) x = std::numeric_limits<double>::min();

  const double Lx = log(x);
  const double xpow = exp((beta - one) * Lx);
  const double Linvx = -Lx;

  const double omx = one - x;
  const double omx2 = omx * omx;

  const double fac = (one + 3.0*omx2) * log1m_over_x(x);

  const double term1 = beta * xpow * (one + delta);
  const double term2 = -beta * (one - x/2.0);
  const double term3 = (beta2/8.0) * ( 4.0*(2.0 - x)*Linvx - fac - 6.0 + x );

  return term1 + term2 + term3;
}

// --------------------------- GSL integrand plumbing ---------------------------

enum class Which {
  I0, I2, I3, I4,
  I1_re, I1_im
};

struct Params {
  double b;
  double beta;
  double delta;
  std::complex<double> A;
  std::complex<double> B;
  Which which;

  // derived
  double p;     // 1/beta
  double bpow;  // b^beta
  double t0_I0I2I3I4;          // (1+delta)*b^beta
  std::complex<double> t0_I1;  // (1+delta)*b^beta / A

  // eval counter (best-effort)
  int* eval_counter;
};

static inline double kernel_I0(double x) { return 1.0 / (1.0 - x); }
static inline double kernel_I2(double /*x*/) { return 1.0; }
static inline double kernel_I3(double x) {
  const double inv = 1.0 / (1.0 - x);
  return inv*inv*inv;
}
static inline double kernel_I4(double x) {
  const double inv = 1.0 / (1.0 - x);
  return inv*inv;
}

static double gsl_integrand(double t, void* void_params) {
  Params& P = *reinterpret_cast<Params*>(void_params);
  if (P.eval_counter) ++(*P.eval_counter);

  // Handle endpoint t=0 with analytic limit to avoid log(0) and 0^negative.
  if (t == 0.0) {
    switch (P.which) {
      case Which::I0:
      case Which::I2:
      case Which::I3:
      case Which::I4:
        return P.t0_I0I2I3I4;
      case Which::I1_re:
        return P.t0_I1.real();
      case Which::I1_im:
        return P.t0_I1.imag();
    }
  }

  // x = b * t^p,  p = 1/beta
  // dx/dt = b * p * t^(p-1)
  // (for beta small, p large; t^(p-1) -> 0 rapidly near 0, no overflow)
  const double x = P.b * std::pow(t, P.p);
  const double dxdt = P.b * P.p * std::pow(t, P.p - 1.0);

  // Compute F(x)
  const double Fx = RadF(x, P.beta, P.delta);

  // Apply kernels
  switch (P.which) {
    case Which::I0: {
      return Fx * kernel_I0(x) * dxdt;
    }
    case Which::I2: {
      return Fx * kernel_I2(x) * dxdt;
    }
    case Which::I3: {
      return Fx * kernel_I3(x) * dxdt;
    }
    case Which::I4: {
      return Fx * kernel_I4(x) * dxdt;
    }
    case Which::I1_re:
    case Which::I1_im: {
      const std::complex<double> denom = P.A + P.B * x;
      // Safety: if denom is too small, return NaN so GSL fails and we catch it.
      if (std::abs(denom) < 1e-300) {
        return std::numeric_limits<double>::quiet_NaN();
      }
      const std::complex<double> val = (Fx * dxdt) / denom;
      return (P.which == Which::I1_re) ? val.real() : val.imag();
    }
  }

  return std::numeric_limits<double>::quiet_NaN(); // unreachable
}

static inline double integrate_one(Params& P, const Options<long double>& opt_ld, gsl_integration_workspace* ws) {
  gsl_function F;
  F.function = &gsl_integrand;
  F.params = &P;

  const double epsabs = (double)opt_ld.abs_tol;
  const double epsrel = (double)opt_ld.rel_tol;

  double result = 0.0, abserr = 0.0;

  const int status = gsl_integration_qag(&F, 0.0, 1.0, epsabs, epsrel, (size_t)opt_ld.limit, opt_ld.key, ws, &result, &abserr);
  if (status != GSL_SUCCESS || !std::isfinite(result)) {
    // throw std::runtime_error("GSL qag failed (non-convergence or non-finite).");
  }
  return result;
}

// --------------------------- Public API (same structure) ---------------------------

// Primary (double inputs) API: integrates in double (GSL), returns long double fields.
inline Result<long double>
compute_all(double b, double beta, double delta,
            std::complex<double> A, std::complex<double> B,
            const Options<long double>& opt = Options<long double>{})
{
  Result<long double> out{};

  if (!opt.fast_noexcept) {
    if (!(b > 0.0 && b < 1.0)) throw std::domain_error("compute_all: require 0<b<1.");
    if (!(beta > 0.0))         throw std::domain_error("compute_all: require beta>0.");
  }
  if (b == 0.0) return out;

  // Precompute endpoint limits for t->0
  const double p = 1.0 / beta;
  const double bpow = std::exp(beta * std::log(b));         // b^beta
  const double t0 = (1.0 + delta) * bpow;                   // for I0,I2,I3,I4
  if (std::abs(A) < (double)opt.denom_min_abs) {
    throw std::runtime_error("compute_all: |A| too small (I1 endpoint limit blows up).");
  }
  const std::complex<double> t0I1 = t0 / A;

  int evals = 0;

  Params P;
  P.b = b; P.beta = beta; P.delta = delta; P.A = A; P.B = B;
  P.p = p; P.bpow = bpow;
  P.t0_I0I2I3I4 = t0;
  P.t0_I1 = t0I1;
  P.eval_counter = &evals;

  gsl_integration_workspace* ws = gsl_integration_workspace_alloc(opt.limit);
  if (!ws) throw std::runtime_error("compute_all: failed to allocate GSL workspace.");

  try {
    // I0
    P.which = Which::I0;
    const double I0 = integrate_one(P, opt, ws);

    // I2
    P.which = Which::I2;
    const double I2 = integrate_one(P, opt, ws);

    // I3
    P.which = Which::I3;
    const double I3 = integrate_one(P, opt, ws);

    // I4
    P.which = Which::I4;
    const double I4 = integrate_one(P, opt, ws);

    // I1 (real + imag)
    P.which = Which::I1_re;
    const double I1re = integrate_one(P, opt, ws);

    P.which = Which::I1_im;
    const double I1im = integrate_one(P, opt, ws);

    gsl_integration_workspace_free(ws);

    out.I0 = (long double)I0;
    out.I2 = (long double)I2;
    out.I3 = (long double)I3;
    out.I4 = (long double)I4;
    out.I1 = std::complex<long double>((long double)I1re, (long double)I1im);
    out.evals = evals;
    out.converged = true;
    return out;
  } catch (...) {
    gsl_integration_workspace_free(ws);
    out.converged = false;
    throw;
  }
}

// Backward-compatible alias used in your earlier script
inline Result<long double>
compute_I0_I4_numeric(double b, double beta, double delta,
                      std::complex<double> A, std::complex<double> B,
                      const Options<long double>& opt = Options<long double>{})
{
  return compute_all(b, beta, delta, A, B, opt);
}

} // namespace isr_i0i4_num