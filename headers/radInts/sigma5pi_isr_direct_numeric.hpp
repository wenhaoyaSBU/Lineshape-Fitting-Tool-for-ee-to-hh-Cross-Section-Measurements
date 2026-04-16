#pragma once
// sigma5pi_isr_direct_numeric.hpp
//
// Direct numerical ISR convolution (no I0..I4 basis) with optional:
//   (1) q_f^3 phase-space factor,
//   (2) vacuum polarization (VP) correction using Vacc(W)=1/|1-Π0(W)|^2,
//       under the “slowly-varying Vacc” approximation (PRD-style):
//       evaluate Vacc at nominal W and pull it outside the ISR integral.
//
// Definition implemented:
//
//   σ_ISR(W) = ∫_0^b dx  F(x;W) · σ_Born(W' = W*sqrt(1-x)),
//   b = 1 - (Wmin/W)^2.
//
// Born model is implemented in “piece” form (continuum/resonance/interference):
//   |1 + R|^2 = 1 + |R|^2 + 2 Re(R)
// so that we can apply VP as in Eq.(14):
//   σ_C -> σ_C * Vacc
//   σ_I -> σ_I * sqrt(Vacc)
//   σ_R unchanged
//
// Numerical integration:
//   - GSL adaptive Gauss–Kronrod: gsl_integration_qag over t ∈ [0,1]
//   - endpoint-regularizing substitution x = b * t^(1/beta), p=1/beta
//     (beta is KF radiator parameter at nominal W)
//
// Link with: -lgsl -lgslcblas
//

#include <complex>
#include <cmath>
#include <limits>
#include <stdexcept>

#include <gsl/gsl_integration.h>

#include "sigma5pi_isr.hpp" // isr_sigma5pi::{Consts,SigmaParts,beta_KF,delta_KF}

namespace isr_sigma5pi {

using cd = std::complex<double>;
static constexpr double PI_LOCAL = 3.141592653589793238462643383279502884;

// -------------------------- Options --------------------------

struct NIOptionsDirect {
  // GSL qag tolerances (double-based)
  double abs_tol = 1e-14;
  double rel_tol = 1e-14;

  // workspace size / recursion limit
  size_t limit = 2000;

  // qag rule key: 1=15pt ... 6=61pt
  int key = 6;

  // safety: if |A + B x| gets too small in I1-like denominators (not used here),
  // kept for uniformity; used below to guard BW denominator too.
  double denom_min_abs = 1e-30;

  // --- Optional corrections ---
  bool enable_qf3 = false; // multiply Born by q_f(W')^3
  bool enable_vp  = false; // apply VP correction using Vacc(W) (slowly-varying approx)

  // Default final-state masses for q_f: omega + pi0 [GeV]
  double m1 = mo;
  double m2 = mp;

  // Vacuum polarization input: Vacc(W) = 1/|1-Π0(W)|^2  (provided by your VP file).
  // If vacc_func is null, the constant Vacc member is used.
  // using VaccFunc = double(*)(double W, void* user);
  // VaccFunc vacc_func = nullptr;
  // void*    vacc_user = nullptr;
  // double   Vacc = 1.0; // used if vacc_func == nullptr
};

// -------------------------- Radiator F(x;W) --------------------------
// (Matches your RadF used in basis-integral work; KF form up to O(beta^2).)

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

static inline double KF_radiator(double x, double beta, double delta) {
  using std::log;
  using std::exp;

  if (x <= 0.0) x = std::numeric_limits<double>::min(); // avoid log(0)

  const double one = 1.0;
  const double beta2 = beta * beta;

  const double Lx    = log(x);
  const double Linvx = -Lx;

  // x^(beta-1)
  const double xpow = exp((beta - one) * Lx);

  const double omx  = one - x;
  const double omx2 = omx * omx;

  // ((1+3(1-x)^2)/x)*log(1-x)
  const double fac = (one + 3.0*omx2) * log1m_over_x(x);

  const double term1 = beta * xpow * (one + delta);
  const double term2 = -beta * (one - x/2.0);
  const double term3 = (beta2/8.0) * ( 4.0*(2.0 - x)*Linvx - fac - 6.0 + x );

  return term1 + term2 + term3;
}

// -------------------------- q_f helper (two-body) --------------------------

static inline double qf_two_body(double W, double m1, double m2) {
  if (W <= 0.0) return 0.0;
  const double s  = W*W;
  const double mp = m1 + m2;
  const double mm = std::abs(m1 - m2);
  const double lam = (s - mp*mp) * (s - mm*mm);
  if (lam <= 0.0) return 0.0;
  return std::sqrt(lam) / (2.0*W);
}

// -------------------------- Born pieces at W' --------------------------
// We compute:
//   σ_BC(W') = pref(W') * 1
//   σ_BR(W') = pref(W') * |R(W')|^2
//   σ_BI(W') = pref(W') * 2 Re R(W')
// where
//   pref(W') = (4π α^2)/(3 W'^2) * (Acal/W'^2)^2 * qf(W')^3
// and
//   R(W') = [3 W'^2 sqrt(Γee Γμμ) C1 e^{iφ}(1 + C2 e^{iΦ})] / [α M (W'^2 - M^2 + i M Γ)]
//
// NOTE: qf^3 is optional; when disabled qf^3=1.
// (VP is NOT applied here; it is applied once at nominal W after the ISR integral.)

static inline void sigma5pi_Born_parts(double Wp,
                                       double M,
                                       double Acal,
                                       double C1,
                                       double C2,
                                       double phi,
                                       double Phi,
                                       const Consts& c,
                                       const NIOptionsDirect& opt,
                                       double& sigBC,
                                       double& sigBR,
                                       double& sigBI)
{
  using std::sqrt;
  using std::cos;
  using std::sin;

  if (Wp <= 0.0) { sigBC = sigBR = sigBI = 0.0; return; }

  const double alpha = c.alpha;
  const double W2 = Wp*Wp;
  const double W4 = W2*W2;

  // optional qf^3
  double qf3 = 1.0;
  if (opt.enable_qf3) {
    const double qf = qf_two_body(Wp, opt.m1, opt.m2);
    qf3 = qf*qf*qf;
    if (qf3 == 0.0) { sigBC = sigBR = sigBI = 0.0; return; }
  }

  // prefactor: (4π α^2)/(3 W^2) * (A/W^2)^2 * qf^3
  const double prefC = (4.0*PI_LOCAL*alpha*alpha)/(3.0*W2) * (Acal*Acal)/(W4) * qf3;

  sigBC = prefC;

  // resonance ratio R(W)
  const cd i(0.0, 1.0);
  const cd phase = std::exp(i * phi);
  const cd strong = (1.0 + C2 * std::exp(i * Phi));

  const cd BW = cd(W2 - M*M, M*c.Gamma);
  const cd denom = alpha * M * BW;

  // guard (rare)
  if (std::abs(denom) < opt.denom_min_abs) {
    // treat as unstable point
    sigBC = sigBR = sigBI = std::numeric_limits<double>::quiet_NaN();
    return;
  }

  const double num_real = 3.0 * W2 * sqrt(c.Gamma_ee * c.Gamma_mumu) * C1;
  const cd R = (num_real * phase * strong) / denom;

  sigBR = prefC * std::norm(R);
  sigBI = prefC * (2.0 * std::real(R));
}

// -------------------------- GSL integrand plumbing --------------------------

enum class WhichPiece { BC, BR, BI };

struct DirectParams {
  double W, M;
  double Acal, C1, C2, phi, Phi;
  Consts c;

  // radiator at nominal W
  double beta, delta;
  double b;

  // substitution exponent p=1/beta
  double p;

  // endpoint limit constant for transformed integrand: b^beta * (1+delta)
  double bpow_1pdelta;

  WhichPiece which;

  // diagnostics
  int* eval_counter;

  // options pointer
  const NIOptionsDirect* opt;
};

static inline double integrand_t(double t, void* vp) {
  DirectParams& P = *reinterpret_cast<DirectParams*>(vp);
  if (P.eval_counter) ++(*P.eval_counter);

  // analytic t->0 limit (finite):
  // g(t) = F(x(t))*σBorn(W*sqrt(1-x(t)))*dx/dt
  // -> (1+delta)*b^beta * σBorn(W) as t->0.
  auto t0_limit = [&]() -> double {
    double bc, br, bi;
    sigma5pi_Born_parts(P.W, P.M, P.Acal, P.C1, P.C2, P.phi, P.Phi, P.c, *P.opt, bc, br, bi);
    const double sig = (P.which == WhichPiece::BC) ? bc : (P.which == WhichPiece::BR) ? br : bi;
    return P.bpow_1pdelta * sig;
  };

  if (t <= 0.0) return t0_limit();

  // x = b * t^p  with p = 1/beta
  // compute t^p via exp(p log t) to control underflow
  const double logt = std::log(t);
  const double y = P.p * logt;

  // if t^p underflows to 0, we are effectively at the limit
  if (y < std::log(std::numeric_limits<double>::min()) - 10.0) {
    return t0_limit();
  }

  const double tp = std::exp(y);
  const double x  = P.b * tp;

  // dx/dt = b * p * t^(p-1) = b*p*exp((p-1)log t)
  const double y2 = (P.p - 1.0) * logt;
  const double dxdt = P.b * P.p * std::exp(y2);

  // W' = W sqrt(1-x)
  const double Wp = P.W * std::sqrt(std::max(0.0, 1.0 - x));

  double bc, br, bi;
  sigma5pi_Born_parts(Wp, P.M, P.Acal, P.C1, P.C2, P.phi, P.Phi, P.c, *P.opt, bc, br, bi);

  const double sig = (P.which == WhichPiece::BC) ? bc : (P.which == WhichPiece::BR) ? br : bi;

  if (!std::isfinite(sig)) return std::numeric_limits<double>::quiet_NaN();

  const double Fx = KF_radiator(x, P.beta, P.delta);
  const double val = Fx * sig * dxdt;

  return val;
}

static inline double integrate_piece(DirectParams& P, const NIOptionsDirect& opt, gsl_integration_workspace* ws) {
  gsl_function F;
  F.function = &integrand_t;
  F.params   = &P;

  double result = 0.0, abserr = 0.0;
  const int status = gsl_integration_qag(&F, 0.0, 1.0,
                                         opt.abs_tol, opt.rel_tol,
                                         opt.limit, opt.key, ws,
                                         &result, &abserr);
  if (status != GSL_SUCCESS || !std::isfinite(result)) {
    throw std::runtime_error("sigma5pi_ISR_direct_num: GSL qag failed (non-convergence or non-finite).");
  }
  return result;
}

// -------------------------- Public API --------------------------
// Direct numerical ISR convolution, with optional qf^3 and VP corrections.

inline SigmaParts sigma5pi_ISR_direct_num(double W,
                                          double M,
                                          double Acal,
                                          double C1,
                                          double C2,
                                          double phi,
                                          double Phi,
                                          const Consts& c = Consts{},
                                          const NIOptionsDirect& opt = NIOptionsDirect{},
                                          const VPOptions& vp = VPOptions{})
{
  SigmaParts out;

  if (!(W > 0.0) || !(M > 0.0)) {
    throw std::invalid_argument("sigma5pi_ISR_direct_num: require W>0 and M>0.");
  }
  if (!(c.Wmin > 0.0 && c.Wmin < W)) {
    throw std::invalid_argument("sigma5pi_ISR_direct_num: require 0 < Wmin < W.");
  }
  if (!(c.Gamma > 0.0 && c.Gamma_ee > 0.0 && c.Gamma_mumu > 0.0)) {
    throw std::invalid_argument("sigma5pi_ISR_direct_num: require positive widths Gamma, Gamma_ee, Gamma_mumu.");
  }

  // KF radiator parameters at nominal W
  const double beta  = beta_KF(W, c);
  const double delta = delta_KF(beta, c);

  // upper limit b = 1 - (Wmin/W)^2
  const double b = 1.0 - (c.Wmin/W)*(c.Wmin/W);

  // substitution exponent p = 1/beta
  const double p = 1.0 / beta;

  // endpoint constant b^beta*(1+delta)
  const double bpow = std::exp(beta * std::log(b));
  const double bpow_1pdelta = bpow * (1.0 + delta);

  int evals = 0;

  DirectParams P;
  P.W = W; P.M = M;
  P.Acal = Acal; P.C1 = C1; P.C2 = C2; P.phi = phi; P.Phi = Phi;
  P.c = c;
  P.beta = beta; P.delta = delta;
  P.b = b;
  P.p = p;
  P.bpow_1pdelta = bpow_1pdelta;
  P.eval_counter = &evals;
  P.opt = &opt;

  gsl_integration_workspace* ws = gsl_integration_workspace_alloc(opt.limit);
  if (!ws) throw std::runtime_error("sigma5pi_ISR_direct_num: failed to allocate GSL workspace.");

  double sigC = 0.0, sigR = 0.0, sigI = 0.0;

  try {
    P.which = WhichPiece::BC;
    sigC = integrate_piece(P, opt, ws);

    P.which = WhichPiece::BR;
    sigR = integrate_piece(P, opt, ws);

    P.which = WhichPiece::BI;
    sigI = integrate_piece(P, opt, ws);

    gsl_integration_workspace_free(ws);
  } catch (...) {
    gsl_integration_workspace_free(ws);
    throw;
  }

  // ---- Vacuum polarization correction (slowly-varying approx) ----
  // Eq.(14): σ_BC/|1-Π0|^2 + σ_BR + σ_BI/|1-Π0|
  // Let Vacc = 1/|1-Π0|^2, so 1/|1-Π0| = sqrt(Vacc).
  double sigC_vp = sigC;
  double sigR_vp = sigR;
  double sigI_vp = sigI;

  if (vp.enable_vp) {
    const double VaccW = (vp.vacc_func) ? vp.vacc_func(W, vp.vacc_user) : vp.Vacc;
    if (!(VaccW > 0.0) || !std::isfinite(VaccW)) {
      throw std::runtime_error("sigma5pi_ISR_direct_num: invalid Vacc(W) (must be finite and >0).");
    }
    const double sqrtVacc = std::sqrt(VaccW);
    sigC_vp *= VaccW;
    sigI_vp *= sqrtVacc;
    // sigR unchanged
    (void)sqrtVacc;
  }

  out.sigmaC = sigC_vp * c.unit_conv;
  out.sigmaR = sigR_vp * c.unit_conv;
  out.sigmaI = sigI_vp * c.unit_conv;
  out.sigma  = (sigC_vp + sigR_vp + sigI_vp) * c.unit_conv;

  out.n_beta = evals;     // reuse to store #integrand calls
  out.n_2f1  = 0;
  out.conv_beta = true;
  out.conv_2f1  = true;

  return out;
}

} // namespace isr_sigma5pi