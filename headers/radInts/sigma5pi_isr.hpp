#pragma once
// sigma5pi_isr.hpp
//
// ISR-corrected cross section for e+e- -> 5π following Liu et al. (PRD 110, 053010, 2024), Appendix B,
// using analytic ISR integrals I0..I4.
//
// This is the “0th-order” core in Wenhao Ya's note (i.e. σ_ISR_noqf), where the q_f^3 factor
// is NOT included.
//
// Patch: add optional vacuum polarization (VP) correction using Vacc(W) = 1/|1-Π0(W)|^2,
// with 1/|1-Π0| = sqrt(Vacc). Under the “slowly-varying Vacc” approximation, we evaluate Vacc at
// nominal W and apply:
//   σ_C -> σ_C * Vacc
//   σ_I -> σ_I * sqrt(Vacc)
//   σ_R unchanged.
//
// Link with GSL: -lgsl -lgslcblas

#include <complex>
#include <cmath>
#include <stdexcept>
#include <algorithm>

#include "isr_integrals_I0_I4.hpp"
#include "isr_integrals_I0_I4_numeric.hpp"
#include "isr_integrals_I5_I7_I8_I9.hpp"   // I5,I7,I8,I9 (namespace isr_int)
#include "isr_integral_I6.hpp"             // I6 (namespace isr6)

namespace isr_sigma5pi {

using cd = std::complex<double>;

struct Consts {
  // QED / radiator
  double alpha = 1.0/137.035999084; // fine structure constant
  double me    = 0.511e-03;         // electron mass [GeV]

  // ISR lower limit Wmin
  double Wmin  = 2.8;               // [GeV]

  // Resonance widths [GeV]
  double Gamma      = 92.6e-6;      // total width Γ
  double Gamma_ee   = 5.55e-6;      // Γee
  double Gamma_mumu = 5.55e-6;      // Γμμ (often set equal to Γee in your tests)

  // GeV^-2 -> nb conversion etc (provided by your physicsFuncs.h)
  double unit_conv = hbar_c_2;
};


struct IsCache {
  // Switch
  bool useCache = false;

  // cache for I0..I4 (shared across orders 1..3)
  isr_i0i4::Result I04;

  // cache for I5,I7,I8,I9 (shared across orders 1..3)
  isr_int::I5I7I8I9 I5789;

  // cache for I6 (shared across orders 1..3)
  isr6::Result I6res;
};

struct SigmaParts {
  double sigma   = 0.0;
  double sigmaC  = 0.0;
  double sigmaR  = 0.0;
  double sigmaI  = 0.0;

  // diagnostics from integral core(s)
  int  n_beta = 0;
  int  n_2f1  = 0;
  bool conv_beta = true;
  bool conv_2f1  = true;
};

// Optional vacuum polarization input (Vacc = 1/|1-Π0|^2)
struct VPOptions {
  bool enable_vp = false;

  using VaccFunc = double(*)(double W, void* user);
  VaccFunc vacc_func = nullptr;
  void*    vacc_user = nullptr;

  // used if vacc_func==nullptr
  double Vacc = 1.0;
};

inline double get_Vacc(double W, const VPOptions& vp) {
  if (!vp.enable_vp) return 1.0;
  const double Vacc = vp.vacc_func ? vp.vacc_func(W, vp.vacc_user) : vp.Vacc;
  if (!(Vacc > 0.0) || !std::isfinite(Vacc)) {
    throw std::runtime_error("VPOptions: invalid Vacc(W). Must be finite and > 0.");
  }
  return Vacc;
}

inline void apply_vp(double W, double& sigmaC, double& sigmaR, double& sigmaI, const VPOptions& vp) {
  if (!vp.enable_vp) return;
  const double Vacc = get_Vacc(W, vp);
  const double sqrtVacc = std::sqrt(Vacc);
  sigmaC *= Vacc;
  sigmaI *= sqrtVacc;
  (void)sigmaR; // unchanged
}

// Kuraev–Fadin radiator parameters
inline double beta_KF(double W, const Consts& c) {
  // β = 2α/π (2 ln(W/me) - 1)
  return 2.0*c.alpha/PI * (2.0*std::log(W/c.me) - 1.0);
}

inline double delta_KF(double beta, const Consts& c) {
  // δ = 3/4 β + α/π (π^2/3 - 1/2) + β^2 (9/32 - π^2/12)
  const double PI2 = PI*PI;
  return 0.75*beta + c.alpha/PI*(PI2/3.0 - 0.5) + beta*beta*(9.0/32.0 - PI2/12.0);
}

// Paper's complex propagator parameters A,B
inline void make_AB_propagator(double W, double M, const Consts& c, cd& Aprop, cd& Bprop) {
  const double x = (W/M)*(W/M);          // (W/M)^2
  Aprop = cd(c.Gamma/M, x - 1.0);        // Γ/M + i[(W/M)^2 - 1]
  Bprop = cd(0.0, -x);                   // -i (W/M)^2
}

// 0th-order ISR cross section (no qf^3), optional VP
inline SigmaParts sigma5pi_ISR(double W,
                               double M,
                               double Acal,
                               double C1,
                               double C2,
                               double phi,
                               double Phi,
                               const Consts& c = Consts{},
                               const isr_i0i4::Options& iopt = isr_i0i4::Options{},
                               const VPOptions& vp = VPOptions{},
                               const IsCache& cache = IsCache{})
{
  SigmaParts out;

  if (!(W > 0.0) || !(M > 0.0)) {
    throw std::invalid_argument("sigma5pi_ISR: require W>0 and M>0.");
  }
  if (!(c.Wmin > 0.0 && c.Wmin < W)) {
    throw std::invalid_argument("sigma5pi_ISR: require 0 < Wmin < W.");
  }
  if (!(c.Gamma > 0.0 && c.Gamma_ee > 0.0 && c.Gamma_mumu > 0.0)) {
    throw std::invalid_argument("sigma5pi_ISR: require positive widths Gamma, Gamma_ee, Gamma_mumu.");
  }

  // paper's complex A,B
  cd Aprop, Bprop;
  make_AB_propagator(W, M, c, Aprop, Bprop);

  isr_i0i4::Result I; // declare here to assign from cache or compute
  if(!cache.useCache) {
  // radiator params
  const double beta  = beta_KF(W, c);
  const double delta = delta_KF(beta, c);

  // b = 1 - (Wmin/W)^2
  const double b = 1.0 - (c.Wmin/W)*(c.Wmin/W);

  // integrals I0..I4
  I = isr_i0i4::compute_all(b, beta, delta, Aprop, Bprop, iopt);
  }
  else {
    I = cache.I04;
  }

  out.n_beta = I.n_used_beta;
  out.n_2f1  = I.n_used_2f1;
  out.conv_beta = I.conv_beta;
  out.conv_2f1  = I.conv_2f1;

  const double I0 = I.I0;
  const double I3 = I.I3;
  const double I4 = I.I4;
  const cd     I1 = I.I1;

  // Prefactors
  const double alpha = c.alpha;
  const double A2 = Acal * Acal;
  const double W2 = W*W;
  const double W4 = W2*W2;
  const double W6 = W4*W2;
  const double M2 = M*M;
  const double M3 = M2*M;
  const double M5 = M3*M2;

  const cd AB = Aprop + Bprop;

  // H1 = cosφ + C2 cos(φ+Φ),  H2 = sinφ + C2 sin(φ+Φ)
  const double H1 = std::cos(phi) + C2 * std::cos(phi + Phi);
  const double H2 = std::sin(phi) + C2 * std::sin(phi + Phi);

  // σ_C (B4)
  out.sigmaC = (4.0*PI*alpha*alpha*A2 / (3.0*W6)) * I3;

  // resonance bracket Rbr = 1/(A+B) I0 + (B/(A+B)) I1
  const cd Rbr = (cd(I0,0.0)/AB) + (Bprop/AB)*I1;

  // σ_R (B5)
  {
    const double fac = (C1*C1) * (1.0 + C2*C2 + 2.0*C2*std::cos(Phi));
    const double prefR = fac * (6.0*PI * A2 * c.Gamma_ee * c.Gamma_mumu) / (c.Gamma * M5 * W2);
    out.sigmaR = prefR * (2.0 * std::real(Rbr));
  }

  // σ_I (B6)
  {
    const double prefI = (4.0*PI*alpha * A2 * std::sqrt(c.Gamma_ee*c.Gamma_mumu)) / M;

    const double termA = C1 * ( H2/(W4*M2) - H1/(c.Gamma*W4*M) );

    const cd termB =
        (Bprop/(AB*AB)) * cd(I0,0.0)
      + (Bprop/AB)*(Bprop/AB) * I1
      + (cd(I4,0.0)/AB);

    const double termC = C1 * H1 / (c.Gamma * W2 * M3);

    const cd inside = cd(termA,0.0)*termB + cd(termC,0.0)*Rbr;

    out.sigmaI = prefI * (2.0 * std::real(inside));
  }

  // optional VP correction (slowly varying Vacc at nominal W)
  apply_vp(W, out.sigmaC, out.sigmaR, out.sigmaI, vp);

  // units
  out.sigma  = (out.sigmaC + out.sigmaR + out.sigmaI) * c.unit_conv;
  out.sigmaC *= c.unit_conv;
  out.sigmaR *= c.unit_conv;
  out.sigmaI *= c.unit_conv;

  return out;
}

} // namespace isr_sigma5pi