#pragma once
// sigma2body_isr.hpp
//
// 0th-order ISR-corrected cross section for a two-body channel with the
// modified Born-level prefactor
//
//   sigma_Born,noqf(W) ~ 4 pi alpha^2 Acal^2 / (3 W^8) * |1 + R(W)|^2,
//
// i.e. compared with the previous W^{-6} implementation, there is one extra
// factor 1/W'^2 after the ISR substitution W' = W sqrt(1-x). Equivalently,
// every ISR kernel in the old formulas acquires one extra factor 1/(1-x).
//
// The resonance-amplitude structure is kept identical to the existing notes:
//
//   R(W) = 3 W^2 Gamma_ee C1 exp(i phi) (1 + C2 exp(i Phi))
//          / [ M alpha (W^2 - M^2 + i M Gamma) ]
//
// This header performs only the 0th-order ("no q_f^3") calculation.
// The q_f^3 polynomial corrections up to 3rd order are implemented in
// sigma2body_isr_qf3.hpp.
//
// Integral mapping relative to the old W^{-6} formulas:
//   I3  -> I10
//   J1  = ∫ F /[(1-x)    (A+Bx)] -> J2 = ∫ F /[(1-x)^2  (A+Bx)]
//   J2  = ∫ F /[(1-x)^2  (A+Bx)] -> J3 = ∫ F /[(1-x)^3  (A+Bx)]
// where
//   J1 = I0 /(A+B) + B/(A+B) I1,
//   J2 = I4 /(A+B) + B/(A+B) J1,
//   J3 = I3 /(A+B) + B/(A+B) J2.

#include <complex>
#include <cmath>
#include <stdexcept>
#include <algorithm>

#include "isr_integrals_I0_I11.hpp"

#include "sigma5pi_isr.hpp" // for VPOptions

namespace isr_sigma2body {

using cd = std::complex<double>;
#ifndef PI
  static constexpr double PI = 3.141592653589793238462643383279502884;
#endif

struct Consts {
  // QED / radiator
  double alpha = 1.0 / 137.035999084;
  double me    = 0.511e-03;   // GeV

  // ISR lower limit Wmin
  double Wmin  = 2.8;         // GeV

  // Resonance widths [GeV]
  double Gamma      = 92.6e-6;
  double Gamma_ee   = 5.55e-6;
  double Gamma_mumu = 5.55e-6;

  // GeV^{-2} -> user units. Set to hbar_c_2 externally if desired.
  // double unit_conv = hbar_c_2;
  // W^-8: * 1/8
  double unit_conv = hbar_c_2/8.;
};

struct IsCache {
  bool useCache = false;
  isr_i0i11::Result Iall;
};

struct SigmaParts {
  double sigma  = 0.0;
  double sigmaC = 0.0;
  double sigmaR = 0.0;
  double sigmaI = 0.0;

  int  n_beta = 0;
  int  n_2f1  = 0;
  bool conv_beta = true;
  bool conv_2f1  = true;
};

// struct VPOptions {
//   bool enable_vp = false;

//   using VaccFunc = double(*)(double W, void* user);
//   VaccFunc vacc_func = nullptr;
//   void*    vacc_user = nullptr;

//   double Vacc = 1.0;
// };

// inline double get_Vacc(double W, const VPOptions& vp) {
//   if (!vp.enable_vp) return 1.0;
//   const double Vacc = vp.vacc_func ? vp.vacc_func(W, vp.vacc_user) : vp.Vacc;
//   if (!(Vacc > 0.0) || !std::isfinite(Vacc)) {
//     throw std::runtime_error("VPOptions: invalid Vacc(W). Must be finite and > 0.");
//   }
//   return Vacc;
// }

// inline void apply_vp(double W, double& sigmaC, double& sigmaR, double& sigmaI, const VPOptions& vp) {
//   if (!vp.enable_vp) return;
//   const double Vacc = get_Vacc(W, vp);
//   const double sqrtVacc = std::sqrt(Vacc);
//   sigmaC *= Vacc;
//   sigmaI *= sqrtVacc;
//   (void)W;
//   (void)sigmaR; // unchanged
// }

using VPOptions = isr_sigma5pi::VPOptions;
using isr_sigma5pi::get_Vacc;
using isr_sigma5pi::apply_vp;

inline double beta_KF(double W, const Consts& c) {
  return 2.0 * c.alpha / PI * (2.0 * std::log(W / c.me) - 1.0);
}

inline double delta_KF(double beta, const Consts& c) {
  const double PI2 = PI * PI;
  return 0.75 * beta
       + c.alpha / PI * (PI2 / 3.0 - 0.5)
       + beta * beta * (9.0 / 32.0 - PI2 / 12.0);
}

inline void make_AB_propagator(double W, double M, const Consts& c, cd& Aprop, cd& Bprop) {
  const double x = (W / M) * (W / M);
  Aprop = cd(c.Gamma / M, x - 1.0);
  Bprop = cd(0.0, -x);
}

struct ShiftedKernels {
  cd J1 = cd(0.0, 0.0); // ∫ F /[(1-x)   (A+Bx)]
  cd J2 = cd(0.0, 0.0); // ∫ F /[(1-x)^2 (A+Bx)]
  cd J3 = cd(0.0, 0.0); // ∫ F /[(1-x)^3 (A+Bx)]
};

inline ShiftedKernels build_shifted_kernels(const isr_i0i11::Result& I, cd Aprop, cd Bprop) {
  ShiftedKernels K;
  const cd S  = Aprop + Bprop;
  const cd BS = Bprop / S;

  K.J1 = cd(I.I0, 0.0) / S + BS * I.I1;
  K.J2 = cd(I.I4, 0.0) / S + BS * K.J1;
  K.J3 = cd(I.I3, 0.0) / S + BS * K.J2;
  return K;
}

inline SigmaParts sigma2body_ISR(double W,
                                 double M,
                                 double Acal,
                                 double C1,
                                 double C2,
                                 double phi,
                                 double Phi,
                                 const Consts& c = Consts{},
                                 const isr_i0i11::Options& iopt = isr_i0i11::Options{},
                                 const VPOptions& vp = VPOptions{},
                                 const IsCache& cache = IsCache{})
{
  SigmaParts out;

  if (!(W > 0.0) || !(M > 0.0)) {
    throw std::invalid_argument("sigma2body_ISR: require W>0 and M>0.");
  }
  if (!(c.Wmin > 0.0 && c.Wmin < W)) {
    throw std::invalid_argument("sigma2body_ISR: require 0 < Wmin < W.");
  }
  if (!(c.Gamma > 0.0 && c.Gamma_ee > 0.0 && c.Gamma_mumu > 0.0)) {
    throw std::invalid_argument("sigma2body_ISR: require positive widths Gamma, Gamma_ee, Gamma_mumu.");
  }

  cd Aprop, Bprop;
  make_AB_propagator(W, M, c, Aprop, Bprop);

  isr_i0i11::Result I;
  if (!cache.useCache) {
    const double beta  = beta_KF(W, c);
    const double delta = delta_KF(beta, c);
    const double b = 1.0 - (c.Wmin / W) * (c.Wmin / W);
    I = isr_i0i11::compute_all(b, beta, delta, Aprop, Bprop, iopt);
  } else {
    I = cache.Iall;
  }

  // Only I1 is used from the complex/2F1 sector at 0th order.
  out.n_beta = I.n_used_beta;
  out.n_2f1  = I.n_used_2f1;
  out.conv_beta = I.conv_beta;
  out.conv_2f1  = I.conv_2f1;

  const auto K = build_shifted_kernels(I, Aprop, Bprop);

  const double alpha = c.alpha;
  const double A2 = Acal * Acal;
  const double W2 = W * W;
  const double W4 = W2 * W2;
  const double W8 = W4 * W4;
  const double M2 = M * M;
  const double M3 = M2 * M;
  const double M5 = M3 * M2;

  const double gsqrt = std::sqrt(c.Gamma_ee * c.Gamma_mumu);

  const double H1 = std::cos(phi) + C2 * std::cos(phi + Phi);
  const double H2 = std::sin(phi) + C2 * std::sin(phi + Phi);

  // Continuum: extra 1/W'^2 shifts I3 -> I10.
  out.sigmaC = (4.0 * PI * alpha * alpha * A2 / (3.0 * W8)) * I.I10;

  // Resonance: extra 1/(1-x) shifts J1 -> J2 and adds one overall 1/W^2.
  {
    const double fac = (C1 * C1) * (1.0 + C2 * C2 + 2.0 * C2 * std::cos(Phi));
    const double prefR = fac * (6.0 * PI * A2 * c.Gamma_ee * c.Gamma_mumu)
                       / (c.Gamma * M5 * W4);
    out.sigmaR = prefR * (2.0 * std::real(K.J2));
  }

  // Interference: extra 1/(1-x) shifts J2 -> J3 and J1 -> J2, and adds one
  // overall 1/W^2 relative to the old W^{-6} implementation.
  {
    const double prefI = (4.0 * PI * alpha * A2 * gsqrt) / (M * W2);

    const double termA = C1 * ( H2 / (W4 * M2) - H1 / (c.Gamma * W4 * M) );
    const double termC = C1 * H1 / (c.Gamma * W2 * M3);

    const cd inside = cd(termA, 0.0) * K.J3 + cd(termC, 0.0) * K.J2;
    out.sigmaI = prefI * (2.0 * std::real(inside));
  }

  apply_vp(W, out.sigmaC, out.sigmaR, out.sigmaI, vp);

  out.sigma  = (out.sigmaC + out.sigmaR + out.sigmaI) * c.unit_conv;
  out.sigmaC *= c.unit_conv;
  out.sigmaR *= c.unit_conv;
  out.sigmaI *= c.unit_conv;

  return out;
}

} // namespace isr_sigma2body
