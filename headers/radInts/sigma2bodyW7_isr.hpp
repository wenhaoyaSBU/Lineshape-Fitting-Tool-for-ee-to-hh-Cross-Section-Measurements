#pragma once
// sigma2bodyW7_isr.hpp
//
// 0th-order ISR-corrected cross section for a two-body channel with the
// modified Born-level prefactor
//
//   sigma_Born,noqf(W) ~ 4 pi alpha^2 Acal^2 / (3 W^7) * |1 + R(W)|^2,
//
// where the resonance-amplitude structure R(W) is kept identical to the
// existing two-body notes/headers.
//
// Relative to the original W^{-6} formulation, the ISR substitution
//   W' = W sqrt(1-x)
// introduces one extra factor 1/sqrt(1-x) in every kernel, together with one
// extra overall factor 1/W.
//
// Therefore, for the 0th-order piece:
//   continuum:   I3 -> I11
//   resonance:   J1 -> K2
//   interference second kernel: J2 -> K3
//
// with
//   J1 = ∫ F /[(1-x)   (A+Bx)],
//   K1 = ∫ F /[(1-x)^(1/2) (A+Bx)],
//   K2 = ∫ F /[(1-x)^(3/2) (A+Bx)],
//   K3 = ∫ F /[(1-x)^(5/2) (A+Bx)].

#include <complex>
#include <cmath>
#include <stdexcept>

#include "sigma2body_isr.hpp"   // shared Consts / SigmaParts / VPOptions / helpers

namespace isr_sigma2bodyW7 {

using cd = std::complex<double>;
#ifndef PI
static constexpr double PI = isr_sigma2body::PI;
#endif

using Consts     = isr_sigma2body::Consts;
using IsCache    = isr_sigma2body::IsCache;
using SigmaParts = isr_sigma2body::SigmaParts;
using VPOptions  = isr_sigma2body::VPOptions;

using isr_sigma2body::beta_KF;
using isr_sigma2body::delta_KF;
using isr_sigma2body::make_AB_propagator;
using isr_sigma2body::apply_vp;
using isr_sigma2body::build_shifted_kernels;

struct HalfShiftedKernels {
  cd K1 = cd(0.0, 0.0); // ∫ F /[(1-x)^(1/2) (A+Bx)]
  cd K2 = cd(0.0, 0.0); // ∫ F /[(1-x)^(3/2) (A+Bx)]
  cd K3 = cd(0.0, 0.0); // ∫ F /[(1-x)^(5/2) (A+Bx)]
};

inline HalfShiftedKernels build_half_shifted_kernels(const isr_i0i11::Result& I,
                                                     cd Aprop,
                                                     cd Bprop) {
  HalfShiftedKernels K;
  const cd S  = Aprop + Bprop;
  const cd BS = Bprop / S;

  K.K1 = cd(I.I7, 0.0) / S + BS * I.I6;
  K.K2 = cd(I.I8, 0.0) / S + BS * K.K1;
  K.K3 = cd(I.I5, 0.0) / S + BS * K.K2;
  return K;
}

inline SigmaParts sigma2bodyW7_ISR(double W,
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
    throw std::invalid_argument("sigma2bodyW7_ISR: require W>0 and M>0.");
  }
  if (!(c.Wmin > 0.0 && c.Wmin < W)) {
    throw std::invalid_argument("sigma2bodyW7_ISR: require 0 < Wmin < W.");
  }
  if (!(c.Gamma > 0.0 && c.Gamma_ee > 0.0 && c.Gamma_mumu > 0.0)) {
    throw std::invalid_argument("sigma2bodyW7_ISR: require positive widths Gamma, Gamma_ee, Gamma_mumu.");
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

  out.n_beta = I.n_used_beta + I.q_used_romI + I.q_used_23;
  out.n_2f1  = I.n_used_2f1;
  out.conv_beta = I.conv_beta && I.converged_romI && I.converged_23;
  out.conv_2f1  = I.conv_2f1;

  const auto K = build_half_shifted_kernels(I, Aprop, Bprop);

  const double alpha = c.alpha;
  const double A2 = Acal * Acal;
  const double W2 = W * W;
  const double W4 = W2 * W2;
  const double W7 = W4 * W2 * W;
  const double M2 = M * M;
  const double M3 = M2 * M;
  const double M5 = M3 * M2;

  const double gsqrt = std::sqrt(c.Gamma_ee * c.Gamma_mumu);

  const double h0 = (C1 * C1) * (1.0 + C2 * C2 + 2.0 * C2 * std::cos(Phi));
  const double H1 = std::cos(phi) + C2 * std::cos(phi + Phi);
  const double H2 = std::sin(phi) + C2 * std::sin(phi + Phi);

  // Continuum: I3 -> I11 and one extra overall 1/W.
  out.sigmaC = (4.0 * PI * alpha * alpha * A2 / (3.0 * W7)) * I.I11;

  // Resonance: J1 -> K2 and one extra overall 1/W.
  {
    const double prefR = h0 * (6.0 * PI * A2 * c.Gamma_ee * c.Gamma_mumu)
                       / (c.Gamma * M5 * W * W * W);
    out.sigmaR = prefR * (2.0 * std::real(K.K2));
  }

  // Interference: J2 -> K3, J1 -> K2, and one extra overall 1/W.
  {
    const double prefI = (4.0 * PI * alpha * A2 * gsqrt) / (M * W);

    const double termA = C1 * ( H2 / (W4 * M2) - H1 / (c.Gamma * W4 * M) );
    const double termC = C1 * H1 / (c.Gamma * W2 * M3);

    const cd inside = cd(termA, 0.0) * K.K3 + cd(termC, 0.0) * K.K2;
    out.sigmaI = prefI * (2.0 * std::real(inside));
  }

  apply_vp(W, out.sigmaC, out.sigmaR, out.sigmaI, vp);

  out.sigma  = (out.sigmaC + out.sigmaR + out.sigmaI) * c.unit_conv;
  out.sigmaC *= c.unit_conv;
  out.sigmaR *= c.unit_conv;
  out.sigmaI *= c.unit_conv;

  return out;
}

} // namespace isr_sigma2bodyW7
