#pragma once
// sigma2bodyW6_isr.hpp
//
// 0th-order ISR-corrected cross section for a two-body channel with the
// Born-level prefactor
//
//   sigma_Born,noqf(W) ~ 4 pi alpha^2 Acal^2 / (3 W^6) * |1 + R(W)|^2,
//
// i.e. the same W dependence as in the original two-body note. Relative to
// the W^{-8} and W^{-7} variants, this is the unshifted construction.
//
// The basic kernels used here are
//   I3 = ∫ F /(1-x)^3,
//   J1 = ∫ F /[(1-x)   (A+Bx)],
//   J2 = ∫ F /[(1-x)^2 (A+Bx)],
// with
//   J1 = I0 /(A+B) + B/(A+B) I1,
//   J2 = I4 /(A+B) + B/(A+B) J1.

#include <complex>
#include <cmath>
#include <stdexcept>

#include "sigma2body_isr.hpp"   // shared Consts / SigmaParts / VPOptions / helpers

namespace isr_sigma2bodyW6 {

using cd = std::complex<double>;
#ifndef PI
static constexpr double PI = isr_sigma2body::PI;
#endif

using Consts     = isr_sigma2body::Consts;
using IsCache    = isr_sigma2body::IsCache;
using SigmaParts = isr_sigma2body::SigmaParts;
using VPOptions  = isr_sigma2body::VPOptions;
using ShiftedKernels = isr_sigma2body::ShiftedKernels;

using isr_sigma2body::beta_KF;
using isr_sigma2body::delta_KF;
using isr_sigma2body::make_AB_propagator;
using isr_sigma2body::apply_vp;
using isr_sigma2body::build_shifted_kernels;

inline SigmaParts sigma2bodyW6_ISR(double W,
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
    throw std::invalid_argument("sigma2bodyW6_ISR: require W>0 and M>0.");
  }
  if (!(c.Wmin > 0.0 && c.Wmin < W)) {
    throw std::invalid_argument("sigma2bodyW6_ISR: require 0 < Wmin < W.");
  }
  if (!(c.Gamma > 0.0 && c.Gamma_ee > 0.0 && c.Gamma_mumu > 0.0)) {
    throw std::invalid_argument("sigma2bodyW6_ISR: require positive widths Gamma, Gamma_ee, Gamma_mumu.");
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

  const auto J = build_shifted_kernels(I, Aprop, Bprop);

  const double alpha = c.alpha;
  const double A2 = Acal * Acal;
  const double W2 = W * W;
  const double W4 = W2 * W2;
  const double W6 = W4 * W2;
  const double M2 = M * M;
  const double M3 = M2 * M;
  const double M5 = M3 * M2;

  const double gsqrt = std::sqrt(c.Gamma_ee * c.Gamma_mumu);

  const double h0 = (C1 * C1) * (1.0 + C2 * C2 + 2.0 * C2 * std::cos(Phi));
  const double H1 = std::cos(phi) + C2 * std::cos(phi + Phi);
  const double H2 = std::sin(phi) + C2 * std::sin(phi + Phi);

  // Continuum
  out.sigmaC = (4.0 * PI * alpha * alpha * A2 / (3.0 * W6)) * I.I3;

  // Resonance
  {
    const double prefR = h0 * (6.0 * PI * A2 * c.Gamma_ee * c.Gamma_mumu)
                       / (c.Gamma * M5 * W2);
    out.sigmaR = prefR * (2.0 * std::real(J.J1));
  }

  // Interference
  {
    const double prefI = (4.0 * PI * alpha * A2 * gsqrt) / M;

    const double termA = C1 * ( H2 / (W4 * M2) - H1 / (c.Gamma * W4 * M) );
    const double termC = C1 * H1 / (c.Gamma * W2 * M3);

    const cd inside = cd(termA, 0.0) * J.J2 + cd(termC, 0.0) * J.J1;
    out.sigmaI = prefI * (2.0 * std::real(inside));
  }

  apply_vp(W, out.sigmaC, out.sigmaR, out.sigmaI, vp);

  out.sigma  = (out.sigmaC + out.sigmaR + out.sigmaI) * c.unit_conv;
  out.sigmaC *= c.unit_conv;
  out.sigmaR *= c.unit_conv;
  out.sigmaI *= c.unit_conv;

  return out;
}

} // namespace isr_sigma2bodyW6
