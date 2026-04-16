#pragma once
// sigma2body_isr_qf3.hpp
//
// Final ISR cross section for the modified two-body Born prefactor
//   sigma_Born,noqf(W) ~ W^{-8}
// including q_f^3 corrections up to 3rd order:
//
//   q_f^3(W) ≈ q0 + q1 W + q2 W^2 + q3 W^3.
//
// Relative to the previous W^{-6} implementation, every ISR kernel receives
// one extra factor 1/(1-x). This leads to the shifts
//
//   I5 -> I11,   I4 -> I3,   I8 -> I5,
//   I1 -> J1,    I6 -> K1,
//   K1 -> K2,    K2 -> K3,
//
// with
//   J1 = ∫ F /[(1-x)   (A+Bx)],
//   J2 = ∫ F /[(1-x)^2 (A+Bx)],
//   J3 = ∫ F /[(1-x)^3 (A+Bx)],
//   K1 = ∫ F /[(1-x)^(1/2) (A+Bx)],
//   K2 = ∫ F /[(1-x)^(3/2) (A+Bx)],
//   K3 = ∫ F /[(1-x)^(5/2) (A+Bx)].
//
// These are generated recursively from the basis integrals in
// isr_integrals_I0_I11.hpp, so no additional special-function core is needed.

#include <complex>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <vector>

#include "sigma2body_isr.hpp"

namespace isr_sigma2body {

struct Corr123Parts {
  SigmaParts ord1;
  SigmaParts ord2;
  SigmaParts ord3;
  SigmaParts sum;
};

struct HalfShiftedKernels {
  cd K1 = cd(0.0, 0.0); // ∫ F /[(1-x)^(1/2) (A+Bx)]
  cd K2 = cd(0.0, 0.0); // ∫ F /[(1-x)^(3/2) (A+Bx)]
  cd K3 = cd(0.0, 0.0); // ∫ F /[(1-x)^(5/2) (A+Bx)]
};

inline HalfShiftedKernels build_half_shifted_kernels(const isr_i0i11::Result& I, cd Aprop, cd Bprop) {
  HalfShiftedKernels K;
  const cd S  = Aprop + Bprop;
  const cd BS = Bprop / S;

  K.K1 = cd(I.I7, 0.0) / S + BS * I.I6;
  K.K2 = cd(I.I8, 0.0) / S + BS * K.K1;
  K.K3 = cd(I.I5, 0.0) / S + BS * K.K2;
  return K;
}

inline Corr123Parts sigma2body_ISR_corr123(double W,
                                           double M,
                                           double Acal,
                                           double C1,
                                           double C2,
                                           double phi,
                                           double Phi,
                                           double q1,
                                           double q2,
                                           double q3,
                                           const Consts& c = Consts{},
                                           const isr_i0i11::Options& iopt = isr_i0i11::Options{},
                                           const VPOptions& vp = VPOptions{},
                                           const IsCache& cache = IsCache{})
{
  Corr123Parts out;

  if (!(W > 0.0) || !(M > 0.0)) {
    throw std::invalid_argument("sigma2body_ISR_corr123: require W>0 and M>0.");
  }
  if (!(c.Wmin > 0.0 && c.Wmin < W)) {
    throw std::invalid_argument("sigma2body_ISR_corr123: require 0 < Wmin < W.");
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

  const auto J = build_shifted_kernels(I, Aprop, Bprop);
  const auto K = build_half_shifted_kernels(I, Aprop, Bprop);

  const int nb = I.n_used_beta + I.q_used_romI + I.q_used_23;
  const int n2 = I.n_used_2f1;
  const bool conv_beta_all = I.conv_beta && I.converged_romI && I.converged_23;
  const bool conv_2f1_all  = I.conv_2f1;

  const double alpha = c.alpha;
  const double A2 = Acal * Acal;

  const double W2 = W * W;
  const double W3 = W2 * W;
  const double W4 = W2 * W2;
  const double W5 = W4 * W;
  const double W6 = W3 * W3;
  const double W7 = W6 * W;

  const double M2 = M * M;
  const double M3 = M2 * M;
  const double M4 = M2 * M2;
  const double M5 = M3 * M2;

  const double gprod = c.Gamma_ee * c.Gamma_mumu;
  const double gsqrt = std::sqrt(gprod);

  const double h0 = (C1 * C1) * (1.0 + C2 * C2 + 2.0 * C2 * std::cos(Phi));
  const double H1 = std::cos(phi) + C2 * std::cos(phi + Phi);
  const double H2 = std::sin(phi) + C2 * std::sin(phi + Phi);
  const double h1 = C1 * H1;
  const double h2 = C1 * H2;

  // ---------------- Order 1 ----------------
  {
    // Continuum: I5 -> I11 and W^{-5} -> W^{-7}.
    out.ord1.sigmaC = q1 * (4.0 * PI * alpha * alpha * A2 / (3.0 * W7)) * I.I11;

    // Resonance: K1 -> K2 and W^{-1} -> W^{-3}.
    const double prefR1 = q1 * h0 * (6.0 * PI * A2 * gprod) / (c.Gamma * M5 * W3);
    out.ord1.sigmaR = prefR1 * (2.0 * std::real(K.K2));

    // Interference: K1 -> K2, K2 -> K3, and one overall extra 1/W^2.
    const double prefI1 = q1 * (4.0 * PI * alpha * A2 * gsqrt) / (c.Gamma * M4 * W);
    const cd inner = cd(h1 / W2, 0.0) * K.K2
                   + cd((M * c.Gamma * h2 - M2 * h1) / W4, 0.0) * K.K3;
    out.ord1.sigmaI = prefI1 * (2.0 * std::real(inner));

    apply_vp(W, out.ord1.sigmaC, out.ord1.sigmaR, out.ord1.sigmaI, vp);

    out.ord1.sigma  = (out.ord1.sigmaC + out.ord1.sigmaR + out.ord1.sigmaI) * c.unit_conv;
    out.ord1.sigmaC *= c.unit_conv;
    out.ord1.sigmaR *= c.unit_conv;
    out.ord1.sigmaI *= c.unit_conv;

    out.ord1.n_beta = nb;
    out.ord1.n_2f1  = n2;
    out.ord1.conv_beta = conv_beta_all;
    out.ord1.conv_2f1  = conv_2f1_all;
  }

  // ---------------- Order 2 ----------------
  {
    // Continuum: I4 -> I3 and W^{-4} -> W^{-6}.
    out.ord2.sigmaC = q2 * (4.0 * PI * alpha * alpha * A2 / (3.0 * W6)) * I.I3;

    // Resonance: I1 -> J1 and W^{0} -> W^{-2}.
    const double prefR2 = q2 * h0 * (6.0 * PI * A2 * gprod) / (c.Gamma * M5 * W2);
    out.ord2.sigmaR = prefR2 * (2.0 * std::real(J.J1));

    // Interference: I1 -> J1, J1 -> J2, and one overall extra 1/W^2.
    const double prefI2 = q2 * (4.0 * PI * alpha * A2 * gsqrt) / (c.Gamma * M4);
    const cd inner = cd(h1 / W2, 0.0) * J.J1
                   + cd((M * c.Gamma * h2 - M2 * h1) / W4, 0.0) * J.J2;
    out.ord2.sigmaI = prefI2 * (2.0 * std::real(inner));

    apply_vp(W, out.ord2.sigmaC, out.ord2.sigmaR, out.ord2.sigmaI, vp);

    out.ord2.sigma  = (out.ord2.sigmaC + out.ord2.sigmaR + out.ord2.sigmaI) * c.unit_conv;
    out.ord2.sigmaC *= c.unit_conv;
    out.ord2.sigmaR *= c.unit_conv;
    out.ord2.sigmaI *= c.unit_conv;

    out.ord2.n_beta = nb;
    out.ord2.n_2f1  = n2;
    out.ord2.conv_beta = conv_beta_all;
    out.ord2.conv_2f1  = conv_2f1_all;
  }

  // ---------------- Order 3 ----------------
  {
    // Continuum: I8 -> I5 and W^{-3} -> W^{-5}.
    out.ord3.sigmaC = q3 * (4.0 * PI * alpha * alpha * A2 / (3.0 * W5)) * I.I5;

    // Resonance: I6 -> K1 and W^{+1} -> W^{-1}.
    const double prefR3 = q3 * h0 * (6.0 * PI * A2 * gprod) / (c.Gamma * M5 * W);
    out.ord3.sigmaR = prefR3 * (2.0 * std::real(K.K1));

    // Interference: I6 -> K1, K1 -> K2, and one overall extra 1/W^2.
    const double prefI3 = q3 * (4.0 * PI * alpha * A2 * gsqrt) * W / (c.Gamma * M4);
    const cd inner = cd(h1 / W2, 0.0) * K.K1
                   + cd((M * c.Gamma * h2 - M2 * h1) / W4, 0.0) * K.K2;
    out.ord3.sigmaI = prefI3 * (2.0 * std::real(inner));

    apply_vp(W, out.ord3.sigmaC, out.ord3.sigmaR, out.ord3.sigmaI, vp);

    out.ord3.sigma  = (out.ord3.sigmaC + out.ord3.sigmaR + out.ord3.sigmaI) * c.unit_conv;
    out.ord3.sigmaC *= c.unit_conv;
    out.ord3.sigmaR *= c.unit_conv;
    out.ord3.sigmaI *= c.unit_conv;

    out.ord3.n_beta = nb;
    out.ord3.n_2f1  = n2;
    out.ord3.conv_beta = conv_beta_all;
    out.ord3.conv_2f1  = conv_2f1_all;
  }

  out.sum.sigmaC = out.ord1.sigmaC + out.ord2.sigmaC + out.ord3.sigmaC;
  out.sum.sigmaR = out.ord1.sigmaR + out.ord2.sigmaR + out.ord3.sigmaR;
  out.sum.sigmaI = out.ord1.sigmaI + out.ord2.sigmaI + out.ord3.sigmaI;
  out.sum.sigma  = out.ord1.sigma  + out.ord2.sigma  + out.ord3.sigma;

  out.sum.n_beta = out.ord1.n_beta + out.ord2.n_beta + out.ord3.n_beta;
  out.sum.n_2f1  = out.ord1.n_2f1  + out.ord2.n_2f1  + out.ord3.n_2f1;
  out.sum.conv_beta = out.ord1.conv_beta && out.ord2.conv_beta && out.ord3.conv_beta;
  out.sum.conv_2f1  = out.ord1.conv_2f1  && out.ord2.conv_2f1  && out.ord3.conv_2f1;

  return out;
}

inline std::vector<double> convert_p2q(double p0, double p1, double p2, double p3, double W0) {
  const double q0 = p0 - p1 * W0 + p2 * W0 * W0 - p3 * W0 * W0 * W0;
  const double q1 = p1 - 2.0 * p2 * W0 + 3.0 * p3 * W0 * W0;
  const double q2 = p2 - 3.0 * p3 * W0;
  const double q3 = p3;
  return {q0, q1, q2, q3};
}

inline SigmaParts sigma2body_ISR_qf3(double W,
                                     double M,
                                     double Acal,
                                     double C1,
                                     double C2,
                                     double phi,
                                     double Phi,
                                     double q0,
                                     double q1,
                                     double q2,
                                     double q3,
                                     const Consts& c = Consts{},
                                     const isr_i0i11::Options& iopt = isr_i0i11::Options{},
                                     const VPOptions& vp = VPOptions{})
{
  if (!(W > 0.0) || !(M > 0.0)) {
    throw std::invalid_argument("sigma2body_ISR_qf3: require W>0 and M>0.");
  }
  if (!(c.Wmin > 0.0 && c.Wmin < W)) {
    throw std::invalid_argument("sigma2body_ISR_qf3: require 0 < Wmin < W.");
  }

  cd Aprop, Bprop;
  make_AB_propagator(W, M, c, Aprop, Bprop);

  const double beta  = beta_KF(W, c);
  const double delta = delta_KF(beta, c);
  const double b = 1.0 - (c.Wmin / W) * (c.Wmin / W);

  const auto Iall = isr_i0i11::compute_all(b, beta, delta, Aprop, Bprop, iopt);

  IsCache cache;
  cache.useCache = true;
  cache.Iall = Iall;

  auto base = sigma2body_ISR(W, M, Acal, C1, C2, phi, Phi, c, iopt, vp, cache);
  base.sigmaC *= q0;
  base.sigmaR *= q0;
  base.sigmaI *= q0;
  base.sigma  *= q0;

  const auto corr = sigma2body_ISR_corr123(W, M, Acal, C1, C2, phi, Phi,
                                           q1, q2, q3, c, iopt, vp, cache);

  SigmaParts out;
  out.sigmaC = base.sigmaC + corr.sum.sigmaC;
  out.sigmaR = base.sigmaR + corr.sum.sigmaR;
  out.sigmaI = base.sigmaI + corr.sum.sigmaI;
  out.sigma  = base.sigma  + corr.sum.sigma;

  out.n_beta = base.n_beta + corr.sum.n_beta;
  out.n_2f1  = base.n_2f1  + corr.sum.n_2f1;
  out.conv_beta = base.conv_beta && corr.sum.conv_beta;
  out.conv_2f1  = base.conv_2f1  && corr.sum.conv_2f1;

  return out;
}

} // namespace isr_sigma2body
