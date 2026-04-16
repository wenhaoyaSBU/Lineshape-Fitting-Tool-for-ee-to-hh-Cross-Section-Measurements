#pragma once
// sigma2bodyW7_isr_qf5.hpp
//
// ISR-corrected two-body cross section with modified 0th-order Born prefactor
//   sigma_Born,noqf(W) ~ W^{-7}
// including q_f^3 corrections up to 5th order:
//
//   q_f^3(W) ≈ q0 + q1 W + q2 W^2 + q3 W^3 + q4 W^4 + q5 W^5.
//
// Basis-closure check before implementation:
//   Relative to the original W^{-6} construction, every kernel acquires one
//   extra factor 1/sqrt(1-x). Through 5th order this requires
//
//     continuum:   I11, I3, I5, I4, I8, I0,
//     resonance:   K2, J1, K1, I1, I6, L1,
//     interference second kernels: K3, J2, K2, J1, K1, I1,
//
//   where the only "new-looking" object is
//     L1 = ∫ F (1-x)/(A+Bx)
//        = ((A+B) I1 - I2) / B.
//
// Therefore everything still closes within I0..I11; no new special-function
// core is required.

#include <complex>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "sigma2bodyW7_isr.hpp"

namespace isr_sigma2bodyW7 {

struct Corr12345Parts {
  SigmaParts ord1;
  SigmaParts ord2;
  SigmaParts ord3;
  SigmaParts ord4;
  SigmaParts ord5;
  SigmaParts sum;
};

inline cd build_linear_positive_kernel(const isr_i0i11::Result& I,
                                       cd Aprop,
                                       cd Bprop) {
  // L1 = ∫ F (1-x)/(A+Bx) = ((A+B) I1 - I2)/B
  return (((Aprop + Bprop) * I.I1) - cd(I.I2, 0.0)) / Bprop;
}

inline Corr12345Parts sigma2bodyW7_ISR_corr12345(double W,
                                                 double M,
                                                 double Acal,
                                                 double C1,
                                                 double C2,
                                                 double phi,
                                                 double Phi,
                                                 double q1,
                                                 double q2,
                                                 double q3,
                                                 double q4,
                                                 double q5,
                                                 const Consts& c = Consts{},
                                                 const isr_i0i11::Options& iopt = isr_i0i11::Options{},
                                                 const VPOptions& vp = VPOptions{},
                                                 const IsCache& cache = IsCache{})
{
  Corr12345Parts out;

  if (!(W > 0.0) || !(M > 0.0)) {
    throw std::invalid_argument("sigma2bodyW7_ISR_corr12345: require W>0 and M>0.");
  }
  if (!(c.Wmin > 0.0 && c.Wmin < W)) {
    throw std::invalid_argument("sigma2bodyW7_ISR_corr12345: require 0 < Wmin < W.");
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
  const cd   L1 = build_linear_positive_kernel(I, Aprop, Bprop);

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
  const double mix = M * c.Gamma * h2 - M2 * h1;

  auto finalize = [&](SigmaParts& s) {
    apply_vp(W, s.sigmaC, s.sigmaR, s.sigmaI, vp);
    s.sigma  = (s.sigmaC + s.sigmaR + s.sigmaI) * c.unit_conv;
    s.sigmaC *= c.unit_conv;
    s.sigmaR *= c.unit_conv;
    s.sigmaI *= c.unit_conv;
    s.n_beta = nb;
    s.n_2f1  = n2;
    s.conv_beta = conv_beta_all;
    s.conv_2f1  = conv_2f1_all;
  };

  // ---------------- Order 1 ----------------
  {
    out.ord1.sigmaC = q1 * (4.0 * PI * alpha * alpha * A2 / (3.0 * W6)) * I.I3;

    const double prefR1 = q1 * h0 * (6.0 * PI * A2 * gprod) / (c.Gamma * M5 * W2);
    out.ord1.sigmaR = prefR1 * (2.0 * std::real(J.J1));

    const double prefI1 = q1 * (4.0 * PI * alpha * A2 * gsqrt) / (c.Gamma * M4);
    const cd inner = cd(h1 / W2, 0.0) * J.J1
                   + cd(mix / W4, 0.0) * J.J2;
    out.ord1.sigmaI = prefI1 * (2.0 * std::real(inner));

    finalize(out.ord1);
  }

  // ---------------- Order 2 ----------------
  {
    out.ord2.sigmaC = q2 * (4.0 * PI * alpha * alpha * A2 / (3.0 * W5)) * I.I5;

    const double prefR2 = q2 * h0 * (6.0 * PI * A2 * gprod) / (c.Gamma * M5 * W);
    out.ord2.sigmaR = prefR2 * (2.0 * std::real(K.K1));

    const double prefI2 = q2 * (4.0 * PI * alpha * A2 * gsqrt) * W / (c.Gamma * M4);
    const cd inner = cd(h1 / W2, 0.0) * K.K1
                   + cd(mix / W4, 0.0) * K.K2;
    out.ord2.sigmaI = prefI2 * (2.0 * std::real(inner));

    finalize(out.ord2);
  }

  // ---------------- Order 3 ----------------
  {
    out.ord3.sigmaC = q3 * (4.0 * PI * alpha * alpha * A2 / (3.0 * W4)) * I.I4;

    const double prefR3 = q3 * h0 * (6.0 * PI * A2 * gprod) / (c.Gamma * M5);
    out.ord3.sigmaR = prefR3 * (2.0 * std::real(I.I1));

    const double prefI3 = q3 * (4.0 * PI * alpha * A2 * gsqrt) * W2 / (c.Gamma * M4);
    const cd inner = cd(h1 / W2, 0.0) * I.I1
                   + cd(mix / W4, 0.0) * J.J1;
    out.ord3.sigmaI = prefI3 * (2.0 * std::real(inner));

    finalize(out.ord3);
  }

  // ---------------- Order 4 ----------------
  {
    out.ord4.sigmaC = q4 * (4.0 * PI * alpha * alpha * A2 / (3.0 * W3)) * I.I8;

    const double prefR4 = q4 * h0 * (6.0 * PI * A2 * gprod) * W / (c.Gamma * M5);
    out.ord4.sigmaR = prefR4 * (2.0 * std::real(I.I6));

    const double prefI4 = q4 * (4.0 * PI * alpha * A2 * gsqrt) * W3 / (c.Gamma * M4);
    const cd inner = cd(h1 / W2, 0.0) * I.I6
                   + cd(mix / W4, 0.0) * K.K1;
    out.ord4.sigmaI = prefI4 * (2.0 * std::real(inner));

    finalize(out.ord4);
  }

  // ---------------- Order 5 ----------------
  {
    out.ord5.sigmaC = q5 * (4.0 * PI * alpha * alpha * A2 / (3.0 * W2)) * I.I0;

    const double prefR5 = q5 * h0 * (6.0 * PI * A2 * gprod) * W2 / (c.Gamma * M5);
    out.ord5.sigmaR = prefR5 * (2.0 * std::real(L1));

    const double prefI5 = q5 * (4.0 * PI * alpha * A2 * gsqrt) * W4 / (c.Gamma * M4);
    const cd inner = cd(h1 / W2, 0.0) * L1
                   + cd(mix / W4, 0.0) * I.I1;
    out.ord5.sigmaI = prefI5 * (2.0 * std::real(inner));

    finalize(out.ord5);
  }

  out.sum.sigmaC = out.ord1.sigmaC + out.ord2.sigmaC + out.ord3.sigmaC + out.ord4.sigmaC + out.ord5.sigmaC;
  out.sum.sigmaR = out.ord1.sigmaR + out.ord2.sigmaR + out.ord3.sigmaR + out.ord4.sigmaR + out.ord5.sigmaR;
  out.sum.sigmaI = out.ord1.sigmaI + out.ord2.sigmaI + out.ord3.sigmaI + out.ord4.sigmaI + out.ord5.sigmaI;
  out.sum.sigma  = out.ord1.sigma  + out.ord2.sigma  + out.ord3.sigma  + out.ord4.sigma  + out.ord5.sigma;

  out.sum.n_beta = out.ord1.n_beta + out.ord2.n_beta + out.ord3.n_beta + out.ord4.n_beta + out.ord5.n_beta;
  out.sum.n_2f1  = out.ord1.n_2f1  + out.ord2.n_2f1  + out.ord3.n_2f1  + out.ord4.n_2f1  + out.ord5.n_2f1;
  out.sum.conv_beta = out.ord1.conv_beta && out.ord2.conv_beta && out.ord3.conv_beta && out.ord4.conv_beta && out.ord5.conv_beta;
  out.sum.conv_2f1  = out.ord1.conv_2f1  && out.ord2.conv_2f1  && out.ord3.conv_2f1  && out.ord4.conv_2f1  && out.ord5.conv_2f1;

  return out;
}

inline std::vector<double> convert_p2q(double p0,
                                       double p1,
                                       double p2,
                                       double p3,
                                       double p4,
                                       double p5,
                                       double W0) {
  const double W02 = W0 * W0;
  const double W03 = W02 * W0;
  const double W04 = W03 * W0;
  const double W05 = W04 * W0;

  const double q0 = p0 - p1 * W0 + p2 * W02 - p3 * W03 + p4 * W04 - p5 * W05;
  const double q1 = p1 - 2.0 * p2 * W0 + 3.0 * p3 * W02 - 4.0 * p4 * W03 + 5.0 * p5 * W04;
  const double q2 = p2 - 3.0 * p3 * W0 + 6.0 * p4 * W02 - 10.0 * p5 * W03;
  const double q3 = p3 - 4.0 * p4 * W0 + 10.0 * p5 * W02;
  const double q4 = p4 - 5.0 * p5 * W0;
  const double q5 = p5;
  return {q0, q1, q2, q3, q4, q5};
}

inline SigmaParts sigma2bodyW7_ISR_qf5(double W,
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
                                       double q4,
                                       double q5,
                                       const Consts& c = Consts{},
                                       const isr_i0i11::Options& iopt = isr_i0i11::Options{},
                                       const VPOptions& vp = VPOptions{})
{
  if (!(W > 0.0) || !(M > 0.0)) {
    throw std::invalid_argument("sigma2bodyW7_ISR_qf5: require W>0 and M>0.");
  }
  if (!(c.Wmin > 0.0 && c.Wmin < W)) {
    throw std::invalid_argument("sigma2bodyW7_ISR_qf5: require 0 < Wmin < W.");
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

  auto base  = sigma2bodyW7_ISR(W, M, Acal, C1, C2, phi, Phi, c, iopt, vp, cache);
  auto corr  = sigma2bodyW7_ISR_corr12345(W, M, Acal, C1, C2, phi, Phi,
                                          q1, q2, q3, q4, q5,
                                          c, iopt, vp, cache);

  SigmaParts out;
  out.sigmaC = q0 * base.sigmaC + corr.sum.sigmaC;
  out.sigmaR = q0 * base.sigmaR + corr.sum.sigmaR;
  out.sigmaI = q0 * base.sigmaI + corr.sum.sigmaI;
  out.sigma  = q0 * base.sigma  + corr.sum.sigma;

  out.n_beta = base.n_beta + corr.sum.n_beta;
  out.n_2f1  = base.n_2f1  + corr.sum.n_2f1;
  out.conv_beta = base.conv_beta && corr.sum.conv_beta;
  out.conv_2f1  = base.conv_2f1  && corr.sum.conv_2f1;

  return out;
}

} // namespace isr_sigma2bodyW7
