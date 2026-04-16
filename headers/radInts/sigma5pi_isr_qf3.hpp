#pragma once
// sigma5pi_isr_qf3.hpp
//
// Final merged ISR cross section with q_f^3 correction via cubic polynomial
// coefficients q0..q3 (as defined in Wenhao Ya note Eq.(6)).
//
// Strategy:
//   - Use 0th-order analytic core sigma5pi_ISR (no-qf), optionally VP corrected.
//   - Multiply the 0th-order result by q0.
//   - Add order-1..3 corrections from sigma5pi_ISR_corr123 (already includes q1..q3).
//
// Return SigmaParts with the same (C,R,I,total) decomposition and diagnostics.

#include <stdexcept>
#include "sigma5pi_isr.hpp"
#include "sigma5pi_isr_qcorr_123.hpp"

namespace isr_sigma5pi {

inline std::vector<double> convert_p2q(double p0, double p1, double p2, double p3, double W0) {
    double q0 = p0 - p1*W0 + p2*W0*W0 - p3*W0*W0*W0;
    double q1 = p1 - 2*p2*W0 + 3*p3*W0*W0;
    double q2 = p2 - 3*p3*W0;
    double q3 = p3;
    return {q0, q1, q2, q3};
}

inline SigmaParts sigma5pi_ISR_qf3(double W,
                                   double M,
                                   double Acal,
                                   double C1,
                                   double C2,
                                   double phi,
                                   double Phi,
                                   double q0 = 0.176892,
                                   double q1 = -0.3246,
                                   double q2 = 0.0201564,
                                   double q3 = 0.123254,
                                   const Consts& c = Consts{},
                                   const isr_i0i4::Options& iopt = isr_i0i4::Options{},
                                   const isr6::Options& i6opt = isr6::Options{},
                                   const VPOptions& vp = VPOptions{})
{
  // radiator params
  const double beta  = beta_KF(W, c);
  const double delta = delta_KF(beta, c);
  const double b = 1.0 - (c.Wmin/W)*(c.Wmin/W);

  // propagator params
  cd Aprop, Bprop;
  make_AB_propagator(W, M, c, Aprop, Bprop);
  const cd AB = Aprop + Bprop;

  // needed integrals
  const auto I04 = isr_i0i4::compute_all(b, beta, delta, Aprop, Bprop, iopt);
  const auto I5789 = isr_int::compute_all(b, beta, delta);
  const auto I6res = isr6::compute_I6(Aprop, Bprop, b, beta, delta, i6opt);

    // cache the Is calculations
    IsCache cache;
    cache.useCache = true;
    cache.I04 = I04;
    cache.I5789 = I5789;
    cache.I6res = I6res;

  // 0th order (no qf^3), optional VP already applied inside
  auto base = sigma5pi_ISR(W, M, Acal, C1, C2, phi, Phi, c, iopt, vp, cache);

  // scale by q0 (note Eq.(10): σ_ISR ≈ Σ_{n=0}^3 σ_ISR_n, and σ_ISR_0 = q0 * σ_ISR_noqf)
  base.sigmaC *= q0;
  base.sigmaR *= q0;
  base.sigmaI *= q0;
  base.sigma  *= q0;

  // corrections (orders 1..3), optional VP applied inside
  const auto corr = sigma5pi_ISR_corr123(W, M, Acal, C1, C2, phi, Phi, q1, q2, q3, c, iopt, i6opt, vp, cache);

  SigmaParts out;
  out.sigmaC = base.sigmaC + corr.sum.sigmaC;
  out.sigmaR = base.sigmaR + corr.sum.sigmaR;
  out.sigmaI = base.sigmaI + corr.sum.sigmaI;
  out.sigma  = base.sigma  + corr.sum.sigma;

  // diagnostics (simple aggregation)
  out.n_beta = base.n_beta + corr.sum.n_beta;
  out.n_2f1  = base.n_2f1  + corr.sum.n_2f1;
  out.conv_beta = base.conv_beta && corr.sum.conv_beta;
  out.conv_2f1  = base.conv_2f1  && corr.sum.conv_2f1;

  return out;
}

} // namespace isr_sigma5pi