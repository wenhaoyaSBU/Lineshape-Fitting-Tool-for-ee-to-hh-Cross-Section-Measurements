#pragma once
// sigma5pi_isr_qcorr_123.hpp
//
// Order 1–3 correction terms for the ISR integral when including q_f^3(W')
// via a cubic polynomial approximation:
//
//   q_f^3(W) ≈ q0 + q1 W + q2 W^2 + q3 W^3   (note Eq.(6))
//
// This header implements σ_ISR,1..σ_ISR,3 in Wenhao Ya's note
// using analytic basis integrals:
//   I0..I4  from isr_integrals_I0_I4.hpp
//   I5,I7,I8,I9 from isr_integrals_I5_I7_I8_I9.hpp
//   I6 from isr_integral_I6.hpp
//
// VP correction (optional) is applied in the same way as the 0th-order core,
// under the slowly-varying Vacc approximation:
//   σ_C -> σ_C * Vacc
//   σ_I -> σ_I * sqrt(Vacc)
//   σ_R unchanged
//
// Units: returns the same SigmaParts fields (already multiplied by c.unit_conv).

#include <complex>
#include <cmath>
#include <stdexcept>
#include <algorithm>

#include "sigma5pi_isr.hpp"                // Consts, SigmaParts, VPOptions, beta_KF, delta_KF, make_AB_propagator, apply_vp
#include "isr_integrals_I0_I4.hpp"         // I0..I4
#include "isr_integrals_I5_I7_I8_I9.hpp"   // I5,I7,I8,I9 (namespace isr_int)
#include "isr_integral_I6.hpp"             // I6 (namespace isr6)

namespace isr_sigma5pi {

struct Corr123Parts {
  SigmaParts ord1;
  SigmaParts ord2;
  SigmaParts ord3;
  SigmaParts sum;  // ord1+ord2+ord3
};

inline Corr123Parts sigma5pi_ISR_corr123(double W,
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
                                         const isr_i0i4::Options& iopt = isr_i0i4::Options{},
                                         const isr6::Options& i6opt = isr6::Options{},
                                         const VPOptions& vp = VPOptions{},
                                         const IsCache& cache = IsCache{})
{
  Corr123Parts out;

  if (!(W > 0.0) || !(M > 0.0)) throw std::invalid_argument("sigma5pi_ISR_corr123: require W>0, M>0.");
  if (!(c.Wmin > 0.0 && c.Wmin < W)) throw std::invalid_argument("sigma5pi_ISR_corr123: require 0 < Wmin < W.");

  // propagator params
  cd Aprop, Bprop;
  make_AB_propagator(W, M, c, Aprop, Bprop);
  const cd AB = Aprop + Bprop;

  isr_i0i4::Result I04; // declare here to assign from cache or compute
  isr_int::I5I7I8I9 I5789;
  isr6::Result I6res;
  if(!cache.useCache) {
  // radiator params
  const double beta  = beta_KF(W, c);
  const double delta = delta_KF(beta, c);
  const double b = 1.0 - (c.Wmin/W)*(c.Wmin/W);

  // needed integrals
  // const auto I04 = isr_i0i4::compute_all(b, beta, delta, Aprop, Bprop, iopt);
  // const auto I5789 = isr_int::compute_all(b, beta, delta);
  // const auto I6res = isr6::compute_I6(Aprop, Bprop, b, beta, delta, i6opt);
  I04 = isr_i0i4::compute_all(b, beta, delta, Aprop, Bprop, iopt);
  I5789 = isr_int::compute_all(b, beta, delta);
  I6res = isr6::compute_I6(Aprop, Bprop, b, beta, delta, i6opt);
  }
  else {
      I04 = cache.I04;
      I5789 = cache.I5789;
      I6res = cache.I6res;
  }

  // diagnostics (rough aggregation)
  const int nb = I04.n_used_beta + I6res.q_used_romI + I6res.q_used_23;
  const int n2 = I04.n_used_2f1;

  // aliases
  const double I0 = I04.I0;
  const cd     I1 = I04.I1;
  const double I2 = I04.I2;
  const double I4 = I04.I4;

  const double I5 = I5789.I5;
  const double I7 = I5789.I7;
  const double I8 = I5789.I8;
  const double I9 = I5789.I9;

  const cd I6 = I6res.I6;

  // constants/powers
  const double alpha = c.alpha;
  const double A2 = Acal*Acal;

  const double W2 = W*W;
  const double W3 = W2*W;
  const double W4 = W2*W2;
  const double W5 = W4*W;
  const double M2 = M*M;
  const double M3 = M2*M;
  const double M4 = M2*M2;
  const double M5 = M3*M2;

  const double gprod  = c.Gamma_ee * c.Gamma_mumu;
  const double gsqrt  = std::sqrt(gprod);

  // h0,h1,h2 (note Eq.(7)), and reuse H1,H2 used in 0th order
  const double h0 = (C1*C1) * (1.0 + C2*C2 + 2.0*C2*std::cos(Phi));
  const double H1 = std::cos(phi) + C2 * std::cos(phi + Phi);
  const double H2 = std::sin(phi) + C2 * std::sin(phi + Phi);
  const double h1 = C1 * H1;
  const double h2 = C1 * H2;

  // ---------------- Order 1 ----------------
  {
    // σ_C,1 = q1 * 4π α^2 A^2 /(3 W^5) * I5
    out.ord1.sigmaC = q1 * (4.0*PI*alpha*alpha*A2 / (3.0*W5)) * I5;

    // T1 = 1/(A+B) I7 + (B/(A+B)) I6
    const cd T1 = cd(I7,0.0)/AB + (Bprop/AB)*I6;

    // σ_R,1 = q1 h0 * 6π A^2 ΓeeΓμμ /(Γ M^5 W) * (T1 + c.c.) = pref * 2Re(T1)
    const double prefR1 = q1 * h0 * (6.0*PI * A2 * gprod) / (c.Gamma * M5 * W);
    out.ord1.sigmaR = prefR1 * (2.0 * std::real(T1));

    // T2 = 1/(A+B) I8 + 1/(A+B)^2 I7 + (B/(A+B))^2 I6
    // const cd T2 = cd(I8,0.0)/AB + cd(I7,0.0)/(AB*AB) + (Bprop/AB)*(Bprop/AB)*I6;
    const cd T2 = cd(I8,0.0)/AB + cd(I7,0.0)*Bprop/(AB*AB) + (Bprop/AB)*(Bprop/AB)*I6;

    // σ_I,1 = q1 * 4π α A^2 sqrt(ΓeeΓμμ) * W /(Γ M^4) * ( [h1/W^2 T1 + (MΓ h2 - M^2 h1)/W^4 T2] + c.c.)
    const double prefI1 = q1 * (4.0*PI * alpha * A2 * gsqrt) * W / (c.Gamma * M4);
    const cd inner = cd(h1/W2,0.0)*T1 + cd((M*c.Gamma*h2 - M2*h1)/W4,0.0)*T2;
    out.ord1.sigmaI = prefI1 * (2.0 * std::real(inner));

    // VP (optional)
    apply_vp(W, out.ord1.sigmaC, out.ord1.sigmaR, out.ord1.sigmaI, vp);

    // units
    out.ord1.sigma  = (out.ord1.sigmaC + out.ord1.sigmaR + out.ord1.sigmaI) * c.unit_conv;
    out.ord1.sigmaC *= c.unit_conv;
    out.ord1.sigmaR *= c.unit_conv;
    out.ord1.sigmaI *= c.unit_conv;

    out.ord1.n_beta = nb; out.ord1.n_2f1 = n2;
    out.ord1.conv_beta = I04.conv_beta && I6res.converged_romI && I6res.converged_23;
    out.ord1.conv_2f1  = I04.conv_2f1;
  }

  // ---------------- Order 2 ----------------
  {
    // σ_C,2 = q2 * 4π α^2 A^2 /(3 W^4) * I4
    out.ord2.sigmaC = q2 * (4.0*PI*alpha*alpha*A2 / (3.0*W4)) * I4;

    // σ_R,2 = q2 h0 * 6π A^2 ΓeeΓμμ /(Γ M^5) * (I1 + c.c.) = pref * 2Re(I1)
    const double prefR2 = q2 * h0 * (6.0*PI * A2 * gprod) / (c.Gamma * M5);
    out.ord2.sigmaR = prefR2 * (2.0 * std::real(I1));

    // σ_I,2 = q2 * 4π α A^2 sqrt(ΓeeΓμμ) * W^2 /(Γ M^4) * ( [h1/W^2 I1 + (MΓ h2 - M^2 h1)/W^4 * X2] + c.c.)
    // X2 = 1 - B/(A+B)^2 I2 + 1/(A+B) I0 + (B/(A+B)) I1
//     const cd X2 =
//     ((cd(1.0,0.0) - Bprop) / (AB*AB)) * cd(I2,0.0)
//   + cd(I0,0.0)/AB
//   + (Bprop/AB) * I1;
    const cd X2 =
    cd(I0,0.0)/AB
  + (Bprop/AB) * I1;

    const double prefI2 = q2 * (4.0*PI * alpha * A2 * gsqrt) * W2 / (c.Gamma * M4);
    const cd inner = cd(h1/W2,0.0)*I1 + cd((M*c.Gamma*h2 - M2*h1)/W4,0.0)*X2;
    out.ord2.sigmaI = prefI2 * (2.0 * std::real(inner));

    apply_vp(W, out.ord2.sigmaC, out.ord2.sigmaR, out.ord2.sigmaI, vp);

    out.ord2.sigma  = (out.ord2.sigmaC + out.ord2.sigmaR + out.ord2.sigmaI) * c.unit_conv;
    out.ord2.sigmaC *= c.unit_conv;
    out.ord2.sigmaR *= c.unit_conv;
    out.ord2.sigmaI *= c.unit_conv;

    out.ord2.n_beta = nb; out.ord2.n_2f1 = n2;
    out.ord2.conv_beta = I04.conv_beta && I6res.converged_romI && I6res.converged_23;
    out.ord2.conv_2f1  = I04.conv_2f1;
  }

  // ---------------- Order 3 ----------------
  {
    // σ_C,3 = q3 * 4π α^2 A^2 /(3 W^3) * I8
    out.ord3.sigmaC = q3 * (4.0*PI*alpha*alpha*A2 / (3.0*W3)) * I8;

    // σ_R,3 = q3 h0 * 6π A^2 ΓeeΓμμ * W /(Γ M^5) * (I6 + c.c.) = pref * 2Re(I6)
    const double prefR3 = q3 * h0 * (6.0*PI * A2 * gprod) * W / (c.Gamma * M5);
    out.ord3.sigmaR = prefR3 * (2.0 * std::real(I6));

    // σ_I,3 = q3 * 4π α A^2 sqrt(ΓeeΓμμ) * W^3 /(Γ M^4) * ( [h1/W^2 I6 + (MΓ h2 - M^2 h1)/W^4 * X3] + c.c.)
    // X3 = 1 - B/(A+B)^2 I9 + 1/(A+B) I7 + (B/(A+B)) I6
//     const cd X3 =
//     ((cd(1.0,0.0) - Bprop) / (AB*AB)) * cd(I9,0.0)
//   + cd(I7,0.0)/AB
//   + (Bprop/AB) * I6;
    const cd X3 =
    cd(I7,0.0)/AB
  + (Bprop/AB) * I6;

    const double prefI3 = q3 * (4.0*PI * alpha * A2 * gsqrt) * W3 / (c.Gamma * M4);
    const cd inner = cd(h1/W2,0.0)*I6 + cd((M*c.Gamma*h2 - M2*h1)/W4,0.0)*X3;
    out.ord3.sigmaI = prefI3 * (2.0 * std::real(inner));

    apply_vp(W, out.ord3.sigmaC, out.ord3.sigmaR, out.ord3.sigmaI, vp);

    out.ord3.sigma  = (out.ord3.sigmaC + out.ord3.sigmaR + out.ord3.sigmaI) * c.unit_conv;
    out.ord3.sigmaC *= c.unit_conv;
    out.ord3.sigmaR *= c.unit_conv;
    out.ord3.sigmaI *= c.unit_conv;

    out.ord3.n_beta = nb; out.ord3.n_2f1 = n2;
    out.ord3.conv_beta = I04.conv_beta && I6res.converged_romI && I6res.converged_23;
    out.ord3.conv_2f1  = I04.conv_2f1;
  }

  // Sum
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

} // namespace isr_sigma5pi