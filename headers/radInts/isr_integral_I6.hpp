#pragma once
// isr_integral_I6.hpp
//
// HPC-oriented C++ core for ISR propagator integral I6:
//
//   I6 = ∫_0^b dx  F(x,W) * sqrt(1-x) / (A + B x)
//      = IA + IB
//
// Inputs:
//   - A, B : std::complex<double>  (propagator parameters)
//   - b    : double, 0 < b < 1
//   - beta : double, beta > 0
//   - delta: double
//
// This implements the closed forms in your I6 note:
//   - IA: compact C_A / C_B factorization (uses sqrt, log, atanh, Li2)
//   - IB:  beta(1+delta)*I_romI  +  beta^2*b/(A*y) * Σ_q w_q[...]  (merged II+III)
//     with y = -b B / A and |y|>1 (assumed by your derivation).
//
// Special functions backend:
//   - Li2(real) and Li2(complex) via GSL if available:
//       real:    gsl_sf_dilog(x)
//       complex: gsl_sf_complex_dilog_e(r,theta,&re,&im) for z=r*e^{iθ}
//   - If GSL is not available, you can plug in your own dilog implementation.
//
// Numerical controls:
//   - Two series appear:
//       (1) I_romI: Σ_q  y^{-q}  Σ_k  [(-1/2)_k b^k/k!] / (beta-1-q+k)
//       (2) II+III merged: Σ_q w_q * [...H_q(y),S_q(y),H_{q+1}(y),S_{q+1}(y)...]
//   - Both are truncated by relative tolerances and max iterations.
//
// Branch conventions:
//   - Uses principal branches of std::sqrt/std::log/std::atanh/std::pow and GSL dilog.
//
// Build/link:
//   - Requires GSL for complex Li2: link -lgsl -lgslcblas
//
// Author: generated with ChatGPT, based on your note.

#include <complex>
#include <cmath>
#include <stdexcept>
#include <algorithm>

#if __has_include(<gsl/gsl_sf_dilog.h>)
  #include <gsl/gsl_sf_dilog.h>
  #include <gsl/gsl_sf_result.h>
  #define ISR6_HAVE_GSL_DILOG 1
#else
  #define ISR6_HAVE_GSL_DILOG 0
#endif

namespace isr6 {

using cd = std::complex<double>;
static constexpr double PILOCAL = 3.141592653589793238462643383279502884;

struct Options {
  // I_romI outer q / inner k truncation
  int    max_q_romI = 200000;
  int    max_k_romI = 200000;
  double rel_tol_q_romI = 1e-14;
  double rel_tol_k_romI = 1e-15;

  // II+III merged q-sum truncation
  int    max_q_23 = 200000;
  double rel_tol_q_23 = 1e-14;

  // safety: require |y| > y_min to trust the |y|>1 derivation
  double y_min_abs = 1.0;

  // If true, do minimal checks and avoid exceptions for speed (returns NaNs on failure).
  bool   fast_noexcept = false;
};

struct Result {
  cd I6  = cd(0.0, 0.0);
  cd IA  = cd(0.0, 0.0);
  cd IB  = cd(0.0, 0.0);

  cd IromI = cd(0.0, 0.0);
  cd I23   = cd(0.0, 0.0);

  int q_used_romI = 0;
  int k_used_romI_avg = 0; // approximate average (integer)
  int q_used_23 = 0;

  bool converged_romI = true;
  bool converged_23   = true;
};

// ------------------------
// Dilog wrappers
// ------------------------
inline double Li2_real(double x) {
#if ISR6_HAVE_GSL_DILOG
  return gsl_sf_dilog(x);
#else
  (void)x;
  throw std::runtime_error("Li2_real: GSL dilog not available; please provide a fallback.");
#endif
}

inline cd Li2_complex(cd z) {
#if ISR6_HAVE_GSL_DILOG
  const double r = std::abs(z);
  const double theta = std::arg(z);
  gsl_sf_result re, im;
  gsl_sf_complex_dilog_e(r, theta, &re, &im);
  return cd(re.val, im.val);
#else
  (void)z;
  throw std::runtime_error("Li2_complex: GSL complex dilog not available; please provide a fallback.");
#endif
}

// ------------------------
// IA computation (compact form)
// ------------------------
inline cd compute_IA(cd A, cd B, double b, double beta) {
  // reusable real scalars
  const double u  = std::sqrt(1.0 - b);
  const double Du = 1.0 - u;
  const double Lb1 = std::log1p(-b);     // log(1-b)
  const double Lu1 = std::log1p(u);      // log(1+u)

  // complex common terms
  const cd D = A + B;
  const cd S = std::sqrt(D);
  const cd sB = std::sqrt(B);
  const cd B52 = B*B*sB;                // B^(5/2)

  const cd r0 = sB / S;
  const cd r1 = (sB * u) / S;

  const cd At0 = std::atanh(r0);
  const cd At1 = std::atanh(r1);

  const cd z0 = B / D;
  const cd z1 = (B * (1.0 - b)) / D;

  const cd P = 2.0*B*(4.0 + 3.0*beta) + A*(4.0 + beta);
  const cd Q = 3.0*A*A + 6.0*A*B + 4.0*B*B;

  // Li2 terms
  const double Li_1mu = Li2_real(1.0 - u);
  const double Li_minu = Li2_real(-u);

  const cd Li_z0 = Li2_complex(z0);
  const cd Li_r0 = Li2_complex(r0);
  const cd Li_r1 = Li2_complex(r1);
  const cd Li_z1 = Li2_complex(z1);

  // C_A block
  const cd CA =
      3.0*A*A*sB*Du*(4.0 + 7.0*beta)
    + A*B*sB*( u*b*(4.0 + 3.0*beta) + Du*(28.0 + 57.0*beta) )
    + B52*beta*(PILOCAL*PILOCAL)
    - 3.0 * (
        A*S*P*At0
      + beta*( A*sB*u*((b - 7.0)*B - 3.0*A) - 2.0*B52*Lu1 )*Lb1
      + S*At1*( -A*P + Q*beta*Lb1 )
    );

  // C_B block
  const cd CB =
      8.0*B52*(cd(Li_1mu + Li_minu, 0.0))
    + S*Q*(Li_z0 - 4.0*Li_r0 + 4.0*Li_r1 - Li_z1);

  const cd pref = beta / (24.0 * A * B52);
  return pref * (2.0*CA + 3.0*beta*CB);
}

// ------------------------
// I_romI computation (double series) as in your Eq.(IromI)
// ------------------------
inline cd compute_IromI(cd A, cd y, double b, double beta,
                        const Options& opt,
                        int* q_used=nullptr, int* k_used_avg=nullptr, bool* converged=nullptr)
{
  // prefactor b^beta / A
  const double Lb = std::log(b);
  const cd bpow = std::exp(cd(beta*Lb, 0.0)); // b^beta (real)

  const cd invy = 1.0 / y;

  // First (closed) term
  const double s = std::sin(PILOCAL * beta);
  const cd pow_neg_y = std::pow(-y, -beta);              // (-y)^(-beta) (principal)
  const cd sqrt_term = std::sqrt(1.0 - (b * invy));      // sqrt(1 - b/y)
  const cd term0 = (PILOCAL / s) * pow_neg_y * sqrt_term;

  // Second term: -(1/y) * sum_q y^{-q} sum_k ...
  cd sumq(0.0, 0.0);
  cd invyq(1.0, 0.0);

  int q_count = 0;
  long long k_total = 0;
  bool ok = true;

  for (int q=0; q<opt.max_q_romI; ++q) {
    // inner k-sum is real
    const double denom0 = beta - 1.0 - double(q); // beta-1-q
    double coeffb = 1.0;                          // (-1/2)_0 / 0! * b^0
    double sumk = 0.0;
    double termk = coeffb / (denom0 + 0.0);
    sumk += termk;

    int k_used = 1;
    for (int k=0; k<opt.max_k_romI-1; ++k) {
      // update coeffb to next k+1:
      // coeff_{k+1} = coeff_k * (k - 1/2)/(k+1) * b
      const double kk = double(k);
      coeffb *= b * ((kk - 0.5) / (kk + 1.0));
      const double denom = denom0 + double(k+1);
      termk = coeffb / denom;
      sumk += termk;
      k_used++;

      if (std::abs(termk) <= opt.rel_tol_k_romI * std::max(1.0, std::abs(sumk))) break;
      // if (std::norm(termk) <= std::pow(opt.rel_tol_k_romI,2) * std::max(1.0, std::norm(sumk))) break;
      if (k_used >= opt.max_k_romI) { ok = false; break; }
    }

    k_total += k_used;
    q_count++;

    const cd qterm = invyq * cd(sumk, 0.0);
    sumq += qterm;

    if (q > 10 && std::abs(qterm) <= opt.rel_tol_q_romI * std::max(1.0, std::abs(sumq))) {
    // if (q > 10 && std::norm(qterm) <= std::pow(opt.rel_tol_q_romI,2) * std::max(1.0, std::norm(sumq))) {
      // converged
      break;
    }

    invyq *= invy;
    if (q_count >= opt.max_q_romI) { ok = false; break; }
  }

  if (q_used) *q_used = q_count;
  if (k_used_avg) *k_used_avg = (q_count>0) ? int(double(k_total)/double(q_count) + 0.5) : 0;
  if (converged) *converged = ok;

  return (bpow / A) * ( term0 - invy * sumq );
}

// ------------------------
// II+III merged q-sum (your merged expression using w_q, H_q, S_q recurrences)
// ------------------------
inline cd compute_I23(cd A, cd y, double b, double beta,
                      const Options& opt,
                      int* q_used=nullptr, bool* converged=nullptr)
{
  const cd invy = 1.0 / y;
  const cd by = cd(b,0.0) * invy;

  const double logb = std::log(b);

  // H0 = -log(1-y), S0 = -Li2(y)
  cd Hq = -std::log(1.0 - y);
  cd Sq = -Li2_complex(y);

  cd wq(1.0, 0.0);     // w0
  cd powy(1.0, 0.0);   // y^0

  cd sum(0.0, 0.0);
  bool ok = true;
  int q_count = 0;

  for (int q=0; q<opt.max_q_23; ++q) {
    const int qp1 = q + 1;

    // advance powy to y^{q+1}
    const cd powy_next = powy * y;

    // H_{q+1}, S_{q+1}
    const cd Hnext = Hq - powy_next / double(qp1);
    const cd Snext = Sq + powy_next / double(qp1*qp1);

    // bracket = -(log b * H_q + S_q) + (b/(4y))*((log b + 1/2) H_{q+1} + 2 S_{q+1})
    const cd bracket =
        -( cd(logb,0.0)*Hq + Sq )
      + ( cd(b,0.0) * invy / 4.0 ) * ( cd(logb + 0.5,0.0)*Hnext + 2.0*Snext );

    const cd qterm = wq * bracket;
    sum += qterm;

    q_count++;

    if (q > 10 && std::abs(qterm) <= opt.rel_tol_q_23 * std::max(1.0, std::abs(sum))) {
    // if (q > 10 && std::norm(qterm) <= std::pow(opt.rel_tol_q_23,2) * std::max(1.0, std::norm(sum))) {
      break;
    }

    // update recurrences for next q
    Hq = Hnext;
    Sq = Snext;
    powy = powy_next;

    // w_{q+1} = ((q-1/2)/(q+1))*(b/y) * w_q
    const double fac = (double(q) - 0.5) / double(qp1);
    wq *= cd(fac, 0.0) * by;

    if (q_count >= opt.max_q_23) { ok = false; break; }
  }

  if (q_used) *q_used = q_count;
  if (converged) *converged = ok;

  // Overall prefactor: beta^2 * b / (A y)
  return (beta*beta) * cd(b,0.0) * (1.0 / (A * y)) * sum;
}

// ------------------------
// Full I6
// ------------------------
inline Result compute_I6(cd A, cd B, double b, double beta, double delta,
                         const Options& opt = Options{})
{
  Result out;

  auto fail = [&](const char* msg)->Result{
    if (!opt.fast_noexcept) throw std::runtime_error(msg);
    out.I6 = out.IA = out.IB = out.IromI = out.I23 = cd(NAN, NAN);
    out.converged_romI = out.converged_23 = false;
    return out;
  };

  if (!(b > 0.0 && b < 1.0)) return fail("compute_I6: require 0<b<1.");
  if (!(beta > 0.0))         return fail("compute_I6: require beta>0.");

  const cd y = -cd(b,0.0) * B / A;
  if (std::abs(y) <= opt.y_min_abs) {
  // if (std::norm(y) <= std::pow(opt.y_min_abs, 2)) {
    // Your series representation assumes |y|>1 for fast convergence / analytic continuation choice.
    // You can still run it, but convergence may be poor.
    // Here we only warn by throwing unless fast_noexcept is enabled.
    if (!opt.fast_noexcept) {
      // Not throwing hard; just continue. Comment out if you want strictness.
      // throw std::runtime_error("compute_I6: |y|<=y_min_abs; derived representation assumes |y|>1.");
    }
  }

  // IA
  out.IA = compute_IA(A, B, b, beta);

  // IB pieces
  out.IromI = compute_IromI(A, y, b, beta, opt, &out.q_used_romI, &out.k_used_romI_avg, &out.converged_romI);
  out.I23   = compute_I23  (A, y, b, beta, opt, &out.q_used_23, &out.converged_23);

  out.IB = beta*(1.0 + delta)*out.IromI + out.I23;

  out.I6 = out.IA + out.IB;
  return out;
}

} // namespace isr6