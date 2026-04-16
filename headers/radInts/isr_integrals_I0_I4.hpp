#pragma once
// isr_integrals_I0_I4.hpp
//
// HPC-oriented C++ core implementing ISR integral basis I0..I4
// using analytic forms from:
//   Baoxin Liu, Zhenyu Zhang, Xiang Zhou, Phys. Rev. D 110, 053010 (2024)
//   Appendix B, Eqs. (B7)–(B11).
//
// Integrals:
//   I0 = ∫_0^b dx  F(x,W)/(1-x)                             (B7)
//   I1 = ∫_0^b dx  F(x,W)/(A + B x)        (A,B complex)     (B8)
//   I2 = ∫_0^b dx  F(x,W)                                     (B9)
//   I3 = ∫_0^b dx  F(x,W)/(1-x)^3                             (B10)
//   I4 = ∫_0^b dx  F(x,W)/(1-x)^2                             (B11)
//
// where
//   Li2(x) is Spence dilog, B(x;a,b) is incomplete beta,
//   and 2F1 is Gauss hypergeometric.
//
// Numerical strategy (HPC-friendly):
//   - Reuse scalars: log(b), log1p(-b), b^beta, Li2(b), Li2(1-b).
//   - Evaluate incomplete beta B(b;β;0/-1/-2) via fast power series in b (0<b<1):
//       B(b;β,0)  = Σ b^{β+n}/(β+n)
//       B(b;β,-1) = Σ (n+1) b^{β+n}/(β+n)              because (1-t)^(-2)=Σ(n+1)t^n
//       B(b;β,-2) = Σ (n+1)(n+2)/2 * b^{β+n}/(β+n)      because (1-t)^(-3)=Σ C(n+2,2)t^n
//   - Evaluate 2F1(1,β;β+1; -t) with t=bB/A using a large-argument connection:
//       2F1(1,β;β+1;-t) = β/(β-1) * (1/t) * 2F1(1,1-β;2-β; -1/t)
//                         + βπ/sin(πβ) * t^{-β}
//     and compute 2F1(1,1-β;2-β;z) by series in z (|z| small).
//
// Special functions backend:
//   - real Li2: gsl_sf_dilog(x)
//   - complex Li2: gsl_sf_complex_dilog_e(r,theta,&re,&im) for z=r e^{iθ}
//   - complex elementary functions: std::log, std::pow, etc.
//
// Link: -lgsl -lgslcblas

#include <complex>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <limits>

// #include <physicsFuncs.h>
#include <variables.h>

#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_errno.h>

namespace isr_i0i4 {

using cd = std::complex<double>;
// static constexpr double PI = 3.141592653589793238462643383279502884;

struct Options {
  // Incomplete beta series truncation (for B(b;β,0/-1/-2))
  int    max_n_beta = 200000;
  double rel_tol_beta = 1e-15;

  // 2F1 series truncation for 2F1(1,1-β;2-β;z) where |z| is expected small
  int    max_n_2f1 = 200000;
  double rel_tol_2f1 = 1e-15;

  bool   fast_noexcept = false;
};

struct Result {
  double I0 = 0.0;
  cd     I1 = cd(0.0, 0.0);
  double I2 = 0.0;
  double I3 = 0.0;
  double I4 = 0.0;

  // diagnostics
  int n_used_beta = 0;
  int n_used_2f1  = 0;
  bool conv_beta  = true;
  bool conv_2f1   = true;
};

// -------------------- Li2 wrappers --------------------
// inline double Li2_real(double x) {
//   return gsl_sf_dilog(1.-x) - PI*PI/6.;
// }
inline double Li2_real(double x) {
  return gsl_sf_dilog(x);
}

inline cd Li2_complex(cd z) {
  const double r = std::abs(z);
  const double theta = std::arg(z);
  gsl_sf_result re, im;
  gsl_sf_complex_dilog_e(r, theta, &re, &im);
  return cd(re.val, im.val);
}

// -------------------- Incomplete beta series: B(b;β,0/-1/-2) --------------------
struct BetaSeriesOut {
  double B0  = 0.0; // B(b;β,0)
  double Bm1 = 0.0; // B(b;β,-1)
  double Bm2 = 0.0; // B(b;β,-2)
  int n_used = 0;
  bool conv = true;
};

inline BetaSeriesOut beta_series_0m1m2(double b, double beta, const Options& opt) {
  BetaSeriesOut out;

  const double Lb = std::log(b);
  double pow_bn = std::exp(beta * Lb); // b^{beta+n} starting at n=0

  double S0=0.0, Sm1=0.0, Sm2=0.0;
  bool ok = true;

  for (int n=0; n<opt.max_n_beta; ++n) {
    const double denom = beta + double(n);

    const double t0  = pow_bn / denom;
    const double tm1 = (double(n)+1.0) * pow_bn / denom;
    const double tm2 = 0.5*(double(n)+1.0)*(double(n)+2.0) * pow_bn / denom;

    S0  += t0;
    Sm1 += tm1;
    Sm2 += tm2;

    out.n_used = n+1;

    const double scale = std::max({1.0, std::abs(S0), std::abs(Sm1), std::abs(Sm2)});
    const double tmax  = std::max({std::abs(t0), std::abs(tm1), std::abs(tm2)});

    if (n > 10 && tmax <= opt.rel_tol_beta * scale) break;

    pow_bn *= b;

    if (n+1 >= opt.max_n_beta) { ok = false; break; }
  }

  out.B0 = S0;
  out.Bm1 = Sm1;
  out.Bm2 = Sm2;
  out.conv = ok;
  return out;
}

// inline BetaSeriesOut beta_series_0m1m2(double b, double beta, const Options& opt) {
//     BetaSeriesOut out;
//     out.B0 = gsl_sf_hyperg_2F1(beta,1,beta+1,b)*pow(b,beta)/beta;
//     out.Bm1 = gsl_sf_hyperg_2F1(beta,2,beta+1,b)*pow(b,beta)/beta;
//     out.Bm2 = gsl_sf_hyperg_2F1(beta,3,beta+1,b)*pow(b,beta)/beta;
//     out.conv = true;
//     return out;
// }

// -------------------- 2F1 series for 2F1(1,1-β;2-β;z) --------------------
// 2F1(1,1-β;2-β;z) = Σ_{n>=0} ((1-β)_n / (2-β)_n) z^n
inline cd hyp2f1_1_1mb_2mb(cd z, double beta, const Options& opt,
                          int* n_used=nullptr, bool* conv=nullptr)
{
  const double bpar = 1.0 - beta;
  const double cpar = 2.0 - beta;

  cd sum(1.0, 0.0);
  cd term(1.0, 0.0);

  bool ok = true;
  int ncount = 0;

  for (int n=0; n<opt.max_n_2f1; ++n) {
    const double num = bpar + double(n);
    const double den = cpar + double(n);
    term *= (num/den) * z;
    sum += term;
    ncount = n+1;

    if (n > 10 && std::abs(term) <= opt.rel_tol_2f1 * std::max(1.0, std::abs(sum))) break;
    if (ncount >= opt.max_n_2f1) { ok = false; break; }
  }

  if (n_used) *n_used = ncount;
  if (conv) *conv = ok;
  return sum;
}

// -------------------- Compute I0, I2, I3, I4 (real) from (B7),(B9),(B10),(B11) --------------------
inline void compute_I0_I2_I3_I4(double b, double beta, double delta,
                                const BetaSeriesOut& Bs,
                                double Li2b, double Li2_1mb,
                                double Lb, double L1m,
                                double bpow,
                                Result& out)
{
  const double beta2 = beta*beta;
  const double PI2 = PI*PI;
  const double cdelta = 1.0 + delta;

  const double inv1mb  = 1.0/(1.0 - b);
  const double inv1mb2 = inv1mb*inv1mb;

  // ---------------- I0 (B7) ----------------
  out.I0 =
      beta*cdelta*Bs.B0
    - 0.5*b*beta
    + L1m*(0.25*beta2 + 0.5*beta + 0.375*beta2*b)
    + (beta2/16.0)*L1m*L1m
    - 0.5*beta2*b*Lb
    - 0.5*beta2*(Li2_1mb - PI2/6.0 - Li2b);

  // ---------------- I2 (B9) ----------------
  out.I2 =
      cdelta*bpow   // since β(1+δ) b^β / β
    + 0.5*beta2*Li2b
    + b*b*(beta/4.0 + beta2/32.0)
    + b*(-beta - 5.0*beta2/16.0)
    + L1m*( -9.0*beta2/16.0 + 0.75*beta2*b - 3.0*beta2*b*b/16.0 )
    + Lb*( beta2*(b*b/4.0 - b) );
  // I2 no longer needed
  // out.I2 = 0.0;

  // ---------------- I3 (B10) ----------------
  // (β^2/8 + β/2) * [-1 + (1-b)^3 + b(b+3-3b)] / [2(1-b)^3]
  // The numerator simplifies to b^2, so this becomes (β^2/8+β/2) * b^2 / [2(1-b)^2]
  const double tA = (beta/2.0 + beta2/8.0) * (b*b/2.0) * inv1mb2;

  const double tB = -0.5*(beta + 0.75*beta2) * (inv1mb2 - 1.0);

  // - (β^2/8) * [1-(1-b)^2 + 2 log(1-b)] / [4(1-b)^2]
  // 1-(1-b)^2 = b(2-b)
  const double tC = -(beta2/32.0) * ( (b*(2.0 - b) + 2.0*L1m) * inv1mb2 );

  const double tD = -0.75*beta2*L1m;

  // - (β^2/4) * b(b - 1 - (b-2)log b)/(1-b)^2
  // Implemented as +(β^2/4)* b( (1-b) + (b-2)log b )/(1-b)^2
  const double tE = (beta2/4.0) * ( b*((1.0 - b) + (b - 2.0)*Lb) * inv1mb2 );

  // + (β^2/2) * b log b/(b-1) = - (β^2/2) * b log b/(1-b)
  const double tF = -(beta2/2.0) * ( b*Lb*inv1mb );

  // - (1/2)β^2 ( -Li2(b) - 1/2 log^2(1-b) ) = + (β^2/2)Li2(b) + (β^2/4)log^2(1-b)
  const double tG = 0.5*beta2*Li2b + 0.25*beta2*L1m*L1m;

  // - (β/8) * (β b + log(1-b)) / (1-b)   [THIS IS IMPORTANT]??
  // -> - (β^2/8) * (b + log(1-b)) / (1-b)
  const double tH = -(beta2/8.0) * ( (b + L1m) * inv1mb );

  out.I3 = beta*cdelta*Bs.Bm2 + tA + tB + tC + tD + tE + tF + tG + tH;

  // ---------------- I4 (B11) ----------------
  out.I4 =
      beta*cdelta*Bs.Bm1
    + 0.5*beta2*(Li2b - Li2_1mb + PI2/6.0)
    + (b*inv1mb)*(-0.75*beta2 - 0.5*beta)
    + L1m*(0.5*beta - 0.375*beta2)
    - (b*Lb*beta2) * (0.5*inv1mb)
    + (beta2/16.0)*L1m*L1m
    - (beta2/8.0) * (L1m*inv1mb);
}

// -------------------- Compute I1 (complex) from (B8), with reused blocks --------------------
inline cd compute_I1(cd A, cd B, double b, double beta, double delta,
                     double Li2b, double L1m, double Lb, double bpow,
                     const Options& opt,
                     int* n_used_2f1=nullptr, bool* conv_2f1=nullptr)
{
  const double beta2 = beta*beta;
  const double cdelta = 1.0 + delta;

  const cd AB  = A + B;
  const cd ABb = A + B*cd(b,0.0);

  // L = log((A+Bb)/A) = log(1 + (B/A)b)
  const cd L = std::log(ABb / A);

  const cd logB_AB = std::log(B / AB);
  const cd logB_A  = std::log(B / A);

  // K = Li2((A+Bb)/(A+B)) - Li2(A/(A+B)) + log(B/(A+B))*log((A+Bb)/A)
  const cd K = (Li2_complex(ABb / AB) - Li2_complex(A / AB)) + logB_AB * L;

  // J = -π^2/6 + 1/2 log^2(1 + (B/A)b) + Li2(A/(A+Bb)) - log(B/A) log(1+(B/A)b)
  const cd J = cd(-PI*PI/6.0, 0.0) + 0.5*L*L + Li2_complex(A / ABb) - logB_A * L;

  // ---- 2F1(1,β;β+1; -bB/A) using connection formula ----
  // t = bB/A, so w = -t
  const cd t = cd(b,0.0) * B / A;
  const cd z_small = -1.0 / t;

  bool ok=true; int n2=0;
  const cd hyp_small = hyp2f1_1_1mb_2mb(z_small, beta, opt, &n2, &ok);

  if (n_used_2f1) *n_used_2f1 = n2;
  if (conv_2f1) *conv_2f1 = ok;

  const cd t_mbeta = std::pow(t, -beta);
  const double s = std::sin(PI*beta);

  // hyp = 2F1(1,β;β+1; -t)
  const cd hyp =
      (beta/(beta-1.0)) * (1.0/t) * hyp_small
    + (beta*PI/s) * t_mbeta;

  // First term of (B8): β(1+δ) * [b^β 2F1/(Aβ)] = (1+δ) b^β/A * 2F1
  const cd term_hyp = (cdelta / A) * cd(bpow,0.0) * hyp;

  const cd invB  = 1.0 / B;
  const cd invB2 = invB * invB;
  const cd invA  = 1.0 / A;

  const cd term1 = (beta2/8.0 + beta/2.0) * ( cd(b,0.0)*invB - A*invB2*L );
  const cd term2 = -(beta + 0.75*beta2) * invB * L;
  const cd term3 = 0.75*beta2 * (-invB) * K;
  const cd term4 = -(3.0/8.0)*beta2 * invB * ( cd(-b,0.0) + cd(b-1.0,0.0)*cd(L1m,0.0) + (A*invB)*K );
  const cd term5 = -beta2 * invB * J;
  const cd term6 = 0.5*beta2 * invB * ( cd(b*(Lb-1.0),0.0) - (A*invB)*J );
  const cd term7 = -0.5*beta2 * invA * ( cd(-Li2b,0.0) + K );

  return term_hyp + term1 + term2 + term3 + term4 + term5 + term6 + term7;
}

// -------------------- Public API: compute_all(I0..I4) --------------------
inline Result compute_all(double b, double beta, double delta, cd A, cd B,
                          const Options& opt = Options{})
{
  Result out;

  auto fail = [&](const char* msg)->Result{
    if (!opt.fast_noexcept) throw std::runtime_error(msg);
    out.I0 = out.I2 = out.I3 = out.I4 = std::numeric_limits<double>::quiet_NaN();
    out.I1 = cd(NAN, NAN);
    out.conv_beta = out.conv_2f1 = false;
    return out;
  };

  if (!(b > 0.0 && b < 1.0)) return fail("compute_all(I0..I4): require 0<b<1.");
  if (!(beta > 0.0))         return fail("compute_all(I0..I4): require beta>0.");

  const double Lb  = std::log(b);
  const double L1m = std::log1p(-b);
  const double bpow = std::exp(beta * Lb);

  // Real dilogs
  const double Li2b   = Li2_real(b);
  const double Li2_1m = Li2_real(1.0 - b);

  // Beta series
  const auto Bs = beta_series_0m1m2(b, beta, opt);
  out.n_used_beta = Bs.n_used;
  out.conv_beta   = Bs.conv;

  // I0/I2/I3/I4
  compute_I0_I2_I3_I4(b, beta, delta, Bs, Li2b, Li2_1m, Lb, L1m, bpow, out);

  // I1
  int n2=0; bool ok2=true;
  out.I1 = compute_I1(A, B, b, beta, delta, Li2b, L1m, Lb, bpow, opt, &n2, &ok2);
  out.n_used_2f1 = n2;
  out.conv_2f1   = ok2;

  return out;
}


// Wrapper to compute all using W,Mv,Gv,Wm instead of A,B,b,beta,delta
inline void make_params(double W, double Mv, double Gv, double Wm,
                    double& b, double& beta, double& delta, cd& A, cd& B)
{
  const double x = (W/Mv)*(W/Mv);    // (W/Mv)^2
  A = cd(Gv/Mv, x - 1.0);            // Gv/Mv + i(x-1)
  B = cd(0.0, -x);                   // -i x
  b = 1.0 - (Wm/W)*(Wm/W);
  beta = 2.0*alpha/PI * (2.0*std::log(W/me) - 1.0);
  delta = 0.75*beta + alpha/PI*(PI*PI/3.0 - 0.5) + beta*beta*(9.0/32.0 - PI*PI/12.0);
}

inline Result compute_all(double W, double Mv, double Gv, double Wm,
                          const Options& opt = Options{})
{
  cd A,B; double b,beta,delta;
  make_params(W, Mv, Gv, Wm, b, beta, delta, A, B);
  return compute_all(b, beta, delta, A, B, opt);
}

} // namespace isr_i0i4