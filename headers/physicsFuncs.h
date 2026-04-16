#ifndef PHYSICSFUNCS_H
#define PHYSICSFUNCS_H

#include <iostream>
#include <fstream>
#include "TComplex.h"
#include <complex>

#include "TPaveText.h"
#include <complex>
#include <iostream>
#include <fstream>
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "Math/QuantFuncMathCore.h"

#include <chrono>

#include <RooMath.h>
#include <iomanip>
#include <iostream>
#include <TTree.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TMinuit.h>
#include <TComplex.h>
#include <TLorentzVector.h>
#include <TGraph.h>

#include <cassert>


#include "Riostream.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TGraphErrors.h"
#include "TDecompChol.h"
#include "TDecompSVD.h"
#include "TF1.h"
#include  "TVirtualPad.h"
#include "TCanvas.h"
// #include "bes3plotstyle.h"
// #include "bes3plotstyle.c"
#include <Math/Interpolator.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_integration.h>

#include "variables.h"
#include "HyperGeometrics.h"

#include "sigma5pi_isr_qf3.hpp"    // namespace isr_sigma5pi
// #include "isr_integrals_I0_I4.hpp"
// #include "isr_integral_I6.hpp"
#include "vaccFunc.hpp"
#include "sigma5pi_isr.hpp"

// W^-8 dependence
// #include "sigma2body_isr_qf3.hpp"
#include "sigma2body_isr_qf5.hpp"

// W^-7 dependence
#include "sigma2bodyW7_isr_qf5.hpp"

// W^-6 dependence
#include "sigma2bodyW6_isr_qf5.hpp"

// Vaccum polarization correction
isr_sigma5pi::VPOptions vpLocal;


// GH integration settings
namespace gh128 {

  static bool inited = false;
  static double nodes[128];
  static double weights[128];

  inline void init_once() {
    if (inited) return;

    constexpr int N = 128;
    // Jacobi matrix for orthonormal Hermite polynomials (weight exp(-x^2)):
    // diagonal a_k = 0
    // off-diagonal b_k = sqrt((k+1)/2), k=0..N-2
    TMatrixDSym J(N);
    for (int i = 0; i < N; ++i) J(i,i) = 0.0;
    for (int k = 0; k < N-1; ++k) {
      const double b = std::sqrt((k+1) / 2.0);
      J(k, k+1) = b;
      J(k+1, k) = b;
    }

    TMatrixDSymEigen eig(J);
    TVectorD eval = eig.GetEigenValues();
    TMatrixD evec = eig.GetEigenVectors(); // columns are eigenvectors

    // Weights for ∫ exp(-x^2) f(x) dx are:
    //   w_i = μ0 * (v0_i)^2,  μ0 = ∫ exp(-x^2) dx = sqrt(pi)
    // where v0_i is the first component of normalized eigenvector i.
    const double mu0 = std::sqrt(M_PI);

    // Sort by eigenvalue increasing (safer across ROOT versions)
    std::vector<int> idx(N);
    for (int i=0;i<N;++i) idx[i]=i;
    std::sort(idx.begin(), idx.end(), [&](int a, int b){ return eval[a] < eval[b]; });

    for (int j=0; j<N; ++j) {
      const int i = idx[j];
      const double x = eval[i];
      const double v0 = evec(0, i);        // first component
      const double w  = mu0 * v0 * v0;
      nodes[j] = x;
      weights[j] = w;
    }

    inited = true;
  }

} // namespace gh128



ROOT::Math::Interpolator interpolator(ROOT::Math::Interpolation::kLINEAR);
using namespace std;

double vacc(double W)
{
	    return interpolator.Eval(W);
}

complex<double> LI2(complex<double> z)
{
    int k=1;
    double eps = 1;
    complex<double> result = {0,0};
    complex<double> temp = {0,0};
    for(int i=0;i<500;i++)
    {
        // temp = result;
        temp = pow(z,k)/pow((double)k,2)/pow((double)k+1.,2)/pow((double)k+2.,2);
        result += temp;
        k++;
        if(fabs(temp)<=TOLS*fabs(result))
        {
            break;
        }
    }
    return (4.*pow(z,2)*result+ 4.*z + 23./4.*z*z + 3.*(1.-pow(z,2))*log(1.-z))/(1.+4.*z+pow(z,2));
}

double Beta(double x, double a, double b)
{
    return gsl_sf_hyperg_2F1(a,1-b,a+1,x)*pow(x,a)/a;

}

double Spencer(double z)
{
    return gsl_sf_dilog(1.-z)-PI*PI/6.;
}
// double Spencer(double z)
// {
//     return gsl_sf_dilog(z);
// }
// double Spencer(double z)
// {
//     return gsl_sf_dilog(1.-z);
// }

// complex<double> Spencec(complex<double> z)
// {
//     return LI2(1.-z)-PI*PI/6.;
// }
complex<double> Spencec(complex<double> z)
{
    const double r = std::abs(z);
    const double theta = std::arg(z);
    gsl_sf_result re, im;
    gsl_sf_complex_dilog_e(r, theta, &re, &im);
    return complex<double>(re.val, im.val);
}


#include "additionalFunctions.h"



double Ana_legacy(double x, double M,double CC1, double CC2, double FF,double phi1,double phi2,double GammaVar,double Gamma_eeVar)
{
    double W = x;
    // double M = par[0];
    // double phi = par[1];
    // double Vacc = vacc(W);
    double Vacc = 1.0;                                  // Turn off vacuum polarization for test
    double fGam = GammaVar/M;
    double bbeta = 2*alpha/PI*(2*log(W/me)-1);
    double ddelta = 3./4.*bbeta+alpha/PI*(PI*PI/3.-1./2.)+pow(bbeta,2)*(9./32.-PI*PI/12.);
    double bb = 1-pow(Wm/W,2);
    double ff = W/M;
    complex<double> AA = {fGam, pow(ff,2)-1};
    complex<double> BB = {0, -pow(ff,2)};
    //Fitting parameter
    double qf = sqrt(pow(W*W-mo*mo-mp*mp,2)-4*mp*mp*mo*mo)/2/W;
    // double qf = sqrt(pow(M*M-mo*mo-mp*mp,2)-4*mp*mp*mo*mo)/2/M;
    // double qf = 1;              // Turn off qf for test
    // double qf = ;
    // double qf = sqrt(1- 4 * mo * mp / pow(M, 2));
    // double qf = 1;
    // double C = pow(FF,2)/3.*4.*PI*alpha*alpha;
    double C = pow(FF,2)/3.*4.*PI*alpha*alpha*pow(qf,3);
    double D = 3*Gamma_eeVar/alpha/M;
    double H1 = CC1*(cos(phi1)+CC2*cos(phi1)*cos(phi2)-CC2*sin(phi1)*sin(phi2));
    double H2 = CC1*(sin(phi1)+CC2*cos(phi1)*sin(phi2)+CC2*sin(phi1)*cos(phi2));
    double C1 = C*D*D*(H1*H1+H2*H2)+2*D*C*H1*sqrt(Vacc);
    double C2 = 2*C*D*( H2*M*GammaVar-H1*M*M )*sqrt(Vacc);
    //hyp2f1 part
    complex<double> hyp = bbeta/(bbeta-1.)/(bb*BB/AA)*hyp2f1(1, 1-bbeta,2-bbeta,1./(-bb*BB/AA)) + bbeta*PI/sin(bbeta*PI)*pow(bb*BB/AA,-bbeta);
    //F1 part
    double p3F1 = Beta(bb,bbeta,-2);
    double p2F1 = Beta(bb,bbeta,-1);
    double p1F1 = Beta(bb,bbeta,0);
    complex<double> p0F1 = pow(bb,(bbeta))*hyp/AA/bbeta;
    //F2 x part
    double p3F21 = pow((1.-bb),(-3))*(-1+pow((1.-bb),3)+bb*(bb+3.-bb*3))/2.;
    double p2F21 = bb/(1.-bb)+log(1.-bb);
    double p1F21 = -bb-log(1.-bb);
    complex<double> p0F21 = bb/BB-AA/BB/BB*log((AA+BB*bb)/AA);
    //F2 constant part
    double p3F22 = 1./2.*pow((1.-bb),(-2))-1./2.;
    double p2F22 = 1./(1.-bb)-1.;
    double p1F22 = -log(1.-bb);
    complex<double> p0F22 = 1./BB*log((AA+BB*bb)/AA);
    //F2 ln(1-x) part
    double p3F23 = pow((1.-bb),(-3))*(1.-pow((1.-bb),3)-bb+2.*(1.-bb)*log(1.-bb))/4.;
    double p2F23 = pow((1.-bb),(-2))*(1.-pow((1.-bb),2)-bb+(1.-bb)*log(1.-bb));
    double p1F23 = -1./2.*pow(log(1.-bb),2);
    complex<double> p0F23 = -1./BB*( Spencec(BB/(AA+BB)*(1.-bb))-Spencec(BB/(AA+BB))+log(BB/(AA+BB))*log((AA+BB*bb)/AA) );
    //F2 xln(1-x) part
    double p3F24 = -pow((1.-bb),(-2))*(1.-pow((1.-bb),2)-bb+(1.-bb)*log(1.-bb)) + p3F23;
    double p2F24 = 1./2.*pow(log(1.-bb),2)+p2F23;
    double p1F24 = -(-bb+(bb-1.)*log(1.-bb)) + p1F23;
    complex<double> p0F24 = 1./BB*( -bb+(bb-1.)*log(1.-bb)-AA*p0F23 );
    //F3 lnx part
    double p3F31 = 1./2.*(log(1.-bb)+(bb*(bb-1.+(2.-bb)*log(bb))/pow((1.-bb),2)));
    double p2F31 = log(1.-bb)-bb*log(bb)/(bb-1.);
    double p1F31 = Spencer(bb);
    complex<double> p0F31 = 1./BB*( -PI*PI/6.+1./2.*pow(log(1.+BB/AA*bb),2)+Spencec(BB*bb/(AA+BB*bb))-Spencer(1)-log(BB/AA)*log(1.+BB/AA*bb));
    //F3 xlnx part
    double p3F32 = -log(1.-bb) + bb*log(bb)/(bb-1.) + 1./2.*(log(1.-bb)+bb*(bb-1.-(bb-2.)*log(bb))/pow(1-bb,2));
    double p2F32 = -Spencer(bb)+log(1.-bb)-bb*log(bb)/(bb-1.);
    double p1F32 = bb*(1.-log(bb))+Spencer(bb);
    complex<double> p0F32 = 1./BB*( bb*(log(bb)-1.) - AA/BB*( -PI*PI/6.+1./2.*pow(log(1.+BB/AA*bb),2)+Spencec(BB*bb/(AA+BB*bb))-Spencer(1)-log(BB/AA)*log(1.+BB/AA*bb))  );
    //F3 ln(1-x)/x part
    double p3F33 = p3F23 + p2F23 + p1F23 + Spencer(1.) - Spencer(1.-bb);
    double p2F33 = p2F23 + p1F23 + Spencer(1.) - Spencer(1.-bb);
    double p1F33 = p1F23 + Spencer(1.) - Spencer(1.-bb);
    complex<double> p0F33 = 1./AA*(Spencer(1)-Spencer(1.-bb)) + 1./AA*(Spencec(BB/(AA+BB)*(1.-bb))-Spencec(BB/(AA+BB))+log(BB/(AA+BB))*log((AA+BB*bb)/AA));

    //All of them
    double anaIII = bbeta*(1+ddelta)*p3F1 + (bbeta*bbeta/8.+bbeta/2.)*p3F21 - (bbeta+3./4.*bbeta*bbeta)*p3F22 + 3./4.*bbeta*bbeta*p3F23 - 3./8.*bbeta*bbeta*p3F24 - bbeta*bbeta*p3F31 + 1./2.*bbeta*bbeta*p3F32 - 1./2.*bbeta*bbeta*p3F33;
    double anaII = bbeta*(1+ddelta)*p2F1 + (bbeta*bbeta/8.+bbeta/2.)*p2F21 - (bbeta+3./4.*bbeta*bbeta)*p2F22 + 3./4.*bbeta*bbeta*p2F23 - 3./8.*bbeta*bbeta*p2F24 - bbeta*bbeta*p2F31 + 1./2.*bbeta*bbeta*p2F32 - 1./2.*bbeta*bbeta*p2F33;
    double anaI = bbeta*(1+ddelta)*p1F1 + (bbeta*bbeta/8.+bbeta/2.)*p1F21 - (bbeta+3./4.*bbeta*bbeta)*p1F22 + 3./4.*bbeta*bbeta*p1F23 - 3./8.*bbeta*bbeta*p1F24 - bbeta*bbeta*p1F31 + 1./2.*bbeta*bbeta*p1F32 - 1./2.*bbeta*bbeta*p1F33;
    complex<double> ana0 = bbeta*(1+ddelta)*p0F1 + (bbeta*bbeta/8.+bbeta/2.)*p0F21 - (bbeta+3./4.*bbeta*bbeta)*p0F22 + 3./4.*bbeta*bbeta*p0F23 - 3./8.*bbeta*bbeta*p0F24 - bbeta*bbeta*p0F31 + 1./2.*bbeta*bbeta*p0F32 - 1./2.*bbeta*bbeta*p0F33;

    //parameters
    complex<double> k21 = 1./(AA+BB);
    complex<double> k22 = BB/(AA+BB);
    complex<double> k41 = 1./(AA+BB);
    complex<double> k42 = BB/pow((AA+BB),2);
    complex<double> k43 = BB*BB/pow(AA+BB,2);
    //All kinds of calculation
    double con = anaIII*Vacc;
    double bw2 = (k21*anaI + k22*ana0).real()/fGam;
    double bw4 = (k41*anaII+ k42*anaI + k43*ana0).real()/fGam;

    return (C*con/pow(ff,6)/pow(M,6) + C1/pow(ff,2)/pow(M,6)*bw2 + C2/pow(M,8)/pow(ff,4)*bw4)*hbar_c_2;
}

double Ana_legacy_component(double x, double M,double CC1, double CC2, double FF,double phi1,double phi2,double GammaVar,double Gamma_eeVar, TString component = "total")
{
    double W = x;
    // double M = par[0];
    // double phi = par[1];
    // double Vacc = vacc(W);
    double Vacc = 1.0;
    double fGam = GammaVar/M;
    double bbeta = 2*alpha/PI*(2*log(W/me)-1);
    double ddelta = 3./4.*bbeta+alpha/PI*(PI*PI/3.-1./2.)+pow(bbeta,2)*(9./32.-PI*PI/12.);
    double bb = 1-pow(Wm/W,2);
    double ff = W/M;
    complex<double> AA = {fGam, pow(ff,2)-1};
    complex<double> BB = {0, -pow(ff,2)};
    //Fitting parameter
    // double qf = sqrt(pow(W*W-mo*mo-mp*mp,2)-4*mp*mp*mo*mo)/2/W;
    // double qf = sqrt(pow(M*M-mo*mo-mp*mp,2)-4*mp*mp*mo*mo)/2/M;
    // double qf = ;
    // double qf = sqrt(1- 4 * mo * mp / pow(M, 2));
    double qf = 1;
    // double C = pow(FF,2)/3.*4.*PI*alpha*alpha;
    double C = pow(FF,2)/3.*4.*PI*alpha*alpha*pow(qf,3);
    double D = 3*Gamma_eeVar/alpha/M;
    double H1 = CC1*(cos(phi1)+CC2*cos(phi1)*cos(phi2)-CC2*sin(phi1)*sin(phi2));
    double H2 = CC1*(sin(phi1)+CC2*cos(phi1)*sin(phi2)+CC2*sin(phi1)*cos(phi2));
    double C1 = C*D*D*(H1*H1+H2*H2)+2*D*C*H1*sqrt(Vacc);
    double C2 = 2*C*D*( H2*M*GammaVar-H1*M*M )*sqrt(Vacc);
    //hyp2f1 part
    complex<double> hyp = bbeta/(bbeta-1.)/(bb*BB/AA)*hyp2f1(1, 1-bbeta,2-bbeta,1./(-bb*BB/AA)) + bbeta*PI/sin(bbeta*PI)*pow(bb*BB/AA,-bbeta);
    //F1 part
    double p3F1 = Beta(bb,bbeta,-2);
    double p2F1 = Beta(bb,bbeta,-1);
    double p1F1 = Beta(bb,bbeta,0);
    complex<double> p0F1 = pow(bb,(bbeta))*hyp/AA/bbeta;
    //F2 x part
    double p3F21 = pow((1.-bb),(-3))*(-1+pow((1.-bb),3)+bb*(bb+3.-bb*3))/2.;
    double p2F21 = bb/(1.-bb)+log(1.-bb);
    double p1F21 = -bb-log(1.-bb);
    complex<double> p0F21 = bb/BB-AA/BB/BB*log((AA+BB*bb)/AA);
    //F2 constant part
    double p3F22 = 1./2.*pow((1.-bb),(-2))-1./2.;
    double p2F22 = 1./(1.-bb)-1.;
    double p1F22 = -log(1.-bb);
    complex<double> p0F22 = 1./BB*log((AA+BB*bb)/AA);
    //F2 ln(1-x) part
    double p3F23 = pow((1.-bb),(-3))*(1.-pow((1.-bb),3)-bb+2.*(1.-bb)*log(1.-bb))/4.;
    double p2F23 = pow((1.-bb),(-2))*(1.-pow((1.-bb),2)-bb+(1.-bb)*log(1.-bb));
    double p1F23 = -1./2.*pow(log(1.-bb),2);
    complex<double> p0F23 = -1./BB*( Spencec(BB/(AA+BB)*(1.-bb))-Spencec(BB/(AA+BB))+log(BB/(AA+BB))*log((AA+BB*bb)/AA) );
    //F2 xln(1-x) part
    double p3F24 = -pow((1.-bb),(-2))*(1.-pow((1.-bb),2)-bb+(1.-bb)*log(1.-bb)) + p3F23;
    double p2F24 = 1./2.*pow(log(1.-bb),2)+p2F23;
    double p1F24 = -(-bb+(bb-1.)*log(1.-bb)) + p1F23;
    complex<double> p0F24 = 1./BB*( -bb+(bb-1.)*log(1.-bb)-AA*p0F23 );
    //F3 lnx part
    double p3F31 = 1./2.*(log(1.-bb)+(bb*(bb-1.+(2.-bb)*log(bb))/pow((1.-bb),2)));
    double p2F31 = log(1.-bb)-bb*log(bb)/(bb-1.);
    double p1F31 = Spencer(bb);
    complex<double> p0F31 = 1./BB*( -PI*PI/6.+1./2.*pow(log(1.+BB/AA*bb),2)+Spencec(BB*bb/(AA+BB*bb))-Spencer(1)-log(BB/AA)*log(1.+BB/AA*bb));
    //F3 xlnx part
    double p3F32 = -log(1.-bb) + bb*log(bb)/(bb-1.) + 1./2.*(log(1.-bb)+bb*(bb-1.-(bb-2.)*log(bb))/pow(1-bb,2));
    double p2F32 = -Spencer(bb)+log(1.-bb)-bb*log(bb)/(bb-1.);
    double p1F32 = bb*(1.-log(bb))+Spencer(bb);
    complex<double> p0F32 = 1./BB*( bb*(log(bb)-1.) - AA/BB*( -PI*PI/6.+1./2.*pow(log(1.+BB/AA*bb),2)+Spencec(BB*bb/(AA+BB*bb))-Spencer(1)-log(BB/AA)*log(1.+BB/AA*bb))  );
    //F3 ln(1-x)/x part
    double p3F33 = p3F23 + p2F23 + p1F23 + Spencer(1.) - Spencer(1.-bb);
    double p2F33 = p2F23 + p1F23 + Spencer(1.) - Spencer(1.-bb);
    double p1F33 = p1F23 + Spencer(1.) - Spencer(1.-bb);
    complex<double> p0F33 = 1./AA*(Spencer(1)-Spencer(1.-bb)) + 1./AA*(Spencec(BB/(AA+BB)*(1.-bb))-Spencec(BB/(AA+BB))+log(BB/(AA+BB))*log((AA+BB*bb)/AA));

    //All of them
    double anaIII = bbeta*(1+ddelta)*p3F1 + (bbeta*bbeta/8.+bbeta/2.)*p3F21 - (bbeta+3./4.*bbeta*bbeta)*p3F22 + 3./4.*bbeta*bbeta*p3F23 - 3./8.*bbeta*bbeta*p3F24 - bbeta*bbeta*p3F31 + 1./2.*bbeta*bbeta*p3F32 - 1./2.*bbeta*bbeta*p3F33;
    double anaII = bbeta*(1+ddelta)*p2F1 + (bbeta*bbeta/8.+bbeta/2.)*p2F21 - (bbeta+3./4.*bbeta*bbeta)*p2F22 + 3./4.*bbeta*bbeta*p2F23 - 3./8.*bbeta*bbeta*p2F24 - bbeta*bbeta*p2F31 + 1./2.*bbeta*bbeta*p2F32 - 1./2.*bbeta*bbeta*p2F33;
    double anaI = bbeta*(1+ddelta)*p1F1 + (bbeta*bbeta/8.+bbeta/2.)*p1F21 - (bbeta+3./4.*bbeta*bbeta)*p1F22 + 3./4.*bbeta*bbeta*p1F23 - 3./8.*bbeta*bbeta*p1F24 - bbeta*bbeta*p1F31 + 1./2.*bbeta*bbeta*p1F32 - 1./2.*bbeta*bbeta*p1F33;
    complex<double> ana0 = bbeta*(1+ddelta)*p0F1 + (bbeta*bbeta/8.+bbeta/2.)*p0F21 - (bbeta+3./4.*bbeta*bbeta)*p0F22 + 3./4.*bbeta*bbeta*p0F23 - 3./8.*bbeta*bbeta*p0F24 - bbeta*bbeta*p0F31 + 1./2.*bbeta*bbeta*p0F32 - 1./2.*bbeta*bbeta*p0F33;

    //parameters
    complex<double> k21 = 1./(AA+BB);
    complex<double> k22 = BB/(AA+BB);
    complex<double> k41 = 1./(AA+BB);
    complex<double> k42 = BB/pow((AA+BB),2);
    complex<double> k43 = BB*BB/pow(AA+BB,2);
    //All kinds of calculation
    double con = anaIII*Vacc;
    double bw2 = (k21*anaI + k22*ana0).real()/fGam;
    double bw4 = (k41*anaII+ k42*anaI + k43*ana0).real()/fGam;

    if(component == "total"){
        return (C*con/pow(ff,6)/pow(M,6) + C1/pow(ff,2)/pow(M,6)*bw2 + C2/pow(M,8)/pow(ff,4)*bw4)*hbar_c_2;
    }
    else if(component == "con"){
        return (C*con/pow(ff,6)/pow(M,6))*hbar_c_2;
    }
    else if(component == "res"){
        return (C1/pow(ff,2)/pow(M,6)*bw2)*hbar_c_2;
    }
    else if(component == "int"){
        return (C2/pow(M,8)/pow(ff,4)*bw4)*hbar_c_2;
    }

    else if(component == "I0"){
        return anaI;
    }
    else if(component == "I4"){
        return anaII;
    }
    else if(component == "I3"){
        return anaIII;
    }
    else if(component == "I1real"){
        return ana0.real();
    }
    else if(component == "I1imag"){
        return ana0.imag();
    }
    else{
        cout<<"No such component!"<<endl;
        return 0;
    }
}


double BrCalc(double CC1, double FF) {

    // use sqrt instead of pow(...,0.5)
    // double num = MJpsi_mean*MJpsi_mean - mo*mo - mp*mp;
    // double discriminant = num*num - 4.0 * mp*mp * mo*mo;
    // double qf = sqrt(discriminant) / (2.0 * MJpsi_mean);
    double qf = sqrt(pow(1 + (mo*mo - mp*mp)/pow(MJpsi_mean, 2), 2) - (4 * mo*mo)/pow(MJpsi_mean, 2));

    // square and cube explicitly, or via std::pow
    double CC1_sq = CC1 * CC1;
    double FF_sq  = FF  * FF;
    double qf_cu  = qf   * qf * qf;

    // build the branching ratio
    // double result = CC1_sq * FF_sq * Gamma_ee_mean * qf_cu
    //               / (pow(MJpsi_mean, 4) * Gamma_mean);
    double result = CC1_sq * FF_sq * Gamma_ee_mean * qf_cu
                  / (8 * pow(MJpsi_mean, 6) * Gamma_mean);              // W^-8

    return result;
}

double BrCalc(double CC1, double FF, double M, double Gamma, double Gamma_ee) {

    // use sqrt instead of pow(...,0.5)
    // double num = M*M - mo*mo - mp*mp;
    // double discriminant = num*num - 4.0 * mp*mp * mo*mo;
    // double qf = sqrt(discriminant) / (2.0 * M);
    double qf = sqrt(pow(1 + (mo*mo - mp*mp)/pow(M, 2), 2) - (4 * mo*mo)/pow(M, 2));

    // square and cube explicitly, or via std::pow
    double CC1_sq = CC1 * CC1;
    double FF_sq  = FF  * FF;
    double qf_cu  = qf   * qf * qf;

    // build the branching ratio
    // double result = CC1_sq * FF_sq * Gamma_ee * qf_cu
    //               / (pow(M, 4) * Gamma);

    // W^-8
    double result = CC1_sq * FF_sq * Gamma_ee * qf_cu
                  / (8 * pow(M, 6) * Gamma);     

    // W^-7
    // double result = CC1_sq * FF_sq * Gamma_ee * qf_cu
    //               / (8 * pow(M, 5) * Gamma);

    // W^-6
    // double result = CC1_sq * FF_sq * Gamma_ee * qf_cu
    //               / (8 * pow(M, 4) * Gamma);  

    return result;
}

double BrCalcEM(double FF, double M, double Gamma, double Gamma_ee) {

    // use sqrt instead of pow(...,0.5)
    // double num = M*M - mo*mo - mp*mp;
    // double discriminant = num*num - 4.0 * mp*mp * mo*mo;
    // double qf = sqrt(discriminant) / (2.0 * M);
    double qf = sqrt(pow(1 + (mo*mo - mp*mp)/pow(M, 2), 2) - (4 * mo*mo)/pow(M, 2));

    // square and cube explicitly, or via std::pow
    double FF_sq  = FF  * FF;
    double qf_cu  = qf   * qf * qf;

    // build the branching ratio
    // double result = FF_sq * Gamma_ee * qf_cu
    //               / (pow(M, 4) * Gamma);

    // W^-8
    double result = FF_sq * Gamma_ee * qf_cu
                  / (8 * pow(M, 6) * Gamma);

    // W^-7
    // double result = FF_sq * Gamma_ee * qf_cu
    //               / (8 * pow(M, 5) * Gamma);

    // W^-6
    // double result = FF_sq * Gamma_ee * qf_cu
    //               / (8 * pow(M, 4) * Gamma);

    return result;
}

double CC1fromCC2(double CC2, double phi2) {
    std::complex<double> Ceiphi = std::complex<double>(1.0, 0.0) + CC2 * std::polar(1.0, phi2);
    return std::abs(Ceiphi);
}

// The ISR function actually used in the fit
double Ana(double W, double M,double CC1, double CC2, double FF,double phi1,double phi2,double GammaVar,double Gamma_eeVar)
{
    isr_sigma5pi::Consts consts;
    isr_sigma2body::Consts consts2body;
    consts.Gamma = GammaVar;
    consts.Gamma_ee = Gamma_eeVar;
    consts2body.Gamma = GammaVar;
    consts2body.Gamma_ee = Gamma_eeVar;
    // auto calcRes = isr_sigma5pi::sigma5pi_ISR_qf3(W, M, FF, CC1, CC2, phi1, phi2, q0Corr, q1Corr, q2Corr, q3Corr, consts, isr_i0i4::Options{}, isr6::Options{}, vpLocal);
    // auto calcRes = isr_sigma2body::sigma2body_ISR_qf3(W, M, FF, CC1, CC2, phi1, phi2, q0Corr, q1Corr, q2Corr, q3Corr, consts2body, isr_i0i11::Options{}, vpLocal);
    
    // W^-8
    auto calcRes = isr_sigma2body::sigma2body_ISR_qf5(W, M, FF, CC1, CC2, phi1, phi2, q0Corr_2body, q1Corr_2body, q2Corr_2body, q3Corr_2body, q4Corr_2body, q5Corr_2body, consts2body, isr_i0i11::Options{}, vpLocal);
    
    // W^-7
    // auto calcRes = isr_sigma2bodyW7::sigma2bodyW7_ISR_qf5(W, M, FF, CC1, CC2, phi1, phi2, q0Corr_2body, q1Corr_2body, q2Corr_2body, q3Corr_2body, q4Corr_2body, q5Corr_2body, consts2body, isr_i0i11::Options{}, vpLocal);

    // W^-6
    // auto calcRes = isr_sigma2bodyW6::sigma2bodyW6_ISR_qf5(W, M, FF, CC1, CC2, phi1, phi2, q0Corr_2body, q1Corr_2body, q2Corr_2body, q3Corr_2body, q4Corr_2body, q5Corr_2body, consts2body, isr_i0i11::Options{}, vpLocal);

    double sigma = calcRes.sigma;
    return sigma;
}

double Ana_component(double W, double M,double CC1, double CC2, double FF,double phi1,double phi2,double GammaVar,double Gamma_eeVar,TString component = "total")
{
    isr_sigma5pi::Consts consts;
    consts.Gamma = GammaVar;
    consts.Gamma_ee = Gamma_eeVar;
    // auto calcRes = isr_sigma5pi::sigma5pi_ISR_qf3(W, M, FF, CC1, CC2, phi1, phi2, q0Corr, q1Corr, q2Corr, q3Corr, consts, isr_i0i4::Options{}, isr6::Options{}, vpLocal);
    isr_sigma2body::Consts consts2body;
    auto calcRes = isr_sigma2body::sigma2body_ISR_qf5(W, M, FF, CC1, CC2, phi1, phi2, q0Corr_2body, q1Corr_2body, q2Corr_2body, q3Corr_2body, q4Corr_2body, q5Corr_2body, consts2body, isr_i0i11::Options{}, vpLocal);
    if(component == "total"){
        return calcRes.sigma;
    }
    else if(component == "con"){
        return calcRes.sigmaC;
    }
    else if(component == "res"){
        return calcRes.sigmaR;
    }
    else if(component == "int"){
        return calcRes.sigmaI;
    }
    else{
        cout<<"[Ana_component]  No such component: " << component << "!" << endl;
        return 0;
    }
}

// gaussian integrands
double pgaus(double *x, double *par)//x,W,M,CC1,CC2,FF,phi1,phi2,GammaVar,Gamma_eeVar
{
    double xx = x[0];
    double W = par[0];
    double M = par[1];
    double CC1 = par[2];
    double CC2 = par[3];
    double FF = par[4];
    double phi1 = par[5];
    double phi2 = par[6];
    double GammaVar = par[7];
    double Gamma_eeVar = par[8];
    double SE = par[9];
    return 1./(sqrt(2*PI)*SE) * exp(-pow((W - xx),2) / (2.*SE*SE)) * Ana(xx, M,CC1,CC2,FF,phi1,phi2, GammaVar,Gamma_eeVar);
}

double pgaus_con(double *x, double *par)//x,W,M,CC1,CC2,FF,phi1,phi2
{
    double xx = x[0];
    double W = par[0];
    double M = par[1];
    double CC1 = par[2];
    double CC2 = par[3];
    double FF = par[4];
    double phi1 = par[5];
    double phi2 = par[6];
    double GammaVar = par[7];
    double Gamma_eeVar = par[8];
    double SE = par[9];
    return 1./(sqrt(2*PI)*SE) * exp(-pow((W - xx),2) / (2.*SE*SE)) * Ana_component(xx, M,CC1,CC2,FF,phi1,phi2, GammaVar,Gamma_eeVar, "con");
}

double pgaus_res(double *x, double *par)//x,W,M,CC1,CC2,FF,phi1,phi2
{
    double xx = x[0];
    double W = par[0];
    double M = par[1];
    double CC1 = par[2];
    double CC2 = par[3];
    double FF = par[4];
    double phi1 = par[5];
    double phi2 = par[6];
    double GammaVar = par[7];
    double Gamma_eeVar = par[8];
    double SE = par[9];
    return 1./(sqrt(2*PI)*SE) * exp(-pow((W - xx),2) / (2.*SE*SE)) * Ana_component(xx, M,CC1,CC2,FF,phi1,phi2, GammaVar,Gamma_eeVar, "res");
}

double pgaus_int(double *x, double *par)//x,W,M,CC1,CC2,FF,phi1,phi2
{
    double xx = x[0];
    double W = par[0];
    double M = par[1];
    double CC1 = par[2];
    double CC2 = par[3];
    double FF = par[4];
    double phi1 = par[5];
    double phi2 = par[6];
    double GammaVar = par[7];
    double Gamma_eeVar = par[8];
    double SE = par[9];
    return 1./(sqrt(2*PI)*SE) * exp(-pow((W - xx),2) / (2.*SE*SE)) * Ana_component(xx, M,CC1,CC2,FF,phi1,phi2, GammaVar,Gamma_eeVar, "int");
}


double Conv(double x,double *par)
{
	TF1* func=new TF1("func",pgaus,-1,1,10);
    func->SetParameters(x,par[0],par[1],par[2],par[3],par[4],par[5],par[6],par[7],par[8]);        // Here x is W
	//return (func->Integral(x-5*0.0013,x+5*0.0013));
    double a=func->Integral(x-10*par[8],x+10*par[8],TOLINTEGRATE);
    delete func;
    return a;
}




double Conv_component(double x,double *par, TString component = "total")
{
	// TF1* func=new TF1("func",pgaus_component,-1,1,8);
	TF1* func = nullptr;
    if(component == "total"){
        func=new TF1("func",pgaus,-1,1,10);
    }
    else if(component == "con"){
        func=new TF1("func",pgaus_con,-1,1,10);
    }
    else if(component == "res"){
        func=new TF1("func",pgaus_res,-1,1,10);
    }
    else if(component == "int"){
        func=new TF1("func",pgaus_int,-1,1,10);
    }
    else{
        cout<<"No such component for convolution!"<<endl;
        return 0;
    }
    func->SetParameters(x,par[0],par[1],par[2],par[3],par[4],par[5],par[6],par[7],par[8]);
	//return (func->Integral(x-5*0.0013,x+5*0.0013));
    double a=func->Integral(x-10*par[8],x+10*par[8],1e-10);
    delete func;
    return a;
}

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{   // arg list: M, CC1, CC2, FF, phi1, phi2, Gamma, Gamma_ee, SE, fCorr
    //           0  1    2    3   4     5     6      7         8   9
    double chisq = 0;
    // double data[Arsize];
    const int total = numPara;
    double arr[total]={0};
    double fCorr = par[total-1];
    double m_M = par[0];
    double m_SE = par[total-2];
    double m_Gamma = par[total-4];
    double m_Gamma_ee = par[total-3];
    double errorCorr = yerrsyscor[0]/100.;

    // double Br = BrCalc(par[1], par[3]);
    double Br = BrCalc(par[1], par[3], par[0], par[6], par[7]);

    for (Int_t i=0;i<total-1;i++)
    {
        arr[i]=par[i];
    }

    for (Int_t i= 0; i < Arsize; i++) 
    {   
        double errorTot = sqrt(yerrsta[i] * yerrsta[i] + yerrsysuncor[i] * yerrsysuncor[i]);
        chisq = chisq+pow((ydata[i]-fCorr*Conv(par[i+total],arr))/errorTot,2); // Considering uncorr. and sta. syst.
    }

    for (Int_t i= 0; i < Arsize; i++) 
    {   
       
            // chisq = chisq+pow((par[i+total]-xdata[i])/(3*dEnergy[i]),2);
            chisq = chisq+pow((par[i+total]-xdata[i])/(dEnergy[i]),2);
        
    }
    
    if(errorCorr != 0) chisq = chisq+pow((1-fCorr)/errorCorr,2); // correlated systematics

    // SE pull term
    chisq = chisq + pow((m_SE - SE_mean)/SE_std,2);

    // Gamma and Gamma_ee pull terms
    chisq = chisq + pow((m_Gamma - Gamma_mean)/Gamma_std,2);
    chisq = chisq + pow((m_Gamma_ee - Gamma_ee_mean)/Gamma_ee_std,2);    

    // Mass pull term
    chisq = chisq + pow((m_M - MJpsi_mean)/MJpsi_std,2);

    // Br pull term
    // chisq = chisq + pow((Br - Br_mean)/Br_std,2);

    f = chisq;
    // res=chisq+pow((1-par[total-1])/0.05,2);

}

void fcn_C2phi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{   // arg list: M, CC1, CC2, FF, phi1, phi2, Gamma, Gamma_ee, SE, fCorr
    //           0  1    2    3   4     5     6      7         8   9
    double chisq = 0;
    // double data[Arsize];
    const int total = numPara;
    double arr[total]={0};
    double fCorr = par[total-1];
    double m_M = par[0];
    double m_SE = par[total-2];
    double m_Gamma = par[total-4];
    double m_Gamma_ee = par[total-3];
    double errorCorr = yerrsyscor[0]/100.;

    double C1 = CC1fromCC2(par[2], par[5]);
    double Br = BrCalc(C1, par[3], par[0], par[6], par[7]);

    for (Int_t i=0;i<total-1;i++)
    {
        arr[i]=par[i];
    }

    for (Int_t i= 0; i < Arsize; i++) 
    {   
        double errorTot = sqrt(yerrsta[i] * yerrsta[i] + yerrsysuncor[i] * yerrsysuncor[i]);
        chisq = chisq+pow((ydata[i]-fCorr*Conv(par[i+total],arr))/errorTot,2); // Considering uncorr. and sta. syst.
    }

    for (Int_t i= 0; i < Arsize; i++) 
    {   
       
            // chisq = chisq+pow((par[i+total]-xdata[i])/(3*dEnergy[i]),2);
            chisq = chisq+pow((par[i+total]-xdata[i])/(dEnergy[i]),2);
        
    }
    
    if(errorCorr != 0) chisq = chisq+pow((1-fCorr)/errorCorr,2); // correlated systematics

    // SE pull term
    chisq = chisq + pow((m_SE - SE_mean)/SE_std,2);

    // Gamma and Gamma_ee pull terms
    chisq = chisq + pow((m_Gamma - Gamma_mean)/Gamma_std,2);
    chisq = chisq + pow((m_Gamma_ee - Gamma_ee_mean)/Gamma_ee_std,2);    

    // Mass pull term
    chisq = chisq + pow((m_M - MJpsi_mean)/MJpsi_std,2);

    // Br pull term
    // chisq = chisq + pow((Br - Br_mean)/Br_std,2);

    f = chisq;
    // res=chisq+pow((1-par[total-1])/0.05,2);

}

void fcn_linearW(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{   // arg list: M, CC1, CC2, FF, phi1, phi2, Gamma, Gamma_ee, SE, fCorr
    //           0  1    2    3   4     5     6      7         8   9
    double chisq = 0;
    // double data[Arsize];
    const int total = numPara;
    double arr[total]={0};
    double fCorr = par[total-1];
    double m_M = par[0];
    double m_SE = par[total-2];
    double m_Gamma = par[total-4];
    double m_Gamma_ee = par[total-3];
    double errorCorr = yerrsyscor[0]/100.;

    // double Br = BrCalc(par[1], par[3]);
    double Br = BrCalc(par[1], par[3], par[0], par[6], par[7]);

    for (Int_t i=0;i<total-1;i++)
    {
        arr[i]=par[i];
    }

    for (Int_t i= 0; i < Arsize; i++) 
    {   
        // double errorTot = sqrt(yerrsta[i] * yerrsta[i] + yerrsysuncor[i] * yerrsysuncor[i]);
        // linear approximation of derivative
        double dW = dEnergy[i];
        double yerr_W = fabs(Conv(par[i+total]+dW, arr) - Conv(par[i+total]-dW, arr)) / 2.0; // numerical derivative
        double errorTot = sqrt(yerr_W * yerr_W + yerrsta[i] * yerrsta[i] + yerrsysuncor[i] * yerrsysuncor[i]);

        chisq = chisq+pow((ydata[i]-fCorr*Conv(par[i+total],arr))/errorTot,2); // Considering uncorr. and sta. syst.
    }

    // for (Int_t i= 0; i < Arsize; i++) 
    // {   
       
    //         // chisq = chisq+pow((par[i+total]-xdata[i])/(3*dEnergy[i]),2);
    //         chisq = chisq+pow((par[i+total]-xdata[i])/(dEnergy[i]),2);
        
    // }
    
    if(errorCorr != 0) chisq = chisq+pow((1-fCorr)/errorCorr,2); // correlated systematics

    // SE pull term
    chisq = chisq + pow((m_SE - SE_mean)/SE_std,2);

    // Gamma and Gamma_ee pull terms
    chisq = chisq + pow((m_Gamma - Gamma_mean)/Gamma_std,2);
    chisq = chisq + pow((m_Gamma_ee - Gamma_ee_mean)/Gamma_ee_std,2);    

    // Mass pull term
    chisq = chisq + pow((m_M - MJpsi_mean)/MJpsi_std,2);

    // Br pull term
    // chisq = chisq + pow((Br - Br_mean)/Br_std,2);

    f = chisq;
    // res=chisq+pow((1-par[total-1])/0.05,2);

}



// Br pieces
double CC1Calc(double Br, double FF) {
    // use sqrt instead of pow(...,0.5)
    // double num = MJpsi_mean*MJpsi_mean - mo*mo - mp*mp;
    // double discriminant = num*num - 4.0 * mp*mp * mo*mo;
    // double qf = sqrt(discriminant) / (2.0 * MJpsi_mean);
    double qf = sqrt(pow(1 + (mo*mo - mp*mp)/pow(MJpsi_mean, 2), 2) - (4 * mo*mo)/pow(MJpsi_mean, 2));

    // square and cube explicitly, or via std::pow
    // double CC1_sq = CC1 * CC1;
    double FF_sq  = FF  * FF;
    double qf_cu  = qf   * qf * qf;

    // build the branching ratio
    // double result = CC1_sq * FF_sq * GAM_ee * qf_cu
    //               / (pow(MV, 4) * GAM_tot);
    // double result = sqrt(pow(MJpsi_mean, 4) * Gamma_mean * Br / (FF_sq * qf_cu * Gamma_ee_mean));
    double result = sqrt(8 * pow(MJpsi_mean, 6) * Gamma_mean * Br / (FF_sq * qf_cu * Gamma_ee_mean));       // W^-8

    return result;
}

double CC1Calc(double Br, double FF, double M, double Gamma, double Gamma_ee) {
    // use sqrt instead of pow(...,0.5)
    // double num = M*M - mo*mo - mp*mp;
    // double discriminant = num*num - 4.0 * mp*mp * mo*mo;
    // double qf = sqrt(discriminant) / (2.0 * M);
    double qf = sqrt(pow(1 + (mo*mo - mp*mp)/pow(M, 2), 2) - (4 * mo*mo)/pow(M, 2));

    // square and cube explicitly, or via std::pow
    // double CC1_sq = CC1 * CC1;
    double FF_sq  = FF  * FF;
    double qf_cu  = qf   * qf * qf;

    // build the branching ratio
    // double result = CC1_sq * FF_sq * GAM_ee * qf_cu
    //               / (pow(MV, 4) * GAM_tot);
    // double result = sqrt(pow(M, 4) * Gamma * Br / (FF_sq * qf_cu * Gamma_ee));
    double result = sqrt(8 * pow(M, 6) * Gamma * Br / (FF_sq * qf_cu * Gamma_ee));                      // W^-8

    return result;
}

double ConvBr(double x,double *par)
{
    // double CC1 = CC1Calc(par[1], par[3]);
    double CC1 = CC1Calc(par[1], par[3], par[0], par[6], par[7]);
	TF1* func=new TF1("func",pgaus,-1,1,10);
    // func->SetParameters(x,par[0],CC1,par[2],par[3],par[4],par[5],par[6]);
    func->SetParameters(x,par[0],CC1,par[2],par[3],par[4],par[5],par[6],par[7],par[8]);        // Here x is W
	//return (func->Integral(x-5*0.0013,x+5*0.0013));
    double a=func->Integral(x-10*par[8],x+10*par[8],1e-10);
    delete func;
    return a;
}

void fcn_Br(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    double chisq = 0;
    // double data[Arsize];
    const int total = numPara;
    double arr[total]={0};
    double fCorr = par[total-1];
    double m_M = par[0];
    double m_SE = par[total-2];
    double m_Gamma = par[total-4];
    double m_Gamma_ee = par[total-3];
    double errorCorr = yerrsyscor[0]/100.;

    double Br = par[1]; // Br is now a fit parameter

    for (Int_t i=0;i<total-1;i++)
    {
        arr[i]=par[i];
    }

    for (Int_t i= 0; i < Arsize; i++) 
    {   
        double errorTot = sqrt(yerrsta[i] * yerrsta[i] + yerrsysuncor[i] * yerrsysuncor[i]);
        chisq = chisq+pow((ydata[i]-fCorr*ConvBr(par[i+total],arr))/errorTot,2); // Considering uncorr. and corr. syst.
    }

    for (Int_t i= 0; i < Arsize; i++) 
    {   
       
            // chisq = chisq+pow((par[i+total]-xdata[i])/(3*dEnergy[i]),2);
            chisq = chisq+pow((par[i+total]-xdata[i])/(dEnergy[i]),2);
        
    }
    
    if(errorCorr != 0) chisq = chisq+pow((1-fCorr)/errorCorr,2); // correlated systematics

    // SE pull term
    chisq = chisq + pow((m_SE - SE_mean)/SE_std,2);

    // Gamma and Gamma_ee pull terms
    chisq = chisq + pow((m_Gamma - Gamma_mean)/Gamma_std,2);
    chisq = chisq + pow((m_Gamma_ee - Gamma_ee_mean)/Gamma_ee_std,2);    

    // Mass pull term
    chisq = chisq + pow((m_M - MJpsi_mean)/MJpsi_std,2);

    // Br pull term
    // chisq = chisq + pow((Br - Br_mean)/Br_std,2);

    f = chisq;
    // res=chisq+pow((1-par[total-1])/0.05,2);

}





#endif // PHYSICSFUNCS_H