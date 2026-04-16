#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <random>  // for Monte Carlo in the tail




double xi2 = 1.64493407;
double xi3 = 1.2020569;
double hbar = 1.05e-34;
double clight = 3e8;
double NItoSI = pow(hbar*clight,2);
double m2Topb = 1e40;
double GeVtoJ = 1.6e-10;





double LI2USER(double x)
{
    // double li2= -x +x*x/4. - x*x*x/9; 
    double li2= +x +x*x/4. + x*x*x/9; // corrected by wangwp
    return li2;
}

double qf_calc(double W){
    return sqrt(pow(W*W-mo*mo-mp*mp,2)-4*mp*mp*mo*mo)/2/W;
}
double qf_fixed = qf_calc(MV);

double qfDL_calc(double W){ // dimensionless definition
    return sqrt(pow(1 + (mo*mo - mp*mp)/(W*W), 2) - (4*mo*mo)/(W*W));
}

// Dressed XS
double FFuser(double W, double A){
    // return A / pow(W, 2);
    return A / pow(W, 4);       // W^-8
}
double sigma_born_dressed(double W, double psi, double A, double C, double vaccV, double MJ, bool useqf=false){
    double qf = qf_fixed;
    if (useqf){
        // qf = qf_calc(W);
        qf = qfDL_calc(W);
    }
    double norm = NItoSI * m2Topb / pow(GeVtoJ, 2);
    // double cont = 4 * PI * pow(alpha, 2) * pow(FFuser(W, A), 2) / (3 * pow(W, 2));
    double cont = PI * pow(alpha, 2) / 6 * pow(FFuser(W, A), 2);       //W^-8
    // double cont = PI * pow(alpha, 2) / 6 * pow(FFuser(W, A), 2) * W;       //W^-7
    // double cont = PI * pow(alpha, 2) / 6 * pow(FFuser(W, A), 2) * W*W;       //W^-6
    complex<double> amp = vaccV + C * 3 * W * W * ee * exp(complex<double>(0, psi)) /
                          (alpha * MJ * (pow(W, 2) - pow(MJ, 2) + complex<double>(0, MJ * Gamma)));
    double ampMod2 = pow(abs(amp), 2);
    return norm * cont * ampMod2 * pow(qf,3);
}

double sigma_born_dressed(double W, double psi, double A, double C, double vaccV, double MJ, TString mode, bool useqf=false){
    // mode 0: total, 1: continuum only, 2: resonance only, 3: interference only
    double qf = qf_fixed;
    if (useqf){
        // qf = qf_calc(W);
        qf = qfDL_calc(W);
    }

    double norm = NItoSI * m2Topb / pow(GeVtoJ, 2);
    // double cont = 4 * PI * pow(alpha, 2) * pow(FFuser(W, A), 2) / (3 * pow(W, 2));
    double cont = PI * pow(alpha, 2) / 6 * pow(FFuser(W, A), 2);       //W^-8
    // double cont = PI * pow(alpha, 2) / 6 * pow(FFuser(W, A), 2) * W;       //W^-7
    // double cont = PI * pow(alpha, 2) / 6 * pow(FFuser(W, A), 2) * W*W;       //W^-6
    complex<double> amp_total = vaccV + C * 3 * W * W * ee * exp(complex<double>(0, psi)) /
                          (alpha * MJ * (pow(W, 2) - pow(MJ, 2) + complex<double>(0, MJ * Gamma)));
    complex<double> amp_cont = vaccV;
    complex<double> amp_res = amp_total - amp_cont;

    double ampMod2_total = pow(abs(amp_total), 2);
    double ampMod2_cont  = pow(abs(amp_cont), 2);
    double ampMod2_res   = pow(abs(amp_res), 2);
    double ampMod2_int = 0.;
    ampMod2_int   = ampMod2_total - ampMod2_cont - ampMod2_res;
    // possibly more robust way to get interference term:
    ampMod2_int = 2.0 * real(amp_cont * conj(amp_res));
    if (mode == "total"){
        return norm * cont * ampMod2_total * pow(qf,3);
    } else if (mode == "con"){
        return norm * cont * ampMod2_cont * pow(qf,3);
    } else if (mode == "res"){
        return norm * cont * ampMod2_res * pow(qf,3);
    } else if (mode == "int"){
        return norm * cont * ampMod2_int * pow(qf,3);
    } else if (mode == "int_alt"){
        return norm * cont * ampMod2_int * pow(qf,3);
    } else {
        // default to total
        return norm * cont * ampMod2_total * pow(qf,3);
    }
}


// ---------- Radiator function etc. definition --------------------------------

double S(double W){
    return pow(W, 2);
}
double X(double s, double m){
    return 1.0 - pow(m, 2)/s;
}
double LL(double s){
    return 2.0 * log( sqrt(s)/me );
}
double BETA(double l){
    return 2.0 * alpha / PI * ( l - 1.0 );
}
double DELTA2(double l){
    // Be careful: 9/8 etc. must be double, not integer division
    return (9.0/8.0 - 2.0*xi2) * pow(l, 2)
         - (45.0/16.0 - 11.0/2.0*xi2 - 3.0*xi3) * l
         - 6.0/5.0*pow(xi2, 2)
         - 9.0/2.0*xi3
         - 6.0*xi2*log(2.0)
         + 3.0/8.0*xi2
         + 57.0/12.0;
}
double DELTA(double delta2, double l){
    return 1.0
         + (alpha/PI) * (1.5*l + (1.0/3.0)*pow(PI, 2) - 2.0)
         + pow((alpha/PI), 2) * delta2;
}

double RADIATOR(double s, double x){
    // Regularise x to avoid log(0) and division by zero
    // const double tiny = 1e-10;
    // if (x < tiny) x = tiny;
    // if (x > 1.0 - tiny) x = 1.0 - tiny;

    double l = LL(s);
    double beta = BETA(l);
    double delta2 = DELTA2(l);
    double delta = DELTA(delta2, l);

    double term1 = delta * beta * pow(x, (beta - 1.0));
    double term2 = -(beta/2.0) * (2.0 - x);
    double term3 = (pow(beta, 2)/8.0) *
                   ((2.0 - x) * (3.0*log(1.0-x) - 4.0*log(x))
                    - 4.0*log(1.0-x)/x
                    - 6.0 + x);
    return term1 + term2 + term3;
}

double narrowCorr(double s){
    double mv = 3.0969;
    double xv = 1.0 - mv*mv/s;
    double contFactor = 12*PI*PI * ee / (mv * s);
    double rad = RADIATOR(s, xv);
    return contFactor * rad;
}

double rho_mJpsi(double m, double W_beam, double psi, double A, double C, double MJ, bool useqf=false){
    double s = S(W_beam);
    double x = X(s, m);
    double vaccValue = vacc(m);
    double radValue = RADIATOR(s, x);
    double dressedValue = sigma_born_dressed(m, psi, A, C, vaccValue, MJ, useqf);
    return (2.0 * m / s) * radValue * dressedValue;
}

double integrand_rho_2(double m, void *p){
    auto *params = static_cast<double *>(p);
    double W_beam = params[0];
    double psi = params[1];
    double A = params[2];
    double C = params[3];
    double MJ = params[4];
    return rho_mJpsi(m, W_beam, psi, A, C, MJ, true);
}
double integrand_rho_2_noqf(double m, void *p){
    auto *params = static_cast<double *>(p);
    double W_beam = params[0];
    double psi = params[1];
    double A = params[2];
    double C = params[3];
    double MJ = params[4];
    return rho_mJpsi(m, W_beam, psi, A, C, MJ, false);
}
double anaIntegral_tail(double b, double W_beam, double psi, double A, double C, double MJ, bool useqf=false){
    double s = S(W_beam);
    double vaccValue = vacc(W_beam);
    double l = LL(s);
    double beta = BETA(l);
    double delta2 = DELTA2(l);
    double delta = DELTA(delta2, l);
    
    double bornValue = sigma_born_dressed(W_beam, psi, A, C, vaccValue, MJ, useqf);
    double ana = delta * pow(b, beta) + 1./32. * (beta * beta * b * b) + 1./4. * beta * b * b - 3./16. * beta * beta * b * b * log(1-b) + 1./4. * beta * beta * b * b * log(b) - 5./16. * beta * beta * b - beta * b + 3./4. * beta * beta * b * log(1-b) - beta * beta * b * log(b) - 9./16. * beta * beta * log(1-b) + 1./2. * beta * beta * LI2USER(b);
    return ana * bornValue;
}

// ---------- main integral with region splitting ------------------------------

double integral_rho_2(double W_beam, double psi, double A, double C, double MJ, bool useqf=false){
    // workspace for GSL
    const size_t limit = 10000;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(limit);

    double params[5] = {W_beam, psi, A, C, MJ};
    gsl_function F;
    if (useqf){
        F.function = &integrand_rho_2;
    } else {
        F.function = &integrand_rho_2_noqf;
    }
    F.params = params;

    // integration limits
    double m_low = Wm;               // physical threshold
    // Avoid exactly hitting m = W_beam (x = 0)
    // double m_up  = W_beam * (1.0 - 1e-8);
    double m_up  = W_beam;

    if (m_up <= m_low) {
        gsl_integration_workspace_free(w);
        return 0.0;
    }

    // Split out a tiny tail region near the upper limit where the integrand is nasty
    // double tail_frac = 0;  // fraction of the interval reserved for ana tail
    // double m_split = m_up - tail_frac * (m_up - m_low);
    // if (m_split <= m_low) m_split = 0.5 * (m_low + m_up);
    double offset = 1e-4; // absolute offset from m_up for the split
    double m_split = m_up - offset;

    // GSL precision settings
    double result_main = 0.0, error_main = 0.0;
    double epsabs = 1e-9;
    double epsrel = 1e-6;

    // Region 1: [m_low, m_split] – handled by qags (no singular endpoint here)
    // int status = gsl_integration_qags(&F, m_low, m_split,
    int status = gsl_integration_qags(&F, m_low, m_split,
                                      epsabs, epsrel, limit,
                                      w, &result_main, &error_main);
    if (status != GSL_SUCCESS) {
        // You can optionally print a warning here
        fprintf(stderr, "qags status %d at W_beam = %g\n", status, W_beam);
    }

    // Region 2: [m_split, m_up] – tiny tail near singular endpoint, use analytical approx.
    double result_tail = 0;
    double b = 1.0 - pow(m_split, 2) / S(W_beam);
    result_tail = anaIntegral_tail(b, W_beam, psi, A, C, MJ, useqf);
    // result_tail = 0;
    

    gsl_integration_workspace_free(w);

    return result_main + result_tail;
}

// Alternative approach: sample 1M points over x = 1 - m^2/s and calculate the integral (Monte Carlo method as in ConExc)
double integral_rho_2_MC(double W_beam, double psi, double A, double C, double MJ, bool useqf=false, bool narrowCorrFlag=false){
    size_t Npoints = 1000000;
    double s = S(W_beam);
    double m_low = Wm;
    double m_up = W_beam;
    double x_low = X(s, m_up); // = 0
    double x_up = X(s, m_low); // = 1 - m_low^2/s

    // Random number generator
    std::mt19937 rng(42); // fixed seed for reproducibility
    std::uniform_real_distribution<double> dist(x_low, x_up);

    double sum = 0.0;
    for (size_t i = 0; i < Npoints; ++i) {
        double x = dist(rng);
        double m = sqrt(s * (1.0 - x));
        double integrand = rho_mJpsi(m, W_beam, psi, A, C, MJ, useqf);
        sum += integrand;
    }

    double integral = (x_up - x_low) * (sum / Npoints);

    double narrowCorrTerm = 0.0;
    if (narrowCorrFlag){
        narrowCorrTerm = narrowCorr(s);
        if(std::isnan(narrowCorrTerm) || std::isinf(narrowCorrTerm)) narrowCorrTerm = 0.0;
    }

    return integral + narrowCorrTerm;
}
void dump_integral_rho_MC(double* fittedParr) {
    const double Wstart = 3.05;
    const double Wend   = 3.12;
    const double dW     = 0.0001;

    const int nPoints = static_cast<int>((Wend - Wstart) / dW + 0.5); // ~700

    double MJ   = fittedParr[0];
    double CC1  = fittedParr[1];
    double FF   = fittedParr[3];
    double phi1 = fittedParr[4];

    std::vector<double> Wvals(nPoints);
    std::vector<double> integrals(nPoints);

    // Precompute W values
    for (int i = 0; i < nPoints; ++i) {
        Wvals[i] = Wstart + i * dW;
    }

    // Parallel loop over Wtest
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < nPoints; ++i) {
        double Wtest = Wvals[i];
        // bool useNarrowCorr = false;
        // if (Wtest > MV - 0.01 && Wtest < MV + 0.01) useNarrowCorr = true;
        double integral = integral_rho_2_MC(Wtest, phi1, FF, CC1, MJ, true, false);
        integrals[i] = integral;
    }

    // Single-threaded file output
    std::ofstream oufFileTest_intgralMC("integral_rhoMC_at_Ws.txt");
    oufFileTest_intgralMC << std::setprecision(10);
    for (int i = 0; i < nPoints; ++i) {
        oufFileTest_intgralMC << Wvals[i] << "\t" << integrals[i] << "\n";
    }
    oufFileTest_intgralMC.close();
}



// -----------------------------------



// Another approach ------------------

double deltad(double beta){
    return 3./4. * beta + alpha/PI * (PI*PI/3. - 1./2.) + beta*beta*(9./32. - PI*PI/12.);
}

double RADIATORF(double x, double W){
    double s = S(W);
    double l = LL(s);
    double beta = BETA(l);
    // double beta = 2*alpha/PI * (2*log(W/me) - 1);
    double deltadd = deltad(beta);

    double term1 = beta * pow(x, beta-1.)*(1.+deltadd);
    double term2 = -1*beta * (1.-x/2.);
    double term3 = 1./8.*beta*beta * ( 4*(2-x) * log(1./x) - (1 + 3*pow(1-x, 2))/x * log(1-x) - 6 + x );
    return term1 + term2 + term3;
}

double sigma_born_dressed_vt(double x, double W, double psi, double A, double C, double MJ, bool useqf = false){
    double m = W * sqrt(1-x);
    double vaccValue = vacc(m);
    return sigma_born_dressed(m, psi, A, C, vaccValue, MJ, useqf);
}

double integrandF(double x, void* p){
    auto *params = static_cast<double *>(p);
    double W_beam = params[0];
    double psi = params[1];
    double A = params[2];
    double C = params[3];
    double MJ = params[4];
    return RADIATORF(x, W_beam) * sigma_born_dressed_vt(x, W_beam, psi, A, C, MJ, true);
}
double integrandF_noqf(double x, void* p){
    auto *params = static_cast<double *>(p);
    double W_beam = params[0];
    double psi = params[1];
    double A = params[2];
    double C = params[3];
    double MJ = params[4];
    return RADIATORF(x, W_beam) * sigma_born_dressed_vt(x, W_beam, psi, A, C, MJ, false);
}

double integralF(double W_beam, double psi, double A, double C, double MJ, bool useqf = false){
    // workspace for GSL
    const size_t limit = 10000;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(limit);

    double params[5] = {W_beam, psi, A, C, MJ};
    gsl_function F;
    if (useqf){
        F.function = &integrandF;
    } else {
        F.function = &integrandF_noqf;
    }
    F.params = params;

    // integration limits
    double xLow = 0;
    double xHigh = 1 - pow(Wm/W_beam, 2);
    if (xHigh <= xLow){
        gsl_integration_workspace_free(w);
        return 0.0;
    }

    // GSL precision settings
    double result_main = 0.0, error_main = 0.0;
    double epsabs = 1e-9;
    double epsrel = 1e-6;

    int status = gsl_integration_qags(&F, xLow, xHigh,
                                      epsabs, epsrel, limit,
                                      w, &result_main, &error_main);
    if (status != GSL_SUCCESS) {
        // You can optionally print a warning here
        fprintf(stderr, "qags status %d at W_beam = %g\n", status, W_beam);
    }

    gsl_integration_workspace_free(w);

    return result_main;

}

// -----------------------------------