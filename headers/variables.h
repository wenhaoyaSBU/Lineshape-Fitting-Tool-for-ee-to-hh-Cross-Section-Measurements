#ifndef VARIABLES_H
#define VARIABLES_H

#include <iostream>
#include <fstream>
#include "TComplex.h"
#include <complex>


#define TOLS 1e-13
#define TOLH 1e-13
#define TOLINTEGRATE 1e-7
#define PI acos(-1)


double MV = 3.097;
double hbar_c_2 = 0.3893793721e9;
double Gamma = 92.6e-6; // width of mother particle
double ee = 0.05971*92.6e-6; // width of mother particle to e+e-
// double mumu = 5.13e-6;
//double delta = 0.0013;
//double precise=delta/100;
// double Wm = 0.782+0.135;
double Wm = 2.8;
// double Wm = 1.01700000;
double me = 0.511e-03;
double mo = 0.78266;
double mp = 0.1349768;
double alpha = 1.0/137.0;
double B1=0.98823;
double B2=0.893;
//double B3=0.3936;

double res;
int dof=8;
int sum=0;

// SE measurement result
double SE_mean = 0.000900;
double SE_std  = 0.000030;

// Gamma measurement result
double Gamma_mean = 92.6e-6;
double Gamma_std  = 1.7e-6;
double Br_ee_mean = 5.971e-2;
double Br_ee_std  = 0.032e-2;
double Gamma_ee_mean = Gamma_mean * Br_ee_mean;
double Gamma_ee_std  = sqrt( pow(Br_ee_mean*Gamma_std,2) + pow(Gamma_mean*Br_ee_std,2) );

// MJpsi measurement result
double MJpsi_mean = 3.096900;
double MJpsi_std  = 0.000006;


// qf correction parameters
double p0 = 2.13171;
double p1 = 2.68721;
double p2 = 1.05549;
double p3 = 0.123254;
// double p4 = 0.00052782;
double q0Corr = 0.176892;
double q1Corr = -0.3246;
double q2Corr = 0.0201564;
double q3Corr = 0.123254;

double q0Corr_2body = -1.95742;
double q1Corr_2body = 3.04027;
double q2Corr_2body = -1.46013;
double q3Corr_2body = 0.383267;
double q4Corr_2body = -0.0532937;
double q5Corr_2body = 0.00308655;


// Br measurement reslt from pipi Jpsi
double Br_mean = 5.13e-4;
double Br_std_sta  = 0.02e-4;
double Br_std_sys = 0.14e-4;
double Br_std = sqrt(pow(Br_std_sta,2) + pow(Br_std_sys,2));



const int Arsize=24; // 能量点个数
// const int numPara = 8;      // FOR GOD SAKE STOP HARDCODING!!!!!
const int numPara = 10;
const int Param=Arsize + numPara; // energy pull terms + fitting params
// const int vaccsize = 53000;

double xdata[Arsize] = {0};
double Nob[Arsize] =            {0};
double epsilon[Arsize] =        {0};
double yerrsta[Arsize]  =          {0};
double yerrsysuncor[Arsize] =     {0};
double yerrsyscor[Arsize] =     {0};
double yerrsyscor_quad[Arsize] = {0};
double L[Arsize] =          {0};
double yerr1[Arsize] = {0};
double yerr2[Arsize] = {0};
double Nerr[Arsize]={0};
double yerr[Arsize] = {0};
double yerrsys[Arsize]={0};
double xerr[Arsize]={0};
double ydata[Arsize]={0};
double delta[Arsize]={0};
double ddelta[Arsize]={0};
double dEnergy[Arsize]={0};
double Lerr[Arsize]={0};
double desEnergy[Arsize]={0};
double vstart[Param]={0};//0-3 fitting params  4-14:M  15-25:SE
double step[Param]={0};
const int vaccsize = 53000;
double Evac[vaccsize]={0};
double CSvac[vaccsize]={0};




#endif // VARIABLES_H