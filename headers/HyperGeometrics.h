#ifndef HYPERGEOMETRICS_H
#define HYPERGEOMETRICS_H

#include <iostream>
#include <complex>
#include <cmath>
#include "variables.h"

using namespace std;

inline double pochhammer(double x, int n)
{
    double result = 1.0;
    for (int i = 0; i < n; ++i) result *= (x + i);
    return result;
}

complex<double> hyp2f1_legacy(double a, double b, double c, complex<double> z)
{
    int n=1;
    double eps = 1;
    double m = 1;
    double atemp = 1;
    double btemp = 1;
    double ctemp = 1;
    complex<double> result = {1, 0};
    complex<double> temp = {1, 0};
    for(int i=0;i<10000;i++)
    {
        m = m*n;
        atemp = atemp *a;
        btemp = btemp *b;
        ctemp = ctemp *c;
        temp = atemp*btemp/ctemp*pow(z,n)/m;
        result += temp;
        a = a+1;
        b = b+1;
        c = c+1;
        n = n+1;
        if(fabs(temp)<=TOLH*fabs(result))
        {
            break;
        }

    }
    return result;
}

// Plain power-series 2F1 around z=0 (ONLY reliable for |z| < 1)
inline std::complex<double> hyp2f1(double a, double b, double c,
                                          std::complex<double> z,
                                          bool *converged = nullptr)
{
    using cd = std::complex<double>;

    const int MAXIT = 10000;
    cd sum(1.0, 0.0);   // n=0 term
    cd term(1.0, 0.0);

    bool ok = false;

    for (int n = 0; n < MAXIT; ++n)
    {
        // recurrence:
        // t_{n+1} = t_n * ((a+n)(b+n)/((c+n)(n+1))) * z
        double nn = static_cast<double>(n);
        cd ratio = ((a + nn) * (b + nn) / ((c + nn) * (nn + 1.0))) * z;
        term *= ratio;
        sum += term;

        if (std::abs(term) <= TOLH * std::abs(sum))
        {
            ok = true;
            break;
        }
    }

    if (converged) *converged = ok;
    return sum;
}

// Appell F1 using the single-sum expansion in x
// F1(a;b1,b2;c;x,y) = sum_n (a)_n (b1)_n / ((c)_n n!) x^n * 2F1(a+n,b2;c+n;y)
inline std::complex<double> AppellF1(double a, double b1, double b2, double c,
                                     std::complex<double> x, std::complex<double> y,
                                     bool printi = false)
{
    using cd = std::complex<double>;

    const int MAXIT = 10000;
    cd sum(0.0, 0.0);

    // coefficient for n=0:
    // (a)_0 (b1)_0 / ((c)_0 0!) * x^0 = 1
    cd coeff(1.0, 0.0);

    bool ok = false;
    int used_terms = 0;

    for (int n = 0; n < MAXIT; ++n)
    {
        bool conv2f1 = false;
        cd inner = hyp2f1(a + n, b2, c + n, y, &conv2f1);

        // If y is outside radius, this will not be reliable.
        // You should replace hyp2f1 with an analytic-continuation version.
        cd term = coeff * inner;
        sum += term;

        if (std::abs(term) <= TOLH * std::abs(sum))
        {
            ok = true;
            break;
        }

        // Update coeff -> next n:
        // coeff_{n+1} = coeff_n * ((a+n)(b1+n)/((c+n)(n+1))) * x
        double nn = static_cast<double>(n);
        coeff *= ((a + nn) * (b1 + nn) / ((c + nn) * (nn + 1.0))) * x;
        
        used_terms = n + 1;
    }

    if (printi)
    {
        std::cout << "Calculated F1: " << sum
                  << " with " << used_terms << " terms"
                  << (ok ? "" : " (MAXIT reached)") << std::endl;
    }

    return sum;
}

#endif // HYPERGEOMETRICS_H