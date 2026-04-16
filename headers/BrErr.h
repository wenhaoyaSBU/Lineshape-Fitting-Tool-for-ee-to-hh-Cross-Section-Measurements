#ifndef BRERR_H
#define BRERR_H

#include <TMinuit.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>
#include <TString.h>

#include <array>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

// Make sure BrCalc and BrCalcEM are declared before using this header.
// Usually they are already declared in physicsFuncs.h.
// If needed, uncomment the next line.
// #include "physicsFuncs.h"

namespace BrErr
{
    struct Result
    {
        double value = 0.0;
        double error = 0.0;
    };

    using EvalFunc = double (*)(const double*);

    // =========================================================
    // Function wrappers using your EXTERNAL parameter numbering
    //
    // 0 -> M
    // 1 -> CC1
    // 3 -> FF
    // 6 -> Gamma
    // 7 -> Gamma_ee
    // =========================================================

    inline double BrFromPars(const double *par)
    {
        return BrCalc(par[1], par[3], par[0], par[6], par[7]);
    }

    inline double BrEMFromPars(const double *par)
    {
        return BrCalcEM(par[3], par[0], par[6], par[7]);
    }

    // =========================================================
    // Build covariance matrix in EXTERNAL numbering
    // Fixed parameters get zero rows/columns
    // =========================================================
    inline TMatrixDSym GetExternalCovMatrix(TMinuit *minuit, int nExtWanted)
    {
        const int nFree = minuit->GetNumFreePars();

        std::vector<double> covFree(nFree * nFree, 0.0);
        if (nFree > 0) {
            minuit->mnemat(covFree.data(), nFree);
        }

        TMatrixDSym covExt(nExtWanted);
        covExt.Zero();

        // map external index -> free parameter index in mnemat matrix
        std::vector<int> extToFree(nExtWanted, -1);

        for (int iext = 0; iext < nExtWanted; ++iext) {
            TString name;
            double val = 0.0, err = 0.0, low = 0.0, up = 0.0;
            int iuint = -1;

            minuit->mnpout(iext, name, val, err, low, up, iuint);

            // iuint > 0 means variable, and is 1-based
            if (iuint > 0) {
                extToFree[iext] = iuint - 1;
            }
        }

        for (int i = 0; i < nExtWanted; ++i) {
            if (extToFree[i] < 0) continue;
            for (int j = 0; j < nExtWanted; ++j) {
                if (extToFree[j] < 0) continue;
                covExt(i, j) = covFree[extToFree[i] * nFree + extToFree[j]];
            }
        }

        return covExt;
    }

    // =========================================================
    // Generic numerical gradient
    // Only parameters listed in activeIdx are differentiated
    // =========================================================
    inline TVectorD GradientNumeric(EvalFunc func,
                                    const double *par,
                                    int nExtPars,
                                    const std::vector<int> &activeIdx)
    {
        TVectorD grad(nExtPars);
        grad.Zero();

        std::vector<double> pplus(nExtPars), pminus(nExtPars);
        for (int i = 0; i < nExtPars; ++i) {
            pplus[i]  = par[i];
            pminus[i] = par[i];
        }

        for (int k : activeIdx) {
            double h = 1e-6 * std::max(1.0, std::abs(par[k]));
            if (h == 0.0) h = 1e-8;

            pplus[k]  = par[k] + h;
            pminus[k] = par[k] - h;

            const double fplus  = func(pplus.data());
            const double fminus = func(pminus.data());

            grad[k] = (fplus - fminus) / (2.0 * h);

            pplus[k]  = par[k];
            pminus[k] = par[k];
        }

        return grad;
    }

    // =========================================================
    // Specific gradients
    // =========================================================
    inline TVectorD BrGradient(const double *par, int nExtPars)
    {
        // BrCalc depends on M, CC1, FF, Gamma, Gamma_ee
        static const std::vector<int> idx = {0, 1, 3, 6, 7};
        return GradientNumeric(&BrFromPars, par, nExtPars, idx);
    }

    inline TVectorD BrEMGradient(const double *par, int nExtPars)
    {
        // BrCalcEM depends on M, FF, Gamma, Gamma_ee
        static const std::vector<int> idx = {0, 3, 6, 7};
        return GradientNumeric(&BrEMFromPars, par, nExtPars, idx);
    }

    // =========================================================
    // Standard covariance propagation
    // =========================================================
    inline double PropagateError(const TVectorD &grad, const TMatrixDSym &cov)
    {
        double var = 0.0;
        const int n = grad.GetNrows();

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                var += grad[i] * cov(i, j) * grad[j];
            }
        }

        if (var < 0.0 && std::abs(var) < 1e-14) var = 0.0;
        return std::sqrt(std::max(0.0, var));
    }

    // =========================================================
    // One-stop wrappers
    // =========================================================
    inline Result CalcBrAndErr(TMinuit *minuit, const double *fittedParr, int nExtWanted)
    {
        Result out;
        out.value = BrFromPars(fittedParr);

        TMatrixDSym covExt = GetExternalCovMatrix(minuit, nExtWanted);
        TVectorD grad = BrGradient(fittedParr, nExtWanted);
        out.error = PropagateError(grad, covExt);

        return out;
    }

    inline Result CalcBrEMAndErr(TMinuit *minuit, const double *fittedParr, int nExtWanted)
    {
        Result out;
        out.value = BrEMFromPars(fittedParr);

        TMatrixDSym covExt = GetExternalCovMatrix(minuit, nExtWanted);
        TVectorD grad = BrEMGradient(fittedParr, nExtWanted);
        out.error = PropagateError(grad, covExt);

        return out;
    }

    // =========================================================
    // Optional debug printouts
    // =========================================================
    inline void PrintBrInputs(const double *par)
    {
        std::cout << "BrCalc inputs:\n"
                  << "  M        = " << par[0] << "\n"
                  << "  CC1      = " << par[1] << "\n"
                  << "  FF       = " << par[3] << "\n"
                  << "  Gamma    = " << par[6] << "\n"
                  << "  Gamma_ee = " << par[7] << std::endl;
    }

    inline void PrintBrEMInputs(const double *par)
    {
        std::cout << "BrCalcEM inputs:\n"
                  << "  M        = " << par[0] << "\n"
                  << "  FF       = " << par[3] << "\n"
                  << "  Gamma    = " << par[6] << "\n"
                  << "  Gamma_ee = " << par[7] << std::endl;
    }

} // namespace BrErr

#endif