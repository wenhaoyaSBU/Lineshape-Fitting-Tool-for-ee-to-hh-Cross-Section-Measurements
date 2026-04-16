#include <Math/Interpolator.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <stdexcept>

namespace vaccFunc {

  static std::unique_ptr<ROOT::Math::Interpolator> interp;
  static double Wmin = 0.0, Wmax = 0.0;

  inline void InitVacc(TString fname = "./vaccFile/vacc_nojpsi.dat") {
    std::ifstream fin(fname);
    if (!fin.is_open()) throw std::runtime_error("Failed to open vacc file.");

    std::vector<double> E, V;
    E.reserve(60000);
    V.reserve(60000);

    double e, v;
    while (fin >> e >> v) {
      E.push_back(e);
      V.push_back(v);
    }
    fin.close();

    if (E.size() < 2) throw std::runtime_error("Vacc file too small.");

    Wmin = E.front();
    Wmax = E.back();

    interp = std::make_unique<ROOT::Math::Interpolator>(
        (unsigned)E.size(), ROOT::Math::Interpolation::kLINEAR);

    interp->SetData((unsigned)E.size(), E.data(), V.data());

    std::cout << "[vaccFunc] Vacc initialized from file: " << fname << "\n"
              << "  W range: [" << Wmin << ", " << Wmax << "]\n"
              << "  Points: " << E.size() << "\n";
  }

  inline double Vacc_from_table(double W, void* user) {
    auto* tab = static_cast<ROOT::Math::Interpolator*>(user);
    // Optional: clamp to range to avoid exceptions
    if (W < Wmin) W = Wmin;
    if (W > Wmax) W = Wmax;
    return tab->Eval(W); // must be Vacc = 1/|1-Pi0|^2
  }
}