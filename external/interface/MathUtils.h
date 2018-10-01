#ifndef DELPHYS_EXTERNAL_MATHUTILS_H_
#define DELPHYS_EXTERNAL_MATHUTILS_H_

#include "TMath.h"

#include <iostream>
#include <cmath>
#include <tuple>

namespace delphys {

Double_t ComputeDeltaR(Double_t eta1, Double_t eta2,
                       Double_t phi1, Double_t phi2);

std::tuple<Double_t, Double_t> ComputeAxes(
    const std::vector<Double_t> & x_values,
    const std::vector<Double_t> & y_values,
    const std::vector<Double_t> & weights);

} // NOTE end of delphys namespace

#endif // DELPHYS_EXTERNAL_MATHUTILS_H_
