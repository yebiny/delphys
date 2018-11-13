#ifndef DELPHYS_EXTERNAL_MATHUTILS_H_
#define DELPHYS_EXTERNAL_MATHUTILS_H_

#include "TMath.h"

#include <iostream>
#include <cmath>
#include <tuple>

namespace delphys {

template<typename T>
T ComputeDeltaR(T eta1, T eta2, T phi1, T phi2);

template<typename Element>
std::tuple<Element, Element, Element> ComputeAxes(
    const std::vector<Element> & x_values,
    const std::vector<Element> & y_values,
    const std::vector<Element> & weights);

} // NOTE end of delphys namespace

#include "delphys/external/src/MathUtils.tpp"
#endif // DELPHYS_EXTERNAL_MATHUTILS_H_
