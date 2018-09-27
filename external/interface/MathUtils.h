#ifndef DELPHYS_EXTERNAL_MATHUTILS_H_
#define DELPHYS_EXTERNAL_MATHUTILS_H_

#include "TMath.h"

#include <cmath>
#include <iostream>

namespace delphys {

static const Double_t kPi = TMath::Pi();
static const Double_t kTwoPi = TMath::TwoPi();

inline Double_t ComputeDeltaPhi(Double_t phi1, Double_t phi2) {
  Double_t dphi = phi1 - phi2;
  if(TMath::IsNaN(dphi)) {
    std::cerr << "ComputeDeltaPhi function called with NaN" << std::endl;
    return dphi;
  }
  while (dphi >= kPi) dphi -= kTwoPi;
  while (dphi < -kPi) dphi += kTwoPi;
  return dphi;
}

inline Double_t ComputeDeltaR(Double_t deta, Double_t dphi) {
  return std::hypot(deta, dphi);
}

inline Double_t ComputeDeltaR(Double_t eta1, Double_t eta2,
                              Double_t phi1, Double_t phi2) {
  Double_t deta = eta1 - eta2;
  Double_t dphi = ComputeDeltaPhi(phi1, phi2);
  return ComputeDeltaR(deta, dphi);
}

} // NOTE end of delphys namespace

#endif // DELPHYS_EXTERNAL_MATHUTILS_H_
