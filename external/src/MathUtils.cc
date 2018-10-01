#include "delphys/external/interface/MathUtils.h"

#include "TMatrixDfwd.h"
#include "TVectorD.h"
#include "TVector2.h"

namespace delphys {


Double_t ComputeDeltaR(Double_t eta1, Double_t eta2,
                       Double_t phi1, Double_t phi2) {
  Double_t deta = eta1 - eta2;
  Double_t dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
  Double_t dr = std::hypot(deta, dphi);
  return dr;
}


std::tuple<Double_t, Double_t> ComputeAxes(
    const std::vector<Double_t> & x_values,
    const std::vector<Double_t> & y_values,
    const std::vector<Double_t> & weights) {

  if (x_values.size() != y_values.size() or x_values.size() != weights.size()) {
    std::cerr << "different size" << std::endl;
    return std::make_tuple(-1.0, -1.0);
  }
  Int_t num_points = x_values.size();

  // sum of squared weight 
  Double_t w2_sum = 0.0;
  // components of a covariance matrix
  Double_t m00 = 0.0, m01 = 0.0, m11 = 0.0;

  for (Int_t idx = 0; idx < num_points; idx++) {
    Double_t x = x_values.at(idx);
    Double_t y = y_values.at(idx);
    Double_t w2 = std::pow(weights.at(idx), 2);

    w2_sum += w2;

    m00 += w2 * std::pow(x, 2);
    m11 += w2 * std::abs(x * y);
    m00 += w2 * std::pow(y, 2);
  } 

  TMatrixDSym covariance_matrix(2);
  covariance_matrix(0, 0) = m00;
  covariance_matrix(0, 1) = m01;
  covariance_matrix(1, 1) = m11;

  TVectorD eigen_values;
  covariance_matrix.EigenVectors(eigen_values);

  // length of
  Double_t major_axis = std::sqrt(eigen_values[0] / w2_sum);
  Double_t minor_axis = std::sqrt(eigen_values[1] / w2_sum);

  return std::make_tuple(major_axis, minor_axis);
}




} // end delphys namespace
