#include "TMatrixTSym.h"
#include "TMatrixDSymfwd.h"
#include "TMatrixFSymfwd.h"
#include "TMatrixDfwd.h"
#include "TMatrixTBase.h" // TVectorT
#include "TVectorF.h"
#include "TVectorD.h"
#include "TVector2.h"

namespace delphys {


template<typename Element>
T ComputeDeltaR(T eta1, T eta2,
                T phi1, T phi2) {
  T deta = eta1 - eta2;
  T dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
  T delta_r = std::hypot(deta, dphi);
  return delta_r;
}


template<typename Element>
std::tuple<Element, Element, Element> ComputeAxes(
    const std::vector<Element> & x_values,
    const std::vector<Element> & y_values,
    const std::vector<Element> & weights) {

  if (x_values.size() != y_values.size() or
      y_values.size() != weights.size() or
      weights.size() != x_values.size()) {

    std::cerr << "ComputeAxes:: different size" << std::endl;
    return std::make_tuple(-1.0, -1.0, -1.0);
  }
  Int_t num_points = x_values.size();

  // sum of weight weight 
  Element w2_sum = 0.0;
  // components of a covariance matrix
  Element m00 = 0.0, m01 = 0.0, m11 = 0.0;

  Element x, y, w2;
  for (Int_t idx = 0; idx < num_points; idx++) {
    x = x_values.at(idx);
    y = y_values.at(idx);
    w2 = std::pow(weights.at(idx), 2);

    w2_sum += w2;

    m00 += w2 * std::pow(x, 2);
    m01 -= w2 * x * y;
    m11 += w2 * std::pow(y, 2);
  } 

  // TMatrixTSm
  // Note that in this implementation both matrix element m[i][j] and m[j][i]
  // are updated and stored in memory . However, when making the object
  // persistent only the upper right triangle is stored .
  TMatrixTSym<Element> covariance_matrix(2);
  covariance_matrix(0, 0) = m00;
  covariance_matrix(0, 1) = m01;
  covariance_matrix(1, 1) = m11;

  TVectorT<Element> eigen_values;
  covariance_matrix.EigenVectors(eigen_values);

  // length of
  Element major_axis = std::sqrt(eigen_values[0] / w2_sum);
  Element minor_axis = std::sqrt(eigen_values[1] / w2_sum);

  Element eccentricity = std::sqrt(1 - (eigen_values[1] / eigen_values[0]));

  return std::make_tuple(major_axis, minor_axis, eccentricity);
}


} // end delphys namespace
