#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include <math.h>

using namespace std;

class MPC {
 public:
  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

static double to_pi_range(double psi) {
  while (true) { // convert to [-pi, pi]
    if (psi > M_PI)
      psi -= 2*M_PI;
    else if (psi < -M_PI)
      psi += 2*M_PI;
    else
      break;
  }
 return psi;
}

static void closest_point(Eigen::VectorXd coeffs, double x, double y,
                          double *dst_x, double *dst_y, double *dst_psi, bool reverse) {
  *dst_x = x;
  *dst_y = 0;
  for (int i=0; i<coeffs.size(); i++) {
    *dst_y += coeffs[i] * pow(*dst_x,i);
  }
  *dst_psi = 0;
  for (int i=1; i<coeffs.size(); i++) {
    *dst_psi += i * coeffs[i] * pow(*dst_psi, i-1);
  }
  *dst_psi = atan(*dst_psi);
  if (reverse)
    *dst_psi += M_PI;
  *dst_psi = to_pi_range(*dst_psi);
}

static void state_transform(double co_x, double co_y, double co_psi,
                       double *x, double *y, double *psi) {
  *x -= co_x;
  *y -= co_y;
  *x = *x * cos(co_psi) + *y * sin(co_psi);
  *y = *y * cos(co_psi) - *x * sin(co_psi);
  if (!psi)
    return;

  *psi -= co_psi;
  *psi = to_pi_range(*psi);
}

#endif /* MPC_H */
