#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;
  int estimations_size = estimations.size();
  int ground_truth_size = ground_truth.size();
  if (estimations_size != ground_truth_size || estimations_size == 0) {
    cout << "Invalid estimations and ground_truth" << endl;
    return rmse;
  }

  for (unsigned int i = 0; i < estimations_size; ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];

    residual = residual.array() * residual.array();
    rmse += residual;    
  }
  // Calculate mean
  rmse = rmse / estimations_size;
  // SQRT 
  rmse = rmse.array().sqrt();

  // return error
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3, 4); 
  if ( x_state.size() != 4 ) {
    cout << "ERROR - CalculateJacobian () - The state vector must have size 4." << endl;
    return Hj;
  }
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float c1 = px *px + py * py;
  float c2 = sqrt(c1);
  float c3 = (c1 * c2);

  if (fabs(c1) < 0.0001) {
    cout << "CalculateJacobian () - Error Division by Zero" << endl;
    return Hj;
  }
  //compute the Jacobian matrix
	Hj << (px/c2), (py/c2), 0, 0,
		  -(py/c1), (px/c1), 0, 0,
		  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

	return Hj;
}
