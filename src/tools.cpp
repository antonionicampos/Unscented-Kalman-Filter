#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	
  VectorXd rmse = VectorXd(4);
  rmse << 0, 0, 0, 0;

  if(estimations.size() == 0 || estimations.size() != ground_truth.size()) {
  	cout << "Invalid estimation or ground truth data." << endl;
  	return rmse;
  }

  // Accumulate squared residuals
  for(unsigned int i = 0; i < estimations.size(); i++) {

  	VectorXd residual = estimations[i] - ground_truth[i];

  	// Coefficient-wise multiplication
  	residual = residual.array()*residual.array();
  	rmse += residual;
  }

  // Calculate the mean
  rmse = rmse/estimations.size();

  // Calculate the squared root
  rmse = rmse.array().sqrt();

  return rmse;
}