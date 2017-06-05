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
  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  // ... your code here
  if(estimations.size() != ground_truth.size() || estimations.size() == 0)
  {
    //std::cout << "Invalid estimation or ground_truth data" << std::endl;
    return rmse;
  }

  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i)
  {
    // ... your code here
    VectorXd residual = estimations[i] - ground_truth[i];
    
    //coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  //calculate the mean
  // ... your code here
  rmse = rmse/estimations.size();

  //calculate the squared root
  // ... your code here
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  
  //TODO: YOUR CODE HERE 
  //pre-compute a set of terms to avoid repeated calculation
  float c1 = px*px+py*py;
  float c2 = sqrt(c1);
  float c3 = (c1*c2);	

  //check division by zero
  if(fabs(c1) < 0.0001)
  {
    //std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
    //To prevent errors, I suggest you initialize this matrix to a zero matrix before returning here. In case this is returned right now, your algorithm will throw an error.
    return Hj;
  }
  
  //compute the Jacobian matrix
  Hj << px/sqrt(px*px + py*py), py/sqrt(px*px + py*py), 0, 0,
       -py/(px*px + py*py), px/(px*px + py*py), 0, 0, 
        py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3,  
                         px/sqrt(px*px + py*py), py/sqrt(px*px + py*py);
  return Hj;
}
