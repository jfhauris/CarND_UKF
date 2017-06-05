#include "kalman_filter.h"
#include <iostream>
#include "tools.h"
Tools tools;
//Tools tools = Tools::Tools();  DOES NOT WORK - WHY ???

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
  return;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  y[1] = atan2(sin(y[1]),cos(y[1]));  // Normalize phi=y[1] to be between -pi and +pi
  //Normalization is not necessary here because, we do not deal with angles in this update function while dealing with LASER
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
  return;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  // convert current state into polar coordinates
  float rho = sqrt(x_[0]*x_[0] + x_[1]*x_[1]);
  
  float phi = 0.0;   // this will be default case if pw=0;
  if ( (fabs(x_[0])>0.0001 ))
  {
    phi = atan2(x_[1], x_[0]);
  }
  
  float rhodot = 0.0;  // make rhodot 0.0 if too small
  if ( fabs(rho)>0.0001 )
  {
    rhodot = (x_[0]*x_[2] + x_[1]*x_[3]) / rho;
  }
  
  VectorXd hx(3);
  hx << rho, phi, rhodot;
  
  VectorXd z_pred = hx;
  VectorXd y = z - z_pred;
  y[1] = atan2(sin(y[1]),cos(y[1]));  // Normalize phi=y[1] to be between -pi and +pi  
  std::cout << "y[1] = " << y[1] << std::endl;
  
  //std::cout << "!z_ = " << z(1) << "!!!" << std::endl;
  //std::cout << "!!hx_ = " <<  hx(1) << "!!!" << std::endl;
  //std::cout << "!!!y_ = " << y(1) << "!!!" << std::endl;
  
  // compute Jacobian here = Hj = H below
  MatrixXd Hj = tools.CalculateJacobian(x_);
  //Hj_ = tools.CalculateJacobian(x_);
  H_ = Hj;
  
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
  return;
}
