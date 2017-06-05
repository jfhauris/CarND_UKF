#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#define SMA 0.001 // Just a small number


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;


/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_.fill(0.0);



  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = .50;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:
  Complete the initialization. See ukf.h for other member properties.
  Hint: one or more values initialized above might be wildly off...
  */

  ///* Augmented state dimension
  is_initialized_ = false;

  // launch the missing values
  n_x_ = 5;
  n_aug_ = n_x_ + 2;
  lambda_ = 3 - n_aug_;
  x_aug = VectorXd(7);
  x_aug.fill(0.0);
  Q_ = MatrixXd(2, 2);
  P_ = MatrixXd(5, 5);
  P_aug = MatrixXd(7, 7);
  P_aug.fill(0.0);

  H_laser_ = MatrixXd(2, 5);

  R_laser_ = MatrixXd(2, 2);


  //fill in the values accordingly
  Q_ << std_a_ * std_a_, 0,
          0, std_yawdd_ * std_yawdd_;

  P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

  H_laser_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0;

  //P_aug.topLeftCorner(n_x_, n_x_) = P_;
  //P_aug.bottomRightCorner(Q_.rows(), Q_.cols()) = Q_;

  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill(0.0);
  weights_.segment(1, 2 * n_aug_).fill(0.5 / (n_aug_ + lambda_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  //make the sigma point matrix
  Xsig_ = MatrixXd::Zero(n_aug_, 2 * n_aug_ + 1);
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  n_z_radar_ = 3;
  Zsig_ = MatrixXd(n_z_radar_, 2 * n_aug_ + 1);
  Zsig_.fill(0.0);

  R_laser_ << std_laspx_, 0,
          0, std_laspy_;

  R_radar_ = MatrixXd(n_z_radar_, n_z_radar_);
  R_radar_ << std_radr_ * std_radr_, 0, 0,
          0, std_radphi_ * std_radphi_, 0,
          0, 0, std_radrd_ * std_radrd_;
}

UKF::~UKF() {}



/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if(!is_initialized_){
    double py=0;
    double px=0;


    if(meas_package.sensor_type_==MeasurementPackage::RADAR){
      float rho=meas_package.raw_measurements_[0];
      float phi=meas_package.raw_measurements_[1];
      float rho_dot=meas_package.raw_measurements_[2];

      //convert to cartisian
      float px=rho*cos(phi);
      float py=rho*sin(phi);
      float vx=rho_dot*cos(phi);
      float vy=rho_dot*sin(phi);
      x_<<px,py,vx,vy;
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;

      //avoid divide by zero
      if (fabs(x_(0)) < SMA and fabs(x_(1)) < SMA){
        x_(0) = SMA;
        x_(1) = SMA;
      }
    }

    //initialize weights
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for (int i = 1; i < weights_.size(); i++) {
      weights_(i) = 0.5 / (n_aug_ + lambda_);
    }
    //get timestamp
    time_us_ = meas_package.timestamp_;
    //now initialized
    is_initialized_=true;
    return;
  }                     //END if(!is_finished) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  //calc the difference between timestamps
  double dt= (meas_package.timestamp_-time_us_);
  //convert to seconds
  dt=dt/1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(dt);



  //update
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    //cout << "Radar " << meas_package.raw_measurements_[0] << " " << meas_package.raw_measurements_[1] << endl;
    UpdateRadar(meas_package);
  }
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    //cout << "Lidar " << meas_package.raw_measurements_[0] << " " << meas_package.raw_measurements_[1] << endl;
    UpdateLidar(meas_package);
  }
}                     //END UKF :: ukf <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<,
              

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */


  //create augmented mean state
  //x_aug << x_.array(), 0, 0;

  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(Q_.rows(), Q_.cols()) = Q_;

  //calculate square root of P
  MatrixXd A = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_.colwise() = x_aug;
  MatrixXd offset = A * sqrt(lambda_ + n_aug_);

  Xsig_.block(0, 1, n_aug_, n_aug_) += offset;
  Xsig_.block(0, n_aug_ + 1, n_aug_, n_aug_) -= offset;


  //predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    //extract values for better readability
    double p_x = Xsig_(0, i);
    double p_y = Xsig_(1, i);
    double v = Xsig_(2, i);
    double yaw = Xsig_(3, i);
    double yawd = Xsig_(4, i);
    double nu_a = Xsig_(5, i);
    double nu_yawdd = Xsig_(6, i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    } else {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);

  }


  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }


}




/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:
  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.
  You'll also need to calculate the lidar NIS.
  */


  VectorXd z = meas_package.raw_measurements_;

  VectorXd z_pred = H_laser_ * x_;
  VectorXd z_diff = z - z_pred;
  MatrixXd Ht = H_laser_.transpose();
  MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * z_diff);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_laser_) * P_;


  NIS_laser_ = z_diff.transpose() * Si * z_diff;

}



/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.
  You'll also need to calculate the radar NIS.
  */
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    // extract values
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double v_y = cos(yaw) * v;
    double v_x = sin(yaw) * v;

    // measurement model
    double rho = sqrt(p_x * p_x + p_y * p_y);
    double phi = atan2(p_y, p_x);
    double rho_dot = (p_x * v_y + p_y * v_x) / rho;

    if (rho != rho) {
      rho = 0;
    }
    if (phi != phi) {
      phi = 0;
    }
    if (rho_dot != rho_dot) {
      rho_dot = 0;
    }

    Zsig_(0, i) = rho;
    Zsig_(1, i) = phi;
    Zsig_(2, i) = rho_dot;
  }
  cout<<"NEW CYLCE------------------------------------"<<endl;
  //cout<<"Zsig"<<Zsig_<<endl;

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_radar_);
  z_pred.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred = z_pred + weights_(i) * Zsig_.col(i);

  }
  //cout<<"z_pred"<<z_pred<<endl;

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_radar_, n_z_radar_);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred;

    //angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();

  }

  cout<<"S"<<S<<endl;

  //add measurement noise covariance matrix
  S = S + R_radar_;

  //cout<<"S_2"<<S<<endl;


  //calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z_radar_);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff_2 = Zsig_.col(i) - z_pred;
    //angle normalization
    while (z_diff_2(1) > M_PI) z_diff_2(1) -= 2. * M_PI;
    while (z_diff_2(1) < -M_PI) z_diff_2(1) += 2. * M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff_2.transpose();
  }

  //cout<<"Tc"<<Tc<<endl;

  //Kalman gain K;

  MatrixXd K = Tc * S.inverse();



  cout<<"K"<<K<<endl;

  //residual
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;
  cout<<"z_diff"<<z_diff<<endl;

  //angle normalization
  while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
  cout<<"X"<<x_<<endl;
  cout<<"P"<<P_<<endl;

  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}
