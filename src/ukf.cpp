#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.25; // consider bicycle acceleration of 0.5 m/s^2

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2*M_PI;

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

  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);

  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for(unsigned int i = 1; i < 2*n_aug_+1; i++) {
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  // initial time
  time_us_ = 0;

  // NIS
  NIS_radar_ = 0;
  NIS_laser_ = 0;
}

UKF::~UKF() {}

/**
 * AugmentedSigmaPoints Generate Sigma points for the augmented state
 */
MatrixXd UKF::AugmentedSigmaPoints() {

  // create augmented mean vector
  VectorXd x_aug = VectorXd(7);
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // create augmented covariance matrix
  MatrixXd P_aug = MatrixXd(7, 7);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);

  // create square root matrix of P_aug
  MatrixXd L = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for(unsigned int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i+1) = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

  return Xsig_aug;
}

/**
 * PredictSigmaPoints Predicts sigma points
 * @param Xsig_aug Matrix with the augmented Sigma points
 */
void UKF::PredictSigmaPoints(MatrixXd Xsig_aug, double delta_t) {

  // for each column of the augmented Sigma points matrix
  for(unsigned int i = 0; i < 2*n_aug_+1; i++) {
    // extract values for better readability
    double px = Xsig_aug(0,i);
    double py = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if(fabs(yawd) > 0.001) {
      px_p = px + v/yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
      py_p = py + v/yawd * (cos(yaw) - cos(yaw+yawd*delta_t));
    } else {
      px_p = px + v*delta_t*cos(yaw);
      py_p = py + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // add noise
    px_p += 0.5*nu_a*delta_t*delta_t*cos(yaw);
    py_p += 0.5*nu_a*delta_t*delta_t*sin(yaw);
    v_p += nu_a*delta_t;
    yaw_p += 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p += nu_yawdd*delta_t;

    // write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
}

/**
 * PredictMeanAndCovariance Predict state mean and covariance matrix
 */
void UKF::PredictMeanAndCovariance() {

  // predicted state mean
  x_.fill(0.0);
  for(unsigned int i = 0; i < 2*n_aug_+1; i++) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  // predicted state covariance matrix
  P_.fill(0.0);
  for(unsigned int i = 0; i < 2*n_aug_+1; i++) {

    //state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // angle normalization
    while(x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while(x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  if(!is_initialized_) {

    if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {

      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double rho_dot = meas_package.raw_measurements_(2);

      // initializing x_ (CAN BE TUNED!!!)
      x_ << rho * cos(phi), rho * sin(phi), 0, rho_dot * cos(phi), rho_dot * sin(phi);

      // initializing P_ (CAN BE TUNED!!!)
      P_ << std_radr_*std_radr_, 0, 0, 0, 0,
            0, std_radr_*std_radr_, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, std_radphi_*std_radphi_, 0,
            0, 0, 0, 0, std_radphi_*std_radphi_;

    } else if(meas_package.sensor_type_ == MeasurementPackage::LASER) {

      double px = meas_package.raw_measurements_(0);
      double py = meas_package.raw_measurements_(1);

      // initializing x_ (CAN BE TUNED!!!)
      x_ << px, py, 5.2, 0.0, 0.0;

      // initializing P_ (CAN BE TUNED!!!)
      P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
            0, std_laspx_*std_laspx_, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;
    }

    // done initializing
    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;

    return;
  }

  // calculate delta_t
  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  // Predict
  Prediction(delta_t);

  // Measurement updates
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  MatrixXd Xsig_aug_ = AugmentedSigmaPoints();
  PredictSigmaPoints(Xsig_aug_, delta_t);
  PredictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  // predict lidar measurement --------------------------------------
  // set measurement dimension (lidar can measure px, py)
  int n_z = 2;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);

  // transform sigma points into measurement space
  for(unsigned int i = 0; i < 2*n_aug_+1; i++) {

    // extract values for better readibility
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);

    // measurement model
    Zsig(0,i) = px; //px
    Zsig(1,i) = py; //py
  }

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for(unsigned int i = 0; i < 2*n_aug_+1; i++) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for(unsigned int i = 0; i < 2*n_aug_+1; i++) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S += weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;
  S += R;

  // Update step ----------------------------------------------------

  // lidar measurement
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_(0),
       meas_package.raw_measurements_(1);

  // create matrix for cross-correlation Tc
  MatrixXd Tc = MatrixXd(n_x_,n_z);
  Tc.fill(0.0);

  // calculate cross-correlation matrix
  for(unsigned int i = 0; i < 2*n_aug_+1; i++) {

    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while(z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while(z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while(x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while(x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;

  // angle normalization
  while(z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while(z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  // calculate NIS
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  // predict radar measurement --------------------------------------
  // set measurement dimension (radar can measure rho, phi, rho_dot)
  int n_z = 3;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);

  // transform sigma points into measurement space
  for(unsigned int i = 0; i < 2*n_aug_+1; i++) {

    // extract values for better readibility
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(px*px + py*py);                    // rho
    Zsig(1,i) = atan2(py,px);                           // phi
    Zsig(2,i) = (px*v1 + py*v2) / sqrt(px*px + py*py);  // rho_dot
  }

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for(unsigned int i = 0; i < 2*n_aug_+1; i++) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for(unsigned int i = 0; i < 2*n_aug_+1; i++) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    while(z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while(z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S += weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_*std_radr_, 0, 0,
       0, std_radphi_*std_radphi_, 0,
       0, 0, std_radrd_*std_radrd_;
  S += R;

  // Update step ----------------------------------------------------

  // radar measurement
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_(0),
       meas_package.raw_measurements_(1),
       meas_package.raw_measurements_(2);

  // create matrix for cross-correlation Tc
  MatrixXd Tc = MatrixXd(n_x_,n_z);
  Tc.fill(0.0);

  // calculate cross-correlation matrix
  for(unsigned int i = 0; i < 2*n_aug_+1; i++) {

    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while(z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while(z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while(x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while(x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;

  // angle normalization
  while(z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while(z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  // calculate NIS
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

}