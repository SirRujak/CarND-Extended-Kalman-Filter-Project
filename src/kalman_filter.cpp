#include "kalman_filter.h"

#include <iostream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &Hj_in, MatrixXd &R_laser_in,
                        MatrixXd &R_radar_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  Hj_ = Hj_in;
  R_laser_ = R_laser_in;
  R_radar_ = R_radar_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_laser_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  float rho = sqrt(x_[0]*x_[0] + x_[1]*x_[1]);
  float phi = atan2(x_[1],x_[0]);
  float rho_v = 0;
  if (rho >= 0.0001) {
    rho_v = (x_[0]*x_[2] + x_[1]*x_[3])/rho;
  }
  VectorXd hj = VectorXd(3);
  hj << rho, phi, rho_v;
  VectorXd y = z - hj;
  y[1] = atan2(sin(y[1]),cos(y[1]));
  MatrixXd Ht = Hj_.transpose();
	MatrixXd S = Hj_ * P_ * Ht + R_radar_;
	MatrixXd Si = S.inverse();
	MatrixXd K = P_ * Ht * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * Hj_) * P_;
}
