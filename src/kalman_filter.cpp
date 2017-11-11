#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

#define eps 0.00001

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
	x_ = F_ * x_;
	P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
	VectorXd y = z - H_ * x_;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K =  P_ * Ht * Si;
	MatrixXd I = MatrixXd::Identity(4, 4);

	//new state
	x_ = x_ + (K * y);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
	float f1 = sqrt(pow(x_[0],2)+pow(x_[1],2));
	if(f1 < eps) {
		f1 = eps;
	}
	VectorXd Hx_ = VectorXd(3);
	Hx_ << f1,
		   atan2(x_[1],x_[0]),
		   (x_[0]*x_[2]+x_[1]*x_[3])/f1;

	VectorXd y = z - Hx_;
	while(y[1]>M_PI) {
		y[1] -= 2*M_PI;
	}
	while (y[1]<-M_PI) {
		y[1] += 2*M_PI;
	}

	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K =  P_ * Ht * Si;
	MatrixXd I = MatrixXd::Identity(4, 4);

	//new state
	x_ = x_ + (K * y);
	P_ = (I - K * H_) * P_;
}
