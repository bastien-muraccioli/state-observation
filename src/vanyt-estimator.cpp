#include "state-observation/tools/definitions.hpp"
#include <state-observation/observer/vanyt-estimator.hpp>

namespace stateObservation
{

VanytEstimator::VanytEstimator(double alpha, double beta, double tau, double dt)
: ZeroDelayObserver(13, 9), alpha_(alpha), beta_(beta), dt_(dt)
{
  setTau(tau);
}

void VanytEstimator::initEstimator(const Vector3 & pos, const Vector3 & x1, const Vector3 & x2_prime, const Vector4 & R)
{
  Eigen::VectorXd initStateVector = Eigen::VectorXd::Zero(getStateSize());

  initStateVector.segment<3>(0) = x1;
  initStateVector.segment<3>(3) = x2_prime;
  initStateVector.segment<3>(6) = pos;
  initStateVector.tail(4) = R;

  setState(initStateVector, 0);

  pos_contacts_ = pos;
  T_hat_.position = pos;
  T_hat_.orientation.fromVector4(R);
  pos_x1_.setZero();
}

void VanytEstimator::startNewIteration()
{
  if(k_est_ == k_data_)
  {
    ++k_data_;
    imuAnchorPos_.reset();
    worldAnchorPos_.reset();

    sigma_.setZero();
    oriCorrFromOriMeas_.setZero();
    posCorrFromContactPos_.setZero();
    oriCorrFromContactPos_.setZero();
  }
}

void VanytEstimator::setMeasurement(const Vector3 & yv_k, const Vector3 & ya_k, const Vector3 & yg_k, TimeIndex k)
{
  startNewIteration();

  ObserverBase::MeasureVector y_k(getMeasureSize());
  y_k << yv_k, ya_k, yg_k;

  ZeroDelayObserver::setMeasurement(y_k, k);
}

void VanytEstimator::addOrientationMeasurement(const Matrix3 & oriMeasurement, double gain)
{
  startNewIteration();

  Matrix3 rot_diff = oriMeasurement * T_hat_.orientation.toMatrix3().transpose();
  Vector3 rot_diff_vec = kine::skewSymmetricToRotationVector(rot_diff - rot_diff.transpose());

  oriCorrFromOriMeas_ -= gain * T_hat_.orientation.toMatrix3().transpose() * Vector3::UnitZ()
                         * (Vector3::UnitZ()).transpose() * rot_diff_vec;
}

void VanytEstimator::addContactPosMeasurement(const Vector3 & posMeasurement,
                                              const Vector3 & imuContactPos,
                                              double gainDelta,
                                              double gainSigma)
{
  startNewIteration();

  oriCorrFromContactPos_ +=
      gainSigma
      * (T_hat_.orientation.toMatrix3().transpose() * (posMeasurement - T_hat_.position())).cross(imuContactPos);

  posCorrFromContactPos_ +=
      gainDelta * (imuContactPos - T_hat_.orientation.toMatrix3().transpose() * (posMeasurement - T_hat_.position()));
}

void VanytEstimator::addPositionMeasurement(const Vector3 & worldAnchorPos, const Vector3 & ImuAnchorPos)
{
  k_contacts_ = getCurrentTime();
  k_data_ = getCurrentTime();

  worldAnchorPos_ = worldAnchorPos;
  imuAnchorPos_ = ImuAnchorPos;
}

ObserverBase::StateVector VanytEstimator::oneStepEstimation_()
{
  TimeIndex k = this->x_.getTime();

  BOOST_ASSERT(this->y_.size() > 0 && this->y_.checkIndex(k + 1) && "ERROR: The measurement vector is not set");

  Vector3 yv = getMeasurement(k + 1).head<3>();
  Vector3 ya = getMeasurement(k + 1).segment<3>(3);
  Vector3 yg = getMeasurement(k + 1).segment<3>(6);

  ObserverBase::StateVector x_hat = getCurrentEstimatedState();
  x1_hat_ = x_hat.segment<3>(0);
  x2_hat_prime_ = x_hat.segment<3>(3);

  Vector dx_hat(getStateSize());
  dx_hat.setZero();
  dx_hat.segment<3>(0) = x1_hat_.cross(yg) - cst::gravityConstant * x2_hat_prime_ + ya + alpha_ * (yv - x1_hat_); // x1
  dx_hat.segment<3>(3) = x2_hat_prime_.cross(yg) - beta_ * (yv - x1_hat_); // x2_prime
  x_hat += dx_hat * dt_;

  Vector3 dt_vl = (x1_hat_ - posCorrFromContactPos_) * dt_; // using p_dot = R(v_l) = R(x1 - delta)

  sigma_ += rho1_ * (T_hat_.orientation.toMatrix3().transpose() * Vector3::UnitZ()).cross(x2_hat_prime_)
            + oriCorrFromContactPos_ + oriCorrFromOriMeas_;
  Vector3 dt_omega = (yg - sigma_) * dt_; // using R_dot = RS(w_l) = RS(yg-sigma)

  T_hat_.SE3_integration(dt_vl, dt_omega);

  x_hat.segment<3>(6) = T_hat_.position();
  x_hat.tail(4) = T_hat_.orientation.toVector4();

  /*
  // once the orientation of the IMU in the world is estimated, we can use it to estimate the position of the IMU in the
  // world
  if(k_contacts_ == k)
  {
    Vector3 worldPosFromContacts = worldAnchorPos_() - T_hat_.orientation.toMatrix3() * imuAnchorPos_();
    pos_contacts_ = expMinDtOverTau_ * pos_contacts_ + (1 - expMinDtOverTau_) * worldPosFromContacts;

    pos_x1_ = expMinDtOverTau_ * pos_x1_
              + (tau_ - expMinDtOverTau_ * tau_) * T_hat_.orientation.toMatrix3() * x_hat.segment<3>(3);

    T_hat_.position = pos_x1_ + pos_contacts_;
    x_hat.segment<3>(0) = T_hat_.position(); // pos
  }
  */

  setState(x_hat, k + 1);

  k_est_++;

  return x_hat;
}

void VanytEstimator::resetImuLocVelHat()
{
  x_().segment<3>(0) = x1_;
}

} // namespace stateObservation
