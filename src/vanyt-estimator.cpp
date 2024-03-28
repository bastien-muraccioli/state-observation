#include "state-observation/tools/definitions.hpp"
#include <state-observation/observer/vanyt-estimator.hpp>

namespace stateObservation
{

VanytEstimator::VanytEstimator(double alpha, double beta, double tau, double dt)
: TiltEstimatorHumanoid(alpha, beta, 0.0, 13, 9, dt)
{
  setTau(tau);
}

void VanytEstimator::initEstimator(const Vector3 & pos, const Vector3 & x1, const Vector3 & x2_prime, const Vector4 & R)
{
  Eigen::VectorXd initStateVector = Eigen::VectorXd::Zero(getStateSize());

  initStateVector.segment<3>(0) = pos;
  initStateVector.segment<3>(3) = x1;
  initStateVector.segment<3>(6) = x2_prime;
  initStateVector.tail(4) = R;

  setState(initStateVector, 0);

  pos_contacts_ = pos;
  pos_x1_.setZero();
}

void VanytEstimator::resetForNextIteration()
{
  imuAnchorPos_.reset();
  worldAnchorPos_.reset();

  sigma_.setZero();
}

void VanytEstimator::setMeasurement(const Vector3 & ya_k, const Vector3 & yg_k, TimeIndex k)
{
  x1_ = R_S_C_.transpose() * (v_C_ + v_S_C_) + (yg_k - R_S_C_.transpose() * w_S_C_).cross(R_S_C_.transpose() * p_S_C_);
  ObserverBase::MeasureVector y_k(getMeasureSize());

  y_k << x1_, ya_k, yg_k;

  ZeroDelayObserver::setMeasurement(y_k, k);
}

void VanytEstimator::setMeasurement(const Vector3 & yv_k, const Vector3 & ya_k, const Vector3 & yg_k, TimeIndex k)
{
  ObserverBase::MeasureVector y_k(getMeasureSize());
  y_k << yv_k, ya_k, yg_k;

  ZeroDelayObserver::setMeasurement(y_k, k);
}

void VanytEstimator::addOrientationMeasurement(const Matrix3 & oriMeasurement, double gain)
{
  Matrix3 rot_diff = R_hat_.toMatrix3() * oriMeasurement.transpose();
  Vector3 rot_diff_vec = kine::skewSymmetricToRotationVector(rot_diff - rot_diff.transpose());

  sigma_ += gain * R_hat_.toMatrix3().transpose() * Vector3::UnitZ() * (Vector3::UnitZ()).transpose() * rot_diff_vec;
}

void VanytEstimator::addPositionMeasurement(const Vector3 & worldAnchorPos, const Vector3 & ImuAnchorPos)
{
  k_contacts_ = getCurrentTime();

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
  x1_hat_ = x_hat.segment<3>(3);
  x2_hat_prime_ = x_hat.segment<3>(6);
  R_hat_.fromVector4(x_hat.tail(4));

  Vector dx_hat(getStateSize());
  dx_hat.setZero();
  dx_hat.segment<3>(3) = x1_hat_.cross(yg) - cst::gravityConstant * x2_hat_prime_ + ya + alpha_ * (yv - x1_hat_); // x1
  dx_hat.segment<3>(6) = x2_hat_prime_.cross(yg) - beta_ * (yv - x1_hat_); // x2_prime
  x_hat += dx_hat * dt_;

  sigma_ += rho1_ * (R_hat_.toMatrix3().transpose() * Vector3::UnitZ()).cross(x2_hat_prime_);
  Vector3 dt_omega = (yg - sigma_) * dt_; // using R_dot = RS(w_l) = RS(yg-sigma)

  R_hat_.integrateRightSide(dt_omega);
  x_hat.tail(4) = R_hat_.toVector4();

  // once the orientation of the IMU in the world is estimated, we can use it to estimate the position of the IMU in the
  // world
  if(k_contacts_ == k)
  {
    Vector3 worldPosFromContacts = worldAnchorPos_() - R_hat_.toMatrix3() * imuAnchorPos_();
    pos_contacts_ = expMinDtOverTau_ * pos_contacts_ + (1 - expMinDtOverTau_) * worldPosFromContacts;

    pos_x1_ = expMinDtOverTau_ * pos_x1_ + (tau_ - expMinDtOverTau_ * tau_) * R_hat_.toMatrix3() * x_hat.segment<3>(3);

    x_hat.segment<3>(0) = pos_x1_ + pos_contacts_; // pos
  }

  setState(x_hat, k + 1);

  resetForNextIteration();

  return x_hat;
}

} // namespace stateObservation
