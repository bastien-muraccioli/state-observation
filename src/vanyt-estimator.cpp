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

void VanytEstimator::startNewIteration()
{
  // Checks if new data was passed to the estimator since the last estimation
  if(k_data_ == getCurrentTime())
  {
    oriCorrection_.setZero();

    imuAnchorPos_.reset();
    worldAnchorPos_.reset();

    ++k_data_;
  }
}

void VanytEstimator::setMeasurement(const Vector3 & ya_k, const Vector3 & yg_k, TimeIndex k)
{
  startNewIteration();
  x1_ = R_S_C_.transpose() * (v_C_ + v_S_C_) + (yg_k - R_S_C_.transpose() * w_S_C_).cross(R_S_C_.transpose() * p_S_C_);
  ObserverBase::MeasureVector y_k(getMeasureSize());

  y_k << x1_, ya_k, yg_k;

  ZeroDelayObserver::setMeasurement(y_k, k);
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
  Matrix3 rot_diff = oriMeasurement * R_hat_.toMatrix3().transpose();
  Vector3 rot_diff_vec = kine::skewSymmetricToRotationVector(rot_diff - rot_diff.transpose());

  oriCorrection_ += gain * R_hat_.toMatrix3().transpose() * Vector3::UnitZ()
                    * (R_hat_.toMatrix3().transpose() * Vector3::UnitZ()).transpose() * rot_diff_vec;
}

void VanytEstimator::addPositionMeasurement(const Vector3 & worldAnchorPos, const Vector3 & ImuAnchorPos)
{
  startNewIteration();

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
  x_hat.segment<3>(6) /= x_hat.segment<3>(6).norm();

  oriCorrection_ += rho1_ * (R_hat_.toMatrix3().transpose() * Vector3::UnitZ()).cross(x2_hat_prime_);
  Matrix3 dR_hat = R_hat_.toMatrix3() * kine::skewSymmetric(yg - oriCorrection_);
  Vector3 dt_omega = kine::skewSymmetricToRotationVector(dR_hat * R_hat_.toMatrix3().transpose()) * dt_;
  R_hat_.integrate(dt_omega);
  x_hat.tail(4) = R_hat_.toVector4();

  // once the orientation of the IMU in the world is estimated, we can use it to estimate the position of the IMU in the
  // world
  if(worldAnchorPos_.isSet())
  {
    Vector3 worldPosFromContacts = worldAnchorPos_() - R_hat_.toMatrix3() * imuAnchorPos_();
    pos_contacts_ = expMinDtOverTau_ * pos_contacts_ + (1 - expMinDtOverTau_) * worldPosFromContacts;

    pos_x1_ = expMinDtOverTau_ * pos_x1_ + (tau_ - expMinDtOverTau_ * tau_) * R_hat_.toMatrix3() * x_hat.segment<3>(3);

    x_hat.segment<3>(0) = pos_x1_ + pos_contacts_; // pos
  }

  setState(x_hat, k + 1);

  return x_hat;
}

} // namespace stateObservation
