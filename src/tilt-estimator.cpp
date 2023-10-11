#include <state-observation/observer/tilt-estimator.hpp>

namespace stateObservation
{

TiltEstimator::TiltEstimator(double alpha, double beta, double gamma)
: ZeroDelayObserver(9, 6), alpha_(alpha), beta_(beta), gamma_(gamma), dt_(0.005), p_S_C_(Vector3::Zero()),
  R_S_C_(Matrix3::Identity()), v_S_C_(Vector3::Zero()), w_S_C_(Vector3::Zero()), v_C_(Vector3::Zero())
{
}

void TiltEstimator::initEstimator(Vector3 x1, Vector3 x2_prime, Vector3 x2)
{
  Eigen::VectorXd initStateVector = Eigen::VectorXd::Zero(getStateSize());

  initStateVector.segment<3>(0) = x1;
  initStateVector.segment<3>(3) = x2_prime;
  initStateVector.segment<3>(6) = x2;

  setState(initStateVector, 0);
}

void TiltEstimator::setMeasurement(const Vector3 ya_k, const Vector3 yg_k, TimeIndex k)
{
  ObserverBase::MeasureVector y_k(6);
  y_k << ya_k, yg_k;

  ZeroDelayObserver::setMeasurement(y_k, k);
}

void TiltEstimator::setExplicitX1(const Vector3 & x1)
{
  withExplicitX1_ = true;
  x1_ = x1;
}

void TiltEstimator::resetX1hat()
{
  resetX1hat_ = true;
}

void TiltEstimator::checkResetX1hat()
{
  if(resetX1hat_)
  {
    x_().segment<3>(0) = x1_;
    resetX1hat_ = false;
  }
}

ObserverBase::StateVector TiltEstimator::oneStepEstimation_()
{
  TimeIndex k = this->x_.getTime();

  BOOST_ASSERT(this->y_.size() > 0 && this->y_.checkIndex(k + 1) && "ERROR: The measurement vector is not set");

  Vector3 ya = getMeasurement(k + 1).head<3>();
  Vector3 yg = getMeasurement(k + 1).tail<3>();

  if(!withExplicitX1_)
  {
    x1_ = R_S_C_.transpose() * (v_C_ + v_S_C_) + (yg - R_S_C_.transpose() * w_S_C_).cross(R_S_C_.transpose() * p_S_C_);
  }
  else
  {
    withExplicitX1_ = false;
  }

  checkResetX1hat();

  ObserverBase::StateVector x_hat = getCurrentEstimatedState();
  x1_hat_ = x_hat.segment<3>(0);
  x2_hat_prime_ = x_hat.segment<3>(3);
  x2_hat_ = x_hat.segment<3>(6);

  Vector dx_hat(9);
  dx_hat.segment<3>(0) = x1_hat_.cross(yg) - cst::gravityConstant * x2_hat_prime_ + ya + alpha_ * (x1_ - x1_hat_);
  dx_hat.segment<3>(3) = x2_hat_prime_.cross(yg) - beta_ * (x1_ - x1_hat_);
  dx_hat.segment<3>(6) = x2_hat_.cross(yg - gamma_ * x2_hat_.cross(x2_hat_prime_));

  x_hat += dx_hat * dt_;

  x_hat.tail<3>() /= x_hat.tail<3>().norm();

  setState(x_hat, k + 1);

  return x_hat;
}

} // namespace stateObservation
