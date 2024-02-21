#include "state-observation/tools/definitions.hpp"
#include "state-observation/tools/rigid-body-kinematics.hpp"
#include "state-observation/tools/rigid-body-kinematics.hxx"
#include <state-observation/observer/tilt-visual-humanoid.hpp>

namespace stateObservation
{

TiltVisualHumanoid::TiltVisualHumanoid(double alpha, double beta) : TiltEstimatorHumanoid(alpha, beta, 0.0, 10, 13) {}

void TiltVisualHumanoid::initEstimator(Vector3 x1, Vector3 x2_prime, Vector4 R)
{
  Eigen::VectorXd initStateVector = Eigen::VectorXd::Zero(getStateSize());

  initStateVector.segment<3>(0) = x1;
  initStateVector.segment<3>(3) = x2_prime;
  initStateVector.tail(4) = R;

  setState(initStateVector, 0);
}

void TiltVisualHumanoid::setMeasurement(const Vector3 & ya_k, const Vector3 & yg_k, const Vector4 & yR_k, TimeIndex k)
{
  x1_ = R_S_C_.transpose() * (v_C_ + v_S_C_) + (yg_k - R_S_C_.transpose() * w_S_C_).cross(R_S_C_.transpose() * p_S_C_);
  ObserverBase::MeasureVector y_k(getMeasureSize());

  y_k << x1_, ya_k, yg_k, yR_k;

  ZeroDelayObserver::setMeasurement(y_k, k);
}

void TiltVisualHumanoid::setMeasurement(const Vector3 & yv_k,
                                        const Vector3 & ya_k,
                                        const Vector3 & yg_k,
                                        const Vector4 & yR_k,
                                        TimeIndex k)
{
  ObserverBase::MeasureVector y_k(getMeasureSize());
  y_k << yv_k, ya_k, yg_k, yR_k;

  ZeroDelayObserver::setMeasurement(y_k, k);
}

void TiltVisualHumanoid::resetImuLocVelHat()
{
  resetX1hat_ = true;
}

void TiltVisualHumanoid::checkResetX1hat()
{
  if(resetX1hat_)
  {
    x_().segment<3>(0) = x1_;
    resetX1hat_ = false;
  }
}

ObserverBase::StateVector TiltVisualHumanoid::oneStepEstimation_()
{
  TimeIndex k = this->x_.getTime();

  BOOST_ASSERT(this->y_.size() > 0 && this->y_.checkIndex(k + 1) && "ERROR: The measurement vector is not set");

  Vector3 yv = getMeasurement(k + 1).head<3>();
  Vector3 ya = getMeasurement(k + 1).segment<3>(3);
  Vector3 yg = getMeasurement(k + 1).segment<3>(6);
  Orientation yR;
  yR.fromVector4(getMeasurement(k + 1).tail(4));

  ObserverBase::StateVector x_hat = getCurrentEstimatedState();
  x1_hat_ = x_hat.segment<3>(0);
  x2_hat_prime_ = x_hat.segment<3>(3);
  R_hat_.fromVector4(x_hat.tail(4));

  Vector dx_hat(getStateSize());
  dx_hat.setZero();
  dx_hat.segment<3>(0) = x1_hat_.cross(yg) - cst::gravityConstant * x2_hat_prime_ + ya + alpha_ * (yv - x1_hat_);
  dx_hat.segment<3>(3) = x2_hat_prime_.cross(yg) - beta_ * (yv - x1_hat_);

  Matrix3 rot_diff = yR.toMatrix3().transpose() * R_hat_.toMatrix3();
  Vector3 rot_diff_vec = kine::skewSymmetricToRotationVector(rot_diff - rot_diff.transpose());

  // Matrix3 rot_diff = R_hat_.toMatrix3() * yR.toMatrix3().transpose();
  Matrix3 Rhat_t = R_hat_.toMatrix3().transpose();

  sigma_part1_ = rho1_ * (Rhat_t * Vector3::UnitZ()).cross(x2_hat_prime_);
  sigma_part2_ = rho2_ * rot_diff_vec;
  sigma_part3_ = mu_ * Rhat_t * Vector3::UnitZ() * (Rhat_t * Vector3::UnitZ()).transpose() * rot_diff_vec;

  Vector3 sigma = sigma_part1_ + sigma_part2_ + sigma_part3_;

  Matrix3 dR_hat = R_hat_.toMatrix3() * kine::skewSymmetric(yg - sigma);

  Vector3 dt_omega = kine::skewSymmetricToRotationVector(dR_hat * R_hat_.toMatrix3().transpose()) * dt_;
  x_hat += dx_hat * dt_;
  x_hat.segment<3>(6) /= x_hat.segment<3>(6).norm();

  R_hat_.integrate(dt_omega);

  x_hat.tail(4) = R_hat_.toVector4();

  setState(x_hat, k + 1);

  return x_hat;
}

} // namespace stateObservation
