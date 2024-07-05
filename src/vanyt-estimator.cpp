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
    resetCorrectionTerms();
  }
}

void VanytEstimator::resetCorrectionTerms()
{
  sigma_.setZero();
  oriCorrFromOriMeas_.setZero();
  posCorrFromContactPos_.setZero();
  oriCorrFromContactPos_.setZero();
}

void VanytEstimator::setMeasurement(const Vector3 & yv_k, const Vector3 & ya_k, const Vector3 & yg_k, TimeIndex k)
{
  startNewIteration();

  ObserverBase::MeasureVector y_k(getMeasureSize());
  y_k << yv_k, ya_k, yg_k;

  ZeroDelayObserver::setMeasurement(y_k, k);
}

void VanytEstimator::setMeasurement(const ObserverBase::MeasureVector & y_k, TimeIndex k)
{
  startNewIteration();

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

ObserverBase::StateVector VanytEstimator::oneStepEstimation_()
{
  TimeIndex k = this->x_.getTime();

  BOOST_ASSERT(this->y_.size() > 0 && this->y_.checkIndex(k + 1) && "ERROR: The measurement vector is not set");

  ObserverBase::StateVector x_hat = getCurrentEstimatedState();

  Eigen::Matrix<double, 12, 1> dx_hat = computeStateDerivatives(x_hat, getMeasurement(k + 1));
  integrateState(x_hat, T_hat_, dx_hat);

  setState(x_hat, k + 1);

  k_est_++;

  return x_hat;
}

ObserverBase::StateVector VanytEstimator::replayBufferedIteration(BufferedIter & bufferedIter,
                                                                  const std::array<double, 3> & gains)
{
  double prevAlpha = alpha_;
  double prevBeta = beta_;
  double prevRho = rho_;

  setAlpha(gains.at(0));
  setBeta(gains.at(1));
  setRho(gains.at(2));

  StateVector integratedState = bufferedIter.x_k_;
  kine::Kinematics T_hat(bufferedIter.x_k_.tail(7), kine::Kinematics::Flags::pose);

  setMeasurement(bufferedIter.y_k_, getCurrentTime() + 1);
  for(auto & contactPosMeas : bufferedIter.contactPosMeasurements_)
  {
    addContactPosMeasurement(contactPosMeas.worldContactRefPos_, contactPosMeas.imuContactPos_, contactPosMeas.lambda_,
                             contactPosMeas.gamma_);
  }
  for(auto & oriMeas : bufferedIter.oriMeasurements_)
  {
    addOrientationMeasurement(oriMeas.measuredOri_, oriMeas.gain_);
  Eigen::Matrix<double, 12, 1> dx_hat = computeStateDerivatives(bufferedIter.x_k_, bufferedIter.y_k_);
  integrateState(integratedState, T_hat, dx_hat);

  resetCorrectionTerms();

  setAlpha(prevAlpha);
  setBeta(prevBeta);
  setRho(prevRho);

  return integratedState;
}

Eigen::Matrix<double, 12, 1> VanytEstimator::computeStateDerivatives(const ObserverBase::StateVector & x_hat,
                                                                     const ObserverBase::MeasureVector & y_k)
{
  const Vector3 & yv = y_k.head<3>();
  const Vector3 & ya = y_k.segment<3>(3);
  const Vector3 & yg = y_k.segment<3>(6);

  x1_hat_ = x_hat.segment<3>(0);
  x2_hat_prime_ = x_hat.segment<3>(3);

  Eigen::Matrix<double, 12, 1> dx_hat;
  dx_hat.segment<3>(0) = x1_hat_.cross(yg) - cst::gravityConstant * x2_hat_prime_ + ya + alpha_ * (yv - x1_hat_); // x1
  dx_hat.segment<3>(3) = x2_hat_prime_.cross(yg) - beta_ * (yv - x1_hat_); // x2_prime

  dx_hat.segment<3>(6) = (x1_hat_ - posCorrFromContactPos_); // using p_dot = R(v_l) = R(x1 - delta)

  sigma_ += rho_ * (T_hat_.orientation.toMatrix3().transpose() * Vector3::UnitZ()).cross(x2_hat_prime_)
            + oriCorrFromContactPos_ + oriCorrFromOriMeas_;
  dx_hat.segment<3>(9) = (yg - sigma_); // using R_dot = RS(w_l) = RS(yg-sigma)

  return dx_hat;
}

void VanytEstimator::integrateState(ObserverBase::StateVector & x_hat,
                                    kine::Kinematics & T_hat,
                                    const Eigen::Matrix<double, 12, 1> & dx_hat)
{
  const Vector3 & vl = dx_hat.segment<3>(6);
  const Vector3 & omega = dx_hat.segment<3>(9);

  // discrete-time integration of x1 and x2
  x_hat.segment<6>(0) += dx_hat.segment<6>(0) * dt_;

  // discrete-time integration of p and R
  T_hat.SE3_integration(vl * dt_, omega * dt_);
  std::cout << std::endl << "Tafter: " << T_hat.orientation.toRotationVector().transpose() << std::endl;
  x_hat.segment<3>(6) = T_hat.position();
  x_hat.tail(4) = T_hat.orientation.toVector4();
}

void VanytEstimator::resetImuLocVelHat()
{
  x_().segment<3>(0) = x1_;
}

} // namespace stateObservation
