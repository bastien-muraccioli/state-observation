#include "state-observation/tools/definitions.hpp"
#include <state-observation/observer/vanyt-estimator.hpp>

namespace stateObservation
{
VanytEstimator::VanytEstimator(double alpha, double beta, double rho, double dt)
: ZeroDelayObserver(13, 9), iterInfos_(alpha, beta, rho, dt)
{
  iterInfos_.alpha_ = alpha;
  iterInfos_.beta_ = beta;
  iterInfos_.rho_ = rho;
  iterInfos_.dt_ = dt;
}

void VanytEstimator::initEstimator(const Vector3 & pos, const Vector3 & x1, const Vector3 & x2_prime, const Vector4 & R)
{
  Eigen::VectorXd initStateVector = Eigen::VectorXd::Zero(getStateSize());

  initStateVector.segment<3>(0) = x1;
  initStateVector.segment<3>(3) = x2_prime;
  initStateVector.segment<3>(6) = pos;
  initStateVector.tail(4) = R;

  setState(initStateVector, 0);

  getCurrentIter().initState_ = initStateVector;
  getCurrentIter().updatedState_ = initStateVector;

  getCurrentIter().currentPose_.position = pos;
  getCurrentIter().currentPose_.orientation.fromVector4(R);
}

void VanytEstimator::IterInfos::startNewIteration()
{
  if(k_est_ == k_data_)
  {
    ++k_data_;

    iter_++;
    initState_ = updatedState_;
    resetCorrectionTerms();
  }
}

void VanytEstimator::IterInfos::resetCorrectionTerms()
{
  sigma_.setZero();
  oriCorrFromOriMeas_.setZero();
  posCorrFromContactPos_.setZero();
  oriCorrFromContactPos_.setZero();
}

void VanytEstimator::setMeasurement(const Vector3 & yv_k,
                                    const Vector3 & ya_k,
                                    const Vector3 & yg_k,
                                    TimeIndex k,
                                    bool resetImuLocVelHat)
{
  getCurrentIter().startNewIteration();

  ObserverBase::MeasureVector y_k(getMeasureSize());
  y_k << yv_k, ya_k, yg_k;

  setMeasurement(y_k, k);

  if(resetImuLocVelHat)
  {
    x_().segment<3>(0) = yv_k;
  }
}

void VanytEstimator::setMeasurement(const ObserverBase::MeasureVector & y_k, TimeIndex k)
{
  getCurrentIter().startNewIteration();

  ZeroDelayObserver::setMeasurement(y_k, k);

  getCurrentIter().saveMeasurement(getMeasurement(getMeasurementTime()));
}

void VanytEstimator::IterInfos::addOrientationMeasurement(const Matrix3 & oriMeasurement, double gain)
{
  startNewIteration();

  Matrix3 rot_diff = oriMeasurement * currentPose_.orientation.toMatrix3().transpose();
  Vector3 rot_diff_vec = kine::skewSymmetricToRotationVector(rot_diff - rot_diff.transpose());
  oriCorrFromOriMeas_ -= gain * currentPose_.orientation.toMatrix3().transpose() * Vector3::UnitZ()
                         * (Vector3::UnitZ()).transpose() * rot_diff_vec;
}

void VanytEstimator::addContactPosMeasurement(const Vector3 & posMeasurement,
                                              const Vector3 & imuContactPos,
                                              double gainDelta,
                                              double gainSigma)
{
  getCurrentIter().addContactPosMeasurement(posMeasurement, imuContactPos, gainDelta, gainSigma);
}

void VanytEstimator::addOrientationMeasurement(const Matrix3 & oriMeasurement, double gain)
{
  getCurrentIter().addOrientationMeasurement(oriMeasurement, gain);
}

void VanytEstimator::IterInfos::addContactPosMeasurement(const Vector3 & posMeasurement,
                                                         const Vector3 & imuContactPos,
                                                         double gainDelta,
                                                         double gainSigma)
{
  startNewIteration();

  oriCorrFromContactPos_ +=
      gainSigma
      * (currentPose_.orientation.toMatrix3().transpose() * (posMeasurement - currentPose_.position()))
            .cross(imuContactPos);

  posCorrFromContactPos_ +=
      gainDelta
      * (imuContactPos - currentPose_.orientation.toMatrix3().transpose() * (posMeasurement - currentPose_.position()));
}

ObserverBase::StateVector VanytEstimator::oneStepEstimation_()
{
  TimeIndex k = this->x_.getTime();
  IterInfos & currentIter = getCurrentIter();

  BOOST_ASSERT(this->y_.size() > 0 && this->y_.checkIndex(k + 1) && "ERROR: The measurement vector is not set");

  ObserverBase::StateVector x_hat = getCurrentEstimatedState();

  Eigen::Matrix<double, 12, 1> dx_hat = currentIter.computeStateDerivatives();
  currentIter.integrateState(dx_hat);

  setState(currentIter.updatedState_, k + 1);

  getCurrentIter().k_est_++;
  if(withDelayedOri_)
  {
    bufferedIters_.push_front(currentIter);
  }

  return x_hat;
}

ObserverBase::StateVector VanytEstimator::IterInfos::replayBufferedIteration()
{
  StateVector integratedState = initState_;
  currentPose_.fromVector(initState_.tail(7), kine::Kinematics::Flags::pose);

  Eigen::Matrix<double, 12, 1> dx_hat = computeStateDerivatives();
  integrateState(dx_hat);

  return integratedState;
}

Eigen::Matrix<double, 12, 1> VanytEstimator::IterInfos::computeStateDerivatives()
{
  const Vector3 & yv = y_k_.head<3>();
  const Vector3 & ya = y_k_.segment<3>(3);
  const Vector3 & yg = y_k_.segment<3>(6);

  Vector3 & oriCorrFromOriMeas = oriCorrFromOriMeas_;
  Vector3 & oriCorrFromContactPos = oriCorrFromContactPos_;
  Vector3 & posCorrFromContactPos = posCorrFromContactPos_;

  const Eigen::Ref<Vector3> x1_hat = initState_.segment<3>(0);
  const Eigen::Ref<Vector3> x2_hat_prime = initState_.segment<3>(3);

  Eigen::Matrix<double, 12, 1> dx_hat;
  dx_hat.segment<3>(0) = x1_hat.cross(yg) - cst::gravityConstant * x2_hat_prime + ya + alpha_ * (yv - x1_hat); // x1
  dx_hat.segment<3>(3) = x2_hat_prime.cross(yg) - beta_ * (yv - x1_hat); // x2_prime

  dx_hat.segment<3>(6) = (x1_hat - posCorrFromContactPos); // using p_dot = R(v_l) = R(x1 - delta)

  sigma_ = rho_ * (currentPose_.orientation.toMatrix3().transpose() * Vector3::UnitZ()).cross(x2_hat_prime)
           + oriCorrFromContactPos + oriCorrFromOriMeas;
  dx_hat.segment<3>(9) = (yg - sigma_); // using R_dot = RS(w_l) = RS(yg-sigma)

  return dx_hat;
}

void VanytEstimator::IterInfos::integrateState(const Eigen::Matrix<double, 12, 1> & dx_hat)
{
  const Vector3 & vl = dx_hat.segment<3>(6);
  const Vector3 & omega = dx_hat.segment<3>(9);

  updatedState_ = initState_;
  kine::Kinematics & T_hat = currentPose_;

  // discrete-time integration of x1 and x2
  updatedState_.segment<6>(0) += dx_hat.segment<6>(0) * dt_;

  // discrete-time integration of p and R
  T_hat.SE3_integration(vl * dt_, omega * dt_);
  std::cout << std::endl << "Tafter: " << T_hat.orientation.toRotationVector().transpose() << std::endl;
  updatedState_.segment<3>(6) = T_hat.position();
  updatedState_.tail(4) = T_hat.orientation.toVector4();
}

ObserverBase::StateVector VanytEstimator::replayIterationWithDelayedOri(unsigned long delay,
                                                                        const Matrix3 & meas,
                                                                        double gain)
{
  BOOST_ASSERT_MSG(withDelayedOri_, "The mode allowing to deal with delayed orientations has not been switched on.");

  IterInfos & bufferedIter = bufferedIters_.at(delay - 1);
  bufferedIter.addOrientationMeasurement(meas, gain);

  StateVector replayedState = bufferedIter.replayBufferedIteration();
  return replayedState;
}

ObserverBase::StateVector VanytEstimator::replayIterationWithDelayedOri(unsigned long delay,
                                                                        const Matrix3 & meas,
                                                                        double gain,
                                                                        kine::Kinematics & poseTransfo)
{
  BOOST_ASSERT_MSG(withDelayedOri_, "The mode allowing to deal with delayed orientations has not been switched on.");

  Eigen::Ref<Eigen::Matrix<double, 7, 1>> latestStatePose = getCurrentEstimatedState().tail(7);

  // extracting the buffered iteration
  IterInfos & bufferedIter = bufferedIters_.at(delay - 1);

  // kinematics that were obtained after the integration without the orientation measurement.
  const kine::Kinematics kineEstWithoutOri = bufferedIter.currentPose_;

  // we add the delayed orientation measurement to the inputs of the buffered iteration and recompute the state update.
  bufferedIter.addOrientationMeasurement(meas, gain);
  StateVector replayedEstimation = bufferedIter.replayBufferedIteration();

  // kinematics that are obtained after the integration with the orientation measurement.
  const kine::Kinematics & kineEstWithOri = bufferedIter.currentPose_;
  // transformation coming from the orientation correction
  poseTransfo = kineEstWithoutOri.getInverse() * kineEstWithOri;

  // we apply the previously computed transformation to the new estimation of the past pose
  kine::Kinematics & latestKineWithOriMeas = getCurrentIter().currentPose_;
  latestKineWithOriMeas = poseTransfo * bufferedIter.currentPose_;

  latestStatePose = latestKineWithOriMeas.toVector(kine::Kinematics::Flags::pose);

  return replayedEstimation;
}

} // namespace stateObservation
