/**
 * \file      tilt-estimator.hpp
 * \author    Arnaud Demont, Mehdi Benallegue
 * \date       2018
 * \brief      Version of the Tilt Estimator that implements all the necessary functions to perform the estimation for
 * humanoid robots.
 *
 * \details
 *
 *
 */

#ifndef VanytEstimatorHPP
#define VanytEstimatorHPP

#include "state-observation/observer/tilt-estimator-humanoid.hpp"
#include <state-observation/api.h>

namespace stateObservation
{

/**
 * \class  VanytEstimator
 * \brief  Version of the Tilt Estimator for humanoid robots.
 *
 */
class STATE_OBSERVATION_DLLAPI VanytEstimator : public ZeroDelayObserver // : public TiltEstimatorHumanoid
{
  typedef kine::Orientation Orientation;

public:
  /// The constructor
  ///  \li alpha : parameter related to the convergence of the linear velocity
  ///              of the IMU expressed in the control frame
  ///  \li beta  : parameter related to the fast convergence of the tilt
  ///  \li rho  : parameter related to the orthogonality
  VanytEstimator(double alpha, double beta, double rho, double dt);

  /// @brief initializes the state vector.
  /// @param x1 The initial local linear velocity of the IMU.
  /// @param x2_p The initial value of the intermediate estimate of the IMU's tilt.
  /// @param x2 The initial tilt of the IMU.
  void initEstimator(const Vector3 & pos = Vector3::Zero(),
                     const Vector3 & x1 = Vector3::Zero(),
                     const Vector3 & x2_prime = Vector3::UnitZ(),
                     const Vector4 & R = Vector4(0, 0, 0, 1));

  /// @brief Resets x1hat (the estimate of the local linear velocity of the IMU in the world)
  /// @details Avoid discontinuities when the computation mode of the anchor point changes
  void resetImuLocVelHat();

  /// sets the measurement
  void setMeasurement(const Vector3 & yv_k, const Vector3 & ya_k, const Vector3 & yg_k, TimeIndex k);

  void setMeasurement(const ObserverBase::MeasureVector & y_k, TimeIndex k) override;

  // add an orientation measurement to the correction
  void addOrientationMeasurement(const Matrix3 & meas, double gain);

  void addContactPosMeasurement(const Vector3 & posMeasurement,
                                const Vector3 & imuContactPos,
                                double gainDelta,
                                double gainSigma);

  StateVector replayBufferedIteration(BufferedIter & bufferedIter, const std::array<double, 3> & gains);

  Vector3 getVirtualLocalVelocityMeasurement()
  {
    return x1_;
  }

  /// set the sampling time of the measurements
  virtual void setSamplingTime(const double dt)
  {
    dt_ = dt;
  }
  double getSamplingTime() const
  {
    return dt_;
  }

  /// set the gain of x1_hat variable
  void setAlpha(const double alpha)
  {
    alpha_ = alpha;
  }
  double getAlpha() const
  {
    return alpha_;
  }

  /// set the gain of x2prime_hat variable
  void setBeta(const double beta)
  {
    beta_ = beta;
  }
  double getBeta() const
  {
    return beta_;
  }

  /// set rho
  void setRho(const double rho)
  {
    rho_ = rho;
  }
  double getRho() const
  {
    return rho_;
  }

  // returns the correction term applied on the estimated orientation
  inline const stateObservation::Vector3 & getOriCorrection() const
  {
    return sigma_;
  }
  // correction of the position coming from the contact positions, passed as a local linear velocity.
  inline const stateObservation::Vector3 & getPosCorrectionFromContactPos() const
  {
    return posCorrFromContactPos_;
  }
  // correction of the orientation coming from the contact positions, passed as a local angular velocity.
  inline const stateObservation::Vector3 & geOriCorrectionFromContactPos() const
  {
    return oriCorrFromContactPos_;
  }
  // correction of the orientation coming from direct orientation measurements, passed as a local angular velocity.
  inline const stateObservation::Vector3 & getOriCorrFromOriMeas() const
  {
    return oriCorrFromOriMeas_;
  }

protected:
  /// Orientation estimator loop
  StateVector oneStepEstimation_() override;
  Eigen::Matrix<double, 12, 1> computeStateDerivatives(const ObserverBase::StateVector & x_hat,
                                                       const ObserverBase::MeasureVector & y_k);

  /// @brief integrates the given dx into the given state.
  /// @param x_hat The state to update
  /// @param T_hat The transformation matrix containing the state position and orientation (passed again under this
  /// shape for convenience)
  /// @param dx_hat The state increment to integrate
  void integrateState(ObserverBase::StateVector & x_hat,
                      kine::Kinematics & T_hat,
                      const Eigen::Matrix<double, 12, 1> & dx_hat);

  void startNewIteration();
  void resetCorrectionTerms();

  TimeIndex k_est_ = 0.0; // time index of the last estimation
  TimeIndex k_data_ = 0.0; // time index of the current measurements
  TimeIndex k_contacts_ = 0.0; // time index of the contact measurements

  /// Sampling time
  double dt_;

  /// variables used for the computation
  Vector3 x1_;
  Vector3 x1_hat_;
  Vector3 x2_hat_prime_;
};

} // namespace stateObservation

#endif // VanytEstimatorHPP
