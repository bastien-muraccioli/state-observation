/**
 * \file      tilt-estimator.hpp
 * \author    Rafael Cisneros, Mehdi Benallegue
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
  VanytEstimator(double alpha, double beta, double tau, double dt);

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

  /// sets ths measurement (accelero and gyro stacked in one vector)
  void setMeasurement(const Vector3 & yv_k, const Vector3 & ya_k, const Vector3 & yg_k, TimeIndex k);

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

  /// set rho1
  void setRho1(const double rho1)
  {
    rho1_ = rho1;
  }
  double getRho1() const
  {
    return rho1_;
  }

  /// set tau
  void setTau(const double tau)
  {
    tau_ = tau;
    expMinDtOverTau_ = exp(-dt_ / tau_);
  }
  double getTau() const
  {
    return tau_;
  }

  /// set the sampling time of the measurements
  virtual void setSamplingTime(const double dt) override
  {
    TiltEstimator::setSamplingTime(dt);
    expMinDtOverTau_ = exp(-dt_ / tau_);
  }

  // returns the correction term applied on the estimated orientation
  inline const stateObservation::Vector3 & getOriCorrection() const
  {
    return sigma_;
  }
  // returns the position pseudo-mmesurement coming from the contacts
  inline const stateObservation::Vector3 & getPosContacts() const
  {
    return pos_contacts_;
  }
  // returns the position pseudo-mmesurement coming from the integration of x1
  inline const stateObservation::Vector3 & getPosX1() const
  {
    return pos_x1_;
  }

  /// sets ths measurement (accelero and gyro stacked in one vector)
  void setMeasurement(const Vector3 & ya_k, const Vector3 & yg_k, TimeIndex k);

  /// sets ths measurement (accelero and gyro stacked in one vector)
  void setMeasurement(const Vector3 & yv_k, const Vector3 & ya_k, const Vector3 & yg_k, TimeIndex k);

  // add an orientation measurement to the correction
  void addOrientationMeasurement(const Matrix3 & meas, double gain);

  // add a position measurement to the correction
  void addPositionMeasurement(const Vector3 & worldAnchorPos, const Vector3 & ImuAnchorPos);

#if defined(__clang__)
#  pragma clang diagnostic pop
#else
#  if defined(__GNUC__)
#    pragma GCC diagnostic pop
#  endif
#endif

public:
protected:
  /// The parameters of the estimator
  double alpha_, beta_;

  /// Estimated pose of the IMU
  kine::Kinematics T_hat_;

  CheckedVector3 worldAnchorPos_;
  CheckedVector3 imuAnchorPos_;
  Vector3 pos_contacts_;
  Vector3 pos_x1_;

  double rho1_ = 0.0;

  double tau_;
  double expMinDtOverTau_;

  Vector3 sigma_ = Vector3::Zero();

  /// Orientation estimator loop
  StateVector oneStepEstimation_() override;

  void resetForNextIteration();

  /// Sampling time
  double dt_;

  /// variables used for the computation
  Vector3 x1_;
  Vector3 x1_hat_;
  Vector3 x2_hat_prime_;
};

} // namespace stateObservation

#endif // VanytEstimatorHPP
