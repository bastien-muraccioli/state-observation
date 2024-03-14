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
class STATE_OBSERVATION_DLLAPI VanytEstimator : public TiltEstimatorHumanoid
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

  /// sets the position of the IMU sensor in the control frame
  void setSensorPositionInC(const Vector3 & p)
  {
    p_S_C_ = p;
  }

  Vector3 getSensorPositionInC()
  {
    return p_S_C_;
  }

  /// sets the oriantation of the IMU sensor in the control frame
  void setSensorOrientationInC(const Matrix3 & R)
  {
    R_S_C_ = R;
  }
  Matrix3 getSensorOrientationInC()
  {
    return R_S_C_;
  }

  Vector3 getVirtualLocalVelocityMeasurement()
  {
    return x1_;
  }

  /// sets teh linear velocity of the IMU sensor in the control frame
  void setSensorLinearVelocityInC(const Vector3 & v)
  {
    v_S_C_ = v;
  }

  Vector3 getSensorLinearVelocityInC()
  {
    return v_S_C_;
  }

  /// sets the angular velocity of the IMU sensor in the control frame
  void setSensorAngularVelocityInC(const Vector3 & w)
  {
    w_S_C_ = w;
  }
  Vector3 getSensorAngularVelocityInC()
  {
    return w_S_C_;
  }

  /// sets the velocity of the control origin in the world frame
  /// this velocity has to be expressed in the control frame.
  void setControlOriginVelocityInW(const Vector3 & v)
  {
    v_C_ = v;
  }
  Vector3 getControlOriginVelocityInW()
  {
    return v_C_;
  }

/// prevent c++ overloaded virtual function warning
#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Woverloaded-virtual"
#else
#  if defined(__GNUC__)
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Woverloaded-virtual"
#  endif
#endif

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

  // returns the correction term applied on the estimated orientation
  inline const stateObservation::Vector3 & getOriCorrection() const
  {
    return oriCorrection_;
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
  /// Position of the IMU in the control frame
  Vector3 p_S_C_;

  /// Orientation of the IMU in the control frame
  Matrix3 R_S_C_;

  /// Linear velocity of the IMU in the control frame
  Vector3 v_S_C_;

  /// Angular velocity of the IMU in the control frame
  Vector3 w_S_C_;

  /// Linear velocity of the control frame
  Vector3 v_C_;

  CheckedVector3 worldAnchorPos_;
  CheckedVector3 imuAnchorPos_;
  Vector3 pos_contacts_;
  Vector3 pos_x1_;

  double rho1_ = 0.0;

  double tau_;
  double expMinDtOverTau_;

  Vector3 oriCorrection_ = Vector3::Zero();

  /// Orientation estimator loop
  StateVector oneStepEstimation_();

  void startNewIteration();

  Index k_data_ = 0.0;
};

} // namespace stateObservation

#endif // VanytEstimatorHPP
