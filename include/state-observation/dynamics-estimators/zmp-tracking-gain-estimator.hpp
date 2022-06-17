
///\file      lipm-dcm-estimator.hpp
///\author    Mehdi Benallegue
///\date      2020
///\brief     Filtering of divergent component of motion (DCM) and estimation of a bias betweeen the DCM
///           and the corresponding zero moment point for a linearized inverted
///           pendulum model
///
///\detail
///
///
#ifndef ZMPTRACKINGGAINESTIMATOR_HPP
#define ZMPTRACKINGGAINESTIMATOR_HPP

#include <state-observation/api.h>
#include <state-observation/observer/linear-kalman-filter.hpp>
#include <state-observation/tools/miscellaneous-algorithms.hpp>
#include <state-observation/tools/rigid-body-kinematics.hpp>

namespace stateObservation
{

/// \class order1-gain-estimator
/// \brief tobedone
///

class STATE_OBSERVATION_DLLAPI ZmpTrackingGainEstimator
{
private:
  constexpr static double defaultDt_ = 0.005;
  constexpr static double defaultGain_ = 1;

public:
  /// default expected drift of the bias every second
  constexpr static double defaultGainDriftSecond = 0.002;

  /// default error in the measurements of the sensor
  constexpr static double defaultZmpErrorStd = 0.005;

  /// default value for gain minimum, the gain has to be a positive definite matrix with the smallest eigenvalue bigger
  /// than this value
  constexpr static double defaultGainMinimum = .01;

  /// default valu for the initial ZMP uncertainty. It should be quite low since this is supposed to be measured
  constexpr static double defaultZmpUncertainty = 0.005;

  /// default value for the uncertainty of the Gain
  constexpr static double defaultGainUncertainty = 20;

  /// @brief Construct a new ZMP Gain Estimator object
  ///
  /// @param dt                     the sampling time in seconds
  /// @param gainDriftPerSecond     the standard deviation of the gain drift (1/s)
  /// @param zmpMeasureErrorStd     the standard deviation of the zmp estimation error (m)
  /// @param minimumGain            minimum, the gain has to be a positive definite matrix with the smallest
  /// eigenvalue bigger than this value
  /// @param initZMP                the initial value of the ZMP (m)
  /// @param initGain               the initial value of the Gain (Vector3 :two diagonals first then the nondiagonal
  /// term)
  /// @param biasLimit              the X and Y (expressed in local frame) largest accepted absolute values of the bias
  ///                               (zero means no limit)
  /// @param initZMPUncertainty     the uncertainty in the ZMP initial value in meters
  /// @param initGainUncertainty    the uncertainty in the Gain initial value in (1/s)
  ZmpTrackingGainEstimator(double dt = defaultDt_,
                           double gainDriftPerSecond = defaultGainDriftSecond,
                           double zmpMeasureErrorStd = defaultZmpErrorStd,
                           double minimumGain = defaultGainMinimum,
                           const Vector2 & initZMP = Vector2::Zero(),
                           const Vector3 & initGain = Vector3::Zero(),
                           const Vector2 & initZMPUncertainty = Vector2::Constant(defaultZmpUncertainty),
                           const Vector3 & initGainUncertainty = Vector3::Constant(defaultGainUncertainty));

  /// @brief Resets the estimator
  ///
  /// @param initZMP                the initial value of the ZMP (m)
  /// @param initGain               the initial value of the Gain (Vector3 :two diagonals first then the nondiagonal
  /// term)
  /// @param yaw                    the yaw expressed as a 2d matrix (a rotation in the X Y plane)
  /// @param zmpMeasureErrorStd     the standard deviation of the zmp estimation error (m)
  /// @param initZMPUncertainty     the uncertainty in the ZMP initial value in meters
  /// @param initGainUncertainty    the uncertainty in the Gain initial value in (1/s)
  void resetWithMeasurements(const Vector2 & initZMP = Vector2::Zero(),
                             const Vector3 & initGain = Vector3::Zero(),
                             const Matrix2 & yaw = Matrix2::Identity(),
                             const Vector2 & initZMPUncertainty = Vector2::Constant(defaultZmpUncertainty),
                             const Vector3 & initGainUncertainty = Vector3::Constant(defaultGainUncertainty));

  /// @brief Resets the estimator
  ///
  /// @param initZMP                the initial value of the ZMP (m)
  /// @param initGain               the initial value of the Gain (Vector3 :two diagonals first then the nondiagonal
  /// term)
  /// @param yaw                    the yaw angle
  /// @param zmpMeasureErrorStd     the standard deviation of the zmp estimation error (m)
  /// @param initZMPUncertainty     the uncertainty in the ZMP initial value in meters
  /// @param initGainUncertainty    the uncertainty in the Gain initial value in (1/s)
  inline void resetWithMeasurements(const Vector2 & initZMP,
                                    const Vector3 & initGain,
                                    double yaw,
                                    const Vector2 & initZMPUncertainty = Vector2::Constant(defaultZmpUncertainty),
                                    const Vector3 & initGainUncertainty = Vector3::Constant(defaultGainUncertainty))
  {
    resetWithMeasurements(initZMP, initGain, Rotation2D(yaw).toRotationMatrix(), initZMPUncertainty,
                          initGainUncertainty);
  }

  /// @brief Resets the estimator with first measurements
  /// @details Use this when initializing with an available DCM (biased / or not) measurement
  ///
  /// @param initZMP                the initial value of the ZMP (m)
  /// @param initGain               the initial value of the Gain (Vector3 :two diagonals first then the nondiagonal
  /// term)
  /// @param rotation                the 3d orientation from which the initial yaw angle will be extracted using the
  /// angle agnostic approach. This orientation is from local to global. i.e. zmp_global == orientation * zmp_local
  /// @param zmpMeasureErrorStd     the standard deviation of the zmp estimation error (m)
  /// @param initZMPUncertainty     the uncertainty in the ZMP initial value in meters
  /// @param initGainUncertainty    the uncertainty in the Gain initial value in (1/s)
  inline void resetWithMeasurements(const Vector2 & initZMP,
                                    const Vector3 & initGain,
                                    const Matrix3 & rotation,
                                    const Vector2 & initZMPUncertainty = Vector2::Constant(defaultZmpUncertainty),
                                    const Vector3 & initGainUncertainty = Vector3::Constant(defaultGainUncertainty))
  {
    resetWithMeasurements(initZMP, initGain, kine::rotationMatrixToYawAxisAgnostic(rotation), initZMPUncertainty,
                          initGainUncertainty);
  }

  ///@brief Destroy the Lipm Dcm Bias Estimator object
  ~ZmpTrackingGainEstimator();

  ///@brief Set the Sampling Time
  ///
  ///@param dt sampling time
  void setSamplingTime(double dt);

  ///@brief Set the Gain from a guess
  ///
  ///@param bias guess
  void setGain(const Vector2 & bias);

  ///@copydoc setGain(double bias)
  ///
  ///@param the uncertainty you have in this guess in meters
  void setGain(const Vector2 & bias, const Vector2 & uncertainty);

  /// @brief Set the Gain Drift Per Second
  ///
  /// @param driftPerSecond the standard deviation of the gain
  void setGainDriftPerSecond(double driftPerSecond);

  /// @brief Set the Gain Limit
  ///
  /// @param minGain the minimal value of the gain
  void setMinimumGain(const double & minGain);

  /// @brief Set the Zmp Measurement Error Stamdard deviation
  ///
  void setZmpMeasureErrorStd(double);

  /// @brief Set the Inputs of the estimator.
  /// @details The yaw will be extracted from the orientation using the axis agnostic
  /// approach.
  ///
  /// @param zmperr      measurement of zmp error (zmp-zmp^{ref})
  /// @param zmp         mesaurement of the ZMP in the world frame
  /// @param orientation the 3d orientation from which the yaw will be extracted. This orientation is from local to
  /// global. i.e. bias_global == orientation * bias*local
  ///
  inline void setInputs(const Vector2 & zmpErr, const Vector2 & zmp, const Matrix3 & orientation)
  {
    setInputs(zmpErr, zmp, kine::rotationMatrixToYawAxisAgnostic(orientation));
  }

  /// @brief Set the Inputs of the estimator.
  ///
  /// @param zmperr      measurement of zmp error (zmp-zmp^{ref})
  /// @param zmp mesaurement of the ZMP in the world frame
  /// @param yaw is the yaw angle to be used. This orientation is from local to global. i.e. bias_global == R *
  /// bias*local
  inline void setInputs(const Vector2 & zmpErr, const Vector2 & zmp, double yaw)
  {
    setInputs(zmpErr, zmp, Rotation2D(yaw).toRotationMatrix());
  }

  /// @brief Set the Inputs of the estimator.
  ///
  /// @param zmperr  measurement of the DCM in the world frame
  /// @param zmp  mesaurement of the ZMP in the world frame
  /// @param R    the 2x2 Matrix'representing the yaw angle i.e. bias_global == R * bias*local
  void setInputs(const Vector2 & zmpErr, const Vector2 & zmp, const Matrix2 & R = Matrix2::Identity());

  /// @brief Runs the estimation. Needs to be called every timestep
  ///
  /// @return Vector2
  void update();

  /// @brief Get the Unbiased DCM filtered by the estimator
  ///
  /// @detailt This is the recommended output to take
  /// @return double
  inline Matrix2 getEstimatedGain() const
  {
    /// the meaning is previousOrientation_ * getLocalGain() * previousOrientation_.transpose();
    return (previousOrientation_ * getEstimatedLocalGain().triangularView<Eigen::Upper>()
            * previousOrientation_.transpose())
        .selfadjointView<Eigen::Upper>();
  }

  /// @brief Get the estimated Bias expressed in the local frame of the robot
  ///
  /// @return double
  Matrix2 getEstimatedLocalGain() const;

  /// @brief Get the Kalman Filter object
  /// This can be used to run specific Advanced Kalman filter related funcions
  /// @return LinearKalmanFilter&
  inline LinearKalmanFilter & getFilter()
  {
    return filter_;
  }

  /// @copydoc getFilter()
  /// const version
  inline const LinearKalmanFilter & getFilter() const
  {
    return filter_;
  }

protected:
  typedef Eigen::Matrix<double, 4, 2> Matrix42;
  typedef Eigen::Matrix<double, 2, 4> Matrix24;

  /// @brief set Matrices: A, B, Q
  void updateMatricesABQ_();

  double omega0_;
  double dt_;
  double biasDriftStd_;
  double zmpErrorStd_;

  Vector2 previousZmp_;

  Vector2 biasLimit_;

  LinearKalmanFilter filter_;
  Matrix4 A_;
  Matrix42 B_;
  /// this needs to be transposed
  Matrix24 C_;
  /// measurement noise
  Matrix2 R_;

  /// process noise
  Matrix4 Q_;

  Matrix2 previousOrientation_;

  /// @brief builds a diagonal out of the square valued of the Vec2
  inline static Matrix2 Vec2ToSqDiag(const Vector2 & v)
  {
    return Vector2(v.array().square()).asDiagonal();
  }

  /// @brief builds a constant 2x2 diagonal from a double
  inline static Matrix2 dblToDiag(const double & d)
  {
    return Vector2::Constant(d).asDiagonal();
  }

  /// @brief builds a constant 2x2 diagonal from a square of a double
  inline static Matrix2 dblToSqDiag(const double & d)
  {
    return dblToDiag(d * d);
  }

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

} // namespace stateObservation

#endif /// ZMPTRACKINGGAINESTIMATOR_HPP
