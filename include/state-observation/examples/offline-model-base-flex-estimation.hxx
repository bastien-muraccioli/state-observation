const double acc_cov_const=5e-3;
const double gyr_cov_const=5e-6;
const double force_sensor_const=1e-10;
const double torque_sensor_const=1e-30;
const double state_fc_const=5e-4;
const double pos_state_cov_const=0;
const double vel_state_cov_const=1e-10;
const double ori_state_cov_const=0;
const double angv_state_cov_const=1e-8;

typedef flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem::state state;

stateObservation::IndexedMatrixArray offlineModelBaseFlexEstimation(
  const stateObservation::IndexedMatrixArray & y,
  const stateObservation::IndexedMatrixArray & u,
  const Matrix & xh0,
  const stateObservation::IndexedMatrixArray numberOfContacts,
  double dt,
  double mass,
  const Matrix3 & kfe,
  const Matrix3 & kfv,
  const Matrix3 & kte,
  const Matrix3 & ktv,
  IndexedMatrixArray * ino,
  IndexedMatrixArray * premea,
  IndexedMatrixArray * simumea,
  int verbose)
{


  flexibilityEstimation::ModelBaseEKFFlexEstimatorIMU estimator;

  estimator.setContactModel(flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem::contactModel::elasticContact);
  estimator.setSamplingPeriod(dt);
  estimator.setRobotMass(mass);

  Matrix R,Q,P;

  int measurementSize=6;
  int stateSize=flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem::state::size;

  double acceleroCovariance=acc_cov_const;
  double gyroCovariance=gyr_cov_const;
  double posStateCov=pos_state_cov_const;
  double oristateCov=ori_state_cov_const;
  double velStateCov=vel_state_cov_const;
  double angvStateCov=angv_state_cov_const;
  double stateForceCov=state_fc_const;

  R.noalias()=Matrix::Identity(measurementSize,measurementSize)*acceleroCovariance;
  R(3,3)=R(4,4)=R(5,5)=gyroCovariance;
  Q.noalias()=Matrix::Identity(stateSize,stateSize)*posStateCov;

  Q.diagonal().segment<3>(state::ori).setConstant(oristateCov);

  Q.diagonal().segment<3>(state::linVel).setConstant(velStateCov);
  Q.diagonal().segment<3>(state::angVel).setConstant(angvStateCov);
  Q.diagonal().segment<12>(state::fc).setConstant(stateForceCov);

  estimator.setProcessNoiseCovariance(Q);
  estimator.setMeasurementNoiseCovariance(R);
  Matrix forcevariance= Matrix::Identity(6,6)*force_sensor_const;
  forcevariance.block(3,3,3,3)=Matrix::Identity(3,3)*torque_sensor_const;
  /// consider the vertical torque as reliable as a force
  forcevariance(2,2)=torque_sensor_const;
  /// consider the vertical force as reliable as torque
  forcevariance(6,6)=force_sensor_const;
  estimator.setForceVariance(forcevariance);

  ///initialize flexibility
  estimator.setFlexibilityCovariance(Q);
  estimator.setFlexibilityGuess(xh0);

  estimator.setWithUnmodeledForces(true);


  bool withForce=true;
  estimator.setWithForcesMeasurements(withForce);


  if (kfe!=Matrix3::Zero())
  {
    estimator.setKfe(kfe);
    estimator.setKfv(kfv);
    estimator.setKte(kte);
    estimator.setKtv(ktv);
  }

  ///the array of the state estimations over time
  stateObservation::IndexedMatrixArray xh;
  if (y.getFirstIndex()>0)
    xh.setValue(xh0,y.getFirstIndex()-1);

  ///the reconstruction of the state
  for (unsigned i=y.getFirstIndex(); i<y.getNextIndex(); ++i)
  {

    estimator.setContactsNumber(numberOfContacts[i](0));

    if (verbose>0)
    {
      if (i%1000==0 || verbose >1)
      {
        std::cout << "iteration: " << i  << std::endl ;

        if (verbose >2)
        {
          std::cout << "number of contacts: "<< numberOfContacts[i]<< std::endl;
          std::cout << "size of the measurement: "<< y[i].size()
                    << ", supposed to be " << estimator.getMeasurementSize() << std::endl;
          std::cout << "size of the input: "<< u[i].size()
                    << ", supposed to be " << estimator.getInputSize() << std::endl;

          if (verbose > 3)
          {

            std::cout << "numberOfContacts: " << numberOfContacts[i].transpose() << std::endl;
            std::cout << "Measurement: " << y[i].transpose() << std::endl;
            std::cout << "Input: " << u[i].transpose() << std::endl;

          }
        }
      }


    }

    ///introduction of the measurement
    estimator.setMeasurement(Vector(y[i]).head(estimator.getMeasurementSize()));

    estimator.setMeasurementInput(u[i]);

    ///get the estimation and give it to the array
    Vector xhk=estimator.getFlexibilityVector();

    xh.pushBack(xhk);

    if (ino != 0)
    {

      ino->setValue(estimator.getInovation(),i);
    }

    if (premea != 0)
    {
      premea->setValue(estimator.getLastPredictedMeasurement(),i);
    }

    if (simumea != 0)
    {
      simumea->setValue(estimator.getSimulatedMeasurement(),i);
    }

    if (verbose >1)
    {
      std::cout << "Success";
      if (verbose >2)
      {
        std::cout << ", rebuilt state "<<std::endl;
        std::cout << xhk.transpose();
      }

      std::cout<<std::endl;

    }

  }

  return xh;
}



