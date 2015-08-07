#ifndef KALMANFILTERVICTORIAPARK_HPP
#define KALMANFILTERVICTORIAPARK_HPP

#include "KalmanFilter.hpp"
#include "MeasurementModel_VictoriaPark.hpp"

namespace rfs{

  /**
   * \class KalmanFilter_VictoriaPark
   * \brief A Kalman filter for updating a 2d landmark position from
   * a single range-bearing measurements. This is derived from the base
   * Kalman Filter to handle innovation involving a rotation (bearing).
   */
  class KalmanFilter_VictoriaPark :
    public KalmanFilter <StaticProcessModel<Landmark3d>, MeasurementModel_VictoriaPark>{

  public:

    typedef MeasurementModel_VictoriaPark::TMeasurement::Vec Vec;

    /** \brief Configuration for this Kalman Filter */
    struct Config{
      /** If positive, the innovation threshold above which an update is not processed for stability reasons. */
      double rangeInnovationThreshold_;
      double bearingInnovationThreshold_;
    }config;

    /**
     * \brief Default constructor.
     */
    KalmanFilter_VictoriaPark(){
      config.rangeInnovationThreshold_ = -1;
      config.bearingInnovationThreshold_ = -1;
    };

    /**
     * \brief Constructor.
     * \param[in] pMeasurementModel Pointer to the measurement model
     * \param[in] pProcessModel Pointer to the process model
     */
    KalmanFilter_VictoriaPark(StaticProcessModel<Landmark3d> *pProcessModel,
			      MeasurementModel_VictoriaPark *pMeasurementModel):
      KalmanFilter<StaticProcessModel<Landmark3d>, MeasurementModel_VictoriaPark>
      (pProcessModel, pMeasurementModel){
      config.rangeInnovationThreshold_ = -1;
      config.bearingInnovationThreshold_ = -1;
    }

    /**
     * \brief Function to calculate the innovation
     * \param[in] z_exp expected measurement predicted using the measurement model
     * \param[in] z_act actual measurement
     * \param[out] z_innov innovation
     */
    bool calculateInnovation(Vec &z_exp, Vec &z_act, Vec &z_innov){

      z_innov = z_act - z_exp;
    
      while(z_innov(1)>PI){
	z_innov(1)-=2*PI;
      }
      while(z_innov(1)<-PI){
	z_innov(1)+=2*PI;
      }
      if(config.rangeInnovationThreshold_ > 0 && fabs(z_innov(0)) > config.rangeInnovationThreshold_){
	return false;
      }
      if(config.bearingInnovationThreshold_ > 0 && fabs(z_innov(1)) > config.bearingInnovationThreshold_){
	return false;
      }

      return true;
    }

  };

}


#endif
