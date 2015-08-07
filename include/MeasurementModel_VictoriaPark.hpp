#ifndef MEASUREMENTMODEL_VICTORIAPARK_HPP
#define MEASUREMENTMODEL_VICTORIAPARK_HPP

#include <boost/shared_ptr.hpp>
#include <vector>

#include "MeasurementModel.hpp"
#include "MeasurementModel_RngBrg.hpp"

namespace rfs{
  
  /**
   * \class MeasurementModel_VictoriaPark
   * \brief A 2d range-bearing measurement model for 2d Circle Features extracted from lidar scans.
   * \author Felipe Inostroza
   */
  class MeasurementModel_VictoriaPark : public MeasurementModel <Pose2d, Landmark3d, Measurement3d>{

  public:

    /** Default constructor */
    MeasurementModel_VictoriaPark();

    /**
     * Constructor that sets the uncertainty (covariance) of the measurement model, \f$\mathbf{R}\f$
     * \param covZ measurement covariance, \f$\mathbf{R}\f$
     * \param Slb laser bearing variance used in calculating the diameter variance
     */
    MeasurementModel_VictoriaPark(Eigen::Matrix3d &covZ, double Slb);


    /**
     * Constructor that sets the uncertainty (covariance) of the measurement model, \f$\mathbf{R}\f$
     * \param covP measurement covariance for the range bearing model, \f$\mathbf{R_p}\f$
     * \param covR measurement covariance for the diameter, \f$\mathbf{R}_R\f$
     * \param Slb laser bearing variance used in calculating the diameter variance
     */
    MeasurementModel_VictoriaPark(Eigen::Matrix2d &covP, double covR, double Slb);

    /**
     * Constructor that sets the uncertainty (covariance) of the measurement model,
     * \f[\mathbf{R} = \begin{bmatrix} \sigma_r^2 & 0 & 0 \\ 0 & \sigma_b^2 & 0 \\ 0 & 0 & \sigma_R^2 \end{bmatrix}\f]
     * radius, range and bearing are assumed to be uncorrelated
     * \param Sr Range variance \f$\sigma_r^2\f$
     * \param Sb Bearing variance \f$\sigma_b^2\f$
     * \param Sd flat diameter variance \f$\sigma_R^2\f$
     * \param Slb laser bearing variance used in calculating the diameter variance
     */
    MeasurementModel_VictoriaPark(double Sr, double Sb, double Sd, double Slb);

    /** Default destructor */
    ~MeasurementModel_VictoriaPark();

    /**
     * Set the zero-mean-white-Gaussian additive noise covariance matrix, \f$\mathbf{R}\f$
     * \param[in] R covariance matrix
     * \param Slb laser bearing variance used in calculating the diameter variance
     */
    void setNoise( Measurement3d::Mat &R , double Slb);

    /**
     * Obtain a measurement from a given robot pose and landmark position
     * \f[ \mathbf{z} = \mathbf{h}(\mathbf{x}, \mathbf{m} ) + \mathbf{e}, \mathbf{e} \sim (\mathbf{0}, \mathbf{R}) \f]
     * where \f$\mathbf{z}\f$ is a measurement, \f$\mathbf{x}\f$ is the robot pose, \f$\mathbf{m}\f$ is a landmark position, \f$\mathbf{e}\f$ is the zero-mean Gaussian noise.
     * \param[in] pose \f$\mathbf{x}\f$, robot pose from which the measurement is made
     * \param[in] landmark \f$\mathbf{m}\f$, the measured landmark
     * \param[out] measurement \f$\mathbf{x}\f$, the measurement
     * \param[out] jacobian_wrt_lmk  if not NULL, the pointed-to matrix is overwritten with the Jacobian w.r.t. the landmark position
     * \param[out] jacobian_wrt_pose if not NULL, the pointed-to matrix is overwritten with the Jacobian w.r.t. the robot pose
     * by the Jacobian of the measurement model, \f$\mathbf{H}\f$, evaluated at \f$\mathbf{x}\f$ and \f$\mathbf{m}\f$
     * \return true if a valid measurement is produced
     */
    bool measure( const Pose2d &pose, const Landmark3d &landmark,
		  Measurement3d &measurement, Eigen::Matrix3d *jacobian_wrt_lmk = NULL, Eigen::Matrix3d *jacobian_wrt_pose = NULL);

    /**
     * \f[ \mathbf{m} = \mathbf{h}^{-1}(\mathbf{x}, \mathbf{z} )\f]
     * where \f$\mathbf{z}\f$ is a measurement, \f$\mathbf{x}\f$ is the robot pose, \f$\mathbf{m}\f$ is a landmark position
     * \param[in] pose \f$\mathbf{x}\f$, robot pose (the uncertainty is not used here because the
     * RBPHDFilter represents robot pose estimates with particles)
     * \param[in] measurement  \f$\mathbf{z}\f$ measurement, for which the uncertainty is \f$\mathbf{R}\f$
     * \param[out] landmark  \f$\mathbf{m}\f$, predicted landmark position with uncertainty
     */
    void inverseMeasure(const Pose2d &pose, const Measurement3d &measurement, Landmark3d &landmark);


    /**
     * \brief This function calculates the probability of detection (using probabilityODetection2) in several places to determine 
     * whether it is both zero and nonzero in the landmark's neighborhood, if it is isCloseToSensingLimit is set to 1.
     * \param[in] pose robot pose
     * \param[in] landmark landmark position
     * \param[out] isCloseToSensingLimit true if landmark is close to the sensing limit
     * \return probability of detection
     */
    double probabilityOfDetection( const Pose2d &pose,
				   const Landmark3d &landmark,
				   bool &isCloseToSensingLimit);

    /**
     * \brief Function to calculate the probability of detection at a known exact landmark position.
     * \param[in] pose robot pose
     * \param[in] landmark landmark position
     * \param[out] isCloseToSensingLimit true if landmark is close to the sensing limit
     * \return probability of detection
     */
    double probabilityOfDetection2( const Pose2d &pose,
				    const Landmark3d &landmark,
				    bool &isCloseToSensingLimit);

    /**
     * Determine the clutter intensity in measurement space.
     * Uniform clutter intensity is assumed
     * \param[in] z measurement point at which clutter intensity will be determined
     * \param[in] nZ the cardinality of Z, of which z is a member.
     * \return clutter intensity
     */
    double clutterIntensity( Measurement3d &z, int nZ );

    /**
     * Determine the clutter intensity integral in measurement space.
     * This is calculated based on the probablity of false alarm,
     * defined as p( NULL | measurement exists)
     * \param[in] nZ the cardinality of Z
     * \return clutter intensity
     */
    double clutterIntensityIntegral( int nZ = 0);


    /**
     * Set the raw lidar scan, which is used for determining the probability of detection of objects
     * \param[in] laserscan The lidar scan as a std vector
     */ 
    void setLaserScan(const std::vector<double> &laserscan);

    /** \brief Configuration for this measurement model */
    struct Config{
      std::vector<double> probabilityOfDetection_; /**< Array containing probabilities of detection for a circle given the amount of points that should show up on the scan */
      double expectedClutterNumber_; /**< Expected number of clutter measurements, the intensity is asummed constant over the field of view of the robot */
      double rangeLimMax_; /**< sensing range limit, beyond which \f$ P_D = 0 \f$*/
      double rangeLimMin_; /**< sensing range limit, below which \f$ P_D = 0 \f$*/
      double bearingLimitMin_; /**< sensing angle limit, below which \f$ P_D = 0 \f$*/
      double bearingLimitMax_; /**< sensing angle limit, beyond which \f$ P_D = 0 \f$*/
      double bufferZonePd_; /**< Pd below which a feature is considered in the buffer zone */

    }config;

  private:

    double Slb_;
    MeasurementModel_RngBrg rangeBearingModel;
    std::vector<double> laserscan_;
    double clutterIntensity_;

  };

}








#endif

