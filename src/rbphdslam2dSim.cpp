/*
 * Software License Agreement (New BSD License)
 *
 * Copyright (c) 2013, Keith Leung, Felipe Inostroza
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Advanced Mining Technology Center (AMTC), the
 *       Universidad de Chile, nor the names of its contributors may be 
 *       used to endorse or promote products derived from this software without 
 *       specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE AMTC, UNIVERSIDAD DE CHILE, OR THE COPYRIGHT 
 * HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE 
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
 * THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <libconfig.h++>
#include "ProcessModel_Odometry2D.hpp"
#include "RBPHDFilter.hpp"
#include "KalmanFilter_RngBrg.hpp"
#include <stdio.h>
#include <string>

using namespace rfs;

/**
 * \class Simulator2d
 * \brief A 2d SLAM Simulator using the RB-PHD Filter
 * \author Keith Leung
 */
class Simulator2d{

public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  Simulator2d(){
    pFilter_ = NULL;
  }
  
  ~Simulator2d(){
    
    if(pFilter_ != NULL){
      delete pFilter_;
    }

  }

  /** Read the simulator configuration file */
  bool readConfigFile(const char* fileName){
    
    cfgFileName_ = fileName;

    try{
      cfg_.readFile( fileName );
    }catch( libconfig::FileIOException &ex){
      printf("\nCannot read file: %s\n\n", fileName);
      return false;
    }catch( libconfig::ParseException &ex){
      const char* error = ex.getError();
      int line = ex.getLine();
      printf("\n%s LINE %d\n\n", error, line);
      return false;
    }
    kMax_ = cfg_.lookup("timesteps");
    dT_ = cfg_.lookup("sec_per_timestep");
    dTimeStamp_ = TimeStamp(dT_);
    
    nSegments_ = cfg_.lookup("Trajectory.nSegments");
    max_dx_ = cfg_.lookup("Trajectory.max_dx_per_sec");
    max_dy_ = cfg_.lookup("Trajectory.max_dy_per_sec");
    max_dz_ = cfg_.lookup("Trajectory.max_dz_per_sec");
    min_dx_ = cfg_.lookup("Trajectory.min_dx_per_sec");
    vardx_ = cfg_.lookup("Trajectory.vardx");
    vardy_ = cfg_.lookup("Trajectory.vardy");
    vardz_ = cfg_.lookup("Trajectory.vardz");

    nLandmarks_ = cfg_.lookup("Landmarks.nLandmarks");
    varlmx_ = cfg_.lookup("Landmarks.varlmx");
    varlmy_ = cfg_.lookup("Landmarks.varlmy");

    rangeLimitMax_ = cfg_.lookup("Measurement.rangeLimitMax");
    rangeLimitMin_ = cfg_.lookup("Measurement.rangeLimitMin");
    rangeLimitBuffer_ = cfg_.lookup("Measurement.rangeLimitBuffer");
    Pd_ = cfg_.lookup("Measurement.probDetection");
    c_ = cfg_.lookup("Measurement.clutterIntensity");
    varzr_ = cfg_.lookup("Measurement.varzr");
    varzb_ = cfg_.lookup("Measurement.varzb");

    nParticles_ = cfg_.lookup("Filter.nParticles");
    pNoiseInflation_ = cfg_.lookup("Filter.processNoiseInflationFactor");
    zNoiseInflation_ = cfg_.lookup("Filter.measurementNoiseInflationFactor");
    innovationRangeThreshold_ = cfg_.lookup("Filter.innovationRangeThreshold");
    innovationBearingThreshold_ = cfg_.lookup("Filter.innovationBearingThreshold");
    birthGaussianWeight_ = cfg_.lookup("Filter.birthGaussianWeight");
    newGaussianCreateInnovMDThreshold_ = cfg_.lookup("Filter.newGaussianCreateInnovMDThreshold");
    importanceWeightingMeasurementLikelihoodMDThreshold_ = cfg_.lookup("Filter.importanceWeightingMeasurementLikelihoodMDThreshold");
    effNParticleThreshold_ = cfg_.lookup("Filter.effectiveNumberOfParticlesThreshold");
    minUpdatesBeforeResample_ = cfg_.lookup("Filter.minUpdatesBeforeResample");
    gaussianMergingThreshold_ = cfg_.lookup("Filter.gaussianMergingThreshold");
    gaussianMergingCovarianceInflationFactor_ = cfg_.lookup("Filter.gaussianMergingCovarianceInflationFactor");
    gaussianPruningThreshold_ = cfg_.lookup("Filter.gaussianPruningThreshold");
    importanceWeightingEvalPointCount_ = cfg_.lookup("Filter.importanceWeightingEvalPointCount");
    importanceWeightingEvalPointGuassianWeight_ = cfg_.lookup("Filter.importanceWeightingEvalPointGuassianWeight");
    useClusterProcess_ = cfg_.lookup("Filter.useClusterProcess");

    logToFile_ = cfg_.lookup("Computation.logToFile");
    logDirPrefix_ = cfg_.lookup("Computation.logDirPrefix");
    
    return true;   
  }

  /** Generate a random trajectory in 2d space */
  void generateTrajectory(int randSeed = 0){

    printf("Generating trajectory with random seed = %d\n", randSeed);
    srand48( randSeed);

    TimeStamp t;
    int seg = 0;
    MotionModel_Odometry2d::TState::Mat Q;
    Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
    MotionModel_Odometry2d motionModel(Q);
    MotionModel_Odometry2d::TInput input_k(0, 0, 0, 0, 0, 0, t);
    MotionModel_Odometry2d::TState pose_k(0, 0, 0, 0, 0, 0, t);
    MotionModel_Odometry2d::TState pose_km(0, 0, 0, 0, 0, 0, t);
    groundtruth_displacement_.reserve( kMax_ );
    groundtruth_pose_.reserve( kMax_ );
    groundtruth_displacement_.push_back(input_k);
    groundtruth_pose_.push_back(pose_k);

    for( int k = 1; k < kMax_; k++ ){

      t += dTimeStamp_;

      if( k <= 50 ){
	double dx = 0;
	double dy = 0;
	double dz = 0;
	input_k = MotionModel_Odometry2d::TInput(dx, dy, dz, 
						0, 0, 0, k);
      }else if( k >= kMax_ / nSegments_ * seg ){
	seg++;
	double dx = drand48() * max_dx_ * dT_;
	while( dx < min_dx_ * dT_ ){
	  dx = drand48() * max_dx_ * dT_;
	}
	double dy = (drand48() * max_dy_ * 2 - max_dy_) * dT_;
	double dz = (drand48() * max_dz_ * 2 - max_dz_) * dT_; 
	input_k = MotionModel_Odometry2d::TInput(dx, dy, dz, 
						Q(0,0), Q(1,1), Q(2,2), k);  
      }

      groundtruth_displacement_.push_back(input_k);
      groundtruth_displacement_.back().setTime(t);

      MotionModel_Odometry2d::TState x_k;
      motionModel.step(x_k, groundtruth_pose_[k-1], input_k, dTimeStamp_);
      groundtruth_pose_.push_back( x_k );
      groundtruth_pose_.back().setTime(t);

    }

  }
  
  /** Generate odometry measurements */
  void generateOdometry(){

    odometry_.reserve( kMax_ );
    MotionModel_Odometry2d::TInput zero;
    MotionModel_Odometry2d::TInput::Vec u0;
    u0.setZero();
    zero.set(u0, 0);
    odometry_.push_back( zero );

    MotionModel_Odometry2d::TState::Mat Q;
    Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
    MotionModel_Odometry2d motionModel(Q);
    deadReckoning_pose_.reserve( kMax_ );
    deadReckoning_pose_.push_back( groundtruth_pose_[0] );

    TimeStamp t;

    for( int k = 1; k < kMax_; k++){
      
      t += dTimeStamp_;
      double dt = dTimeStamp_.getTimeAsDouble();

      MotionModel_Odometry2d::TInput in = groundtruth_displacement_[k];
      MotionModel_Odometry2d::TState::Mat Qk = Q * dt * dt;
      in.setCov(Qk);
      MotionModel_Odometry2d::TInput out;
      in.sample(out);
      
      odometry_.push_back( out );

      MotionModel_Odometry2d::TState p;
      motionModel.step(p, deadReckoning_pose_[k-1], odometry_[k], dTimeStamp_);
      p.setTime(t);
      deadReckoning_pose_.push_back( p );
    }

  }

  /** Generate landmarks */
  void generateLandmarks(){

    MeasurementModel_RngBrg measurementModel( varzr_, varzb_);
    MeasurementModel_RngBrg::TPose pose;

    groundtruth_landmark_.reserve(nLandmarks_);

    int nLandmarksCreated = 0;
    for( int k = 1; k < kMax_; k++ ){

      if( k >= kMax_ / nLandmarks_ * nLandmarksCreated){

	MeasurementModel_RngBrg::TPose pose;
	MeasurementModel_RngBrg::TMeasurement measurementToCreateLandmark;
	MeasurementModel_RngBrg::TMeasurement::Vec z;
	double r = drand48() * rangeLimitMax_;
	double b = drand48() * 2 * PI;
	z << r, b;
	measurementToCreateLandmark.set(z);
	MeasurementModel_RngBrg::TLandmark lm;
	
	measurementModel.inverseMeasure( groundtruth_pose_[k], 
					 measurementToCreateLandmark, 
					 lm);

	groundtruth_landmark_.push_back(lm);

	nLandmarksCreated++;
	
      }

    }

  }

  /** Generate landmark measurements */
  void generateMeasurements(){

    MeasurementModel_RngBrg measurementModel( varzr_, varzb_ );
    MeasurementModel_RngBrg::TMeasurement::Mat R;
    measurementModel.getNoise(R);
    measurementModel.config.rangeLimMax_ = rangeLimitMax_;
    measurementModel.config.rangeLimMin_ = rangeLimitMin_;
    measurementModel.config.probabilityOfDetection_ = Pd_;
    measurementModel.config.uniformClutterIntensity_ = c_;
    double meanClutter = measurementModel.clutterIntensityIntegral();
    
    double expNegMeanClutter = exp( -meanClutter );
    double poissonPmf[100];
    double poissonCmf[100];
    double mean_pow_i = 1;
    double i_factorial = 1;
    poissonPmf[0] = expNegMeanClutter;
    poissonCmf[0] = poissonPmf[0];
    for( int i = 1; i < 100; i++){
      mean_pow_i *= meanClutter;
      i_factorial *= i;
      poissonPmf[i] = mean_pow_i / i_factorial * expNegMeanClutter;
      poissonCmf[i] = poissonCmf[i-1] + poissonPmf[i]; 
    }

    lmkFirstObsTime_.resize( groundtruth_landmark_.size());
    for( int m = 0; m < lmkFirstObsTime_.size(); m++ ){
      lmkFirstObsTime_[m] = -1;
    }

    TimeStamp t;

    for( int k = 1; k < kMax_; k++ ){
      
      t += dTimeStamp_;

      groundtruth_pose_[k];
      
      // Real detections
      for( int m = 0; m < groundtruth_landmark_.size(); m++){
	
	bool success;
	MeasurementModel_RngBrg::TMeasurement z_m_k;
	success = measurementModel.sample( groundtruth_pose_[k],
					   groundtruth_landmark_[m],
					   z_m_k);
	if(success){
   
	  if(z_m_k.get(0) <= rangeLimitMax_ && z_m_k.get(0) >= rangeLimitMin_ && drand48() <= Pd_){
	    z_m_k.setTime(t);
	    // z_m_k.setCov(R);
	    measurements_.push_back( z_m_k );
	  }

	  if(lmkFirstObsTime_[m] == -1){
	    lmkFirstObsTime_[m] = t.getTimeAsDouble();
	  }
	}

      }

      // False alarms
      double randomNum = drand48();
      int nClutterToGen = 0;
      while( randomNum > poissonCmf[ nClutterToGen ] ){
	nClutterToGen++;
      }
      for( int i = 0; i < nClutterToGen; i++ ){
	
	double r = drand48() * rangeLimitMax_;
	while(r < rangeLimitMin_)
	  r = drand48() * rangeLimitMax_;
	double b = drand48() * 2 * PI - PI;
	MeasurementModel_RngBrg::TMeasurement z_clutter;
	MeasurementModel_RngBrg::TMeasurement::Vec z;
	z << r, b;
	z_clutter.set(z, t);
	measurements_.push_back(z_clutter);
	
      }
      
    }
    
  }

  /** Data Logging */
  void exportSimData(){

    if(!logToFile_)
      return;

    boost::filesystem::path dir(logDirPrefix_);
    boost::filesystem::create_directory(dir);

    boost::filesystem::path cfgFilePathSrc( cfgFileName_ );
    std::string cfgFileDst( logDirPrefix_ );
    cfgFileDst += "simSettings.cfg";
    boost::filesystem::path cfgFilePathDst( cfgFileDst.data() );
    boost::filesystem::copy_file( cfgFilePathSrc, cfgFilePathDst, boost::filesystem::copy_option::overwrite_if_exists);

    TimeStamp t;

    FILE* pGTPoseFile;
    std::string filenameGTPose( logDirPrefix_ );
    filenameGTPose += "gtPose.dat";
    pGTPoseFile = fopen(filenameGTPose.data(), "w");
    MotionModel_Odometry2d::TState::Vec x;
    for(int i = 0; i < groundtruth_pose_.size(); i++){
      groundtruth_pose_[i].get(x, t);
      fprintf( pGTPoseFile, "%f   %f   %f   %f\n", t.getTimeAsDouble(), x(0), x(1), x(2));
    }
    fclose(pGTPoseFile);

    FILE* pGTLandmarkFile;
    std::string filenameGTLandmark( logDirPrefix_ );
    filenameGTLandmark += "gtLandmark.dat";
    pGTLandmarkFile = fopen(filenameGTLandmark.data(), "w");
    MeasurementModel_RngBrg::TLandmark::Vec m;
    for(int i = 0; i < groundtruth_landmark_.size(); i++){
      groundtruth_landmark_[i].get(m);
      fprintf( pGTLandmarkFile, "%f   %f   %f\n", m(0), m(1), lmkFirstObsTime_[i]);
    }
    fclose(pGTLandmarkFile);

    FILE* pOdomFile;
    std::string filenameOdom( logDirPrefix_ );
    filenameOdom += "odometry.dat";
    pOdomFile = fopen(filenameOdom.data(),"w");
    MotionModel_Odometry2d::TInput::Vec u;
    for(int i = 0; i < odometry_.size(); i++){
      odometry_[i].get(u, t);
      fprintf( pOdomFile, "%f   %f   %f   %f\n", t.getTimeAsDouble(), u(0), u(1), u(2));
    }
    fclose(pOdomFile);

    FILE* pMeasurementFile;
    std::string filenameMeasurement( logDirPrefix_ );
    filenameMeasurement += "measurement.dat";
    pMeasurementFile = fopen(filenameMeasurement.data(), "w");
    MeasurementModel_RngBrg::TMeasurement::Vec z;
    for(int i = 0; i < measurements_.size(); i++){
      measurements_[i].get(z, t);
      fprintf( pMeasurementFile, "%f   %f   %f\n", t.getTimeAsDouble(), z(0), z(1) );
    }
    fclose(pMeasurementFile);

    FILE* pDeadReckoningFile;
    std::string filenameDeadReckoning( logDirPrefix_ );
    filenameDeadReckoning += "deadReckoning.dat";
    pDeadReckoningFile = fopen(filenameDeadReckoning.data(), "w");
    MotionModel_Odometry2d::TState::Vec odo;
    for(int i = 0; i < deadReckoning_pose_.size(); i++){
      deadReckoning_pose_[i].get(odo, t);
      fprintf( pDeadReckoningFile, "%f   %f   %f   %f\n", t.getTimeAsDouble(), odo(0), odo(1), odo(2));
    }
    fclose(pDeadReckoningFile);

  }

  /** RB-PHD Filter Setup */
  void setupRBPHDFilter(){
    
    pFilter_ = new RBPHDFilter<MotionModel_Odometry2d,
			       StaticProcessModel<Landmark2d>,
			       MeasurementModel_RngBrg,
			       KalmanFilter_RngBrg>( nParticles_ );

    double dt = dTimeStamp_.getTimeAsDouble();

    // configure robot motion model (only need to set once since timesteps are constant)
    MotionModel_Odometry2d::TState::Mat Q;
    Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
    Q *= (pNoiseInflation_ * dt * dt);
    pFilter_->getProcessModel()->setNoise(Q);

    // configure landmark process model (only need to set once since timesteps are constant)
    Landmark2d::Mat Q_lm;
    Q_lm << varlmx_, 0, 0, varlmy_;
    Q_lm = Q_lm * dt * dt;
    pFilter_->getLmkProcessModel()->setNoise(Q_lm); 

    // configure measurement model
    MeasurementModel_RngBrg::TMeasurement::Mat R;
    R << varzr_, 0, 0, varzb_;
    R *= zNoiseInflation_;
    pFilter_->getMeasurementModel()->setNoise(R);
    pFilter_->getMeasurementModel()->config.probabilityOfDetection_ = Pd_;
    pFilter_->getMeasurementModel()->config.uniformClutterIntensity_ = c_;
    pFilter_->getMeasurementModel()->config.rangeLimMax_ = rangeLimitMax_;
    pFilter_->getMeasurementModel()->config.rangeLimMin_ = rangeLimitMin_;
    pFilter_->getMeasurementModel()->config.rangeLimBuffer_ = rangeLimitBuffer_;

    // configure the Kalman filter for landmark updates
    pFilter_->getKalmanFilter()->config.rangeInnovationThreshold_ = innovationRangeThreshold_;
    pFilter_->getKalmanFilter()->config.bearingInnovationThreshold_ = innovationBearingThreshold_;

    // configure the filter
    pFilter_->config.birthGaussianWeight_ = birthGaussianWeight_;
    pFilter_->setEffectiveParticleCountThreshold(effNParticleThreshold_);
    pFilter_->config.minUpdatesBeforeResample_ = minUpdatesBeforeResample_;
    pFilter_->config.newGaussianCreateInnovMDThreshold_ = newGaussianCreateInnovMDThreshold_;
    pFilter_->config.importanceWeightingMeasurementLikelihoodMDThreshold_ = importanceWeightingMeasurementLikelihoodMDThreshold_;
    pFilter_->config.importanceWeightingEvalPointCount_ = importanceWeightingEvalPointCount_;
    pFilter_->config.importanceWeightingEvalPointGuassianWeight_ = importanceWeightingEvalPointGuassianWeight_;
    pFilter_->config.gaussianMergingThreshold_ = gaussianMergingThreshold_;
    pFilter_->config.gaussianMergingCovarianceInflationFactor_ = gaussianMergingCovarianceInflationFactor_;
    pFilter_->config.gaussianPruningThreshold_ = gaussianPruningThreshold_;
    pFilter_->config.useClusterProcess_ = useClusterProcess_;
  }

  /** Run the simulator */
  void run(){
    
    printf("Running simulation\n\n");

    //////// Initialization at first timestep //////////

    FILE* pParticlePoseFile;
    if(logToFile_){
      std::string filenameParticlePoseFile( logDirPrefix_ );
      filenameParticlePoseFile += "particlePose.dat";
      pParticlePoseFile = fopen(filenameParticlePoseFile.data(), "w");
    }
    FILE* pLandmarkEstFile;
    if(logToFile_){
      std::string filenameLandmarkEstFile( logDirPrefix_ );
      filenameLandmarkEstFile += "landmarkEst.dat";
      pLandmarkEstFile = fopen(filenameLandmarkEstFile.data(), "w");
    }
    MotionModel_Odometry2d::TState x_i;
    int zIdx = 0;

    if(logToFile_){
      for(int i = 0; i < pFilter_->getParticleCount(); i++){
	pFilter_->getParticleSet()->at(i)->getPose(x_i);
	fprintf( pParticlePoseFile, "%f   %d   %f   %f   %f   1.0\n", 0.0, i, x_i.get(0), x_i.get(1), x_i.get(2));
      }   
    }

    boost::timer::auto_cpu_timer *stepTimer = NULL;
    boost::timer::auto_cpu_timer *processTimer = NULL;
    
    processTimer = new boost::timer::auto_cpu_timer(6, "Total run time: %ws\n");
    /////////// Run simulator from k = 1 to kMax_ /////////

    TimeStamp time;

    for(int k = 1; k < kMax_; k++){

      time += dTimeStamp_;
      
      if( k % 100 == 0)
	printf("k = %d\n", k);
      
      ////////// Prediction Step //////////

      // configure robot motion model ( not necessary since in simulation, timesteps are constant)
      // MotionModel_Odometry2d::TState::Mat Q;
      // Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
      // Q *= (pNoiseInflation_ * dt * dt);
      // pFilter_->getProcessModel()->setNoise(Q);

      // configure landmark process model ( not necessary since in simulation, timesteps are constant)
      // Landmark2d::Mat Q_lm;
      // Q_lm << varlmx_, 0, 0, varlmy_;
      // Q_lm = Q_lm * dt * dt;
      // pFilter_->getLmkProcessModel()->setNoise(Q_lm);

      pFilter_->predict( odometry_[k], dTimeStamp_ );
      
      if( k <= 100){
	for( int i = 0; i < nParticles_; i++)
	  pFilter_->setParticlePose(i, groundtruth_pose_[k]);
      }

      // Prepare measurement vector for update
      std::vector<MeasurementModel_RngBrg::TMeasurement> Z;
      TimeStamp kz = measurements_[ zIdx ].getTime();
      while( kz == time ){ 
	Z.push_back( measurements_[zIdx] );
	zIdx++;
	if(zIdx >= measurements_.size())
	  break;
	kz = measurements_[ zIdx ].getTime();
      }

      ////////// Update Step //////////
      pFilter_->update(Z);

      // Log particle poses
      if(logToFile_){
	for(int i = 0; i < pFilter_->getParticleCount(); i++){
	  pFilter_->getParticleSet()->at(i)->getPose(x_i);
	  double w = pFilter_->getParticleSet()->at(i)->getWeight();
	  fprintf( pParticlePoseFile, "%f   %d   %f   %f   %f   %f\n", time.getTimeAsDouble(), i, x_i.get(0), x_i.get(1), x_i.get(2), w);
	}
	fprintf( pParticlePoseFile, "\n");
      }

      // Log landmark estimates
      if(logToFile_){
	for(int i = 0; i < pFilter_->getParticleCount(); i++){
	  int mapSize = pFilter_->getGMSize(i);
	  for( int m = 0; m < mapSize; m++ ){
	    MeasurementModel_RngBrg::TLandmark::Vec u;
	    MeasurementModel_RngBrg::TLandmark::Mat S;
	    double w;
	    pFilter_->getLandmark(i, m, u, S, w);
	    
	    fprintf( pLandmarkEstFile, "%f   %d   ", time.getTimeAsDouble(), i);
	    fprintf( pLandmarkEstFile, "%f   %f      ", u(0), u(1));
	    fprintf( pLandmarkEstFile, "%f   %f   %f", S(0,0), S(0,1), S(1,1));
	    fprintf( pLandmarkEstFile, "   %f\n", w );
	  }
	}
      }

    }

    delete processTimer;

    
    printf("Elapsed Timing Information [nsec]\n");
    printf("Prediction    -- wall: %lld   cpu: %lld\n", 
	   pFilter_->getTimingInfo()->predict_wall, pFilter_->getTimingInfo()->predict_cpu);
    printf("Map Update    -- wall: %lld   cpu: %lld\n", 
	   pFilter_->getTimingInfo()->mapUpdate_wall, pFilter_->getTimingInfo()->mapUpdate_cpu);
    printf("Map Update KF -- wall: %lld   cpu: %lld\n", 
	   pFilter_->getTimingInfo()->mapUpdate_kf_wall, pFilter_->getTimingInfo()->mapUpdate_kf_cpu);
    printf("Weighting     -- wall: %lld   cpu: %lld\n", 
	   pFilter_->getTimingInfo()->particleWeighting_wall, pFilter_->getTimingInfo()->particleWeighting_cpu);
    printf("Map Merge     -- wall: %lld   cpu: %lld\n", 
	   pFilter_->getTimingInfo()->mapMerge_wall, pFilter_->getTimingInfo()->mapMerge_cpu);
    printf("Map Prune     -- wall: %lld   cpu: %lld\n", 
	   pFilter_->getTimingInfo()->mapPrune_wall, pFilter_->getTimingInfo()->mapPrune_cpu);
    printf("Resampling    -- wall: %lld   cpu: %lld\n", 
	   pFilter_->getTimingInfo()->particleResample_wall, pFilter_->getTimingInfo()->particleResample_cpu);
    printf("Total         -- wall: %lld   cpu: %lld\n",
	   pFilter_->getTimingInfo()->predict_wall +
	   pFilter_->getTimingInfo()->mapUpdate_wall +
	   pFilter_->getTimingInfo()->particleWeighting_wall +
	   pFilter_->getTimingInfo()->mapMerge_wall +
	   pFilter_->getTimingInfo()->mapPrune_wall +
	   pFilter_->getTimingInfo()->particleResample_wall,
	   pFilter_->getTimingInfo()->predict_cpu +
	   pFilter_->getTimingInfo()->mapUpdate_cpu +
	   pFilter_->getTimingInfo()->particleWeighting_cpu +
	   pFilter_->getTimingInfo()->mapMerge_cpu +
	   pFilter_->getTimingInfo()->mapPrune_cpu + 
	   pFilter_->getTimingInfo()->particleResample_cpu);
    printf("\n");

    if(logToFile_){
      fclose(pParticlePoseFile);
      fclose(pLandmarkEstFile);
    }
  }

private:

  libconfig::Config cfg_;
  const char* cfgFileName_;

  int kMax_; /**< number of timesteps */
  double dT_; /**< duration of timestep in seconds */
  TimeStamp dTimeStamp_; /**< duration of timestep in timestamp */

  // Trajectory
  int nSegments_;
  double max_dx_;
  double max_dy_;
  double max_dz_;
  double min_dx_;
  double vardx_;
  double vardy_;
  double vardz_;
  std::vector<MotionModel_Odometry2d::TInput> groundtruth_displacement_;
  std::vector<MotionModel_Odometry2d::TState> groundtruth_pose_;
  std::vector<MotionModel_Odometry2d::TInput> odometry_;
  std::vector<MotionModel_Odometry2d::TState> deadReckoning_pose_;

  // Landmarks 
  int nLandmarks_;
  std::vector<MeasurementModel_RngBrg::TLandmark> groundtruth_landmark_;
  double varlmx_;
  double varlmy_;
  std::vector<double> lmkFirstObsTime_;

  // Range-Bearing Measurements
  double rangeLimitMax_;
  double rangeLimitMin_;
  double rangeLimitBuffer_;
  double Pd_;
  double c_;
  double varzr_;
  double varzb_;
  std::vector<MeasurementModel_RngBrg::TMeasurement> measurements_;

  // Filters
  KalmanFilter_RngBrg kf_;
  RBPHDFilter<MotionModel_Odometry2d, 
	      StaticProcessModel<Landmark2d>,
	      MeasurementModel_RngBrg,  
	      KalmanFilter_RngBrg> *pFilter_; 
  int nParticles_;
  double pNoiseInflation_;
  double zNoiseInflation_;
  double innovationRangeThreshold_;
  double innovationBearingThreshold_;
  double birthGaussianWeight_;
  double newGaussianCreateInnovMDThreshold_;
  double importanceWeightingMeasurementLikelihoodMDThreshold_;
  double importanceWeightingEvalPointGuassianWeight_;
  double effNParticleThreshold_;
  int minUpdatesBeforeResample_;
  double gaussianMergingThreshold_;
  double gaussianMergingCovarianceInflationFactor_;
  double gaussianPruningThreshold_;
  int importanceWeightingEvalPointCount_;
  bool useClusterProcess_;

  bool logToFile_;
  const char *logDirPrefix_;
};











int main(int argc, char* argv[]){

  int initRandSeed = 0;
  const char* logFileName = "cfg/rbphdslam2dSim.cfg";
  if( argc >= 2 ){
    initRandSeed = boost::lexical_cast<int>(argv[1]);
  }
  if( argc >= 3 ){
    logFileName = argv[2];
  }

  Simulator2d sim;
  if( !sim.readConfigFile( logFileName ) ){
    return -1;
  }
  sim.generateTrajectory( initRandSeed );
  sim.generateLandmarks();

  sim.generateOdometry();
  sim.generateMeasurements();
  sim.exportSimData();
 
  sim.setupRBPHDFilter();

  srand48( time(NULL) );

  boost::timer::auto_cpu_timer *timer = new boost::timer::auto_cpu_timer(6, "Simulation run time: %ws\n");

  sim.run(); 

  delete timer;

  return 0;

}
