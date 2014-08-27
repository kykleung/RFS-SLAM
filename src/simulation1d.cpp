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

// 1D Simulator for testing the RBPHDFilter
// Keith Leung 2013
// Requires libconfig to be installed


#include <boost/lexical_cast.hpp>
#include <libconfig.h++>
#include "RBPHDFilter.hpp"
#include <stdio.h>

class Simulator1d{

public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  Simulator1d(){

  }

  ~Simulator1d(){

  }

  bool readConfigFile(){

    libconfig::Config cfg;
    const char* fileName= "cfg/simulator1d.cfg";
    try{
      cfg.readFile( fileName );
    }catch( libconfig::FileIOException &ex){
      printf("\nCannot read file: %s\n\n", fileName);
      return false;
    }catch( libconfig::ParseException &ex){
      const char* error = ex.getError();
      int line = ex.getLine();
      printf("\n%s LINE %d\n\n", error, line);
      return false;
    }
    kMax_ = cfg.lookup("timesteps");
    dT_ = cfg.lookup("sec_per_timestep");

    nSegments_ = cfg.lookup("Trajectory.nSegments");
    max_dx_ = cfg.lookup("Trajectory.max_dx_per_sec");
    min_dx_ = cfg.lookup("Trajectory.min_dx_per_sec");
    vardx_ = cfg.lookup("Trajectory.vardx");

    nLandmarks_ = cfg.lookup("Landmarks.nLandmarks");
    varlmx_ = cfg.lookup("Landmarks.varlmx");

    rangeLimitMax_ = cfg.lookup("Measurement.rangeLimitMax");
    rangeLimitMin_ = cfg.lookup("Measurement.rangeLimitMin");
    rangeLimitBuffer_ = cfg.lookup("Measurement.rangeLimitBuffer");
    Pd_ = cfg.lookup("Measurement.probDetection");
    c_ = cfg.lookup("Measurement.clutterIntensity");
    varzr_ = cfg.lookup("Measurement.varzr");

    nParticles_ = cfg.lookup("Filter.nParticles");
    zNoiseInflation_ = cfg.lookup("Filter.measurementNoiseInflationFactor");
    birthGaussianWeight_ = cfg.lookup("Filter.birthGaussianWeight");
    newGaussianCreateInnovMDThreshold_ = cfg.lookup("Filter.newGaussianCreateInnovMDThreshold");
    importanceWeightingMeasurementLikelihoodMDThreshold_ = cfg.lookup("Filter.importanceWeightingMeasurementLikelihoodMDThreshold");
    effNParticleThreshold_ = cfg.lookup("Filter.effectiveNumberOfParticlesThreshold");
    minInterSampleTimesteps_ = cfg.lookup("Filter.minInterSampleTimesteps");
    gaussianMergingThreshold_ = cfg.lookup("Filter.gaussianMergingThreshold");
    gaussianMergingCovarianceInflationFactor_ = cfg.lookup("Filter.gaussianMergingCovarianceInflationFactor");
    gaussianPruningThreshold_ = cfg.lookup("Filter.gaussianPruningThreshold");
    importanceWeightingEvalPointCount_ = cfg.lookup("Filter.importanceWeightingEvalPointCount");
    reportTimingInfo_ = cfg.lookup("Filter.reportTimingInfo");

    nThreadsPropagationStep_ = cfg.lookup("Computation.nThreadsPropagationStep");
    logToFile_ = cfg.lookup("Computation.logToFile");

    return true;

  }


  /** Generate a random trajectory */
  void generateTrajectory(int randSeed = 0){

    printf("Generating trajectory with random seed = %d\n", randSeed);
    srand48( randSeed);

    OdometryMotionModel1d::TState::Mat Q;
    Q << vardx_;
    Q = dT_ * Q * dT_;
    // std::cout << "\n\n" << Q << "\n\n";
    OdometryMotionModel1d motionModel(Q);

    OdometryMotionModel1d::TInput::Vec u_k;
    OdometryMotionModel1d::TInput::Mat Su_k;
    OdometryMotionModel1d::TInput input_k;
    u_k.setZero();
    Su_k.setZero();
    input_k.set(u_k, Su_k, 0);
    OdometryMotionModel1d::TState pose_k(0, 0, 0);
    OdometryMotionModel1d::TState pose_km(0, 0, 0);

    groundtruth_displacement_.reserve( kMax_ );
    groundtruth_pose_.reserve( kMax_ );
    groundtruth_displacement_.push_back(input_k);
    groundtruth_pose_.push_back(pose_k);

    int seg = 0;
    for( int k = 1; k < kMax_; k++ ){

      if( k >= kMax_ / nSegments_ * seg ){
	seg++;
	double dx = drand48() * max_dx_ * dT_;
	while( dx < min_dx_ * dT_ ){
	  dx = drand48() * max_dx_ * dT_;
	}
	if( drand48() < 0.5 )
	  dx *= -1;
	u_k << dx;
	input_k.set(u_k, Su_k, k);
      }

      groundtruth_displacement_.push_back(input_k);

      OdometryMotionModel1d::TState x_k;
      motionModel.step(x_k, groundtruth_pose_[k-1], input_k);
      x_k.setTime(k);
      groundtruth_pose_.push_back( x_k );

    }

  }

  /* Generate odometry measurements */ 
  void generateOdometry(){

    odometry_.reserve( kMax_ );
    OdometryMotionModel1d::TInput zero;
    OdometryMotionModel1d::TInput::Vec u0;
    u0.setZero();
    zero.set(u0, 0);
    odometry_.push_back( zero );

    OdometryMotionModel1d::TState::Mat Q;
    Q << vardx_;
    OdometryMotionModel1d motionModel(Q);
    deadReckoning_pose_.reserve( kMax_ );
    deadReckoning_pose_.push_back( groundtruth_pose_[0] );

    for( int k = 1; k < kMax_; k++){
      
      OdometryMotionModel1d::TInput in = groundtruth_displacement_[k];
      in.setCov(Q);
      OdometryMotionModel1d::TInput out;
      RandomVecMathTools<OdometryMotionModel1d::TInput>::sample(in, out);
      odometry_.push_back( out );

      OdometryMotionModel1d::TState p;
      motionModel.step(p, deadReckoning_pose_[k-1], odometry_[k]);
      p.setTime(k);
      deadReckoning_pose_.push_back( p );
    }

  }

  /* Generate landmarks */
  void generateLandmarks(){

    MeasurementModel1d measurementModel( varzr_ );
    groundtruth_landmark_.reserve(nLandmarks_);

    int nLandmarksCreated = 0;
    for( int k = 1; k < kMax_; k++ ){

      if( k >= kMax_ / nLandmarks_ * nLandmarksCreated){

	MeasurementModel1d::TMeasurement measurementToCreateLandmark;
	MeasurementModel1d::TMeasurement::Vec z;
	double r = drand48() * rangeLimitMax_;
	z << r;
	measurementToCreateLandmark.set(z);
	MeasurementModel1d::TLandmark lm;
	
	measurementModel.inverseMeasure( groundtruth_pose_[k], 
					 measurementToCreateLandmark, 
					 lm);

	groundtruth_landmark_.push_back(lm);

	nLandmarksCreated++;
	
      }

    }

  }

  /* Generate measurements */
  void generateMeasurements(){

    MeasurementModel1d measurementModel( varzr_ );
    MeasurementModel1d::TMeasurement::Mat R;
    measurementModel.getNoise(R);
    measurementModel.config.rangeLimMax_ = rangeLimitMax_;
    measurementModel.config.rangeLimMin_ = rangeLimitMin_;
    measurementModel.config.probabilityOfDetection_ = Pd_;
    measurementModel.config.uniformClutterIntensity_ = c_;
    double meanClutter = measurementModel.clutterIntensityIntegral();
    
    printf("Range lim max = %f\n", measurementModel.config.rangeLimMax_);
    printf("Range lim min = %f\n", measurementModel.config.rangeLimMin_);
    printf("Var_z = %f\n", varzr_);

    // lookup table for generating clutter
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

    for( int k = 1; k < kMax_; k++ ){
      
      // Real detections
      for( int m = 0; m < groundtruth_landmark_.size(); m++){
	
	bool success;
	MeasurementModel1d::TMeasurement z_m_k;
	MeasurementModel1d::TMeasurement::Vec zv;
	success = measurementModel.sample( groundtruth_pose_[k],
					   groundtruth_landmark_[m],
					   z_m_k);
	/*success = measurementModel.measure( groundtruth_pose_[k],
					    groundtruth_landmark_[m],
					    z_m_k);*/
	z_m_k.get(zv);
	if(success){
	  if(drand48() <= Pd_){
	    /*printf("Measurement[%d] = [%f %f]\n", int(measurements_.size()),
	      z_m_k.get(0), z_m_k.get(1)); */
	    z_m_k.setTime(k);
	    z_m_k.setCov(R);
	    measurements_.push_back( z_m_k );
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
	
	MeasurementModel1d::TMeasurement z_clutter;
	MeasurementModel1d::TMeasurement::Vec z;
	double r = drand48() * rangeLimitMax_;
	z << r;
	z_clutter.set(z, k);
	measurements_.push_back(z_clutter);
	
      }
      
    }
    
  }

  /* Export simulation data for plotting */
  void exportSimData(){

    double t;

    FILE* pGTPoseFile;
    pGTPoseFile = fopen("data/gtPose.dat", "w");
    OdometryMotionModel1d::TState::Vec x;
    for(int i = 0; i < groundtruth_pose_.size(); i++){
      groundtruth_pose_[i].get(x, t);
      fprintf( pGTPoseFile, "%f   %f\n", t, x(0));
    }
    fclose(pGTPoseFile);

    FILE* pGTLandmarkFile;
    pGTLandmarkFile = fopen("data/gtLandmark.dat", "w");
    MeasurementModel1d::TLandmark::Vec m;
    for(int i = 0; i < groundtruth_landmark_.size(); i++){
      groundtruth_landmark_[i].get(m);
      fprintf( pGTLandmarkFile, "%f\n", m(0));
    }
    fclose(pGTLandmarkFile);

    FILE* pOdomFile;
    pOdomFile = fopen("data/odometry.dat","w");
    OdometryMotionModel1d::TInput::Vec u;
    for(int i = 0; i < odometry_.size(); i++){
      odometry_[i].get(u, t);
      fprintf( pOdomFile, "%f   %f\n", t, u(0));
    }
    fclose(pOdomFile);

    FILE* pMeasurementFile;
    pMeasurementFile = fopen("data/measurement.dat", "w");
    MeasurementModel1d::TMeasurement::Vec z;
    for(int i = 0; i < measurements_.size(); i++){
      measurements_[i].get(z, t);
      fprintf( pMeasurementFile, "%f   %f\n", t, z(0));
    }
    fclose(pMeasurementFile);

    FILE* pDeadReckoningFile;
    pDeadReckoningFile = fopen("data/deadReckoning.dat", "w");
    OdometryMotionModel1d::TState::Vec odo;
    for(int i = 0; i < deadReckoning_pose_.size(); i++){
      deadReckoning_pose_[i].get(odo, t);
      fprintf( pDeadReckoningFile, "%f   %f\n", t, odo(0));
    }
    fclose(pDeadReckoningFile);

  }

  /* Setup the RB-PHD-SLAM filter */
  void setupRBPHDFilter(){
    
    pFilter_ = new RBPHDFilter<OdometryMotionModel1d,
			       StaticProcessModel<Landmark1d>,
			       MeasurementModel1d,
			       KalmanFilter< StaticProcessModel<Landmark1d>, MeasurementModel1d > >( nParticles_ );

    // configure robot motion model
    OdometryMotionModel1d::TState::Mat Q;
    Q << vardx_;
    pFilter_->getProcessModel()->setNoise(Q);

    // configure landmark process model
    Landmark1d::Mat Q_lm;
    Q_lm << varlmx_;
    pFilter_->getLmkProcessModel()->setNoise(Q_lm);

    // configure measurement model
    MeasurementModel1d::TMeasurement::Mat R;
    R << varzr_;
    R *= zNoiseInflation_;
    pFilter_->getMeasurementModel()->setNoise(R);
    pFilter_->getMeasurementModel()->config.probabilityOfDetection_ = Pd_;
    pFilter_->getMeasurementModel()->config.uniformClutterIntensity_ = c_;
    pFilter_->getMeasurementModel()->config.rangeLimMax_ = rangeLimitMax_;
    pFilter_->getMeasurementModel()->config.rangeLimMin_ = rangeLimitMin_;
    pFilter_->getMeasurementModel()->config.rangeLimBuffer_ = rangeLimitBuffer_;

    // configure the RB-PHD-SLAM filter
    pFilter_->config.birthGaussianWeight_ = birthGaussianWeight_;
    pFilter_->setEffectiveParticleCountThreshold(effNParticleThreshold_);
    pFilter_->config.minInterSampleTimesteps_ = minInterSampleTimesteps_;
    pFilter_->config.newGaussianCreateInnovMDThreshold_ = newGaussianCreateInnovMDThreshold_;
    pFilter_->config.importanceWeightingMeasurementLikelihoodMDThreshold_ = importanceWeightingMeasurementLikelihoodMDThreshold_;
    pFilter_->config.importanceWeightingEvalPointCount_ = importanceWeightingEvalPointCount_;
    pFilter_->config.gaussianMergingThreshold_ = gaussianMergingThreshold_;
    pFilter_->config.gaussianMergingCovarianceInflationFactor_ = gaussianMergingCovarianceInflationFactor_;
    pFilter_->config.gaussianPruningThreshold_ = gaussianPruningThreshold_;
    pFilter_->config.reportTimingInfo_ = reportTimingInfo_;

    // set multi-threading parameter
    pFilter_->PFconfig.nThreadsPropagationStep_ = nThreadsPropagationStep_;

  }


  void run(){

    printf("Running simulation\n\n");
    
    FILE* pParticlePoseFile;
    if(logToFile_){
      pParticlePoseFile = fopen("data/particlePose.dat", "w");
      fprintf( pParticlePoseFile, "Timesteps: %d\n", kMax_);
      fprintf( pParticlePoseFile, "nParticles: %d\n\n", pFilter_->getParticleCount());
    }
    FILE* pLandmarkEstFile;
    if(logToFile_){
      pLandmarkEstFile = fopen("data/landmarkEst.dat", "w");
      fprintf( pLandmarkEstFile, "Timesteps: %d\n", kMax_);
      fprintf( pLandmarkEstFile, "nParticles: %d\n\n", pFilter_->getParticleCount());
    }
    OdometryMotionModel1d::TState x_i;
    int zIdx = 0;
    if(logToFile_){
      for(int i = 0; i < pFilter_->getParticleCount(); i++){
	pFilter_->getParticleSet()->at(i)->getPose(x_i);
	fprintf( pParticlePoseFile, "%f   1.0\n", x_i.get(0));
      }
    }

    for(int k = 1; k < kMax_; k++){

      if( k % 100 == 0)
	printf("k = %d\n", k);
      if(logToFile_){
	fprintf( pParticlePoseFile, "k = %d\n", k);
      }
      pFilter_->predict( odometry_[k], k );

      // Prepare measurement vector for update
      std::vector<MeasurementModel1d::TMeasurement> Z;
      double kz = measurements_[ zIdx ].getTime();
      while( kz == k ){
	Z.push_back( measurements_[zIdx] );
	zIdx++;
	if(zIdx >= measurements_.size())
	  break;
	kz = measurements_[ zIdx ].getTime();
      }

      pFilter_->update(Z, k);
      
      if(logToFile_){
      // Log particle pose
	for(int i = 0; i < pFilter_->getParticleCount(); i++){
	  pFilter_->getParticleSet()->at(i)->getPose(x_i);
	  double w = pFilter_->getParticleSet()->at(i)->getWeight();
	  fprintf( pParticlePoseFile, "%f   %f\n", x_i.get(0), w);
	}
	fprintf( pParticlePoseFile, "\n");
      
	// Log landmark estimates
	for(int i = 0; i < pFilter_->getParticleCount(); i++){
	  int mapSize = pFilter_->getGMSize(i);
	  fprintf( pLandmarkEstFile, "Timestep: %d\tParticle: %d\tMap Size: %d\n", k, i, mapSize);
	  for( int m = 0; m < mapSize; m++ ){
	    MeasurementModel1d::TLandmark::Vec u;
	    MeasurementModel1d::TLandmark::Mat S;
	    double w;
	    pFilter_->getLandmark(i, m, u, S, w);
	  
	    fprintf( pLandmarkEstFile, "%f      ", u(0));
	    fprintf( pLandmarkEstFile, "%f      ", S(0,0));
	    fprintf( pLandmarkEstFile, "%f\n", w );
	  }
	  fprintf( pLandmarkEstFile, "\n");
	}
	fprintf( pLandmarkEstFile, "\n");
      }
    }
    if(logToFile_){
      fclose(pParticlePoseFile);
      fclose(pLandmarkEstFile);
    }
  }


private:

  int kMax_; /**< number of timesteps */
  double dT_;

  /* Trajectory-related variables */
  int nSegments_; 
  double max_dx_;
  double min_dx_; 
  double vardx_; 
  std::vector<OdometryMotionModel1d::TInput> groundtruth_displacement_;
  std::vector<OdometryMotionModel1d::TState> groundtruth_pose_;
  std::vector<OdometryMotionModel1d::TInput> odometry_;
  std::vector<OdometryMotionModel1d::TState> deadReckoning_pose_;

  /* Landmark-related variables */
  int nLandmarks_;
  double varlmx_;
  std::vector<MeasurementModel1d::TLandmark> groundtruth_landmark_;

  /* Measurement-related variables */
  double rangeLimitMax_;
  double rangeLimitMin_;
  double rangeLimitBuffer_;
  double Pd_;
  double c_;
  double varzr_;
  std::vector<MeasurementModel1d::TMeasurement> measurements_;

  /* RB-PHD-SLAM-Filter-related variables */
  RangeBearingKalmanFilter kf_;
  RBPHDFilter<OdometryMotionModel1d,
	      StaticProcessModel<Landmark1d>,
	      MeasurementModel1d,
	      KalmanFilter< StaticProcessModel<Landmark1d>, MeasurementModel1d > > *pFilter_; 
  int nParticles_;
  double zNoiseInflation_;
  double birthGaussianWeight_;
  double newGaussianCreateInnovMDThreshold_;
  double importanceWeightingMeasurementLikelihoodMDThreshold_;
  double effNParticleThreshold_;
  int minInterSampleTimesteps_;
  double gaussianMergingThreshold_;
  double gaussianMergingCovarianceInflationFactor_;
  double gaussianPruningThreshold_;
  int importanceWeightingEvalPointCount_;
  bool reportTimingInfo_;

  unsigned int nThreadsPropagationStep_;
  bool logToFile_;

};




int main(int argc, char* argv[]){

  int initRandSeed = 0;
  if( argc == 2 ){
    initRandSeed = boost::lexical_cast<int>(argv[1]);
  }

  Simulator1d sim;
  if( !sim.readConfigFile() ){
    return -1;
  }
  sim.generateTrajectory( initRandSeed );
  sim.generateOdometry();
  sim.generateLandmarks();
  sim.generateMeasurements();
  sim.exportSimData();
 
  sim.setupRBPHDFilter();
  sim.run();

  return 0;

}
