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

#include <assert.h>
#define BOOST_NO_CXX11_SCOPED_ENUMS // required for boost/filesystem to work with C++11
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "ProcessModel_Ackerman2D.hpp"
#include "RBPHDFilter.hpp"
#include "KalmanFilter_VictoriaPark.hpp"
#include <stdio.h>
#include <string>
#include <sstream>
#include <sys/ioctl.h>

#ifdef _PERFTOOLS_CPU
#include <gperftools/profiler.h>
#endif
#ifdef _PERFTOOLS_HEAP
#include <gperftools/heap-profiler.h>
#endif

using namespace rfs;

/**
 * \class RBPHDSLAM_VictoriaPark
 * \brief Run the RB-PHD-SLAM Algorithm on the Victoria Park Dataset.
 * \author Keith Leung
 */
class RBPHDSLAM_VictoriaPark{

public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  typedef RBPHDFilter<MotionModel_Ackerman2d, 
		      StaticProcessModel<Landmark3d>,
		      MeasurementModel_VictoriaPark,  
		      KalmanFilter_VictoriaPark> SLAM_Filter;

  RBPHDSLAM_VictoriaPark(){
    pFilter_ = NULL;
  }
  
  ~RBPHDSLAM_VictoriaPark(){
    
    if(pFilter_ != NULL){
      delete pFilter_;
    }

  }

  /** Read the simulator configuration file */
  bool readConfigFile(const char* fileName){
    
    cfgFileName_ = fileName;

    boost::property_tree::ptree pt;
    boost::property_tree::xml_parser::read_xml(fileName, pt);

    std::string dataDir = pt.get<std::string>("config.dataset.directory");
    if(*(dataDir.rbegin()) != '/'){
      dataDir += '/';
    }
    dataFileGPS_ = dataDir + pt.get<std::string>("config.dataset.filename.gps");
    dataFileDetection_ = dataDir + pt.get<std::string>("config.dataset.filename.detection");
    dataFileLidar_ = dataDir + pt.get<std::string>("config.dataset.filename.lidar");
    dataFileInput_ = dataDir + pt.get<std::string>("config.dataset.filename.input");
    dataFileSensorManager_ = dataDir + pt.get<std::string>("config.dataset.filename.manager");

    logResultsToFile_ = false;
    if( pt.get("config.logging.logResultsToFile", 0) == 1 )
      logResultsToFile_ = true;
       logTimingToFile_ = false;
    if( pt.get("config.logging.logTimingToFile", 0) == 1 )
      logTimingToFile_ = true;
    logDirPrefix_ = pt.get<std::string>("config.logging.logDirPrefix", "./");
    if( *logDirPrefix_.rbegin() != '/')
      logDirPrefix_ += '/';

    ackerman_h_ = pt.get<double>("config.process.AckermanModel.rearWheelOffset");
    ackerman_l_ = pt.get<double>("config.process.AckermanModel.frontToRearDist");
    ackerman_dx_ = pt.get<double>("config.process.AckermanModel.sensorOffset_x");
    ackerman_dy_ = pt.get<double>("config.process.AckermanModel.sensorOffset_y");
    var_uv_ = pt.get<double>("config.process.varuv");
    var_ur_ = pt.get<double>("config.process.varur");
    scale_ur_ = pt.get<double>("config.process.ur_scale");

    varlmx_ = pt.get<double>("config.landmarks.varlmx");
    varlmy_ = pt.get<double>("config.landmarks.varlmy");
    varlmd_ = pt.get<double>("config.landmarks.varlmd");

    rangeLimitMax_ = pt.get<double>("config.measurements.rangeLimitMax");
    rangeLimitMin_ = pt.get<double>("config.measurements.rangeLimitMin");
    bearingLimitMax_ = pt.get<double>("config.measurements.bearingLimitMax");
    bearingLimitMin_ = pt.get<double>("config.measurements.bearingLimitMin");
    clutterExpected_ = pt.get<double>("config.measurements.expectedNClutter");
    clutterAdded_ = pt.get("config.measurements.addedClutter", 0);
    bufferZonePd_ = pt.get<double>("config.measurements.bufferZonePd");
    varzr_ = pt.get<double>("config.measurements.varzr");
    varzb_ = pt.get<double>("config.measurements.varzb");
    varzd_ = pt.get<double>("config.measurements.varzd");
    varza_ = pt.get<double>("config.measurements.varza");
    BOOST_FOREACH(const boost::property_tree::ptree::value_type& v, 
		  pt.get_child("config.measurements.Pd")){
      Pd_.push_back( boost::lexical_cast<double>(v.second.data()) );
    }

    nMessageToProcess_ = pt.get<int>("config.filter.nMsgToProcess");

    nParticles_ = pt.get("config.filter.nParticles", 200);

    pNoiseInflation_ = pt.get("config.filter.predict.processNoiseInflationFactor", 1.0);
    birthGaussianWeight_ = pt.get("config.filter.predict.birthGaussian.Weight", 0.01);
    birthGaussianSupportMeasurementDist_ = pt.get<double>("config.filter.predict.birthGaussian.SupportMeasurementDist");
    birthGaussianSupportMeasurementThreshold_ = pt.get<int>("config.filter.predict.birthGaussian.SupportMeasurementThreshold");
    birthGaussianCheckCountThreshold_ = pt.get<int>("config.filter.predict.birthGaussian.CheckCountThreshold");
    birthGaussianCurrentMeasurementCountThreshold_ = pt.get<double>("config.filter.predict.birthGaussian.CurrentMeasurementCountThreshold");

    zNoiseInflation_ = pt.get("config.filter.update.measurementNoiseInflationFactor", 1.0);
    innovationRangeThreshold_ = pt.get<double>("config.filter.update.KalmanFilter.innovationThreshold.range");
    innovationBearingThreshold_ = pt.get<double>("config.filter.update.KalmanFilter.innovationThreshold.bearing");
    newGaussianCreateInnovMDThreshold_ = pt.get<double>("config.filter.update.GaussianCreateInnovMDThreshold");

    importanceWeightingEvalPointCount_ = pt.get("config.filter.weighting.nEvalPt", 15);
    importanceWeightingEvalPointGuassianWeight_ = pt.get("config.filter.weighting.minWeight", 0.75);
    importanceWeightingMeasurementLikelihoodMDThreshold_ = pt.get("config.filter.weighting.threshold", 3.0);
    useClusterProcess_ = false;
    if( pt.get("config.filter.weighting.useClusterProcess", 0) == 1 )
      useClusterProcess_ = true;

    effNParticleThreshold_ = pt.get("config.filter.resampling.effNParticle", nParticles_);
    minUpdatesBeforeResample_ = pt.get("config.filter.resampling.minTimesteps", 1);
    minMeasurementsBeforeResample_ = pt.get("config.filter.resampling.minMeasurements", 0);
    
    gaussianMergingThreshold_ = pt.get<double>("config.filter.merge.threshold");
    gaussianMergingCovarianceInflationFactor_ = pt.get("config.filter.merge.covInflationFactor", 1.0);
    
    gaussianPruningThreshold_ = pt.get("config.filter.prune.threshold", birthGaussianWeight_);

    // Copy config file to logDir
    if(logResultsToFile_ || logTimingToFile_){
      boost::filesystem::path dir(logDirPrefix_);
      boost::filesystem::create_directories(dir);
      boost::filesystem::path cfgFilePathSrc( cfgFileName_ );
      std::string cfgFileDst( logDirPrefix_ );
      cfgFileDst += "settings.xml";
      boost::filesystem::path cfgFilePathDst( cfgFileDst.data() );
      boost::filesystem::copy_file( cfgFilePathSrc, cfgFilePathDst, boost::filesystem::copy_option::overwrite_if_exists);
    }

    return true;   
  }

  struct SensorManagerMsg{
    TimeStamp t;
    enum Type {GPS=1, Input=2, Lidar=3};
    Type sensorType;
    uint idx;
  };

  struct LidarScanMsg{
    TimeStamp t;
    std::vector<double> scan;
  };

  /** \brief Import dataset from files */
  void readData(){

    std::string msgLine;

    // Read sensor manager log
    std::cout << "Reading input file: " << dataFileSensorManager_ << std::endl;
    std::ifstream file_sensorManager( dataFileSensorManager_.c_str() );
    assert( file_sensorManager.is_open() );
    while( std::getline( file_sensorManager, msgLine ) ){
      SensorManagerMsg msg;
      double time;
      int type;
      std::stringstream ss( msgLine );
      ss >> time >> type >> msg.idx;
      msg.t.setTime(time);
      msg.sensorType = (SensorManagerMsg::Type)type;
      msg.idx--;
      sensorManagerMsgs_.push_back(msg);
      //std::cout << std::setw(10) << std::fixed << std::setprecision(3) << msg.t.getTimeAsDouble() 
      //	<< std::setw(10) << (int)(msg.sensorType) 
      //	<< std::setw(10) << msg.idx << std::endl;
    }
    file_sensorManager.close();

    // Read Process Model Inputs    
    std::cout << "Reading input file: " << dataFileInput_ << std::endl;
    std::ifstream file_input( dataFileInput_.c_str() );
    assert( file_input.is_open() );
    while( std::getline( file_input, msgLine ) ){
      double time, vel, steer;
      std::stringstream ss( msgLine );
      ss >> time >> vel >> steer;
      SLAM_Filter::TInput::Vec uVec;
      uVec << vel, steer;
      SLAM_Filter::TInput u(uVec, TimeStamp(time)); 
      motionInputs_.push_back( u );
      //std::cout << std::setw(10) << std::fixed << std::setprecision(3) << u.getTime().getTimeAsDouble() 
      //	<< std::setw(10) << u[0] 
      //	<< std::setw(10) << std::fixed << std::setprecision(4) << u[1] << std::endl;
    }
    file_input.close();
    if(logResultsToFile_){
      boost::filesystem::path src( dataFileInput_.c_str() );
      std::string dst( logDirPrefix_ );
      dst += "inputs.dat";
      boost::filesystem::path cfgFilePathDst( dst.c_str() );
      boost::filesystem::copy_file( src, dst, boost::filesystem::copy_option::overwrite_if_exists);
    }

    // Read Lidar detections
    std::cout << "Reading input file: " << dataFileDetection_ << std::endl;
    std::ifstream file_measurements( dataFileDetection_.c_str() );
    assert( file_measurements.is_open() );
    while( std::getline( file_measurements, msgLine ) ){
      double time, range, bearing, diameter;
      std::stringstream ss(msgLine);
      ss >> time >> range >> bearing >> diameter;
      SLAM_Filter::TMeasurement::Vec zVec;
      zVec << range, bearing, diameter;
      SLAM_Filter::TMeasurement z(zVec, TimeStamp(time));
      measurements_.push_back(z);
      //std::cout << std::setw(10) << std::fixed << std::setprecision(3) << z.getTime().getTimeAsDouble() 
      //	<< std::setw(10) << std::fixed << std::setprecision(5) << z[0] 
      //<< std::setw(10) << std::fixed << std::setprecision(5) << z[1]
      //<< std::setw(10) << std::fixed << std::setprecision(5) << z[2] << std::endl;
    }
    file_measurements.close();
    if(logResultsToFile_){
      boost::filesystem::path src( dataFileDetection_.c_str() );
      std::string dst( logDirPrefix_ );
      dst += "measurements.dat";
      boost::filesystem::path cfgFilePathDst( dst.c_str() );
      boost::filesystem::copy_file( src, dst, boost::filesystem::copy_option::overwrite_if_exists);
    }

    // Read Lidar raw scans
    std::cout << "Reading input file: " << dataFileLidar_ << std::endl;
    std::ifstream file_lidar( dataFileLidar_.c_str() );
    assert( file_lidar.is_open() );
    double t = -1;
    LidarScanMsg msg;
    msg.scan.resize(361);
    file_lidar >> t; 
    while(t >= 0){
      msg.t.setTime(t);
      for(int i = 0; i < 361; i++){
	file_lidar >> msg.scan[i];
      }
      lidarScans_.push_back(msg);
      t = -1;
      file_lidar >> t;
      //std::cout << std::setw(10) << std::fixed << std::setprecision(3) << msg.t.getTimeAsDouble() 
      //		<< std::setw(10) << std::fixed << std::setprecision(5) << msg.scan[0] 
      //		<< std::setw(10) << std::fixed << std::setprecision(5) << msg.scan[360] << std::endl;
    }
    file_lidar.close();

    
    // Ground truth GPS
    //groundtruth_pose_.push_back( x_k );
    //groundtruth_pose_.back().setTime(t);
    std::cout << "Reading input file: " << dataFileGPS_ << std::endl;
    std::ifstream file_gps( dataFileGPS_.c_str() );
    assert( file_gps.is_open() );
    while( std::getline( file_gps, msgLine ) ){
      double time, x, y;
      std::stringstream ss(msgLine);
      ss >> time >> x >> y;
      Position2d::Vec pVec;
      pVec << x, y;
      Position2d p(pVec, TimeStamp(time));
      groundtruth_pos_.push_back(p);
      //std::cout << std::setw(10) << std::fixed << std::setprecision(3) << p.getTime().getTimeAsDouble() 
      //	<< std::setw(10) << std::fixed << std::setprecision(2) << p[0] 
      //	<< std::setw(10) << std::fixed << std::setprecision(2) << p[1] << std::endl;
    }
    file_gps.close();
    if(logResultsToFile_){
      boost::filesystem::path src( dataFileGPS_.c_str() );
      std::string dst( logDirPrefix_ );
      dst += "gps.dat";
      boost::filesystem::path cfgFilePathDst( dst.c_str() );
      boost::filesystem::copy_file( src, dst, boost::filesystem::copy_option::overwrite_if_exists);
    }

  }

  /** \brief Peform dead reckoning with process input data */
  void deadReckoning(){

    if(logResultsToFile_){
      std::ofstream drPoseFile( (logDirPrefix_ + "deadReckoning.dat").c_str() );
      std::cout << "Calculating dead reckoning estimate\n";
      std::cout << "Writing to: " << logDirPrefix_ + "deadReckoning.dat\n";
      SLAM_Filter::TPose x;
      TimeStamp t_k  = sensorManagerMsgs_[0].t;
      TimeStamp t_km = sensorManagerMsgs_[0].t;
      SLAM_Filter::TInput u_km; // last process input
      for(uint k = 0; k < sensorManagerMsgs_.size() ; k++ ){  
	if(sensorManagerMsgs_[k].sensorType == SensorManagerMsg::Input){
	  t_k = sensorManagerMsgs_[k].t;
	  TimeStamp dt = t_k - t_km;
	  SLAM_Filter::TInput::Vec u_km_vec = u_km.get();
	  u_km_vec[1] *= scale_ur_;
	  u_km.set(u_km_vec);
	  pFilter_->getProcessModel()->step( x, x, u_km, dt);
	  u_km = motionInputs_[ sensorManagerMsgs_[k].idx ];
	  
	  drPoseFile << std::fixed << std::setprecision(3) 
		     << std::setw(10) << sensorManagerMsgs_[k].t.getTimeAsDouble()
		     << std::setw(10) << x[0] 
		     << std::setw(10) << x[1] 
		     << std::setw(10) << x[2] << std::endl;
	  t_km = t_k;
	}
      }
      drPoseFile.close();
    }
  }

  /** RB-PHD Filter Setup */
  void setupRBPHDFilter(){

    pFilter_ = new SLAM_Filter( nParticles_ );

    // Configure process model
    pFilter_->getProcessModel()->setAckermanParams( ackerman_h_, ackerman_l_, ackerman_dx_, ackerman_dy_);
  
    // configure measurement model
    SLAM_Filter::TMeasurement::Cov R;
    R << varzr_, 0, 0, 0, varzb_, 0, 0, 0, varzd_;
    R *= zNoiseInflation_;
    pFilter_->getMeasurementModel()->setNoise(R, varza_);
    pFilter_->getMeasurementModel()->config.probabilityOfDetection_ = Pd_;
    pFilter_->getMeasurementModel()->config.expectedClutterNumber_ = clutterExpected_;
    pFilter_->getMeasurementModel()->config.rangeLimMax_ = rangeLimitMax_;
    pFilter_->getMeasurementModel()->config.rangeLimMin_ = rangeLimitMin_;
    pFilter_->getMeasurementModel()->config.bearingLimitMax_ = bearingLimitMax_ * PI / 180;
    pFilter_->getMeasurementModel()->config.bearingLimitMin_ = bearingLimitMin_ * PI / 180;
    pFilter_->getMeasurementModel()->config.bufferZonePd_ = bufferZonePd_;
    
    // configure the filter
    pFilter_->getKalmanFilter()->config.rangeInnovationThreshold_ = innovationRangeThreshold_;
    pFilter_->getKalmanFilter()->config.bearingInnovationThreshold_ = innovationBearingThreshold_;
    pFilter_->config.birthGaussianWeight_ = birthGaussianWeight_;
    pFilter_->config.birthGaussianMeasurementSupportDist_ = birthGaussianSupportMeasurementDist_;
    pFilter_->config.birthGaussianMeasurementCountThreshold_ = birthGaussianSupportMeasurementThreshold_;
    pFilter_->config.birthGaussianMeasurementCheckThreshold_ = birthGaussianCheckCountThreshold_;
    pFilter_->config.birthGaussianCurrentMeasurementCountThreshold_ = birthGaussianCurrentMeasurementCountThreshold_;
    pFilter_->setEffectiveParticleCountThreshold(effNParticleThreshold_);
    pFilter_->config.minUpdatesBeforeResample_ = minUpdatesBeforeResample_;
    pFilter_->config.minMeasurementsBeforeResample_ = minMeasurementsBeforeResample_;
    pFilter_->config.newGaussianCreateInnovMDThreshold_ = newGaussianCreateInnovMDThreshold_;
    pFilter_->config.importanceWeightingMeasurementLikelihoodMDThreshold_ = importanceWeightingMeasurementLikelihoodMDThreshold_;
    pFilter_->config.importanceWeightingEvalPointCount_ = importanceWeightingEvalPointCount_;
    pFilter_->config.importanceWeightingEvalPointGuassianWeight_ = importanceWeightingEvalPointGuassianWeight_;
    pFilter_->config.gaussianMergingThreshold_ = gaussianMergingThreshold_;
    pFilter_->config.gaussianMergingCovarianceInflationFactor_ = gaussianMergingCovarianceInflationFactor_;
    pFilter_->config.gaussianPruningThreshold_ = gaussianPruningThreshold_;
    pFilter_->config.useClusterProcess_ = useClusterProcess_;

  }

  /** \brief Process the data and peform SLAM */
  void run(){

#ifdef _PERFTOOLS_CPU
    std::string perfCPU_file = logDirPrefix_ + "rbphdslam_VictoriaPark_cpu.prof";
    ProfilerStart(perfCPU_file.data());
#endif
#ifdef _PERFTOOLS_HEAP
    std::string perfHEAP_file = logDirPrefix_ + "rbphdslam_VictoriaPark_heap.prof";
    HeapProfilerStart(perfHEAP_file.data());
#endif

    // Initialization at first timestep
    if(!logResultsToFile_){
      std::cout << "Note: results are NOT being logged to file (see config xml file)\n";
    }

    std::ofstream particlePoseFile;
    std::ofstream landmarkEstFile;
    std::ofstream clutterFile;
    if(logResultsToFile_){
      particlePoseFile.open( (logDirPrefix_ + "particlePose.dat").c_str() );
      landmarkEstFile.open( (logDirPrefix_ + "landmarkEst.dat").c_str() );
      clutterFile.open( (logDirPrefix_ + "clutter.dat").c_str() );
    }
    if(logResultsToFile_){
      MotionModel_Ackerman2d::TState x_i;
      for(int i = 0; i < pFilter_->getParticleCount(); i++){
	SLAM_Filter::TPose x = *(pFilter_->getParticleSet()->at(i));
	double w = pFilter_->getParticleSet()->at(i)->getWeight();
	particlePoseFile << std::fixed << std::setprecision(3) 
			 << std::setw(10) << sensorManagerMsgs_[0].t.getTimeAsDouble()
			 << std::setw(5) << i
			 << std::setw(10) << x[0] 
			 << std::setw(10) << x[1] 
			 << std::setw(10) << x[2] 
			 << std::setw(10) << w << std::endl;
      }
    }
    double expNegMeanClutter = exp( -clutterAdded_ );
    double poissonPmf[100];
    double poissonCmf[100];
    double mean_pow_i = 1;
    double i_factorial = 1;
    poissonPmf[0] = expNegMeanClutter;
    poissonCmf[0] = poissonPmf[0];
    for( int i = 1; i < 100; i++){
      mean_pow_i *= clutterAdded_;
      i_factorial *= i;
      poissonPmf[i] = mean_pow_i / i_factorial * expNegMeanClutter;
      poissonCmf[i] = poissonCmf[i-1] + poissonPmf[i]; 
    }

    // Process all sensor messages sequentially
    bool isInInitialStationaryState = true;
    TimeStamp t_km(0); // last sensor msg time
    SLAM_Filter::TInput u_km; // last process input
    SLAM_Filter::TInput::Cov u_km_cov;
    u_km_cov << var_uv_, 0, 0, var_ur_;
    u_km.setCov( u_km_cov );
    u_km.setTime( t_km );
    Landmark3d::Cov Q_m_k; // landmark process model additive noise
    int zIdx = 0;
    bool birthGaussianCheck = true;
    if(nMessageToProcess_ <= 0){
      nMessageToProcess_ = sensorManagerMsgs_.size();
    }else if(nMessageToProcess_ > sensorManagerMsgs_.size()){
      nMessageToProcess_ = sensorManagerMsgs_.size();
    }
    for(uint k = 0; k < nMessageToProcess_ ; k++ ){ 

      if( k % 500 == 0 || k == nMessageToProcess_ - 1){
	float progressPercent = float(k+1) / float(nMessageToProcess_);
	int progressBarW = 50;
	struct winsize ws;
	if(ioctl(1, TIOCGWINSZ, &ws) >= 0)
	  progressBarW = ws.ws_col - 30;
	int progressPos = progressPercent * progressBarW;
	if(progressBarW >= 50){
	  std::cout << "["; 
	  for(int i = 0; i < progressBarW; i++){
	    if(i < progressPos)
	      std::cout << "=";
	    else if(i == progressPos)
	      std::cout << ">";
	    else
	      std::cout << " ";
	  }
	  std::cout << "] ";
	}
	std::cout << "k = " << k << " (" << int(progressPercent * 100.0) << " %)\r"; 
	std::cout.flush();
      }
      if(k == nMessageToProcess_ - 1)
	std::cout << std::endl << std::endl;

#ifdef _PERFTOOLS_HEAP
      if( k % 50 == 0)
	HeapProfilerDump("Timestep interval dump");
#endif

      if(sensorManagerMsgs_[k].sensorType == SensorManagerMsg::Input){

	TimeStamp t_k = sensorManagerMsgs_[k].t;
	TimeStamp dt = t_k - t_km;

	Q_m_k << varlmx_, 0, 0, 0, varlmy_, 0, 0, 0, varlmd_; 
	Q_m_k = Q_m_k * dt.getTimeAsDouble() * dt.getTimeAsDouble();
	pFilter_->getLmkProcessModel()->setNoise(Q_m_k);

	if(isInInitialStationaryState){
	  pFilter_->predict( u_km, dt, false, false, birthGaussianCheck ); // this basically makes all initial particles sit still
	}else{
	  pFilter_->predict( u_km, dt, false, true, birthGaussianCheck); // true for use noise from u_km
	}
	birthGaussianCheck = false;
	
	u_km = motionInputs_[ sensorManagerMsgs_[k].idx ];
	SLAM_Filter::TInput::Vec u_km_vec = u_km.get();
	u_km_vec[1] *= scale_ur_;
	u_km.set(u_km_vec, u_km_cov);
	if(u_km[0] != 0){
	  isInInitialStationaryState = false;
	}

	t_km = t_k;

      }else if(sensorManagerMsgs_[k].sensorType == SensorManagerMsg::Lidar){

	TimeStamp t_k = sensorManagerMsgs_[k].t;
	TimeStamp dt = t_k - t_km;

	// Propagate particles up to lidar scan time

	Q_m_k << varlmx_, 0, 0, 0, varlmy_, 0, 0, 0, varlmd_; 
	Q_m_k = Q_m_k * dt.getTimeAsDouble() * dt.getTimeAsDouble();
	pFilter_->getLmkProcessModel()->setNoise(Q_m_k);

	if(isInInitialStationaryState){
	  pFilter_->predict( u_km, dt, false, false, birthGaussianCheck ); // this basically makes all initial particles sit still
	}else{
	  pFilter_->predict( u_km, dt, false, true, birthGaussianCheck ); // true for use noise from u_km
	}	
	birthGaussianCheck = false;
	  
	// Update particles with lidar scan data

	std::vector<SLAM_Filter::TMeasurement> Z;
	while( zIdx < measurements_.size() && measurements_[zIdx].getTime() == t_k ){
	  Z.push_back(measurements_[zIdx]);
	  zIdx++;
	}

	if(clutterAdded_ > 0 ){
	  double randomNum = drand48();
	  int nClutterToGen = 0;
	  while( randomNum > poissonCmf[ nClutterToGen ] ){
	    nClutterToGen++;
	  }
	  for( int i = 0; i < nClutterToGen; i++ ){
	
	    double r = drand48() * (rangeLimitMax_ - rangeLimitMin_) + rangeLimitMin_;
	    double b = drand48() * ((bearingLimitMax_ - bearingLimitMin_) + bearingLimitMin_) * PI / 180;
	    double d = 1.0;
	    SLAM_Filter::TMeasurement z_clutter;
	    SLAM_Filter::TMeasurement::Vec z;
	    z << r, b, d;
	    z_clutter.set(z, t_k);
	    Z.push_back(z_clutter);
	    if(logResultsToFile_){
	      clutterFile << std::fixed << std::setprecision(3)
			  << std::setw(10) << t_k.getTimeAsDouble()
			  << std::setprecision(5)
			  << std::setw(11) << r
			  << std::setw(11) << b 
			  << std::setw(11) << d << std::endl;
	    }
	  }
	}

	pFilter_->getMeasurementModel()->setLaserScan( lidarScans_[sensorManagerMsgs_[k].idx].scan );
	pFilter_->update(Z);
	birthGaussianCheck = true;

	// Log data
	if (logResultsToFile_){
    
	  double w_max = 0;
	  uint i_w_max = 0;
	  for(int i = 0; i < pFilter_->getParticleCount(); i++){
	    SLAM_Filter::TPose x = *(pFilter_->getParticleSet()->at(i));
	    double w = pFilter_->getParticleSet()->at(i)->getWeight();
	    particlePoseFile << std::fixed << std::setprecision(3) 
			     << std::setw(10) << sensorManagerMsgs_[k].t.getTimeAsDouble()
			     << std::setw(5) << i
			     << std::setw(10) << x[0] 
			     << std::setw(10) << x[1] 
			     << std::setw(10) << x[2] 
			     << std::setw(10) << w << std::endl;
	    if(w > w_max){
	      w_max = w;
	      i_w_max = i;
	    }
	  }
	  for( int m = 0; m < pFilter_->getGMSize(i_w_max); m++ ){
	    MeasurementModel_VictoriaPark::TLandmark::Vec u;
	    MeasurementModel_VictoriaPark::TLandmark::Mat S;
	    double w;
	    pFilter_->getLandmark(i_w_max, m, u, S, w);
	    landmarkEstFile << std::fixed << std::setprecision(3) 
			    << std::setw(10) << sensorManagerMsgs_[k].t.getTimeAsDouble()
			    << std::setw(5) << i_w_max
			    << std::setw(10) << u(0) 
			    << std::setw(10) << u(1)
			    << std::setw(10) << S(0,0) 
			    << std::setw(10) << S(0,1)
			    << std::setw(10) << S(1,1) 
			    << std::setw(10) << w << std::endl;
	  }

	}

	t_km = t_k;
	
      }

    }

    // Use the highest-weight particle and get its trajectory
    if(logResultsToFile_){
      std::ofstream bestTrajFile;
      bestTrajFile.open( (logDirPrefix_ + "trajectory.dat").c_str() );

      double w_max = 0;
      uint i_w_max = 0;
      for(int i = 0; i < pFilter_->getParticleCount(); i++){
	double w = pFilter_->getParticle(i)->getWeight();
	if(w > w_max){
	  w_max = w;
	  i_w_max = i;
	}
      }
      std::vector<SLAM_Filter::TTrajectory> traj;
      traj.push_back( *(pFilter_->getParticle(i_w_max)) );
      boost::shared_ptr<SLAM_Filter::TTrajectory> x_prev = traj[0].prev;
      while(x_prev != NULL){
	traj.push_back( *x_prev );
	x_prev = x_prev->prev;
      }
      for(int k = traj.size() - 1; k >= 0; k--){
	bestTrajFile << std::fixed << std::setprecision(3) 
		     << std::setw(10) << traj[k].getTime().getTimeAsDouble()
		     << std::setw(10) << traj[k][0] 
		     << std::setw(10) << traj[k][1] 
		     << std::setw(10) << traj[k][2] << std::endl;
      }

      bestTrajFile.close();
    }
    
#ifdef _PERFTOOLS_HEAP
    HeapProfilerStop();
#endif
#ifdef _PERFTOOLS_CPU
    ProfilerStop();
#endif

    
    std::cout << "Elapsed Timing Information [nsec]\n";
    std::cout << std::setw(15) << std::left << "Prediction" << std::setw(15)
	      << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->predict_wall
	      << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->predict_cpu << std::endl;
    std::cout << std::setw(15) << std::left << "Map Update" << std::setw(15)
	      << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->mapUpdate_wall
	      << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->mapUpdate_cpu << std::endl;
    std::cout << std::setw(15) << std::left << "Map Update (KF)" << std::setw(15)
	      << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->mapUpdate_kf_wall
	      << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->mapUpdate_kf_cpu << std::endl;
    std::cout << std::setw(15) << std::left << "Weighting" << std::setw(15)
	      << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->particleWeighting_wall
	      << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->particleWeighting_cpu << std::endl;
    std::cout << std::setw(15) << std::left << "Map Merge" << std::setw(15)
	      << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->mapMerge_wall
	      << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->mapMerge_cpu << std::endl;
    std::cout << std::setw(15) << std::left << "Map Prune" << std::setw(15)
	      << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->mapPrune_wall
	      << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->mapPrune_cpu << std::endl;
    std::cout << std::setw(15) << std::left << "Resampling" << std::setw(15)
	      << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->particleResample_wall
	      << std::setw(6) << std::left << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->particleResample_cpu << std::endl;
    std::cout << std::setw(15) << std::left << "Total" << std::setw(15)
	      << std::setw(6) << std::right << "wall:" << std::setw(15)
	      << pFilter_->getTimingInfo()->predict_wall +
                 pFilter_->getTimingInfo()->mapUpdate_wall +
                 pFilter_->getTimingInfo()->particleWeighting_wall +
                 pFilter_->getTimingInfo()->mapMerge_wall +
                 pFilter_->getTimingInfo()->mapPrune_wall +
                 pFilter_->getTimingInfo()->particleResample_wall 
	      << std::setw(6) << std::right << "cpu:" << std::setw(15)
	      << pFilter_->getTimingInfo()->predict_cpu +
                 pFilter_->getTimingInfo()->mapUpdate_cpu +
                 pFilter_->getTimingInfo()->particleWeighting_cpu +
                 pFilter_->getTimingInfo()->mapMerge_cpu +
                 pFilter_->getTimingInfo()->mapPrune_cpu + 
                 pFilter_->getTimingInfo()->particleResample_cpu << std::endl;

    if(logTimingToFile_){
      std::ofstream timingFile( (logDirPrefix_ + "timing.dat").data() );
      timingFile << "Elapsed Timing Information [nsec]\n";
      timingFile << std::setw(15) << std::left << "Prediction" << std::setw(15)
		 << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->predict_wall
		 << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->predict_cpu << std::endl;
      timingFile << std::setw(15) << std::left << "Map Update" << std::setw(15)
		 << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->mapUpdate_wall
		 << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->mapUpdate_cpu << std::endl;
      timingFile << std::setw(15) << std::left << "Map Update (KF)" << std::setw(15)
		 << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->mapUpdate_kf_wall
		 << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->mapUpdate_kf_cpu << std::endl;
      timingFile << std::setw(15) << std::left << "Weighting" << std::setw(15)
		 << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->particleWeighting_wall
		 << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->particleWeighting_cpu << std::endl;
      timingFile << std::setw(15) << std::left << "Map Merge" << std::setw(15)
		 << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->mapMerge_wall
		 << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->mapMerge_cpu << std::endl;
      timingFile << std::setw(15) << std::left << "Map Prune" << std::setw(15)
		 << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->mapPrune_wall
		 << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->mapPrune_cpu << std::endl;
      timingFile << std::setw(15) << std::left << "Resampling" << std::setw(15)
		 << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->particleResample_wall
		 << std::setw(6) << std::left << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->particleResample_cpu << std::endl;
      timingFile << std::setw(15) << std::left << "Total" << std::setw(15)
		 << std::setw(6) << std::right << "wall:" << std::setw(15)
		 << pFilter_->getTimingInfo()->predict_wall +
	pFilter_->getTimingInfo()->mapUpdate_wall +
	pFilter_->getTimingInfo()->particleWeighting_wall +
	pFilter_->getTimingInfo()->mapMerge_wall +
	pFilter_->getTimingInfo()->mapPrune_wall +
	pFilter_->getTimingInfo()->particleResample_wall 
		 << std::setw(6) << std::right << "cpu:" << std::setw(15)
		 << pFilter_->getTimingInfo()->predict_cpu +
	pFilter_->getTimingInfo()->mapUpdate_cpu +
	pFilter_->getTimingInfo()->particleWeighting_cpu +
	pFilter_->getTimingInfo()->mapMerge_cpu +
	pFilter_->getTimingInfo()->mapPrune_cpu + 
	pFilter_->getTimingInfo()->particleResample_cpu << std::endl;
      timingFile.close();
    }
    
    if(logResultsToFile_){
      particlePoseFile.close();
      landmarkEstFile.close();
      clutterFile.close();
    }
  }

private:

  const char* cfgFileName_;

  // Dataset
  std::string dataFileGPS_;
  std::string dataFileDetection_;
  std::string dataFileLidar_;
  std::string dataFileInput_;
  std::string dataFileSensorManager_;
  std::vector<SensorManagerMsg> sensorManagerMsgs_;

  // Process model
  double ackerman_h_;
  double ackerman_l_;
  double ackerman_dx_;
  double ackerman_dy_;
  double var_uv_;
  double var_ur_;
  double scale_ur_;
  std::vector<Position2d> groundtruth_pos_;
  std::vector<MotionModel_Ackerman2d::TInput> motionInputs_;
  std::vector<MotionModel_Ackerman2d::TState> deadReckoning_pose_;

  // Landmarks 
  double varlmx_;
  double varlmy_;
  double varlmd_;

  // Range-Bearing Measurements
  double rangeLimitMax_;
  double rangeLimitMin_;
  double bearingLimitMax_;
  double bearingLimitMin_;
  double clutterExpected_;
  double clutterAdded_;
  double bufferZonePd_;
  double varzr_;
  double varzb_;
  double varzd_;
  double varza_;
  std::vector<SLAM_Filter::TMeasurement> measurements_;
  std::vector<double> Pd_;
  std::vector<LidarScanMsg> lidarScans_;

  // Filters
  KalmanFilter_VictoriaPark kf_;
  SLAM_Filter *pFilter_; 
  int nParticles_;
  double pNoiseInflation_;
  double zNoiseInflation_;
  double innovationRangeThreshold_;
  double innovationBearingThreshold_;
  double birthGaussianWeight_;
  double birthGaussianSupportMeasurementDist_;
  int birthGaussianSupportMeasurementThreshold_;
  int birthGaussianCheckThreshold_;
  int birthGaussianCheckCountThreshold_;
  double birthGaussianCurrentMeasurementCountThreshold_;
  double newGaussianCreateInnovMDThreshold_;
  double importanceWeightingMeasurementLikelihoodMDThreshold_;
  double importanceWeightingEvalPointGuassianWeight_;
  double effNParticleThreshold_;
  int minUpdatesBeforeResample_;
  int minMeasurementsBeforeResample_;
  double gaussianMergingThreshold_;
  double gaussianMergingCovarianceInflationFactor_;
  double gaussianPruningThreshold_;
  int importanceWeightingEvalPointCount_;
  bool useClusterProcess_;

  bool logResultsToFile_;
  bool logTimingToFile_;
  int nMessageToProcess_;

public:
  std::string logDirPrefix_;
};



int main(int argc, char* argv[]){

  RBPHDSLAM_VictoriaPark slam;

  int seed = time(NULL);
  srand(seed);
  std::string cfgFileName;
  boost::program_options::options_description desc("Options");
  desc.add_options()
    ("help,h", "produce this help message")
    ("cfg,c", boost::program_options::value<std::string>(&cfgFileName)->default_value("cfg/rbphdslam_VictoriaPark.xml"), "configuration xml file")
    ("seed,s", boost::program_options::value<int>(&seed), "random seed for running the simulation (default: based on current system time)");
  boost::program_options::variables_map vm;
  boost::program_options::store( boost::program_options::parse_command_line(argc, argv, desc), vm);
  boost::program_options::notify(vm);

  if( vm.count("help") ){
    std::cout << desc << "\n";
    return 1;
  }

  if( vm.count("cfg") ){
    cfgFileName = vm["cfg"].as<std::string>();
  }
  std::cout << "Configuration file: " << cfgFileName << std::endl;
  if( !slam.readConfigFile( cfgFileName.data() ) ){
    std::cout << "[Error] Unable to read config file: " << cfgFileName << std::endl;
    return -1;
  }

  slam.readData();
  slam.setupRBPHDFilter();
  slam.deadReckoning();

  if( vm.count("seed") ){
    seed = vm["seed"].as<int>();
    std::cout << "Simulation random seed manually set to: " << seed << std::endl;
  }
  srand48( seed );

  //boost::timer::auto_cpu_timer *timer = new boost::timer::auto_cpu_timer(6, "Run time: %ws\n");

  slam.run(); 

  //delete timer;
 
  return 0;

}
