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

// Convert log file from old format to new format

#include <boost/filesystem.hpp>
#include <stdio.h>
#include <string>

int main(int argc, char* argv[]){

  if( argc != 2 ){
    printf("Change 2d simulation logs from the old format to the new format\n");
    printf("Usage: convertLogFiles DATA_DIR/\n");
    return 0;
  }
  const char* logDir = argv[1];
  printf("Log directory: %s\n", logDir);

  boost::filesystem::path dir(logDir);
  if(!exists(dir)){
    printf("Log directory %s does not exist\n", logDir);
    return 0;
  }

  std::string filenameLandmarkEstOld( logDir );
  std::string filenamePoseEstOld( logDir );
  std::string filenameLandmarkEstNew( logDir );
  std::string filenamePoseEstNew( logDir );
  filenameLandmarkEstNew += "landmarkEst.dat";
  filenamePoseEstNew += "particlePose.dat";
  filenameLandmarkEstOld += "landmarkEst.bak";
  filenamePoseEstOld += "particlePose.bak";
  boost::filesystem::rename( boost::filesystem::path(filenameLandmarkEstNew),
			     boost::filesystem::path(filenameLandmarkEstOld) );
  boost::filesystem::rename( boost::filesystem::path(filenamePoseEstNew),
			     boost::filesystem::path(filenamePoseEstOld) );
  FILE* pLandmarkEstFileNew = fopen(filenameLandmarkEstNew.data(), "w");
  FILE* pLandmarkEstFileOld = fopen(filenameLandmarkEstOld.data(), "r");
  FILE* pPoseEstFileNew = fopen(filenamePoseEstNew.data(), "w");
  FILE* pPoseEstFileOld = fopen(filenamePoseEstOld.data(), "r");


  // Process particle pose estimate file
  printf("Processing: %s\n", filenamePoseEstNew.data());
  int nTimesteps = 0;
  if (fscanf(pPoseEstFileOld, "Timesteps: %d\n", &nTimesteps) != 1)
    return 1;

  double time;
  int nParticles;
  double x, y, r, w;
  for(int k = 0; k < nTimesteps; k++){
    if(fscanf(pPoseEstFileOld, "k = %lf\n", &time) != 1)
      return 2;
    if(fscanf(pPoseEstFileOld, "nParticles = %d\n", &nParticles) != 1)
      return 3;
    for(int i = 0; i < nParticles; i++){
      if(fscanf(pPoseEstFileOld, "%lf %lf %lf %lf\n", &x, &y, &r, &w) != 4)
	return 4;
      fprintf(pPoseEstFileNew, "%lf %d %lf %lf %lf %lf\n", time, i, x, y, r, w); 
    }
  }
  fclose(pPoseEstFileOld);
  fclose(pPoseEstFileNew);

  // Process map estimate file
  printf("Processing: %s\n", filenameLandmarkEstNew.data());
  double Sxx, Sxy, Syx, Syy;
  if( fscanf(pLandmarkEstFileOld, "Timesteps: %d\n", &nTimesteps) != 1)
    return 5;
  if( fscanf(pLandmarkEstFileOld, "nParticles: %d\n", &nParticles) != 1)
    return 6;
  int nM, pid;
  while (fscanf(pLandmarkEstFileOld, "Timestep: %lf   Particle: %d   Map Size: %d\n", &time, &pid, &nM) == 3){
    for(int m = 0; m < nM; m++ ){
      if (fscanf(pLandmarkEstFileOld, "%lf %lf %lf %lf %lf %lf %lf\n", &x, &y, &Sxx, &Sxy, &Syx, &Syy, &w) != 7)
	return 7;
      fprintf(pLandmarkEstFileNew, "%lf %d %lf %lf %lf %lf %lf %lf\n", time, pid, x, y, Sxx, Sxy, Syy, w);
    }
  }
  fclose(pLandmarkEstFileOld);
  fclose(pLandmarkEstFileNew);
  return 0;
  
}

