#!/usr/bin/python

 #
 # Software License Agreement (New BSD License)
 #
 # Copyright (c) 2013, Keith Leung
 # All rights reserved.
 # 
 # Redistribution and use in source and binary forms, with or without
 # modification, are permitted provided that the following conditions are met:
 #     * Redistributions of source code must retain the above copyright
 #       notice, this list of conditions and the following disclaimer.
 #     * Redistributions in binary form must reproduce the above copyright
 #       notice, this list of conditions and the following disclaimer in the
 #       documentation and/or other materials provided with the distribution.
 #     * Neither the name of the Advanced Mining Technology Center (AMTC), the
 #       Universidad de Chile, nor the names of its contributors may be 
 #       used to endorse or promote products derived from this software without 
 #       specific prior written permission.
 # 
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 # ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 # WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 # DISCLAIMED. IN NO EVENT SHALL THE AMTC, UNIVERSIDAD DE CHILE, OR THE COPYRIGHT 
 # HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 # CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE 
 # GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
 # HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
 # LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
 # THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 #

import sys
import os.path
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

# Necessary to generate Type 1 fonts for pdf figures as required by IEEE for paper submissions
#matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
#matplotlib.rcParams['text.usetex'] = True
#plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 22})

if len(sys.argv) < 2:
    print "Usage: plotError2dSim DATA_PARENT_DIR\n";
    print "DATA_PARENT_DIR should contain rbphdslam, fastslam, mhfastslam\n"
    sys.exit(0);

# Open files
dataDir = sys.argv[1];
if dataDir[-1] != '/':
    dataDir += '/'
poseEstFile_rbphdslam = 'rbphdslam/poseEstError.dat';
poseEstFile_rbphdslam = dataDir + poseEstFile_rbphdslam;
poseEstFile_fastslam = 'fastslam/poseEstError.dat';
poseEstFile_fastslam = dataDir + poseEstFile_fastslam;
poseEstFile_mhfastslam = 'mhfastslam/poseEstError.dat';
poseEstFile_mhfastslam = dataDir + poseEstFile_mhfastslam;
if os.path.exists(poseEstFile_rbphdslam):
    print('Opening ' + poseEstFile_rbphdslam);
else:
    print(poseEstFile_rbphdslam + ' does not exist')
    sys.exit(0);
if os.path.exists(poseEstFile_fastslam):
    print('Opening ' + poseEstFile_fastslam);
else:
    print(poseEstFile_fastslam + ' does not exist')
    sys.exit(0);
if os.path.exists(poseEstFile_mhfastslam):
    print('Opening ' + poseEstFile_mhfastslam);
else:
    print(poseEstFile_mhfastslam + ' does not exist')
    sys.exit(0);


mapEstFile_rbphdslam = 'rbphdslam/landmarkEstError.dat';
mapEstFile_rbphdslam = dataDir + mapEstFile_rbphdslam;
if os.path.exists(mapEstFile_rbphdslam):
    print('Opening ' + mapEstFile_rbphdslam);
else:
    print(mapEstFile_rbphdslam + ' does not exist')
    sys.exit(0);
mapEstFile_fastslam = 'fastslam/landmarkEstError.dat';
mapEstFile_fastslam = dataDir + mapEstFile_fastslam;
if os.path.exists(mapEstFile_fastslam):
    print('Opening ' + mapEstFile_fastslam);
else:
    print(mapEstFile_fastslam + ' does not exist')
    sys.exit(0);
mapEstFile_mhfastslam = 'mhfastslam/landmarkEstError.dat';
mapEstFile_mhfastslam = dataDir + mapEstFile_mhfastslam;
if os.path.exists(mapEstFile_mhfastslam):
    print('Opening ' + mapEstFile_mhfastslam);
else:
    print(mapEstFile_mhfastslam + ' does not exist')
    sys.exit(0);


# Read Files
print('Reading ' + poseEstFile_rbphdslam);
poseEst = np.genfromtxt(poseEstFile_rbphdslam);
poseTimesteps = poseEst[:,0];
poseErr_d_rbphdslam = poseEst[:,4]*10;
#print poseErr_d_rbphdslam[::10];

print('Reading ' + poseEstFile_fastslam);
poseEst = np.genfromtxt(poseEstFile_fastslam);
poseTimesteps = poseEst[:,0];
poseErr_d_fastslam = poseEst[:,4]*10;
#print poseErr_d_fastslam[::10];

print('Reading ' + poseEstFile_mhfastslam);
poseEst = np.genfromtxt(poseEstFile_mhfastslam);
poseTimesteps = poseEst[:,0];
poseErr_d_mhfastslam = poseEst[:,4]*10;
#print poseErr_d_mhfastslam[::10];

poseEst = None;

print('Reading ' + mapEstFile_rbphdslam);
mapEst = np.genfromtxt(mapEstFile_rbphdslam);
mapTimesteps = mapEst[:,0];
landmarksMeasured = mapEst[:,1];
landmarksEstimated_rbphdslam = mapEst[:,2];
errorOSPA_rbphdslam = mapEst[:,3];

print('Reading ' + mapEstFile_fastslam);
mapEst = np.genfromtxt(mapEstFile_fastslam);
mapTimesteps = mapEst[:,0];
landmarksMeasured = mapEst[:,1];
landmarksEstimated_fastslam = mapEst[:,2];
errorOSPA_fastslam = mapEst[:,3];

print('Reading ' + mapEstFile_mhfastslam);
mapEst = np.genfromtxt(mapEstFile_mhfastslam);
mapTimesteps = mapEst[:,0];
landmarksMeasured = mapEst[:,1];
landmarksEstimated_mhfastslam = mapEst[:,2];
errorOSPA_mhfastslam = mapEst[:,3];

mapEst = None;

# Generate plots
plt.figure(1);
p1, = plt.plot(poseTimesteps[::10], poseErr_d_rbphdslam[::10], 'r-');
p2, = plt.plot(poseTimesteps[::10], poseErr_d_mhfastslam[::10], 'b-');
p3, = plt.plot(poseTimesteps[::10], poseErr_d_fastslam[::10], 'g-');
plt.legend([p1, p2, p3], [r"RB-PHD-SLAM", r"MH-FastSLAM", r"FastSLAM"], loc=2);
plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')
plt.xlabel(r'Time [s]');
plt.ylabel(r'Distance error [m]');
plt.grid(True);
plt.ylim(ymax = 3);
#plt.ylim(ymax = 15);

plt.figure(2, figsize=(12,6), facecolor='w');
p1, = plt.plot(mapTimesteps[::10], errorOSPA_rbphdslam[::10], 'r-');
p2, = plt.plot(mapTimesteps[::10], errorOSPA_mhfastslam[::10], 'b-');
p3, = plt.plot(mapTimesteps[::10], errorOSPA_fastslam[::10], 'g-');
plt.legend([p1, p2, p3], [r"RB-PHD-SLAM", r"MH-FastSLAM", r"FastSLAM"], loc=2);
plt.setp(plt.gca().get_legend().get_texts(), fontsize='18')
plt.xlabel(r'Time [s]');
plt.ylabel(r'OSPA error [m]');
plt.grid(True);
plt.ylim(ymax = 12);

plt.figure(3, figsize=(12,6), facecolor='w');
p4, = plt.plot(mapTimesteps[::10], landmarksMeasured[::10], 'k-', linewidth=5.0);
p1, = plt.plot(mapTimesteps[::10], landmarksEstimated_rbphdslam[::10], 'r-');
p2, = plt.plot(mapTimesteps[::10], landmarksEstimated_mhfastslam[::10], 'b-');
p3, = plt.plot(mapTimesteps[::10], landmarksEstimated_fastslam[::10], 'g-');
plt.legend([p4, p1, p2, p3], [r"Actual", r"RB-PHD-SLAM", r"MH-FastSLAM", r"FastSLAM"], loc=4);
plt.setp(plt.gca().get_legend().get_texts(), fontsize='18')
plt.xlabel(r'Time [s]');
plt.ylabel(r'Number of landmarks');
plt.grid(True);
plt.ylim(ymax = 70);

#plt.show();

# Save plots
errorPoseDstFile = dataDir + 'errorPoseDst.pdf';
errorOSPAFile = dataDir + 'errorOSPA.pdf';
errorCardinalityFile = dataDir + 'errorCardinality.pdf';
plt.figure(1);
print('Saving  ' + errorPoseDstFile);
plt.savefig(errorPoseDstFile, format='pdf', bbox_inches='tight');
plt.figure(2);
print('Saving  ' + errorOSPAFile); 
plt.savefig(errorOSPAFile, format='pdf', bbox_inches='tight');
plt.figure(3);
print('Saving  ' + errorCardinalityFile); 
plt.savefig(errorCardinalityFile, format='pdf', bbox_inches='tight');

