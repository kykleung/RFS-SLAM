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
import matplotlib
import matplotlib.pyplot as plt
#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

# Necessary to generate Type 1 fonts for pdf figures as required by IEEE for paper submissions
#matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['axes.unicode_minus'] = False
#matplotlib.rcParams['text.usetex'] = True
#plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rcParams.update({'font.size': 20})

if len(sys.argv) < 2:
    print "Usage: plotError2dSim DATA_DIR\n";
    sys.exit(0);

# Open files
dataDir = sys.argv[1];
if dataDir[-1] != '/':
    dataDir += '/'
poseEstFile = 'poseEstError.dat';
poseEstFile = dataDir + poseEstFile;
if os.path.exists(poseEstFile):
    print('Opening ' + poseEstFile);
else:
    print(poseEstFile + ' does not exist')
    sys.exit(0);
mapEstFile = 'landmarkEstError.dat';
mapEstFile = dataDir + mapEstFile;
if os.path.exists(mapEstFile):
    print('Opening ' + mapEstFile);
else:
    print(mapEstFile + ' does not exist')
    sys.exit(0);


# Read File
print('Reading ' + poseEstFile);
poseEst = np.genfromtxt(poseEstFile);
poseTimesteps = poseEst[:,0];
poseErr_x = poseEst[:,1]*10;
poseErr_y = poseEst[:,2]*10;
poseErr_r = poseEst[:,3] * 180.0 / np.pi ;
poseErr_d = poseEst[:,4]*10;
print('Reading ' + mapEstFile);
mapEst = np.genfromtxt(mapEstFile);
mapTimesteps = mapEst[:,0];
landmarksMeasured = mapEst[:,1];
landmarksEstimated = mapEst[:,2];
errorOSPA = mapEst[:,3];

# Generate plots
plt.figure(1);
p1, = plt.plot(poseTimesteps[::10], poseErr_x[::10], 'r-');
p2, = plt.plot(poseTimesteps[::10], poseErr_y[::10], 'b-');
plt.legend([p1, p2], [r"$x$", r"$y$"], loc=4);
#plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')
plt.xlabel(r'Time [s]');
plt.ylabel(r'Position error [m]');
plt.grid(True);
plt.ylim(-3, 3);
#plt.ylim(-15, 15);

plt.figure(2);
plt.plot(poseTimesteps[::10], poseErr_r[::10], 'r-');
plt.xlabel(r'Time [s]');
plt.ylabel(r'Rotation error [deg]');
plt.grid(True);
plt.ylim(-6, 6);
#plt.ylim(-15, 15);

plt.figure(3);
plt.plot(poseTimesteps[::10], poseErr_d[::10], 'r-');
plt.xlabel(r'Time [s]');
plt.ylabel(r'Position error [m]');
plt.grid(True);
plt.ylim(ymax = 3);
#plt.ylim(ymax = 15);

plt.figure(4);
plt.plot(mapTimesteps[::10], errorOSPA[::10], 'r-');
plt.xlabel(r'Time [s]');
plt.ylabel(r'OSPA error');
plt.grid(True);
plt.ylim(ymax = 4);

plt.figure(5);
p1, = plt.plot(mapTimesteps[::10], landmarksMeasured[::10], 'k-', linewidth=3.0);
p2, = plt.plot(mapTimesteps[::10], landmarksEstimated[::10], 'r-');
plt.legend([p1, p2], [r"Actual", r"Estimated"], loc=4);
#plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')
plt.xlabel(r'Time [s]');
plt.ylabel(r'Number of landmarks observed');
plt.grid(True);
plt.ylim(ymax = 70);

fig, ax1 = plt.subplots()
p1, = ax1.plot(poseTimesteps[::10], poseErr_x[::10], 'r-');
p2, = ax1.plot(poseTimesteps[::10], poseErr_y[::10], 'b-');
ax2 = ax1.twinx();
p3, = ax2.plot(poseTimesteps[::10], poseErr_r[::10], 'g-')

plt.legend([p1, p2, p3], ["x-position", "y-position", r"orientation"], loc=3);
plt.setp(plt.gca().get_legend().get_texts(), fontsize='20')
ax1.set_xlabel("Time [s]", fontsize=24);
ax1.set_ylabel("Position error [m]", fontsize=24);
ax2.set_ylabel("Orientation error [deg]", fontsize=24);
ax1.grid(True);
ax1.set_ylim(-3, 3);
ax2.set_ylim(-6, 6);
#ax1.set_ylim(-15, 15);
#ax2.set_ylim(-15, 15);
#fig.show();

# Save plots
errorPosePosFile = dataDir + 'errorPosePos.pdf';
errorPoseRotFile = dataDir + 'errorPoseRot.pdf';
errorPoseFile = dataDir + 'errorPose.pdf';
errorPoseDstFile = dataDir + 'errorPoseDst.pdf';
errorOSPAFile = dataDir + 'errorOSPA.pdf';
errorCardinalityFile = dataDir + 'errorCardinality.pdf';
plt.figure(1);
print('Saving  ' + errorPosePosFile);
plt.savefig(errorPosePosFile, format='pdf', bbox_inches='tight');
plt.figure(2);
print('Saving  ' + errorPoseRotFile);
plt.savefig(errorPoseRotFile, format='pdf', bbox_inches='tight');
plt.figure(3);
print('Saving  ' + errorPoseDstFile);
plt.savefig(errorPoseDstFile, format='pdf', bbox_inches='tight');
plt.figure(4);
print('Saving  ' + errorOSPAFile); 
plt.savefig(errorOSPAFile, format='pdf', bbox_inches='tight');
plt.figure(5);
print('Saving  ' + errorCardinalityFile); 
plt.savefig(errorCardinalityFile, format='pdf', bbox_inches='tight');
plt.figure(6);
print('Saving  ' + errorPoseFile);
plt.savefig(errorPoseFile, format='pdf', bbox_inches='tight');

