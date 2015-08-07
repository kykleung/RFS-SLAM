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
matplotlib.use("TkAgg");

#print matplotlib.__version__

import matplotlib.pyplot as plt
import matplotlib.animation as anim
from matplotlib.patches import Ellipse, Circle
from matplotlib import transforms
import matplotlib.ticker as ticker   

matplotlib.rcParams.update({'font.size': 20})

saveMovie = False;
saveFig = True
timestepStart = 2950

nLandmarksDrawMax = 500;
nMeasurementsDrawMax = 500;

if len(sys.argv) < 2:
    print "Usage: animate2dSim DATA_DIR\n";
    sys.exit(0);

# Setting for file names

dataDir = sys.argv[1];
if dataDir[-1] != '/':
    dataDir += '/'

gtPoseFile = 'gtPose.dat';
gtPoseFile = dataDir + gtPoseFile;
if os.path.exists(gtPoseFile):
    print('Opening ' + gtPoseFile);
else:
    print(gtPoseFile + ' does not exist')
    sys.exit(0);

gtMapFile = 'gtLandmark.dat';
gtMapFile = dataDir + gtMapFile;
if os.path.exists(gtMapFile):
    print('Opening ' + gtMapFile);
else:
    print(gtMapFile + ' does not exist')
    sys.exit(0);

estPoseFile = 'particlePose.dat';
estPoseFile = dataDir + estPoseFile;
if os.path.exists(estPoseFile):
    print('Opening ' + estPoseFile);
else:
    print(estPoseFile + ' does not exist');
    sys.exit(0);
estPoseFileHandle = open(estPoseFile, "r");

estMapFile = 'landmarkEst.dat';
estMapFile = dataDir + estMapFile;
if os.path.exists(estMapFile):
    print('Opening ' + estMapFile);
else:
    print(estMapFile + ' does not exist')
    sys.exit(0);
estMapFileHandle = open(estMapFile, "r");

measurementFile = 'measurement.dat';
measurementFile = dataDir + measurementFile;
if os.path.exists(measurementFile):
    print('Opening ' + measurementFile);
else:
    print(measurementFile + ' does not exist')
    sys.exit(0);
measurementFileHandle = open(measurementFile, "r");

estimateImageFile = 'estimate.pdf';
estimateImageFile = dataDir + estimateImageFile;
estimateMovieFile = 'estimate.mp4';
estimateMovieFile = dataDir + estimateMovieFile;

# Reading files

print('Reading ' + gtPoseFile);
gtPose = np.genfromtxt(gtPoseFile);
gtPose_t = gtPose[:,0];
gtPose_x = gtPose[:,1];
gtPose_y = gtPose[:,2];
gtPose_r = gtPose[:,3];
print('Number of poses: ' + str(len(gtPose)));

print('Reading ' + gtMapFile);
gtMap = np.genfromtxt(gtMapFile);
gtMap_x = gtMap[:,0];
gtMap_y = gtMap[:,1];


p = np.fromfile(estPoseFileHandle, dtype=float, count=6, sep=" ");
p_idx = 0;
p_idx_maxWeight = 0;
p_maxWeight = p[5];
p_x = [];
p_y = [];
p_w = [];
p_t = p[0];
while p[0] == p_t:
    p_x.append(p[2]);
    p_y.append(p[3]);
    p_w.append(p[5]);
    if p[5] > p_maxWeight :
        p_maxWeight = p[5];
        p_idx_maxWeight = p_idx;
    p = np.fromfile(estPoseFileHandle, dtype=float, count=6, sep=" ");
    p_idx += 1;
px_best = gtPose_x * 0;
py_best = gtPose_y * 0;

m = np.fromfile(estMapFileHandle, count=8, sep=" ", dtype=float);

z = np.fromfile(measurementFileHandle, count=3, sep=" ", dtype=float);

# Plotting 

fig = plt.figure( figsize=(12,10), facecolor='w')
gtMapHandle, = plt.plot(gtMap_x, gtMap_y, 'k.');
ax = plt.gca();
plt.axis('equal');
plt.grid(True);

gtPoseHandle, = plt.plot(gtPose_x, gtPose_y, 'r--');

gtPoseCurrentPos, = plt.plot([], [], 'ro');
gtPoseCurrentDir, = plt.plot([], [], 'r-');

measurements = [];
for i in range(0, nMeasurementsDrawMax) : 
    measurement_line, = plt.plot([], [], 'b-');
    measurements.append( measurement_line );

landmarks = [];
for i in range(0, nLandmarksDrawMax) : 
    landmark_ellipse = Ellipse(xy=(0,0), width=0, height=0, angle=0);
    landmarks.append(landmark_ellipse); 
    ax.add_patch(landmarks[i]);

particles, = plt.plot([], [], 'b.');

xLim = plt.getp(ax, 'xlim');
yLim = plt.getp(ax, 'ylim');
txt = plt.text(xLim[1]-1, yLim[1]-1, " ");

def animateInit():

    txt.set_text("Time: ");
    gtPoseCurrentPos.set_data([],[]);
    gtPoseCurrentDir.set_data([],[]);
    gtPoseHandle.set_data([],[]);
    gtMapHandle.set_data([],[]);
    particles.set_data([],[]);
    for i in range(0, nMeasurementsDrawMax) :
        measurements[i].set_data([],[]);
        measurements[i].set_color([0.5,0.5,0.9]);
    for i in range(0, nLandmarksDrawMax):
        landmarks[i].center = (0,0);
        landmarks[i].width = 0;
        landmarks[i].height = 0;
        landmarks[i].set_facecolor([0.2,0.2,0.8])
    return [];

def animate(i):
    
    global p;
    global m;
    global z;

    currentTime = i * 0.1;
    drawnObjects = [];

    # Time
    txt.set_text("Time: {0}".format(currentTime));
    drawnObjects.append(txt);
    
    timeStart = timestepStart * 0.1

    p_idx_maxWeight = 0
    p_maxWeight = 0
    while p.any() and p[0] < timeStart:
        p = np.fromfile(estPoseFileHandle, dtype=float, count=6, sep=" ");
        currentTime = p[0]
        k = int(np.ceil(currentTime / 0.1))
        if p[1] == 0:
            p_idx_maxWeight = p[1]
            p_maxWeight = p[5]
            px_best[k] = p[2];
            py_best[k] = p[3];
        elif p[5] > p_maxWeight :
            p_maxWeight = p[5];
            p_idx_maxWeight = p[1];
            px_best[k] = p[2];
            py_best[k] = p[3];
    
    while m.any() and m[0] < timeStart:
        m = np.fromfile(estMapFileHandle, count=8, sep=" ", dtype=float);
    while z.any() and z[0] < timeStart:
        z = np.fromfile(measurementFileHandle, count=3, sep=" ", dtype=float);


    # Groundtruth 
    gtPoseCurrentPos.set_xdata(gtPose_x[i]);
    gtPoseCurrentPos.set_ydata(gtPose_y[i]); # Append to drawn objects later so it shows above measurements
    gtPoseHead = [gtPose_x[i] + 0.3*np.cos(gtPose_r[i]), gtPose_y[i] + 0.3*np.sin(gtPose_r[i]) ];
    gtPoseCurrentDir.set_xdata([gtPose_x[i], gtPoseHead[0]]);
    gtPoseCurrentDir.set_ydata([gtPose_y[i], gtPoseHead[1]]); # Append to drawn objects later so it shows above measurements

    gtPoseHandle.set_data(gtPose_x, gtPose_y);
    gtMapHandle.set_data(gtMap_x, gtMap_y);

    # Particles
    p_idx = 0;
    p_idx_maxWeight = 0;
    p_maxWeight = p[5];
    px_best[i] = p[2];
    py_best[i] = p[3];
    p_x = [];
    p_y = [];
    p_w = [];
    while p.any() and abs(p[0] - currentTime) < 1e-12:
        p_x.append(p[2]);
        p_y.append(p[3]);
        p_w.append(p[5]);
        if p[5] > p_maxWeight :
            p_maxWeight = p[5];
            p_idx_maxWeight = p_idx;
            px_best[i] = p[2];
            py_best[i] = p[3];
        p = np.fromfile(estPoseFileHandle, dtype=float, count=6, sep=" ");
        p_idx += 1;
    particles.set_data(p_x, p_y);

    # Landmarks
    m_idx = 0;
    while m.any() and abs(m[0] - currentTime) < 1e-12:
        if round(m[1]) == round(p_idx_maxWeight):

            cov = np.array([ [ m[4], m[5] ], [ m[5], m[6] ] ]);
            w = m[7];
            eVal, eVec = np.linalg.eig(cov);
            eVal = eVal.real;
            a1 = 5*np.sqrt(eVal[0]); # Assume this is semi-major axis first
            a2 = 5*np.sqrt(eVal[1]); 
            semiMajorAxis = eVec[:,0];
            if a2 > a1:
                aTmp = a1
                a1 = a2
                a2 = aTmp
                semiMajorAxis = eVec[:,1];
            a1Angle = np.arctan2(semiMajorAxis[1], semiMajorAxis[0]);

            landmarks[m_idx].set_alpha(min(w, 0.75));
            landmarks[m_idx].center = (m[2], m[3]);
            landmarks[m_idx].height = a2;
            landmarks[m_idx].width = a1;
            t_start = ax.transData;
            t_rot = transforms.Affine2D().rotate_around(m[2], m[3], a1Angle);
            t_compound = t_rot + t_start;
            landmarks[m_idx].set_transform(t_compound);
            drawnObjects.append(landmarks[m_idx]);
            m_idx += 1;

        m = np.fromfile(estMapFileHandle, count=8, sep=" ", dtype=float);

    while landmarks[m_idx].height != 0:
        landmarks[m_idx].set_alpha(0);
        landmarks[m_idx].center = (0, 0);
        landmarks[m_idx].height = 0;
        landmarks[m_idx].width = 0;
        m_idx += 1;

    # Measurements
    nZ = 0;
    while z.any() and abs(z[0] -  currentTime) < 1e-12:
        z_dir = gtPose_r[i] + z[2];
        z_end = [gtPose_x[i] + z[1]*np.cos(z_dir), gtPose_y[i] + z[1]*np.sin(z_dir) ];
        measurements[nZ].set_data([gtPose_x[i], z_end[0]], [gtPose_y[i], z_end[1]]);
        drawnObjects.append(measurements[nZ]);
        z = np.fromfile(measurementFileHandle, count=3, sep=" ", dtype=float);
        nZ += 1;
    while measurements[nZ].get_xdata() != []:
        measurements[nZ].set_data([], []);
        nZ += 1;
    
    drawnObjects.append(gtPoseCurrentPos);
    drawnObjects.append(gtPoseCurrentDir);
    drawnObjects.append(gtPoseHandle);
    drawnObjects.append(gtMapHandle);
    drawnObjects.append(particles);

    return drawnObjects;

animation = anim.FuncAnimation(plt.figure(1), animate, np.arange(timestepStart, len(gtPose_t)), interval=1, 
                               init_func=animateInit, blit=True, repeat=False);
if saveMovie:
    FFMpegWriter = matplotlib.animation.writers['ffmpeg']
    animation.save(estimateMovieFile, writer=FFMpegWriter(fps = 30))
else:
    plt.show(block=False)

if saveFig:
    # Necessary to generate Type 1 fonts for pdf figures as required by IEEE for paper submissions
    #matplotlib.rcParams['pdf.fonttype'] = 42
    #matplotlib.rcParams['ps.fonttype'] = 42
    matplotlib.rcParams['ps.useafm'] = True
    matplotlib.rcParams['pdf.use14corefonts'] = True
    #matplotlib.rcParams['text.usetex'] = True
    #plt.rc('text', usetex=True)
    estPoseHandle, = plt.plot(px_best, py_best, 'b-')
    for i in range(0, nMeasurementsDrawMax) : 
        measurements[i].remove();
    plt.setp(gtPoseHandle, linewidth=2.0)
    txt.set_text(" ");
    plt.legend([gtPoseHandle, estPoseHandle, gtMapHandle, landmarks[0]], ["Ground-truth trajectory", "Estimated trajectory", "Ground-truth landmark", "Estimated landmark" ], loc=4);
    plt.setp(plt.gca().get_legend().get_texts(), fontsize='18')
    scale = 10;
    ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*scale))
    plt.gca().xaxis.set_major_formatter(ticks)
    ticks = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y*scale))
    plt.gca().yaxis.set_major_formatter(ticks)
    plt.savefig(estimateImageFile, format='pdf', bbox_inches='tight')
    #plt.savefig('estimate.eps', format='eps', bbox_inches='tight')

measurementFileHandle.close();


