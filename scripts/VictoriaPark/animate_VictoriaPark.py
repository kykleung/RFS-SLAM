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
from scipy import ndimage

import matplotlib
matplotlib.use("TKAgg");
#print matplotlib.__version__

# Necessary to generate Type 1 fonts for pdf figures as required by IEEE for paper submissions
#matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
#matplotlib.rcParams['text.usetex'] = True
#matplotlib.rcParams['axes.unicode_minus'] = False


import matplotlib.pyplot as plt
import matplotlib.animation as anim
from matplotlib.patches import Ellipse, Circle
from matplotlib import transforms
import matplotlib.ticker as ticker   

import argparse

######### DEFAULT OPTION VALUES ####################
#saveMovie = False
#saveFig = True
#useSatImg = True
#timeStepStart = 0 # max 7230
#########################################

parser = argparse.ArgumentParser(description="Victoria Park data / results animation")
parser.add_argument("-v", "--verbosity", help="increase output verbosity", action="store_true")
parser.add_argument("--satImg", help="use satellite image in plot", action="store_true")
parser.add_argument("--satImgLowRes", help="use the low resolution satellite image in the plot to reduce figure file size", action="store_true")
parser.add_argument("--allMeasure", help="If --saveFig is used, all measurements will be shown on the saved file. This will also force option --startTimeStep 0 to be set", action="store_true")
parser.add_argument("--noLandmarkEst", help="Do not show landmark estimates", action="store_true")
parser.add_argument("--noGPS", help="Do not show GPS trajectory", action="store_true")
parser.add_argument("--saveFig", help="save last animation frame as a pdf file", action="store_true")
parser.add_argument("--saveMovie", help="save animation as movie (mp4) file", action="store_true")
parser.add_argument("-s", "--startTimeStep", help="starting timestep (max 7230)", type=int, default=0, metavar='NUM')
parser.add_argument("--noEllipse", help="Uncertainty ellipses are not plotted", action="store_true")
parser.add_argument("resultDir", help="path to result data")
args = parser.parse_args()
if args.allMeasure:
    args.startTimeStep = 0
if args.startTimeStep > 7220:
    args.startTimeStep = 7220

nLandmarksDrawMax = 1000;
nMeasurementsDrawMax = 100;
nClutterDrawMax = 30

#sys.exit(0);

# Setting for file names

dataDir = args.resultDir
if dataDir[-1] != '/':
    dataDir += '/'

estPoseFile = 'particlePose.dat';
estPoseFile = dataDir + estPoseFile;
if os.path.exists(estPoseFile):
    print('Opening ' + estPoseFile);
else:
    print(estPoseFile + ' does not exist');
    sys.exit(0);
estPoseFileHandle = open(estPoseFile, "r");

bestTrajFile = 'trajectory.dat'
bestTrajFile = dataDir + bestTrajFile
if os.path.exists(bestTrajFile):
    print('Opening ' + bestTrajFile)
else:
    print(bestTrajFile + ' does not exist')
bestTrajFileHandle = open(bestTrajFile, "r")

estMapFile = 'landmarkEst.dat';
estMapFile = dataDir + estMapFile;
if os.path.exists(estMapFile):
    print('Opening ' + estMapFile);
else:
    print(estMapFile + ' does not exist')
    sys.exit(0);
estMapFileHandle = open(estMapFile, "r");

measurementFile = 'measurements.dat';
measurementFile = dataDir + measurementFile;
if os.path.exists(measurementFile):
    print('Opening ' + measurementFile);
else:
    print(measurementFile + ' does not exist')
    sys.exit(0);
measurementFileHandle = open(measurementFile, "r");

clutterFile = 'clutter.dat'
clutterFile = dataDir + clutterFile
clutterFileHandle = 0
if os.path.exists(clutterFile):
    print('Opening ' + clutterFile);
    clutterFileHandle = open(clutterFile, "r")

gpsFile = 'gps.dat'
gpsFile = dataDir + gpsFile
gpsFileHandle = 0
if os.path.exists(gpsFile):
    print('Opening ' + gpsFile)
    gpsFileHandle = open(gpsFile, "r")

estimateImageFile = 'estimate.pdf';
estimateImageFile = dataDir + estimateImageFile;
estimateMovieFile = 'estimate.mp4';
estimateMovieFile = dataDir + estimateMovieFile;

# Reading files

p = np.fromfile(estPoseFileHandle, dtype=float, count=6, sep=" ");
p_idx = 0;
p_idx_maxWeight = 0;
p_maxWeight = p[5];
p_x = [];
p_y = [];
p_r = [];
p_w = [];
p_t = p[0]
while p[0] == p_t:
    p_x.append(p[2]);
    p_y.append(p[3]);
    p_r.append(p[4]);
    p_w.append(p[5]);
    if p[5] > p_maxWeight :
        p_maxWeight = p[5];
        p_idx_maxWeight = p_idx;
    p = np.fromfile(estPoseFileHandle, dtype=float, count=6, sep=" ");
    p_idx += 1;
px_best = []
py_best = []
pr_best = []
px_best.append(p[2]);
py_best.append(p[3]);
pr_best.append(p[4]);
p_maxWeight = p[5];

bt = np.fromfile(bestTrajFileHandle, dtype=float, count=4, sep=" ")
while bt[0] < p[0]:
    bt = np.fromfile(bestTrajFileHandle, dtype=float, count=4, sep=" ")
bt_x = []
bt_y = []
bt_r = []

m = np.fromfile(estMapFileHandle, count=8, sep=" ", dtype=float);

z = np.fromfile(measurementFileHandle, count=4, sep=" ", dtype=float);
allZ_x = []
allZ_y = []

c = np.array([])
if(clutterFileHandle != 0):
    c = np.fromfile(clutterFileHandle, count=4, sep=" ", dtype=float)

g = np.array([])
if(gpsFileHandle != 0):
    g = np.fromfile(gpsFileHandle, count=3, sep=" ", dtype=float) # first line is all 0
    g = np.fromfile(gpsFileHandle, count=3, sep=" ", dtype=float)
gps_x = []
gps_y = []

# Plotting 

fig = plt.figure( figsize=(12,10),facecolor='w')

particles, = plt.plot(p_x, p_y, 'b.');
ax = plt.gca()
plt.axis('equal');
plt.grid(True);
plt.xlim([-125, 225])
plt.ylim([-50, 250])


if args.satImg and os.path.exists(dataDir + "VictoriaParkSatellite.png"):
    satImg = np.flipud(plt.imread(dataDir + "VictoriaParkSatellite.png"))
    satImg = ndimage.rotate(satImg, -7, reshape=False)
    img_x_min = -195
    img_y_min = -80
    img_x_max = img_x_min + 475
    img_y_max = img_y_min + 375
    plt.imshow(satImg, zorder=0, origin='lower', extent=[img_x_min, img_x_max, img_y_min, img_y_max])
elif args.satImgLowRes and os.path.exists(dataDir + "VictoriaParkSatellite_lowRes.png"):
    satImg = np.flipud(plt.imread(dataDir + "VictoriaParkSatellite_lowRes.png"))
    satImg = ndimage.rotate(satImg, -7, reshape=False)
    img_x_min = -195
    img_y_min = -80
    img_x_max = img_x_min + 475
    img_y_max = img_y_min + 375
    plt.imshow(satImg, zorder=0, origin='lower', extent=[img_x_min, img_x_max, img_y_min, img_y_max])
else:
    args.satImg = False;

measurements = [];
for i in range(0, nMeasurementsDrawMax) : 
    measurement_line, = plt.plot([], [], 'r-');
    measurements.append( measurement_line );
allZ, = plt.plot([], [], 'k.')
if( not args.satImg and not args.satImgLowRes ):
    allZ.set_color([0.2, 0.2 , 0.8])

clutter = []
for i in range(0, nClutterDrawMax) : 
    clutter_line, = plt.plot([],[], 'g-')
    clutter.append( clutter_line )

landmarks = [];
for i in range(0, nLandmarksDrawMax) : 
    landmark_ellipse = Ellipse(xy=(0,0), width=0, height=0, angle=0);
    landmarks.append(landmark_ellipse); 
    ax.add_patch(landmarks[i]);

landmarkCenters, = plt.plot([], [], '+')
landmarkCenters.set_color([1,0.6,0])
if( not args.satImg and not args.satImgLowRes ):
    landmarkCenters.set_color([0.2, 0.2 , 0.8])
trajectory, = plt.plot(0, 0, 'b-')
trajectory_best, = plt.plot(0, 0, 'y-')
if( not args.satImg and not args.satImgLowRes ):
    trajectory_best.set_color('k')
gps, = plt.plot(0, 0, 'r.')
gps_legend, = plt.plot(0, 0, 'r-')
plt.setp(gps, ms=0.75)


txt = plt.text(150, -40, " ");


timeStepCurrent = 0


def animateInit():

    txt.set_text("Time: ");
    particles.set_data([],[]);
    for i in range(0, nMeasurementsDrawMax) :
        measurements[i].set_data([],[]);
        measurements[i].set_color([1.0, 0.2 ,0.2]);
    for i in range(0, nClutterDrawMax):
        clutter[i].set_data([],[])
        clutter[i].set_color([0.2, 1, 0.2])
    for i in range(0, nLandmarksDrawMax):
        landmarks[i].center = (0,0);
        landmarks[i].width = 0;
        landmarks[i].height = 0;
        landmarks[i].set_facecolor([1,0.6,0])
        if( not args.satImg and not args.satImgLowRes):
            landmarks[i].set_facecolor([0.2, 0.2, 0.8])
    return [];

def animate(i):
    
    global p;
    global bt;
    global m;
    global z;
    global c
    global g
    global timeStepCurrent
    global p_maxWeight

    if not p.any():
        return []

    currentTime = p[0];
    drawnObjects = [];

    # Time
    txt.set_text("Time: {0}".format(currentTime));
    drawnObjects.append(txt);

    if args.verbosity and i % 100 == 0:
        print(i)

    if timeStepCurrent < i:
        while p.any() and timeStepCurrent < i:
            if currentTime != p[0]:                
                px_best.append(p[2]);
                py_best.append(p[3]);
                pr_best.append(p[4]);
                p_maxWeight = p[5];
                timeStepCurrent = timeStepCurrent + 1
                currentTime = p[0];
            elif p[5] > p_maxWeight :
                p_maxWeight = p[5];
                px_best[-1] = p[2];
                py_best[-1] = p[3];
                pr_best[-1] = p[4];
            p = np.fromfile(estPoseFileHandle, dtype=float, count=6, sep=" ")
                
        while bt.any() and bt[0] < currentTime:
            bt_x.append(bt[1])
            bt_y.append(bt[2])
            bt_r.append(bt[3])
            bt = np.fromfile(bestTrajFileHandle, dtype=float, count=4, sep=" ")
        while m.any() and m[0] < currentTime:
            m = np.fromfile(estMapFileHandle, count=8, sep=" ", dtype=float);
        while z.any() and z[0] < currentTime:
            z = np.fromfile(measurementFileHandle, count=4, sep=" ", dtype=float);
        while c.any() and c[0] < currentTime:
            c = np.fromfile(clutterFileHandle, count=4, sep=" ", dtype=float);
        while g.any() and g[0] < currentTime:
            gps_x.append(g[1])
            gps_y.append(g[2])
            g = np.fromfile(gpsFileHandle, count=3, sep=" ", dtype=float)
            

    # Particles
    p_idx = 0;
    p_idx_maxWeight = 0;
    p_maxWeight = p[5];
    px_best.append(p[2]);
    py_best.append(p[3]);
    pr_best.append(p[4]);
    p_x = [];
    p_y = [];
    p_w = [];
    while p.any() and abs(p[0] - currentTime) < 1e-4:
        p_x.append(p[2]);
        p_y.append(p[3]);
        p_r.append(p[4]);
        p_w.append(p[5]);
        if p[5] > p_maxWeight :
            p_maxWeight = p[5];
            p_idx_maxWeight = p_idx;
            px_best[i] = p[2];
            py_best[i] = p[3];
            pr_best[i] = p[4];
        p = np.fromfile(estPoseFileHandle, dtype=float, count=6, sep=" ");
        p_idx += 1;
    particles.set_data(p_x, p_y)
    trajectory.set_data(px_best, py_best)
    timeStepCurrent = timeStepCurrent + 1

    # Best trajectory
    while bt[0] < currentTime:
        bt_x.append(bt[1])
        bt_y.append(bt[2])
        bt_r.append(bt[3])
        bt = np.fromfile(bestTrajFileHandle, dtype=float, count=4, sep=" ")
    trajectory_best.set_data(bt_x, bt_y)
    

    # Landmarks
    m_idx = 0;
    m_x = []
    m_y = []

    while m.any() and abs(m[0] - currentTime) < 1e-12:

        if not args.noLandmarkEst:

            m_x.append(m[2])
            m_y.append(m[3])

            if not args.noEllipse:
                cov = np.array([ [ m[4], m[5] ], [ m[5], m[6] ] ]);
                w = m[7];
                eVal, eVec = np.linalg.eig(cov);
                eVal = eVal.real;
                a1 = 4*np.sqrt(eVal[0]); # Assume this is semi-major axis first
                a2 = 4*np.sqrt(eVal[1]); # 3 dof, 4 stdev is roughly prob = 0.997
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

    landmarkCenters.set_data(m_x, m_y)
    drawnObjects.append(landmarkCenters)

    # Measurements
    nZ = 0;
    while z.any() and abs(z[0] -  currentTime) < 1e-12:
        z_dir = pr_best[i] + z[2] - np.pi / 2;
        z_end = [px_best[i] + z[1]*np.cos(z_dir), py_best[i] + z[1]*np.sin(z_dir) ];
        measurements[nZ].set_data([px_best[i], z_end[0]], [py_best[i], z_end[1]]);
        if args.allMeasure:
            allZ_x.append(z_end[0])
            allZ_y.append(z_end[1])
        drawnObjects.append(measurements[nZ]);
        z = np.fromfile(measurementFileHandle, count=4, sep=" ", dtype=float);
        nZ += 1;
    while measurements[nZ].get_xdata() != []:
        measurements[nZ].set_data([], []);
        nZ += 1;
    nZ = 0
    while c.any() and abs(c[0] -  currentTime) < 1e-12:
        c_dir = pr_best[i] + c[2] - np.pi / 2;
        c_end = [px_best[i] + c[1]*np.cos(c_dir), py_best[i] + c[1]*np.sin(c_dir) ];
        clutter[nZ].set_data([px_best[i], c_end[0]], [py_best[i], c_end[1]]);
        if args.allMeasure:
            allZ_x.append(c_end[0])
            allZ_y.append(c_end[1])
        drawnObjects.append(clutter[nZ]);
        c = np.fromfile(clutterFileHandle, count=4, sep=" ", dtype=float);
        nZ += 1;
    while clutter[nZ].get_xdata() != []:
        clutter[nZ].set_data([], []);
        nZ += 1;


    # GPS
    while g.any() and g[0] < currentTime:
        if not args.noGPS:
            gps_x.append(g[1])
            gps_y.append(g[2])
        g = np.fromfile(gpsFileHandle, count=3, sep=" ", dtype=float)
    if not args.noGPS:
        gps.set_data(gps_x, gps_y)
    
    drawnObjects.append(gps)
    drawnObjects.append(particles)
    #drawnObjects.append(trajectory)
    drawnObjects.append(trajectory_best)

    return drawnObjects;

animation = anim.FuncAnimation(plt.figure(1), animate, np.arange(args.startTimeStep, 7230), interval=1, init_func=animateInit, blit=True, repeat=False);

if args.saveMovie:
    FFMpegWriter = matplotlib.animation.writers['ffmpeg']
    animation.save(estimateMovieFile, writer=FFMpegWriter(fps=30, codec='mpeg4', bitrate=300000))
    # run ffmpeg -codecs on shell to see available video encoding codecs
else:
    plt.show(block=False)

if args.saveFig:
    
    for i in range(0, nMeasurementsDrawMax) : 
        measurements[i].remove();
    for i in range(0, nClutterDrawMax) :
        clutter[i].remove()
    plt.setp(gps, ms=0.75)
    plt.setp(gps_legend, linewidth=1)
    txt.set_text(" ");
    trajectory.set_data([], [])
    particles.set_data([], [])
    trajectory_best.set_data(bt_x[0::20], bt_y[0::20])
    
    plt.setp(landmarkCenters, markersize=4)
    matplotlib.rcParams.update({'font.size': 6})
    plt.gcf().set_size_inches(4, 3.3)

    legend_list = [trajectory_best]
    legend_name = ["Trajectory Est."]

    if not args.noLandmarkEst:
        legend_list.append(landmarkCenters)
        legend_name.append("Landmark Position Est.")
        if not args.noEllipse:
            legend_list.append(landmarks[0])
            legend_name.append("Landmark Est. Uncertainty")

    if args.allMeasure:
        allZ.set_data(allZ_x, allZ_y)
        plt.setp(allZ, markersize=2)
        legend_list.append(allZ)
        legend_name.append("Measurements")

    if not args.noGPS:
        legend_list.append(gps_legend)
        legend_name.append("GPS Vehicle Trajectory")

    plt.legend(legend_list, legend_name, loc='upper right', labelspacing=0.5);

    plt.setp(plt.gca().get_legend().get_texts(), fontsize='6')

    scale = 1;
    ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*scale))
    plt.gca().xaxis.set_major_formatter(ticks)
    ticks = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y*scale))
    plt.gca().yaxis.set_major_formatter(ticks)
    
    plt.savefig(estimateImageFile, format='pdf', bbox_inches='tight')

measurementFileHandle.close();


