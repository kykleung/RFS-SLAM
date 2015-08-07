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

# Necessary to generate Type 1 fonts for pdf figures as required by IEEE for paper submissions
#matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['axes.unicode_minus'] = False
#matplotlib.rcParams['text.usetex'] = True
#plt.rc('text', usetex=True)

#matplotlib.rcParams.update({'font.size': 14})

import matplotlib.pyplot as plt
import matplotlib.animation as anim
from matplotlib.patches import Ellipse, Circle
from matplotlib import transforms
import matplotlib.ticker as ticker   
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

import argparse

parser = argparse.ArgumentParser(description="Batch Simulation Results Plotter")
parser.add_argument("-v", "--verbosity", help="increase output verbosity", action="store_true")
parser.add_argument("--saveFig", help="save last animation frame as a pdf file", action="store_true")
parser.add_argument("resultFile", help="path to result log file")
args = parser.parse_args()

# Open and read result log
if not os.path.exists(args.resultFile):
    print(args.resultFile + ' does not exist');
    sys.exit(0);
print('Reading ' + args.resultFile);
data = np.genfromtxt(args.resultFile)
np.set_printoptions(threshold=np.nan)
data = data[ np.lexsort((data[:,1], data[:,0])) ]

data_stat = np.empty((0,6), float)
Pd_current = data[0,0]
c_current = data[0,1]
trajErrors = [data[0,2]]
mapColas = [data[0,3]]
c_unique = []
for i in range(1, data.shape[0]):

    if data[i,0] == Pd_current and data[i,1] == c_current :
        trajErrors.append(data[i,2])
        mapColas.append(data[i,3])
    else:
        # Calculate statistics
        if(Pd_current > 0.3):
            data_add = np.array([Pd_current, np.log10(c_current)-2, np.mean(trajErrors), np.std(trajErrors), np.mean(mapColas), np.std(mapColas) ])
            data_stat = np.vstack( [data_stat, data_add])
            c_unique.append(np.log10(c_current)-2)

        #print trajErrors

        Pd_current = data[i,0]
        c_current = data[i,1]
        trajErrors = [data[i,2]]
        mapColas = [data[i,3]]

if(Pd_current > 0.3):
    data_add = np.array([Pd_current, np.log10(c_current)-2, np.mean(trajErrors), np.std(trajErrors), np.mean(mapColas), np.std(mapColas) ])
    data_stat = np.vstack( [data_stat, data_add])

    data_stat = data_stat[ np.lexsort((data_stat[:,0], data_stat[:,1])) ]

c_unique = np.unique(c_unique)

fig = plt.figure(1)
ax = fig.gca(projection='3d')
X = data_stat[:,0]
Y = data_stat[:,1]
Z = data_stat[:,2]
ax.plot_trisurf(X, Y, Z, cmap=cm.coolwarm, linewidth=0, vmin=0, vmax=1.5)
ax.set_zlim3d(0, 1.5)
ax.set_xlabel(r"Prob. of detection")
ax.set_ylabel(r"log clutter intensity")
ax.set_zlabel(r"Avg. robot position error [m]")
if(args.saveFig):
    plt.savefig("batch_robot_error.pdf", format='pdf', bbox_inches='tight')

fig = plt.figure(2)
ax = fig.gca(projection='3d')
X = data_stat[:,0]
Y = data_stat[:,1]
Z = data_stat[:,4]
ax.plot_trisurf(X, Y, Z, cmap=cm.coolwarm, linewidth=0, vmin=0, vmax=50)
ax.set_zlim3d(0, 50)
ax.set_xlabel("Prob. of detection")
ax.set_ylabel("log clutter intensity")
ax.set_zlabel(r"Avg. landmark position error [m]")
if(args.saveFig):
    plt.savefig("batch_map_error.pdf", format='pdf', bbox_inches='tight')

if( not args.saveFig ):
    plt.show()
