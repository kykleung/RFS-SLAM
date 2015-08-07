#!/usr/bin/python

import sys
import os.path
import matplotlib
import matplotlib.pyplot as plt
import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Map Octree Visualizer")
parser.add_argument("dataFile", help="path to tree data")
args = parser.parse_args()

if os.path.exists(args.dataFile):
    print('Opening ' + args.dataFile)
else:
    print(args.dataFile + ' does not exist');
    sys.exit(0);
treeData = np.genfromtxt(args.dataFile)

pt = np.empty([0,2])
box = np.empty([0,5])
for i in range(0, treeData.shape[0]):
    if( treeData[i][1] == treeData[i][3] and treeData[i][2] == treeData[i][4] ): # point
        pt = np.vstack( [pt, treeData[i,1:3]] )
    else: #box
        box = np.vstack( [box, treeData[i,1:6]] )

print pt
print box

fig = plt.figure()
ax1 = fig.add_subplot(111, aspect='equal')
for i in range (0, box.shape[0]):
    ax1.add_patch( matplotlib.patches.Rectangle( box[i,:], box[i,2]-box[i,0], box[i,3]-box[i,1], fill=False ) )
    
plt.plot(pt[:,0], pt[:,1], marker='+', color='r', linewidth=0)

plt.show()
