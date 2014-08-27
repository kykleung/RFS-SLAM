/*
 * Software License Agreement (New BSD License)
 *
 * Copyright (c) 2014, Keith Leung, Felipe Inostroza
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

The main purpose of this RB-PHD Filter library is promote research in random
finite set (RFS) estimation methods for the problem of simultaneous localization
and mapping (SLAM). The intention of the authors is to keep the library general 
in the sense that users can define system models that are relevant to their 
specific problem, while not having to program and test their own implementation 
of RFS filters. This library is an on-going project, and we intend to update it 
when any related work is published. Any feedback will be appreciated. 

============================ VERSION HISTORY ==============================

1.0 - Initial release

1.1
  - RB-PHD SLAM algorithm updates:
    - included multi-feature particle importance weighting strategy
    - included single-cluster (SC)-PHD SLAM weighting strategy
  - new FastSLAM and MH-FastSLAM algorithms included 
  - updated 2-D simulations
  - new Optimal Sub-Pattern Assignment (OSPA) error metric class 
  - visualization tools are now in Python instead of Matlab
  - introduced namespace rfs to the library
  - some updates made to naming convention of classes in the library
  - cmake now generates rfsslam-config.cmake to enable find(rfsslam)
    from other projects

============================= INSTALLATION =================================

This library has been tested under Ubuntu 13.04, 13.10, and 14.04.

Prerequisite Libraries:
  Boost (version 1.53 minimum) with components:
    math_c99, timer, system, thread, filesystem, graph 
  Eigen (version 3.0.0 minimum)

Optional Libraries:
  libconfig (for the 2-D SLAM simulators)
  python-numpy (for visualizing 2-D SLAM results)
  python-matplotlib (for visualizing 2-D SLAM results)
    
Compile the library and the 2d Simulator using cmake. Run:

  cmake .
  make
  make install (option, and produces install_manifest.txt)
     
============================== DOCUMENTATION =================================

Documentations can be generated using Doxygen, and we have provided a Doxyfile
(configuration file). To generate the html and pdf documentations in the doc
directory, run:

  doxygen

============================= 2-D SLAM EXAMPLES ===============================

* SLAM Filters

Rao-Blackwellized Probability Hypothesis Density (RB-PHD) SLAM 2-D Simulation
  Source: src/rbphdslam2dSim.cpp
  Config: cfg/rbphdslam2dSim.cfg
  Run: bin/rbphdslam2dSim [random_seed=0] [cfg_file=cfg/rbphdslam2dSim.cfg]

Factored Solution to SLAM (FastSLAM 1.0) 2-D Simulation
  Source: src/fastslam2dSim.cpp 
  Config: cfg/fastslam2dSim.cfg
  Run: bin/fastslam2dSim [random_seed=0] [cfg_file=cfg/fastslam2dSim.cfg]

Multi-Hypothesis Factored Solution to SLAM (MH FastSLAM) 2-D Simulation
  Source: src/fastslam2dSim.cpp 
  Config: cfg/mhfastslam2dSim.cfg
  Run: bin/fastslam2dSim [random_seed=0] [cfg_file=cfg/mhfastslam2dSim.cfg]

* Visualization Tools

Animation
  Source: bin/animate2dSim.py
  Config: edit bin/animate2dSim.py directly
  Run: bin/animate2dSim.py [results_dir]


================================ FUTURE WORK ===================================

- Remove dependency on libconfig, and switch to boost::program_option
- Multi-threaded versions of SLAM algorithms
- Cardinalized Probability Hypothesis Density (CPHD) filter
- Multi-Bernoulli filter




