README
===============

### Summary ###

The main purpose of this library is promote research in random
finite set (RFS) estimation methods for the problem of simultaneous localization
and mapping (SLAM). The intention of the authors is to keep the library general 
in the sense that users can define system models that are relevant to their 
specific problem, while not having to program and test their own implementation 
of RFS filters. This library is an on-going project, and we intend to update it 
when any related work is published. Any feedback will be appreciated. 

* License: New BSD
* Version: 1.1.0
* Compiles with gcc on Linux (Ubuntu 13.04, 13.10, 14.04)

### Installation ###

#### Source ####

Obtain from git repository: `https://github.com/kykleung/RFS-SLAM.git`

#### C++ Library Dependencies ####

* Boost (version 1.53 minimum) with components:
    * math_c99
    * timer
    * system
    * thread
    * filesystem
    * graph 
* Eigen (version 3.0.0 minimum)
* gtest 
  
    If using Ubuntu apt-get to install: 

    1. `sudo apt-get install libgtest-dev`
    2. `cd /usr/src/gtest`
    3. `cmake .`
    4. `make`
    5. `sudo mv libgtest_main.a /usr/lib/`
    6. `sudo mv libgtest.a /usr/lib/`

* libconfig (optional, for the 2-D SLAM simulators )

#### Other Dependencies ####

For visualizing 2-D SLAM results

* python-numpy (for visualizing 2-D SLAM results)
* python-matplotlib (for visualizing 2-D SLAM results)

#### Compiling ####

For out-of-source build of the library and the 2D simulator:

* `mkdir build`
* `cd build`
* `cmake ..`
* `make`
* `make install` (option, and produces install_manifest.txt)

#### Documentation ####

Documentations can be generated using *Doxygen*, and we have provided a Doxyfile
(configuration file). To generate the html and pdf documentations **Doxygen needs
to be installed**. Run `doxygen` in project root directory.
Documentation will be generated in `doc/` directory.

### Usage Example ###

#### SLAM Filters ####

* Rao-Blackwellized Probability Hypothesis Density (RB-PHD) SLAM 2-D Simulation
    - Source: src/rbphdslam2dSim.cpp
    - Config: cfg/rbphdslam2dSim.cfg
    - Run: `bin/rbphdslam2dSim [random_seed=0] [cfg_file=cfg/rbphdslam2dSim.cfg]`

* Factored Solution to SLAM (FastSLAM 1.0) 2-D Simulation
    - Source: src/fastslam2dSim.cpp 
    - Config: cfg/fastslam2dSim.cfg
    - Run: `bin/fastslam2dSim [random_seed=0] [cfg_file=cfg/fastslam2dSim.cfg]`

* Multi-Hypothesis Factored Solution to SLAM (MH FastSLAM) 2-D Simulation
    - Source: src/fastslam2dSim.cpp 
    - Config: cfg/mhfastslam2dSim.cfg
    - Run: `bin/fastslam2dSim [random_seed=0] [cfg_file=cfg/mhfastslam2dSim.cfg]`

#### Visualization Tools ####

For animating 2D SLAM results, run: `bin/animate2dSim.py [results_dir]`. 
Edit the python script and set `saveMoive=False` to see animation.
Set `saveMoive=True` to generate a mp4 file.

### Version History ###

* 1.0.0 
    - Initial release

* 1.1.0
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

### Future Work ###

- Remove dependency on libconfig, and switch to boost::program_option
- Multi-threaded versions of SLAM algorithms
- Cardinalized Probability Hypothesis Density (CPHD) filter
- Multi-Bernoulli filter

### Contact ###

* Maintainers:
    * Keith Leung (kykleung[at]gmail.com)
    * Felipe Inostroza

* Contributors:
    * Daniel Luhr


### License ###

 >
 > Software License Agreement (New BSD License)
 >
 > Copyright (c) 2014, Keith Leung, Felipe Inostroza
 > All rights reserved.
 > 
 > Redistribution and use in source and binary forms, with or without
 > modification, are permitted provided that the following conditions are met:
 >     * Redistributions of source code must retain the above copyright
 >       notice, this list of conditions and the following disclaimer.
 >     * Redistributions in binary form must reproduce the above copyright
 >       notice, this list of conditions and the following disclaimer in the
 >       documentation and/or other materials provided with the distribution.
 >     * Neither the name of the Advanced Mining Technology Center (AMTC), the
 >       Universidad de Chile, nor the names of its contributors may be 
 >       used to endorse or promote products derived from this software without 
 >       specific prior written permission.
 > 
 > THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 > ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 > WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 > DISCLAIMED. IN NO EVENT SHALL THE AMTC, UNIVERSIDAD DE CHILE, OR THE COPYRIGHT 
 > HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 > CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE 
 > GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
 > HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
 > LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
 > THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 >
