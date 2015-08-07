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
* Compiles with llvm on Mac (OSX 10.10)

### Installation ###

#### Source ####

Obtain from git repository: `https://kykleung@bitbucket.org/kykleung/phdfilter.git`

#### C++ Library Dependencies ####

* Boost (version 1.53 minimum) with components:
    * timer
    * chrono
    * system
    * filesystem
    * graph
    * program_options
* Eigen (version 3.0.0 minimum)
* gtest (optional)
  
    If using Ubuntu apt-get to install: 

    1. `sudo apt-get install libgtest-dev`
    2. `cd /usr/src/gtest`
    3. `cmake .`
    4. `make`
    5. `sudo mv libgtest_main.a /usr/lib/`
    6. `sudo mv libgtest.a /usr/lib/`

#### Other Dependencies ####

For visualizing 2-D SLAM results

* python-numpy (for visualizing 2-D SLAM results)
* python-matplotlib (for visualizing 2-D SLAM results)

#### Compiling ####

For out-of-source build of the library and the 2D simulator:

* `mkdir build`
* `cd build`
* `cmake ..`  or `ccmake ..`
* `make`
* `make install` (optional, and produces install_manifest.txt)

#### Documentation ####

Documentations can be generated using *Doxygen*, and we have provided a Doxyfile
(configuration file). To generate the html and pdf documentations **Doxygen needs
to be installed**. Run `doxygen` in project root directory.
Documentation will be generated in `doc/` directory.

### Usage Example ###

#### SLAM Filters ####

* Rao-Blackwellized Probability Hypothesis Density (RB-PHD) SLAM 2-D Simulation
    - Source: src/rbphdslam2dSim.cpp
    - Config: cfg/rbphdslam2dSim.xml
    - Run: `bin/rbphdslam2dSim`

* Factored Solution to SLAM (FastSLAM 1.0) 2-D Simulation
    - Source: src/fastslam2dSim.cpp 
    - Config: cfg/fastslam2dSim.xml
    - Run: `bin/fastslam2dSim`

* Multi-Hypothesis Factored Solution to SLAM (MH FastSLAM) 2-D Simulation
    - Source: src/fastslam2dSim.cpp 
    - Config: cfg/mhfastslam2dSim.xml
    - Run: `bin/fastslam2dSim`

* RB-PHD SLAM on the Victoria Park dataset
    - Source: src/rbphdslam_VictoriaPark.cpp
    - Config: cfg/rbphdslam_VictoriaPark.xml and cfg/rbphdslam_VictoriaPark_artificialClutter.xml
    - Run: `bin/rbphdslam_VictoriaPark`

* FastSLAM on the Victoria Park dataset
    - Source: src/fastslam_VictoriaPark.cpp
    - Config: cfg/fastslam_VictoriaPark.xml and cfg/fastslam_VictoriaPark_artificialClutter.xml
    - Run: `bin/fastslam_VictoriaPark`

* MH FastSLAM on the Victoria Park dataset
    - Source: src/mhfastslam_VictoriaPark.cpp
    - Config: cfg/mhfastslam_VictoriaPark.xml and cfg/mhfastslam_VictoriaPark_artificialClutter.xml
    - Run: `bin/fastslam_VictoriaPark`

#### Analysis Tools ####

For calculating errors for 2d simulations, run: `bin/analysis2dSim [results_dir]`.

For plotting the errors after running the analysis executable, use: 
    - `scripts/sim/plotError2dSim.py`
    - `scripts/sim/plotErrorCompare2dSim.py`

#### Visualization Tools ####

For animating 2D SLAM simulation results, run: `scripts/sim/animate2dSim.py [results_dir]`. 
Edit the python script and set `saveMoive=False` to see animation.
Set `saveMoive=True` to generate a mp4 file.

For animating Victoria Park dataset results, run: `scripts/VictoriaPark/animate_VictoriaPark.py [results_dir]`.
Use `-h` or `--help` to see options.

#### Performance Profiling Tools ####

Performance profiling is currently available for:
  - `bin/rbphdslam2dSim`
  - `bin/rbphdslam_VictoriaPark`
  - `bin/fastslam2dSim`
  - `bin/fastslam_VictoriaPark`

Use `ccmake` to turn on `USE_CPU_PROFILER` and or `USE_HEAP_PROFILER`.
Performance profiles are recorded in `.prof` files in the current directory.
Use `google-pprof` to parse the profiles.
At the moment, profiling does not provide meaningful results on OS X machines due to Address space layout randomization (ASLR).

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

* Latest (1.2.0)
    - Implementation of joint compatibility branch and bound (JCBB) 
      for data association in vector-based methods
    - Config files now use xml format to removed dependency on the libconfig library
    - OSX compatible when compiling with Clang/LLVM
    - Multi-threaded versions of SLAM algorithms using OpenMP
        - multithreading with OpenMP is currently not supported by Clang/LLVM 
        - An OpenMP-supported LLVM compiler is available at: [http://clang-omp.github.io/](http://clang-omp.github.io/)	
    - Inclusion of the Victoria Park dataset and the code for processing it using various SLAM filters.
    - Performance profiling option using Google Perftools
    - Updates to CMakeLists.txt to make options more operating system specific
    - Executables now use the Boost program_options library to handle arguments

### Future Work ###

- Cardinalized Probability Hypothesis Density (CPHD) filter
- Cardinality Balanced Multi-Bernoulli filter
- Reimplementation of the data stuctures for Gaussian mixtures (map/landmarks)

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




