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

The main purpose of this RB-PHD Filter library is promote research in random
finite set (RFS) estimation methods for the problem of simultaneous localization
and mapping (SLAM). The intention of the authors, Keith and Felipe, is to keep
the library general in the sense that users can define system models that are
relevant to their specific problem, while not having to program and test their
own implementation of the RB-PHD filter. We have also provided some generic vehicle
motion models and measurement models. This library is an on-going project, 
and we intend to update it with new improvements whenever possible. Of course, 
any feedback will be appreciated. 

While we have completed the documentation for all the individual classes and
functions in this library, we are still lacking a document to explain how 
everything fits together (although we plan to have this ready in the near
future). For those who wants to use the current version of the library, the
2d simulation that is included will provide some guidance of the necessary steps
to get things to work. We have also included some Matlab scripts for visualizing
the results from the simulator. 

Please note that we have only tested the library in Ubuntu 13.04.

To get started, make sure that you have the following libraries installed:

Boost::math_c99
Boost::timer
Boost::system
Boost::thread
Eigen
libconfig (This is only necessary for the simulator)

Compiling for parts of the library and the 2d Simulator can be done with cmake.
We have provided a CMakeLists.txt already. To compile, run:

cmake .
make

Documentations can be generated using Doxygen, and we have provided a Doxyfile
(configuration file). To generate the html and pdf documentations, run:

doxygen

The executable for the 2d simulator should be located in bin/sim2d. There is 
a configuration file for the simulator in cfg/simulator2d.cfg, where many
settings such as the length of the simulation, the number of landmarks, 
measurement noise, etc, can be modified. To start the simulator, run:

bin/sim2d <NUM>

where <NUM> is a seed for the random number generators that are used. Logging
of the results can be turned on in simulation2d.cfg. Under the tools directory,
there is a Matlab script for viewing the results. 





