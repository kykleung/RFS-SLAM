/*
 * Software License Agreement (New BSD License)
 *
 * Copyright (c) 2013, Keith Leung
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

#ifndef TIMER_HPP
#define TIMER_HPP

#include<boost/timer/timer.hpp> 

namespace rfs{

/** 
 * \class Timer
 * \brief A wrapper for the boost cpu timer
 */
class Timer
{
public:

  /** Default constructor */
  Timer();
  /** Default destructor */
  ~Timer();
  /** Reset and start the timer */
  void start();
  /** Stop the timer */
  void stop();
  /** Resume the timer after stopping */
  void resume();

  /** 
   * Get the elapsed time from the timer since the start (as a string)
   * \param[out] t_wall wall time in nsec
   * \param[out] t_cpu cpu (system + user) time in nsec
   */
  void elapsed(std::string &t_wall, std::string &t_cpu);

  /** 
   * Get the elapsed time from the timer since the start
   * \param[out] t_wall wall time in nsec
   * \param[out] t_cpu cpu (system + user) time in nsec
   */
  void elapsed(long long &t_wall, long long &t_cpu);

private:

  ::boost::timer::cpu_timer timer_;

};

}

#endif
