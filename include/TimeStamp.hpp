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

#ifndef TIMESTAMP_HPP
#define TIMESTAMP_HPP

#include <stdint.h>

namespace rfs{

/**
 * \class TimeStamp
 * \brief A timestamp
 */
class TimeStamp{

public:

  /** Default constructor */
  TimeStamp();

  /** Constructor 
   *  \param[in] s number of seconds
   *  \param[in] ns number of nanoseconds
   */
  TimeStamp(int32_t s, int32_t ns);

  /** Constructor 
   *  \param[in] t time in seconds
   */
  TimeStamp(double t);

  /** Destructor */
  ~TimeStamp();

  /** Comparison operator == */
  bool operator== (const TimeStamp &other);

  /** Binary arithmetic operator += */ 
  TimeStamp& operator+= (const TimeStamp &rhs);

  /** Binary arithmetic operator + */
  TimeStamp operator+ (const TimeStamp &rhs) const;

  /** Binary arithmetic operator -= */ 
  TimeStamp& operator-= (const TimeStamp &rhs);

  /** Binary arithmetic operator - */
  TimeStamp operator- (const TimeStamp &rhs) const;

  /** Get the time expressed as a double 
   *  \return time
   */
  double getTimeAsDouble() const;

  /** Set the time
   *  \param[in] t time in seconds
   */
  void setTime(double const t);

  /** Set the time
   *  \param[in] s seconds
   *  \param[in] ns nanoseconds
   */
  void setTime(int32_t const s, int32_t const ns);

  int32_t sec; /**< number of seconds */
  int32_t nsec; /**< number of nanoseconds */

private:

  /** Normalize so that nsec < 1e9 */
  void normalize();

};

// Inline functions

inline bool TimeStamp::operator== (const TimeStamp &other){
  if(this->sec == other.sec && this->nsec == other.nsec){
    return true;
  }
  return false;
}

inline TimeStamp& TimeStamp::operator+= (const TimeStamp &rhs){
  sec += rhs.sec;
  nsec += rhs.nsec;
  normalize();
  return *this;
}

inline TimeStamp TimeStamp::operator+ (const TimeStamp &rhs) const{

  TimeStamp lhs(this->sec + rhs.sec, this->nsec + rhs.nsec);
  lhs.normalize();
  return lhs;
}

inline TimeStamp& TimeStamp::operator-= (const TimeStamp &rhs){
  sec -= rhs.sec;
  nsec -= rhs.nsec;
  normalize();
  return *this;
}

inline TimeStamp TimeStamp::operator- (const TimeStamp &rhs) const{
  
  TimeStamp lhs(this->sec - rhs.sec, this->nsec - rhs.nsec);
  lhs.normalize();
  return lhs;
}

inline void TimeStamp::normalize(){
  
  int32_t nsec_temp = nsec % 1000000000L;
  int32_t sec_temp = sec + nsec / 1000000000L;
  if(nsec_temp < 0){
    nsec_temp += 1000000000L;
    sec_temp--;
  }
  nsec = nsec_temp;
  sec = sec_temp;
}

}
#endif
