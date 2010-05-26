/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

//  ctimer.h
//  timer class


#ifndef __DMTK_CTIMER_H__
#define __DMTK_CTIMER_H__

#include <boost/date_time/posix_time/posix_time.hpp>
//#include <sys/resource.h>
#include <stdio.h>
#include <string>
#include <sstream>

namespace dmtk
{

enum
{
  CTIMER_HOURS,
  CTIMER_SECONDS,
};

typedef double cpu_time;

void mytime(cpu_time *my_timer) 
{ 
  static boost::posix_time::ptime start = boost::posix_time::second_clock::local_time();
        *my_timer = (boost::posix_time::second_clock::local_time()-start).seconds(); 
}

class CTimer
{
private:
        cpu_time _start;
        cpu_time _lap;
        size_t _format;


public:
        //  constructor
        CTimer (): _start(0), _lap(0), _format(CTIMER_HOURS) {}

        inline CTimer& Start () { mytime(&_start); _lap = _start; return *this; }
        inline CTimer& Lap () { mytime(&_lap); return *this; }

        CTimer& SetFormat(size_t format) 
          { _format = format; return *this; }
        size_t Format() const { return _format; }

        std::string LapTime () const
        {
                cpu_time lap;
                mytime(&lap);
                return Hours(lap-_lap, _format);
        }

        std::string TotalTime () const
        {
                cpu_time lap;
                mytime(&lap);
                return Hours(lap-_start, _format);
        }

        std::string Hours (const cpu_time& t, size_t f = CTIMER_HOURS) const
        {
                char s[200];
                int hh, mm, ss, dd;
                  char cm[30];
                char cs[30];
                char cd[30];

                switch(f){
                  case CTIMER_HOURS:
                    hh = int(t) / 3600;
                        mm = int(t - hh * 3600) / 60;
                    ss = int(t - hh*3600 - mm*60);
                    dd = int((t - hh*3600 - mm*60 - ss)*100);
                    
                    if(ss < 10) 
                      sprintf(cs,"0%i",ss);
                    else
                      sprintf(cs,"%i",ss);
                    if(mm < 10) 
                      sprintf(cm,"0%i",mm);
                    else
                      sprintf(cm,"%i",mm);
                    if(dd < 10) 
                      sprintf(cd,"0%i",dd);
                    else
                      sprintf(cd,"%i",dd);
                    sprintf(s,"%i:%s:%s.%s",hh,cm,cs,cd);
                    break;
                  case CTIMER_SECONDS:
                  default:
                    sprintf(s,"%f sec.",t);
                    break;
                }

                return std::string(s);
        }

};
        
} // namespace dmtk

#endif // __DMTK_CTIMER_H__
