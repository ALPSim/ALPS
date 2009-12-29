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
	  	char cm[3];
		char cs[3];
		char cd[3];

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
