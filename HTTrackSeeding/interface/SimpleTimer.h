#ifndef HTTrackSeeding_SimpleTimer_H
#define HTTrackSeeding_SimpleTimer_H

/*** \class  SimpleTimer
  *
  *  A simple timer. Just for now, it will soon be removed.
  *
  *  \author Maurizio Lo Vetere
  */

#include <sys/time.h>


class SimpleTimer {
  public:
    SimpleTimer      ( );
    double     lapse ( )  const;
    void       reset ( );
    void       start ( );
    void       stop  ( );
  private:
    timeval _start;
    double  _timer;
};


inline SimpleTimer::SimpleTimer()
{
  _timer = 0;
  gettimeofday( &_start, 0 ); 
}


inline double  SimpleTimer::lapse()  const
{ 
  return _timer; 
}  


inline void  SimpleTimer::reset() 
{
  _timer = 0; 
}

inline void  SimpleTimer::start()
{
  gettimeofday( &_start, 0 );
}


inline void  SimpleTimer::stop()
{
  timeval now;
  gettimeofday( &now , 0 );  
  double  time1 = ((double)(_start.tv_sec) + (double)(_start.tv_usec)/1000000.);
  double  time2 = ((double)(   now.tv_sec) + (double)(   now.tv_usec)/1000000.);
  _timer += (time2 - time1);
}


#endif // HTTrackSeeding_SimpleTimer_H
