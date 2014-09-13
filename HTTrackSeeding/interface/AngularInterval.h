#ifndef HTTrackSeeding_AngularInterval_H
#define HTTrackSeeding_AngularInterval_H

/*** \class  AngularInterval
  *
  *  Only sets of non null meausure are relevant in this context therefore we consider open intervals  
  *
  *  \author Maurizio Lo Vetere
  */

#include "DataFormats/GeometryVector/interface/Phi.h"
#include "MLoVetere/HTTrackSeeding/interface/Interval.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>
#include <vector>


template <typename T>
class AngularIntervalBase {
  public:
    AngularIntervalBase  ( );
    AngularIntervalBase  ( T inf, T sup );
    T                         center               ( )                               const;
    IntervalBase<T>           cosRange             ( )                               const;
    bool                      include              ( T angle )                       const;
    bool                      isEmpty              ( )                               const;
    bool                      isFull               ( )                               const;
    T                         length               ( )                               const;
    T                         lower                ( )                               const;
    T                         normalize            ( T angle )                       const;
    std::pair<IntervalBase<T>,IntervalBase<T>>               
                              normalize            ( AngularIntervalBase<T> other )  const;
    AngularIntervalBase<T>    operator+            ( AngularIntervalBase<T> other )  const;
    AngularIntervalBase<T>    operator-            ( AngularIntervalBase<T> other )  const;
    AngularIntervalBase<T>    operator-            ( )                               const;
    bool                      overlaps             ( AngularIntervalBase<T> other )  const;
    IntervalBase<T>           sinRange             ( )                               const;
    AngularIntervalBase<T>    smallestCovering     ( AngularIntervalBase<T> other )  const;
    T                         upper                ( )                               const;
    AngularIntervalBase<T> &  turnClockWise        ( T angle );
    AngularIntervalBase<T> &  turnClockWise        ( AngularIntervalBase<T> other );
    AngularIntervalBase<T> &  turnCounterClockWise ( T angle );
    AngularIntervalBase<T> &  turnCounterClockWise ( AngularIntervalBase<T> other );
  private:
    T   _mean;
    T  _delta;
};



template <typename T>
inline static AngularIntervalBase<T>  smallestCovering( AngularIntervalBase<T> first, AngularIntervalBase<T> second ) 
{
  return first.smallestCovering(second);
}


template <typename T>
inline  AngularIntervalBase<T>::AngularIntervalBase ( )
  : _mean(0.), _delta(0.)
{ }


template <typename T>
inline  AngularIntervalBase<T>::AngularIntervalBase ( T inf, T sup )
{ 
  if ( sup<inf ) {
    _delta = 0.;
    _mean  = 0.;
  } else if ( sup>inf+2*M_PI ) {
    _delta = M_PI;
    _mean  = 0.;
  } else {
    _delta = (sup-inf)/2;
    _mean  = Geom::Phi<T>((sup+inf)/2);
  }
}


template <typename T>
inline T  AngularIntervalBase<T>::center ( )  const
{
  return _mean;
}


template <typename T>
inline IntervalBase<T>  AngularIntervalBase<T>::cosRange()  const
{
  if ( _delta==0. ) return IntervalBase<T>();
  T max = std::max(cos(_mean-_delta),cos(_mean+_delta));
  T min = std::min(cos(_mean-_delta),cos(_mean+_delta));
  if ( include( 0   ) ) max= 1;
  if ( include(-M_PI) ) min=-1;
  return IntervalBase<T>(min,max);
}


template <typename T>
inline bool  AngularIntervalBase<T>::include ( T angle )  const
{
  if ( isEmpty() ) return false;
  if ( isFull () ) return true;
  angle = Geom::Phi<T>(angle);
  if ( _mean>angle )
    return ( _mean-_delta<angle ) || ( _mean+_delta>angle+2*M_PI ); 
  else
    return ( angle<_mean+_delta ) || ( angle>_mean-_delta+2*M_PI ); 
}


template <typename T>
inline bool  AngularIntervalBase<T>::isEmpty ( )  const
{ 
  return _delta==0.; 
} 


template <typename T>
inline bool  AngularIntervalBase<T>::isFull( )  const
{ 
  return _delta==M_PI; 
} 


template <typename T>
inline T  AngularIntervalBase<T>::length ( )  const
{
  return 2.*_delta;
}


template <typename T>
inline T  AngularIntervalBase<T>::lower ( )  const
{
  return _mean-_delta;
}


template <typename T>
inline T  AngularIntervalBase<T>::normalize ( T angle )  const
{ 
  if ( !include(angle) ) return 0.;
  angle = Geom::Phi<T>(angle);
  if ( _mean>angle )
    return    std::max(0.,std::max(angle-_mean+_delta,_mean+_delta-angle-2*M_PI))/(2*_delta); 
  else
    return 1.-std::max(0.,std::max(_mean+_delta-angle,angle-_mean+_delta-2*M_PI))/(2*_delta); 
}


template <typename T>
inline std::pair<IntervalBase<T>,IntervalBase<T>>  AngularIntervalBase<T>::normalize ( AngularIntervalBase<T> other )  const
{ 
  if ( isEmpty() || other.isEmpty() )
    return std::pair<IntervalBase<T>,IntervalBase<T>>(IntervalBase<T>()     ,IntervalBase<T>());
  if ( other.isFull()  )
    return std::pair<IntervalBase<T>,IntervalBase<T>>(IntervalBase<T>(0.,1.),IntervalBase<T>());
  IntervalBase<T>  first ( std::min(1.,std::max(0.,(other.lower()-lower()       )/length())) , 
                           std::min(1.,std::max(0.,(other.upper()-lower()       )/length())) ); 
  IntervalBase<T>  second( std::min(1.,std::max(0.,(other.lower()-lower()+2*M_PI)/length())) , 
                           std::min(1.,std::max(0.,(other.upper()-lower()+2*M_PI)/length())) ); 
  IntervalBase<T>  third ( std::min(1.,std::max(0.,(other.lower()-lower()-2*M_PI)/length())) , 
                           std::min(1.,std::max(0.,(other.upper()-lower()-2*M_PI)/length())) ); 
  assert ( second.isEmpty() || third.isEmpty() );
  if ( third.isEmpty() )
    return std::pair<IntervalBase<T>,IntervalBase<T>>(first,second);
  else
    return std::pair<IntervalBase<T>,IntervalBase<T>>(first,third );
}


template <typename T>
inline AngularIntervalBase<T> AngularIntervalBase<T>::operator+ ( AngularIntervalBase<T> other ) const
{
  if ( other.isEmpty() ) return *this;
  if ( this->isEmpty() ) return other;
  return AngularIntervalBase<T>( this->lower()+other.lower(), this->upper()+other.upper() );
}


template <typename T>
inline AngularIntervalBase<T> AngularIntervalBase<T>::operator- ( AngularIntervalBase<T> other ) const
{
  if ( other.isEmpty() ) return *this;
  if ( this->isEmpty() ) return AngularIntervalBase<T>( -other.upper(), -other.lower() ) ;
  return AngularIntervalBase<T>( this->lower()-other.upper(), this->upper()-other.lower() );
}


template <typename T>
inline AngularIntervalBase<T> AngularIntervalBase<T>::operator- ( ) const
{
  if ( isEmpty() || isFull() ) return *this;
  return AngularIntervalBase<T>( -this->upper(), -this->lower() );
}


template <typename T>
inline bool  AngularIntervalBase<T>::overlaps ( AngularIntervalBase<T> other ) const { 
  if ( isEmpty() || other.isEmpty() ) return false;
  if ( isFull () || other.isFull () ) return true;
  if ( _mean>other._mean )
    return ( _mean-_delta<other._mean+other._delta ) || ( _mean+_delta>other._mean-other._delta+2*M_PI ); 
  else
    return ( other._mean-other._delta<_mean+_delta ) || ( other._mean+other._delta>_mean-_delta+2*M_PI ); 
}


template <typename T>
inline IntervalBase<T>  AngularIntervalBase<T>::sinRange()  const
{
  if ( _delta==0. ) return IntervalBase<T>();
  T max = std::max(sin(_mean-_delta),sin(_mean+_delta));
  T min = std::min(sin(_mean-_delta),sin(_mean+_delta));
  if ( include( M_PI/2.) ) max= 1;
  if ( include(-M_PI/2.) ) min=-1;
  return IntervalBase<T>(min,max);
}


template <typename T>
inline AngularIntervalBase<T>  AngularIntervalBase<T>::smallestCovering ( AngularIntervalBase<T> other )  const
{
  if ( isEmpty()       ) return other;
  if ( other.isEmpty() ) return *this;
  if ( center()>other.center() ) {
    if ( center()-other.center()<M_PI )
      return AngularIntervalBase<T>( std::min(lower(),other.lower()), std::max(upper(),other.upper()) );
    else
      return AngularIntervalBase<T>( std::min(lower(),other.lower()+2*M_PI), std::max(upper(),other.upper()+2*M_PI) );							    
  } else {
    if ( other.center()-center()<M_PI )
      return AngularIntervalBase<T>( std::min(lower(),other.lower()), std::max(upper(),other.upper()) );
    else
      return AngularIntervalBase<T>( std::min(lower()+2*M_PI,other.lower()), std::max(upper()+2*M_PI,other.upper()) );							    
  }
}


template <typename T>
inline T  AngularIntervalBase<T>::upper ( )  const
{
  return _mean+_delta;
}


template <typename T>
inline AngularIntervalBase<T> &  AngularIntervalBase<T>::turnClockWise ( T angle )
{ 
  if ( isEmpty() || isFull() ) return *this;
  _mean = Geom::Phi<T>(_mean-angle);
  return *this;
}


template <typename T>
inline AngularIntervalBase<T> &  AngularIntervalBase<T>::turnClockWise ( AngularIntervalBase<T> other )
{ 
  if ( isEmpty() || isFull() ) return *this;
  if ( other.isEmpty() || other.isFull() ) {
    *this = other;
  } else {
    *this = *this-other;
  }
  return *this;
}


template <typename T>
inline AngularIntervalBase<T> &  AngularIntervalBase<T>::turnCounterClockWise ( T angle )
{ 
  if ( isEmpty() || isFull() ) return *this;
  _mean = Geom::Phi<T>(_mean+angle);
  return *this;
}


template <typename T>
inline AngularIntervalBase<T> &  AngularIntervalBase<T>::turnCounterClockWise ( AngularIntervalBase<T> other )
{ 
  if ( isEmpty() || isFull() ) return *this;
  if ( other.isEmpty() || other.isFull() ) {
    *this = other;
  } else {
    *this = *this+other;
  }
  return *this;
}


typedef  AngularIntervalBase<double>  AngularInterval;


#endif //  HTTrackSeeding_AngularInterval_H
