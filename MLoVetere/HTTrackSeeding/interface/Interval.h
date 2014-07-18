#ifndef HTTrackSeeding_Interval_H
#define HTTrackSeeding_Interval_H

/*** \class  Interval
  *
  *  Only sets of non null meausure are relevant in this context therefore we consider open intervals  
  *
  *  \author Maurizio Lo Vetere
  */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

template <typename T>
class IntervalBase {
  public:
    IntervalBase ( )              : _inf( 0.), _sup( 0.) { }
    IntervalBase ( T inf, T sup ) : _inf(inf), _sup(sup) { }
    T                   center            ( )                         const  { return (_sup+_inf)/2.; }
    bool                include           ( T  value )                const;
    bool                isEmpty           ( )                         const  { return _sup<=_inf; }
    T                   length            ( )                         const  { return std::max(_sup-_inf,0.); }
    T                   lower             ( )                         const  { return _inf; }
    T                   normalize         ( T  value )                const;
    IntervalBase<T>     normalize         ( IntervalBase<T> value  )  const;
    bool                operator<         ( IntervalBase<T> other  )  const;
    IntervalBase<T>     operator+         ( IntervalBase<T> other  )  const;
    IntervalBase<T>     operator-         ( IntervalBase<T> other  )  const;
    IntervalBase<T>     operator-         ( )                         const;
    bool                overlaps          ( IntervalBase<T> other  )  const;
    IntervalBase<T>     smallestCovering  ( IntervalBase<T> other  )  const;
    T                   upper             ( )                         const  { return _sup; }
    IntervalBase<T> &   downScale         ( IntervalBase<T> factor );
    IntervalBase<T> &   scale             ( T factor );
    IntervalBase<T> &   setBound          ( IntervalBase range     );
    IntervalBase<T> &   setLowerBound     ( T  value );
    IntervalBase<T> &   setUpperBound     ( T  value );
    IntervalBase<T> &   shift             ( T  value );
    IntervalBase<T> &   upScale           ( IntervalBase<T> factor );
  private:
    T _inf;
    T _sup;
};


template <typename T>
inline static bool overlaps( IntervalBase<T> first, IntervalBase<T> second ) 
{
  return first.overlaps(second);
}


// da verificare con attenzione;
template <typename T>
inline static IntervalBase<T> intersection( IntervalBase<T> first, IntervalBase<T> second ) 
{
  return first.setBound(second);
}


template <typename T>
inline static IntervalBase<T>  smallestCovering( IntervalBase<T> first, IntervalBase<T> second ) 
{
  return first.smallestCovering(second);
}


template <typename T>
inline bool  IntervalBase<T>::include ( T value )  const
{
  if ( isEmpty() ) return false;
  return ( _inf<value ) && ( value<_sup ); 
}


template <typename T>
inline T  IntervalBase<T>::normalize ( T value )  const
{ 
  if ( isEmpty()     ) return 0.;
  if ( value <= _inf ) return 0.;
  if ( value >= _sup ) return 1.;
  return (value-_inf)/(_sup-_inf);
}


template <typename T>
inline IntervalBase<T>  IntervalBase<T>::normalize ( IntervalBase<T> value )  const
{ 
  value.setBound(*this);
  if ( isEmpty() || value.isEmpty() ) return IntervalBase<T>();
  return value.shift(-_inf).scale(1./length());
}


template <typename T>
inline bool  IntervalBase<T>::operator< ( IntervalBase<T> other ) const  
{ 
  if ( other.isEmpty() || this->isEmpty() ) return false;
  return ( _sup<=other._inf ); 
}


template <typename T>
inline IntervalBase<T> IntervalBase<T>::operator+ ( IntervalBase<T> other ) const
{
  if ( other.isEmpty() ) return *this;
  if ( this->isEmpty() ) return other;
  return IntervalBase<T>( this->lower()+other.lower(), this->upper()+other.upper() );
}


template <typename T>
inline IntervalBase<T> IntervalBase<T>::operator- ( IntervalBase<T> other ) const
{
  if ( other.isEmpty() ) return *this;
  if ( this->isEmpty() ) return IntervalBase<T>( -other.upper(), -other.lower() ) ;
  return IntervalBase<T>( this->lower()-other.upper(), this->upper()-other.lower() );
}


template <typename T>
inline IntervalBase<T> IntervalBase<T>::operator- ( ) const
{
  return IntervalBase<T>( -this->upper(), -this->lower() );
}


template <typename T>
inline bool  IntervalBase<T>::overlaps ( IntervalBase<T> other ) const  
{ 
  if ( other.isEmpty() || this->isEmpty() ) return false;
  return ! ( *this<other || other<*this ); 
}


template <typename T>
inline IntervalBase<T>  IntervalBase<T>::smallestCovering ( IntervalBase<T> other  )  const
{
  if ( isEmpty()       ) return other;
  if ( other.isEmpty() ) return *this;
  return IntervalBase<T>( std::min(lower(),other.lower()), std::max(upper(),other.upper()) ); 
}


/*
 *  This function assume other.min()>0
 */

template <typename T>
inline IntervalBase<T> &  IntervalBase<T>::downScale ( IntervalBase<T> other )
{ 
  assert ( !other.isEmpty() );
  assert ( other._inf>0. );
  _sup /= other._inf; 
  _inf /= other._sup; 
  return *this; 
}


template <typename T>
inline IntervalBase<T> &  IntervalBase<T>::scale ( T factor )
{ 
  _sup *= factor; 
  _inf *= factor; 
  if ( factor<0 ) std::swap(_inf,_sup);
  return *this; 
}


template <typename T>
inline IntervalBase<T> &  IntervalBase<T>::setBound ( IntervalBase<T> range )
{
  setUpperBound( range.upper() );
  setLowerBound( range.lower() );
  return *this; 
}


template <typename T>
inline IntervalBase<T> &  IntervalBase<T>::setLowerBound ( T value )
{
  _sup = std::max( _sup, value );
  _inf = std::max( _inf, value );
  return *this; 
}


template <typename T>
inline IntervalBase<T> &  IntervalBase<T>::setUpperBound ( T value )
{
  _sup = std::min( _sup, value );
  _inf = std::min( _inf, value );
  return *this; 
}

template <typename T>
inline IntervalBase<T> &  IntervalBase<T>::shift ( T value )
{ 
  _sup += value; 
  _inf += value; 
  return *this; 
}


/*
 *  This function assume other.min()>0
 */

template <typename T>
inline IntervalBase<T> &  IntervalBase<T>::upScale ( IntervalBase<T> other )
{ 
  assert ( !other.isEmpty() );
  assert ( other._inf>0. );
  _sup *= other._sup; 
  _inf *= other._inf; 
  return *this; 
}


/*
 *  Add element to the list of intervals. If the interval overlaps with others intervals in the
 *  list, collapse all the overlapping intervals into one equivalent interval.
 */

template <typename T>
inline static void  addToPairWiseDisjointIntervalSet ( std::vector<IntervalBase<T> > & theList, IntervalBase<T> item )
{
  if ( item.isEmpty() ) return;
  std::vector<IntervalBase<T> > alist;
  typedef typename std::vector<IntervalBase<T> >::iterator IT;
  for ( IT iter=theList.begin(); iter!=theList.end(); ++iter ) {
    if ( iter->overlaps(item) ) 
      item = IntervalBase<T>( std::min(item.lower(),iter->lower()), std::max(item.upper(),iter->upper()) );        
    else
      alist.push_back(*iter);
  }
  alist.push_back(item);
  alist.swap(theList);
}  


/*
 *  Return the first bin in the bins list that overlaps with val, or -1 if it doesn't belong to any
 *  It assumes that bins is an ordered list.
 */

// Da sistemare!
// Non funziona perche' bin(val,val) e' un insieme vuoto non ha ordine. 

template <typename T>
inline static int seed_bin ( std::vector<IntervalBase<T> > & bins, T val )
{
  IntervalBase<T> bin(val,val);  
  if ( ( bin<bins[0] ) || ( bins.back()<bin ) ) return -1;
  typedef typename std::vector<IntervalBase<T> >::iterator IT;
  std::pair<IT,IT> bounds;
  bounds = equal_range( bins.begin(), bins.end(), bin );
  if ( bounds.first!=bounds.second )
    return ((int) ( bounds.first - bins.begin() ));
  else 
    return -1;
}


typedef IntervalBase<double>  Interval;


#endif // HTTrackSeeding_Interval_H
