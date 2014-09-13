#include "MLoVetere/HTTrackSeeding/interface/RangeFinderNExact.h"

#include "MLoVetere/HTTrackSeeding/interface/Interval.h"

#include <algorithm>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <utility>
#include <vector>


/*
 *  This function set tip and curv intervals and compute and caches some values in order to speed ranges calculations
 *  other computations returns the arc length and phi range combinations compatible with given tip and curv ranges.
 *  The case hfturns>1 and hbturns>1 could be improved... but for the moment I don't have time.
 *  In general we can improve computation caching values...
 */

RangeFinderNExact::RangeFinderNExact ( Interval curv, Interval tip )
  : NoCachePhi(std::numeric_limits<int>::max()),NoCacheArcLength (NoCachePhi-1), _curv(curv), _tip(tip), _forwPhiRange(), _backPhiRange(), _hfturn(NoCacheArcLength), _arcLengthOneHalf(), _arcLength()
{
  tip.setBound(Interval(-1.,1.));
  cachePhiRange();
}


/*
 *  This function returns the phi range compatible with given tip and curv ranges and turns number.
 */

AngularInterval  RangeFinderNExact::dPhi ( int hfturn )  const
{
  return ( hfturn<0 ) ? _backPhiRange : _forwPhiRange;
}


/*
 *  This function returns the (transverst) arc length compatible with given tip and curv  tip and curv ranges and turns number.
 */

Interval  RangeFinderNExact::dTal ( int hfturn )  const
{
  if ( _hfturn==NoCacheArcLength )  cacheInitArcLengthOneHalfRange();
  if ( _hfturn==NoCacheArcLength || _hfturn!=hfturn ) {
    _arcLength = Interval();
    if ( !_forwPhiRange.isEmpty() && !_arcLengthOneHalf.isEmpty() ) {
      div_t res = std::div(hfturn+1,2);
      assert( res.rem==0 || res.rem==1     );   // make sure it doesn't do sylly things on negative numbers
      assert( hfturn+1==2*res.quot+res.rem );   // make sure it doesn't do sylly things on negative numbers
      if ( res.rem==1 ) {
        _arcLength =   _arcLengthOneHalf; 
      } else {
        _arcLength = - _arcLengthOneHalf;
      }
      if ( res.quot!=0 ) {
        double sup = ( _curv.lower()*_curv.upper()>0 )
                   ? 1./std::min( std::abs(_curv.lower()), std::abs(_curv.upper()) )
                   : std::numeric_limits<float>::max();
        double inf = 1./std::max( std::abs(_curv.lower()), std::abs(_curv.upper()) );
        _arcLength = _arcLength + Interval(inf,sup).scale(2*M_PI*res.quot);

      }
      hfturn<0 ? _arcLength.setUpperBound(0.) : _arcLength.setLowerBound(0.);
      _hfturn    = hfturn;
    }
  } 
  return _arcLength;
}


/*
 *  This function returns the eta range compatible with given tip, curv, lip ranges and turns number.
 */

Interval  RangeFinderNExact::dEta ( Interval lip, int hfturn )  const
{
  std::vector<Interval> alist = etaRangeGivenLipArcLength( lip, dTal(hfturn) );
  assert( alist.size()<2 );
  Interval avalue;
  if ( alist.size()==1 ) avalue = alist[0];
  return   avalue;
}


/*
 *  This function returns the theta combinations compatible with given tip, curv, lip ranges and turns number.
 */

Interval  RangeFinderNExact::dTheta ( Interval lip, int hfturn )  const
{
  std::vector<Interval> alist = thetaRangeGivenLipArcLength( lip, dTal(hfturn) );
  assert( alist.size()<2 );
  Interval avalue;
  if ( alist.size()==1 ) avalue = alist[0];
  return   avalue;
}


/*  **************************************************************************************************************************  */

/*
 *  This function return the arc length ranges in the transverse plane between the poca and the hit as a function of
 *  tip and curv range. The calculation ot the range is non exact everywhere. In some region the arc length range is 
 *  sligthty approximated to best value. It returns a list of ranges. Null measure intervals are discarded.
 */

void  RangeFinderNExact::cacheInitArcLengthOneHalfRange ( )  const
{
  const double epsilon = 1e-15;
  double max_t = std::max( std::min(_tip.upper(),1.),-1.);
  double min_t = std::max( std::min(_tip.lower(),1.),-1.);
  double max_c = _curv.upper();
  double min_c = _curv.lower();
  std::vector<double> sValues;
  sValues.reserve(12);
  std::vector<Interval> alist;
  if ( ( max_t==1. && max_c>=1. && min_c<=1. ) || ( max_t==-1. && max_c>=-1. && min_c<=-1. ) ) {
     _arcLengthOneHalf = Interval(0,M_PI);
  }
  // the function doesn't have any minima or maxima inside the tip, curv range minima and maxima
  // on the horizontal boundaries
  for ( int i=0; i<2; i++ ) {
    double curv = (i==0) ? min_c : max_c;
    if        ( curv>1. ) {
      double tip = 2./curv-1.-epsilon;
      if        ( max_t<=tip ) {
         sValues.push_back( arcLengthGivenNormTipCurv(max_t,curv) ); 
         sValues.push_back( arcLengthGivenNormTipCurv(min_t,curv) ); 
      } else if ( min_t<=tip ) {
         sValues.push_back( arcLengthGivenNormTipCurv(tip  ,curv) ); 
         sValues.push_back( arcLengthGivenNormTipCurv(min_t,curv) ); 
      }        
    } else if ( curv>-1. ) {
      double tip = (fabs(curv)<epsilon) ? curv/2.: (1.-sqrt(1.-curv*curv))/curv;
      sValues.push_back( arcLengthGivenNormTipCurv(min_t,curv) );
      sValues.push_back( arcLengthGivenNormTipCurv(max_t,curv) );
      if ( (min_t<tip) && (tip<max_t) ) sValues.push_back( arcLengthGivenNormTipCurv(tip,curv) ); 
    } else {
      double tip = 1.+2./curv+epsilon;
      if        ( tip<min_t ) {
        sValues.push_back( arcLengthGivenNormTipCurv(max_t,curv) ); 
        sValues.push_back( arcLengthGivenNormTipCurv(min_t,curv) ); 
      } else if ( tip<max_t ) {
        sValues.push_back( arcLengthGivenNormTipCurv(max_t,curv) ); 
	sValues.push_back( arcLengthGivenNormTipCurv(tip  ,curv) ); 
      } 
    }
  }
  // minima and maxima on the vertical boundaries
  for ( int i=0; i<2; i++ ) {
    double tip = (i==0) ? min_t : max_t; 
    if   ( tip==1. ) {
      if ( min_c<= 1. )               sValues.push_back( 0.   );
      if ( max_c>= 1. && min_c<= 1. ) sValues.push_back( M_PI );
    } if ( tip==-1.) {
      if ( max_c>=-1. )               sValues.push_back( 0.   );
      if ( max_c>=-1. && min_c<=-1. ) sValues.push_back( M_PI );
    } else {
      double sup_c =  2./(1.+tip)-epsilon;
      sup_c = std::min(sup_c,max_c);
      double inf_c = -2./(1.-tip)+epsilon;
      inf_c = std::max(inf_c,min_c);
      if ( sup_c<inf_c ) continue;
      sValues.push_back( arcLengthGivenNormTipCurv(tip,sup_c) ); 
      sValues.push_back( arcLengthGivenNormTipCurv(tip,inf_c) );
      if ( arcLengthDerCurvSignGivenNormTipCurv(tip,max_c) != arcLengthDerCurvSignGivenNormTipCurv(tip,min_c) )
        sValues.push_back( arcLengthMinimumEstimateGivenNormTip(tip) );
    }
  }
  if ( !sValues.empty() ) { 
    double min_s = *min_element( sValues.begin(), sValues.end() );
    double max_s = *max_element( sValues.begin(), sValues.end() );
    if ( min_s<max_s ) _arcLengthOneHalf = Interval(min_s,max_s);
  }
}


void  RangeFinderNExact::cachePhiRange ( )
{
  double min_t =  _tip.lower();
  double max_t =  _tip.upper();
  double min_c = _curv.lower();
  double max_c = _curv.upper();
  min_t = std::max(std::min(min_t,1.),-1.);
  max_t = std::max(std::min(max_t,1.),-1.);
  // the function doesn't have any minima or maxima inside the tip,curv ranges
  // the domain is complex but we don't have to take this into account
  // add minima and maxima on the corners 
  // we don't have to take into account if the function is evaluated outside domain
  double max_s = sinPhiGivenNormTipCurv(max_t,min_c);
  if ( fabs(min_c)>1. && max_s< 1. ) max_s = std::max( max_s, sinPhiGivenNormTipCurv(min_t,min_c) );
  double min_s = sinPhiGivenNormTipCurv(min_t,max_c);
  if ( fabs(max_c)>1. && min_s>-1. ) min_s = std::min( min_s, sinPhiGivenNormTipCurv(max_t,max_c) );
  // the function doesn't have any minima or maxima beloging to the border tip = min_t or max_t
  // add minima or maxima beloging to the border curv = min_c or max_c
  if ( min_c>= 1. ) {
    double inv_c = 1./min_c;
    double tip   = inv_c-sqrt(1.-inv_c*inv_c);
    if ( (min_t<tip) && (tip<max_t) ) max_s = std::max( max_s, sinPhiGivenNormTipCurv(tip,min_c) );
  }
  if ( max_c<=-1. ) {
    double inv_c = 1./max_c;
    double tip   = inv_c+sqrt(1.-inv_c*inv_c);
    if ( (min_t<tip) && (tip<max_t) ) min_s = std::min( min_s, sinPhiGivenNormTipCurv(tip,max_c) );
  }
  if ( min_s<max_s ) {
    double  phi_l = asin(min_s);
    double  phi_h = asin(max_s);
    _forwPhiRange = AngularInterval( phi_l, phi_h );
    _backPhiRange = AngularInterval( M_PI-phi_h, M_PI-phi_l );
  }
}


/*
 *  This function returns a list of eta ranges as a function of range of normalized lip of the poca and the allowed arc lengths
 *  in the transverse plane between the poca and the hit.
 */

std::vector<Interval>  RangeFinderNExact::etaRangeGivenLipArcLength ( Interval lip, Interval arcl )
{
  std::vector<Interval> alist;
  Interval parcl = arcl;
  parcl.setLowerBound(0.);
  if ( !parcl.isEmpty() ) {
    double min_eta = (   lip.upper() <0. ) ? etaGivenLipArcLength(lip.upper(),parcl.upper()) : 
                     ( parcl.lower()!=0. ) ? etaGivenLipArcLength(lip.upper(),parcl.lower()) : -std::numeric_limits<float>::max();
    double max_eta = (   lip.lower() >0. ) ? etaGivenLipArcLength(lip.lower(),parcl.upper()) : 
                     ( parcl.lower()!=0. ) ? etaGivenLipArcLength(lip.lower(),parcl.lower()) :  std::numeric_limits<float>::max();
    if ( min_eta<max_eta ) addToPairWiseDisjointIntervalSet(alist,Interval(min_eta,max_eta));
  }
  Interval narcl = arcl;
  narcl.setUpperBound(0.);
  if ( !narcl.isEmpty() ) {
    double min_eta = (   lip.lower() >0. ) ? etaGivenLipArcLength(lip.lower(),narcl.lower()) :
                     ( narcl.upper()!=0. ) ? etaGivenLipArcLength(lip.lower(),narcl.upper()) : -std::numeric_limits<float>::max();
    double max_eta = (   lip.upper() <0. ) ? etaGivenLipArcLength(lip.upper(),narcl.lower()) : 
                     ( narcl.upper()!=0. ) ? etaGivenLipArcLength(lip.upper(),narcl.upper()) :  std::numeric_limits<float>::max();
    if ( min_eta<max_eta ) addToPairWiseDisjointIntervalSet(alist,Interval(min_eta,max_eta));
  }
  return alist;
}


/*
 *  This function returns a list of theta ranges as a function of range of normalized lip of the poca and the allowed arc lengths
 *  in the transverse plane between the poca and the hit.
 */

std::vector<Interval>  RangeFinderNExact::thetaRangeGivenLipArcLength ( Interval lip, Interval arcl )
{
  std::vector<Interval> alist;
  Interval parcl = arcl;
  parcl.setLowerBound(0.);
  if ( !parcl.isEmpty() ) {
    double max_theta = (  lip.upper() <0. ) ? thetaGivenLipArcLength(lip.upper(),parcl.upper()) : thetaGivenLipArcLength(lip.upper(),parcl.lower());
    double min_theta = (  lip.lower() >0. ) ? thetaGivenLipArcLength(lip.lower(),parcl.upper()) : thetaGivenLipArcLength(lip.lower(),parcl.lower());
    if ( min_theta<max_theta ) addToPairWiseDisjointIntervalSet(alist,Interval(min_theta,max_theta));
  }
  Interval narcl = arcl;
  narcl.setUpperBound(0.);
  if ( !narcl.isEmpty() ) {
    double max_theta = (   lip.lower() >0. ) ? thetaGivenLipArcLength(lip.lower(),narcl.lower()) :
                       ( narcl.upper()!=0. ) ? thetaGivenLipArcLength(lip.lower(),narcl.upper()) : -std::numeric_limits<float>::max();
    double min_theta = (   lip.upper() <0. ) ? thetaGivenLipArcLength(lip.upper(),narcl.lower()) : 
                       ( narcl.upper()!=0. ) ? thetaGivenLipArcLength(lip.upper(),narcl.upper()) :  std::numeric_limits<float>::max();
    if ( min_theta<max_theta ) addToPairWiseDisjointIntervalSet(alist,Interval(min_theta,max_theta));
  }
  return alist;
}
