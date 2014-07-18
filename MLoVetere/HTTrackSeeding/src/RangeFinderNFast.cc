#include "MLoVetere/HTTrackSeeding/interface/RangeFinderNFast.h"

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

RangeFinderNFast::RangeFinderNFast ( Interval curv, Interval tip )
  : NoCacheArcLengthOneHalf (std::numeric_limits<int>::max()), NoCacheArcLength (NoCacheArcLengthOneHalf-1), _curv(curv), _tip(tip), _hfturn(NoCacheArcLengthOneHalf)
{
  cachePhiRange();
}


/*
 *  This function returns the phi range compatible with given tip and curv ranges and turns number.
 */

AngularInterval  RangeFinderNFast::dPhi ( int hfturn )  const
{
  return ( hfturn<0 ) ? _backPhiRange : _forwPhiRange;
}


/*
 *  This function returns the (trasverse) arc length compatible with given tip and curv  tip and curv ranges and turns number.
 */

Interval  RangeFinderNFast::dTal ( int hfturn )  const
{
  if ( _forwPhiRange.isEmpty() )     return Interval();
  cacheArcLengthOneHalfRange();
  if ( _arcLengthOneHalf.isEmpty() ) return Interval();
  if ( _hfturn==NoCacheArcLength || _hfturn!=hfturn ) {
    _arcLength = Interval();
    const double  epsilon = 1e-15;
    double  curv = _curv.center();
    div_t res = std::div(hfturn+1,2);
    assert( res.rem==0 || res.rem==1     );   // make sure it doesn't do sylly things on negative numbers
    assert( hfturn+1==2*res.quot+res.rem );   // make sure it doesn't do sylly things on negative numbers
    if ( res.rem==1 ) {
      _arcLength =   _arcLengthOneHalf; 
    } else {
      _arcLength = - _arcLengthOneHalf;
    }
    if ( res.quot!=0 ) {
      double length = ( curv!=0 ) ? 1./std::abs(curv) : std::numeric_limits<float>::max();
      _arcLength = _arcLength + Interval(length-epsilon,length+epsilon).scale(2*M_PI*res.quot);
    }
    hfturn<0 ? _arcLength.setUpperBound(0.) : _arcLength.setLowerBound(0.);
    //_hfturn    = hfturn;
  } 
  return _arcLength;
}


/*
 *  This function returns the eta range compatible with given tip, curv, lip ranges and turns number.
 */

Interval  RangeFinderNFast::dEta ( Interval lip, int hfturn )  const
{
  if ( lip.isEmpty() )  return Interval();
  Interval  arclRange = dTal(hfturn);
  if ( arclRange.isEmpty() ) return Interval();
  double  arcl = arclRange.center();
  if ( arcl==0 ) return Interval();
  const double  epsilon = 1e-15;
  double eta = etaGivenLipArcLength(lip.center(),arcl);
  return Interval(eta-epsilon,eta+epsilon);
}


/*
 *  This function returns the theta combinations compatible with given tip, curv, lip ranges and turns number.
 */

Interval  RangeFinderNFast::dTheta ( Interval lip, int hfturn )  const
{
  if ( lip.isEmpty() )  return Interval();
  Interval  arclRange = dTal(hfturn);
  if ( arclRange.isEmpty() ) return Interval();
  double  arcl = arclRange.center();
  if ( arcl==0 ) return Interval();
  const double  epsilon = 1e-15;
  double  theta = thetaGivenLipArcLength(lip.center(),arcl);
  return  Interval(theta-epsilon,theta+epsilon);
}


/*  **************************************************************************************************************************  */

/*
 *  This function return the arc length ranges in the transverse plane between the poca and the hit as a function of
 *  tip and curv range. The calculation ot the range is non exact everywhere. In some region the arc length range is 
 *  sligthty approximated to best value. It returns a list of ranges. Null measure intervals are discarded.
 */

void  RangeFinderNFast::cacheArcLengthOneHalfRange ( )  const
{
  //if ( _hfturn!=NoCacheArcLengthOneHalf ) return;
  const double  epsilon = 1e-15;
  double  tip  =  _tip.center();
  double  curv = _curv.center();
  if ( (tip>=1.) || (tip<=-1.) ) return;
  double sup_c =  2./(1.+tip)-epsilon;
  double inf_c = -2./(1.-tip)+epsilon;
  if ( (inf_c<curv) && (curv<sup_c) ) {
    double length = arcLengthGivenNormTipCurv(tip,curv);
    _arcLengthOneHalf = Interval( length-epsilon, length+epsilon ); 
  } else
    _arcLengthOneHalf = Interval();
  _arcLength = Interval();
  _hfturn = NoCacheArcLength;
}


void  RangeFinderNFast::cachePhiRange ( )
{
  const double  epsilon = 1e-15;
  double  tip  =  _tip.center();
  double  curv = _curv.center();
  if ( (tip>=1.) || (tip<=-1.) ) return;
  double  sinPhi = sinPhiGivenNormTipCurv(tip,curv);
  if ( (sinPhi>=1.) || (sinPhi<=-1.) ) return; 
  double  phi = asin(sinPhi);
  _forwPhiRange = AngularInterval( phi     -epsilon, phi     +epsilon );
  _backPhiRange = AngularInterval( M_PI-phi-epsilon, M_PI-phi+epsilon );
}
