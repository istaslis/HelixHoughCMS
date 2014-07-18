#include "MLoVetere/HTTrackSeeding/interface/RangeFinderNorm.h"

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
 *  This function is a utility used by arcLengthRangeGivenNormTipCurv. It returns the sign of the derivative of the
 *  arc length with respect to the curvature at fixed tip.
 */

bool  RangeFinderNorm::arcLengthDerCurvSignGivenNormTipCurv ( double tip, double curv )
{
  const double arg0 = tip*curv;
  bool  pos = true;
  if ( arg0<0 ) {
    const double arg1 = curv/2.*sqrt((1-tip*tip)/(1-arg0));
    pos = fabs(arg1)<1. ? (2.-arg0)/(2.-2.*arg0)>asin(arg1)/arg1*sqrt(1.-arg1*arg1) : true;
  }
  return pos==(curv>0);
}


/*
 *  This function return the adimensional transverse plane arc length parameter between poca and the hit, given tip and curv.
 *  It assume tip and curv in adimensional units (d/D, kD). Unphysical situation must be handled before calling this function.
 *  In other words it must be granted: -1<=tip<=1<=|tip-2./curv|, 1-tip*curv>=0.
 */

double  RangeFinderNorm::arcLengthGivenNormTipCurv ( double tip, double curv )
{
  const double epsilon = 1e-15;
  assert ( tip>=-1. );
  assert ( tip<= 1. );
  assert ( tip== 1. || curv >= -2./(1.-tip) );
  assert ( tip==-1. || curv <=  2./(1.+tip) );
  if ( fabs(curv)>epsilon ) {
    double arg = (1.-tip*tip)/(1.-tip*curv);
    arg = ( arg>0. ) ? curv/2.*sqrt(arg) : 0. ;
    arg = std::max( std::min(arg,1.), -1. );
    return 2./curv*asin(arg);
  } else {
    double arg = (1.-tip*tip);   
    arg = ( arg>0. ) ? sqrt(arg) : 0. ;
    return (1.+tip*curv/2.)*arg;
  }
}


/*
 *  This function is a utility used by arcLengthRangeGivenNormTipCurv. It returns a lower bound to the minimum
 *  of the arc length at fixed tip. The estimate is better than 1e-8.
 */

double   RangeFinderNorm::arcLengthMinimumEstimateGivenNormTip ( double tip )
{
  double x = fabs(tip);
  assert ( x<=1. );
  double z;
  if ( x<0.1882937000875663 ) { 
    if ( x<0.06596858831822221 ) {
      //  0. < x < 0.06596858831822221,
      const double x2 = x *x ;
      const double x3 = x *x2;
      z = 3*x - 69*x3/5.;
    } else {
      //  0.06596858831822221 < x < 0.1882937000875663
      x -= 0.1450174226103326;
      const double x2 = x *x ;
      const double x3 = x *x2;
      const double x4 = x2*x2;
      z = 0.4 +  2.3521810135995493  * x
              -  3.3451901084911166  * x2
              -  0.12636832648132904 * x3
              + 18.48769908449483    * x4;
    }
  } else {
    if ( x<0.5560316177508693 ) {
      if ( x<0.31665222267769044 ) {
        //  0.1882937000875663 < x < 0.31665222267769044
        x -= 0.24329023161205376;
        const double x2 = x *x ;
        const double x3 = x *x2;
        const double x4 = x2*x2;
        z = 0.6 + 1.7387525391545398 * x
                - 2.7494479454055574 * x2
                + 3.011527656186959  * x3
                + 0.855126803290202  * x4;
      } else { 
        //  0.31665222267769044 < x < 0.5560316177508693
        x -= 0.4327119935763125;
        const double x2 = x *x ;
        const double x3 = x *x2;
        const double x4 = x2*x2;
        z = 0.85 + 0.990539897185676  * x
                 - 1.3404236294472376 * x2
                 + 1.7050504762167324 * x3
                 - 1.8883762869891854 * x4;
      }
    } else {
      if ( x<0.8354729073528635 ) {
        //  0.5560316177508693 < x < 0.8354729073528635
        x -= 0.7118687685952387;
        const double x2 = x *x ;
        const double x3 = x *x2;
        const double x4 = x2*x2;
        z =  1.05 + 0.5168880139126358 * x
                  - 0.5178690915880463 * x2
                  + 0.5164375050182545 * x3
                  - 0.5048724399397361 * x4;
      } else {
        //  0.8354729073528635 < x < 1.
        x -= 1.0000000000000002;
        const double x2 = x *x ;
        const double x3 = x *x2;
        const double x4 = x2*x2;
        z = 1.1655611852072114 + 0.3112334710573324  * x
                               - 0.2417459512623898  * x2
                               + 0.18868564575636362 * x3
                               - 0.14704546013373293 * x4;
      }
    }
  }
  const double epsilon = 1e-8;
  const double b = ( z>epsilon ) ? sin(z)*sin(z)*tan(z)/(z-tan(z)) : -3.0; 
  const double curv = 2.*tip/(1.-tip*tip)*b;
  return std::max(0.,2./curv*asin(curv/2.*sqrt((1.-tip*tip)/(1.-tip*curv)))-1.e-8);
}


/*
 *  This function return the pseudorapidity of the track given lip, the distance of hit from the poca in z, and arcl, the 
 *  transverse plane arc length between poca and the hit. When arc length is positive, the track is assumed to go from the 
 *  poca to the hit. When arc length is negative, the track is assumed to go from the hit to the poca.
 */

double  RangeFinderNorm::etaGivenLipArcLength ( double lip, double arcl )
{
  double ds = sqrt(arcl*arcl+lip*lip); 
  if ( arcl<0. ) lip = -lip;
  return 0.5*log((ds-lip)/(ds+lip));
}


/*
 *  This function return the sine of the angle between the direction of the charged particle transverse momentum and the direction 
 *  of the hit in the transverse plane. It assume tip and curv in adimensional units ( d/D, kD ) and assume -1<=tip<=1. 
 *  When there are no solutions, this function returns -1. for positive curvature and 1. for negative curvature.
 */

double  RangeFinderNorm::sinPhiGivenNormTipCurv ( double tip, double curv )
{
  assert( tip>=-1. );
  assert( tip<= 1. );
  if ( tip>-1. && curv>= 2./(1.+tip) ) return -1.;
  if ( tip< 1. && curv<=-2./(1.-tip) ) return  1.;
  double ak = (2.*tip-tip*tip*curv-curv)/2.;
  double sinphi = ak/(1.-tip*curv);
  if ( sinphi> 1. ) sinphi= 1.;
  if ( sinphi<-1. ) sinphi=-1.;
  return sinphi;
}


/*
 *  This function return the theta angle of the track given lip, the distance of hit from the poca in z, and arcl, the transverse 
 *  plane arc length between poca and the hit. When arc length is positive, the track is assumed to go from the poca to the hit. 
 *  When arc length is negative, the track is assumed to go from the hit to the poca.
 */

double  RangeFinderNorm::thetaGivenLipArcLength ( double lip, double arcl )
{
  double ds = sqrt(arcl*arcl+lip*lip); 
  if ( arcl<0. ) lip = -lip;
  return acos(-lip/ds);
}
