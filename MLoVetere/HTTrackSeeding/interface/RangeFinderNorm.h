#ifndef HTTrackSeeding_RangeFinderNorm_H
#define HTTrackSeeding_RangeFinderNorm_H

/*** \class  RangeFinderNorm
  *
  *  Under development.
  *
  *  \author Maurizio Lo Vetere
  */

#include "MLoVetere/HTTrackSeeding/interface/AngularInterval.h"
#include "MLoVetere/HTTrackSeeding/interface/Interval.h"


class RangeFinderNorm
{
  public:
    RangeFinderNorm          ( )  { };
    virtual ~RangeFinderNorm ( )  { };
    virtual Interval         dTal   ( int hfturn =0 )                  const =0;
    virtual AngularInterval  dPhi   ( int hfturn =0 )                  const =0;
    virtual Interval         dEta   ( Interval lip  , int hfturn =0 )  const =0;
    virtual Interval         dTheta ( Interval lip  , int hfturn =0 )  const =0;
  protected:
    static bool    arcLengthDerCurvSignGivenNormTipCurv ( double tip, double curv );
    static double  arcLengthGivenNormTipCurv            ( double tip, double curv );
    static double  arcLengthMinimumEstimateGivenNormTip ( double tip  );
    static double  etaGivenLipArcLength                 ( double lip, double arcl );
    static double  sinPhiGivenNormTipCurv               ( double tip, double curv );
    static double  thetaGivenLipArcLength               ( double lip, double arcl );
};


#endif // HTTrackSeeding_RangeFinderNorm_H
