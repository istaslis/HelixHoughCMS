#ifndef HTTrackSeeding_RangeFinderNExact_H
#define HTTrackSeeding_RangeFinderNExact_H

/*** \class  RangeFinderNExact
  *
  *  Under development.
  *
  *  \author Maurizio Lo Vetere
  */

#include "MLoVetere/HTTrackSeeding/interface/AngularInterval.h"
#include "MLoVetere/HTTrackSeeding/interface/Interval.h"
#include "MLoVetere/HTTrackSeeding/interface/RangeFinderNorm.h"


class RangeFinderNExact : public RangeFinderNorm
{
  public:
    RangeFinderNExact ( Interval curv, Interval tip );
    virtual AngularInterval  dPhi   ( int hfturn =0 )                const;
    virtual Interval         dTal   ( int hfturn =0 )                const;
    virtual Interval         dEta   ( Interval lip, int hfturn =0 )  const;
    virtual Interval         dTheta ( Interval lip, int hfturn =0 )  const;
  private:
    virtual void      cacheInitArcLengthOneHalfRange ( )  const;
    virtual void      cachePhiRange                  ( );
  private:
    static std::vector<Interval>  etaRangeGivenLipArcLength   ( Interval lip, Interval arcl );
    static std::vector<Interval>  thetaRangeGivenLipArcLength ( Interval lip, Interval arcl );
  private:
    const int         NoCachePhi ;
    const int         NoCacheArcLength;
    const Interval    _curv;
    const Interval    _tip;
    AngularInterval   _forwPhiRange;
    AngularInterval   _backPhiRange;
    mutable int       _hfturn;
    mutable Interval  _arcLengthOneHalf;
    mutable Interval  _arcLength;
};


#endif // HTTrackSeeding_RangeFinderNExact_H
