#ifndef HTTrackSeeding_RangeFinderNFast_H
#define HTTrackSeeding_RangeFinderNFast_H

/*** \class  RangeFinderNFast
  *
  *  Under development.
  *
  *  \author Maurizio Lo Vetere
  */

#include "MLoVetere/HTTrackSeeding/interface/AngularInterval.h"
#include "MLoVetere/HTTrackSeeding/interface/Interval.h"
#include "MLoVetere/HTTrackSeeding/interface/RangeFinderNorm.h"


class RangeFinderNFast : public RangeFinderNorm
{
  public:
    RangeFinderNFast  ( Interval curv, Interval tip );
    virtual AngularInterval  dPhi   ( int hfturn =0 )                const;
    virtual Interval         dTal   ( int hfturn =0 )                const;
    virtual Interval         dEta   ( Interval lip, int hfturn =0 )  const;
    virtual Interval         dTheta ( Interval lip, int hfturn =0 )  const;
  private:
    virtual void      cacheArcLengthOneHalfRange ( )  const;
    virtual void      cachePhiRange              ( );
  private:
    const int         NoCacheArcLengthOneHalf; 
    const int         NoCacheArcLength;
    const Interval    _curv;
    const Interval    _tip;
    AngularInterval   _forwPhiRange;
    AngularInterval   _backPhiRange;
    mutable int       _hfturn;
    mutable Interval  _arcLengthOneHalf;
    mutable Interval  _arcLength;
};


#endif // HTTrackSeeding_RangeFinderNFast_H
