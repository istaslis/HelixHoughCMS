#ifndef HTTrackSeeding_SimpleTrack3D_H
#define HTTrackSeeding_SimpleTrack3D_H

/*** \class  SimpleTrack3D
  *
  *  WARNING. This class has to go away quickly!
  *
  *  \author Maurizio Lo Vetere
  */
#include "TObject.h"

#include "MLoVetere/HTTrackSeeding/interface/AngularInterval.h"
#include "MLoVetere/HTTrackSeeding/interface/Interval.h"
#include "MLoVetere/HTTrackSeeding/interface/SimpleHit3D.h"


class SimpleTrack3D
{
  public:
  SimpleTrack3D ( ) : nlayers(0), nhits(0), shared(0)  { };
    SimpleTrack3D ( Interval curv_0, Interval eta_0, Interval lip_0, AngularInterval phi_0, Interval tip_0 ) 
      :  curv(curv_0), eta(eta_0), lip(lip_0), phi(phi_0), tip(tip_0), nlayers(0), nhits(0), shared(0)  { }
   ~SimpleTrack3D ( ) { }
  public: 
    std::vector<SimpleHit3D>  hits;
    Interval         curv;
    Interval         eta;
    Interval         lip;
    AngularInterval  phi;
    Interval         tip;
    int              nlayers;
    int              nhits;
    int              shared;
};


#endif // HTTrackSeeding_SimpleTrack3D_H
