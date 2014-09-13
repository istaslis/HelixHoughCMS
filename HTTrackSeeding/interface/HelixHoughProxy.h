#ifndef HTTrackSeeding_HelixHoughProxy_H
#define HTTrackSeeding_HelixHoughProxy_H

/*** \class  HelixHoughProxy
  *
  *  This class is a proxy to and/or wrapper of HelixHough class
  *
  *  \author Maurizio Lo Vetere
  */

#include "MLoVetere/HTTrackSeeding/interface/HelixHoughInterface.h"
#include "MLoVetere/HTTrackSeeding/interface/SimpleHit3D.h"
#include "MLoVetere/HTTrackSeeding/interface/SimpleTimer.h"
#include "MLoVetere/HTTrackSeeding/interface/SimpleTrack3D.h"

#include <vector>

class HelixHoughProxy
{
  public:
    HelixHoughProxy ( HelixHough      & context ) : _context(context       ) { }
    HelixHoughProxy ( HelixHoughProxy & other   ) : _context(other._context) { }

    /*
     *  Stop recursion
     */
    bool   breakRecursion    ( const std::vector<SimpleHit3D> &  hits  ,
                               const HelixParRange            &  range )  const { return _context.breakRecursion(hits,range); }

    float  decreasePerZoom   ( )  const  { return _context.decreasePerZoom();   } 

   
    /*
     *  Target resolution on helix parameters
     */
    const HelixParResolution &  maximumResolution ( )  const  { return _context.maximumResolution(); }
    const HelixParResolution &  minimumResolution ( )  const  { return _context.minimumResolution(); }

    /*
     *  Number of half turns that define the helix segment
     *  -  in stands for  incoming (wrt to tip) half turns
     *  - out stands for outcoming (wrt to tip) half turns
     */
    float   inHalfTurns  ( )  const  { return _context.negHalfTurns(); }
    float  outHalfTurns  ( )  const  { return _context.posHalfTurns(); }
 
    /*
     *  Additional errors to take in account the multiple scattering
     */
    float  etaError ( const SimpleHit3D &  hit    ,
                      float                minCurv,
                      float                maxCurv,
                      float                minEta ,
                      float                maxEta )  const  { return _context.etaError(hit,minCurv,maxCurv,minEta,maxEta); }
    float  phiError ( const SimpleHit3D &  hit    ,
                      float                minCurv,
                      float                maxCurv,
                      float                minEta ,
                      float                maxEta )  const  { return _context.phiError(hit,minCurv,maxCurv,minEta,maxEta); }

    const HelixParRange &  topRange             ( )  const  { return _context.range  ();        }

    /*
     *  Simple timers to compute voting in 
     */ 
    SimpleTimer &  voteTime   ( )  const  { return _context.voteTime  (); }
    SimpleTimer &  voteTimeXY ( )  const  { return _context.voteTimeXY(); }
    SimpleTimer &  voteTimeZ  ( )  const  { return _context.voteTimeZ (); }

    /*
     *  Find tracks
     */
    void  findTracks       ( const std::vector<SimpleHit3D>    &  hits   ,
                             std::vector<SimpleTrack3D>        &  tracks ,
                             const HelixParRange               &  range  )  { return _context.findTracks(hits,tracks,range); }
    void  findSeededTracks ( const std::vector<SimpleTrack3D>  &  seeds  ,
                             const std::vector<SimpleHit3D>    &  hits   ,
                             std::vector<SimpleTrack3D>        &  tracks ,
                             const HelixParRange               &  range  )  { return _context.findSeededTracks(seeds,hits,tracks,range); }

  private:
    HelixHough &  _context;
};


#endif // HTTrackSeeding_HelixHoughProxy_H
