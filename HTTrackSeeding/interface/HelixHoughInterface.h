#ifndef HTTrackSeeding_HelixHough_H
#define HTTrackSeeding_HelixHough_H

/*** \class  HelixHough
  *
  *  This class define the external interface and provided context for HelixHoughEngine
  *
  *  \author Maurizio Lo Vetere
  */

#include "MLoVetere/HTTrackSeeding/interface/HelixParNBins.h"
#include "MLoVetere/HTTrackSeeding/interface/HelixParRange.h"
#include "MLoVetere/HTTrackSeeding/interface/HelixParResolution.h"
#include "MLoVetere/HTTrackSeeding/interface/SimpleHit3DInterface.h"
#include "MLoVetere/HTTrackSeeding/interface/SimpleTimer.h"
#include "MLoVetere/HTTrackSeeding/interface/SimpleTrack3DInterface.h"
#include "RecoPixelVertexing/PixelTriplets/interface/OrderedHitTriplets.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "RecoTracker/TkTrackingRegions/interface/TrackingRegion.h"

#include "MLoVetere/HTTrackSeeding/interface/SimpleHit3D.h"
#include "MLoVetere/HTTrackSeeding/interface/SimpleTrack3D.h"



#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include <vector>


class HelixHough
{
  public:
    HelixHough               ( GlobalPoint        origin         ,
                               HelixParRange      range          ,
                               HelixParNBins      nbins          ,
                               HelixParResolution minResolution  ,
                               HelixParResolution maxResolution  ,
                               unsigned int       requiredLayers );
    virtual  ~HelixHough     ( );
    void  findHelices        ( const TrackingRegion::Hits & hits         ,
			       std::vector<unsigned int> &HitsInLayers, 
                               unsigned int                 min_hits     ,
                               unsigned int                 max_hits     ,
                               OrderedHitTriplets         & tracks       ,
			       std::vector<std::vector<unsigned int> > &zoomlevels ,
			       double phimin,
			       double phimax,
			       double SV,
			       double dzdlmin,
			       double dzdlmax,
			       double pTmin,
			       double debugOutput,
			       double HoughTransformSeeding,
			       double SeedRefitterOnly,
                               unsigned int                 maxtracks =0 );
    void  findSeededHelices  ( std::vector<SimpleTrack3D> & seeds        ,
                               const TrackingRegion::Hits & hits         ,
                               unsigned int                 min_hits     ,
                               unsigned int                 max_hits     ,
                               OrderedHitTriplets         & tracks       ,
                               unsigned int                 maxtracks =0 );
    void  setDecreasePerZoom ( float value )  { _decreasePerZoom = value; }    
  public:
    virtual bool                breakRecursion    ( const std::vector<SimpleHit3D>   & hits   ,
                                                    const HelixParRange              & range  )  const;
    float                       decreasePerZoom   ( )                                            const  { return _decreasePerZoom;   }
    const HelixParResolution &  maximumResolution ( )                                            const  { return _maximumResolution; }
    const HelixParResolution &  minimumResolution ( )                                            const  { return _minimumResolution; }
    unsigned int                negHalfTurns      ( )                                            const  { return 0;      }
    unsigned int                posHalfTurns      ( )                                            const  { return 1;      }
    const HelixParRange      &  range             ( )                                            const  { return _range; } 
    virtual float               etaError          ( const SimpleHit3D & hit,
                                                    float min_curv ,
                                                    float max_curv ,
                                                    float min_eta  ,
                                                    float max_eta  )                             const  { return 0.;     }
    virtual float               phiError          ( const SimpleHit3D & hit,
                                                    float min_curv ,
                                                    float max_curv ,
                                                    float min_eta  ,
                                                    float max_eta  )                             const  { return 0.;     }
    virtual void   findTracks       ( const std::vector<SimpleHit3D>   & hits   ,
                                      std::vector<SimpleTrack3D>       & tracks ,
                                      const HelixParRange              & range  );
    virtual void   findSeededTracks ( const std::vector<SimpleTrack3D> & seeds  ,
                                      const std::vector<SimpleHit3D>   & hits   ,
                                      std::vector<SimpleTrack3D>       & tracks ,
                                      const HelixParRange              & range  )                { };
    SimpleTimer &  voteTime         ( )                                                          const { return * _voteTime;   } 
    SimpleTimer &  voteTimeXY       ( )                                                          const { return * _voteTimeXY; }
    SimpleTimer &  voteTimeZ        ( )                                                          const { return * _voteTimeZ;  }
    
    unsigned int nEv;
    
    void SetVertex(float xvtx, float xvtxErr,float  yvtx,float  yvtxErr,float  zvtx, float zvtxErr)
    {
        xVtx = xvtx;
        xVtxErr = xvtxErr;
        yVtx = yvtx;
        yVtxErr =yvtxErr;
        zVtx = zvtx;
        zVtxErr = zvtxErr;
        
    }

    void SetSeeds(std::vector<PHENIXHough::SimpleTrack3D> seedTracks,std::vector<TransientTrackingRecHitBuilder::RecHitPointer> &hits) {seeds = seedTracks;seedHits = hits;}
  private:
    virtual void  initEvent   ( std::vector<SimpleHit3D> & hits, unsigned int min_hits )                   { }
    virtual void  initSeeding ( )                                                                          { }    
    virtual void  finalize    ( const TrackingRegion::Hits & hits, const std::vector<SimpleTrack3D> & input, OrderedHitTriplets & output );
  private:
    GlobalPoint         _origin;
    HelixParRange       _range;
    HelixParNBins       _nBins;
    HelixParResolution  _minimumResolution;
    HelixParResolution  _maximumResolution;
    unsigned int        _requiredLayers;
    float               _decreasePerZoom;
    SimpleTimer      *  _voteTime;
    SimpleTimer      *  _voteTimeXY;
    SimpleTimer      *  _voteTimeZ;
    
    
    float xVtx, xVtxErr,  yVtx,  yVtxErr,  zVtx,  zVtxErr;
    std::vector<PHENIXHough::SimpleTrack3D> seeds;
    std::vector<TransientTrackingRecHitBuilder::RecHitPointer> seedHits;
};


#endif // HTTrackSeeding_HelixHough_H
