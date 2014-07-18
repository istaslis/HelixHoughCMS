#ifndef HTTrackSeeding_HTTripletGenerator_H
#define HTTrackSeeding_HTTripletGenerator_H

/*** \class  HTTripletGenerator
  *
  *  This plugin is used to generate hit triples seeds using Hough transform method
  *
  *  \author Maurizio Lo Vetere
  */

#include "FWCore/ParameterSet/interface/ParameterSet.h" 
#include "MLoVetere/HTTrackSeeding/interface/HelixParNBins.h"
#include "MLoVetere/HTTrackSeeding/interface/HelixParRange.h"
#include "MLoVetere/HTTrackSeeding/interface/HelixParResolution.h"
#include "RecoTracker/TkHitPairs/interface/LayerHitMapCache.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedingLayerSets.h"
#include "RecoPixelVertexing/PixelTriplets/interface/HitTripletGenerator.h"

#include <string>


class HTTripletGenerator : public HitTripletGenerator {
  public:
    HTTripletGenerator          ( const edm::ParameterSet & cfg );
    virtual ~HTTripletGenerator ( ) { };
    virtual void  hitTriplets   ( const TrackingRegion & reg, OrderedHitTriplets & prs, const edm::Event & ev, const edm::EventSetup & es );
    void          init          ( const edm::EventSetup & es, const TrackingRegion & reg );
  private:
    edm::ParameterSet             theConfig;
    HelixParNBins                 theNumBins;
    HelixParResolution            theMinRes;
    HelixParResolution            theMaxRes;
    unsigned int                  thePHTurns;
    unsigned int                  theNHTurns;
    unsigned int                  theRequiredLayers;
    std::string                   theLayerBuilderName;
    ctfseeding::SeedingLayerSets  theLayerSets;
    float                         theField;
    GlobalPoint                   theRefPoint; 
    HelixParRange                 theRange;
    std::string			  vertexSrc;
    std::string           seedSrc;
    bool                  seedfromHits;
};


#endif // HTTrackSeeding_HTTripletGenerator_H
