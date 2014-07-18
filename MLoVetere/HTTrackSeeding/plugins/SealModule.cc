#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoTracker/TkTrackingRegions/interface/OrderedHitsGeneratorFactory.h"
#include "MLoVetere/HTTrackSeeding/plugins/HTTripletGenerator.h"

DEFINE_EDM_PLUGIN(OrderedHitsGeneratorFactory, HTTripletGenerator, "HTTripletGenerator");

