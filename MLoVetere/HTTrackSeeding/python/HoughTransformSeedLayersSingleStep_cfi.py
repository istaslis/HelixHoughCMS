import FWCore.ParameterSet.Config as cms

from RecoLocalTracker.Configuration.RecoLocalTracker_cff import *

htSeedLayers = cms.Sequence(siPixelRecHits*siStripMatchedRecHits)

from MLoVetere.HTTrackSeeding.HoughTransformSeedLayersNC_cfi import *

