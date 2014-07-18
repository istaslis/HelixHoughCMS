import FWCore.ParameterSet.Config as cms

houghTransformStepClusters = cms.EDProducer("TrackClusterRemover",
    clusterLessSolution    = cms.bool    (True),
    oldClusterRemovalInfo  = cms.InputTag("detachedTripletStepClusters"),
    trajectories           = cms.InputTag("detachedTripletStepTracks"),
    overrideTrkQuals       = cms.InputTag('detachedTripletStep'),
    TrackQuality           = cms.string  ('highPurity'),
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    pixelClusters          = cms.InputTag("siPixelClusters"),
    stripClusters          = cms.InputTag("siStripClusters"),
    Common = cms.PSet(
        maxChi2 = cms.double(9.0)
    )
)

from RecoLocalTracker.Configuration.RecoLocalTracker_cff import *
from RecoTracker.IterativeTracking.iterativeTk_cff import *

htSeedLayers = cms.Sequence(
    siPixelRecHits*siStripMatchedRecHits*
    InitialStep*LowPtTripletStep*PixelPairStep*DetachedTripletStep*
    houghTransformStepClusters
)

from MLoVetere.HTTrackSeeding.HoughTransformSeedLayers_cfi import *

