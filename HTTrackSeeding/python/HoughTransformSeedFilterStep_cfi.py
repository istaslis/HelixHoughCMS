import FWCore.ParameterSet.Config as cms
from math import pi
###############################################################
# Hough Transform Tracking                                    #
###############################################################

#houghTransformStepClusters = cms.EDProducer("TrackClusterRemover",
#    clusterLessSolution    = cms.bool    (True),
#    oldClusterRemovalInfo  = cms.InputTag("detachedTripletStepClusters"),
#    trajectories           = cms.InputTag("detachedTripletStepTracks"),
#    overrideTrkQuals       = cms.InputTag('detachedTripletStep'),
#    TrackQuality           = cms.string  ('highPurity'),
#    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
#    pixelClusters          = cms.InputTag("siPixelClusters"),
#    stripClusters          = cms.InputTag("siStripClusters"),
#    Common = cms.PSet(
#        maxChi2 = cms.double(9.0)
#    )
#)

# SEEDING LAYERS
from MLoVetere.HTTrackSeeding.HoughTransformSeedLayersNC_cfi import *

# SEEDS
HoughTransformSeedGeneratorPset = cms.PSet(
    DivPSet = cms.PSet(
        nBinsCurv = cms.uint32(4),
        nBinsEta  = cms.uint32(4),
        nBinsLip  = cms.uint32(4),
        nBinsPhi  = cms.uint32(4),
        nBinsTip  = cms.uint32(4)
    ),
    MinResPSet = cms.PSet(
        resCurv = cms.double(3e-3),
        resEta  = cms.double(1e-2),
        resLip  = cms.double(1.00),
        resPhi  = cms.double(0.31), 
        resTip  = cms.double(0.30)
    ),
    MaxResPSet = cms.PSet(
        resCurv = cms.double(3e-4),
        resEta  = cms.double(1e-3),
        resLip  = cms.double(0.10),
        resPhi  = cms.double(0.03),
        resTip  = cms.double(0.30)
    ),
    HalfTurns = cms.PSet(
        positive = cms.uint32(1),
        negative = cms.uint32(0)
    ),
    RequiredLayers = cms.uint32(3)
)

HoughTransformZoomLevelsv1_barrelonly = cms.PSet(
    #phi,d,kappa,dzdl,z0  
    ZoomLevels = cms.VPSet(
         cms.PSet (levels=cms.vuint32(8, 1, 1, 1, 1)),    
         cms.PSet (levels=cms.vuint32(8, 1, 5, 5, 1)),      
         cms.PSet (levels=cms.vuint32(8, 1, 6, 6, 1)),      
         cms.PSet (levels=cms.vuint32(8, 1, 5, 5, 1)),      
         cms.PSet (levels=cms.vuint32(4, 1, 2, 6, 1))
    )
)

HoughTransformZoomLevelsv2 = cms.PSet(
    #phi,d,kappa,dzdl,z0  
    ZoomLevels = cms.VPSet(
        cms.PSet (levels=cms.vuint32(8, 1, 1, 1, 1)),
        cms.PSet (levels=cms.vuint32(5, 1, 5, 5, 1)),
        cms.PSet (levels=cms.vuint32(1, 1, 1, 3, 1)),
        cms.PSet (levels=cms.vuint32(1, 1, 1, 2, 1)),
        cms.PSet (levels=cms.vuint32(3, 1, 5, 5, 1)),
        cms.PSet (levels=cms.vuint32(1, 1, 1, 3, 1)),
        cms.PSet (levels=cms.vuint32(1, 1, 1, 2, 1)),
        cms.PSet (levels=cms.vuint32(3, 1, 2, 5, 1)),
        cms.PSet (levels=cms.vuint32(1, 1, 1, 3, 1)),
        cms.PSet (levels=cms.vuint32(1, 1, 1, 2, 1)),
        cms.PSet (levels=cms.vuint32(3, 1, 3, 5, 1)),
        cms.PSet (levels=cms.vuint32(3, 1, 3, 5, 1)),
        cms.PSet (levels=cms.vuint32(3, 1, 3, 5, 1))
    )
)

HoughTransformFullRegion = cms.PSet(
    HoughTransformRegion = cms.PSet(
        pTmin = cms.double(0.8),
        SV = cms.double(0.5),
        phimin = cms.double(0),
        phimax = cms.double(2*pi),
        dzdlmin = cms.double(-1),
        dzdlmax = cms.double(1)
    )
)


import RecoTracker.TkSeedGenerator.GlobalSeedsFromTriplets_cff
import RecoTracker.TkSeedingLayers.PixelLayerTriplets_cfi
houghTransformStepSeeds = RecoTracker.TkSeedGenerator.GlobalSeedsFromTriplets_cff.globalSeedsFromTriplets.clone()
houghTransformStepSeeds.RegionFactoryPSet.RegionPSet.ptMin = 0.6
houghTransformStepSeeds.RegionFactoryPSet.RegionPSet.originHalfLength = 25.0
houghTransformStepSeeds.RegionFactoryPSet.RegionPSet.originRadius = 2.5
houghTransformStepSeeds.OrderedHitsFactoryPSet = cms.PSet(
        HoughTransformSeedGeneratorPset,
        HoughTransformZoomLevelsv2,
        HoughTransformFullRegion,
        ComponentName = cms.string('HTTripletGenerator'),
        SeedingLayers = cms.string('houghTransformSeedLayersPixelAndOuterMatched'),
        maxElement = cms.uint32(1000000),#100000
        VertexSrc = cms.string('hiSelectedVertex'),
        SeedsFromHits = cms.bool(True), #False - from tracks
        SeedSrc = cms.string('hiSecondPixelTripletSeeds'),
        SeedRefitterOnly = cms.bool(False),
        HoughTransformSeeding = cms.bool(False),
        DebugOutput = cms.bool(False)
    )
houghTransformStepSeeds.SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('PixelClusterShapeSeedComparitor'),
        FilterAtHelixStage = cms.bool(False),
        FilterPixelHits = cms.bool(False),
        FilterStripHits = cms.bool(False),
        ClusterShapeHitFilterName = cms.string('ClusterShapeHitFilter')
    )
houghTransformStepSeeds.SeedCreatorPSet.ComponentName = 'SeedFromConsecutiveHitsTripletOnlyCreator'


#otherwise it is broken on high multiplicity events
#houghTransformStepSeeds.ClusterCheckPSet.doClusterCheck = cms.bool(False)
                            

