import FWCore.ParameterSet.Config as cms

houghTransformSeedLayersAllHitsTwoSets = cms.ESProducer("SeedingLayersESProducer",
    ComponentName = cms.string('HoughTransformSeedLayersAllHitsTwoSets'),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3+'
        'FPix1_pos+FPix2_pos+'
        'TIB1+TIB2+TIB3+TIB4+'
        'TID1_pos+TID2_pos+TID3_pos+'
        'TOB1+TOB2+TOB3+TOB4+TOB5+TOB6+'
        'TEC1_pos+TEC2_pos+TEC3_pos+TEC4_pos+TEC5_pos+TEC6_pos+TEC7_pos+TEC8_pos+TEC9_pos'
      , 'BPix1+BPix2+BPix3+'
        'FPix1_neg+FPix2_neg+'
        'TIB1+TIB2+TIB3+TIB4+'
        'TID1_neg+TID2_neg+TID3_neg+'
        'TOB1+TOB2+TOB3+TOB4+TOB5+TOB6+'
        'TEC1_neg+TEC2_neg+TEC3_neg+TEC4_neg+TEC5_neg+TEC6_neg+TEC7_neg+TEC8_neg+TEC9_neg'
    ),
    BPix = cms.PSet(
        HitProducer     = cms.string  ('siPixelRecHits'),
        useErrorsFromParam = cms.bool (True),
        hitErrorRZ      = cms.double  (0.0060),
        hitErrorRPhi    = cms.double  (0.0027),
        TTRHBuilder     = cms.string  ('WithTrackAngle'),
        skipClusters    = cms.InputTag('houghTransformStepClusters' )
    ),
    FPix = cms.PSet(
        HitProducer     = cms.string  ('siPixelRecHits'),
        useErrorsFromParam = cms.bool (True),
        hitErrorRPhi    = cms.double  (0.0051),
        hitErrorRZ      = cms.double  (0.0036),
        TTRHBuilder     = cms.string  ('WithTrackAngle'),
        skipClusters    = cms.InputTag('houghTransformStepClusters' )
    ),
    TIB = cms.PSet(
        matchedRecHits  = cms.InputTag('siStripMatchedRecHits','matchedRecHit'        ),
        rphiRecHits     = cms.InputTag('siStripMatchedRecHits','rphiRecHitUnmatched'  ),
        stereoRecHits   = cms.InputTag('siStripMatchedRecHits','stereoRecHitUnmatched'),
        TTRHBuilder     = cms.string  ('WithTrackAngle'),
        skipClusters    = cms.InputTag('houghTransformStepClusters' )
    ),
    TID = cms.PSet(
        matchedRecHits  = cms.InputTag('siStripMatchedRecHits','matchedRecHit'        ),
        rphiRecHits     = cms.InputTag('siStripMatchedRecHits','rphiRecHitUnmatched'  ),
        stereoRecHits   = cms.InputTag('siStripMatchedRecHits','stereoRecHitUnmatched'),
        TTRHBuilder     = cms.string  ('WithTrackAngle'),
        skipClusters    = cms.InputTag('houghTransformStepClusters' )
    ),
    TOB = cms.PSet(
        matchedRecHits  = cms.InputTag('siStripMatchedRecHits','matchedRecHit'        ),
        rphiRecHits     = cms.InputTag('siStripMatchedRecHits','rphiRecHitUnmatched'  ),
        stereoRecHits   = cms.InputTag('siStripMatchedRecHits','stereoRecHitUnmatched'),
        TTRHBuilder     = cms.string  ('WithTrackAngle'),
        skipClusters    = cms.InputTag('houghTransformStepClusters' )
    ),
    TEC = cms.PSet(
        matchedRecHits  = cms.InputTag('siStripMatchedRecHits','matchedRecHit'        ),
        rphiRecHits     = cms.InputTag('siStripMatchedRecHits','rphiRecHitUnmatched'  ),
        stereoRecHits   = cms.InputTag('siStripMatchedRecHits','stereoRecHitUnmatched'),
        TTRHBuilder     = cms.string  ('WithTrackAngle'),
        skipClusters    = cms.InputTag('houghTransformStepClusters' )
    )
)

houghTransformSeedLayersAllHitsOneSet = houghTransformSeedLayersAllHitsTwoSets.clone( 
    ComponentName = cms.string('HoughTransformSeedLayersAllHitsOneSet'),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3+'
        'FPix1_pos+FPix2_pos+'
        'FPix1_neg+FPix2_neg+'
        'TIB1+TIB2+TIB3+TIB4+'
        'TID1_neg+TID2_neg+TID3_neg+'
        'TID1_pos+TID2_pos+TID3_pos+'
        'TOB1+TOB2+TOB3+TOB4+TOB5+TOB6+'
        'TEC1_pos+TEC2_pos+TEC3_pos+TEC4_pos+TEC5_pos+TEC6_pos+TEC7_pos+TEC8_pos+TEC9_pos+'
        'TEC1_neg+TEC2_neg+TEC3_neg+TEC4_neg+TEC5_neg+TEC6_neg+TEC7_neg+TEC8_neg+TEC9_neg'
    ),
)

houghTransformSeedLayersAllHitsBarrelOnly = houghTransformSeedLayersAllHitsTwoSets.clone(
    ComponentName = cms.string('HoughTransformSeedLayersAllHitsBarrelOnly'),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3+'
        'TIB1+TIB2+TIB3+TIB4+'
        'TOB1+TOB2+TOB3+TOB4+TOB5+TOB6'
    ),
)

houghTransformSeedLayersPixelAndMatchedHitsTwoSets = cms.ESProducer("SeedingLayersESProducer",
    ComponentName = cms.string('HoughTransformSeedLayersPixelAndMatchedHitsTwoSets'),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3+'
        'FPix1_pos+FPix2_pos+'
        'TIB1+TIB2+'
        'TID1_pos+TID2_pos+TID3_pos+'
        'TOB1+TOB2+'
        'TEC1_pos+TEC2_pos+TEC3_pos+TEC4_pos+TEC5_pos+TEC6_pos+TEC7_pos+TEC8_pos+TEC9_pos'
      , 'BPix1+BPix2+BPix3+'
        'FPix1_neg+FPix2_neg+'
        'TIB1+TIB2+'
        'TID1_neg+TID2_neg+TID3_neg+'
        'TOB1+TOB2+'
        'TEC1_neg+TEC2_neg+TEC3_neg+TEC4_neg+TEC5_neg+TEC6_neg+TEC7_neg+TEC8_neg+TEC9_neg'
    ),
    BPix = cms.PSet(
        HitProducer     = cms.string  ('siPixelRecHits'),
        useErrorsFromParam = cms.bool (True),
        hitErrorRZ      = cms.double  (0.0060),
        hitErrorRPhi    = cms.double  (0.0027),
        TTRHBuilder     = cms.string  ('WithTrackAngle'),
        skipClusters    = cms.InputTag('houghTransformStepClusters' )
    ),
    FPix = cms.PSet(
        HitProducer     = cms.string  ('siPixelRecHits'),
        useErrorsFromParam = cms.bool (True),
        hitErrorRPhi    = cms.double  (0.0051),
        hitErrorRZ      = cms.double  (0.0036),
        TTRHBuilder     = cms.string  ('WithTrackAngle'),
        skipClusters    = cms.InputTag('houghTransformStepClusters' )
    ),
    TIB = cms.PSet(
        matchedRecHits  = cms.InputTag('siStripMatchedRecHits','matchedRecHit'),
        TTRHBuilder     = cms.string  ('WithTrackAngle'),
        skipClusters    = cms.InputTag('houghTransformStepClusters' )
    ),
    TID = cms.PSet(
        matchedRecHits  = cms.InputTag('siStripMatchedRecHits','matchedRecHit'),
        useRingSlector  = cms.bool    (True),
        minRing         = cms.int32   (1),
        maxRing         = cms.int32   (2),
        TTRHBuilder     = cms.string  ('WithTrackAngle'),
        skipClusters    = cms.InputTag('houghTransformStepClusters' )
    ),
    TOB = cms.PSet(
        matchedRecHits  = cms.InputTag('siStripMatchedRecHits','matchedRecHit'),
        TTRHBuilder     = cms.string  ('WithTrackAngle'),
        skipClusters    = cms.InputTag('houghTransformStepClusters' )
    ),
    TEC = cms.PSet(
        matchedRecHits  = cms.InputTag('siStripMatchedRecHits','matchedRecHit'),
        useRingSlector  = cms.bool    (True),
        minRing         = cms.int32   (1),
        maxRing         = cms.int32   (5),
        TTRHBuilder     = cms.string  ('WithTrackAngle'),
        skipClusters    = cms.InputTag('houghTransformStepClusters' )
    )
)
 
houghTransformSeedLayersPixelAndMatchedHitsOneSet = houghTransformSeedLayersPixelAndMatchedHitsTwoSets.clone(
    ComponentName = cms.string('HoughTransformSeedLayersPixelAndMatchedHitsOneSet'),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3+'
        'FPix1_pos+FPix2_pos+'
        'FPix1_neg+FPix2_neg+'
        'TIB1+TIB2+'
        'TID1_pos+TID2_pos+TID3_pos+'
        'TID1_neg+TID2_neg+TID3_neg+'
        'TOB1+TOB2+'
        'TEC1_pos+TEC2_pos+TEC3_pos+TEC4_pos+TEC5_pos+TEC6_pos+TEC7_pos+TEC8_pos+TEC9_pos+' 
        'TEC1_neg+TEC2_neg+TEC3_neg+TEC4_neg+TEC5_neg+TEC6_neg+TEC7_neg+TEC8_neg+TEC9_neg'
    )
)

houghTransformSeedLayersPixelAndMatchedHitsBarrelOnly = houghTransformSeedLayersPixelAndMatchedHitsTwoSets.clone(
    ComponentName = cms.string('HoughTransformSeedLayersPixelAndMatchedHitsBarrelOnly'),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3+'
        'TIB1+TIB2+'
        'TOB1+TOB2'
    )
)

