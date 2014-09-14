from step2_RAW2DIGI_L1Reco_RECO import *

process.Timing = cms.Service("Timing")
process.SimpleMemoryCheck=cms.Service("SimpleMemoryCheck", oncePerEventMode=cms.untracked.bool(False))
#input

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.source.fileNames = cms.untracked.vstring('file:HIN-HiFall13DR53X-00001_step1_10_1_Glb.root') 
process.source.skipEvents = cms.untracked.uint32(0)

#output
outputfilename='BJet_STRECO.root'

process.RECODEBUGoutput.outputCommands = cms.untracked.vstring('keep *')
process.RECODEBUGoutput.fileName = cms.untracked.string(outputfilename)


#change low-pt tracking to high-displacement tracking

sv = 0.5

process.hiSecondPixelTripletSeeds.RegionFactoryPSet.RegionPSet.ptMin = cms.double(0.9)
process.hiSecondPixelTripletSeeds.RegionFactoryPSet.RegionPSet.originRadius = cms.double(sv)
#?process.hiSecondPixelTripletSeeds.RegionFactoryPSet.RegionPSet.originHalfLength = 1.0                                                                                            
process.hiSecondPixelTripletSeeds.RegionFactoryPSet.RegionPSet.fixedError = cms.double(sv)
process.hiSecondPixelTripletStepSelector.trackSelectors[0].keepAllTracks = cms.bool(True)

#doesn't work well, there are tracks with such a big momentum                                                                                                                                             
#process.hiPixelPairSeeds.RegionFactoryPSet.RegionPSet.ptMin = 999                                                                                                                                        

#check if the following brings back the highPurity!                                                                                                                                                       
process.hiTracking.remove(process.hiPixelPairClusters)
process.hiTracking.remove(process.hiPixelPairSeeds)
process.hiTracking.remove(process.hiPixelPairTrackCandidates)
process.hiTracking.remove(process.hiPixelPairGlobalPrimTracks)
process.hiTracking.remove(process.hiPixelPairStepSelector)
process.hiGeneralTracks.selectedTrackQuals = cms.VInputTag(cms.InputTag("hiInitialStepSelector","hiInitialStep"), cms.InputTag("hiSecondPixelTripletStepSelector","hiSecondPixelTripletStep"))
process.hiGeneralTracks.TrackProducers = cms.VInputTag(cms.InputTag("hiGlobalPrimTracks"), cms.InputTag("hiSecondPixelTripletGlobalPrimTracks"))

process.hiSecondPixelTripletStepSelector.trackSelectors[0].d0_par2=cms.vdouble(99999.0,99999.0)
process.hiSecondPixelTripletStepSelector.trackSelectors[1].d0_par2=cms.vdouble(99999.0,99999.0)
process.hiSecondPixelTripletStepSelector.trackSelectors[2].d0_par2=cms.vdouble(99999.0,99999.0)
