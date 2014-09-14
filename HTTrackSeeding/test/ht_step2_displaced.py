import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from hi_step2_displaced import *

outputfilename = 'BJet_HTRECO.root'

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.source.fileNames = cms.untracked.vstring('file:HIN-HiFall13DR53X-00001_step1_10_1_Glb.root')
process.source.skipEvents = cms.untracked.uint32(0)

process.RECODEBUGoutput.fileName = outputfilename
process.RECODEBUGoutput.outputCommands=cms.untracked.vstring('keep *')


#HT setup

process.load('MLoVetere.HTTrackSeeding.HoughTransformSeedFilterStep_cfi')
from MLoVetere.HTTrackSeeding.HoughTransformAloneStep_cfi import *

process.houghTransformStepSeeds.OrderedHitsFactoryPSet.SeedsFromHits = True
process.houghTransformStepSeeds.OrderedHitsFactoryPSet.SeedSrc = cms.string('hiSecondPixelTripletSeeds')
process.houghTransformStepSeeds.OrderedHitsFactoryPSet.VertexSrc = 'hiSelectedVertex'
process.hiSecondPixelTripletTrackCandidates.src = cms.InputTag("houghTransformStepSeeds")
process.houghTransformStepSeeds.ClusterCheckPSet.doClusterCheck = cms.bool(False)

process.hiTracking.replace(process.hiSecondPixelTripletSeeds,process.hiSecondPixelTripletSeeds+process.houghTransformStepSeeds)




