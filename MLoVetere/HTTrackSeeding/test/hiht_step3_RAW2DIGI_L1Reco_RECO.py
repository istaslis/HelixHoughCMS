import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from hi_step3_RAW2DIGI_L1Reco_RECO import *

process.load('MLoVetere.HTTrackSeeding.HoughTransformAloneStep_cfi')

outputfilename = 'Pythia_BJet_Mix_Pt80_TuneZ2_2760GeV_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_L1Reco_HIHTRECO.root'

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('file:Pythia_BJet_Mix_Pt80_TuneZ2_2760GeV_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_L1Reco.root'),#Pythia_BJet_Mix_Pt80_TuneZ2_2760GeV.root')
    skipEvents = cms.untracked.uint32(0)#9#7#32)
)

process.RECODEBUGoutput.fileName = outputfilename
#process.RECODEBUGoutput.outputCommands=cms.untracked.vstring('drop *')


from MLoVetere.HTTrackSeeding.HoughTransformAloneStep_cfi import *

houghTransformStepSeeds.OrderedHitsFactoryPSet.SeedsFromHits = True
houghTransformStepSeeds.OrderedHitsFactoryPSet.SeedSrc = cms.string('hiSecondPixelTripletSeeds')
houghTransformStepSeeds.OrderedHitsFactoryPSet.VertexSrc = 'hiSelectedVertex'

process.hiTracking.replace(process.hiSecondPixelTripletSeeds,process.hiSecondPixelTripletSeeds+process.houghTransformStepSeeds)


process.hiSecondPixelTripletTrackCandidates.src = cms.InputTag("houghTransformStepSeeds")
houghTransformStepSeeds.ClusterCheckPSet.doClusterCheck = cms.bool(False)


