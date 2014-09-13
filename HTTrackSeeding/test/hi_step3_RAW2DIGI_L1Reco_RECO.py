# Auto generated configuration file
# using: 
# Revision: 1.381.2.28 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: step3 --conditions STARTHI53_LV1::All -s RAW2DIGI,L1Reco,RECO --scenario HeavyIons --datatier GEN-SIM-RECO --eventcontent RECODEBUG --no_exec --filein file:Pythia_BJet_Mix_Pt80_TuneZ2_2760GeV.root
import FWCore.ParameterSet.Config as cms

process = cms.Process('HIRECO')

import glob

#files = ['file:'+x for x in glob.glob('ppPythia/*.root')]
#print files

import sys


rawData = True


process.Timing = cms.Service("Timing")
process.SimpleMemoryCheck=cms.Service("SimpleMemoryCheck", oncePerEventMode=cms.untracked.bool(False))

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("SimGeneral.TrackingAnalysis.Playback_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

#from bJet_2760GeV_RAW_100_INPUTFILES import *
#process.source = source

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
        #'file:Pythia_BJet_Mix_Pt80_TuneZ2_2760GeV.root') #53x file
        #'file:FCFA630E-5110-E311-80EC-00266CFAE818.root') #44x file
        'file:Pythia_BJet_Mix_Pt80_TuneZ2_2760GeV_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_L1Reco.root') #53x L1Reco file

)

outputfilename='Pythia_BJet_Mix_Pt80_TuneZ2_2760GeV_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_L1Reco_HIRECO.root'

process.options = cms.untracked.PSet()

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.381.2.28 $'),
    annotation = cms.untracked.string('step3 nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)

sv = 0.5

process.hiSecondPixelTripletSeeds.RegionFactoryPSet.RegionPSet.ptMin = cms.double(0.9)
process.hiSecondPixelTripletSeeds.RegionFactoryPSet.RegionPSet.originRadius = cms.double(sv)
#process.hiSecondPixelTripletSeeds.RegionFactoryPSet.RegionPSet.originHalfLength = 1.0
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


#for standard seeding use pixel seeding only (no EC)                                                                                                                                                    
#process.hiPixel3PrimTracks.OrderedHitsFactoryPSet.SeedingLayers = cms.string("houghTransformSeedLayersPixelBarrelOnly")                                                     
#process.hiSecondPixelTripletSeeds.OrderedHitsFactoryPSet.SeedingLayers = cms.string("houghTransformSeedLayersPixelBarrelOnly")                                                                            



# Output definition

process.RECODEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = cms.untracked.vstring('keep *'),#('drop *', 'keep *_*_*_SIM', 'keep *_mergedtruth_*_*', 'keep *Si*Cluster*_*_*_*', 'keep *Centrality*_*_*_*', 'keep *_*_*_HIRECO'),#process.RECODEBUGEventContent.outputCommands,#
    fileName = cms.untracked.string(outputfilename),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
    )
)


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'STARTHI53_LV1::All', '')

process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECODEBUGoutput_step = cms.EndPath(process.RECODEBUGoutput)

process.RawReco = cms.Path(process.RawToDigi+process.L1Reco+
                           process.reconstructionHeavyIons)
process.ReReco = cms.Path(process.siPixelRecHits+process.siStripMatchedRecHits+process.hiTracking)

process.reco = process.RawReco if rawData else process.ReReco

# Schedule definition
process.schedule = cms.Schedule(process.reco,process.endjob_step,process.RECODEBUGoutput_step)
