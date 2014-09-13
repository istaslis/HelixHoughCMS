import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import sys

#import glob
#files = ['file:'+x for x in glob.glob('ppPythia/*.root')]
#for x in files: print x

options = VarParsing.VarParsing ()

options.register ('hi',
                  True, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.bool,          # string, int, or float
                  "HI ana (True, default) or pp analyzer (False)")
#if len(sys.argv)>1:
if __name__ == '__main__':
    options.parseArguments()
    print "Arguments parsed"

process = cms.Process("Analyzer")

hiReco = options.hi
HT=False

print ("HI" if hiReco else "pp") + " analysis"

outputfilename = 'Pythia_BJet_Mix_Pt80_TuneZ2_2760GeV_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_L1Reco_HIHTRECO_ANA.root'

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# for writing TFiles
process.load("CommonTools.UtilAlgos.TFileService_cfi")

### standard includes
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('MLoVetere.HTTrackSeeding.HoughTransformAloneStep_cfi')

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        'file:Pythia_BJet_Mix_Pt80_TuneZ2_2760GeV_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_L1Reco_HIHTRECO.root' if hiReco else  'file:pp_step3_RAW2DIGI_L1Reco_RECO.root'
        ),
                            skipEvents = cms.untracked.uint32(0),

 			    dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
                            inputCommands=cms.untracked.vstring(
                                                                'keep *',
                                                                'drop *_*voronoi*_*_*',
#								'drop *_hi*_*_*'
                                                                )
                            )


process.TFileService = cms.Service("TFileService",fileName=cms.string(
                                        outputfilename if hiReco else 'trackAnalysis_pp.root'
   				  ))



#GLOBALTAG
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'STARTHI53_LV1::All' # if hiReco else 'START53_V27::All'

 
from RecoLocalTracker.Configuration.RecoLocalTracker_cff import *
from RecoTracker.IterativeTracking.iterativeTk_cff import *
from MLoVetere.HTTrackSeeding.HoughTransformAloneStep_cfi import *

#process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger = cms.Service("MessageLogger", 
#       destinations   = cms.untracked.vstring('cout'),
#       categories     = cms.untracked.vstring('HTTrackSeeding'),
#       debugModules   = cms.untracked.vstring('houghTransformStepSeeds'),
#       cout           = cms.untracked.PSet(
#                        threshold  = cms.untracked.string('INFO'),
#                        INFO           = cms.untracked.PSet ( limit = cms.untracked.int32(0) ),
#                        DEBUG          = cms.untracked.PSet ( limit = cms.untracked.int32(0) ),
#                       HTTrackSeeding = cms.untracked.PSet ( limit = cms.untracked.int32(10000000) ),
#       )
#)

# standard geometry stuff
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')

# needed to re-do tracker rechits (not stored in reco event content)
process.load('RecoHI.HiTracking.LowPtTracking_PbPb_cff')

# some stuff to associate simulated tracks
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi") #Chi2
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")

# track analyzer
process.load("MitHig.PixelTrackletAnalyzer.trackAnalyzer_cff")
process.load("MitHig.PixelTrackletAnalyzer.recHitAnalyzer_cff")

#process.load("CmsHi/HiHLTAlgos.hievtanalyzer_cfi")
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")

#load centrality
from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import *
overrideCentrality(process)
process.HeavyIonGlobalParameters = cms.PSet(
  centralityVariable = cms.string("HFtowers"),
  nonDefaultGlauberModel = cms.string("Hydjet_Drum"),
  centralitySrc = cms.InputTag("hiCentrality")
  )


process.anaTrack.tpFakeSrc = cms.untracked.InputTag('mergedtruth','MergedTrackTruth')
process.anaTrack.tpEffSrc = cms.untracked.InputTag('mergedtruth','MergedTrackTruth')

process.anaTrack.doSimTrack = True
process.anaTrack.fillSimTrack = True
process.anaTrack.doTrackHits = True
process.anaTrack.useCentrality = hiReco #True if hiReco

process.anaTrack.vertexSrc = ['hiSelectedVertex' if hiReco else 'offlinePrimaryVertices']
process.anaTrack.SeedSrc = ''#'hiPixelTrackSeeds' if hiReco else 'initialStepSeeds'
process.anaTrack.trackSrc = 'hiSelectedTracks' if hiReco else 'generalTracks'

process.anaTrackSecond = process.anaTrack.clone()
process.TrackRefitterSecond = process.TrackRefitter.clone()
process.TrackRefitterSecond.src = 'hiSecondQual'
process.anaTrackSecond.trackSrc = 'TrackRefitterSecond'


process.TrackRefitter.src = 'hiGeneralTracks' if hiReco else 'generalTracks'#cms.InputTag("hiSelectedTracks" if hiReco else "generalTracks") #"houghTransformStepTracks"
process.anaTrackGeneral = process.anaTrack.clone()
process.anaTrackGeneral.SeedSrc = 'houghTransformStepSeeds' if HT else 'hiSecondPixelTripletSeeds'
process.anaTrackGeneral.trackSrc = 'TrackRefitter'


from SimTracker.TrackHistory.CategorySelectors_cff import *
process.hiTrackingCategorySelector = TrackingParticleCategorySelector(
src = cms.InputTag('mergedtruth','MergedTrackTruth'),
cut = cms.string("is('Bottom')")
)

genTag = "hiSignal"#"generator"

process.hiTrackingCategorySelector.trackProducer = 'hiGeneralTracks' if hiReco else 'generalTracks'
process.hiTrackingCategorySelector.hepMC = cms.untracked.InputTag(genTag)

process.bWeakDecaySelector = process.hiTrackingCategorySelector.clone(
    cut = cms.string("is('BWeakDecay')")
)

process.bWeakDecayTracks = process.anaTrackGeneral.clone(
    tpFakeSrc = 'bWeakDecaySelector',
    tpEffSrc = 'bWeakDecaySelector'
)


process.rechits = cms.Sequence(
    process.siPixelRecHits * process.siStripMatchedRecHits * 
    process.anaRecHit)

process.Path = cms.Path(process.rechits
#                        *process.anaTrack
#                        *process.TrackRefitter
#                        *process.anaTrackGeneral
                        *process.TrackRefitterSecond
                        *process.anaTrackSecond
#                        *process.bWeakDecaySelector
#                        *process.bWeakDecayTracks
                        )



process.Timing = cms.Service("Timing")
process.SimpleMemoryCheck=cms.Service("SimpleMemoryCheck",
                                      oncePerEventMode=cms.untracked.bool(False))

