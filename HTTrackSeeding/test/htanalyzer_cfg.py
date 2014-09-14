from analyzer_cfg import *

outputfilename = 'BJet_HTRECO_ANA.root'
process.source.fileNames = cms.untracked.vstring('file:BJet_HTRECO.root')
process.TFileService.fileName=cms.string(outputfilename)
