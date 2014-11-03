import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('MetFilter')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
		    limit = cms.untracked.int32(-1)
				)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        # 'file:/afs/cern.ch/user/y/yoshin/work/public/MetTesting/2014-10-13v2/outStdMET_CHS.root'
	 "file:/afs/cern.ch/user/y/yoshin/work/public/relval/7E674AE9-933F-E411-9C63-0026189438B1.root"
    )
)

process.dummy = cms.EDAnalyzer("Dummy")

process.demo = cms.EDFilter('MetFilter'
)

process.out = cms.OutputModule(
    "PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = cms.untracked.vstring('drop *_*_*_*', 'keep pfMet*',),
    fileName = cms.untracked.string('test.root')
)

process.p = cms.Path(process.demo)
