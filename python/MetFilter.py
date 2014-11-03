import FWCore.ParameterSet.Config as cms

process = cms.Process('SKIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')


## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")



### =========== Global configuration ==========================

## MET flavor configuration -> may be time consuming
## NoPu and MVA met need to be generated with some refeernce objects (e;g. leptons)
## look for the corresponding area in the config file to set your own definition


### ===========================================================
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

ZMM_720p6_PU50 = cms.untracked.vstring(
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre6/RelValZMM_13/GEN-SIM-RECO/PU50ns_PRE_LS172_V12-v1/00000/025E473F-6442-E411-BDFE-0025905A60BC.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre6/RelValZMM_13/GEN-SIM-RECO/PU50ns_PRE_LS172_V12-v1/00000/2684640B-6342-E411-83F8-003048FFD76E.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre6/RelValZMM_13/GEN-SIM-RECO/PU50ns_PRE_LS172_V12-v1/00000/28AC76E7-7242-E411-9BEF-0026189438DF.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre6/RelValZMM_13/GEN-SIM-RECO/PU50ns_PRE_LS172_V12-v1/00000/905DBDE6-7242-E411-9EC6-002618943950.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre6/RelValZMM_13/GEN-SIM-RECO/PU50ns_PRE_LS172_V12-v1/00000/E44FB840-6442-E411-8AFC-0025905A60DA.root",
)

ZMM_720p4_PU25 = cms.untracked.vstring(
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre4/RelValZMM_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3-v3/00000/221E98F3-C627-E411-B9FD-0025905B8572.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre4/RelValZMM_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3-v3/00000/B8FDC5D5-C327-E411-BCDD-003048FFCB9E.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre4/RelValZMM_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3-v3/00000/DE123EAA-C827-E411-ACA6-002590596484.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre4/RelValZMM_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3-v3/00000/E878E0E5-C527-E411-9FEF-0025905A60E0.root",
)

ZEE_720p4_PU25 = cms.untracked.vstring(
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre4/RelValZEE_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3-v3/00000/22186B06-C427-E411-8A04-002618943940.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre4/RelValZEE_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3-v3/00000/96625DF3-C627-E411-B609-0025905B8572.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre4/RelValZEE_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3-v3/00000/ACECF575-C527-E411-A2F0-002590596486.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre4/RelValZEE_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3-v3/00000/BE29946E-C127-E411-9C4C-00259059391E.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre4/RelValZEE_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3-v3/00000/E4497BCC-C227-E411-9112-0026189437FD.root",
)

process.source = cms.Source("PoolSource",
    fileNames = ZMM_720p4_PU25,
    skipEvents = cms.untracked.uint32(0)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

##
## To remove the "begin job processing line " from printing
##
#process.MessageLogger = cms.Service("MessageLogger",
#   destinations = cms.untracked.vstring('cout',),
#   cout         = cms.untracked.PSet(threshold = cms.untracked.string('ERROR'))
#)
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.default.limit = 10
#process.options.wantSummary = False


#
# for debugging
# Needed by PAT which want an output module
process.dummy = cms.EDAnalyzer("Dummy")

#
# out module for tests
#
# process.out = cms.OutputModule(
#     "PoolOutputModule",
#     splitLevel = cms.untracked.int32(0),
#     outputCommands = cms.untracked.vstring('keep *_*_*_*',),
#     fileName = cms.untracked.string('output/outStdMET_CHS.root')
# )

# PAT
# Use PF2PAT to use PF as input to PAT instead of standard RECO
# from PhysicsTools.PatAlgos.patTemplate_cfg import *
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.pfTools import *


process.load("JetMETCorrections.Type1MET.correctedMet_cff")
process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType1Type2_cff")



process.load("JetMETCorrections.Configuration.JetCorrectionServices_cff")
process.ak4PFCHSL1Fastjet.algorithm = cms.string('AK4PFchs')
process.ak4PFCHSL2Relative.algorithm = cms.string('AK4PFchs')
process.ak4PFCHSL3Absolute.algorithm = cms.string('AK4PFchs')

process.corrPfMetType1CHS = process.corrPfMetType1.clone()
process.corrPfMetType1CHS.src = cms.InputTag('ak4PFJetsCHS')
process.corrPfMetType1CHS.offsetCorrLabel = cms.string("ak4PFCHSL1Fastjet")
process.corrPfMetType1CHS.jetCorrLabel = cms.string("ak4PFCHSL1FastL2L3")
process.pfMetT1CHS = process.pfMetT1.clone(
		srcCorrections = cms.VInputTag(
				cms.InputTag('corrPfMetType1CHS', 'type1')
		)
)

#
# PATH definition , define which filter and sequence will be used before the NtupleMaker
#

process.metFilter = cms.EDFilter('MetFilter'
)

process.p = cms.Path(
    process.correctionTermsPfMetType1Type2*
    process.corrPfMetType1CHS*
    process.pfMetT1*
    process.pfMetT1CHS*
		process.patDefaultSequence*
		process.metFilter
   #pat sequence
#   process.pat_sequence
)

process.out = cms.OutputModule(
    "PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = cms.untracked.vstring(
			# 'drop *_*_*_*',
			# 'keep recoGenMETs_*_*_*', 
			# 'keep recoCaloMETs_*_*_*', 
			'keep recoPFMETs_*_*_*', 
			'keep CorrMETData_*_*_*', 
			'keep *_cleanPatMuons_*_*', 
			'keep *_ak4PFJets_*_*', 
			'keep *_ak4PFJetsCHS_*_*', 
			'keep *_offlinePrimaryVertices_*_*', 
			# 'keep CorrMETData_*_*_*', 
			# Keep anything produced by metFilter
			'keep *_metFilter_*_*',
			),
    # outputCommands = cms.untracked.vstring('keep *',),
    fileName = cms.untracked.string('output/MetFilter_ZMM_720p4_PU25.root')
)


# storage
process.outpath = cms.EndPath(process.out) #dummy

  
