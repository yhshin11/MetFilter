import FWCore.ParameterSet.Config as cms
# Command line argument parsing
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process('SKIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')


## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

## Options and Output Report
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False),
																			SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# setup any defaults you want
options.maxEvents = -1 # -1 means all events
options.inputFiles= 'root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre6/RelValZMM_13/GEN-SIM-RECO/PU50ns_PRE_LS172_V12-v1/00000/E44FB840-6442-E411-8AFC-0025905A60DA.root'
options.outputFile = 'output/test.root'
options.register ('test',
                  0, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Testing option. Set to 1 to enable testing options.")
options.test = 0
options.register ('dataset',
                  'ZMM_720p6_PU50', # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "Label for dataset. See config file for options.")


# get and parse the command line arguments
options.parseArguments()

if options.test: print """
options.test set to True. Using test settings.
Running over 100 events.
Output file set to output/MetFilter_basicSkim_test.root
"""

# # Unscheduled mode
process.options.allowUnscheduled = cms.untracked.bool(True)



### =========== Global configuration ==========================

## MET flavor configuration -> may be time consuming
## NoPu and MVA met need to be generated with some refeernce objects (e;g. leptons)
## look for the corresponding area in the config file to set your own definition


### ===========================================================
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

# Define dictionary of dataset->filenames
datasets={}
datasets['ZMM_720p6_PU50'] = cms.untracked.vstring(
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre6/RelValZMM_13/GEN-SIM-RECO/PU50ns_PRE_LS172_V12-v1/00000/025E473F-6442-E411-BDFE-0025905A60BC.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre6/RelValZMM_13/GEN-SIM-RECO/PU50ns_PRE_LS172_V12-v1/00000/2684640B-6342-E411-83F8-003048FFD76E.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre6/RelValZMM_13/GEN-SIM-RECO/PU50ns_PRE_LS172_V12-v1/00000/28AC76E7-7242-E411-9BEF-0026189438DF.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre6/RelValZMM_13/GEN-SIM-RECO/PU50ns_PRE_LS172_V12-v1/00000/905DBDE6-7242-E411-9EC6-002618943950.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre6/RelValZMM_13/GEN-SIM-RECO/PU50ns_PRE_LS172_V12-v1/00000/E44FB840-6442-E411-8AFC-0025905A60DA.root",
)
datasets['ZMM_720p6_PU25'] = cms.untracked.vstring(
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre6/RelValZMM_13/GEN-SIM-RECO/PU25ns_PRE_LS172_V11-v1/00000/20D22466-5F42-E411-AA3A-0026189438B1.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre6/RelValZMM_13/GEN-SIM-RECO/PU25ns_PRE_LS172_V11-v1/00000/72FE953D-6142-E411-9511-002618943886.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre6/RelValZMM_13/GEN-SIM-RECO/PU25ns_PRE_LS172_V11-v1/00000/7C26F64C-6F42-E411-9978-002618943981.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre6/RelValZMM_13/GEN-SIM-RECO/PU25ns_PRE_LS172_V11-v1/00000/C65D2806-7142-E411-8352-0025905A610A.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre6/RelValZMM_13/GEN-SIM-RECO/PU25ns_PRE_LS172_V11-v1/00000/D4310C16-5F42-E411-AD54-002618943925.root",
)
datasets['ZMM_720p6_PU0'] = cms.untracked.vstring(
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre6/RelValZMM_13/GEN-SIM-RECO/PRE_LS172_V11-v1/00000/7E674AE9-933F-E411-9C63-0026189438B1.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre6/RelValZMM_13/GEN-SIM-RECO/PRE_LS172_V11-v1/00000/D0DE9EE9-933F-E411-8F7F-0025905A60DE.root",
)

datasets['ZMM_720p4_PU25'] = cms.untracked.vstring(
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre4/RelValZMM_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3-v3/00000/221E98F3-C627-E411-B9FD-0025905B8572.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre4/RelValZMM_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3-v3/00000/B8FDC5D5-C327-E411-BCDD-003048FFCB9E.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre4/RelValZMM_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3-v3/00000/DE123EAA-C827-E411-ACA6-002590596484.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre4/RelValZMM_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3-v3/00000/E878E0E5-C527-E411-9FEF-0025905A60E0.root",
)

datasets['ZEE_720p4_PU25'] = cms.untracked.vstring(
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre4/RelValZEE_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3-v3/00000/22186B06-C427-E411-8A04-002618943940.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre4/RelValZEE_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3-v3/00000/96625DF3-C627-E411-B609-0025905B8572.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre4/RelValZEE_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3-v3/00000/ACECF575-C527-E411-A2F0-002590596486.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre4/RelValZEE_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3-v3/00000/BE29946E-C127-E411-9C4C-00259059391E.root",
"root://eoscms//eos/cms/store/relval/CMSSW_7_2_0_pre4/RelValZEE_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3-v3/00000/E4497BCC-C227-E411-9112-0026189437FD.root",
)

process.source = cms.Source("PoolSource",
    fileNames = datasets[options.dataset],
    skipEvents = cms.untracked.uint32(0)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )
if options.test: process.maxEvents.input = cms.untracked.int32(100)

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

# PAT
# Use PF2PAT to use PF as input to PAT instead of standard RECO
# from PhysicsTools.PatAlgos.patTemplate_cfg import *
process.load("PhysicsTools.PatAlgos.patSequences_cff")
# Copied from patTuple_PF2PAT_cfg.py
from PhysicsTools.PatAlgos.tools.pfTools import *

# First PF2PAT
postfix = "PFlowNoChs"
jetAlgo="AK4"
pvCollection = cms.InputTag('offlinePrimaryVertices')
# jetCorrections=('AK4PF', ['L1FastJet','L2Relative','L3Absolute'],'None')
usePF2PAT(process, runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=True, postfix=postfix, pvCollection=pvCollection)
switchJetCollection(process,
		postfix = postfix,
		jetSource = cms.InputTag('ak4PFJets'),
		pvSource = pvCollection,
		algo = jetAlgo,
		jetCorrections = ('AK4PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], ''),
		# outputModules = 'out'
)
# Second PF2PAT
postfix = "PFlowChs"
jetAlgo="AK4"
pvCollection = cms.InputTag('offlinePrimaryVertices')
# jetCorrections=('AK4PFchs', ['L1FastJet','L2Relative','L3Absolute'],'None')
usePF2PAT(process, runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=True, postfix=postfix, pvCollection=pvCollection)
switchJetCollection(process,
		postfix = postfix,
		jetSource = cms.InputTag('ak4PFJetsCHS'),
		pvSource = pvCollection,
		algo = jetAlgo,
		jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], ''),
		# outputModules = 'out'
)
# addJetCollection(process,
# 		labelName = 'CHS',
# 		postfix = postfix,
# 		jetSource = cms.InputTag('ak4PFJetsCHS'),
# 		pvSource = pvCollection,
# 		algo = jetAlgo,
# 		jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], ''),
# 		# outputModules = 'out'
# )
# process.cleanPatJetsCHSPFlow = process.cleanPatJetsPFlow.clone(src = 'selectedPatJetsCHSPFlow')
# # top projections in PF2PAT:
# getattr(process,"pfNoPileUpJME"+postfix).enable = True
# getattr(process,"pfNoMuonJME"+postfix).enable = True
# getattr(process,"pfNoElectronJME"+postfix).enable = True
# getattr(process,"pfNoTau"+postfix).enable = False
# getattr(process,"pfNoJet"+postfix).enable = True
# # to use tau-cleaned jet collection uncomment the following:
# #getattr(process,"pfNoTau"+postfix).enable = True
#
# # verbose flags for the PF2PAT modules
# getattr(process,"pfNoMuonJME"+postfix).verbose = False
#
# # enable delta beta correction for muon selection in PF2PAT?
# getattr(process,"pfIsolatedMuons"+postfix).doDeltaBetaCorrection = cms.bool(False)

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
		# process.patDefaultSequence
		# process.cleanPatJetsCHSPFlow
		process.metFilter
   #pat sequence
#   process.pat_sequence
)

# Add PF2PAT output to the created file
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

process.out = cms.OutputModule(
    "PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    # outputCommands = cms.untracked.vstring('keep *',),
    fileName = cms.untracked.string('output/basicSkim_ZMM_720p6_PU0.root')
    # fileName = cms.untracked.string('output/basicSkim_test.root')
)
if options.test: process.out.fileName = cms.untracked.string('output/MetFilter_basicSkim_test.root')
process.out.outputCommands = cms.untracked.vstring(
		# 'keep *',
	'drop *',
	# 'keep recoGenMETs_*_*_*',
	# 'keep recoCaloMETs_*_*_*',
	# 'keep recoPFMETs_*_*_*',
	# 'keep CorrMETData_*_*_*',
	'keep *_pfMetT1_*_*',
	'keep *_corrPfMetType1_*_*',
	'keep *_pfMetT1CHS_*_*',
	'keep *_corrPfMetType1CHS_*_*',
	########
	# Muons
	'keep *_cleanPatMuons_*_*',
	########
	# Jets
	'keep *_selectedPatJets_*_*',
	# 'keep *_selectedPatJetsPFlow_*_*',
	# 'keep *_selectedPatJetsCHSPFlow_*_*',
	'keep *_selectedPatJetsPFlowNoChs_*_*',
	'keep *_selectedPatJetsPFlowChs_*_*',
	'keep *_cleanPatJets_*_*',
	# 'keep *_cleanPatJetsPFlow_*_*',
	# 'keep *_cleanPatJetsCHSPFlow_*_*',
	'keep *_cleanPatJetsPFlowNoChs_*_*',
	'keep *_cleanPatJetsPFlowChs_*_*',
	'keep *_ak4PFJets_*_*',
	'keep *_ak4PFJetsCHS_*_*',
	# 'keep *_offlinePrimaryVertices_*_*',
	# 'keep CorrMETData_*_*_*',
	# Keep anything produced by metFilter
	'keep *_metFilter_*_*',
	# *patEventContent)
)


# load the coreTools of PAT
from PhysicsTools.PatAlgos.tools.coreTools import removeMCMatching
removeMCMatching(process, ['All'], 'PFlowNoChs', ['out'])
removeMCMatching(process, ['All'], 'PFlowChs', ['out'])


# storage
process.outpath = cms.EndPath(process.out) #dummy


