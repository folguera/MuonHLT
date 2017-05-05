# MuonHLT

Some example and ready-to-run config files are on the test direcory: 

## Getting a new configuration: 
Follow information here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGlobalHLT#Preparing_a_working_area

hltGetConfiguration /dev/CMSSW_8_0_0/GRun --full --offline --mc --unprescale --process SFHLT --globaltag auto:run2_mc_GRun --paths HLTriggerFirstPath,HLT_Mu50_v*,HLTriggerFinalPath --input file:/afs/cern.ch/work/f/folguera/trees/2CE425FC-3D11-E611-BE83-02163E01356D.root --output minimal [--l1-emulator 'Full'] > myHLT_MC.py

# At the end of the config file put this: 
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *

process.muonDebugger =cms.EDAnalyzer("MuonHLTDebugger",
                                     MuonServiceProxy,
                                     triggerResults  = cms.untracked.InputTag("TriggerResults::REHLT"),
                                     triggerSummary  = cms.untracked.InputTag("hltTriggerSummaryAOD::REHLT"),
                                     L3Candidates    = cms.untracked.InputTag("hltIterL3MuonCandidates"), 
                                     L2Candidates    = cms.untracked.InputTag("hltL2MuonCandidates"), 
                                     L1Candidates    = cms.untracked.InputTag("hltGtStage2Digis", "Muon"), 
                                     MuonLinksTag    = cms.untracked.InputTag("hltIterL3MuonsLinksCombination","","REHLT"),
                                     genParticlesTag = cms.untracked.InputTag("genParticles"),
                                     muonTag         = cms.untracked.InputTag("muons"),
                                     triggerProcess  = cms.string("REHLT"),
                                     triggerName     = cms.string("HLT_Mu50_v5"),
                                     l1filterLabel   = cms.string("hltL1fL1sMu22Or25L1Filtered0"),
                                     l2filterLabel   = cms.string("hltL2fL1sMu22Or25L1f0L2Filtered10Q"),
                                     l3filterLabel   = cms.string("hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q"),
                                     debuglevel      = cms.untracked.uint32(0),
                                     isMC            = cms.untracked.bool(True),
                                     UseGenInfo      = cms.untracked.bool(True),
                                     runSharedHits   = cms.untracked.bool(False),
                                     )   

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outfileName),
                                   closeFileFast = cms.untracked.bool(False)
                                   )


process.validation = cms.EndPath(
    process.muonDebugger
)
