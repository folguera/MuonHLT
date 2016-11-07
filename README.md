# MuonHLT

Some example and ready-to-run config files are on the test direcory: 

## Getting a new configuration: 
Follow information here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGlobalHLT#Preparing_a_working_area

hltGetConfiguration /dev/CMSSW_8_0_0/GRun --full --offline --mc --unprescale --process SFHLT --globaltag auto:run2_mc_GRun --paths HLTriggerFirstPath,HLT_Mu50_v*,HLTriggerFinalPath --input file:/afs/cern.ch/work/f/folguera/trees/2CE425FC-3D11-E611-BE83-02163E01356D.root --output minimal [--l1-emulator 'Full'] > myHLT_MC.py


