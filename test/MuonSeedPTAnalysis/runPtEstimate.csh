#! /bin/csh
cd /afs/cern.ch/user/d/dfortin/scratch0/test/CMSSW_1_2_0/src/RecoMuon/SeedGenerator/test/MuonSeedPTAnalysis
eval `scramv1 runtime -csh`
cmsRun -p estimatePt.cfg
rfcp /tmp/pt_estimate_5-200.root /castor/cern.ch/user/d/dfortin/pt_estimate_5-200.root
