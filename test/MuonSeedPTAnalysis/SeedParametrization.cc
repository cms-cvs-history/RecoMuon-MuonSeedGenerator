// Class Header
#include "SeedParametrization.h"
#include "RecoMuon/MuonSeedGenerator/test/MuonSeedPTAnalysis/SegSelector.h"

// for MuonSeedBuilder
/*
#include "RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h"
#include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
*/

// Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "PluginManager/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "TFile.h"
#include "TVector3.h"

#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <string>
#include <stdio.h>
#include <algorithm>


DEFINE_FWK_MODULE(SeedParametrization);
using namespace std;
using namespace edm;

// constructors
SeedParametrization::SeedParametrization(const ParameterSet& pset){ 

  debug             = pset.getUntrackedParameter<bool>("debug");
  rootFileName      = pset.getUntrackedParameter<string>("rootFileName");
  recHitLabel       = pset.getUntrackedParameter<string>("recHitLabel");
  cscSegmentLabel   = pset.getUntrackedParameter<string>("cscSegmentLabel");
  dtrecHitLabel     = pset.getUntrackedParameter<string>("dtrecHitLabel");
  dtSegmentLabel    = pset.getUntrackedParameter<string>("dtSegmentLabel");
  simHitLabel       = pset.getUntrackedParameter<string>("simHitLabel");
  simTrackLabel     = pset.getUntrackedParameter<string>("simTrackLabel");
  muonseedLabel     = pset.getUntrackedParameter<string>("muonseedLabel");

  recsegSelector    = new SegSelector(pset); 
  //muonSeedBuilder_  = new MuonSeedBuilder( pset );
  //SeedBuilder_      = new SeedBuilder( pset );

  // Create the root file
  theFile = new TFile(rootFileName.c_str(), "RECREATE");
  theFile->mkdir("AllMuonSys");
  theFile->cd();
  theFile->mkdir("CSC_All");
  theFile->cd();
  theFile->mkdir("DT_All");
  theFile->cd();
  theFile->mkdir("ME_All");
  theFile->cd();
  theFile->mkdir("MB_All");
  theFile->cd();
  theFile->mkdir("OL_All");
  // TTree test
  tr_muon = new TNtuple1();

  // The dphi between chamber by chamber
  // All possible segment pair in CSC 
  //                  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 
  int csc1[2][15]= {{11,11,12,12,13,11,11,12,13,11,21,21,22,21,31}
                   ,{12,21,21,22,22,31,32,32,32,41,31,32,32,41,41} }; 
  char ME_nu1[8];
  for (int i =0; i<15; i++) {
      sprintf(ME_nu1,"ME_%d-%d", csc1[0][i] , csc1[1][i] );
      hME1[i] = new H2DRecHit4(ME_nu1);
      cout <<"hME1_"<<i<<" = "<< ME_nu1<<endl;
  }
  // All possible segment pair in DT 
  //                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25  
  int dt1[2][26]={{10,10,10,10,10,11,11,11,11,11,11,12,12,20,20,20,21,21,21,21,22,22,30,31,31,32}
                 ,{20,30,31,40,41,21,22,31,32,41,42,22,32,30,40,41,31,32,41,42,32,42,40,41,42,42}};
  char MB_nu1[8];
  for (int i =0; i<26; i++) {
      sprintf(MB_nu1,"MB_%d-%d", dt1[0][i] , dt1[1][i] );
      hMB1[i] = new H2DRecHit5(MB_nu1);
      cout <<"hMB1_"<<i<<" = "<< MB_nu1<<endl;
  }
  // All possible segment pair in Overlap region between DT and CSC
  //               0  1  2  3  4  5 
  int olp[2][6]={{12,12,12,22,22,32} 
                ,{13,22,32,13,22,13}};
  char OL_nu1[8];
  for (int i =0; i<6; i++) {
      sprintf(OL_nu1,"OL_%d-%d", olp[0][i] , olp[1][i] );
      hOL1[i] = new H2DRecHit10(OL_nu1);
      cout <<"hOL_"<<i<<" = "<< OL_nu1<<endl;
  }

  // All single chamber segment in CSC
  //             0  1  2  3  4  5  6  7
  int csc2[8] ={11,12,13,21,22,31,32,41};
  char ME_nu2[6];
  for (int i =0; i<8; i++) {
      sprintf(ME_nu2,"SME_%d", csc2[i] );
      hME2[i] = new H2DRecHit6(ME_nu2);
      cout <<"hME2_"<<i<<" = "<< ME_nu2<<endl;
  }

  // All single chamber segment in DT
  //           0  1  2  3  4  5  6  7  8  9 10 11 
  int dt2[12] ={10,11,12,20,21,22,30,31,32,40,41,42};
  char MB_nu2[6];
  for (int i =0; i<12; i++) {
      sprintf(MB_nu2,"SMB_%d", dt2[i] );
      hMB2[i] = new H2DRecHit7(MB_nu2);
      cout <<"hMB2_"<<i<<" = "<< MB_nu2<<endl;
  }
  h_all  = new H2DRecHit1("AllMu_");
  h_csc  = new H2DRecHit2("CSC_A");
  h_dt   = new H2DRecHit3("DT_A");

}

// destructor
SeedParametrization::~SeedParametrization(){

  if (debug) cout << "[SeedQualityAnalysis] Destructor called" << endl;
  delete recsegSelector;
  //delete muonSeedBuilder_; 
  //delete SeedBuilder_; 
  // Write the histos to file
  theFile->cd();
  theFile->cd("AllMuonSys");
  h_all->Write();
  theFile->cd();

  theFile->cd("CSC_All");
  h_csc->Write();
  theFile->cd();

  theFile->cd("DT_All");
  h_dt->Write();
  theFile->cd();

  theFile->cd("ME_All");
  for (int i=0; i<15; i++) {
      hME1[i]->Write();
  }
  for (int i=0; i<8; i++) {
      hME2[i]->Write();
  }
  theFile->cd();

  theFile->cd("MB_All");
  for (int i=0; i<26; i++) {
      hMB1[i]->Write();
  }
  for (int i=0; i<12; i++) {
      hMB2[i]->Write();
  }
  theFile->cd();

  theFile->cd("OL_All");
  for (int i=0; i<6; i++) {
      hOL1[i]->Write();
  }
  // for tree
  theFile->cd();
  tr_muon->Write();

  // Release the memory ...
  for (int i=0; i<15; i++) {
      delete hME1[i];
  }
  for (int i=0; i<26; i++) {
      delete hMB1[i];
  }
  for (int i=0; i<6; i++) {
      delete hOL1[i];
  }
  for (int i=0; i<8; i++) {
      delete hME2[i];
  }
  for (int i=0; i<12; i++) {
      delete hMB2[i];
  }
  delete h_all;
  delete h_csc;
  delete h_dt;
  delete tr_muon;

  theFile->Close();
  if (debug) cout << "************* Finished writing histograms to file" << endl;

}

// The Main...Aanlysis...

void SeedParametrization::analyze(const Event& event, const EventSetup& eventSetup)
{
  //Get the CSC Geometry :
  ESHandle<CSCGeometry> cscGeom;
  eventSetup.get<MuonGeometryRecord>().get(cscGeom);

  //Get the DT Geometry :
  ESHandle<DTGeometry> dtGeom;
  eventSetup.get<MuonGeometryRecord>().get(dtGeom);

  // Get the RecHits collection :
  Handle<CSCRecHit2DCollection> csc2DRecHits; 
  event.getByLabel(recHitLabel, csc2DRecHits);

  // Get the CSC Segments collection :
  Handle<CSCSegmentCollection> cscSegments; 
  event.getByLabel(cscSegmentLabel, cscSegments);

  // Get the DT RecHits collection :
  Handle<DTRecHitCollection> dt1DRecHits; 
  event.getByLabel(dtrecHitLabel, dt1DRecHits);

  // Get the DT Segments collection :
  Handle<DTRecSegment4DCollection> dt4DSegments;
  event.getByLabel(dtSegmentLabel, dt4DSegments);

  // Get the SimHit collection :
  Handle<PSimHitContainer> csimHits;
  event.getByLabel(simHitLabel,"MuonCSCHits", csimHits);  
  Handle<PSimHitContainer> dsimHits;
  event.getByLabel(simHitLabel,"MuonDTHits", dsimHits);  
 
  // Get the SimTrack
  Handle<SimTrackContainer> simTracks;
  event.getByLabel(simTrackLabel, simTracks);

  H2DRecHit1 *histo1 = 0;   
  H2DRecHit2 *histo2 = 0;   
  H2DRecHit3 *histo3 = 0;   
  H2DRecHit4 *histo4 = 0;   
  H2DRecHit5 *histo5 = 0;   
  H2DRecHit6 *histo6 = 0;   
  H2DRecHit7 *histo7 = 0;   
  H2DRecHit10 *histo10 = 0;   
  //  TNtuple1 *tt = 0;
 
  // 0. Run the class SegSelector
  //SegSelector recsegSelector();
  std::vector<SimSegment> sCSC_v      = recsegSelector->Sim_CSCSegments(1,csimHits,cscGeom);
  std::vector<SimSegment> sDT_v       = recsegSelector->Sim_DTSegments(1,dsimHits,dtGeom);
  std::vector<CSCSegment> cscseg_V    = recsegSelector->Select_CSCSeg(cscSegments,cscGeom, sCSC_v);
  std::vector<DTRecSegment4D> dtseg_V = recsegSelector->Select_DTSeg(dt4DSegments,dtGeom, sDT_v);

  // 1. Get the track information  
  SimInfo(simTracks,dsimHits,csimHits,dtGeom,cscGeom);

  // 2. Check # of segments in each chambers
  CSCsegment_stat(cscSegments);
  DTsegment_stat(dt4DSegments);
  if ( (cscseg_stat[5] < 2)&&(dtseg_stat[5] < 2) ) {
     /// check # of rechit for those 1 or 0 segment event
     CSCRecHit_stat(csc2DRecHits,cscGeom);
     DTRecHit_stat(dt1DRecHits,dtGeom);
  }
  // 3. Event statistic
  // number of segments or rechits w.r.t. eta in each event 
  int allmu_stat = cscseg_stat[5] + dtseg_stat[5];
  int allrh_stat =   cscrh_sum[5] + dtrh_sum[5];
  double allmu_eta = -9.0;
  if      ((cscseg_stat[5]==0)&&(dtseg_stat[5]!=0)) { allmu_eta = eta_d; }
  else if ((cscseg_stat[5]!=0)&&(dtseg_stat[5]==0)) { allmu_eta = eta_c; }
  else if ((cscseg_stat[5]!=0)&&(dtseg_stat[5]!=0)) { allmu_eta = (eta_c + eta_d)/2 ; }
  else {
       // look up the rechits vs eta in events with no segments
       if      ((eta_d == -9.0)&&(eta_c != -9.0)) {allmu_eta = eta_c ; }
       else if ((eta_d != -9.0)&&(eta_c == -9.0)) {allmu_eta = eta_d ; }
       else    { allmu_eta = (eta_c + eta_d)/2 ; }

       histo1 = h_all;
       histo1->Fill1a( allmu_eta, allrh_stat, cscrh_sum[0], dtrh_sum[0]);
  }
  histo1 = h_all;
  histo1->Fill1(cscseg_stat[5],dtseg_stat[5],allmu_stat,eta_c,eta_d,allmu_eta,eta_trk);

  // look up the # of segments in each station
  if (cscseg_stat[0] != 0){
     histo2 = h_csc;
     histo2->Fill3(pt1[0],pa[0],cscseg_stat[0]);
     for (int k=1; k<5; k++) {
         histo2->Fill4(k,cscseg_stat[k],cscseg_stat1[k]);
     }
  }
  if (dtseg_stat[0] != 0) {
     histo3 = h_dt;
     histo3->Fill3a(pt1[0],pa[0],dtseg_stat[0]); 
     for (int k=1; k<5; k++) {
         histo3->Fill4a(k,dtseg_stat[k],dtseg_stat1[k]);
     }
  }

  // flag the overlap case
  bool overlap = false;
  if ( (dtseg_stat[0] != 0)&&(cscseg_stat[0] != 0) ){
     overlap = true;
  }


  // 4. the simulated segments statistic
  int simcscseg[6] = {0};
  double ns1 =0.0;
  double eta_sim1 =0;
  for (std::vector<SimSegment>::const_iterator it = sCSC_v.begin(); it != sCSC_v.end(); it++) {
      int st = ((*it).csc_DetId).station();
      eta_sim1 += fabs(((*it).sGlobalOrg).eta());
      simcscseg[st]++;
      ns1++;
  }
  simcscseg[0]=simcscseg[1]+simcscseg[2]+simcscseg[3]+simcscseg[4];
  for (int i=1; i<5; i++) {
      if (simcscseg[i]!=0) { simcscseg[5]++; }
  }
  eta_sim1 = eta_sim1/ns1 ;

  int simdtseg[6] = {0};
  double ns2 =0.0;
  double eta_sim2 =0;
  for (std::vector<SimSegment>::const_iterator it = sDT_v.begin(); it != sDT_v.end(); it++) {
      int st = ((*it).dt_DetId).station();
      eta_sim2 += fabs(((*it).sGlobalOrg).eta());
      simdtseg[st]++;
      ns2++;
  }
  simdtseg[0]=simdtseg[1]+simdtseg[2]+simdtseg[3]+simdtseg[4];
  for (int i=1; i<5; i++) {
      if (simdtseg[i]!=0) { simdtseg[5]++; }
  }
  eta_sim2 = eta_sim2/ns2 ;

  int allmu_stat1 = simcscseg[5]+ simdtseg[5];
  double allmu_eta1 = -9.0;
  if      ((simcscseg[0]==0)&&( simdtseg[0]!=0)) { allmu_eta1 = eta_sim2; }
  else if (( simdtseg[0]==0)&&(simcscseg[0]!=0)) { allmu_eta1 = eta_sim1; }
  else { allmu_eta1 = (eta_sim1 + eta_sim2)/2 ; }
    
  histo1 = h_all;
  histo1->Fill1b(allmu_stat1,allmu_eta1);


  //5.  look at different Bxdl btw. stations
  FromCSCSeg(cscseg_V,cscGeom,sCSC_v);
  FromDTSeg(dtseg_V,dtGeom,sDT_v);
  FromOverlap();

  FromCSCSingleSeg(cscseg_V,cscGeom,sCSC_v);
  FromDTSingleSeg(dtseg_V,dtGeom,sDT_v);

  /// chi2 distribution
  for (int i=0; i<5; i++) {
      if (chi2_dof1[i]< 0.0) continue;     
      histo2 = h_csc;
      histo2->Fill3b( chi2_dof1[i] );
  }
  for (int i=1; i<5; i++) {
      if (chi2_dof3[i]< 0.0) continue;
      histo3 = h_dt;
      histo3->Fill3c( chi2_dof3[i] );
  }


  if (dPhiV1[1][0][1]!=99.0 && chi2_dof1[0] < 2000.0 && chi2_dof1[1] < 2000.0){
     histo2 = h_csc;
     histo2->Fill5_0( dPhiV1[0][0][1],dPhiV1[1][0][1],pt1[1]*dPhiV1[0][0][1],pt1[1]*dPhiV1[1][0][1],
		     fabs(EtaP1[0][1]),fabs(EtaP1[1][1]) );
  }
  if (dPhiV1[1][0][2]!=99.0 && chi2_dof1[0] < 2000.0 && chi2_dof1[2] < 2000.0){
     histo2 = h_csc;
     histo2->Fill5_1( dPhiV1[0][0][2],dPhiV1[1][0][2],pt1[1]*dPhiV1[0][0][2],pt1[1]*dPhiV1[1][0][2],
		     fabs(EtaP1[0][2]),fabs(EtaP1[1][2]) );
  }
  if (dPhiV1[1][0][3]!=99.0 && chi2_dof1[0] < 2000.0 && chi2_dof1[3] < 2000.0){
     histo2 = h_csc;
     histo2->Fill5_2( dPhiV1[0][0][3],dPhiV1[1][0][3],pt1[1]*dPhiV1[0][0][3],pt1[1]*dPhiV1[1][0][3],
		     fabs(EtaP1[0][3]),fabs(EtaP1[1][3]) );
  }
  if (dPhiV1[1][0][4]!=99.0 && chi2_dof1[0] < 2000.0 && chi2_dof1[4] < 2000.0){
     histo2 = h_csc;
     histo2->Fill5_3( dPhiV1[0][0][4],dPhiV1[1][0][4],pt1[1]*dPhiV1[0][0][4],pt1[1]*dPhiV1[1][0][4],
		     fabs(EtaP1[0][4]),fabs(EtaP1[1][4]) );
  }
  if (dPhiV1[1][1][2]!=99.0 && chi2_dof1[1] < 2000.0 && chi2_dof1[2] < 2000.0){
     histo2 = h_csc;
     histo2->Fill5_1( dPhiV1[0][1][2],dPhiV1[1][1][2],pt1[1]*dPhiV1[0][1][2],pt1[1]*dPhiV1[1][1][2],
		     fabs(EtaP1[0][2]),fabs(EtaP1[1][2]) );
  }
  if (dPhiV1[1][1][3]!=99.0 && chi2_dof1[1] < 2000.0 && chi2_dof1[3] < 2000.0){
     histo2 = h_csc;
     histo2->Fill5_2( dPhiV1[0][1][3],dPhiV1[1][1][3],pt1[1]*dPhiV1[0][1][3],pt1[1]*dPhiV1[1][1][3],
		     fabs(EtaP1[0][3]),fabs(EtaP1[1][3]) );
  }
  if (dPhiV1[1][1][4]!=99.0 && chi2_dof1[1] < 2000.0 && chi2_dof1[4] < 2000.0){
     histo2 = h_csc;
     histo2->Fill5_3( dPhiV1[0][1][4],dPhiV1[1][1][4],pt1[1]*dPhiV1[0][1][4],pt1[1]*dPhiV1[1][1][4],
		     fabs(EtaP1[0][4]),fabs(EtaP1[1][4]) );
  }
  if (dPhiV1[1][2][3]!=99.0 && chi2_dof1[2] < 2000.0 && chi2_dof1[3] < 2000.0){
     histo2 = h_csc;
     histo2->Fill5_4( dPhiV1[0][2][3],dPhiV1[1][2][3],pt1[2]*dPhiV1[0][2][3],pt1[2]*dPhiV1[1][2][3],
		     fabs(EtaP1[0][3]),fabs(EtaP1[1][3]) );
  }
  if (dPhiV1[1][2][4]!=99.0 && chi2_dof1[2] < 2000.0 && chi2_dof1[4] < 2000.0){
     histo2 = h_csc;
     histo2->Fill5_5( dPhiV1[0][2][4],dPhiV1[1][2][4],pt1[2]*dPhiV1[0][2][4],pt1[2]*dPhiV1[1][2][4],
		     fabs(EtaP1[0][4]),fabs(EtaP1[1][4]) );
  }
  if (dPhiV1[1][3][4]!=99.0 && chi2_dof1[3] < 2000.0 && chi2_dof1[4] < 2000.0){
     histo2 = h_csc;
     histo2->Fill5_6( dPhiV1[0][3][4],dPhiV1[1][3][4],pt1[3]*dPhiV1[0][3][4],pt1[3]*dPhiV1[1][3][4],
		     fabs(EtaP1[0][4]),fabs(EtaP1[1][4]) );
  }


  /// For DT
  if (dPhiV3[1][1][2]!=99.0 && chi2_dof3[1] < 2000.0 && chi2_dof3[2] < 2000.0){
     histo3 = h_dt;
     histo3->Fill6_1( dPhiV3[0][1][2],dPhiV3[1][1][2],pt1[1]*dPhiV3[0][1][2],pt1[1]*dPhiV3[1][1][2],
		     fabs(EtaP3[0][2]),fabs(EtaP3[1][2]) );
  }
  if (dPhiV3[1][1][3]!=99.0 && chi2_dof3[1] < 2000.0 && chi2_dof3[3] < 2000.0){
     histo3 = h_dt;
     histo3->Fill6_2( dPhiV3[0][1][3],dPhiV3[1][1][3],pt1[1]*dPhiV3[0][1][3],pt1[1]*dPhiV3[1][1][3],
		     fabs(EtaP3[0][3]),fabs(EtaP3[1][3]) );
  }
  if (dPhiV3[1][1][4]!=99.0 && chi2_dof3[1] < 2000.0 && chi2_dof3[4] < 2000.0){
     histo3 = h_dt;
     histo3->Fill6_3( dPhiV3[0][1][4],dPhiV3[1][1][4],pt1[1]*dPhiV3[0][1][4],pt1[1]*dPhiV3[1][1][4],
		     fabs(EtaP3[0][1]),fabs(EtaP3[1][1]) );
  }
  if (dPhiV3[1][2][3]!=99.0 && chi2_dof3[2] < 2000.0 && chi2_dof3[3] < 2000.0){
     histo3 = h_dt;
     histo3->Fill6_4( dPhiV3[0][2][3],dPhiV3[1][2][3],pt1[2]*dPhiV3[0][2][3],pt1[2]*dPhiV3[1][2][3],
		     fabs(EtaP3[0][3]),fabs(EtaP3[1][3]) );
  }
  if (dPhiV3[1][2][4]!=99.0 && chi2_dof3[2] < 2000.0 && chi2_dof3[4] < 2000.0){
     histo3 = h_dt;
     histo3->Fill6_5( dPhiV3[0][2][4],dPhiV3[1][2][4],pt1[2]*dPhiV3[0][2][4],pt1[2]*dPhiV3[1][2][4],
		     fabs(EtaP3[0][2]),fabs(EtaP3[1][2]) );
  }
  if (dPhiV3[1][3][4]!=99.0 && chi2_dof3[3] < 2000.0 && chi2_dof3[4] < 2000.0){
     histo3 = h_dt;
     histo3->Fill6_6( dPhiV3[0][3][4],dPhiV3[1][3][4],pt1[3]*dPhiV3[0][3][4],pt1[3]*dPhiV3[1][3][4],
		     fabs(EtaP3[0][3]),fabs(EtaP3[1][3]) );
  }


  //  Look at different Bxdl btw. stations & rings
  /// All possible segment pairs in CSC 
  ///                 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   
  int csc1[2][15]= {{ 1, 1,12,12,13, 1, 1,12,13, 1,21,21,22,21,31},
                    {12,21,21,22,22,31,32,32,32,41,31,32,32,41,41}}; 
  for (int l=0; l<15; l++) {
      int s1 = csc1[0][l]/10; // 0 0 1 1 1 0 0 1 1 0 2 2 2 2 3
      int r1 = csc1[0][l]%10; // 1 1 2 2 3 1 1 2 3 1 1 1 2 1 1
      int s2 = csc1[1][l]/10; 
      int r2 = csc1[1][l]%10;
      if ( MEPath[1][s1][r1] && MEPath[1][s2][r2] && MEPath[0][s1][r1] && MEPath[0][s2][r2] ) {
         double ME_Resol = dPhiP1[0][s1][s2]-dPhiP1[1][s1][s2];
         histo4 = hME1[l];
         histo4->Fill8(  (pt1[s1]*dPhiP1[0][s1][s2]) ,dPhiP1[0][s1][s2] ,dEtaP1[0][s1][s2], fabs(EtaP1[0][s2]), pt1[s1] );
         histo4->Fill8a( (pt1[s1]*dPhiP1[1][s1][s2]) ,dPhiP1[1][s1][s2] ,dEtaP1[1][s1][s2], fabs(EtaP1[1][s2]), pt1[s1], ME_Resol );
      }
  }
  
  //  Look at different Bxdl btw. stations & rings
  /// All possible segment pair in DT 
  ///               0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 
  int dt1[2][26]={{10,10,10,10,10,11,11,11,11,11,11,12,12,20,20,20,21,21,21,21,22,22,30,31,31,32}
                 ,{20,30,31,40,41,21,22,31,32,41,42,22,32,30,40,41,31,32,41,42,32,42,40,41,42,42}};
  for (int l=0; l<26; l++) {
      int s1 = dt1[0][l]/10;
      int w1 = dt1[0][l]%10;
      int s2 = dt1[1][l]/10;
      int w2 = dt1[1][l]%10;
       
      if ( MBPath[1][s1][w1] && MBPath[1][s2][w2] && MBPath[0][s1][w1] && MBPath[0][s2][w2] ) {
         double MB_Resol = dPhiP3[0][s1][s2]-dPhiP3[1][s1][s2];
         if (s2!=4 ){
            histo5 = hMB1[l];
            histo5->Fill9( (pt1[s1]*dPhiP3[0][s1][s2]), dPhiP3[0][s1][s2], dEtaP3[0][s1][s2], fabs(EtaP3[0][s2]), pt1[s1] );
            histo5->Fill9a((pt1[s1]*dPhiP3[1][s1][s2]), dPhiP3[1][s1][s2], dEtaP3[1][s1][s2], fabs(EtaP3[1][s2]), pt1[s1], MB_Resol );
         }
         if (s2==4 ){
            histo5 = hMB1[l];
            histo5->Fill9( (pt1[s1]*dPhiP3[0][s1][s2]), dPhiP3[0][s1][s2], dEtaP3[0][s1][s2], fabs(EtaP3[0][s1]), pt1[s1] );
            histo5->Fill9a((pt1[s1]*dPhiP3[1][s1][s2]), dPhiP3[1][s1][s2], dEtaP3[1][s1][s2], fabs(EtaP3[1][s1]), pt1[s1], MB_Resol );
         }
      }
  }

  //  Look at different Bxdl in overlap region
  /// All possible segment pairs in overlap region 
  ///              0  1  2  3  4  5   
  int olp[2][6]={{12,12,12,22,22,32} 
                ,{13,22,32,13,22,13}};
  for (int l=0; l<6; l++) {
      int s1 = olp[0][l]/10;
      int w1 = olp[0][l]%10;
      int s2 = olp[1][l]/10;
      int w2 = olp[1][l]%10; 
      if ( MBPath[1][s1][w1] && MEPath[1][s2][w2] ) {
         double OL_Resol = dPhiP2[s1][s2]-dPhiP2[s1][s2];
         histo10 = hOL1[l];
         histo10->Fill12( (pt1[s1]*dPhiP2[0][s1][s2]), dPhiP2[0][s1][s2], dEtaP2[0][s1][s2], fabs(EtaP3[0][s2]), pt1[s1] );
         histo10->Fill12a((pt1[s1]*dPhiP2[1][s1][s2]), dPhiP2[1][s1][s2], dEtaP2[1][s1][s2], fabs(EtaP3[1][s2]), pt1[s1], OL_Resol );
      }
  }
  
  // Fill  the 1 segment case for CSC and DT
  int csc2[8] ={ 1,12,13,21,22,31,32,41};
  for (int l=0; l<8; l++) {
      int s1 = csc2[l]/10;
      int r1 = csc2[l]%10;
      if ( MEPath[1][s1][r1] && MEPath[0][s1][r1]) {
         double dME_phi = ME_phi[1][s1][r1] - ME_phi[0][s1][r1];
         double dME_eta = ME_eta[1][s1][r1] - ME_eta[0][s1][r1];
         histo6 = hME2[l];
         histo6->Fill8b( (pt1[0]*ME_phi[0][s1][r1]) ,ME_phi[0][s1][r1], ME_eta[0][s1][r1], pt1[0]);
         histo6->Fill8c( (pt1[0]*ME_phi[1][s1][r1]) ,ME_phi[1][s1][r1], dME_phi , dME_eta, ME_eta[1][s1][r1] ,pt1[0]);
      
      }
  }

  int dt2[12] ={10,11,12,20,21,22,30,31,32,40,41,42};
  for (int l=0; l<12; l++) {
      int s1 = dt2[l]/10;
      int w1 = dt2[l]%10;
      if ( MBPath[1][s1][w1] && MBPath[0][s1][w1] ) {
         double dMB_phi = MB_phi[1][s1][w1] - MB_phi[0][s1][w1];
         double dMB_eta = MB_eta[1][s1][w1] - MB_eta[0][s1][w1];
         histo7 = hMB2[l];
         histo7->Fill9b( (pt1[0]*MB_phi[0][s1][w1]) , MB_phi[0][s1][w1], MB_eta[0][s1][w1], pt1[0]);
         histo7->Fill9c( (pt1[0]*MB_phi[1][s1][w1]) , MB_phi[1][s1][w1], dMB_phi , dMB_eta, MB_eta[1][s1][w1], pt1[0]);

      }
  }

  /// For reco-segment treea
  /*
  tt = tr_muon;
  if ( MEPath[1][1] && MEPath[1][2] && MEPath[1][3] ) {
      tt->Fill_b1(fabs(EtaP1[1][1]),fabs(EtaP1[1][2]),fabs(EtaP1[1][3]),fabs(EtaP1[1][4]), 
                  EtaP1[1][1], EtaP1[1][2], EtaP1[1][3], EtaP1[1][4], pt1[0]);
      tt->Fill_l1(pa[0]);
  }
  if ( MBPath[1][1] && MBPath[1][2] && MBPath[1][3] ) {
      tt->Fill_b2(fabs(EtaP3[1][1]),fabs(EtaP3[1][2]),fabs(EtaP3[1][3]),fabs(EtaP3[1][4]), 
                  EtaP3[1][1], EtaP3[1][2], EtaP3[1][3], EtaP3[1][4], pt1[0]);
      tt->Fill_l1(pa[0]);
  }
  tt->FillTree();
  */ 
  
}

// ********************************************
// ***********  Utility functions  ************
// ********************************************

// number of csc segments in one chamber for each station
// cscseg_stat[0] = total segments in all stations
// cscseg_stat[i] = the number of segments in station i
// cscseg_stat[5] = the number of stations which have segments
// cscseg_stat1 is the statistic for segments w/ more than 4 rechits
void SeedParametrization::CSCsegment_stat( Handle<CSCSegmentCollection> cscSeg ) {

     for (int i=0; i<6; i++) {
         cscseg_stat[i]=0;
         cscseg_stat1[i]=0;
     }
     for(CSCSegmentCollection::const_iterator seg_It = cscSeg->begin(); seg_It != cscSeg->end(); seg_It++)
     { 
        CSCDetId DetId = (CSCDetId)(*seg_It).cscDetId();
        cscseg_stat[DetId.station()] += 1;
        if ((*seg_It).nRecHits() > 4 ) {
           cscseg_stat1[DetId.station()] += 1;
        }
     }
     cscseg_stat[0] = cscseg_stat[1]+cscseg_stat[2]+cscseg_stat[3]+cscseg_stat[4];
     cscseg_stat1[0] = cscseg_stat1[1]+cscseg_stat1[2]+cscseg_stat1[3]+cscseg_stat1[4];
     for (int i =1; i<5; i++){
         if(cscseg_stat[i]!=0) {
            cscseg_stat[5]++ ;
            cscseg_stat1[5]++ ;
         }
     }
     
}
// the same as CSCsegment_stat
void SeedParametrization::DTsegment_stat( Handle<DTRecSegment4DCollection> dtSeg ) {

     for (int i=0; i<6; i++) {
         dtseg_stat[i]=0;
         dtseg_stat1[i]=0;
     }
     for(DTRecSegment4DCollection::const_iterator seg_It = dtSeg->begin(); seg_It != dtSeg->end(); seg_It++)
     { 
        DTChamberId DetId = (*seg_It).chamberId();
        dtseg_stat[DetId.station()] += 1;
        int n_phiHits = ((*seg_It).phiSegment())->specificRecHits().size();
        if ( (*seg_It).hasZed() && (n_phiHits > 4) ) {
           dtseg_stat1[DetId.station()] += 1;
        }
     }
     dtseg_stat[0] = dtseg_stat[1]+dtseg_stat[2]+dtseg_stat[3]+dtseg_stat[4];
     dtseg_stat1[0] = dtseg_stat1[1]+dtseg_stat1[2]+dtseg_stat1[3]+dtseg_stat1[4];
     for (int i =1; i<5; i++){
         if(dtseg_stat[i]!=0) {
            if((i==4)&&(dtseg_stat[5]==0)) {
              dtseg_stat[5]=0;
            } else {
              dtseg_stat[5]++ ;
              dtseg_stat1[5]++ ;
            } 
         }
     }

}
// number rechit in each station, basically only for those station cannot form a segment
void SeedParametrization::CSCRecHit_stat(Handle<CSCRecHit2DCollection> cscrechit, ESHandle<CSCGeometry> cscGeom){
     for (int i=0; i <6; i++) {
         cscrh_sum[i]=0;
     }
     for(CSCRecHit2DCollection::const_iterator r_it = cscrechit->begin(); r_it != cscrechit->end(); r_it++)
     { 
        CSCDetId det_id = (CSCDetId)(*r_it).cscDetId();
        //const CSCLayer* csclayer = cscGeom->layer( det_id );
        //const CSCChamber* cscchamber = cscGeom->chamber( det_id );
        //LocalPoint lrh = (*r_it).localPosition();
        //GlobalPoint grh = csclayer->toGlobal(lrh);
        cscrh_sum[det_id.station()]++;
     }
     cscrh_sum[0] = cscrh_sum[1]+cscrh_sum[2]+cscrh_sum[3]+cscrh_sum[4];
     for (int i =1; i<5; i++){
         if(cscrh_sum[i]!=0) {
            cscrh_sum[5]++ ;
         }
     }
}

void SeedParametrization::DTRecHit_stat(Handle<DTRecHitCollection> dtrechit, ESHandle<DTGeometry> dtGeom){

     //double phi[4]={999.0};
     for (int i=0; i <6; i++) {
         dtrh_sum[i]=0;
     }

     double eta=-9.0;
     double nn=0.0;
     for (DTRecHitCollection::const_iterator r_it = dtrechit->begin(); r_it != dtrechit->end(); r_it++){
         DTWireId det_id = (*r_it).wireId();
         const DTChamber* dtchamber = dtGeom->chamber( det_id );
         LocalPoint lrh = (*r_it).localPosition();
         GlobalPoint grh = dtchamber->toGlobal( lrh );
         dtrh_sum[det_id.station()]++;
         eta += grh.eta();
         nn += 1.0;
     }
     eta = eta/nn ;

     dtrh_sum[0] = dtrh_sum[1]+dtrh_sum[2]+dtrh_sum[3]+dtrh_sum[4];
     for (int i =1; i<5; i++){
         if (dtrh_sum[i]!=0) {
            dtrh_sum[5]++ ;
         }
     }
}

// find the simHits which is corresponding to the segment
bool SeedParametrization::SameChamber(CSCDetId SimDetId, CSCDetId SegDetId){
  
     if ( SimDetId.endcap()== SegDetId.endcap() && SimDetId.station()== SegDetId.station() &&
          SimDetId.ring()  == SegDetId.ring()   && SimDetId.chamber()== SegDetId.chamber() ){
          return true;
     }
     else {
          return false;
     }

}

void SeedParametrization::SimInfo(Handle<edm::SimTrackContainer> simTracks,
                            Handle<edm::PSimHitContainer> dsimHits, Handle<edm::PSimHitContainer> csimHits,
                            ESHandle<DTGeometry> dtGeom, ESHandle<CSCGeometry> cscGeom){

  // pt1 -> pt in different station pt1[0] is the track pt
  for (int i=0; i<5; i++) {
      pt1[i] = 0.0;
  }
  for (int i=0; i<5; i++) {
      pa[i] = 0.0;
  }

  // eta_c -> ave.eta from all csc stations ; cta_d -> ave.eta from dt stations
  eta_c = -9.0;
  eta_d = -9.0;
  eta_trk = -9.0;

  for (SimTrackContainer::const_iterator simTk_It = simTracks->begin(); simTk_It != simTracks->end(); simTk_It++)
  {
      if (abs((*simTk_It).type())!=13) continue;

      if ((*simTk_It).type()==13) {
         theQ = -1.0;
      }else {
         theQ = 1.0;
      }

      float px = ((*simTk_It).momentum()).x();
      float py = ((*simTk_It).momentum()).y();
      float pz = ((*simTk_It).momentum()).z();
      pt1[0] = sqrt((px*px) + (py*py));
      pa[0] = sqrt((px*px) + (py*py) + (pz*pz));
      double theta = acos( pz/pa[0] );
      eta_trk = fabs((-1.0)*log( tan(theta/2.0) ));

      double eta_c1=0.0;
      double enu1 = 0.0;
      for (PSimHitContainer::const_iterator cs_It = csimHits->begin(); cs_It != csimHits->end(); cs_It++)
      {
          CSCDetId C_Id = CSCDetId((*cs_It).detUnitId());
          const CSCChamber* cscchamber = cscGeom->chamber( C_Id );
          GlobalVector m1 = cscchamber->toGlobal((*cs_It).momentumAtEntry() );
          Local3DPoint lp = (*cs_It).localPosition();
          GlobalPoint gp = cscchamber->toGlobal(lp );
          if ( ( abs((*cs_It).particleType())==13 ) && ( (*cs_It).trackId()==1 )) {
             pt1[C_Id.station()] = sqrt( (m1.x()*m1.x()) + (m1.y()*m1.y()) );
             pa[C_Id.station()] = sqrt( (m1.x()*m1.x()) + (m1.y()*m1.y()) + (m1.z()*m1.z()));
             eta_c1 += fabs(gp.eta());
             enu1   += 1.0;
          }
      }
      if (enu1 !=0.0 ) {
          eta_c = eta_c1/enu1 ;
      }else {
          eta_c = -9.0 ;
      }

      double eta_d1=0.0;
      double enu2 = 0.0;
      for (PSimHitContainer::const_iterator ds_It = dsimHits->begin(); ds_It != dsimHits->end(); ds_It++)
      {
          Local3DPoint lp = (*ds_It).localPosition();

          DTLayerId D_Id = DTLayerId( (*ds_It).detUnitId() );
          const DTLayer* dtlayer = dtGeom->layer(D_Id);
          GlobalVector m2 = dtlayer->toGlobal((*ds_It).momentumAtEntry() );
          GlobalPoint gp = dtlayer->toGlobal(lp );

          if ( ( abs((*ds_It).particleType())==13 ) && ( (*ds_It).trackId()==1 )) {
             pt1[D_Id.station()] = sqrt( (m2.x()*m2.x()) + (m2.y()*m2.y()) );
             pa[D_Id.station()] = sqrt( (m2.x()*m2.x()) + (m2.y()*m2.y()) + (m2.z()*m2.z()));
             eta_d1 += fabs(gp.eta());
             enu2  += 1.0;
          }
      }
      if (enu2 !=0.0 ) {
         eta_d = eta_d1/enu2 ;
      }else {
         eta_d = -9.0 ;
      }
  }
}

// Fill the phi and eta information for CSC
void SeedParametrization::FromCSCSeg( std::vector<CSCSegment> cscSeg, ESHandle<CSCGeometry> cscGeom,
                                           std::vector<SimSegment> seg) {
     
     /// get the global dphi from recHits for CSC
     for (int l=0; l<10; l++) {
        int i = l/2; // 0 0 1 1 2 2 3 3 4 4 
        int k = l%2; // 0 1 0 1 0 1 0 1 0 1
        PhiV1[k][i]=99.;
	EtaV1[k][i]=99.;
	PhiP1[k][i]=99.;
	EtaP1[k][i]=99.;
        chi2_dof1[i]= -1.0;
        for (int j=0; j<5; j++) {
            dPhiV1[k][i][j]=99.;
            dEtaV1[k][i][j]=99.;
            dPhiP1[k][i][j]=99.;
            dEtaP1[k][i][j]=99.;
        }
     }
     bool  layer[5]= { false };
     // Fill the phi and eta of segment direction in different station
     for (std::vector<CSCSegment>::const_iterator it = cscSeg.begin(); it != cscSeg.end(); it++) {

         CSCDetId DetId = (CSCDetId)(*it).cscDetId();
         const CSCChamber* cscchamber = cscGeom->chamber( DetId );
         GlobalPoint  gp = cscchamber->toGlobal( (*it).localPosition()  );
	 GlobalVector gv = cscchamber->toGlobal( (*it).localDirection() );
	 int st = DetId.station();
	 int rg = DetId.ring();
  
         if (st==1 && (rg ==1 || rg==4) ) {
            PhiV1[1][0] = gv.phi();
            EtaV1[1][0] = gv.eta();
            PhiP1[1][0] = gp.phi();
            EtaP1[1][0] = gp.eta();
            layer[0] =true;
            chi2_dof1[st] = (*it).chi2() / static_cast<double>((*it).degreesOfFreedom()); 
         }
         else {
            PhiV1[1][st] = gv.phi();
            EtaV1[1][st] = gv.eta();
            PhiP1[1][st] = gp.phi();
            EtaP1[1][st] = gp.eta();
            layer[st]=true;
            chi2_dof1[st] = (*it).chi2() / static_cast<double>((*it).degreesOfFreedom()); 
         }
     }
     for (std::vector<SimSegment>::const_iterator it = seg.begin(); it != seg.end(); it++) {
         if ( (*it).chamber_type != 1) continue; 
         GlobalPoint  gp = (*it).sGlobalOrg;
         GlobalVector gv = (*it).sGlobalVec;
         int st = (*it).csc_DetId.station();
         int rg = (*it).csc_DetId.ring();
         if (st==1 && (rg ==1 || rg==4) ) {
            PhiV1[0][0] = gv.phi();
            EtaV1[0][0] = gv.eta();
            PhiP1[0][0] = gp.phi();
            EtaP1[0][0] = gp.eta();
         }
         else {
            PhiV1[0][st] = gv.phi();
            EtaV1[0][st] = gv.eta();
            PhiP1[0][st] = gp.phi();
            EtaP1[0][st] = gp.eta();
         }
     }
     
     // Get the dPhi and dEta for different combination
     for (int m=0; m<2; m++){
        for (int l=0; l<16; l++) {
           int s1 = (l%4);    // 0  0,1  0,1,2  0,1,2,3
           int s2 = (l/4)+1;  // 1  2,2  3,3,3  4,4,4,4
           if ( layer[s1] && layer[s2] && (s1<s2) ) {
  	      dPhiV1[m][s1][s2] = PhiV1[m][s1]-PhiV1[m][s2];  
	      dEtaV1[m][s1][s2] = EtaV1[m][s1]-EtaV1[m][s2];
	      dPhiP1[m][s1][s2] = PhiP1[m][s1]-PhiP1[m][s2];  
	      dEtaP1[m][s1][s2] = EtaP1[m][s1]-EtaP1[m][s2];
           }
        }
     }
}

// Fill the phi and eta information for DT
void SeedParametrization::FromDTSeg( std::vector<DTRecSegment4D> dtSeg, ESHandle<DTGeometry> dtGeom,
                                          std::vector<SimSegment> seg) {
     
     /// get the global dphi from recHits for DT
     for (int l=0; l<10; l++) {
        int i = l/2; // 0 0 1 1 2 2 3 3 4 4
        int k = l%2; // 0 1 0 1 0 1 0 1 0 1
        PhiV3[k][i]=99.;
	EtaV3[k][i]=99.;
	PhiP3[k][i]=99.;
	EtaP3[k][i]=99.;
        chi2_dof3[i]= -1.0;
        for (int j=0; j<5; j++) {
            dPhiV3[k][i][j]=99.;
            dEtaV3[k][i][j]=99.;
            dPhiP3[k][i][j]=99.;
            dEtaP3[k][i][j]=99.;
        }
     }
     bool  layer[5]= { false };
     // Fill the phi and eta of segment direction in different station
     for (std::vector<DTRecSegment4D>::const_iterator it = dtSeg.begin(); it != dtSeg.end(); it++) {

         DTChamberId DetId = (*it).chamberId();
         const DTChamber* dtchamber = dtGeom->chamber( DetId );
         GlobalPoint  gp = dtchamber->toGlobal( (*it).localPosition() );
         GlobalVector gv = dtchamber->toGlobal( (*it).localDirection() );
	 int st = DetId.station();
	 PhiV3[1][st] = gv.phi();
	 EtaV3[1][st] = gv.eta();
	 PhiP3[1][st] = gp.phi();
	 EtaP3[1][st] = gp.eta();
	 layer[st] = true;
         chi2_dof3[st] = (*it).chi2() / static_cast<double>((*it).degreesOfFreedom()); 
     }
     for (std::vector<SimSegment>::const_iterator it = seg.begin(); it != seg.end(); it++) {
         if ( (*it).chamber_type != 2) continue; 
         GlobalPoint  gp = (*it).sGlobalOrg;
         GlobalVector gv = (*it).sGlobalVec;
         int st = (*it).dt_DetId.station();
	 PhiV3[0][st] = gv.phi();
	 EtaV3[0][st] = gv.eta();
	 PhiP3[0][st] = gp.phi();
	 EtaP3[0][st] = gp.eta();
     }
     // Get the dPhi and dEta for different combination
     for (int m=0; m<2; m++) {
         for (int l=0; l<9; l++) {
            int s1 = (l%3)+1;  // 1  1,2  1,2,3  
	    int s2 = (l/3)+2;  // 2  3,3  4,4,4  
	    if ( layer[s1] && layer[s2] && (s1<s2) ) {
	       dPhiV3[m][s1][s2] = PhiV3[m][s1]-PhiV3[m][s2];  
	       dEtaV3[m][s1][s2] = EtaV3[m][s1]-EtaV3[m][s2];
	       dPhiP3[m][s1][s2] = PhiP3[m][s1]-PhiP3[m][s2];  
	       dEtaP3[m][s1][s2] = EtaP3[m][s1]-EtaP3[m][s2];
               //cout <<" ---------------------- DT  ------------------------------------------------"<<endl;
               //cout <<"dPhi from P = "<< PhiP3[m][s1] <<" - "<<PhiP3[m][s2]<<" = "<<dPhiP3[m][s1][s2]<<endl;
               //cout <<"dPhi from V = "<< PhiV3[m][s1] <<" - "<<PhiV3[m][s2]<<" = "<<dPhiV3[m][s1][s2]<<endl;
            }
         }
     }
}

void SeedParametrization::FromOverlap() {
 
     for (int l=0; l<10; l++) {
        int i = l/2;
        int k = l%2;
        for (int j=0; j<5; j++) {
            dPhiV2[k][i][j]=99.;
            dEtaV2[k][i][j]=99.;
            dPhiP2[k][i][j]=99.;
            dEtaP2[k][i][j]=99.;
        }
     }
     for (int m=0; m<2; m++) {
         for (int l=0; l<9; l++) {
            int s1 = (l%3)+1;  // 1,2,3  1,2,3  1,2,3  
	    int s2 = (l/3)+1;  // 1,1,1  2,2,2  3,3,3 

	    if ( (PhiV3[m][s1]==99.0) || (PhiV1[m][s2]==99.0 )) continue;
	    dPhiV2[m][s1][s2] = PhiV3[m][s1]-PhiV1[m][s2];  
	    dEtaV2[m][s1][s2] = EtaV3[m][s1]-EtaV1[m][s2];
	    dPhiP2[m][s1][s2] = PhiP3[m][s1]-PhiP1[m][s2];  
	    dEtaP2[m][s1][s2] = EtaP3[m][s1]-EtaP1[m][s2];
         }
     }

}
     
void SeedParametrization::FromCSCSingleSeg( std::vector<CSCSegment> cscSeg, ESHandle<CSCGeometry> cscGeom,
                                           std::vector<SimSegment> seg) {

  for (int i=0; i<2; i++) {
    for (int j=0; j<5; j++) {
      for (int k=0; k<4; k++) {
         MEPath[i][j][k]=false;
	 ME_phi[i][j][k]=99.;
	 ME_eta[i][j][k]=99.;
      }
    }
  }
  for (std::vector<CSCSegment>::const_iterator it = cscSeg.begin(); it != cscSeg.end(); it++) {
      CSCDetId DetId = (CSCDetId)(*it).cscDetId();
      const CSCChamber* cscchamber = cscGeom->chamber( DetId );
      GlobalPoint  gp = cscchamber->toGlobal((*it).localPosition() );
      GlobalVector gv = cscchamber->toGlobal((*it).localDirection() );
      int st = DetId.station();
      int rg = DetId.ring();
      if (rg==4){ rg=1; }
      if ( st==1 && rg ==1 ) { st = 0; }

      // for single-chamber segment in csc
      double ab = ( gv.x()*gp.x() ) + ( gv.y()*gp.y() ) ;
      double al = sqrt(  (gp.x()*gp.x()) + (gp.y()*gp.y()) );
      double bl = sqrt(  (gv.x()*gv.x()) + (gv.y()*gv.y()) );
      double axb = ( gp.x()*gv.y() ) - ( gp.y()*gv.x() ) ;
      double cc = (axb < 0.) ? 1.0 : -1.0;

      ME_phi[1][st][rg] = cc*acos(ab /(al*bl));
      if ( ME_phi[1][st][rg] > 1.570796 ) {
         ME_phi[1][st][rg] = 3.141592 - ME_phi[1][st][rg];
      }
      ME_eta[1][st][rg] = fabs( gp.eta() );
      MEPath[1][st][rg] = true;
  }
  for (std::vector<SimSegment>::const_iterator it = seg.begin(); it != seg.end(); it++) {
      int st = ((*it).csc_DetId).station();
      int rg = ((*it).csc_DetId).ring();
      if (rg==4){ rg=1; }
      if ( st==1 && rg ==1 ) { st = 0; }
 
      // for single-chamber segment in csc
      double ab = ( ((*it).sGlobalVec).x()*((*it).sGlobalOrg).x() ) +
                  ( ((*it).sGlobalVec).y()*((*it).sGlobalOrg).y() ) ;
      double al = sqrt(  (((*it).sGlobalOrg).x()*((*it).sGlobalOrg).x()) + 
                         (((*it).sGlobalOrg).y()*((*it).sGlobalOrg).y()) );
      double bl = sqrt(  (((*it).sGlobalVec).x()*((*it).sGlobalVec).x()) + 
                         (((*it).sGlobalVec).y()*((*it).sGlobalVec).y()) );
      double axb = ( ((*it).sGlobalOrg).x()*((*it).sGlobalVec).y() ) -
                   ( ((*it).sGlobalOrg).y()*((*it).sGlobalVec).x() ) ;
      double cc = (axb < 0.) ? 1.0 : -1.0;
      ME_phi[0][st][rg] = cc*acos(ab /(al*bl));

      if ( ME_phi[0][st][rg] > 1.570796 ) {
         ME_phi[0][st][rg] = 3.141592 - ME_phi[0][st][rg];
      }
      ME_eta[0][st][rg] = fabs(((*it).sGlobalOrg).eta());
      MEPath[0][st][rg] = true;
  }
}

void SeedParametrization::FromDTSingleSeg( std::vector<DTRecSegment4D> dtSeg, ESHandle<DTGeometry> dtGeom,
                                          std::vector<SimSegment> seg) {

  for (int i=0; i<2; i++) {
    for (int j=0; j<5; j++) {
      for (int k=0; k<3; k++) {
         MBPath[i][j][k]=false;
	 MB_phi[i][j][k]=99.;
	 MB_eta[i][j][k]=99.;
      }
    }
  }
  for (std::vector<DTRecSegment4D>::const_iterator it = dtSeg.begin(); it != dtSeg.end(); it++) {

      DTChamberId DetId = (*it).chamberId();
      const DTChamber* dtchamber = dtGeom->chamber( DetId );
      int st = DetId.station();
      int wl = abs(DetId.wheel());

      MBPath[1][st][wl]=true;
      if (st==4) continue;

      GlobalPoint g_seg_o = dtchamber->toGlobal((*it).localPosition() );
      GlobalVector g_seg_v = dtchamber->toGlobal((*it).localDirection() );

      // for single-chamber segment in dt
      double ab = ( g_seg_v.x()*g_seg_o.x() ) + ( g_seg_v.y()*g_seg_o.y() ) ;
      double al = sqrt(  (g_seg_o.x()*g_seg_o.x()) + (g_seg_o.y()*g_seg_o.y()) );
      double bl = sqrt(  (g_seg_v.x()*g_seg_v.x()) + (g_seg_v.y()*g_seg_v.y()) );
      double axb = ( g_seg_o.x()*g_seg_v.y() ) - ( g_seg_o.y()*g_seg_v.x() ) ;
      double cc = (axb < 0.) ? 1.0 : -1.0;

      MB_phi[1][st][wl] = cc*acos(ab /(al*bl));
      MB_eta[1][st][wl] = fabs( g_seg_o.eta() );
  }
  /// Single segments dphi
  for (std::vector<SimSegment>::const_iterator it = seg.begin(); it != seg.end(); it++) {
      int st = ((*it).dt_DetId).station();
      int wl = abs( ((*it).dt_DetId).wheel());

      MBPath[0][st][wl]=true;
      if (st==4) continue;
   
      // for single-chamber segment in dt
      double ab = ( ((*it).sGlobalVec).x()*((*it).sGlobalOrg).x() ) +
                  ( ((*it).sGlobalVec).y()*((*it).sGlobalOrg).y() ) ;
      double al = sqrt(  (((*it).sGlobalOrg).x()*((*it).sGlobalOrg).x()) + 
                         (((*it).sGlobalOrg).y()*((*it).sGlobalOrg).y()) );
      double bl = sqrt(  (((*it).sGlobalVec).x()*((*it).sGlobalVec).x()) + 
                         (((*it).sGlobalVec).y()*((*it).sGlobalVec).y()) );
      double axb = ( ((*it).sGlobalOrg).x()*((*it).sGlobalVec).y() ) -
                   ( ((*it).sGlobalOrg).y()*((*it).sGlobalVec).x() ) ;
      double cc = (axb < 0.) ? 1.0 : -1.0;
      MB_phi[0][st][wl] = cc*acos(ab /(al*bl));
      MB_eta[0][st][wl] = fabs( ((*it).sGlobalOrg).eta() );
  }
}     