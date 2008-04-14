// Class Header
#include "SeedGeneratorTest.h"
#include "RecoMuon/SeedGenerator/test/MuonSeedPTAnalysis/SegSelector.h"
#include "RecoMuon/SeedGenerator/test/MuonSeedPTAnalysis/MuonSeedBuilder.h"
#include "RecoMuon/SeedGenerator/test/MuonSeedPTAnalysis/SeedBuilder.h"

// for MuonSeedBuilder
#include "RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h"
#include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"

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


DEFINE_FWK_MODULE(SeedGeneratorTest);
using namespace std;
using namespace edm;

// constructors
SeedGeneratorTest::SeedGeneratorTest(const ParameterSet& pset){ 

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
  muonSeedBuilder_  = new MuonSeedBuilder( pset );
  SeedBuilder_      = new SeedBuilder( pset );

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

  // All possible segment pair in CSC 
  //                  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 
  int csc1[2][17]= {{11,11,12,12,13,11,11,12,12,13,11,12,21,21,22,21,31} ,
                    {21,22,21,22,22,31,32,31,32,32,41,41,31,32,32,41,41} }; 
  char ME_nu1[8];
  for (int i =0; i<17; i++) {
      sprintf(ME_nu1,"ME_%d-%d", csc1[0][i] , csc1[1][i] );
      hME1[i] = new H2DRecHit4(ME_nu1);
      cout <<"hME1_"<<i<<" = "<< ME_nu1<<endl;
  }
  // All possible segment pair in DT 
  //                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
  int dt1[2][30]={{10,10,10,10,10,10,11,11,11,11,11,11,12,12,12,20,20,20,20,21,21,21,21,22,22,30,30,31,31,32}
                 ,{20,21,30,31,40,41,21,22,31,32,41,42,22,32,42,30,31,40,41,31,32,41,42,32,42,40,41,41,42,42}};
  char MB_nu1[8];
  for (int i =0; i<30; i++) {
      sprintf(MB_nu1,"MB_%d-%d", dt1[0][i] , dt1[1][i] );
      hMB1[i] = new H2DRecHit5(MB_nu1);
      cout <<"hMB1_"<<i<<" = "<< MB_nu1<<endl;
  }
  // All possible segment pair in Overlap region between DT and CSC
  //                0  1  2  3  4  5  6  7  8  9 10 
  int olp[2][11]={{11,12,12,12,12,21,22,22,22,32,32} 
                 ,{13,12,13,22,32,13,13,22,32,13,22}};
  char OL_nu1[8];
  for (int i =0; i<11; i++) {
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
  //             0  1  2  3  4  5 
  int csc3[6] ={12,13,14,22,23,34};
  char ME_nu3[6];
  for (int i =0; i<6; i++) {
      sprintf(ME_nu3,"RME_%d", csc3[i] );
      hME3[i] = new H2DRecHit8(ME_nu3);
      cout <<"hME3_"<<i<<" = "<< ME_nu3<<endl;
  }

  //            0  1  2  3  4  5 
  int dt3[6] ={12,13,14,23,24,34};
  char MB_nu3[6];
  for (int i =0; i<6; i++) {
      sprintf(MB_nu3,"RMB_%d", dt3[i] );
      hMB3[i] = new H2DRecHit9(MB_nu3);
      cout <<"hMB3_"<<i<<" = "<< MB_nu3<<endl;
  }

  h_all  = new H2DRecHit1("AllMu_");
  h_csc  = new H2DRecHit2("CSC_A");
  h_dt   = new H2DRecHit3("DT_A");

}

// destructor
SeedGeneratorTest::~SeedGeneratorTest(){

  if (debug) cout << "[SeedQualityAnalysis] Destructor called" << endl;
  delete recsegSelector;
  delete muonSeedBuilder_; 
  delete SeedBuilder_; 
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

  theFile->cd();
  theFile->cd("ME_All");
  for (int i=0; i<17; i++) {
      hME1[i]->Write();
  }
  for (int i=0; i<8; i++) {
      hME2[i]->Write();
  }
  for (int i=0; i<6; i++) {
      hME3[i]->Write();
  }
  theFile->cd();
  theFile->cd("MB_All");
  for (int i=0; i<30; i++) {
      hMB1[i]->Write();
  }
  for (int i=0; i<12; i++) {
      hMB2[i]->Write();
  }
  for (int i=0; i<6; i++) {
      hMB3[i]->Write();
  }
  theFile->cd();
  theFile->cd("OL_All");
  for (int i=0; i<11; i++) {
      hOL1[i]->Write();
  }
  // for tree
  theFile->cd();
  tr_muon->Write();

  // Release the memory ...
  for (int i=0; i<17; i++) {
      delete hME1[i];
  }
  for (int i=0; i<30; i++) {
      delete hMB1[i];
  }
  for (int i=0; i<11; i++) {
      delete hOL1[i];
  }
  for (int i=0; i<8; i++) {
      delete hME2[i];
  }
  for (int i=0; i<12; i++) {
      delete hMB2[i];
  }
  for (int i=0; i<6; i++) {
      delete hME3[i];
  }
  for (int i=0; i<6; i++) {
      delete hMB3[i];
  }
  delete h_all;
  delete h_csc;
  delete h_dt;
  delete tr_muon;

  theFile->Close();
  if (debug) cout << "************* Finished writing histograms to file" << endl;

}

// The Main...Aanlysis...

void SeedGeneratorTest::analyze(const Event& event, const EventSetup& eventSetup)
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

  // Get the muon seeds
  Handle<TrajectorySeedCollection> muonseeds;
  event.getByLabel(muonseedLabel, muonseeds);

  // Muon Geometry for testing DF seedGenerator
  ESHandle<MuonDetLayerGeometry> muonLayers;
  eventSetup.get<MuonRecoGeometryRecord>().get(muonLayers);
  const MuonDetLayerGeometry* lgeom = &*muonLayers;
  muonSeedBuilder_->setGeometry( lgeom );
  SeedBuilder_->setGeometry( lgeom );

  // Magnetic field
  ESHandle<MagneticField> field;
  eventSetup.get<IdealMagneticFieldRecord>().get(field);
  const MagneticField* theField = &*field;
  muonSeedBuilder_->setBField( theField );
  SeedBuilder_->setBField( theField );
 
  H2DRecHit1 *histo1 = 0;   
  H2DRecHit2 *histo2 = 0;   
  H2DRecHit3 *histo3 = 0;   
  H2DRecHit4 *histo4 = 0;   
  H2DRecHit5 *histo5 = 0;   
  H2DRecHit6 *histo6 = 0;   
  H2DRecHit7 *histo7 = 0;   
  H2DRecHit8 *histo8 = 0;   
  H2DRecHit9 *histo9 = 0;   
  H2DRecHit10 *histo10 = 0;   
  TNtuple1 *tt = 0;
 
  // 0. Run the class SegSelector
  //SegSelector recsegSelector();
  std::vector<SimSegment> sCSC_v      = recsegSelector->Sim_CSCSegments(csimHits,cscGeom);
  std::vector<SimSegment> sDT_v       = recsegSelector->Sim_DTSegments(dsimHits,dtGeom);
  std::vector<CSCSegment> cscseg_V    = recsegSelector->Select_CSCSeg(cscSegments,cscGeom, sCSC_v);
  std::vector<DTRecSegment4D> dtseg_V = recsegSelector->Select_DTSeg(dt4DSegments,dtGeom, sDT_v);
  
  // MuonSeedBuilder... Create pointer to the seed container
  
  std::auto_ptr<TrajectorySeedCollection> output(new TrajectorySeedCollection() );
  int nSeeds = 0;
  nSeeds = muonSeedBuilder_->build( event, *output);
  std::vector<int> nSegFromSeed  = muonSeedBuilder_->get_nSegOnSeed();
  std::vector<float> etaFromSeed = muonSeedBuilder_->get_etaOfSeed();
  
  // Modified DFSeedBuilder... Create pointer to the seed container
  std::auto_ptr<TrajectorySeedCollection> output1(new TrajectorySeedCollection() );
  int nSeeds1 = 0;
  nSeeds1 = SeedBuilder_->build( event, eventSetup, *output1);
  std::vector<int> nSegFromSeed1  = SeedBuilder_->get_nSegOnSeed();
  std::vector<float> etaFromSeed1 = SeedBuilder_->get_etaOfSeed();
  // reading the information from the damn seeds
  nSegInSeed = 0;
  RecSeedReader(muonseeds,dtGeom,cscGeom);

  histo1 = h_all;
  if (nSegFromSeed.size()> 0) {
     histo1->Fill1e(nSegFromSeed[0], fabs(etaFromSeed[0]));
     histo1->Fill1h(nSegFromSeed[0], fabs(seed_gp.eta()) );
  }
  if (nSegFromSeed1.size()> 0) {
     histo1->Fill1g(nSegFromSeed1[0],fabs(etaFromSeed1[0]));
  }

  // 1. Get the track information  
  double pt =0.0; 
  double pt1[5] = {0.0}; 
  double pa =0; 
  int tk_size=0;
  double eta_c = -9.0;
  double eta_d = -9.0;
  double eta_trk = -9.0;
  for (SimTrackContainer::const_iterator simTk_It = simTracks->begin(); simTk_It != simTracks->end(); simTk_It++)
  {
      tk_size++ ;
      if (abs((*simTk_It).type())!=13) continue;
      float px = ((*simTk_It).momentum()).x();
      float py = ((*simTk_It).momentum()).y();
      float pz = ((*simTk_It).momentum()).z();
      pt = sqrt((px*px) + (py*py));
      pa = sqrt((px*px) + (py*py) + (pz*pz));
      double theta = acos( pz/pa );
      eta_trk = fabs((-1.0)*log( tan(theta/2.0) ));
      pt1[0] = pt;

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
          //const DTChamber* dtchamber = dtGeom->chamber( D_Id );
          GlobalVector m2 = dtlayer->toGlobal((*ds_It).momentumAtEntry() );
          GlobalPoint gp = dtlayer->toGlobal(lp );
          //Local3DPoint  lp1 = dtchamber->toLocal(gp);

          //DTSuperLayerId S_Id = DTSuperLayerId( (*ds_It).detUnitId() );
          //const DTSuperLayer* dtsuperlayer = dtGeom->superLayer(S_Id);
          //GlobalPoint gp2 = dtsuperlayer->toGlobal(lp );
          //Local3DPoint  lp2 = dtsuperlayer->toLocal(gp2);

          //const DTChamber* dtchamber = dtGeom->chamber( D_Id );
          //GlobalPoint gp = dtchamber->toGlobal(lp );
          //LocalPoint  lpa = dtchamber->toLocal(gp);

          if ( ( abs((*ds_It).particleType())==13 ) && ( (*ds_It).trackId()==1 )) {      
             pt1[D_Id.station()] = sqrt( (m2.x()*m2.x()) + (m2.y()*m2.y()) );
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
  /// eta test !!!!
  Eta_Test(simTracks,dsimHits,dtGeom);
  histo1 = h_all;
  histo1->Fill1c(h_trk,h_sl1[1],h_sl2[1],h_sl3[1],h_sl1[2],h_sl2[2],h_sl3[2],h_sl1[3],h_sl2[3],h_sl3[3]);

  // 2. Check # of segments in each chambers
  CSCsegment_stat(cscSegments);
  DTsegment_stat(dt4DSegments);
  if ( (cscseg_stat[5] < 2)&&(dtseg_stat[5] < 2) ) {
     SeedFromCSCRecHit(csc2DRecHits,cscGeom);
     SeedFromDTRecHit(dt1DRecHits,dtGeom);
  }
  // 3. Event statistic
  int allmu_stat = cscseg_stat[5]+dtseg_stat[5];
  int allrh_stat = cscrh_sum[5]+dtrh_sum[5];
  double allmu_eta = -9.0;
  if      ((cscseg_stat[5]==0)&&(dtseg_stat[5]!=0)) { allmu_eta = eta_d; }
  else if ((dtseg_stat[5]==0)&&(cscseg_stat[5]!=0)) { allmu_eta = eta_c; }
  else if ((cscseg_stat[5]!=0)&&(dtseg_stat[5]!=0)) { allmu_eta = (eta_c + eta_d)/2 ; }
  else {
       if ((eta_d == -9.0)&&(eta_c != -9.0)) {allmu_eta = eta_c ; }
       else if ((eta_d != -9.0)&&(eta_c == -9.0)) {allmu_eta = eta_d ; }
       else { allmu_eta = (eta_c + eta_d)/2 ; }
       histo1 = h_all;
       histo1->Fill1a(allrh_stat,allmu_eta,cscrh_sum[0],dtrh_sum[0]);
  }
    
  histo1 = h_all;
  histo1->Fill1(cscseg_stat[5],dtseg_stat[5],eta_c,eta_d,allmu_stat,allmu_eta,eta_trk);
  
  //cout <<" csc eta= "<< eta_c <<"  #= " <<cscseg_stat[0];
  //cout <<" dt eta= "<< eta_d <<"  #= " << dtseg_stat[0]<<endl;

  if (cscseg_stat[0] != 0){
     histo2 = h_csc;
     histo2->Fill3(pt,pa,cscseg_stat[0]);
     for (int k=1; k<5; k++) {
         histo2->Fill4(k,cscseg_stat[k],cscseg_stat1[k]);
     }
  }
  if (dtseg_stat[0] != 0) {
     histo3 = h_dt;
     histo3->Fill3a(pt,pa,dtseg_stat[0]); 
     for (int k=1; k<5; k++) {
         histo3->Fill4a(k,dtseg_stat[k],dtseg_stat1[k]);
     }
  }
  bool overlap = false;
  if ( (dtseg_stat[0] != 0)&&(cscseg_stat[0] != 0) ){
     overlap = true;
  }

  // 4. build the simulated segments
  //    output sCSC_v & sDT_v 

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
  if      ((simcscseg[0]==0)&&(simdtseg[0]!=0)) { allmu_eta1 = eta_sim2; }
  else if ((simdtseg[0]==0)&&(simcscseg[0]!=0)) { allmu_eta1 = eta_sim1; }
  else { allmu_eta1 = (eta_sim1 + eta_sim2)/2 ; }
    
  histo1 = h_all;
  histo1->Fill1b(allmu_stat1,allmu_eta1);

  // 5. match the sim-segment and rec-segment
  //    output cscseg_V & dtseg_V

  // 6. get the global phi & eta from csc simHits
  ///  For CSC
  std::vector<double> ephi(5, 999.);
  std::vector<double> eeta(5, 999.);
  GlobalVector ini_gv = GlobalVector( 0., 0., 0. );
  std::vector<GlobalVector> egv(5, ini_gv);
  bool path[5][4]={{false},{false}};
  bool path1[5]={false};
  ///  Single segment dphi
  double ME_phi[5][4]={{999.},{999.}};
  double ME_eta[5][4]={{999.},{999.}};
  for (std::vector<SimSegment>::const_iterator it = sCSC_v.begin(); it != sCSC_v.end(); it++) {
      int st = ((*it).csc_DetId).station();
      int rg = ((*it).csc_DetId).ring();
      if (rg==4){
         rg=1;
      }
 
      if (st==1){
         if (!path1[1]) {
            ephi[st]= ((*it).sGlobalOrg).phi(); 
            eeta[st]= ((*it).sGlobalOrg).eta();
            path[st][rg]=true;
            path1[st]=true;
            egv[st] = (*it).sGlobalVec ;
         }
      }

      if (st!=1) {
            ephi[st]= ((*it).sGlobalOrg).phi(); 
            eeta[st]= ((*it).sGlobalOrg).eta();
            path[st][rg]=true;
            path1[st]=true;
            egv[st] = (*it).sGlobalVec ;
      }
 
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
      ME_phi[st][rg] = cc*acos(ab /(al*bl));
      if ( ME_phi[st][rg] > 1.570796 ) {
         ME_phi[st][rg] = 3.141592 - ME_phi[st][rg];
      }
      ME_eta[st][rg] = fabs(((*it).sGlobalOrg).eta());

  }
  // For DT
  /// Single segments dphi
  std::vector<double> bphi(5, 999.);
  std::vector<double> beta(5, 999.);
  std::vector<GlobalVector> bgv(5, ini_gv);
  bool bpath[5][4]={{false},{false}};
  bool bpath1[5]={false};
  double MB_phi[5][4]={{999.},{999.}};
  double MB_eta[5][4]={{999.},{999.}};
  for (std::vector<SimSegment>::const_iterator it = sDT_v.begin(); it != sDT_v.end(); it++) {
      int st = ((*it).dt_DetId).station();
      int wl = abs( ((*it).dt_DetId).wheel());
      bphi[st]= ((*it).sGlobalOrg).phi(); 
      beta[st]= ((*it).sGlobalOrg).eta();
      bpath[st][wl]=true;
      bpath1[st]=true;
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
      MB_phi[st][wl] = cc*acos(ab /(al*bl));
      MB_eta[st][wl] = fabs( ((*it).sGlobalOrg).eta() );

      bgv[st] = (*it).sGlobalVec ;
  }

  /// get the d_phi and d_eta between each csc station
  /// segment pairs
  double d_ephi[4][5]={{999.},{999.}};
  double d_eeta[4][5]={{999.},{999.}};
  for (int l=0; l<12; l++) {
      int s1 = (l%3)+1; // 1,2,3  1,2,3  1,2,3  1,2,3
      int s2 = (l/3)+1; // 1,1,1  2,2,2  3,3,3  4,4,4
      d_ephi[s1][s2] = ephi[s1]-ephi[s2];  
      d_eeta[s1][s2] = eeta[s1]-eeta[s2];
  }

  // get the d_phi and d_eta between each dt station
  /// segment pairs
  double d_bphi[4][5]={{999.},{999.}};
  double d_beta[4][5]={{999.},{999.}};
  for (int l=0; l<12; l++) {
      int s1 = (l%3)+1; // 1,2,3  1,2,3  1,2,3  1,2,3
      int s2 = (l/3)+1; // 1,1,1  2,2,2  3,3,3  4,4,4
      d_bphi[s1][s2] = bphi[s1]-bphi[s2];
      d_beta[s1][s2] = beta[s1]-beta[s2];
  }

  // get the d_phi and d_eta in overlap region
  /// segment pairs
  double d_ophi[4][5]={{999.},{999.}};
  double d_oeta[4][5]={{999.},{999.}};
  for (int l=0; l<12; l++) {
      int s1 = (l%3)+1; // 1,2,3  1,2,3  1,2,3  1,2,3
      int s2 = (l/3)+1; // 1,1,1  2,2,2  3,3,3  4,4,4
                        // ^ ^ ^  ^ ^ ^  ^ ^   marking the effective elements!
      d_ophi[s1][s2] = bphi[s1]-ephi[s2];
      d_oeta[s1][s2] = beta[s1]-eeta[s2];
  }

  
  /// get the global dphi from recHits for CSC
  std::vector<double> ephi_r(5, 999.);
  std::vector<double> eeta_r(5, 999.);
  std::vector<GlobalVector> egv_r(5, ini_gv);
  bool path_r[5][4]={{false},{false}};
  bool path_r1[5]={false};
  double ME_phi_r[5][4]={{999.},{999.}};
  double ME_eta_r[5][4]={{999.},{999.}};
  double eta_csc2 = 0.0;
  for (std::vector<CSCSegment>::const_iterator it = cscseg_V.begin(); it != cscseg_V.end(); it++) {
      CSCDetId DetId = (CSCDetId)(*it).cscDetId();
      const CSCChamber* cscchamber = cscGeom->chamber( DetId );
      GlobalPoint g_seg_o = cscchamber->toGlobal((*it).localPosition() );
      GlobalVector g_seg_v = cscchamber->toGlobal((*it).localDirection() );
      int st = DetId.station();
      int rg = DetId.ring();
      eta_csc2 += fabs(g_seg_o.eta());
      if (rg==4){
         rg=1;
      }
      if (st==1) {
         if (!path_r1[1]) {
            ephi_r[st]= g_seg_o.phi(); 
            eeta_r[st]= g_seg_o.eta();
            path_r[st][rg]=true;
            path_r1[st]=true;
            egv_r[st] = g_seg_v ;
         }
      }
      if (st!=1) {
         ephi_r[st]= g_seg_o.phi(); 
         eeta_r[st]= g_seg_o.eta();
         path_r[st][rg]=true;
         path_r1[st]=true;
         egv_r[st] = g_seg_v ;
      }     

      // for single-chamber segment in csc
      double ab = ( g_seg_v.x()*g_seg_o.x() ) + ( g_seg_v.y()*g_seg_o.y() ) ;
      double al = sqrt(  (g_seg_o.x()*g_seg_o.x()) + (g_seg_o.y()*g_seg_o.y()) );
      double bl = sqrt(  (g_seg_v.x()*g_seg_v.x()) + (g_seg_v.y()*g_seg_v.y()) );
      double axb = ( g_seg_o.x()*g_seg_v.y() ) - ( g_seg_o.y()*g_seg_v.x() ) ;
      double cc = (axb < 0.) ? 1.0 : -1.0;
      ME_phi_r[st][rg] = cc*acos(ab /(al*bl));
      if ( ME_phi_r[st][rg] > 1.57 ) {
         ME_phi_r[st][rg] = 3.141592 - ME_phi_r[st][rg];
      }
      ME_eta_r[st][rg] = fabs( g_seg_o.eta() );
  }
  if (cscseg_V.size()!=0){ 
     eta_csc2 =  eta_csc2/cscseg_V.size() ;
     histo1 = h_all;
     histo1->Fill1f(cscseg_V.size(),eta_csc2,eta_c);
  }

  /// get the global dphi from recHits for DT
  std::vector<double> bphi_r(5, 999.);
  std::vector<double> beta_r(5, 999.);
  std::vector<GlobalVector> bgv_r(5, ini_gv);
  bool bpath_r[5][4]={{false},{false}};
  bool bpath_r1[5]={false};
  double MB_phi_r[5][4]={{999.},{999.}};
  double MB_eta_r[5][4]={{999.},{999.}};
  double eta_dt2 = 0.0;
  double nu_eta_dt2 = 0.0;
  for (std::vector<DTRecSegment4D>::const_iterator it = dtseg_V.begin(); it != dtseg_V.end(); it++) {

      DTChamberId DetId = (*it).chamberId();
      const DTChamber* dtchamber = dtGeom->chamber( DetId );

      GlobalPoint g_seg_o = dtchamber->toGlobal((*it).localPosition() );
      GlobalVector g_seg_v = dtchamber->toGlobal((*it).localDirection() );

      /*DTLayerId DetId = DTLayerId((*it).chamberId());
      const DTLayer* dtlayer = dtGeom->layer(DetId);

      GlobalPoint g_seg_o = dtlayer->toGlobal((*it).localPosition() );
      GlobalVector g_seg_v = dtlayer->toGlobal((*it).localDirection() );
      */

      int st = DetId.station();
      int wl = abs(DetId.wheel());
      bphi_r[st]= g_seg_o.phi(); 
      beta_r[st]= g_seg_o.eta();
      if (DetId.station()!= 4) {
        eta_dt2 += fabs(g_seg_o.eta());
        nu_eta_dt2 += 1.0;
      }

      bpath_r[st][wl]=true;
      bpath_r1[st]=true; 
      bgv_r[st] = g_seg_v ;
      // for single-chamber segment in dt
      double ab = ( g_seg_v.x()*g_seg_o.x() ) + ( g_seg_v.y()*g_seg_o.y() ) ;
      double al = sqrt(  (g_seg_o.x()*g_seg_o.x()) + (g_seg_o.y()*g_seg_o.y()) );
      double bl = sqrt(  (g_seg_v.x()*g_seg_v.x()) + (g_seg_v.y()*g_seg_v.y()) );
      double axb = ( g_seg_o.x()*g_seg_v.y() ) - ( g_seg_o.y()*g_seg_v.x() ) ;
      double cc = (axb < 0.) ? 1.0 : -1.0;
      MB_phi_r[st][wl] = cc*acos(ab /(al*bl));
      MB_eta_r[st][wl] = fabs( g_seg_o.eta() );
  }

  if (nu_eta_dt2!=0) {
     eta_dt2 = eta_dt2/nu_eta_dt2 ;
     histo1 = h_all;
     histo1->Fill1d(dtseg_V.size(),eta_dt2,eta_d);
  }
  else {
     eta_dt2 = -9.0;
  }

  /// For reco-segment tree
  tt = tr_muon;
  if ( path_r1[1]&&path_r1[2]&&path_r1[3] ) {
      tt->Fill_b1(fabs(eeta_r[1]),fabs(eeta_r[2]),fabs(eeta_r[3]),fabs(eeta_r[4]), 
                  ephi_r[1], ephi_r[2], ephi_r[3], ephi_r[4],pt);
      tt->Fill_l1(pa);
  }
  if ( bpath_r1[1]&&bpath_r1[2]&&bpath_r1[3] ) {
      tt->Fill_b2(fabs(beta_r[1]),fabs(beta_r[2]),fabs(beta_r[3]),fabs(beta_r[4]), 
                  bphi_r[1], bphi_r[2], bphi_r[3], bphi_r[4],pt);
      tt->Fill_l1(pa);
  }
  tt->FillTree();

   
  /// reco-dphi for CSC
  double d_ephi_r[4][5]={{999.},{999.}};
  double d_eeta_r[4][5]={{999.},{999.}};
  for (int l=0; l<12; l++) {
      int s1 = (l%3)+1; // 1,2,3,1,2,3,1,2,3,1,2,3
      int s2 = (l/3)+1; // 1,1,1,2,2,2,3,3,3,4,4,4
      d_ephi_r[s1][s2] = ephi_r[s1]-ephi_r[s2];  
      d_eeta_r[s1][s2] = eeta_r[s1]-eeta_r[s2];  
  }
  /// reco-dphi for DT
  double d_bphi_r[4][5]={{999.},{999.}};
  double d_beta_r[4][5]={{999.},{999.}};
  for (int l=0; l<12; l++) {
      int s1 = (l%3)+1; // 1,2,3,1,2,3,1,2,3,1,2,3
      int s2 = (l/3)+1; // 1,1,1,2,2,2,3,3,3,4,4,4
      d_bphi_r[s1][s2] = bphi_r[s1]-bphi_r[s2];  
      d_beta_r[s1][s2] = beta_r[s1]-beta_r[s2];  
  }
  /// reco-dphi for Overlap
  double d_ophi_r[4][5]={{999.},{999.}};
  double d_oeta_r[4][5]={{999.},{999.}};
  for (int l=0; l<12; l++) {
      int s1 = (l%3)+1; // 1,2,3,1,2,3,1,2,3,1,2,3
      int s2 = (l/3)+1; // 1,1,1,2,2,2,3,3,3,4,4,4
      d_ophi_r[s1][s2] = bphi_r[s1]-ephi_r[s2];  
      d_oeta_r[s1][s2] = beta_r[s1]-eeta_r[s2];  
  }
  

  //  look at different Bxdl btw. stations
  /// For CSC
  for (int l=0; l<12; l++) {
      int s1 = (l%3)+1; // 1,2,3  1,2,3  1,2,3  1,2,3
      int s2 = (l/3)+1; // 1,1,1  2,2,2  3,3,3  4,4,4
      
      if (path1[1] && path1[2] &&  path_r1[1] && path_r1[2] && (s1==1) && (s2==2) ){
         histo2 = h_csc;
         histo2->Fill5_1( d_ephi[1][2],d_ephi_r[1][2],pt1[1]*d_ephi[1][2],pt1[1]*d_ephi_r[1][2],
                        fabs(eeta[2]),fabs(eeta_r[2]) );
      }
      if (path1[1] && path1[3] &&  path_r1[1] && path_r1[3] && (s1==1) && (s2==3) ){
         histo2 = h_csc;
         histo2->Fill5_2( d_ephi[1][3],d_ephi_r[1][3],pt1[1]*d_ephi[1][3],pt1[1]*d_ephi_r[1][3],
                        fabs(eeta[3]),fabs(eeta_r[3]) );
      }
      if (path1[1] && path1[4] &&  path_r1[1] && path_r1[4] && (s1==1) && (s2==4) ){
         histo2 = h_csc;
         histo2->Fill5_3( d_ephi[1][4],d_ephi_r[1][4],pt1[1]*d_ephi[1][4],pt1[1]*d_ephi_r[1][4],
                        fabs(eeta[4]),fabs(eeta_r[4]) );
      }
      if (path1[2] && path1[3] &&  path_r1[2] && path_r1[3] && (s1==2) && (s2==3) ){
         histo2 = h_csc;
         histo2->Fill5_4( d_ephi[2][3],d_ephi_r[2][3],pt1[2]*d_ephi[2][3],pt1[2]*d_ephi_r[2][3],
                        fabs(eeta[3]),fabs(eeta_r[3]) );
      }
      if (path1[2] && path1[4] &&  path_r1[2] && path_r1[4] && (s1==2) && (s2==4) ){
         histo2 = h_csc;
         histo2->Fill5_5( d_ephi[2][4],d_ephi_r[2][4],pt1[2]*d_ephi[2][4],pt1[2]*d_ephi_r[2][4],
                        fabs(eeta[4]),fabs(eeta_r[4]) );
      }
      if (path1[3] && path1[4] &&  path_r1[3] && path_r1[4] && (s1==3) && (s2==4) ){
         histo2 = h_csc;
         histo2->Fill5_6( d_ephi[3][4],d_ephi_r[3][4],pt1[3]*d_ephi[3][4],pt1[3]*d_ephi_r[3][4],
                        fabs(eeta[4]),fabs(eeta_r[4]) );
      }



  }
  if (( ephi[1]!=999.)&&(ephi[2]!=999.)&&( ephi[3]!=999.) &&
      ( ephi_r[1]!=999.)&&(ephi_r[2]!=999.)&&( ephi_r[3]!=999.) ){
     histo2 = h_csc;
     histo2->Fill5b(fabs(eeta[2]),(d_ephi[2][3]/d_ephi[1][2]),fabs(eeta_r[2]),(d_ephi_r[2][3]/d_ephi_r[1][2]),pt);
  }
  if (( ephi[4]!=999.)&&(ephi[2]!=999.)&&( ephi[3]!=999.) &&
      ( ephi_r[4]!=999.)&&(ephi_r[2]!=999.)&&( ephi_r[3]!=999.) ){
     histo2 = h_csc;
     histo2->Fill5b1(fabs(eeta[2]),(d_ephi[3][4]/d_ephi[2][3]),fabs(eeta_r[2]),(d_ephi_r[3][4]/d_ephi_r[2][3]),pt);
  }

  /// For CSC
  //// CSC dphi ratio from all segment pairs btw station 12 and station 23
  for(CSCSegmentCollection::const_iterator i1 = cscSegments->begin(); i1 != cscSegments->end(); i1++) {
     double dphi_12 = 0.0 ;
     double dphi_23 = 0.0 ;
     histo2 = h_csc;
     CSCDetId DetId1 = (CSCDetId)(*i1).cscDetId();
     const CSCChamber* cscchamber1 = cscGeom->chamber( DetId1 );
     GlobalPoint  gp1 = cscchamber1->toGlobal( (*i1).localPosition() );
     if(DetId1.station()==1){
       for(CSCSegmentCollection::const_iterator i2 = cscSegments->begin(); i2 != cscSegments->end(); i2++) {
          CSCDetId DetId2 = (CSCDetId)(*i2).cscDetId();
          const CSCChamber* cscchamber2 = cscGeom->chamber( DetId2 );
          GlobalPoint  gp2 = cscchamber2->toGlobal( (*i2).localPosition() );
          if (DetId2.station()==2) { 
             dphi_12 = gp1.phi() - gp2.phi() ;
             for(CSCSegmentCollection::const_iterator i3 = cscSegments->begin(); i3 != cscSegments->end(); i3++) {
                CSCDetId DetId3 = (CSCDetId)(*i3).cscDetId();
                const CSCChamber* cscchamber3 = cscGeom->chamber( DetId3 );
                GlobalPoint  gp3 = cscchamber3->toGlobal( (*i3).localPosition() );
                if (DetId3.station()==3) { 
                   dphi_23 = gp2.phi() - gp3.phi() ;
                   histo2->Fill5c( fabs(gp1.eta()) , (dphi_23/dphi_12), pt );
                }
             }
          }
       }
     }
  }  
  //// CSC dphi ratio from all Segment pairs btw station 23 and station 34
  for(CSCSegmentCollection::const_iterator i1 = cscSegments->begin(); i1 != cscSegments->end(); i1++) {
     double dphi_23 = 0.0 ;
     double dphi_34 = 0.0 ;
     histo2 = h_csc;
     CSCDetId DetId1 = (CSCDetId)(*i1).cscDetId();
     const CSCChamber* cscchamber1 = cscGeom->chamber( DetId1 );
     GlobalPoint  gp1 = cscchamber1->toGlobal( (*i1).localPosition() );
     if(DetId1.station()==2){
       for(CSCSegmentCollection::const_iterator i2 = cscSegments->begin(); i2 != cscSegments->end(); i2++) {
          CSCDetId DetId2 = (CSCDetId)(*i2).cscDetId();
          const CSCChamber* cscchamber2 = cscGeom->chamber( DetId2 );
          GlobalPoint  gp2 = cscchamber2->toGlobal( (*i2).localPosition() );
          if (DetId2.station()==3) { 
             dphi_23 = gp1.phi() - gp2.phi() ;
             for(CSCSegmentCollection::const_iterator i3 = cscSegments->begin(); i3 != cscSegments->end(); i3++) {
                CSCDetId DetId3 = (CSCDetId)(*i3).cscDetId();
                const CSCChamber* cscchamber3 = cscGeom->chamber( DetId3 );
                GlobalPoint  gp3 = cscchamber3->toGlobal( (*i3).localPosition() );
                if (DetId3.station()==4) { 
                   dphi_34 = gp2.phi() - gp3.phi() ;
                   histo2->Fill5c1( fabs(gp1.eta()) , (dphi_34/dphi_23), pt );
                }
             }
          }
       }
     }
  }  
  //// CSC dphi ratio from single segment btw different station 
  


  /// For DT
  for (int l=0; l<12; l++) {
      int s1 = (l%3)+1; // 1,2,3  1,2,3  1,2,3  1,2,3
      int s2 = (l/3)+1; // 1,1,1  2,2,2  3,3,3  4,4,4

      if (bpath1[1] && bpath1[2] &&  bpath_r1[1] && bpath_r1[2] && (s1==1) && (s2==2) ){
         histo3 = h_dt;
         histo3->Fill6_1( d_bphi[1][2],d_bphi_r[1][2],pt1[1]*d_bphi[1][2],pt1[1]*d_bphi_r[1][2],
                        fabs(beta[2]),fabs(beta_r[2]) );
      }
      if (bpath1[1] && bpath1[3] &&  bpath_r1[1] && bpath_r1[3] && (s1==1) && (s2==3) ){
         histo3 = h_dt;
         histo3->Fill6_2( d_bphi[1][3],d_bphi_r[1][3],pt1[1]*d_bphi[1][3],pt1[1]*d_bphi_r[1][3],
                        fabs(beta[3]),fabs(beta_r[3]) );
      }
      if (bpath1[1] && bpath1[4] &&  bpath_r1[1] && bpath_r1[4] && (s1==1) && (s2==4) ){
         histo3 = h_dt;
         histo3->Fill6_3( d_bphi[1][4],d_bphi_r[1][4],pt1[1]*d_bphi[1][4],pt1[1]*d_bphi_r[1][4],
                        fabs(beta[1]),fabs(beta_r[1]) );
      }
      if (bpath1[2] && bpath1[3] &&  bpath_r1[2] && bpath_r1[3] && (s1==2) && (s2==3) ){
         histo3 = h_dt;
         histo3->Fill6_4( d_bphi[2][3],d_bphi_r[2][3],pt1[2]*d_bphi[2][3],pt1[2]*d_bphi_r[2][3],
                        fabs(beta[3]),fabs(beta_r[3]) );
      }
      if (bpath1[2] && bpath1[4] &&  bpath_r1[2] && bpath_r1[4] && (s1==2) && (s2==4) ){
         histo3 = h_dt;
         histo3->Fill6_5( d_bphi[2][4],d_bphi_r[2][4],pt1[2]*d_bphi[2][4],pt1[2]*d_bphi_r[2][4],
                        fabs(beta[1]),fabs(beta_r[1]) );
      }
      if (bpath1[3] && bpath1[4] &&  bpath_r1[3] && bpath_r1[4] && (s1==3) && (s2==4) ){
         histo3 = h_dt;
         histo3->Fill6_6( d_bphi[3][4],d_bphi_r[3][4],pt1[3]*d_bphi[3][4],pt1[3]*d_bphi_r[3][4],
                        fabs(beta[1]),fabs(beta_r[1]) );
      }


  }
  if (( bphi[1]!=999.)&&(bphi[2]!=999.)&&( bphi[3]!=999.) &&
      ( bphi_r[1]!=999.)&&(bphi_r[2]!=999.)&&( bphi_r[3]!=999.) ){
     histo3 = h_dt;
     histo3->Fill7b( fabs(beta[2]),(d_bphi[2][3]/d_bphi[1][2]),fabs(beta_r[2]),(d_bphi_r[2][3]/d_bphi_r[1][2]),pt);
  }

  for(DTRecSegment4DCollection::const_iterator i1 = dt4DSegments->begin(); i1 != dt4DSegments->end(); i1++){
     double dphi_12 = 0.0 ;
     double dphi_23 = 0.0 ;
     histo3 = h_dt;
     DTChamberId DetId1 = (*i1).chamberId();
     const DTChamber* dtchamber1 = dtGeom->chamber( DetId1 );
     GlobalPoint gp1 = dtchamber1->toGlobal((*i1).localPosition() );
     if(DetId1.station()==1){
       for(DTRecSegment4DCollection::const_iterator i2 = dt4DSegments->begin(); i2 != dt4DSegments->end(); i2++){
          DTChamberId DetId2 = (*i2).chamberId();
          const DTChamber* dtchamber2 = dtGeom->chamber( DetId2 );
          GlobalPoint gp2 = dtchamber2->toGlobal((*i2).localPosition() );
          if (DetId2.station()==2) { 
             dphi_12 = gp1.phi() - gp2.phi() ;
             for(DTRecSegment4DCollection::const_iterator i3 = dt4DSegments->begin(); i3 != dt4DSegments->end(); i3++){
                DTChamberId DetId3 = (*i3).chamberId();
                const DTChamber* dtchamber3 = dtGeom->chamber( DetId3 );
                GlobalPoint gp3 = dtchamber3->toGlobal((*i3).localPosition() );
                if (DetId3.station()==3) { 
                   dphi_23 = gp2.phi() - gp3.phi() ;
                   histo3->Fill7c( fabs(gp1.eta()) , (dphi_23/dphi_12), pt );
                }
             }
          }
       }
     }
  }  


  //  Look at different Bxdl btw. stations & rings
  /// All possible segment pair in CSC 
  ///                 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 
  int csc1[2][17]= {{11,11,12,12,13,11,11,12,12,13,11,12,21,21,22,21,31},
                    {21,22,21,22,22,31,32,31,32,32,41,41,31,32,32,41,41}}; 
  for (int l=0; l<17; l++) {
      int s1 = csc1[0][l]/10;
      int r1 = csc1[0][l]%10;
      int s2 = csc1[1][l]/10;
      int r2 = csc1[1][l]%10;
      if ( path[s1][r1] && path[s2][r2] && path_r[s1][r1] && path_r[s2][r2] ) {
         double ME_Resol = d_ephi[s1][s2]-d_ephi_r[s1][s2];
         histo4 = hME1[l];
         histo4->Fill8( (pt1[s1]*d_ephi[s1][s2]) ,d_ephi[s1][s2] ,d_eeta[s1][s2], fabs(eeta[s2]),pt1[s1] );
         histo4->Fill8a( (pt1[s1]*d_ephi_r[s1][s2]) ,d_ephi_r[s1][s2] ,d_eeta_r[s1][s2], fabs(eeta_r[s2]),pt1[s1],ME_Resol );
      }
  }
  
  //  Look at different Bxdl btw. stations & rings
  /// All possible segment pair in DT 
  ///               0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
  int dt1[2][30]={{10,10,10,10,10,10,11,11,11,11,11,11,12,12,12,20,20,20,20,21,21,21,21,22,22,30,30,31,31,32}
                 ,{20,21,30,31,40,41,21,22,31,32,41,42,22,32,42,30,31,40,41,31,32,41,42,32,42,40,41,41,42,42}};
  for (int l=0; l<30; l++) {
      int s1 = dt1[0][l]/10;
      int w1 = dt1[0][l]%10;
      int s2 = dt1[1][l]/10;
      int w2 = dt1[1][l]%10;
       
      if ( bpath[s1][w1] && bpath[s2][w2] && bpath_r[s1][w1] && bpath_r[s2][w2] ) {

         double MB_Resol = d_bphi[s1][s2]-d_bphi_r[s1][s2];
         if (s2!=4 ){
            histo5 = hMB1[l];
            histo5->Fill9( (pt1[s1]*d_bphi[s1][s2]), d_bphi[s1][s2], d_beta[s1][s2], fabs(beta[s2]), pt1[s1] );
            histo5->Fill9a( (pt1[s1]*d_bphi_r[s1][s2]), d_bphi_r[s1][s2], d_beta_r[s1][s2], fabs(beta_r[s2]), pt1[s1], MB_Resol);
         }
         if (s2==4 ){
            histo5 = hMB1[l];
            histo5->Fill9( (pt1[s1]*d_bphi[s1][s2]), d_bphi[s1][s2], d_beta[s1][s2], fabs(beta[s1]), pt1[s1] );
            histo5->Fill9a( (pt1[s1]*d_bphi_r[s1][s2]), d_bphi_r[s1][s2], d_beta_r[s1][s2], fabs(beta_r[s1]), pt1[s1], MB_Resol);
         }
      }
  }

  //  Look at different Bxdl in overlap region
  /// All possible segment pairs in overlap region 
  ///               0  1  2  3  4  5  6  7  8  9 10  
  int olp[2][11]={{11,12,12,12,12,21,22,22,22,32,32} 
                 ,{13,12,13,22,32,13,13,22,32,13,22}};
  for (int l=0; l<11; l++) {
      int s1 = olp[0][l]/10;
      int w1 = olp[0][l]%10;
      int s2 = olp[1][l]/10;
      int w2 = olp[1][l]%10; 
      if ( bpath[s1][w1] && path[s2][w2] && bpath_r[s1][w1] && path_r[s2][w2] && overlap ) {

         double OL_Resol = d_ophi[s1][s2]-d_ophi_r[s1][s2];
         histo10 = hOL1[l];
         histo10->Fill12( (pt1[s1]*d_ophi[s1][s2]), d_ophi[s1][s2], d_oeta[s1][s2], fabs(eeta[s2]), pt1[s1] );
         histo10->Fill12a( (pt1[s1]*d_ophi_r[s1][s2]), d_ophi_r[s1][s2], d_oeta_r[s1][s2], fabs(eeta_r[s2]), pt1[s1], OL_Resol);
      }
  }
  
  // Fill  the 1 segment case for CSC and DT
  int csc2[8] ={11,12,13,21,22,31,32,41};
  double ME_phiRT[4][5] = {{999.0},{999.0}};
  double ME_phiRT_r[4][5] = {{999.0},{999.0}};
  for (int l=0; l<8; l++) {
      int s1 = csc2[l]/10;
      int r1 = csc2[l]%10;
      if ( path[s1][r1] && path_r[s1][r1] ) {
         double dME_phi = ME_phi_r[s1][r1] - ME_phi[s1][r1];
         double dME_eta = ME_eta_r[s1][r1] - ME_eta[s1][r1];
         histo6 = hME2[l];
         histo6->Fill8b( (pt*ME_phi[s1][r1]) , ME_phi[s1][r1], ME_eta[s1][r1], pt);
         histo6->Fill8c( (pt*ME_phi_r[s1][r1]) ,ME_phi_r[s1][r1], dME_phi , dME_eta, ME_eta_r[s1][r1] ,pt);
      
         for (int h=0; h<8; h++) {
              int s2 = csc2[h]/10;
              int r2 = csc2[h]%10;
              if ( (s2 > s1) && (path[s2][r2]) && (path_r[s2][r2]) ) {
                 ME_phiRT[s1][s2] = ME_phi[s2][r2] / ME_phi[s1][r1] ;
                 ME_phiRT_r[s1][s2] = ME_phi_r[s2][r2] / ME_phi_r[s1][r1] ;
                 if ( (s1==1) && (s2==2) ){ histo8 = hME3[0]; }
                 if ( (s1==1) && (s2==3) ){ histo8 = hME3[1]; }
                 if ( (s1==1) && (s2==4) ){ histo8 = hME3[2]; }
                 if ( (s1==2) && (s2==3) ){ histo8 = hME3[3]; }
                 if ( (s1==2) && (s2==4) ){ histo8 = hME3[4]; }
                 if ( (s1==3) && (s2==4) ){ histo8 = hME3[5]; }
                 histo8->Fill10(ME_phiRT[s1][s2],ME_phiRT_r[s1][s2],ME_eta[s1][r1],ME_eta_r[s1][r1],pt);
              }
         }
      }
  }


  int dt2[12] ={10,11,12,20,21,22,30,31,32,40,41,42};
  double MB_phiRT[4][5] = {{999.0},{999.0}};
  double MB_phiRT_r[4][5] = {{999.0},{999.0}};
  for (int l=0; l<12; l++) {
      int s1 = dt2[l]/10;
      int w1 = dt2[l]%10;
      if ( bpath[s1][w1] && bpath_r[s1][w1] ) {
         double dMB_phi = MB_phi_r[s1][w1] - MB_phi[s1][w1];
         double dMB_eta = MB_eta_r[s1][w1] - MB_eta[s1][w1];
         histo7 = hMB2[l];
         histo7->Fill9b( (pt*MB_phi[s1][w1]) , MB_phi[s1][w1], MB_eta[s1][w1], pt);
         histo7->Fill9c( (pt*MB_phi_r[s1][w1]) , MB_phi_r[s1][w1], dMB_phi , dMB_eta, MB_eta_r[s1][w1], pt);

         for (int h=0; h<8; h++) {
              int s2 = dt2[h]/10;
              int w2 = dt2[h]%10;
              if ( (s2 > s1) && (bpath[s2][w2]) && (bpath_r[s2][w2]) ) {
                 MB_phiRT[s1][s2] = MB_phi[s2][w2] / MB_phi[s1][w1] ;
                 MB_phiRT_r[s1][s2] = MB_phi_r[s2][w2] / MB_phi_r[s1][w1] ;
                 if ( (s1==1) && (s2==2) ){ histo9 = hMB3[0]; }
                 if ( (s1==1) && (s2==3) ){ histo9 = hMB3[1]; }
                 if ( (s1==1) && (s2==4) ){ histo9 = hMB3[2]; }
                 if ( (s1==2) && (s2==3) ){ histo9 = hMB3[3]; }
                 if ( (s1==2) && (s2==4) ){ histo9 = hMB3[4]; }
                 if ( (s1==3) && (s2==4) ){ histo9 = hMB3[5]; }
                 histo9->Fill11(MB_phiRT[s1][s2],MB_phiRT_r[s1][s2],MB_eta[s1][w1],MB_eta_r[s1][w1],pt);
              }
         }

      }
  }
  // Charge Assignment Testing !
  int q_csc[6] ={12,13,14,23,24,34};
  int eqv[6] = {0};
  int eq[6]  = {0};
  int q_edphi = 0;
  int eqv_r[6] = {0};
  int eq_r[6]  = {0};
  int q_edphi_r = 0;
  int bqv[6] = {0};
  int bq[6]  = {0};
  int q_bdphi = 0;
  int bqv_r[6] = {0};
  int bq_r[6]  = {0};
  int q_bdphi_r = 0;
  for (int l=0; l<6; l++) {
      int s1 = q_csc[l]/10;
      int s2 = q_csc[l]%10;
      if ((ephi[s1]!=999. )||(ephi[s2]!=999. )) {
         q_edphi = ((ephi[s1] - ephi[s2])> 0.0 ? 1:-1) ;
         q_edphi_r = ((ephi_r[s1] - ephi_r[s2])> 0.0 ? 1:-1) ;
      }
      if ((bphi[s1]!=999. )||(bphi[s2]!=999. )) {
         q_bdphi = ((bphi[s1] - bphi[s2])> 0.0 ? 1:-1) ;
         q_bdphi_r = ((bphi_r[s1] - bphi_r[s2])> 0.0 ? 1:-1) ;
      }
      eqv[l] = (l+1)*ChargeAssignment(egv[s1],egv[s2]);
      eq[l]  = (l+1)*q_edphi;
      eqv_r[l] = (l+1)*ChargeAssignment(egv_r[s1],egv_r[s2]);
      eq_r[l]  = (l+1)*q_edphi_r;
      bqv[l] = (l+1)*ChargeAssignment(bgv[s1],bgv[s2]);
      bq[l]  = (l+1)*q_bdphi;
      bqv_r[l] = (l+1)*ChargeAssignment(bgv_r[s1],bgv_r[s2]);
      bq_r[l]  = (l+1)*q_bdphi_r;
      if (  (eqv[l]!=0)&&(eq[l]!=0)  ) {
         histo1 = h_all;
         histo1->Fill2(eqv[l],eq[l]);
      }
      if (  (eqv_r[l]!=0)&&(eq_r[l]!=0)  ) {
         histo1 = h_all;
         histo1->Fill2a(eqv_r[l],eq_r[l]);
      }
      if (  (bqv[l]!=0)&&(bq[l]!=0)  ) {
         histo1 = h_all;
         histo1->Fill2b(bqv[l],bq[l]);
      }
      if (  (bqv_r[l]!=0)&&(bq_r[l]!=0)  ) {
         histo1 = h_all;
         histo1->Fill2c(bqv_r[l],bq_r[l]);
      }
  }
  
}

// ********************************************
// ***********  Utility functions  ************
// ********************************************

double pT_estimation(double p0, double p1, double the_eta, double the_phi) {
       double pT_est = (p0*the_eta + p1)/the_phi ;
       return pT_est; 
}

// number of csc segments in one chamber for each station
// cscseg_stat[0] = total segments in all stations
// cscseg_stat[5] = the number of stations which have segments
void SeedGeneratorTest::CSCsegment_stat( Handle<CSCSegmentCollection> cscSeg ) {

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
// number of dt segments in one chamber for each station
void SeedGeneratorTest::DTsegment_stat( Handle<DTRecSegment4DCollection> dtSeg ) {

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

// find the simHits which is corresponding to the segment
bool SeedGeneratorTest::SameChamber(CSCDetId SimDetId, CSCDetId SegDetId){
  
     if ( SimDetId.endcap()== SegDetId.endcap() && SimDetId.station()== SegDetId.station() &&
          SimDetId.ring()  == SegDetId.ring()   && SimDetId.chamber()== SegDetId.chamber() ){
          return true;
     }
     else {
          return false;
     }

}


void SeedGeneratorTest::SeedFromCSCRecHit(Handle<CSCRecHit2DCollection> cscrechit, ESHandle<CSCGeometry> cscGeom){
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

void SeedGeneratorTest::SeedFromDTRecHit(Handle<DTRecHitCollection> dtrechit, ESHandle<DTGeometry> dtGeom){

     //double phi[4]={999.0};
     for (int i=0; i <6; i++) {
         dtrh_sum[i]=0;
     }

     /*int id0[2]={0};
     int id1[2]={0};
     std::vector<DTRecHit1DPair> sp_f;
     std::vector<DTRecHit1DPair> sp_h;
     for (DTRecHitCollection::const_iterator r_it = dtrechit->begin(); r_it != dtrechit->end(); r_it++){
          DTWireId det_id = (*r_it).wireId();
          id0[0] = (det_id.wheel()/fabs(det_id.wheel()))*
                (det_id.station()*10000 + fabs(det_id.wheel())*1000 + det_id.chamber()) ;
          id0[1] = det_id.superLayer();
          if (id1[0]==0){ 
              id1[0] = id0[0] ; 
              id1[1] = id0[1] ; 
          }
          else if ( (id1[0]=id0[0])&&((id1[1]==1)||(id1[1]==3)) ){
                  
          }
          
     }*/

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
void SeedGeneratorTest::Eta_Test(Handle<edm::SimTrackContainer> simTracks,
                                 Handle<edm::PSimHitContainer> dsimHits,
                                 ESHandle<DTGeometry> dtGeom) {

  h_trk = -9.0;
  for (int j=0; j<5; j++) {
      h_sl1[j] = -9.0;
      h_sl2[j] = -9.0;
      h_sl3[j] = -9.0;
  }
  for (SimTrackContainer::const_iterator simTk_It = simTracks->begin(); simTk_It != simTracks->end(); simTk_It++)
  {
      if (abs((*simTk_It).type())!=13) continue;
      float px = ((*simTk_It).momentum()).x();
      float py = ((*simTk_It).momentum()).y();
      float pz = ((*simTk_It).momentum()).z();
      float pa = sqrt((px*px) + (py*py) + (pz*pz));
      double theta = acos( pz/pa );
      h_trk = fabs((-1.0)*log( tan(theta/2.0) ));

      for (PSimHitContainer::const_iterator ds_It = dsimHits->begin(); ds_It != dsimHits->end(); ds_It++)
      {          
          if ( ( abs((*ds_It).particleType())!=13 ) || ( (*ds_It).trackId()!=1 )) continue; 
          Local3DPoint lp = (*ds_It).localPosition(); 

          DTLayerId D_Id = DTLayerId( (*ds_It).detUnitId() );
          const DTLayer* dtlayer = dtGeom->layer(D_Id);
          //const DTChamber* dtchamber = dtGeom->chamber( D_Id );

          GlobalPoint gp = dtlayer->toGlobal(lp );
          //GlobalPoint gp = dtchamber->toGlobal(lp );

          if (D_Id.superLayer()==1) {
             h_sl1[D_Id.station()] = fabs(gp.eta());
          }
          if (D_Id.superLayer()==2) {
             h_sl2[D_Id.station()] = fabs(gp.eta());
          }
          if (D_Id.superLayer()==3) {
             h_sl3[D_Id.station()] = fabs(gp.eta());
          }
      }

  }

}

int SeedGeneratorTest::ChargeAssignment(GlobalVector Va, GlobalVector Vb){
     int charge = 0;
     float axb = ( Va.x()*Vb.y() ) - ( Vb.x()*Va.y() );
     if (axb != 0.0) {
        charge = ( (axb > 0.0) ?  1:-1 ) ;
     }
     return charge;
}

void SeedGeneratorTest::RecSeedReader( Handle<TrajectorySeedCollection> rec_seeds, ESHandle<DTGeometry> dtGeom, ESHandle<CSCGeometry> cscGeom){

     TrajectorySeedCollection::const_iterator seed_it;
     for (seed_it = rec_seeds->begin(); seed_it !=  rec_seeds->end(); seed_it++) {
         PTrajectoryStateOnDet pTSOD = (*seed_it).startingState();

         LocalTrajectoryParameters seed_para = pTSOD.parameters();
         /*
         double seed_mx = (seed_para.momentum()).x();
         double seed_my = (seed_para.momentum()).y();
         double seed_mT = sqrt((seed_mx*seed_mx)+(seed_my*seed_my));
         cout <<" (1) seed_Q= "<< seed_para.charge() <<endl;
         cout <<" (2) seed_m= "<< seed_para.momentum() <<endl;
         cout <<" (3) seed_p= "<< seed_para.position() <<endl;
         cout <<" (4) seed_mT= "<< seed_mT << endl;
         */ 

         nSegInSeed = (*seed_it).nHits();
         //cout <<"n hits = "<<(*seed_it).nHits() << endl;

         if ( (pTSOD.detId()/10000000) == 57) {
            DTChamberId MB_Id = DTChamberId( pTSOD.detId() );
            const DTChamber* dtchamber = dtGeom->chamber( MB_Id );
            seed_gp = dtchamber->toGlobal( seed_para.position() );
            //cout <<" *** DetID = "<< MB_Id <<" ***" <<endl;
            //cout <<" eta = "<<seed_gp.eta() <<endl; 
            
         }
         if ( (pTSOD.detId()/10000000) == 60) {
            CSCDetId  ME_Id = CSCDetId( pTSOD.detId() );
            const CSCChamber* cscchamber = cscGeom->chamber( ME_Id );
            seed_gp = cscchamber->toGlobal( seed_para.position() );
            //cout <<" *** DetID = "<< ME_Id <<" ***" <<endl;
            //cout <<" eta = "<<seed_gp.eta() <<endl; 
         }

     }
}

