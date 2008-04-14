// Class Header
#include "SeedValidator.h"
#include "RecoMuon/MuonSeedGenerator/test/MuonSeedPTAnalysis/SegSelector.h"

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
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include "TFile.h"
#include "TVector3.h"

#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <string>
#include <stdio.h>
#include <algorithm>

DEFINE_FWK_MODULE(SeedValidator);
using namespace std;
using namespace edm;
using namespace reco;


// constructors
SeedValidator::SeedValidator(const ParameterSet& pset){ 

  debug             = pset.getUntrackedParameter<bool>("debug");
  rootFileName      = pset.getUntrackedParameter<string>("rootFileName");
  recHitLabel       = pset.getUntrackedParameter<string>("recHitLabel");
  cscSegmentLabel   = pset.getUntrackedParameter<string>("cscSegmentLabel");
  dtrecHitLabel     = pset.getUntrackedParameter<string>("dtrecHitLabel");
  dtSegmentLabel    = pset.getUntrackedParameter<string>("dtSegmentLabel");
  simHitLabel       = pset.getUntrackedParameter<string>("simHitLabel");
  simTrackLabel     = pset.getUntrackedParameter<string>("simTrackLabel");
  muonseedLabel     = pset.getUntrackedParameter<string>("muonseedLabel");
  staTrackLabel     = pset.getParameter<InputTag>("staTrackLabel");
  glbTrackLabel     = pset.getParameter<InputTag>("glbTrackLabel");
  expectedPT        = pset.getUntrackedParameter<double>("expectedPT");
  scope             = pset.getUntrackedParameter<bool>("scope");
  pTCutMax          = pset.getUntrackedParameter<double>("pTCutMax");
  pTCutMin          = pset.getUntrackedParameter<double>("pTCutMin");
  eta_Low           = pset.getUntrackedParameter<double>("eta_Low");
  eta_High          = pset.getUntrackedParameter<double>("eta_High");

  ParameterSet serviceParameters = pset.getParameter<ParameterSet>("ServiceParameters");
  theService        = new MuonServiceProxy(serviceParameters);

  recsegSelector    = new SegSelector(pset);

  // Create the root file
  theFile = new TFile(rootFileName.c_str(), "RECREATE");
  theFile->mkdir("AllMuonSys");
  theFile->cd();
  theFile->mkdir("No_Seed");
  theFile->cd();
  theFile->mkdir("No_STA");
  theFile->cd();
  // TTree test
  tr_muon = new TNtuple1();

  h_all  = new H2DRecHit1("AllMu_");
  h_NoSeed  = new H2DRecHit2("NoSeed");
  h_NoSta   = new H2DRecHit3("NoSta");

}

// destructor
SeedValidator::~SeedValidator(){

  if (debug) cout << "[Seed Validation] Destructor called" << endl;
  delete recsegSelector;
  //delete muonSeedBuilder_; 
  // Write the histos to file
  theFile->cd();
  theFile->cd("AllMuonSys");
  h_all->Write();
  theFile->cd();
  theFile->cd("No_Seed");
  h_NoSeed->Write();
  theFile->cd();
  theFile->cd("No_STA");
  h_NoSta->Write();
  theFile->cd();

  // for tree
  theFile->cd();
  tr_muon->Write();

  // Release the memory ...
  delete h_all;
  delete h_NoSeed;
  delete h_NoSta;
  delete tr_muon;

  theFile->Close();
  if (debug) cout << "************* Finished writing histograms to file" << endl;
  if (theService) delete theService;

}

// The Main...Aanlysis...

void SeedValidator::analyze(const Event& event, const EventSetup& eventSetup)
{
  //Get the CSC Geometry :
  theService->update(eventSetup);

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
  // Get sta muon tracks
  Handle<TrackCollection> standAloneMuons;
  event.getByLabel(staTrackLabel, standAloneMuons);
  // Get global muon tracks
  Handle<TrackCollection> globalMuons;
  event.getByLabel(glbTrackLabel, globalMuons);

  // Magnetic field
  ESHandle<MagneticField> field;
  eventSetup.get<IdealMagneticFieldRecord>().get(field);

  H2DRecHit1 *histo1 = 0;   
  H2DRecHit2 *histo2 = 0;   
  H2DRecHit3 *histo3 = 0;   
  //  TNtuple1 *tt = 0;
 
  // Get sim track information  
  // return  eta_c, eta_d, phi_c, phi_d :  average eta and phi of simhits from CSC and DT  
  //         theQ, eta_trk, phi_trk, pt_trk, trackId : Q,eta,phi,pt and trackID of simtrack
  //         pt_layer,palayer : the pt or pa in each layer
  SimInfo(simTracks,dsimHits,csimHits,dtGeom,cscGeom);

  /// A. statistic information for segment, seed and sta muon w.r.t eta
  // 2. Reading the damn seeds
  RecSeedReader(muonseeds);

  //if (nu_seed > 9) nu_seed = 9;
  // 3. Get muon track information ; sta=0, glb=1
  StaTrackReader(standAloneMuons, 0);
  
  //// associate seeds with sim tracks 
  //// seed_trk => the simtrack number j associate with seed i 
  std::vector<int> seeds_trk(seed_gp.size(), -1);
  cout <<""<<endl; 
  cout<<" --------------------------------------------- "<<endl;
  cout<<" # of seed: "<<nu_seed<<"  # of sim: "<<eta_trk.size()<< endl;
  for (unsigned int i=0; i < seed_gp.size(); i++) {
    double dR1= 99.0 ;
    int preferTrack = -1;
    cout<<"  seed eta = "<<seed_gp[i].eta()<<"  seed phi = "<<seed_gp[i].phi()<<endl; 
    for (unsigned int j=0; j < eta_trk.size(); j++) {
      cout<<"  trj eta = "<<eta_trk[j]<<"  trj phi = "<<phi_trk[j]<<endl; 
      double dh = fabs(seed_gp[i].theta() - theta_trk[j]) ; 
      double df = fabs(seed_gp[i].phi() - phi_trk[j]) ;
      double dR2 = sqrt( dh*dh + df*df ) ;
      if ( dR2 < dR1  && dh < 0.35 && df < 0.5 ) {
         preferTrack = static_cast<int>(j) ;
         dR1 = dR2;
      }
    }
    seeds_trk[i] = preferTrack ;
    cout<<" seed("<<i<<") <-->  trk("<<preferTrack<<") / " <<eta_trk.size() <<endl;
  }

  //// associate sta with sim tracks  
  //// sta_simtrk => the simtrack number j associate with sta  i 
  std::vector<int> sta_simtrk(sta_theta.size(), -1);
  for (unsigned int i=0; i < sta_theta.size(); i++) {
    double dR1= 99.0 ;
    int preferTrack = -1;
    for (unsigned int j=0; j < eta_trk.size(); j++) {
      double dh = fabs(sta_theta[i] - theta_trk[j]) ; 
      double df = fabs(sta_phi[i] - phi_trk[j]) ;
      double dR2 = sqrt( dh*dh + df*df ) ;
      if ( dR2 < dR1  && dh < 0.35 && df < 0.5 ) {
         preferTrack = static_cast<int>(j) ;
         dR1 = dR2;
      }
    }
    sta_simtrk[i] = preferTrack ;
  }


  // look at the number of reconstructed seed/sta for each sim-track
  histo1 = h_all;
  int sta_Evt  =0;
  for (unsigned int i=0; i < eta_trk.size(); i++) {

      int nu_seeds_trk =0;
      for(unsigned int j=0; j < seeds_trk.size(); j++){
         if ( seeds_trk[j]== static_cast<int>(i) ) nu_seeds_trk++;
      }

      int nu_sta_trk =0;
      for(unsigned int j=0; j < sta_simtrk.size(); j++){
         if ( sta_simtrk[j]== static_cast<int>(i) ) nu_sta_trk++;
      }
      std::vector<double> pa_tmp = palayer[i] ; 
      histo1->Fill1b( eta_trk[i], nu_seeds_trk, nu_sta_trk, pt_trk[i],
                      theQ[i]/pt_trk[i], theQ[i]/pa_tmp[0] );

      if (nu_sta_trk > 0) {  
         sta_Evt++ ;
      } else {
         sta_Evt = 0 ;
      }
  }
 
  int reco_sim = sta_Evt - eta_trk.size() ;
  histo1->Fill1o(sta_mT.size(), seed_mT.size(), reco_sim );

  // look at those un-associated seeds and sta
  for (unsigned int i=0; i < seeds_trk.size(); i++ ) {
      if ( seeds_trk[i]== -1 )   histo1->Fill1d1( seed_gp[i].eta() );
  }
  for (unsigned int i=0; i < sta_simtrk.size(); i++ ) {
      if ( sta_simtrk[i]== -1 )  histo1->Fill1d2( getEta(sta_theta[i]) );
  }

  /// B. seeds  Q/PT pulls
  if (nu_seed > 0) {
    int bestSeed = -1;
    double dSeedPT = 99999.9 ;
    // look the seed for every sta tracks
    int trkId = 0;
    for (unsigned int i=0; i < eta_trk.size(); i++) {

        // find the best seed whose pt is closest to a simTrack which associate with 
        for(unsigned int j=0; j < seeds_trk.size(); j++){
           if ( ( seeds_trk[j] == trkId ) && ( fabs(seed_mT[j] - pt_trk[trkId]) < dSeedPT ) ){
              dSeedPT = fabs( seed_mT[j] - pt_trk[i] );
              bestSeed = static_cast<int>(j);
           }
        }

        for(unsigned int j=0; j < seeds_trk.size(); j++){

           if (  seeds_trk[j]!= trkId ) continue;

           // put the rest seeds in the minus side 
           double bestSeedPt = 99999.9;
           if (bestSeed == static_cast<int>(j) ) {
              bestSeedPt = seed_mT[j];
           }else {
              bestSeedPt = -1.0*seed_mT[j];
           }

           std::vector<double> pa1 = palayer[i];
           std::vector<double> pt1 = ptlayer[i];
	   double pull_qbp  = ( qbp[j]  - (theQ[j]/pa1[ seed_layer[j] ]) ) / err_qbp[j] ;
	   double pull_qbpt = ( qbpt[j] - (theQ[j]/pt1[ seed_layer[j] ]) ) / err_qbpt[j] ;
	   histo1->Fill1g( seed_mT[j], seed_mA[j], bestSeedPt, seed_gp[j].eta() ) ;
	   histo1->Fill1i( pull_qbp, seed_gp[j].eta(), qbpt[j], pull_qbpt, err_qbp[j], err_qbpt[j] );
	   histo1->Fill1f( seed_gp[j].eta(), err_dx[j], err_dy[j], err_x[j], err_y[j]);

        }
        trkId++;
    }

  }

  /// C. sta track information
  if (nu_sta > 0 ) {
    double dPT = 9999999.9 ;
    int best = 0;
    for (unsigned int i=0; i < eta_trk.size(); i++) {
        // find the best sta whose pt is closest to a simTrack which associate with 
        for(unsigned int j=0; j < sta_simtrk.size(); j++){
           if ( ( sta_simtrk[j]==static_cast<int>(i) ) && ( fabs(sta_mT[j] - pt_trk[i]) < dPT ) ){
              dPT = fabs( seed_mT[j] - pt_trk[i] );
              best = j;
           }
        }
       // looking for the sta muon which is closest to simTrack pt 
       double dh_trk = fabs(sta_theta[best] - theta_trk[i]);
       double df_trk = fabs(sta_phi[best] - phi_trk[i]);
       if (dh_trk < 0.35 && df_trk < 0.35 ) {
          histo1->Fill1j(eta_trk[i], sta_qbp[best], sta_qbpt[best], sta_mT[best], sta_mA[best]);
       }
    }


    // look at the dPhi, dEta, dx, dy from sim-segment and reco-segment of a seed
    for (unsigned int i=0; i<nSegInSeed.size(); i++) {

        int ii = seeds_trk[i]  ;
        std::vector<SimSegment> sCSC_v = recsegSelector->Sim_CSCSegments( trackID[ii], csimHits, cscGeom);
        std::vector<SimSegment> sDT_v  = recsegSelector->Sim_DTSegments( trackID[ii], dsimHits, dtGeom);

        SegOfRecSeed(muonseeds,i,sCSC_v,sDT_v); 
        for ( unsigned int j=0; j<geoID.size(); j++  ) {
            if (geoID[j].subdetId() ==  MuonSubdetId::DT ) {
               DTChamberId MB_Id = DTChamberId( geoID[j]  );
               histo1->Fill1e(-1*MB_Id.station() ,d_h[j],d_f[j],d_x[j],d_y[j]); 
            }
            if (geoID[j].subdetId() ==  MuonSubdetId::CSC ) {
               CSCDetId  ME_Id = CSCDetId( geoID[j] );
               histo1->Fill1e(ME_Id.station(),d_h[j],d_f[j],d_x[j],d_y[j]); 
            }
        }
    }
    // look at simhits information !!! 
    /*
    if (debug) {
       cout <<"   "<<endl;
       cout <<" ^^^^^^^^^^^^^^^^ Sim information ^^^^^^^^^^^^^^^^^ " <<endl;
       for (std::vector<SimSegment>::const_iterator i1 = sDT_v.begin(); i1!=sDT_v.end(); i1++ ) {
	    cout <<"DId: "<<(*i1).dt_DetId <<endl;
	    cout<<"  h0= "<< (*i1).sGlobalOrg.eta() <<"  f0= "<<(*i1).sGlobalOrg.phi();
	    cout<<"  gp= "<< (*i1).sGlobalOrg <<endl;
	    cout<<"  h1= "<< (*i1).sGlobalVec.eta() <<"  f1= "<<(*i1).sGlobalVec.phi();
	    cout<<"  gv= "<< (*i1).sGlobalVec <<endl;
            cout<<" -------------------------------------------- "<<endl; 
       }
       for (std::vector<SimSegment>::const_iterator i1 = sCSC_v.begin(); i1!=sCSC_v.end(); i1++ ) {
	    cout <<"DId: "<<(*i1).csc_DetId <<endl;
	    cout<<"  h0= "<< (*i1).sGlobalOrg.eta() <<"  f0= "<<(*i1).sGlobalOrg.phi();
	    cout<<"  gp= "<< (*i1).sGlobalOrg <<endl;
	    cout<<"  h1= "<< (*i1).sGlobalVec.eta() <<"  f1= "<<(*i1).sGlobalVec.phi();
	    cout<<"  gv= "<< (*i1).sGlobalVec <<endl;
            cout<<" -------------------------------------------- "<<endl; 
       }
       cout <<" "<<endl;
       cout <<"      an event finished !!!    "<<endl;
       cout <<" "<<endl;
    }
    */

    /// open a scope to check all information
     if ( fabs( getEta( sta_theta[best]) ) > eta_Low  && fabs(getEta( sta_theta[best]) ) < eta_High 
          &&    sta_mT[best] > pTCutMin &&        sta_mT[best] < pTCutMax   ){

        if ( scope ) { cout <<" ************ pT-eta scope *************** " <<endl; }

        // check the seed and sta phi dirstribition
        if (getEta(sta_theta[best]) < 0 ) {     
           histo1->Fill1l0( seed_gp[best].phi(), sta_phi[best] );
        }else {
           histo1->Fill1l1( seed_gp[best].phi(), sta_phi[best] );
        }

        // look at the dPhi, dEta, dx, dy from sim-segment and reco-segment of a seed
	for (unsigned int i=0; i<nSegInSeed.size(); i++) {

            //if ( nSegInSeed.size() == 1 ) continue ;
            if (scope) {
               cout <<" "<<i<<" seed w/ "<<nSegInSeed[i]<<" segs";
	       cout <<" & pt= "<<seed_mT[i]<<" +/- "<<err_qbp[i]*seed_mT[i]*seed_mA[i];
	       cout <<" q/p= "<<qbp[i]<<" @ "<<seed_layer[i];
	       cout <<" dx= "<<err_dx[i]<<" dy= "<<err_dy[i]<<" x= "<<err_x[i]<<" y= "<<err_y[i]<<endl;
            }
    
            // check the segments of the seeds in the scope
            bool flip_debug = false;
            if (!debug)  {
               debug = true;
               flip_debug = true;
            }

            int ii = seeds_trk[i]  ;
            std::vector<SimSegment> sCSC_v = recsegSelector->Sim_CSCSegments( trackID[ii], csimHits, cscGeom);
            std::vector<SimSegment> sDT_v  = recsegSelector->Sim_DTSegments( trackID[ii], dsimHits, dtGeom);
	    SegOfRecSeed(muonseeds,i,sCSC_v,sDT_v);
            if (flip_debug) {
               debug = false;
            }

	    for ( unsigned int j=0; j<geoID.size(); j++  ) {
		if (geoID[j].subdetId() ==  MuonSubdetId::DT ) {
		   DTChamberId MB_Id = DTChamberId( geoID[j]  );
		   histo1->Fill1k(-1*MB_Id.station() ,d_h[j],d_f[j],d_x[j],d_y[j]); 
		}
		if (geoID[j].subdetId() ==  MuonSubdetId::CSC ) {
	           CSCDetId  ME_Id = CSCDetId( geoID[j] );
		   histo1->Fill1k(ME_Id.station(),d_h[j],d_f[j],d_x[j],d_y[j]); 
		}
	   }
	}


        // look at sta information
        if (scope) {
           cout <<"  "<<endl;
	   for (unsigned int x=0; x<sta_mT.size(); x++) {
	       cout <<" sta_pt= "<<sta_mT[x]<<" q/p= "<<sta_qbp[x]<<" w/"<<sta_nHits[x] <<endl;
	   }
           cout <<"************************************************"<<endl;
           cout <<"  "<<endl;
        }


       // look at the sta pt vs. #_of_hits and chi2
       int elected = -1 ;
       int length  = -1 ;
       double X2   = -1. ;

       for (unsigned int x=0; x<sta_mT.size(); x++) {

           double ndf = static_cast<double>(2*sta_nHits[x]-4);
           //double ndf = 1.0;
           double chi2_ndf = sta_chi2[x]/ndf;
 
           histo1->Fill1m(sta_mT[x], sta_nHits[x], chi2_ndf);
           // pass the sta cleaning criteria
           
           if ( length < 0 ) {
              length = sta_nHits[x];
              X2= chi2_ndf;
              elected = x;
           } else if ( sta_nHits[x] -  length > 3){
           //} else if ( sta_nHits[x] -  length > 2){
              elected = x;
              length = sta_nHits[x];
              X2= chi2_ndf;
           } else if ( sta_nHits[x] - length > 0 && sta_nHits[x] - length < 3){
           //} else if ( sta_nHits[x] == length ){
              if ( chi2_ndf < X2 ) {
                 elected = x;
                 length = sta_nHits[x];
                 X2= chi2_ndf;
              }
           } else {
              continue;
           }
           
           /*
           if ( X2 < 0 ) {
              length = sta_nHits[x];
              X2= chi2_ndf;
              elected = x;
           } else if ( chi2_ndf < X2 ){
              elected = x;
              length = sta_nHits[x];
              X2= chi2_ndf;
           } else {
             continue;
           }
           */
       }
       if (elected >= 0)  histo1->Fill1n(sta_mT[elected], sta_nHits[elected], sta_chi2[elected]/(2.0*sta_nHits[elected]-4.0));

     }

  }

  // D. fail sta cases
  if (nu_sta == 0 && nu_seed!=0) {
     double sim_eta = -9.0 ; 
     for (int i=0; i < nu_seed; i++) {

         histo3 = h_NoSta;
         double pull_qbpt = ( qbpt[i] - (theQ[0]/pt_trk[0]) ) / err_qbpt[i] ;
         double pt_err = err_qbpt[i]*seed_mT[i]*seed_mT[i];
         if (seeds_trk[i] != -1 )  sim_eta = eta_trk[ seeds_trk[i] ] ; 
         histo3->Fill3a( sim_eta , seed_gp[i].eta(), nu_seed, seed_mT[i] , pt_err, pull_qbpt);
    
         CSCsegment_stat(cscSegments, cscGeom, seed_gp[i].eta(), seed_gp[i].phi() );
         DTsegment_stat(dt4DSegments, dtGeom, seed_gp[i].eta(), seed_gp[i].phi());

         int allseg1 = cscseg_stat1[5]+dtseg_stat1[5];
         int allseg  = cscseg_stat[0]+dtseg_stat[0];
         histo3->Fill3b(sim_eta,allseg1,allseg );

         int types = RecSegReader(cscSegments,dt4DSegments,cscGeom,dtGeom, seed_gp[i].eta(), seed_gp[i].phi());
         if (debug) cout<<" seed type for fail STA= "<<types<<endl;
         for (unsigned i=0; i< phi_resid.size(); i++) {
             histo3->Fill3c(phi_resid[i],eta_resid[i]);
         }
     }
     
  }
 
  cout <<" ================================================== "<<endl; 
  cout <<" # of seeds: "<<nu_seed<<" # of STA: "<<nu_sta<<endl;
  // Basic simulation and reco information
  int idx=0;
  for (SimTrackContainer::const_iterator stk = simTracks->begin(); stk != simTracks->end(); stk++)
  {
      //if( abs((*stk).type())==13 ) cout<<" vtxId : "<<(*stk).vertIndex()<<endl;
      //if (abs((*stk).type())!=13 || (*stk).vertIndex() != 0 ) continue;

      bool rechitSize = (dsimHits->size() < 8 && csimHits->size() < 4 )? true:false ;
      if (abs((*stk).type())!=13 || rechitSize || (*stk).vertIndex() != 0 ) continue;

      // 1. Run the class SegSelector
      int trkId = static_cast<int>( (*stk).trackId() );
      std::vector<SimSegment>     sCSC_v = recsegSelector->Sim_CSCSegments( trkId, csimHits, cscGeom);
      std::vector<SimSegment>     sDT_v = recsegSelector->Sim_DTSegments( trkId, dsimHits, dtGeom);
      std::vector<CSCSegment>     cscseg_V = recsegSelector->Select_CSCSeg(cscSegments,cscGeom, sCSC_v);
      std::vector<DTRecSegment4D> dtseg_V = recsegSelector->Select_DTSeg(dt4DSegments,dtGeom, sDT_v);

      cout<<" DT  Size: "<< dtseg_V.size()<<" CSC  Size: "<< cscseg_V.size()<<endl;
      cout<<" sDT Size: "<<  sDT_v.size() <<" sCSC Size: "<< sCSC_v.size()  <<endl;
      cout<<" rDT Size: "<<  dt4DSegments->size() <<" rCSC Size: "<< cscSegments->size()  <<endl;
      // 2. Reading the reco segmenst
      //   get the appropriate eta and phi to select reco-segments
      double px = ((*stk).momentum()).x();
      double py = ((*stk).momentum()).y();
      double pz = ((*stk).momentum()).z();
      double pt = sqrt(px*px + py*py);
      double trkEta = 0.0;
      double trkPhi = 0.0;
      double theta = acos( pz/sqrt(pz*pz + pt*pt));
       
      cout<<" <<<<< Idx : "<<idx<<" >>>>>>>>>>>>>>>>>>>>> "<<endl;
      cout<<"  eta from SimTrack= "<<getEta(px,py,pz)<<"  theta = "<<theta<<"("<<(180.0*theta)/3.1415 <<")"<<endl;
      if( eta_d[idx] != -99.0 && eta_c[idx] == -99.0 ) {
          trkEta = eta_d[idx] ;
          trkPhi = phi_d[idx] ; }
      if( eta_c[idx] != -99.0 && eta_d[idx] == -99.0) {
          trkEta = eta_c[idx] ;
          trkPhi = phi_c[idx] ; }
      if( eta_c[idx] != -99.0 && eta_d[idx] != -99.0 ) {
          trkEta = (eta_c[idx] + eta_d[idx])/2.0;
          trkPhi = (phi_c[idx] + phi_d[idx])/2.0;       }
      if( eta_c[idx] == -99.0 && eta_d[idx] == -99.0 ) {
          //cout<<" no simhit !"<<endl;
          //trkEta = getEta(px,py,pz) ; 
          trkEta = -99.0 ; 
      }
      double trkTheta = 2.*atan( exp(-1.*trkEta) ) ;
      cout<<"  eta in MuonSystem=  "<<trkEta<<"  theta= "<<trkTheta<<"("<<(180.0*trkTheta)/3.1415 <<")"<< endl; 
      cout<<"  phi from SimTrack= "<<atan2(py,px)<<endl;
      cout<<"  phi in MuonSystem= "<<trkPhi<<endl;
 
      // 3. Read out reco segments
      //   return : the ave_phi, ave_eta, phi_resid, eta_resid, dx_error, dy_error, x_error, y_error 
      int types = RecSegReader(cscSegments,dt4DSegments,cscGeom,dtGeom,trkTheta,trkPhi);

      // 4. Check # of segments and rechits in each chambers for this track
      CSCsegment_stat(cscSegments, cscGeom, trkTheta, trkPhi);
      DTsegment_stat(dt4DSegments, dtGeom, trkTheta, trkPhi);
      Simsegment_stat(sCSC_v,sDT_v);
      // 5. if less than 1 segment, check the # of rechits in each chambers for this track
      if ( (cscseg_stat[5] < 2)&&(dtseg_stat[5] < 2) ) {
         CSCRecHit_Stat(csc2DRecHits, cscGeom, trkEta, trkPhi);
         DTRecHit_Stat(dt1DRecHits, dtGeom, trkEta, trkPhi);
      }

      // seg_stat[0] = total segments in all stations
      // seg_stat[5] = the number of stations which have segments
      // seg_stat1[x] = the number of stations/segments which only count segments w/ more than 4 rechits
      int layer_sum  = cscseg_stat1[5] + dtseg_stat1[5];
      int seg_sum    = cscseg_stat[0]  + dtseg_stat[0];
      int leftrh_sum = cscrh_sum[5]    + dtrh_sum[5];
      cout <<"  sims:"<<simseg_sum<<" seg:"<<seg_sum<<" layer:"<<layer_sum<<endl;

      // look at those events without any reco segments
      histo1 = h_all;
      // look at all information
      histo1->Fill1(cscseg_stat[5], dtseg_stat[5], layer_sum, simseg_sum, getEta(px,py,pz) );
      // 6. pt vs. # of segments in a event
      histo1->Fill1c(pt, (*stk).momentum().mag(), cscseg_stat[0]+dtseg_stat[0]);
      if ((cscseg_stat[5]==0)&&(dtseg_stat[5]==0)) { 

	 histo1->Fill1a(leftrh_sum, trkEta);
         if (cscrh_sum[0] !=0 ) {
            histo2 = h_NoSeed;
            histo2 -> Fill2b(trkEta, cscrh_sum[0]);
         }
         if (dtrh_sum[0] !=0 ) {
            histo2 = h_NoSeed;
            histo2 -> Fill2c(trkEta, dtrh_sum[0]);
         }
      }

      if ( nu_seed == 0 ) {
         cout<<" seed type for no seed : "<<types<<" h: "<< getEta(px,py,pz) <<" w/ seg# "<<layer_sum<<endl; 
         cout<<" idx : "<<idx<<endl;
         histo2 = h_NoSeed;
         histo2->Fill2a( getEta(px,py,pz), layer_sum, seg_sum , simseg_sum, 
                        simseg_sum - seg_sum, simseg_sum - layer_sum );
      }
      idx++;
  }


}

// ********************************************
// ***********  Utility functions  ************
// ********************************************

// number of csc segments in one chamber for each station
// cscseg_stat[0] = total segments in all stations
// cscseg_stat[5] = the number of stations which have segments
void SeedValidator::CSCsegment_stat( Handle<CSCSegmentCollection> cscSeg , ESHandle<CSCGeometry> cscGeom, double trkTheta, double trkPhi) {

     for (int i=0; i<6; i++) {
         cscseg_stat[i]=0;
         cscseg_stat1[i]=0;
     }
     for(CSCSegmentCollection::const_iterator seg_It = cscSeg->begin(); seg_It != cscSeg->end(); seg_It++)
     { 
        CSCDetId DetId = (CSCDetId)(*seg_It).cscDetId();
	const CSCChamber* cscchamber = cscGeom->chamber( DetId );
	GlobalPoint gp = cscchamber->toGlobal((*seg_It).localPosition() );
        if (( fabs(gp.theta()- trkTheta) > 0.35  ) || ( fabs(gp.phi()- trkPhi) > 0.35 ) ) continue;

        cscseg_stat[DetId.station()] += 1;
        if ((*seg_It).nRecHits() > 3 ) {
           cscseg_stat1[DetId.station()] += 1;
        }
     }
     cscseg_stat[0] = cscseg_stat[1]+cscseg_stat[2]+cscseg_stat[3]+cscseg_stat[4];
     cscseg_stat1[0] = cscseg_stat1[1]+cscseg_stat1[2]+cscseg_stat1[3]+cscseg_stat1[4];
     for (int i =1; i<5; i++){
         if(cscseg_stat[i]!=0)  { cscseg_stat[5]++  ;}
         if(cscseg_stat1[i]!=0) { cscseg_stat1[5]++ ;}
     }
     
}
// number of dt segments in one chamber for each station
void SeedValidator::DTsegment_stat( Handle<DTRecSegment4DCollection> dtSeg, ESHandle<DTGeometry> dtGeom, double trkTheta, double trkPhi)  {

     for (int i=0; i<6; i++) {
         dtseg_stat[i]=0;
         dtseg_stat1[i]=0;
     }
     for(DTRecSegment4DCollection::const_iterator seg_It = dtSeg->begin(); seg_It != dtSeg->end(); seg_It++)
     { 
        DTChamberId DetId = (*seg_It).chamberId();
        const DTChamber* dtchamber = dtGeom->chamber( DetId );
        GlobalPoint  gp = dtchamber->toGlobal( (*seg_It).localPosition() );
        if ( ( fabs(gp.theta()- trkTheta) > 0.35  ) || ( fabs(gp.phi()- trkPhi) > 0.35 ) ) continue;

        dtseg_stat[DetId.station()] += 1;
        int n_phiHits = ((*seg_It).phiSegment())->specificRecHits().size();
        if ( (*seg_It).hasZed() && (n_phiHits > 4) ) {
           dtseg_stat1[DetId.station()] += 1;
        }
     }
     dtseg_stat[0]  = dtseg_stat[1] +dtseg_stat[2] +dtseg_stat[3] +dtseg_stat[4];
     dtseg_stat1[0] = dtseg_stat1[1]+dtseg_stat1[2]+dtseg_stat1[3]+dtseg_stat1[4];
     for (int i =1; i<5; i++){
         if(dtseg_stat[i]!=0)  { dtseg_stat[5]++ ;}
         if(dtseg_stat1[i]!=0) {
            if ( i !=4 ) { dtseg_stat1[5]++ ;}
            // because no eta/Z measurement at station 4
            if ((i==4)&&(dtseg_stat[4]!=0)) { dtseg_stat1[5]++ ;} 
         }
     }
}

// number of sim segments in one chamber for each station
void SeedValidator::Simsegment_stat( std::vector<SimSegment> sCSC, std::vector<SimSegment> sDT ) {

     for (int i=0; i<6; i++) {
         simcscseg[i] = 0;
         simdtseg[i]=0;
     }

     double ns1 =0.0;
     double eta_sim1 =0;
     for (std::vector<SimSegment>::const_iterator it = sCSC.begin(); it != sCSC.end(); it++) {
	     int st = ((*it).csc_DetId).station();
	     eta_sim1 += ((*it).sGlobalOrg).eta();
	     simcscseg[st]++;
	     ns1++;
     }
     simcscseg[0]=simcscseg[1]+simcscseg[2]+simcscseg[3]+simcscseg[4];
     for (int i=1; i<5; i++) {
	 if (simcscseg[i]!=0)  simcscseg[5]++; 
     }
     eta_sim1 = eta_sim1/ns1 ;

     double ns2 =0.0;
     double eta_sim2 =0;
     for (std::vector<SimSegment>::const_iterator it = sDT.begin(); it != sDT.end(); it++) {
	     int st = ((*it).dt_DetId).station();
	     eta_sim2 += ((*it).sGlobalOrg).eta();
	     simdtseg[st]++;
	     ns2++;
     }
     simdtseg[0]=simdtseg[1]+simdtseg[2]+simdtseg[3]+simdtseg[4];
     for (int i=1; i<5; i++) {
	 if (simdtseg[i]!=0) simdtseg[5]++; 
     }
     eta_sim2 = eta_sim2/ns2 ;

     simseg_sum = simcscseg[5]+ simdtseg[5];
     simseg_eta = -9.0;
     if      ((simcscseg[0]==0)&&(simdtseg[0]!=0)) { simseg_eta = eta_sim2; }
     else if ((simdtseg[0]==0)&&(simcscseg[0]!=0)) { simseg_eta = eta_sim1; }
     else { simseg_eta = (eta_sim1 + eta_sim2)/2.0 ; }
}

void SeedValidator::CSCRecHit_Stat(Handle<CSCRecHit2DCollection> cscrechit, ESHandle<CSCGeometry> cscGeom, double trkEta, double trkPhi){
     for (int i=0; i <6; i++) {
         cscrh_sum[i]=0;
     }
     for(CSCRecHit2DCollection::const_iterator r_it = cscrechit->begin(); r_it != cscrechit->end(); r_it++)
     { 
        CSCDetId det_id = (CSCDetId)(*r_it).cscDetId();
	const CSCChamber* cscchamber = cscGeom->chamber( det_id );
	GlobalPoint gp = cscchamber->toGlobal((*r_it).localPosition() );
        if (( fabs(gp.eta()- trkEta) > 0.2  ) || ( fabs(gp.phi()- trkPhi) > 0.5 ) ) continue;

        cscrh_sum[det_id.station()]++;
     }
     cscrh_sum[0] = cscrh_sum[1]+cscrh_sum[2]+cscrh_sum[3]+cscrh_sum[4];
     for (int i =1; i<5; i++){
         if(cscrh_sum[i]!=0) {
            cscrh_sum[5]++ ;
         }
     }
}

void SeedValidator::DTRecHit_Stat(Handle<DTRecHitCollection> dtrechit, ESHandle<DTGeometry> dtGeom, double trkEta, double trkPhi){

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
         if ( ( fabs(grh.eta()- trkEta) > 0.2  ) || ( fabs(grh.phi()- trkPhi) > 0.5 ) ) continue;

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

int SeedValidator::ChargeAssignment(GlobalVector Va, GlobalVector Vb){
     int charge = 0;
     float axb = ( Va.x()*Vb.y() ) - ( Vb.x()*Va.y() );
     if (axb != 0.0) {
        charge = ( (axb > 0.0) ?  1:-1 ) ;
     }
     return charge;
}

void SeedValidator::RecSeedReader( Handle<TrajectorySeedCollection> rec_seeds ){

     nu_seed = 0;
     seed_gp.clear();
     seed_gm.clear();
     seed_lp.clear();
     seed_lv.clear();
     qbp.clear();
     qbpt.clear();
     err_qbp.clear();
     err_qbpt.clear();
     err_dx.clear();
     err_dy.clear();
     err_x.clear();
     err_y.clear();
     seed_mT.clear();
     seed_mA.clear();
     seed_layer.clear();
     nSegInSeed.clear();

     TrajectorySeedCollection::const_iterator seed_it;
     for (seed_it = rec_seeds->begin(); seed_it !=  rec_seeds->end(); seed_it++) {
         PTrajectoryStateOnDet pTSOD = (*seed_it).startingState();

         nSegInSeed.push_back( (*seed_it).nHits() );
         // Get the tsos of the seed
         TrajectoryStateTransform tsTransform;
         DetId seedDetId(pTSOD.detId());
         const GeomDet* gdet = theService->trackingGeometry()->idToDet( seedDetId );
         TrajectoryStateOnSurface seedTSOS = tsTransform.transientState(pTSOD, &(gdet->surface()), &*theService->magneticField());
         // seed global position and momentum(direction)
         seed_gp.push_back( seedTSOS.globalPosition() );
         seed_gm.push_back( seedTSOS.globalMomentum() );
         seed_lp.push_back( seedTSOS.localPosition()  );
         seed_lv.push_back( seedTSOS.localDirection() );
        
         LocalTrajectoryParameters seed_para = pTSOD.parameters();

         // the error_vector v[15]-> [0] , [2] , [5] , [9] , [14]
         //                          q/p   dx    dy     x      y
         std::vector<float> err_mx = pTSOD.errorMatrix();
         err_qbp.push_back( sqrt(err_mx[0]) );
         err_dx.push_back( sqrt(err_mx[2])  );
         err_dy.push_back( sqrt(err_mx[5])  );
         err_x.push_back(  sqrt(err_mx[9])  );
         err_y.push_back(  sqrt(err_mx[14]) );
         //for (unsigned i=0; i< err_mx.size(); i++) {
         //    cout <<"Err"<<i<<" = "<<err_mx[i]<<"  -> "<<sqrt(err_mx[i])<<endl;
         //}

         // seed layer
         DetId pdid(pTSOD.detId()); 
         if ( pdid.subdetId()  == MuonSubdetId::DT ) {
            DTChamberId MB_Id = DTChamberId( pTSOD.detId() );
            seed_layer.push_back( MB_Id.station() );
         }
         if ( pdid.subdetId()  == MuonSubdetId::CSC ) {
            CSCDetId  ME_Id = CSCDetId( pTSOD.detId() );
            seed_layer.push_back( ME_Id.station() );
         }
         double seed_mx = seed_gm[nu_seed].x();
         double seed_my = seed_gm[nu_seed].y();
         double seed_mz = seed_gm[nu_seed].z();
         double seed_q  = 1.0*seed_para.charge();
         double seed_sin = sqrt((seed_mx*seed_mx)+(seed_my*seed_my)) / sqrt((seed_mx*seed_mx)+(seed_my*seed_my)+(seed_mz*seed_mz));
         // seed pt, pa, and q/pt or q/pa
         seed_mT.push_back( sqrt((seed_mx*seed_mx)+(seed_my*seed_my)) );
         seed_mA.push_back(  sqrt((seed_mx*seed_mx)+(seed_my*seed_my)+(seed_mz*seed_mz)) );
         qbpt.push_back( seed_q / sqrt((seed_mx*seed_mx)+(seed_my*seed_my)) );
         qbp.push_back( seed_para.signedInverseMomentum() );
         err_qbpt.push_back( err_qbp[nu_seed]/seed_sin );
         nu_seed++;
         if (debug) { cout <<" seed pt: "<<sqrt((seed_mx*seed_mx)+(seed_my*seed_my)) <<endl; }

     }
}

// read the segments associated with the seed
void SeedValidator::SegOfRecSeed( Handle<TrajectorySeedCollection> rec_seeds, int seed_idx,
                                  std::vector<SimSegment> sCSC, std::vector<SimSegment> sDT ){
 
     int idx = 0;
     d_h.clear();
     d_f.clear();
     d_x.clear();
     d_y.clear();
     geoID.clear();

     TrajectorySeedCollection::const_iterator seed_it;
     for (seed_it = rec_seeds->begin(); seed_it !=  rec_seeds->end(); seed_it++) {
 
         idx++; 
         if (seed_idx != (idx-1) ) continue;
         if (debug && scope) {  cout<<" "<<endl; }
 
         for (edm::OwnVector<TrackingRecHit>::const_iterator rh_it = seed_it->recHits().first; rh_it != seed_it->recHits().second; rh_it++) {

             const GeomDet* gdet = theService->trackingGeometry()->idToDet( (*rh_it).geographicalId() );
             GlobalPoint gp = gdet->toGlobal( (*rh_it).localPosition() );
             LocalPoint  lp = (*rh_it).localPosition();
             DetId pdid = (*rh_it).geographicalId() ; 
             

             // for parameters [1]:dx/dz, [2]:dy/dz, [3]:x, [4]:y
             double dxz = (*rh_it).parameters()[0] ;
             double dyz = (*rh_it).parameters()[1] ;
             double dz  = 1.0 / sqrt( (dxz*dxz) + (dyz*dyz) + 1.0 );
             if ( pdid.subdetId()  == MuonSubdetId::DT ) {
                dz = -1.0*dz;
             }
             double dx  = dxz*dz;
             double dy  = dyz*dz;
             LocalVector lv = LocalVector(dx,dy,dz);
             GlobalVector gv = gdet->toGlobal( lv );

	     if ( pdid.subdetId()  == MuonSubdetId::CSC ) {
                double directionSign = gp.z() * gv.z();
		lv =  (directionSign*lv).unit();
		gv = gdet->toGlobal( lv );
             }

             if (debug && scope) { 
                cout<<"============= segs from seed ============== "<<endl;
             }

	     if ( pdid.subdetId()  == MuonSubdetId::DT ) {

	        DTChamberId MB_Id = DTChamberId( pdid );

                if (debug && scope) { cout <<"DId: "<< MB_Id <<endl; }

                // find the sim-reco match case and store the difference
                if ( sDT.size() ==0 ) {
                   geoID.push_back( pdid );
                   d_h.push_back( 999.0 );
                   d_f.push_back( 999.0 );
                   d_x.push_back( 999.0 );
                   d_y.push_back( 999.0 );
                }
                else {
                   bool match=false;
                   for (std::vector<SimSegment>::const_iterator it = sDT.begin(); it != sDT.end(); it++) {
                       if ( (*it).dt_DetId == MB_Id ) {
			  geoID.push_back( pdid );
			  d_h.push_back( gv.eta() - (*it).sGlobalVec.eta() );
			  d_f.push_back( gv.phi() - (*it).sGlobalVec.phi() );
			  d_x.push_back( lp.x() - (*it).sLocalOrg.x() );
			  d_y.push_back( lp.y() - (*it).sLocalOrg.y() );
                          match = true;
                       }
                   }
                   if (!match) {
                      geoID.push_back( pdid );
		      d_h.push_back( 999.0 );
		      d_f.push_back( 999.0 );
		      d_x.push_back( 999.0 );
		      d_y.push_back( 999.0 );
                   }
                }
	     }
	     if ( pdid.subdetId()  == MuonSubdetId::CSC ) {

		CSCDetId  ME_Id = CSCDetId( pdid ) ;
                if (debug && scope) { cout <<"DId: "<< ME_Id <<endl; }

                /*/ flip the z-sign for station 1 and 2 because of the chamber orientation
                if ( (ME_Id.station()==1)||(ME_Id.station()==2) ) {
                   gv = GlobalVector( gv.x(), gv.y(), (-1.*gv.z()) );
                }*/

                // find the sim-reco match case and store the difference
                if ( sCSC.size() ==0 ) {
                   geoID.push_back( pdid );
                   d_h.push_back( 999.0 );
                   d_f.push_back( 999.0 );
                   d_x.push_back( 999.0 );
                   d_y.push_back( 999.0 );
                }
                else {
                   bool match=false;
                   for (std::vector<SimSegment>::const_iterator it = sCSC.begin(); it != sCSC.end(); it++) {
                       if ( (*it).csc_DetId == ME_Id ) {
			  geoID.push_back( pdid );
			  d_h.push_back( gv.eta() - (*it).sGlobalVec.eta() );
			  d_f.push_back( gv.phi() - (*it).sGlobalVec.phi() );
			  d_x.push_back( lp.x() - (*it).sLocalOrg.x() );
			  d_y.push_back( lp.y() - (*it).sLocalOrg.y() );
                          match = true;
                       }
                   }
                   if (!match) {
                      geoID.push_back( pdid );
		      d_h.push_back( 999.0 );
		      d_f.push_back( 999.0 );
		      d_x.push_back( 999.0 );
		      d_y.push_back( 999.0 );
                   }
                }
	     }

             if (debug && scope) {
                cout<<"  h0= "<< gp.eta() <<"  f0= "<< gp.phi()<<"   gp= "<< gp << endl;
                cout<<"  h1= "<< gv.eta() <<"  f1= "<< gv.phi()<<"   gv= "<< gv << endl;
             }
         }

     }
}

void SeedValidator::StaTrackReader( Handle<reco::TrackCollection> sta_trk, int sta_glb){

     // look at the inner most momentum and position
     nu_sta=0;
     sta_mT.clear();
     sta_mA.clear();
     sta_theta.clear();
     sta_phi.clear();
     sta_qbp.clear();
     sta_qbpt.clear();
     sta_chi2.clear();
     sta_nHits.clear();   

     TrackCollection::const_iterator iTrk;
     for ( iTrk = sta_trk->begin(); iTrk !=  sta_trk->end(); iTrk++) {
         nu_sta++;
         sta_mA.push_back( (*iTrk).p() );
         sta_mT.push_back( (*iTrk).pt() );
         sta_theta.push_back( (*iTrk).theta() );
         sta_phi.push_back( (*iTrk).phi() );
         sta_qbp.push_back( (*iTrk).qoverp() );
         sta_qbpt.push_back( ( (*iTrk).qoverp()/(*iTrk).pt() )*(*iTrk).p() );
         sta_chi2.push_back( (*iTrk).chi2() );
         sta_nHits.push_back( (*iTrk).recHitsSize() );
         if (debug) { 
            if (sta_glb==0) cout<<"sta";
            if (sta_glb==1) cout<<"glb";

            cout <<" track  pt: "<< (*iTrk).pt() <<endl; 
         }
     }
}


void SeedValidator::SimInfo(Handle<edm::SimTrackContainer> simTracks,
                            Handle<edm::PSimHitContainer> dsimHits, Handle<edm::PSimHitContainer> csimHits,
                            ESHandle<DTGeometry> dtGeom, ESHandle<CSCGeometry> cscGeom){

  // eta_c -> ave.eta from all csc stations ; cta_d -> ave.eta from dt stations
  eta_c.clear();
  eta_d.clear();
  phi_c.clear();
  phi_d.clear();
  eta_trk.clear();
  theta_trk.clear();
  phi_trk.clear();
  theQ.clear();
  pt_trk.clear();
  ptlayer.clear();
  palayer.clear();
  trackID.clear();

  for (SimTrackContainer::const_iterator simTk_It = simTracks->begin(); simTk_It != simTracks->end(); simTk_It++)
  {

      //if (abs((*simTk_It).type())!=13 || (*simTk_It).vertIndex() != 0 ) continue;
      bool rechitSize = (dsimHits->size() <8 && csimHits->size() <4) ? true:false ;
      if (abs((*simTk_It).type())!=13 || rechitSize || (*simTk_It).vertIndex() != 0 ) continue;
    
      trackID.push_back( static_cast<int>((*simTk_It).trackId())  );
      if ((*simTk_It).type()==13) {
         theQ.push_back( -1.0 );
      }else {
         theQ.push_back(  1.0 );
      }

      std::vector<double> pt1(5,0.0);
      std::vector<double> pa1(5,0.0);

      double px = ((*simTk_It).momentum()).x();
      double py = ((*simTk_It).momentum()).y();
      double pz = ((*simTk_It).momentum()).z();
      pa1[0] = sqrt( px*px + py*py + pz*pz );
      pt1[0] = sqrt( px*px + py*py );

      eta_trk.push_back( getEta(px,py,pz)  );
      theta_trk.push_back( acos(pz/pa1[0])  );
      phi_trk.push_back( atan2(py,px) );
      pt_trk.push_back( pt1[0] );
   
      double eta_d1 = 0.0;
      double phi_d1 = 0.0;
      double enu2   = 0.0;
      for (PSimHitContainer::const_iterator ds_It = dsimHits->begin(); ds_It != dsimHits->end(); ds_It++)
      {          
          Local3DPoint lp = (*ds_It).localPosition(); 

          DTLayerId D_Id = DTLayerId( (*ds_It).detUnitId() );
          const DTLayer* dtlayer = dtGeom->layer(D_Id);
          GlobalVector m2 = dtlayer->toGlobal((*ds_It).momentumAtEntry() );
          GlobalPoint gp = dtlayer->toGlobal(lp );

          if ( ( abs((*ds_It).particleType())==13 ) && ( (*ds_It).trackId()==(*simTk_It).trackId() )) {
 
             pt1[ D_Id.station() ] = sqrt( (m2.x()*m2.x()) + (m2.y()*m2.y()) );
             pa1[ D_Id.station() ] = sqrt( (m2.x()*m2.x()) + (m2.y()*m2.y()) + (m2.z()*m2.z()) );

             eta_d1 += gp.eta();
             phi_d1 += gp.phi();
             enu2  += 1.0;
          }
      }
      if (enu2 !=0.0 ) {
         eta_d.push_back( eta_d1/enu2 );
         phi_d.push_back( phi_d1/enu2 );
      } else {
         eta_d.push_back( -99.0 );
         phi_d.push_back( -99.0 );
      }

      double eta_c1 = 0.0;
      double phi_c1 = 0.0;
      double enu1   = 0.0;
      for (PSimHitContainer::const_iterator cs_It = csimHits->begin(); cs_It != csimHits->end(); cs_It++)
      {
          CSCDetId C_Id = CSCDetId((*cs_It).detUnitId());
          const CSCChamber* cscchamber = cscGeom->chamber( C_Id );
          GlobalVector m1 = cscchamber->toGlobal((*cs_It).momentumAtEntry() );
          Local3DPoint lp = (*cs_It).localPosition(); 
          GlobalPoint gp = cscchamber->toGlobal(lp );

          if ( ( abs((*cs_It).particleType())==13 ) && ( (*cs_It).trackId()==(*simTk_It).trackId() )) {

             if (enu2 == 0.0) {
                pt1[C_Id.station()] = sqrt( (m1.x()*m1.x()) + (m1.y()*m1.y()) ) ; 
                pa1[C_Id.station()] = sqrt( (m1.x()*m1.x()) + (m1.y()*m1.y()) + (m1.z()*m1.z()) );
             }

             eta_c1 += gp.eta();
             phi_c1 += gp.phi();
             enu1   += 1.0;
          }
      }
      if (enu1 !=0.0 ) {
          eta_c.push_back( eta_c1/enu1 );
          phi_c.push_back( phi_c1/enu1 );
      } else {
         eta_c.push_back( -99.0 );
         phi_c.push_back( -99.0 );
      }

      ptlayer.push_back(pt1);
      palayer.push_back(pa1);
  }

}

// Look up what segments we have in a event
int SeedValidator::RecSegReader( Handle<CSCSegmentCollection> cscSeg, Handle<DTRecSegment4DCollection> dtSeg                                , ESHandle<CSCGeometry> cscGeom, ESHandle<DTGeometry> dtGeom, double trkTheta, double trkPhi) {

     // Calculate the ave. eta & phi
     ave_phi = 0.0;
     ave_eta = 0.0;
     phi_resid.clear();
     eta_resid.clear();

     double n=0.0;
     double m=0.0;
     for(CSCSegmentCollection::const_iterator it = cscSeg->begin(); it != cscSeg->end(); it++)
     {
        if ( (*it).nRecHits() < 4) continue;
        CSCDetId DetId = (CSCDetId)(*it).cscDetId();
	const CSCChamber* cscchamber = cscGeom->chamber( DetId );
	GlobalPoint gp = cscchamber->toGlobal((*it).localPosition() );
        if (( fabs(gp.theta()- trkTheta) > 0.35  ) || ( fabs(gp.phi()- trkPhi) > 0.35)  ) continue;
	GlobalVector gv = cscchamber->toGlobal((*it).localDirection() );


        ave_phi += gp.phi();
        ave_eta += gp.eta();
        dx_error.push_back( (*it).localDirectionError().xx() );
        dy_error.push_back( (*it).localDirectionError().yy() );
        x_error.push_back( (*it).localPositionError().xx() );
        y_error.push_back( (*it).localPositionError().yy() );
        n++;
        if (debug) {
           cout <<"~~~~~~~~~~~~~~~~  reco segs  ~~~~~~~~~~~~~~~~  " <<endl;
	   cout <<"DId: "<<DetId<<endl;
	   cout <<"  h0= "<<gp.eta()<<"  f0= "<<gp.phi()<<"   gp= "<< gp <<endl;
	   cout <<"  h1= "<<gv.eta()<<"  f1= "<<gv.phi()<<"   gv= "<< gv <<endl;
        }
     }
     for(DTRecSegment4DCollection::const_iterator it = dtSeg->begin(); it != dtSeg->end(); it++)
     {
        if ( !(*it).hasPhi() || !(*it).hasZed()  ) continue;
        DTChamberId DetId = (*it).chamberId();
        const DTChamber* dtchamber = dtGeom->chamber( DetId );
        GlobalPoint  gp = dtchamber->toGlobal( (*it).localPosition() );
        if (( fabs(gp.eta()- trkTheta) > 0.35  ) || ( fabs(gp.phi()- trkPhi) > 0.35 ) ) continue;
	GlobalVector gv = dtchamber->toGlobal((*it).localDirection() );

        ave_phi += gp.phi();
        ave_eta += gp.eta();
        m++;
        if (debug) {
           cout <<"~~~~~~~~~~~~~~~~  reco segs  ~~~~~~~~~~~~~~~~  " <<endl;
	   cout <<"DId: "<<DetId<<endl;
	   cout <<"  h0= "<<gp.eta()<<"  f0= "<<gp.phi()<<"   gp= "<< gp <<endl;
	   cout <<"  h1= "<<gv.eta()<<"  f1= "<<gv.phi()<<"   gv= "<< gv <<endl;
        }
     }
     ave_phi = ave_phi / (n+m) ;
     ave_eta = ave_eta / (n+m) ;

     // Calculate the residual of phi and eta
     for(CSCSegmentCollection::const_iterator it = cscSeg->begin(); it != cscSeg->end(); it++)
     {
        if ( (*it).nRecHits() < 4) continue;
        CSCDetId DetId = (CSCDetId)(*it).cscDetId();
        const CSCChamber* cscchamber = cscGeom->chamber( DetId );
        GlobalPoint gp = cscchamber->toGlobal((*it).localPosition() );
        if ( (fabs(gp.theta()- trkTheta) > 0.35  ) || ( fabs(gp.phi()- trkPhi) > 0.35 ) ) continue;
        phi_resid.push_back( gp.phi()- ave_phi );
        eta_resid.push_back( gp.eta()- ave_eta );
     }
     for(DTRecSegment4DCollection::const_iterator it = dtSeg->begin(); it != dtSeg->end(); it++)
     {
        if ( !(*it).hasPhi() || !(*it).hasZed()  ) continue;
        DTChamberId DetId = (*it).chamberId();
        const DTChamber* dtchamber = dtGeom->chamber( DetId );
        GlobalPoint  gp = dtchamber->toGlobal( (*it).localPosition() );
        if ( (fabs(gp.eta()- trkTheta) > 0.35  ) || ( fabs(gp.phi()- trkPhi) > 0.35 ) ) continue;
        phi_resid.push_back( gp.phi()- ave_phi );
        eta_resid.push_back( gp.eta()- ave_eta );
     }

     if (n!=0 && m== 0) {
        return 1; // csc 
     } else if (n==0 && m!=0 ) {
        return 3; // dt
     } else if (n!=0 && m!=0 ) {
        return 2; // overlap
     } else {
        return 0; // fail
     }

}

double SeedValidator::getEta(double vx, double vy, double vz ) {

      double va = sqrt( vx*vx + vy*vy + vz*vz );

      double theta = acos( vz/va );
      double eta = (-1.0)*log( tan(theta/2.0) )  ;
      return eta;
}
double SeedValidator::getEta(double theta ) {

      double eta = (-1.0)*log( tan(theta/2.0) )  ;
      return eta;
}
