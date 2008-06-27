// Class Header
#include "RecoMuon/SeedGenerator/test/MuonSeedPTAnalysis/MuonSeedPtEstimate.h"

// Data Formats 
#include <SimDataFormats/Track/interface/SimTrackContainer.h>
#include <SimDataFormats/TrackingHit/interface/PSimHitContainer.h>
#include <DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h>
#include <DataFormats/CSCRecHit/interface/CSCRecHit2D.h>
#include <DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h>
#include <DataFormats/DTRecHit/interface/DTRecSegment4D.h>

// Geometry
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "TrackingTools/DetLayers/interface/DetLayer.h"
#include <RecoMuon/MeasurementDet/interface/MuonDetLayerMeasurements.h>
#include <RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h>
#include <RecoMuon/Records/interface/MuonRecoGeometryRecord.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/CSCGeometry/interface/CSCChamber.h>

// Framework
#include "Geometry/Vector/interface/GlobalPoint.h"
#include "Geometry/Vector/interface/GlobalVector.h"
#include "Geometry/Vector/interface/LocalPoint.h"
#include "Geometry/Vector/interface/LocalVector.h"
 
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "TFile.h"
#include "TVector3.h"

#include <iostream>
#include <string>
#include <map>

using namespace std;
using namespace edm;


// constructors
MuonSeedPtEstimate::MuonSeedPtEstimate(const ParameterSet& pset){ 


  // Local Debug flag
  debug                = pset.getParameter<bool>("DebugMuonSeed");

  // Output file for histograms
  rootFileName         = pset.getUntrackedParameter<string>("rootFileName");

  // enable the DT chamber
  enableDTMeasurement  = pset.getParameter<bool>("EnableDTMeasurement");
  theDTSegmentLabel    = pset.getParameter<edm::InputTag>("DTSegmentLabel");

  // enable the CSC chamber
  enableCSCMeasurement = pset.getParameter<bool>("EnableCSCMeasurement");
  theCSCSegmentLabel   = pset.getParameter<edm::InputTag>("CSCSegmentLabel");

  // Parameters for seed creation in endcap region
  minCSCHitsPerSegment = pset.getParameter<int>("minCSCHitsPerSegment");

  // Create the root file
  theFile = new TFile(rootFileName.c_str(), "RECREATE");
  theFile->cd();
  
  // Book the histograms

  hME[0] = new HSeedPt1("ME11_1");  
  hME[1] = new HSeedPt1("ME11_2");  
  hME[2] = new HSeedPt1("ME11_3");
  hME[3] = new HSeedPt1("ME11_4");  
  hME[4] = new HSeedPt1("ME12_2");
  hME[5] = new HSeedPt1("ME12_3");
  hME[6] = new HSeedPt1("ME12_4");
  hME[7] = new HSeedPt1("ME2_3");
  hME[8] = new HSeedPt1("ME2_4");
  hME[9] = new HSeedPt1("ME3_4");

	  
  hMB[0] = new HSeedPt2("MB1_2");  
  hMB[1] = new HSeedPt2("MB1_3");
  hMB[2] = new HSeedPt2("MB1_4");
  hMB[3] = new HSeedPt2("MB2_3");
  hMB[4] = new HSeedPt2("MB2_4");
  hMB[5] = new HSeedPt2("MB3_4");
}


// destructor
MuonSeedPtEstimate::~MuonSeedPtEstimate(){

  // Write the histos to file
  theFile->cd();
  for (int i = 0; i < 10; ++i) hME[i]->Write(); 
  for (int i = 0; i <  6; ++i) hMB[i]->Write(); 
  theFile->Close();
}


// The Main
void MuonSeedPtEstimate::analyze(const Event& event, const EventSetup& eSetup)
{
  
  // Muon Geometry
  edm::ESHandle<MuonDetLayerGeometry> muonLayers;
  eSetup.get<MuonRecoGeometryRecord>().get(muonLayers);
  
  edm::ESHandle<CSCGeometry> h;
  eSetup.get<MuonGeometryRecord>().get(h);
  cscgeom = &*h;
  
  // Simtracks
  edm::Handle<edm::SimTrackContainer> simTracks;
  event.getByLabel("g4SimHits",simTracks);
  
  // Get the SimHit collection for CSC because of ME1/a strips which are ganged...
  Handle<PSimHitContainer> cscSimHits;
  event.getByLabel("g4SimHits","MuonCSCHits", cscSimHits);
  
  
  float SimEta = 999.;
  float SimPt = 0.;
  
  // Get the true pt of track
  for (edm::SimTrackContainer::const_iterator it = simTracks->begin(); it != simTracks->end(); it++) {
    if (abs((*it).type()) != 13) continue;
    float px = (*it).momentum().x();
    float py = (*it).momentum().y();
    SimEta   = (*it).momentum().eta();
    SimPt    = sqrt(px*px + py*py);
  }    
  
  
  //  if (debug) std::cout << "Generated track has eta: " << SimEta << " and Pt: " << SimPt << std::endl;
  
  // Now obtain the reconstructed segments
  
  // Instantiate the accessor (get the segments: DT + CSC but not RPC=false)
  MuonDetLayerMeasurements muonMeasurements(enableDTMeasurement,enableCSCMeasurement,false,
					    theDTSegmentLabel.label(),theCSCSegmentLabel.label());
  
  
  // Get the various stations and store segments in containers for each station (layers)
  
  // Deal with Endcap pairs first:
  SegmentContainer CSClist4;
  SegmentContainer CSClist3; 
  SegmentContainer CSClist2; 
  SegmentContainer CSClist1; 
  SegmentContainer CSClist0; 
  
  std::vector<DetLayer*> cscAllLayers; 
  
  if (SimEta < 0.) {
    cscAllLayers = muonLayers->backwardCSCLayers();
  } else { 
    cscAllLayers = muonLayers->forwardCSCLayers();
  }
  
  CSClist4 = muonMeasurements.recHits( cscAllLayers[4], event );
  CSClist3 = muonMeasurements.recHits( cscAllLayers[3], event );
  CSClist2 = muonMeasurements.recHits( cscAllLayers[2], event );
  CSClist1 = muonMeasurements.recHits( cscAllLayers[1], event ); // ME1/2 and 1/3
  CSClist0 = muonMeasurements.recHits( cscAllLayers[0], event ); // ME11
     

  GlobalPoint gp1, gp2;
  float dpsi[10]   = { 99., 99., 99., 99., 99., 99., 99., 99., 99., 99.};
  float dphi[10]   = { 99., 99., 99., 99., 99., 99., 99., 99., 99., 99.};
  float RecEta[10] = { 99., 99., 99., 99., 99., 99., 99., 99., 99., 99.};
  int idx = 0;

  // Loop over ME1/1 segments first:
  for (SegmentContainer::iterator it = CSClist0.begin(); it != CSClist0.end(); ++it ){
    if ( int ((*it)->recHits().size()) < minCSCHitsPerSegment ) continue; 
    gp1 = (*it)->globalPosition();
    if ( !cscSimMatch(cscSimHits, CSClist0, idx) ) continue;
    idx++;


    GlobalVector gv = (*it)->globalDirection();
    
    // Psi is angle between the segment origin and segment direction
    // Use dot product between two vectors to get Psi in global x-y plane
    double cosDpsi  = (gv.x()*gp1.x() + gv.y()*gp1.y());
    cosDpsi /= sqrt(gp1.x()*gp1.x() + gp1.y()*gp1.y());
    cosDpsi /= sqrt(gv.x()*gv.x() + gv.y()*gv.y());
    
    double axb = ( gp1.x()*gv.y() ) - ( gp1.y()*gv.x() ) ;
    double sign = (axb < 0.) ? 1.0 : -1.0;
    
    dpsi[0] = acos(cosDpsi) * sign;


    // Loop over ME1/2 segments: 
    for (SegmentContainer::iterator it2 = CSClist1.begin(); it2 != CSClist1.end(); ++it2 ){
      if ( int ((*it2)->recHits().size()) < minCSCHitsPerSegment ) continue; 
      gp2 = (*it2)->globalPosition();
      RecEta[0] = gp2.eta();
      dphi[0] = gp1.phi() - gp2.phi();
    }

    // Loop over ME2 segments: 
    for (SegmentContainer::iterator it2 = CSClist2.begin(); it2 != CSClist2.end(); ++it2 ){
      if ( int ((*it2)->recHits().size()) < minCSCHitsPerSegment ) continue; 
      gp2 = (*it2)->globalPosition();
      RecEta[1] = gp2.eta();
      dphi[1] = gp1.phi() - gp2.phi();
    }

    // Loop over ME3 segments: 
    for (SegmentContainer::iterator it2 = CSClist3.begin(); it2 != CSClist3.end(); ++it2 ){
      if ( int ((*it2)->recHits().size()) < minCSCHitsPerSegment ) continue; 
      gp2 = (*it2)->globalPosition();
      RecEta[2] = gp2.eta();
      dphi[2] = gp1.phi() - gp2.phi();
    }

    // Loop over ME4 segments: 
    for (SegmentContainer::iterator it2 = CSClist4.begin(); it2 != CSClist4.end(); ++it2 ){
      if ( int ((*it2)->recHits().size()) < minCSCHitsPerSegment ) continue; 
      gp2 = (*it2)->globalPosition();
      RecEta[3] = gp2.eta();
      dphi[3] = gp1.phi() - gp2.phi();
    }
  }


  // Loop over ME1/2 and 1/3 segments first:

  for (SegmentContainer::iterator it = CSClist1.begin(); it != CSClist1.end(); ++it ){
    if ( int ((*it)->recHits().size()) < minCSCHitsPerSegment ) continue; 
    gp1 = (*it)->globalPosition();

    GlobalVector gv = (*it)->globalDirection();
    
    // Psi is angle between the segment origin and segment direction
    // Use dot product between two vectors to get Psi in global x-y plane
    double cosDpsi  = (gv.x()*gp1.x() + gv.y()*gp1.y());
    cosDpsi /= sqrt(gp1.x()*gp1.x() + gp1.y()*gp1.y());
    cosDpsi /= sqrt(gv.x()*gv.x() + gv.y()*gv.y());
    
    double axb = ( gp1.x()*gv.y() ) - ( gp1.y()*gv.x() ) ;
    double sign = (axb < 0.) ? 1.0 : -1.0;
    
    dpsi[1] = acos(cosDpsi) * sign;


    // Loop over ME2 segments: 
    for (SegmentContainer::iterator it2 = CSClist2.begin(); it2 != CSClist2.end(); ++it2 ){
      if ( int ((*it2)->recHits().size()) < minCSCHitsPerSegment ) continue; 
      gp2 = (*it2)->globalPosition();
      RecEta[4] = gp2.eta();
      dphi[4] = gp1.phi() - gp2.phi();
    }


    // Loop over ME3 segments: 
    for (SegmentContainer::iterator it2 = CSClist3.begin(); it2 != CSClist3.end(); ++it2 ){
      if ( int ((*it2)->recHits().size()) < minCSCHitsPerSegment ) continue; 
      gp2 = (*it2)->globalPosition();
      RecEta[5] = gp2.eta();
      dphi[5] = gp1.phi() - gp2.phi();
    }

    // Loop over ME4 segments: 
    for (SegmentContainer::iterator it2 = CSClist4.begin(); it2 != CSClist4.end(); ++it2 ){
      if ( int ((*it2)->recHits().size()) < minCSCHitsPerSegment ) continue; 
      gp2 = (*it2)->globalPosition();
      RecEta[6] = gp2.eta();
      dphi[6] = gp1.phi() - gp2.phi();
    }
  }


  // Loop over ME2 segments first:

  for (SegmentContainer::iterator it = CSClist2.begin(); it != CSClist2.end(); ++it ){
    if ( int ((*it)->recHits().size()) < minCSCHitsPerSegment ) continue; 
    gp1 = (*it)->globalPosition();


    GlobalVector gv = (*it)->globalDirection();
    
    // Psi is angle between the segment origin and segment direction
    // Use dot product between two vectors to get Psi in global x-y plane
    double cosDpsi  = (gv.x()*gp1.x() + gv.y()*gp1.y());
    cosDpsi /= sqrt(gp1.x()*gp1.x() + gp1.y()*gp1.y());
    cosDpsi /= sqrt(gv.x()*gv.x() + gv.y()*gv.y());
    
    double axb = ( gp1.x()*gv.y() ) - ( gp1.y()*gv.x() ) ;
    double sign = (axb < 0.) ? 1.0 : -1.0;
    
    dpsi[3] = acos(cosDpsi) * sign;

    // Loop over ME3 segments: 
    for (SegmentContainer::iterator it2 = CSClist3.begin(); it2 != CSClist3.end(); ++it2 ){
      if ( int ((*it2)->recHits().size()) < minCSCHitsPerSegment ) continue; 
      gp2 = (*it2)->globalPosition();
      RecEta[7] = gp2.eta();
      dphi[7] = gp1.phi() - gp2.phi();
    }

    // Loop over ME4 segments: 
    for (SegmentContainer::iterator it2 = CSClist4.begin(); it2 != CSClist4.end(); ++it2 ){
      if ( int ((*it2)->recHits().size()) < minCSCHitsPerSegment ) continue; 
      gp2 = (*it2)->globalPosition();
      RecEta[8] = gp2.eta();
      dphi[8] = gp1.phi() - gp2.phi();
    }
  }

  // Loop over ME3 segments first:
  for (SegmentContainer::iterator it = CSClist3.begin(); it != CSClist3.end(); ++it ){
    if ( int ((*it)->recHits().size()) < minCSCHitsPerSegment ) continue; 
    gp1 = (*it)->globalPosition();
 
    GlobalVector gv = (*it)->globalDirection();
    
    // Psi is angle between the segment origin and segment direction
    // Use dot product between two vectors to get Psi in global x-y plane
    double cosDpsi  = (gv.x()*gp1.x() + gv.y()*gp1.y());
    cosDpsi /= sqrt(gp1.x()*gp1.x() + gp1.y()*gp1.y());
    cosDpsi /= sqrt(gv.x()*gv.x() + gv.y()*gv.y());
    
    double axb = ( gp1.x()*gv.y() ) - ( gp1.y()*gv.x() ) ;
    double sign = (axb < 0.) ? 1.0 : -1.0;
    
    dpsi[4] = acos(cosDpsi) * sign;

    // Loop over ME4 segments: 
    for (SegmentContainer::iterator it2 = CSClist4.begin(); it2 != CSClist4.end(); ++it2 ){
      if ( int ((*it2)->recHits().size()) < minCSCHitsPerSegment ) continue; 
      gp2 = (*it2)->globalPosition();
      RecEta[9] = gp2.eta();
      dphi[9] = gp1.phi() - gp2.phi();
    }
  }
  
  // Fill Encap Histograms:
  for (int i = 0; i < 10; ++i) hME[i]->FillCSC(SimPt, RecEta[i], dphi[i], dpsi[i]);



  // Now look over barrel

  SegmentContainer DTlist4; 
  SegmentContainer DTlist3;
  SegmentContainer DTlist2; 
  SegmentContainer DTlist1; 
  
  
  if ( fabs(SimEta) < 1.1) {
    // Get the DT segments by stations (layers):
    std::vector<DetLayer*> dtLayers = muonLayers->allDTLayers();
    DTlist4= muonMeasurements.recHits( dtLayers[3], event );
    DTlist3= muonMeasurements.recHits( dtLayers[2], event );
    DTlist2= muonMeasurements.recHits( dtLayers[1], event );
    DTlist1= muonMeasurements.recHits( dtLayers[0], event );
  }

  float dphiB[6]   = { 99., 99., 99., 99., 99., 99.};
  float dpsiB[6]   = { 99., 99., 99., 99., 99., 99.};
  float RecEtaB[6] = { 99., 99., 99., 99., 99., 99.};

  // Loop over MB1 segments first:
  for (SegmentContainer::iterator it = DTlist1.begin(); it != DTlist1.end(); ++it ){
    gp1 = (*it)->globalPosition();

    GlobalVector gv = (*it)->globalDirection();
    
    // Psi is angle between the segment origin and segment direction
    // Use dot product between two vectors to get Psi in global x-y plane
    double cosDpsi  = (gv.x()*gp1.x() + gv.y()*gp1.y());
    cosDpsi /= sqrt(gp1.x()*gp1.x() + gp1.y()*gp1.y());
    cosDpsi /= sqrt(gv.x()*gv.x() + gv.y()*gv.y());
    
    double axb = ( gp1.x()*gv.y() ) - ( gp1.y()*gv.x() ) ;
    double sign = (axb < 0.) ? 1.0 : -1.0;
    
    dpsiB[0] = acos(cosDpsi) * sign;
    
    // Loop over MB2 segments: 
    for (SegmentContainer::iterator it2 = DTlist2.begin(); it2 != DTlist2.end(); ++it2 ){
      gp2 = (*it2)->globalPosition();
      RecEtaB[0] = gp2.eta();
      dphiB[0] = gp1.phi() - gp2.phi();
    }

    // Loop over MB3 segments: 
    for (SegmentContainer::iterator it2 = DTlist3.begin(); it2 != DTlist3.end(); ++it2 ){
      gp2 = (*it2)->globalPosition();
      RecEtaB[1] = gp2.eta();
      dphiB[1] = gp1.phi() - gp2.phi();
    }

    // Loop over MB4 segments: 
    for (SegmentContainer::iterator it2 = DTlist4.begin(); it2 != DTlist4.end(); ++it2 ){
      gp2 = (*it2)->globalPosition();
      RecEtaB[2] = gp2.eta();
      dphiB[2] = gp1.phi() - gp2.phi();
    }
  }

  // Loop over MB2 segments first:
  for (SegmentContainer::iterator it = DTlist2.begin(); it != DTlist2.end(); ++it ){
    gp1 = (*it)->globalPosition();

    GlobalVector gv = (*it)->globalDirection();
    
    // Psi is angle between the segment origin and segment direction
    // Use dot product between two vectors to get Psi in global x-y plane
    double cosDpsi  = (gv.x()*gp1.x() + gv.y()*gp1.y());
    cosDpsi /= sqrt(gp1.x()*gp1.x() + gp1.y()*gp1.y());
    cosDpsi /= sqrt(gv.x()*gv.x() + gv.y()*gv.y());
    
    double axb = ( gp1.x()*gv.y() ) - ( gp1.y()*gv.x() ) ;
    double sign = (axb < 0.) ? 1.0 : -1.0;
    
    dpsiB[1] = acos(cosDpsi) * sign;
    
    // Loop over MB3 segments: 
    for (SegmentContainer::iterator it2 = DTlist3.begin(); it2 != DTlist3.end(); ++it2 ){
      gp2 = (*it2)->globalPosition();
      RecEtaB[3] = gp2.eta();
      dphiB[3] = gp1.phi() - gp2.phi();
    }

    // Loop over MB4 segments: 
    for (SegmentContainer::iterator it2 = DTlist4.begin(); it2 != DTlist4.end(); ++it2 ){
      gp2 = (*it2)->globalPosition();
      RecEtaB[4] = gp2.eta();
      dphiB[4] = gp1.phi() - gp2.phi();
    }
  }

  // Loop over MB3 segments first:
  for (SegmentContainer::iterator it = DTlist3.begin(); it != DTlist3.end(); ++it ){
    gp1 = (*it)->globalPosition();

    GlobalVector gv = (*it)->globalDirection();
    
    // Psi is angle between the segment origin and segment direction
    // Use dot product between two vectors to get Psi in global x-y plane
    double cosDpsi  = (gv.x()*gp1.x() + gv.y()*gp1.y());
    cosDpsi /= sqrt(gp1.x()*gp1.x() + gp1.y()*gp1.y());
    cosDpsi /= sqrt(gv.x()*gv.x() + gv.y()*gv.y());
    
    double axb = ( gp1.x()*gv.y() ) - ( gp1.y()*gv.x() ) ;
    double sign = (axb < 0.) ? 1.0 : -1.0;
    
    dpsiB[2] = acos(cosDpsi) * sign;
    
    // Loop over MB4 segments: 
    for (SegmentContainer::iterator it2 = DTlist4.begin(); it2 != DTlist4.end(); ++it2 ){
      gp2 = (*it2)->globalPosition();
      RecEtaB[5] = gp2.eta();
      dphiB[5] = gp1.phi() - gp2.phi();
    }
  }
  
  // Fill Barrel Histograms:
  for (int i = 0; i < 6; ++i) hMB[i]->FillDT(SimPt, RecEtaB[i], dphiB[i], dpsiB[i]);


}



// Find matching ME1/a segment to Simhits
bool MuonSeedPtEstimate::cscSimMatch(const edm::Handle<edm::PSimHitContainer> simHits , SegmentContainer CSClist0, int idx) {

  // Nightmare:  don't have access to DetId in these pseudo segments, so we don't know if ME1/a or ME1/b...   
  // ... but if only have two or less reconstructed segment --> ME1/b
  if (CSClist0.size() < 3) return true;
  // ... if have more than 3 segments, reject
  if (CSClist0.size() > 3) return false;

  // Hence, have 3 segments --> most likely to be ME1/a...

  int idx2 = 0;
  int MatchIdx = idx;
  float DeltaPhiMin = 999.;


  // Loop over ME1/1 segments
  for (SegmentContainer::iterator its = CSClist0.begin(); its != CSClist0.end(); ++its ){
    float gphi = (*its)->globalPosition().phi();

    float avgSimPhi = 0.;
    int size = 0;

    // loop over simhits and compute average global phi for ME1/a  --> loose match with reco segment
    for (PSimHitContainer::const_iterator ith = simHits->begin(); ith != simHits->end(); ith++) {      
      CSCDetId simId = (CSCDetId)(*ith).detUnitId();
      if ((simId.station() == 1) &&
	  (simId.ring()    == 4)) {  
	const CSCChamber* ch = cscgeom->chamber(simId);
	GlobalPoint gp = ch->toGlobal((*ith).localPosition());
	float simphi = gp.phi();
	avgSimPhi += simphi;
	size++;
      }
    }
    if (size > 0) {
      avgSimPhi = avgSimPhi/size;
    } else {
      return false;
    }
    float DeltaPhi = fabs(avgSimPhi - gphi);
    if (DeltaPhi < DeltaPhiMin) {
      DeltaPhiMin = DeltaPhi;
      MatchIdx = idx2;
    }
    
    idx2++;
  }

  if (MatchIdx == idx) return true;

  return false;
}

DEFINE_FWK_MODULE(MuonSeedPtEstimate);
