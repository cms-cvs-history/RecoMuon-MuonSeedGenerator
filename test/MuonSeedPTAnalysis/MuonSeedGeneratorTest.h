#ifndef RecoMuon_MuonSeedGeneratorTest_H
#define RecoMuon_MuonSeedGeneratorTest_H

/** \class MuonSeedGeneratorTest
 *
 *  Author: S.C. Kao  - UC Riverside
 */

//#include "RecoMuon/MuonSeedGenerator/test/MuonSeedPTAnalysis/SegSelector.h"
#include "MuonSeedGeneratorHistograms.h"
#include "MuonSeedGeneratorNtuple.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include <DataFormats/Common/interface/Handle.h>

#include <DataFormats/MuonDetId/interface/CSCDetId.h>
#include <DataFormats/MuonDetId/interface/DTChamberId.h>
#include <DataFormats/CSCRecHit/interface/CSCRecHit2D.h>
#include <DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h>
#include <DataFormats/CSCRecHit/interface/CSCSegment.h>
#include <DataFormats/CSCRecHit/interface/CSCSegmentCollection.h>
#include <DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h>
#include <DataFormats/DTRecHit/interface/DTRecSegment4D.h>
#include <DataFormats/DTRecHit/interface/DTRecSegment2DCollection.h>
#include <DataFormats/DTRecHit/interface/DTRecSegment2D.h>
#include <DataFormats/DTRecHit/interface/DTRecHitCollection.h>
#include <DataFormats/DTRecHit/interface/DTRecHit1D.h>

#include <SimDataFormats/Track/interface/SimTrackContainer.h>
#include <SimDataFormats/TrackingHit/interface/PSimHitContainer.h>

#include <Geometry/CSCGeometry/interface/CSCGeometry.h>
#include <Geometry/CSCGeometry/interface/CSCChamber.h>
#include <Geometry/CSCGeometry/interface/CSCLayer.h>
#include <Geometry/CSCGeometry/interface/CSCLayerGeometry.h>
#include <Geometry/DTGeometry/interface/DTGeometry.h>
#include <Geometry/DTGeometry/interface/DTChamber.h>
#include <Geometry/DTGeometry/interface/DTLayer.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
	
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include <DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h>
#include <DataFormats/TrajectorySeed/interface/TrajectorySeed.h>

#include <vector>
#include <map>
#include <string>
#include <utility> 

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

//class PSimHit;
class TFile;
class CSCLayer;
class CSCDetId;
class DTLayerId;
class DTSuperLayerId;
class DTChamberId;
class SegSelector;
class MuonSeedBuilder;
class SeedBuilder;

class MuonSeedGeneratorTest : public edm::EDAnalyzer {
public:

  /// Constructor
  MuonSeedGeneratorTest(const edm::ParameterSet& pset);

  /// Destructor
  virtual ~MuonSeedGeneratorTest();

  // Operations
  /// Perform the real analysis
  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);


protected:

private:
 
  SegSelector* recsegSelector;
  /// Builder where seeds are formed
  MuonSeedBuilder* muonSeedBuilder_;
  SeedBuilder* SeedBuilder_;

  // Utility functions
  void CSCsegment_stat(edm::Handle<CSCSegmentCollection> cscSeg);
  void DTsegment_stat(edm::Handle<DTRecSegment4DCollection> dtSeg);

  void SeedFromCSCRecHit(edm::Handle<CSCRecHit2DCollection> cscrechit, edm::ESHandle<CSCGeometry> cscGeom);
  void SeedFromDTRecHit(edm::Handle<DTRecHitCollection> dtrechit, edm::ESHandle<DTGeometry> dtGeom);
  
  void Eta_Test(const edm::Handle<edm::SimTrackContainer> simTracks,
                const edm::Handle<edm::PSimHitContainer> simHits,
                      edm::ESHandle<DTGeometry> dtGeom);
  int ChargeAssignment(GlobalVector Va, GlobalVector Vb);
 
  double pT_estimation(double p0, double p1, double the_eta, double the_phi);
  bool SameChamber(CSCDetId SimDetId, CSCDetId SegDetId);
  void RecSeedReader(edm::Handle<TrajectorySeedCollection> rec_seeds, edm::ESHandle<DTGeometry> dtGeom, edm::ESHandle<CSCGeometry> cscGeom);

  // Histograms
  H2DRecHit1 *h_all;
  H2DRecHit2 *h_csc;
  H2DRecHit3 *h_dt;
  H2DRecHit4 *hME1[17];
  H2DRecHit5 *hMB1[30];
  H2DRecHit6 *hME2[8];
  H2DRecHit7 *hMB2[12];
  H2DRecHit8 *hME3[6];
  H2DRecHit9 *hMB3[6];
  H2DRecHit10 *hOL1[11];
  TNtuple1 *tr_muon;

  // The file which will store the histos
  TFile *theFile;


  //cscsegment_stat output
  int cscseg_stat[6];
  int cscseg_stat1[6];
  //dtsegment_stat output
  int dtseg_stat[6];
  int dtseg_stat1[6];

  // SeedfromRecHit
  //std::vector<CSCRecHit2D> csc_rh;
  int cscrh_sum[6];
  int dtrh_sum[6];
  // Eta testing
  double h_trk;
  double h_sl1[5];
  double h_sl2[5];
  double h_sl3[5];
  // reco-seeding reader
  GlobalPoint seed_gp;
  int nSegInSeed;

  // Switch for debug output
  bool debug;

  std::string rootFileName;
  std::string cscSegmentLabel;
  std::string recHitLabel;
  std::string dtSegmentLabel;
  std::string dtrecHitLabel;
  std::string simHitLabel;
  std::string simTrackLabel;
  std::string muonseedLabel;

};


#endif

