#ifndef RecoMuon_MuonSeedPtEstimate_H
#define RecoMuon_MuonSeedPtEstimate_H

/** \class SeedPtEstimate
 *
 *  Author: Dominique Fortin  - UC Riverside
 */

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Handle.h"
#include <FWCore/ParameterSet/interface/ParameterSet.h>

#include <SimDataFormats/TrackingHit/interface/PSimHitContainer.h>
#include <RecoMuon/MeasurementDet/interface/MuonDetLayerMeasurements.h>
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>

#include "RecoMuon/SeedGenerator/test/MuonSeedPTAnalysis/SeedPtEstimateHisto.h"

#include <vector>
#include <map>
#include <string>

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

//class PSimHit;
class TFile;
class DetLayer; 
class MuonDetLayerGeometry;


class MuonSeedPtEstimate: public edm::EDAnalyzer {

 public:

  typedef MuonTransientTrackingRecHit::MuonRecHitContainer SegmentContainer;

  /// Constructor
  MuonSeedPtEstimate(const edm::ParameterSet& pset);

  /// Destructor
  virtual ~MuonSeedPtEstimate();

  // Operations

  /// Perform the real analysis
  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);


protected:

private: 

  /// Find proper ME1/a segment out of 3
  bool cscSimMatch(const edm::Handle<edm::PSimHitContainer> cscsimHits , SegmentContainer cscSeg, int idx);

  // Histograms
  HSeedPt1 *hME[10]; 	  
  HSeedPt2 *hMB[6]; 


  // The file which will store the histos
  TFile *theFile;

 // Switch for debug output
  bool debug;
  
  std::string rootFileName;
  std::string simHitLabel;
  std::string simTrackLabel;


  // Enable the DT measurement
  bool enableDTMeasurement;

  // Enable the CSC measurement
  bool enableCSCMeasurement;

  // Minimum # of hits to consider a CSC Segment;
  int minCSCHitsPerSegment;

  /// Name of the DT segment collection
  edm::InputTag theDTSegmentLabel;

  /// Name of the CSC segment collection
  edm::InputTag theCSCSegmentLabel;

  const CSCGeometry* cscgeom;

};


#endif

