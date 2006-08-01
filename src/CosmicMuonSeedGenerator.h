#ifndef MuonSeedGenerator_CosmicMuonSeedGenerator_H
#define MuonSeedGenerator_CosmicMuonSeedGenerator_H

/** \class CosmicMuonSeedGenerator
 *  SeedGenerator for Cosmic Muon
 *
 *  $Date: 2006/07/15 18:34:24 $
 *  $Revision: 1.3 $
 *  \author Chang Liu - Purdue University 
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"

#include <vector>

namespace edm {class ParameterSet; class Event; class EventSetup;}


class CosmicMuonSeedGenerator: public edm::EDProducer {
 public:

  /// Constructor
  CosmicMuonSeedGenerator(const edm::ParameterSet&);
  
  /// Destructor
  virtual ~CosmicMuonSeedGenerator();
  
  // Operations

  /// reconstruct muon's seeds
  virtual void produce(edm::Event&, const edm::EventSetup&);

 private:

  /// generate TrajectorySeeds and put them into results
  void createSeeds(TrajectorySeedCollection& results,
                   const MuonTransientTrackingRecHit::MuonRecHitContainer& hits,
                   const edm::EventSetup& eSetup) const;

  /// determine if a MuonTransientTrackingRecHit is qualified to build seed
  bool checkQuality(MuonTransientTrackingRecHit::MuonRecHitPointer) const;

  /// create TrajectorySeed from MuonTransientTrackingRecHit 
  std::vector<TrajectorySeed> createSeed(MuonTransientTrackingRecHit::MuonRecHitPointer,
                                         const edm::EventSetup&) const;

 private: 
  /// the name of the DT rec hits collection
  std::string theDTRecSegmentLabel;

  /// the name of the CSC rec hits collection
  std::string theCSCRecSegmentLabel;

  /// the maximum number of Seeds
  unsigned int theMaxSeeds;
  
  /// the maximum chi2 required for dt and csc rechits
  double theMaxDTChi2;
  double theMaxCSCChi2;
 
};
#endif

