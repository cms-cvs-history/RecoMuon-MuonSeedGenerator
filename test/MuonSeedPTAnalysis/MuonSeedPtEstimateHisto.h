#ifndef RecoMuon_MuonSeedPtEstimateHisto_H
#define RecoMuon_MuonSeedPtEstimateHisto_H

/** \class SeedPtEstimateHisto
 *
 * Define histograms and filling methods
 *
 * Author: D. Fortin  - UC Riverside
 */

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TString.h"
#include <string>
#include <iostream>


// Class for CSC segments
class HSeedPt1 {
public:

  /// Constructor from collection name
  HSeedPt1(std::string name_) {
    TString N = name_.c_str();
    name=N;
    
    // Pt Delta Phi 
    hDPhiPt = new TH2F(N+"_hDPhiPt", "Pt DeltaPhi", 60, 0.9, 2.4, 300, -0.25, 1.25);
    hDPsiPt = new TH2F(N+"_hDPsiPt", "Pt DeltaPhi", 60, 0.9, 2.4, 300, 0.50, 2.0);
  }


  /// Constructor from collection name and TFile.
  HSeedPt1(TString name_, TFile* file) {
    name=name_;

    hDPhiPt = (TH2F *) file->Get(name+"_hDPhiPt");
    hDPsiPt = (TH2F *) file->Get(name+"_hDPsiPt");
  }


  virtual ~HSeedPt1() {
    delete hDPhiPt; 
    delete hDPsiPt; 
  }


  /// Fill all the reconstructed histos
  void FillCSC(float simpt, float receta, float dphi, float dpsi) {
    if (fabs(receta) < 3.) {
      hDPhiPt->Fill(fabs(receta), fabs(simpt)*dphi);
      hDPsiPt->Fill(fabs(receta), fabs(simpt)*dpsi);
    }
  }


  void Write() {
    hDPhiPt->Write();
    hDPsiPt->Write();

  }

  TH2F* hDPhiPt;
  TH2F* hDPsiPt;

  TString name;
};


// Class for DT segments
class HSeedPt2 {
public:

  /// Constructor from collection name
  HSeedPt2(std::string name_) {
    TString N = name_.c_str();
    name=N;
    
    // Pt Delta Phi 
    hDPhiPt = new TH2F(N+"_hDPhiPt", "Pt DeltaPhi", 55, 0.0, 1.1, 300, -0.25, 1.25);
    hDPsiPt = new TH2F(N+"_hDPsiPt", "Pt DeltaPsi", 55, 0.0, 1.1, 300, 0.5, 2.0);

  }


  /// Constructor from collection name and TFile.
  HSeedPt2(TString name_, TFile* file) {
    name=name_;

    hDPhiPt = (TH2F *) file->Get(name+"_hDPhiPt");
    hDPsiPt = (TH2F *) file->Get(name+"_hDPsiPt");
  }


  virtual ~HSeedPt2() {
    delete hDPhiPt; 
    delete hDPsiPt; 
  }


  /// Fill all the reconstructed histos
  void FillDT(float simpt, float receta, float dphi, float dpsi) {
    if (fabs(receta) < 3.) {
      hDPhiPt->Fill(fabs(receta), fabs(simpt)*dphi);
      hDPsiPt->Fill(fabs(receta), fabs(simpt)*dpsi);
    }
  }

  void Write() {
    hDPhiPt->Write();
    hDPsiPt->Write();

  }

  TH2F* hDPhiPt;
  TH2F* hDPsiPt;

  TString name;
};



#endif
