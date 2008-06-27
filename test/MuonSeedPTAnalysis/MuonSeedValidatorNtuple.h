#ifndef RecoMuon_MuonSeedValidatorNtuple_H
#define RecoMuon_MuonSeedValidatorNtuple_H

/** \class SeedValidatorNtuple
 *  Collection of Ntuples for SeedGenerator test.
 *
 * Author: S.C. Kao  - UC Riverside
 */

#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include <string>
#include <iostream>

using namespace std;

class TNtuple1 {
public:
 

 TNtuple1() {
    t1 = new TTree("csc_T","eta,phi related information");
    t1->Branch("tcsc",&tcsc.eta_r1,"eta_r1/D:eta_r2:eta_r3:eta_r4:phi_r1:phi_r2:phi_r3:phi_r4:pT"); 
    t1->Branch("tdt" , &tdt.eta_r1,"eta_r1/D:eta_r2:eta_r3:eta_r4:phi_r1:phi_r2:phi_r3:phi_r4:pT"); 
    t1->Branch("tP", &tP.total_p,"total_p/D");
 } 


 /// Destructor
 virtual ~TNtuple1() {
    delete t1;
 }


 void Fill_b1( double eta_r1, double eta_r2, double eta_r3, double eta_r4,
               double phi_r1, double phi_r2, double phi_r3, double phi_r4, double pT ) 
 {
      tcsc.eta_r1 = eta_r1;
      tcsc.eta_r2 = eta_r2;
      tcsc.eta_r3 = eta_r3;
      tcsc.eta_r4 = eta_r4;
      tcsc.phi_r1 = phi_r1;
      tcsc.phi_r2 = phi_r2;
      tcsc.phi_r3 = phi_r3;
      tcsc.phi_r4 = phi_r4;
      tcsc.pT     = pT;
 }
 
 void Fill_b2( double eta_r1, double eta_r2, double eta_r3, double eta_r4,
               double phi_r1, double phi_r2, double phi_r3, double phi_r4, double pT ) 
 {
      tdt.eta_r1 = eta_r1;
      tdt.eta_r2 = eta_r2;
      tdt.eta_r3 = eta_r3;
      tdt.eta_r4 = eta_r4;
      tdt.phi_r1 = phi_r1;
      tdt.phi_r2 = phi_r2;
      tdt.phi_r3 = phi_r3;
      tdt.phi_r4 = phi_r4;
      tdt.pT     = pT;
 }
 
 void Fill_l1(double total_p)
 {
      tP.total_p = total_p;
 }
 
 void FillTree() {
      t1->Fill();
 }

 void Write() {
      t1->Write();
 }

 TTree *t1;

private:

 struct TreeCSC {
        double eta_r1;
        double eta_r2;
        double eta_r3;
        double eta_r4;
        double phi_r1;
        double phi_r2;
        double phi_r3;
        double phi_r4;
        double pT;
 } tcsc;

 struct TreeDT {
        double eta_r1;
        double eta_r2;
        double eta_r3;
        double eta_r4;
        double phi_r1;
        double phi_r2;
        double phi_r3;
        double phi_r4;
        double pT;
 } tdt;

 struct TreeP {
        double total_p;
 } tP;

};

#endif
