#ifndef RecoMuon_MuonSeedGeneratorHistograms_H
#define RecoMuon_MuonSeedGeneratorHistograms_H

/** \class SeedGeneratorHistograms
 *  Collection of histograms for SeedGenerator test.
 *
 * Author: S.C. Kao  - UC Riverside
 */

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TString.h"
#include <string>
#include <iostream>

using namespace std;

class H2DRecHit1 {
public:
 
 H2DRecHit1(std::string name_) {
    TString N1 = name_.c_str();
    name=N1;
    hcsc_dt  = new TH2F(N1+"_hcsc_dt", " CSC Seg vs DT Seg ", 20, -0.25, 9.75, 20, -0.25, 9.75);
    heta_csc = new TH2F(N1+"_heta_csc", " CSC Seg vs eta ", 180, 0.7, 2.5, 20, -0.25, 9.75);
    heta_dt  = new TH2F(N1+"_heta_dt",  " DT  Seg vs eta ", 130, 0.0, 1.3, 20, -0.25, 9.75);
    heta_mu  = new TH2F(N1+"_heta_mu",  " All Seg vs eta ", 300, 0.0, 3.0, 20, -0.25, 9.75);
    heta_mu1 = new TH2F(N1+"_heta_mu1", " All Seg vs eta from track", 300, 0.0, 3.0, 20, -0.25, 9.75);
    heta_mu2 = new TH2F(N1+"_heta_mu2", " All Seg vs eta from simseg", 300, 0.0, 3.0, 20, -0.25, 9.75);

    heta_eta = new TH2F(N1+"_heta_eta", " eta_simhits vs eta_track", 300, 0.0, 3.0, 300, 0.0, 3.0);
    heta_S1L1 = new TH2F(N1+"_heta_S1L1", " eta_S1L1 vs eta_track", 200, 0.0, 2.0, 200, 0.0, 2.0);
    heta_S1L2 = new TH2F(N1+"_heta_S1L2", " eta_S1L2 vs eta_track", 200, 0.0, 2.0, 200, 0.0, 2.0);
    heta_S1L3 = new TH2F(N1+"_heta_S1L3", " eta_S1L3 vs eta_track", 200, 0.0, 2.0, 200, 0.0, 2.0);
    heta_S2L1 = new TH2F(N1+"_heta_S2L1", " eta_S2L1 vs eta_track", 200, 0.0, 2.0, 200, 0.0, 2.0);
    heta_S2L2 = new TH2F(N1+"_heta_S2L2", " eta_S2L2 vs eta_track", 200, 0.0, 2.0, 200, 0.0, 2.0);
    heta_S2L3 = new TH2F(N1+"_heta_S2L3", " eta_S2L3 vs eta_track", 200, 0.0, 2.0, 200, 0.0, 2.0);
    heta_S3L1 = new TH2F(N1+"_heta_S3L1", " eta_S3L1 vs eta_track", 200, 0.0, 2.0, 200, 0.0, 2.0);
    heta_S3L2 = new TH2F(N1+"_heta_S3L2", " eta_S3L2 vs eta_track", 200, 0.0, 2.0, 200, 0.0, 2.0);
    heta_S3L3 = new TH2F(N1+"_heta_S3L3", " eta_S3L3 vs eta_track", 200, 0.0, 2.0, 200, 0.0, 2.0);

    heta_rh  = new TH2F(N1+"_heta_rh",  " All rechits vs eta ", 250, 0.0, 2.5, 20, -0.25, 9.75);
    heta_cscrh = new TH2F(N1+"_heta_cscrh", " csc rechits vs eta ", 250, 0.0, 2.5, 20, -0.25, 9.75);
    heta_dtrh  = new TH2F(N1+"_heta_dtrh",  " dt  rechits vs eta ", 250, 0.0, 2.5, 20, -0.25, 9.75);
    heta_trk = new TH1F(N1+"_heta_trk"," eta from track",300, 0.0, 3.0);

    heta_csc_reco = new TH2F(N1+"_heta_csc_reco", " CSC rec Seg vs rec eta ", 170, 0.8, 2.5, 20, -0.25, 9.75);
    heta_dt_reco = new TH2F(N1+"_heta_dt_reco", " DT rec Seg vs rec eta ", 130, 0.0, 1.3, 20, -0.25, 9.75);
    heta_sim_reco = new TH2F(N1+"_heta_sim_reco", " sim eta vs rec eta ", 130, 0.0, 1.3, 130, 0.0, 1.3);
    heta_sim_reco_csc = new TH2F(N1+"_heta_sim_reco_csc", " sim eta vs rec eta ", 170, 0.8, 2.5, 170, 0.8, 2.5);
    heta_nSeg_DF = new TH2F(N1+"_heta_nSeg_DF",  " eta vs nSeg by DF ", 300, 0.0, 3., 20, -0.25, 9.75);
    heta_nSeg_SC = new TH2F(N1+"_heta_nSeg_SC",  " eta vs nSeg by DF ", 300, 0.0, 3., 20, -0.25, 9.75);
    heta_nSeg_seed = new TH2F(N1+"_heta_nSeg_seed",  " eta vs nSeg from rec seed", 300, 0.0, 3., 20, -0.25, 9.75);

    hcsc_q  = new TH2F(N1+"_hcsc_q", " CSC Q from vec - Q from dphi ", 29, -7.25, 7.25, 29, -7.25, 7.25);
    hcsc_qr = new TH2F(N1+"_hcsc_qr", " CSC Q from vec - Q from dphi R", 29, -7.25, 7.25, 29, -7.25, 7.25);
    hdt_q   = new TH2F(N1+"_hdt_q", " DT Q from vec - Q from dphi ", 29, -7.25, 7.25, 29, -7.25, 7.25);
    hdt_qr  = new TH2F(N1+"_hdt_qr", " DT Q from vec - Q from dphi R", 29, -7.25, 7.25, 29, -7.25, 7.25);
 } 

 H2DRecHit1(TString name_, TFile* file) {
    name=name_;
    hcsc_dt  = (TH2F *) file->Get(name+"_hcsc_dt");
    heta_csc = (TH2F *) file->Get(name+"_heta_csc");
    heta_dt  = (TH2F *) file->Get(name+"_heta_dt");
    heta_mu  = (TH2F *) file->Get(name+"_heta_mu");
    heta_mu1 = (TH2F *) file->Get(name+"_heta_mu1");
    heta_mu2 = (TH2F *) file->Get(name+"_heta_mu2");

    heta_eta = (TH2F *) file->Get(name+"_heta_eta");
    heta_S1L1= (TH2F *) file->Get(name+"_heta_S1L1");
    heta_S1L2= (TH2F *) file->Get(name+"_heta_S1L2");
    heta_S1L3= (TH2F *) file->Get(name+"_heta_S1L3");
    heta_S2L1= (TH2F *) file->Get(name+"_heta_S2L1");
    heta_S2L2= (TH2F *) file->Get(name+"_heta_S2L2");
    heta_S2L3= (TH2F *) file->Get(name+"_heta_S2L3");
    heta_S3L1= (TH2F *) file->Get(name+"_heta_S3L1");
    heta_S3L2= (TH2F *) file->Get(name+"_heta_S3L2");
    heta_S3L3= (TH2F *) file->Get(name+"_heta_S3L3");

    heta_rh  = (TH2F *) file->Get(name+"_heta_rh");
    heta_cscrh = (TH2F *) file->Get(name+"_heta_cscrh");
    heta_dtrh  = (TH2F *) file->Get(name+"_heta_dtrh");
    heta_trk = (TH1F *) file->Get(name+"_heta_trk");

    heta_csc_reco = (TH2F *) file->Get(name+"_heta_csc_reco");
    heta_dt_reco  = (TH2F *) file->Get(name+"_heta_dt_reco");
    heta_sim_reco = (TH2F *) file->Get(name+"_heta_sim_reco");
    heta_sim_reco_csc = (TH2F *) file->Get(name+"_heta_sim_reco_csc");
    heta_nSeg_DF = (TH2F *) file->Get(name+"_heta_nSeg_DF"); 
    heta_nSeg_SC = (TH2F *) file->Get(name+"_heta_nSeg_SC"); 
    heta_nSeg_seed = (TH2F *) file->Get(name+"_heta_nSeg_seed"); 

    hcsc_q  = (TH2F *) file->Get(name+"_hcsc_q");
    hcsc_qr = (TH2F *) file->Get(name+"_hcsc_qr");
    hdt_q   = (TH2F *) file->Get(name+"_hdt_q");
    hdt_qr  = (TH2F *) file->Get(name+"_hdt_qr");
 }

 /// Destructor
 virtual ~H2DRecHit1() {
    delete hcsc_dt;
    delete heta_csc;
    delete heta_dt;
    delete heta_mu;
    delete heta_mu1;
    delete heta_mu2;

    delete heta_eta;
    delete heta_S1L1;
    delete heta_S1L2;
    delete heta_S1L3;
    delete heta_S2L1;
    delete heta_S2L2;
    delete heta_S2L3;
    delete heta_S3L1;
    delete heta_S3L2;
    delete heta_S3L3;

    delete heta_rh;
    delete heta_cscrh;
    delete heta_dtrh;
    delete heta_trk;

    delete heta_csc_reco;
    delete heta_dt_reco;
    delete heta_sim_reco;
    delete heta_sim_reco_csc;
    delete heta_nSeg_DF;
    delete heta_nSeg_SC;
    delete heta_nSeg_seed;
  
    delete hcsc_q;
    delete hcsc_qr;
    delete hdt_q;
    delete hdt_qr;
 }

 void Fill1(int csc_nu,int dt_nu,double eta_c, double eta_d, int all_nu, double eta_a, double eta_trk) {
      hcsc_dt->Fill(csc_nu,dt_nu);
      heta_csc->Fill(eta_c,csc_nu);
      heta_dt->Fill(eta_d,dt_nu);
      heta_mu->Fill(eta_a,all_nu);
      heta_mu1->Fill(eta_trk,all_nu);
      heta_trk->Fill(eta_trk);
      heta_eta->Fill(eta_a,eta_trk);
 }
 void Fill1a(int rh_nu, double eta_a, int cscrh_nu, int dtrh_nu) {
      heta_rh->Fill(eta_a,rh_nu);
      heta_cscrh->Fill(eta_a,cscrh_nu);
      heta_dtrh->Fill(eta_a,dtrh_nu);
 }
 void Fill1b(int sim_nu, double eta_sim ) {
      heta_mu2->Fill(eta_sim,sim_nu);
 }
 void Fill1c(double h_trk, double h11,double h12,double h13,double h21,double h22,double h23,
                           double h31,double h32,double h33) {
      heta_S1L1->Fill(h_trk,h11);
      heta_S1L2->Fill(h_trk,h12);
      heta_S1L3->Fill(h_trk,h13);
      heta_S2L1->Fill(h_trk,h21);
      heta_S2L2->Fill(h_trk,h22);
      heta_S2L3->Fill(h_trk,h23);
      heta_S3L1->Fill(h_trk,h31);
      heta_S3L2->Fill(h_trk,h32);
      heta_S3L3->Fill(h_trk,h33);
 }
 void Fill1d(int nu_rec_seg, double eta_rec, double eta_sim ) {
      heta_dt_reco->Fill(eta_rec,nu_rec_seg);
      heta_sim_reco->Fill(eta_sim,eta_rec);
 }
 void Fill1f(int nu_rec_seg, double eta_rec, double eta_sim ) {
      heta_csc_reco->Fill(eta_rec,nu_rec_seg);
      heta_sim_reco_csc->Fill(eta_sim,eta_rec);
 }
 void Fill1e(int nSeg_DF, float eta_DF) {
      heta_nSeg_DF->Fill(eta_DF,nSeg_DF);
 }
 void Fill1g(int nSeg_SC, float eta_SC) {
      heta_nSeg_SC->Fill(eta_SC,nSeg_SC);
 }
 void Fill1h(int nSeg_seed, float eta_seed) {
      heta_nSeg_seed->Fill(eta_seed,nSeg_seed);
 }
 
 void Fill2(int q_v, int q_f) {
      hcsc_q->Fill(q_v,q_f);
 }
 void Fill2a(int q_v, int q_f) {
      hcsc_qr->Fill(q_v,q_f);
 }
 void Fill2b(int q_v, int q_f) {
      hdt_q->Fill(q_v,q_f);
 }
 void Fill2c(int q_v, int q_f) {
      hdt_qr->Fill(q_v,q_f);
 }

 void Write() {
      hcsc_dt->Write();
      heta_csc->Write();
      heta_dt->Write();
      heta_mu->Write();
      heta_mu1->Write();
      heta_mu2->Write();

      heta_eta->Write();
      heta_S1L1->Write();
      heta_S1L2->Write();
      heta_S1L3->Write();
      heta_S2L1->Write();
      heta_S2L2->Write();
      heta_S2L3->Write();
      heta_S3L1->Write();
      heta_S3L2->Write();
      heta_S3L3->Write();

      heta_rh->Write();
      heta_cscrh->Write();
      heta_dtrh->Write();
      heta_trk->Write();

      heta_csc_reco->Write();
      heta_dt_reco->Write();
      heta_sim_reco->Write();
      heta_sim_reco_csc->Write();
      heta_nSeg_DF->Write();
      heta_nSeg_SC->Write();
      heta_nSeg_seed->Write();

      hcsc_q->Write();
      hcsc_qr->Write();
      hdt_q->Write();
      hdt_qr->Write();
 }

 TH2F *hcsc_dt;
 TH2F *heta_csc;
 TH2F *heta_dt;
 TH2F *heta_mu;
 TH2F *heta_mu1;
 TH2F *heta_mu2;

 TH2F *heta_eta;
 TH2F *heta_S1L1;
 TH2F *heta_S1L2;
 TH2F *heta_S1L3;
 TH2F *heta_S2L1;
 TH2F *heta_S2L2;
 TH2F *heta_S2L3;
 TH2F *heta_S3L1;
 TH2F *heta_S3L2;
 TH2F *heta_S3L3;

 TH2F *heta_rh;
 TH2F *heta_cscrh;
 TH2F *heta_dtrh;
 TH1F *heta_trk;

 TH2F *heta_csc_reco;
 TH2F *heta_dt_reco;
 TH2F *heta_sim_reco;
 TH2F *heta_sim_reco_csc;
 TH2F *heta_nSeg_DF;
 TH2F *heta_nSeg_SC;
 TH2F *heta_nSeg_seed;

 TH2F *hcsc_q;
 TH2F *hcsc_qr;
 TH2F *hdt_q;
 TH2F *hdt_qr;

 TString name;

};


class H2DRecHit2 {
public:

 H2DRecHit2(std::string name_) {
    TString N2 = name_.c_str();
    name=N2;

    hPt = new TH1F(N2+"_hPt", " Pt of Tracks ", 100,5.,205.);
    hPa_Pt     = new TH2F(N2+"_hPa_Pt", "P vs Pt", 50, 5., 205., 100, 0., 1000.);
    hPt_CSCSeg = new TH2F(N2+"_hPt_CSCSeg", "Pt vs total CSC segments number", 40, 5., 205., 40, -0.25, 19.75);
    hst_csc_nu = new TH2F(N2+"_hst_csc_nu", "station vs csc seg number", 6, 0, 6, 20, -0.25, 9.75);
    hst_csc_nu1= new TH2F(N2+"_hst_csc_nu1","station vs csc seg number w/ > 5 hits seg", 6, 0, 6, 20, -0.25, 9.75);

    heta_dphi12 = new TH2F(N2+"_heta_dphi12", "eta vs dphi12", 150, 1., 2.5, 200, -0.03, 0.03);
    heta_dphi13 = new TH2F(N2+"_heta_dphi13", "eta vs dphi13", 150, 1., 2.5, 200, -0.03, 0.03);
    heta_dphi14 = new TH2F(N2+"_heta_dphi14", "eta vs dphi14", 150, 1., 2.5, 200, -0.03, 0.03);
    heta_dphi23 = new TH2F(N2+"_heta_dphi23", "eta vs dphi23", 150, 1., 2.5, 200, -0.03, 0.03);
    heta_dphi24 = new TH2F(N2+"_heta_dphi24", "eta vs dphi24", 150, 1., 2.5, 200, -0.03, 0.03);
    heta_dphi34 = new TH2F(N2+"_heta_dphi34", "eta vs dphi34", 150, 1., 2.5, 200, -0.03, 0.03);

    heta_rdphi12 = new TH2F(N2+"_heta_rdphi12", "eta vs rdphi12", 150, 1., 2.5, 200, -0.03, 0.03);
    heta_rdphi13 = new TH2F(N2+"_heta_rdphi13", "eta vs rdphi13", 150, 1., 2.5, 200, -0.03, 0.03);
    heta_rdphi14 = new TH2F(N2+"_heta_rdphi14", "eta vs rdphi14", 150, 1., 2.5, 200, -0.03, 0.03);
    heta_rdphi23 = new TH2F(N2+"_heta_rdphi23", "eta vs rdphi23", 150, 1., 2.5, 200, -0.03, 0.03);
    heta_rdphi24 = new TH2F(N2+"_heta_rdphi24", "eta vs rdphi24", 150, 1., 2.5, 200, -0.03, 0.03);
    heta_rdphi34 = new TH2F(N2+"_heta_rdphi34", "eta vs rdphi34", 150, 1., 2.5, 200, -0.03, 0.03);

    heta_dphiPt12 = new TH2F(N2+"_heta_dphiPt12", "eta vs dphi*Pt12", 150, 1., 2.5, 200, -1., 1.);
    heta_dphiPt13 = new TH2F(N2+"_heta_dphiPt13", "eta vs dphi*Pt13", 150, 1., 2.5, 200, -1., 1.);
    heta_dphiPt14 = new TH2F(N2+"_heta_dphiPt14", "eta vs dphi*Pt14", 150, 1., 2.5, 200, -1., 1.);
    heta_dphiPt23 = new TH2F(N2+"_heta_dphiPt23", "eta vs dphi*Pt23", 150, 1., 2.5, 200, -0.25, 0.25);
    heta_dphiPt24 = new TH2F(N2+"_heta_dphiPt24", "eta vs dphi*Pt24", 150, 1., 2.5, 200, -0.25, 0.25);
    heta_dphiPt34 = new TH2F(N2+"_heta_dphiPt34", "eta vs dphi*Pt34", 150, 1., 2.5, 200, -0.25, 0.25);

    heta_rdphiPt12 = new TH2F(N2+"_heta_rdphiPt12", "eta vs rdphi*Pt12", 150, 1., 2.5, 200, -1., 1.);
    heta_rdphiPt13 = new TH2F(N2+"_heta_rdphiPt13", "eta vs rdphi*Pt13", 150, 1., 2.5, 200, -1., 1.);
    heta_rdphiPt14 = new TH2F(N2+"_heta_rdphiPt14", "eta vs rdphi*Pt14", 150, 1., 2.5, 200, -1., 1.);
    heta_rdphiPt23 = new TH2F(N2+"_heta_rdphiPt23", "eta vs rdphi*Pt23", 150, 1., 2.5, 200, -0.25, 0.25);
    heta_rdphiPt24 = new TH2F(N2+"_heta_rdphiPt24", "eta vs rdphi*Pt24", 150, 1., 2.5, 200, -0.25, 0.25);
    heta_rdphiPt34 = new TH2F(N2+"_heta_rdphiPt34", "eta vs rdphi*Pt34", 150, 1., 2.5, 200, -0.25, 0.25);

    heta_dphiRatio = new TH2F(N2+"_heta_dphiRatio", "eta vs dphi Ratio", 150, 1., 2.5, 500, -0.5, 1.5);
    heta_rdphiRatio = new TH2F(N2+"_heta_rdphiRatio", "eta vs rdphi Ratio", 150, 1., 2.5, 500, -0.5, 1.5);
    heta_rdphiRatio2 = new TH2F(N2+"_heta_rdphiRatio2", "eta vs rdphi Ratio2", 150, 1., 2.5, 500, -0.5, 1.5);
    hpt_dphiRatio = new TH2F(N2+"_hpt_dphiRatio", "pt vs dphi Ratio", 2000, 0., 200.0, 500, -0.5, 1.5);
    hpt_rdphiRatio = new TH2F(N2+"_hpt_rdphiRatio", "pt vs rdphi Ratio", 200, 0., 200.0, 500, -0.5, 1.5);
    hpt_rdphiRatio2 = new TH2F(N2+"_hpt_rdphiRatio2", "pt vs rdphi Ratio2", 200, 0., 200.0, 500, -0.5, 1.5);

    heta_dphiRatiof = new TH2F(N2+"_heta_dphiRatiof", "eta vs dphi Ratio", 150, 1., 2.5, 500, -0.5, 1.5);
    heta_rdphiRatiof = new TH2F(N2+"_heta_rdphiRatiof", "eta vs rdphi Ratio", 150, 1., 2.5, 500, -0.5, 1.5);
    heta_rdphiRatio2f = new TH2F(N2+"_heta_rdphiRatio2f", "eta vs rdphi Ratio2", 150, 1., 2.5, 500, -0.5, 1.5);
    hpt_dphiRatiof = new TH2F(N2+"_hpt_dphiRatiof", "pt vs dphi Ratio", 2000, 0., 200.0, 500, -0.5, 1.5);
    hpt_rdphiRatiof = new TH2F(N2+"_hpt_rdphiRatiof", "pt vs rdphi Ratio", 200, 0., 200.0, 500, -0.5, 1.5);
    hpt_rdphiRatio2f = new TH2F(N2+"_hpt_rdphiRatio2f", "pt vs rdphi Ratio2", 200, 0., 200.0, 500, -0.5, 1.5);
  }

 H2DRecHit2(TString name_, TFile* file) {
    name=name_;

    hPt = (TH1F *) file->Get(name+"_hPt");
    hPa_Pt     = (TH2F *) file->Get(name+"_hPa_Pt");
    hPt_CSCSeg = (TH2F *) file->Get(name+"_hPt_CSCSeg");
    hst_csc_nu = (TH2F *) file->Get(name+"_hst_csc_nu");
    hst_csc_nu1= (TH2F *) file->Get(name+"_hst_csc_nu1");

    heta_dphi12 = (TH2F *) file->Get(name+"_heta_dphi12");
    heta_dphi13 = (TH2F *) file->Get(name+"_heta_dphi13");
    heta_dphi14 = (TH2F *) file->Get(name+"_heta_dphi14");
    heta_dphi23 = (TH2F *) file->Get(name+"_heta_dphi23");
    heta_dphi24 = (TH2F *) file->Get(name+"_heta_dphi24");
    heta_dphi34 = (TH2F *) file->Get(name+"_heta_dphi34");

    heta_rdphi12 = (TH2F *) file->Get(name+"_heta_rdphi12");
    heta_rdphi13 = (TH2F *) file->Get(name+"_heta_rdphi13");
    heta_rdphi14 = (TH2F *) file->Get(name+"_heta_rdphi14");
    heta_rdphi23 = (TH2F *) file->Get(name+"_heta_rdphi23");
    heta_rdphi24 = (TH2F *) file->Get(name+"_heta_rdphi24");
    heta_rdphi34 = (TH2F *) file->Get(name+"_heta_rdphi34");
   
    heta_dphiPt12 = (TH2F *) file->Get(name+"_heta_dphiPt12");
    heta_dphiPt13 = (TH2F *) file->Get(name+"_heta_dphiPt13");
    heta_dphiPt14 = (TH2F *) file->Get(name+"_heta_dphiPt14");
    heta_dphiPt23 = (TH2F *) file->Get(name+"_heta_dphiPt23");
    heta_dphiPt24 = (TH2F *) file->Get(name+"_heta_dphiPt24");
    heta_dphiPt34 = (TH2F *) file->Get(name+"_heta_dphiPt34");

    heta_rdphiPt12 = (TH2F *) file->Get(name+"_heta_rdphiPt12");
    heta_rdphiPt13 = (TH2F *) file->Get(name+"_heta_rdphiPt13");
    heta_rdphiPt14 = (TH2F *) file->Get(name+"_heta_rdphiPt14");
    heta_rdphiPt23 = (TH2F *) file->Get(name+"_heta_rdphiPt23");
    heta_rdphiPt24 = (TH2F *) file->Get(name+"_heta_rdphiPt24");
    heta_rdphiPt34 = (TH2F *) file->Get(name+"_heta_rdphiPt34");

    heta_dphiRatio = (TH2F *) file->Get(name+"_heta_dphiRatio");
    heta_rdphiRatio = (TH2F *) file->Get(name+"_heta_rdphiRatio");
    heta_rdphiRatio2 = (TH2F *) file->Get(name+"_heta_rdphiRatio2");
    hpt_dphiRatio = (TH2F *) file->Get(name+"_hpt_dphiRatio");
    hpt_rdphiRatio = (TH2F *) file->Get(name+"_hpt_rdphiRatio");
    hpt_rdphiRatio2 = (TH2F *) file->Get(name+"_hpt_rdphiRatio2");

    heta_dphiRatiof = (TH2F *) file->Get(name+"_heta_dphiRatiof");
    heta_rdphiRatiof = (TH2F *) file->Get(name+"_heta_rdphiRatiof");
    heta_rdphiRatio2f = (TH2F *) file->Get(name+"_heta_rdphiRatio2f");
    hpt_dphiRatiof = (TH2F *) file->Get(name+"_hpt_dphiRatiof");
    hpt_rdphiRatiof = (TH2F *) file->Get(name+"_hpt_rdphiRatiof");
    hpt_rdphiRatio2f = (TH2F *) file->Get(name+"_hpt_rdphiRatio2f");
  } 

  /// Destructor
  virtual ~H2DRecHit2() {

    delete hPt;
    delete hPa_Pt;
    delete hPt_CSCSeg;
    delete hst_csc_nu;
    delete hst_csc_nu1;

    delete heta_dphi12;
    delete heta_dphi13;
    delete heta_dphi14;
    delete heta_dphi23;
    delete heta_dphi24;
    delete heta_dphi34;

    delete heta_rdphi12;
    delete heta_rdphi13;
    delete heta_rdphi14;
    delete heta_rdphi23;
    delete heta_rdphi24;
    delete heta_rdphi34;

    delete heta_dphiPt12;
    delete heta_dphiPt13;
    delete heta_dphiPt14;
    delete heta_dphiPt23;
    delete heta_dphiPt24;
    delete heta_dphiPt34;

    delete heta_rdphiPt12;
    delete heta_rdphiPt13;
    delete heta_rdphiPt14;
    delete heta_rdphiPt23;
    delete heta_rdphiPt24;
    delete heta_rdphiPt34;
 
    delete heta_dphiRatio;
    delete heta_rdphiRatio;
    delete heta_rdphiRatio2;
    delete hpt_dphiRatio;
    delete hpt_rdphiRatio;
    delete hpt_rdphiRatio2;

    delete heta_dphiRatiof;
    delete heta_rdphiRatiof;
    delete heta_rdphiRatio2f;
    delete hpt_dphiRatiof;
    delete hpt_rdphiRatiof;
    delete hpt_rdphiRatio2f;
  }

  void Fill3(float Pt, float Pa, int csc_nu)
  {
       hPt->Fill(Pt);
       hPa_Pt->Fill(Pt,Pa);
       hPt_CSCSeg->Fill(Pt,csc_nu);
  }

  void Fill4(int k, int csc_nu, int csc_nu1 )
  {
       hst_csc_nu->Fill(k,csc_nu);
       hst_csc_nu1->Fill(k,csc_nu1);
  }
  /*
  void Fill5(double dphi12, double dphi13, double dphi14,
             double dphi23, double dphi24, double dphi34,
             double eta2, double eta3, double eta4 )
  {
     heta_dphi12->Fill(eta2,dphi12);
     heta_dphi13->Fill(eta3,dphi13);
     heta_dphi14->Fill(eta4,dphi14);
     heta_dphi23->Fill(eta3,dphi23);
     heta_dphi24->Fill(eta4,dphi24);
     heta_dphi34->Fill(eta4,dphi34);

  }
  void Fill5a(double dphi12, double dphi13, double dphi14,
              double dphi23, double dphi24, double dphi34,
              double eta2, double eta3, double eta4 )
  {
     heta_rdphi12->Fill(eta2,dphi12);
     heta_rdphi13->Fill(eta3,dphi13);
     heta_rdphi14->Fill(eta4,dphi14);
     heta_rdphi23->Fill(eta3,dphi23);
     heta_rdphi24->Fill(eta4,dphi24);
     heta_rdphi34->Fill(eta4,dphi34);
  }*/

  void Fill5_1(double dphi,double dphir,double ptxdphi,double ptxdphir,double eta,double eta_r)
  {
    heta_dphi12->Fill(eta,dphi);
    heta_rdphi12->Fill(eta_r,dphir);
    heta_dphiPt12->Fill(eta,ptxdphi);
    heta_rdphiPt12->Fill(eta_r,ptxdphir);
  }
  void Fill5_2(double dphi,double dphir,double ptxdphi,double ptxdphir,double eta,double eta_r)
  {
    heta_dphi13->Fill(eta,dphi);
    heta_rdphi13->Fill(eta_r,dphir);
    heta_dphiPt13->Fill(eta,ptxdphi);
    heta_rdphiPt13->Fill(eta_r,ptxdphir);
  }
  void Fill5_3(double dphi,double dphir,double ptxdphi,double ptxdphir,double eta,double eta_r)
  {
    heta_dphi14->Fill(eta,dphi);
    heta_rdphi14->Fill(eta_r,dphir);
    heta_dphiPt14->Fill(eta,ptxdphi);
    heta_rdphiPt14->Fill(eta_r,ptxdphir);
  }
  void Fill5_4(double dphi,double dphir,double ptxdphi,double ptxdphir,double eta,double eta_r)
  {
    heta_dphi23->Fill(eta,dphi);
    heta_rdphi23->Fill(eta_r,dphir);
    heta_dphiPt23->Fill(eta,ptxdphi);
    heta_rdphiPt23->Fill(eta_r,ptxdphir);
  }
  void Fill5_5(double dphi,double dphir,double ptxdphi,double ptxdphir,double eta,double eta_r)
  {
    heta_dphi24->Fill(eta,dphi);
    heta_rdphi24->Fill(eta_r,dphir);
    heta_dphiPt24->Fill(eta,ptxdphi);
    heta_rdphiPt24->Fill(eta_r,ptxdphir);
  }
  void Fill5_6(double dphi,double dphir,double ptxdphi,double ptxdphir,double eta,double eta_r)
  {
    heta_dphi34->Fill(eta,dphi);
    heta_rdphi34->Fill(eta_r,dphir);
    heta_dphiPt34->Fill(eta,ptxdphi);
    heta_rdphiPt34->Fill(eta_r,ptxdphir);
  }

  void Fill5b(double eta2, double dphi_R,double reta2, double rdphi_R, double pt )
  {
     heta_dphiRatio->Fill(eta2,dphi_R);
     heta_rdphiRatio->Fill(reta2,rdphi_R);

     hpt_dphiRatio->Fill(pt,dphi_R);
     hpt_rdphiRatio->Fill(pt,rdphi_R);
  }

  void Fill5b1(double eta2, double dphi_R,double reta2, double rdphi_R, double pt )
  {
     heta_dphiRatiof->Fill(eta2,dphi_R);
     heta_rdphiRatiof->Fill(reta2,rdphi_R);

     hpt_dphiRatiof->Fill(pt,dphi_R);
     hpt_rdphiRatiof->Fill(pt,rdphi_R);
  }

  void Fill5c(double reta2, double rdphi_R, double pt )
  {
     heta_rdphiRatio2->Fill(reta2,rdphi_R);
     hpt_rdphiRatio2->Fill(pt,rdphi_R);
  }

  void Fill5c1(double reta2, double rdphi_R, double pt )
  {
     heta_rdphiRatio2f->Fill(reta2,rdphi_R);
     hpt_rdphiRatio2f->Fill(pt,rdphi_R);
  }

  /*
  void Fill6(double ptxdphi12, double ptxdphi13, double ptxdphi14,
             double ptxdphi23, double ptxdphi24, double ptxdphi34,
             double eta2, double eta3, double eta4 )
  {
     heta_dphiPt12->Fill(eta2,ptxdphi12);
     heta_dphiPt13->Fill(eta3,ptxdphi13);
     heta_dphiPt14->Fill(eta4,ptxdphi14);
     heta_dphiPt23->Fill(eta3,ptxdphi23);
     heta_dphiPt24->Fill(eta4,ptxdphi24);
     heta_dphiPt34->Fill(eta4,ptxdphi34);
  }
  void Fill6a(double ptxdphi12, double ptxdphi13, double ptxdphi14,
              double ptxdphi23, double ptxdphi24, double ptxdphi34,
              double eta2, double eta3, double eta4 )
  {
     heta_rdphiPt12->Fill(eta2,ptxdphi12);
     heta_rdphiPt13->Fill(eta3,ptxdphi13);
     heta_rdphiPt14->Fill(eta4,ptxdphi14);
     heta_rdphiPt23->Fill(eta3,ptxdphi23);
     heta_rdphiPt24->Fill(eta4,ptxdphi24);
     heta_rdphiPt34->Fill(eta4,ptxdphi34);
  }*/

  void Write() {

       hPt->Write();
       hPa_Pt->Write();
       hPt_CSCSeg->Write();
       hst_csc_nu->Write();
       hst_csc_nu1->Write();
 
       heta_dphi12->Write();
       heta_dphi13->Write();
       heta_dphi14->Write();
       heta_dphi23->Write();
       heta_dphi24->Write();
       heta_dphi34->Write();

       heta_rdphi12->Write();
       heta_rdphi13->Write();
       heta_rdphi14->Write();
       heta_rdphi23->Write();
       heta_rdphi24->Write();
       heta_rdphi34->Write();

       heta_dphiPt12->Write();
       heta_dphiPt13->Write();
       heta_dphiPt14->Write();
       heta_dphiPt23->Write();
       heta_dphiPt24->Write();
       heta_dphiPt34->Write();

       heta_rdphiPt12->Write();
       heta_rdphiPt13->Write();
       heta_rdphiPt14->Write();
       heta_rdphiPt23->Write();
       heta_rdphiPt24->Write();
       heta_rdphiPt34->Write();
  
       heta_dphiRatio->Write();
       heta_rdphiRatio->Write();
       heta_rdphiRatio2->Write();
       hpt_dphiRatio->Write();
       hpt_rdphiRatio->Write();
       hpt_rdphiRatio2->Write();

       heta_dphiRatiof->Write();
       heta_rdphiRatiof->Write();
       heta_rdphiRatio2f->Write();
       hpt_dphiRatiof->Write();
       hpt_rdphiRatiof->Write();
       hpt_rdphiRatio2f->Write();
  }

  TH1F *hPt;
  TH2F *hPa_Pt;
  TH2F *hPt_CSCSeg;
  TH2F *hst_csc_nu;
  TH2F *hst_csc_nu1;

  TH2F *heta_dphi12;
  TH2F *heta_dphi13;
  TH2F *heta_dphi14;
  TH2F *heta_dphi23;
  TH2F *heta_dphi24;
  TH2F *heta_dphi34;

  TH2F *heta_rdphi12;
  TH2F *heta_rdphi13;
  TH2F *heta_rdphi14;
  TH2F *heta_rdphi23;
  TH2F *heta_rdphi24;
  TH2F *heta_rdphi34;

  TH2F *heta_dphiPt12;
  TH2F *heta_dphiPt13;
  TH2F *heta_dphiPt14;
  TH2F *heta_dphiPt23;
  TH2F *heta_dphiPt24;
  TH2F *heta_dphiPt34;

  TH2F *heta_rdphiPt12;
  TH2F *heta_rdphiPt13;
  TH2F *heta_rdphiPt14;
  TH2F *heta_rdphiPt23;
  TH2F *heta_rdphiPt24;
  TH2F *heta_rdphiPt34;

  TH2F *heta_dphiRatio;
  TH2F *heta_rdphiRatio;
  TH2F *heta_rdphiRatio2;
  TH2F *hpt_dphiRatio;
  TH2F *hpt_rdphiRatio;
  TH2F *hpt_rdphiRatio2;

  TH2F *heta_dphiRatiof;
  TH2F *heta_rdphiRatiof;
  TH2F *heta_rdphiRatio2f;
  TH2F *hpt_dphiRatiof;
  TH2F *hpt_rdphiRatiof;
  TH2F *hpt_rdphiRatio2f;

  TString name;
};


class H2DRecHit3 {
public:

 H2DRecHit3(std::string name_) {
    TString N3 = name_.c_str();
    name=N3;

    hPt = new TH1F(N3+"_hPt", " Pt of Tracks ", 100,5.,205.);
    hPa_Pt     = new TH2F(N3+"_hPa_Pt", "P vs Pt", 50, 5., 205., 100, 0., 1000.);
    hPt_DTSeg  = new TH2F(N3+"_hPt_DTSeg",  "Pt vs total DT segments number", 40, 5., 205., 40, -0.25, 19.75);
    hst_dt_nu  = new TH2F(N3+"_hst_dt_nu", "station vs dt seg number", 6, 0, 6, 20, -0.25, 9.75);
    hst_dt_nu1 = new TH2F(N3+"_hst_dt_nu1","station vs dt seg number w/ > 5 hits seg", 6, 0, 6, 20, -0.25, 9.75);

    hrdphi_Pt   = new TH2F(N3+"_hrdphi_Pt", "rec dphi vs Pt", 50, 5., 205., 200, -0.06, 0.06);
    hrdeta_Pt   = new TH2F(N3+"_hrdeta_Pt", "rec deta vs Pt", 50, 5., 205., 202, -0.01, 0.01);

    hbeta_dphi12 = new TH2F(N3+"_hbeta_dphi12", "eta vs dphi12", 110, 0., 1.1, 400, -0.01, 0.01);
    hbeta_dphi13 = new TH2F(N3+"_hbeta_dphi13", "eta vs dphi13", 110, 0., 1.1, 400, -0.01, 0.01);
    hbeta_dphi14 = new TH2F(N3+"_hbeta_dphi14", "eta vs dphi14", 110, 0., 1.1, 400, -0.01, 0.01);
    hbeta_dphi23 = new TH2F(N3+"_hbeta_dphi23", "eta vs dphi23", 110, 0., 1.1, 400, -0.01, 0.01);
    hbeta_dphi24 = new TH2F(N3+"_hbeta_dphi24", "eta vs dphi24", 110, 0., 1.1, 400, -0.01, 0.01);
    hbeta_dphi34 = new TH2F(N3+"_hbeta_dphi34", "eta vs dphi34", 110, 0., 1.1, 400, -0.01, 0.01);

    hbeta_rdphi12 = new TH2F(N3+"_hbeta_rdphi12", "eta vs rdphi12", 110, 0., 1.1, 400, -0.01, 0.01);
    hbeta_rdphi13 = new TH2F(N3+"_hbeta_rdphi13", "eta vs rdphi13", 110, 0., 1.1, 400, -0.01, 0.01);
    hbeta_rdphi14 = new TH2F(N3+"_hbeta_rdphi14", "eta vs rdphi14", 110, 0., 1.1, 400, -0.01, 0.01);
    hbeta_rdphi23 = new TH2F(N3+"_hbeta_rdphi23", "eta vs rdphi23", 110, 0., 1.1, 400, -0.01, 0.01);
    hbeta_rdphi24 = new TH2F(N3+"_hbeta_rdphi24", "eta vs rdphi24", 110, 0., 1.1, 400, -0.01, 0.01);
    hbeta_rdphi34 = new TH2F(N3+"_hbeta_rdphi34", "eta vs rdphi34", 110, 0., 1.1, 400, -0.01, 0.01);

    hbeta_dphiPt12 = new TH2F(N3+"_hbeta_dphiPt12", "eta vs dphi*Pt12", 110, 0., 1.1, 200, -1., 1.);
    hbeta_dphiPt13 = new TH2F(N3+"_hbeta_dphiPt13", "eta vs dphi*Pt13", 110, 0., 1.1, 200, -1., 1.);
    hbeta_dphiPt14 = new TH2F(N3+"_hbeta_dphiPt14", "eta vs dphi*Pt14", 110, 0., 1.1, 200, -1., 1.);
    hbeta_dphiPt23 = new TH2F(N3+"_hbeta_dphiPt23", "eta vs dphi*Pt23", 110, 0., 1.1, 200, -1., 1.);
    hbeta_dphiPt24 = new TH2F(N3+"_hbeta_dphiPt24", "eta vs dphi*Pt24", 110, 0., 1.1, 200, -1., 1.);
    hbeta_dphiPt34 = new TH2F(N3+"_hbeta_dphiPt34", "eta vs dphi*Pt34", 110, 0., 1.1, 200, -1., 1.);

    hbeta_rdphiPt12 = new TH2F(N3+"_hbeta_rdphiPt12", "eta vs rdphi*Pt12", 110, 0., 1.1, 200, -1., 1.);
    hbeta_rdphiPt13 = new TH2F(N3+"_hbeta_rdphiPt13", "eta vs rdphi*Pt13", 110, 0., 1.1, 200, -1., 1.);
    hbeta_rdphiPt14 = new TH2F(N3+"_hbeta_rdphiPt14", "eta vs rdphi*Pt14", 110, 0., 1.1, 200, -1., 1.);
    hbeta_rdphiPt23 = new TH2F(N3+"_hbeta_rdphiPt23", "eta vs rdphi*Pt23", 110, 0., 1.1, 200, -1., 1.);
    hbeta_rdphiPt24 = new TH2F(N3+"_hbeta_rdphiPt24", "eta vs rdphi*Pt24", 110, 0., 1.1, 200, -1., 1.);
    hbeta_rdphiPt34 = new TH2F(N3+"_hbeta_rdphiPt34", "eta vs rdphi*Pt34", 110, 0., 1.1, 200, -1., 1.);

    hbeta_dphiRatio = new TH2F(N3+"_hbeta_dphiRatio", "eta vs dphi Ratio", 110, 0., 1.1, 500, -0.5, 1.5);
    hbeta_rdphiRatio = new TH2F(N3+"_hbeta_rdphiRatio", "eta vs rdphi Ratio", 110, 0., 1.1, 500, -0.5, 1.5);
    hbeta_rdphiRatio2 = new TH2F(N3+"_hbeta_rdphiRatio2", "eta vs rdphi Ratio2", 110, 0., 1.1, 500, -0.5, 1.5);
    hpt_dphiRatio = new TH2F(N3+"_hpt_dphiRatio", "pt vs dphi Ratio", 200, 0., 200.0, 500, -0.5, 1.5);
    hpt_rdphiRatio = new TH2F(N3+"_hpt_rdphiRatio", "pt vs rdphi Ratio", 200, 0., 200.0, 500, -0.5, 1.5);
    hpt_rdphiRatio2 = new TH2F(N3+"_hpt_rdphiRatio2", "pt vs rdphi Ratio2", 200, 0., 200.0, 500, -0.5, 1.5);
  }

 H2DRecHit3(TString name_, TFile* file) {
    name=name_;

    hPt = (TH1F *) file->Get(name+"_hPt");
    hPa_Pt     = (TH2F *) file->Get(name+"_hPa_Pt");
    hPt_DTSeg  = (TH2F *) file->Get(name+"_hPt_DTSeg");
    hst_dt_nu  = (TH2F *) file->Get(name+"_hst_dt_nu");
    hst_dt_nu1 = (TH2F *) file->Get(name+"_hst_dt_nu1");

    hrdeta_Pt   = (TH2F *) file->Get(name+"_hrdeta_Pt");
    hrdphi_Pt   = (TH2F *) file->Get(name+"_hrdphi_Pt");

    hbeta_dphi12 = (TH2F *) file->Get(name+"_hbeta_dphi12");
    hbeta_dphi13 = (TH2F *) file->Get(name+"_hbeta_dphi13");
    hbeta_dphi14 = (TH2F *) file->Get(name+"_hbeta_dphi14");
    hbeta_dphi23 = (TH2F *) file->Get(name+"_hbeta_dphi23");
    hbeta_dphi24 = (TH2F *) file->Get(name+"_hbeta_dphi24");
    hbeta_dphi34 = (TH2F *) file->Get(name+"_hbeta_dphi34");

    hbeta_rdphi12 = (TH2F *) file->Get(name+"_hbeta_rdphi12");
    hbeta_rdphi13 = (TH2F *) file->Get(name+"_hbeta_rdphi13");
    hbeta_rdphi14 = (TH2F *) file->Get(name+"_hbeta_rdphi14");
    hbeta_rdphi23 = (TH2F *) file->Get(name+"_hbeta_rdphi23");
    hbeta_rdphi24 = (TH2F *) file->Get(name+"_hbeta_rdphi24");
    hbeta_rdphi34 = (TH2F *) file->Get(name+"_hbeta_rdphi34");

    hbeta_dphiPt12 = (TH2F *) file->Get(name+"_hbeta_dphiPt12");
    hbeta_dphiPt13 = (TH2F *) file->Get(name+"_hbeta_dphiPt13");
    hbeta_dphiPt14 = (TH2F *) file->Get(name+"_hbeta_dphiPt14");
    hbeta_dphiPt23 = (TH2F *) file->Get(name+"_hbeta_dphiPt23");
    hbeta_dphiPt24 = (TH2F *) file->Get(name+"_hbeta_dphiPt24");
    hbeta_dphiPt34 = (TH2F *) file->Get(name+"_hbeta_dphiPt34");

    hbeta_rdphiPt12 = (TH2F *) file->Get(name+"_hbeta_rdphiPt12");
    hbeta_rdphiPt13 = (TH2F *) file->Get(name+"_hbeta_rdphiPt13");
    hbeta_rdphiPt14 = (TH2F *) file->Get(name+"_hbeta_rdphiPt14");
    hbeta_rdphiPt23 = (TH2F *) file->Get(name+"_hbeta_rdphiPt23");
    hbeta_rdphiPt24 = (TH2F *) file->Get(name+"_hbeta_rdphiPt24");
    hbeta_rdphiPt34 = (TH2F *) file->Get(name+"_hbeta_rdphiPt34");

    hbeta_dphiRatio = (TH2F *) file->Get(name+"_hbeta_dphiRatio");
    hbeta_rdphiRatio = (TH2F *) file->Get(name+"_hbeta_rdphiRatio");
    hbeta_rdphiRatio2 = (TH2F *) file->Get(name+"_hbeta_rdphiRatio2");
    hpt_dphiRatio = (TH2F *) file->Get(name+"_hpt_dphiRatio");
    hpt_rdphiRatio = (TH2F *) file->Get(name+"_hpt_rdphiRatio");
    hpt_rdphiRatio2 = (TH2F *) file->Get(name+"_hpt_rdphiRatio2");
  } 

  /// Destructor
  virtual ~H2DRecHit3() {

    delete hPt;
    delete hPa_Pt;
    delete hPt_DTSeg;
    delete hst_dt_nu;
    delete hst_dt_nu1;

    delete hrdphi_Pt;
    delete hrdeta_Pt;

    delete hbeta_dphi12;
    delete hbeta_dphi13;
    delete hbeta_dphi14;
    delete hbeta_dphi23;
    delete hbeta_dphi24;
    delete hbeta_dphi34;

    delete hbeta_rdphi12;
    delete hbeta_rdphi13;
    delete hbeta_rdphi14;
    delete hbeta_rdphi23;
    delete hbeta_rdphi24;
    delete hbeta_rdphi34;

    delete hbeta_dphiPt12;
    delete hbeta_dphiPt13;
    delete hbeta_dphiPt14;
    delete hbeta_dphiPt23;
    delete hbeta_dphiPt24;
    delete hbeta_dphiPt34;

    delete hbeta_rdphiPt12;
    delete hbeta_rdphiPt13;
    delete hbeta_rdphiPt14;
    delete hbeta_rdphiPt23;
    delete hbeta_rdphiPt24;
    delete hbeta_rdphiPt34;

    delete hbeta_dphiRatio;
    delete hbeta_rdphiRatio;
    delete hbeta_rdphiRatio2;
    delete hpt_dphiRatio;
    delete hpt_rdphiRatio;
    delete hpt_rdphiRatio2;
  }

  void Fill3a(float Pt, float Pa, int dt_nu)
  {
       hPt->Fill(Pt);
       hPa_Pt->Fill(Pt,Pa);
       hPt_DTSeg->Fill(Pt,dt_nu);
  }

  void Fill4a(int k, int dt_nu, int dt_nu1)
  {
       hst_dt_nu->Fill(k,dt_nu);
       hst_dt_nu1->Fill(k,dt_nu1);
  }
  /*
  void Fill7(double ptxdphi12, double ptxdphi13, double ptxdphi14,
             double ptxdphi23, double ptxdphi24, double ptxdphi34,
             double eta1, double eta2, double eta3)
  {
     hbeta_dphiPt12->Fill(eta1,ptxdphi12);
     hbeta_dphiPt13->Fill(eta1,ptxdphi13);
     hbeta_dphiPt14->Fill(eta1,ptxdphi14);
     hbeta_dphiPt23->Fill(eta2,ptxdphi23);
     hbeta_dphiPt24->Fill(eta2,ptxdphi24);
     hbeta_dphiPt34->Fill(eta3,ptxdphi34);

  }
  void Fill7a(double ptxdphi12, double ptxdphi13, double ptxdphi14,
              double ptxdphi23, double ptxdphi24, double ptxdphi34,
              double eta1, double eta2, double eta3, double dphi13, double pt)
  {
     hbeta_rdphiPt12->Fill(eta1,ptxdphi12);
     hbeta_rdphiPt13->Fill(eta1,ptxdphi13);
     hbeta_rdphiPt14->Fill(eta1,ptxdphi14);
     hbeta_rdphiPt23->Fill(eta2,ptxdphi23);
     hbeta_rdphiPt24->Fill(eta2,ptxdphi24);
     hbeta_rdphiPt34->Fill(eta3,ptxdphi34);
     hrdphi_Pt->Fill(pt, dphi13);
     hrdeta_Pt->Fill(pt, (eta3-eta1) );
  }*/

  void Fill6_1(double dphi,double dphir,double ptxdphi,double ptxdphir,double eta,double eta_r)
  {
    hbeta_dphi12->Fill(eta,dphi);
    hbeta_rdphi12->Fill(eta_r,dphir);
    hbeta_dphiPt12->Fill(eta,ptxdphi);
    hbeta_rdphiPt12->Fill(eta_r,ptxdphir);
  }
  void Fill6_2(double dphi,double dphir,double ptxdphi,double ptxdphir,double eta,double eta_r)
  {
    hbeta_dphi13->Fill(eta,dphi);
    hbeta_rdphi13->Fill(eta_r,dphir);
    hbeta_dphiPt13->Fill(eta,ptxdphi);
    hbeta_rdphiPt13->Fill(eta_r,ptxdphir);
  }
  void Fill6_3(double dphi,double dphir,double ptxdphi,double ptxdphir,double eta,double eta_r)
  {
    hbeta_dphi14->Fill(eta,dphi);
    hbeta_rdphi14->Fill(eta_r,dphir);
    hbeta_dphiPt14->Fill(eta,ptxdphi);
    hbeta_rdphiPt14->Fill(eta_r,ptxdphir);
  }
  void Fill6_4(double dphi,double dphir,double ptxdphi,double ptxdphir,double eta,double eta_r)
  {
    hbeta_dphi23->Fill(eta,dphi);
    hbeta_rdphi23->Fill(eta_r,dphir);
    hbeta_dphiPt23->Fill(eta,ptxdphi);
    hbeta_rdphiPt23->Fill(eta_r,ptxdphir);
  }
  void Fill6_5(double dphi,double dphir,double ptxdphi,double ptxdphir,double eta,double eta_r)
  {
    hbeta_dphi24->Fill(eta,dphi);
    hbeta_rdphi24->Fill(eta_r,dphir);
    hbeta_dphiPt24->Fill(eta,ptxdphi);
    hbeta_rdphiPt24->Fill(eta_r,ptxdphir);
  }
  void Fill6_6(double dphi,double dphir,double ptxdphi,double ptxdphir,double eta,double eta_r)
  {
    hbeta_dphi34->Fill(eta,dphi);
    hbeta_rdphi34->Fill(eta_r,dphir);
    hbeta_dphiPt34->Fill(eta,ptxdphi);
    hbeta_rdphiPt34->Fill(eta_r,ptxdphir);
  }

  void Fill7b(double eta2, double dphi_R, double reta2, double rdphi_R, double pt)
  {
       hbeta_dphiRatio->Fill(eta2,dphi_R);
       hbeta_rdphiRatio->Fill(reta2,rdphi_R);
       hpt_dphiRatio->Fill(pt,dphi_R);
       hpt_rdphiRatio->Fill(pt,rdphi_R);
  }
  void Fill7c(double reta2, double rdphi_R, double pt )
  {
     hbeta_rdphiRatio2->Fill(reta2,rdphi_R);
     hpt_rdphiRatio2->Fill(pt,rdphi_R);
  }

  void Write() {

       hPt->Write();
       hPa_Pt->Write();
       hPt_DTSeg->Write();
       hst_dt_nu->Write();
       hst_dt_nu1->Write();

       hrdphi_Pt->Write();
       hrdeta_Pt->Write();

       hbeta_dphi12->Write();
       hbeta_dphi13->Write();
       hbeta_dphi14->Write();
       hbeta_dphi23->Write();
       hbeta_dphi24->Write();
       hbeta_dphi34->Write();

       hbeta_rdphi12->Write();
       hbeta_rdphi13->Write();
       hbeta_rdphi14->Write();
       hbeta_rdphi23->Write();
       hbeta_rdphi24->Write();
       hbeta_rdphi34->Write();

       hbeta_dphiPt12->Write();
       hbeta_dphiPt13->Write();
       hbeta_dphiPt14->Write();
       hbeta_dphiPt23->Write();
       hbeta_dphiPt24->Write();
       hbeta_dphiPt34->Write();

       hbeta_rdphiPt12->Write();
       hbeta_rdphiPt13->Write();
       hbeta_rdphiPt14->Write();
       hbeta_rdphiPt23->Write();
       hbeta_rdphiPt24->Write();
       hbeta_rdphiPt34->Write();

       hbeta_dphiRatio->Write();
       hbeta_rdphiRatio->Write();
       hbeta_rdphiRatio2->Write();
       hpt_dphiRatio->Write();
       hpt_rdphiRatio->Write();
       hpt_rdphiRatio2->Write();
  }

  TH1F *hPt;
  TH2F *hPa_Pt;
  TH2F *hPt_DTSeg;
  TH2F *hst_dt_nu;
  TH2F *hst_dt_nu1;

  TH2F *hrdphi_Pt;
  TH2F *hrdeta_Pt;

  TH2F *hbeta_dphi12;
  TH2F *hbeta_dphi13;
  TH2F *hbeta_dphi14;
  TH2F *hbeta_dphi23;
  TH2F *hbeta_dphi24;
  TH2F *hbeta_dphi34;

  TH2F *hbeta_rdphi12;
  TH2F *hbeta_rdphi13;
  TH2F *hbeta_rdphi14;
  TH2F *hbeta_rdphi23;
  TH2F *hbeta_rdphi24;
  TH2F *hbeta_rdphi34;

  TH2F *hbeta_dphiPt12;
  TH2F *hbeta_dphiPt13;
  TH2F *hbeta_dphiPt14;
  TH2F *hbeta_dphiPt23;
  TH2F *hbeta_dphiPt24;
  TH2F *hbeta_dphiPt34;

  TH2F *hbeta_rdphiPt12;
  TH2F *hbeta_rdphiPt13;
  TH2F *hbeta_rdphiPt14;
  TH2F *hbeta_rdphiPt23;
  TH2F *hbeta_rdphiPt24;
  TH2F *hbeta_rdphiPt34;

  TH2F *hbeta_dphiRatio;
  TH2F *hbeta_rdphiRatio;
  TH2F *hbeta_rdphiRatio2;
  TH2F *hpt_dphiRatio;
  TH2F *hpt_rdphiRatio;
  TH2F *hpt_rdphiRatio2;

  TString name;
};

class H2DRecHit4 {
public:

 H2DRecHit4(std::string name_) {
    TString N4 = name_.c_str();
    name=N4;

    heta_dphiA  = new TH2F(N4+"_heta_dphiA", "eta vs dphit @ME", 150, 1., 2.5, 200, -0.03, 0.03);
    heta_rdphiA = new TH2F(N4+"_heta_rdphiA", "eta vs rdphi @ME", 150, 1., 2.5, 200, -0.03, 0.03);
    heta_detaA  = new TH2F(N4+"_heta_detaA", "eta vs deta  @ME", 150, 1., 2.5, 101, -0.01, 0.01);
    heta_rdetaA  = new TH2F(N4+"_heta_rdetaA", "eta vs rdeta  @ME", 150, 1., 2.5, 101, -0.01, 0.01);

    heta_dphiPtA = new TH2F(N4+"_heta_dphiPtA", "eta vs dphi*Pt @ME", 150, 1., 2.5, 200, -1.0, 1.0);
    heta_rdphiPtA = new TH2F(N4+"_heta_rdphiPtA", "eta vs rdphi*Pt @ME", 150, 1., 2.5, 200, -1.0, 1.0);

    heta_Pt    = new TH2F(N4+"_heta_Pt", "eta vs Pt @ME", 150, 1., 2.5, 50, 5.0, 205.0);
    heta_dphir = new TH2F(N4+"_heta_dphir", "eta vs dphi resolution @ME", 150, 1., 2.5, 101, -0.0005, 0.0005);
    hpt_dphiA  = new TH2F(N4+"_hpt_dphiA", "pt vs dphi @ME", 50, 5., 205, 200, 0.0, 0.1);
    hpt_rdphiA = new TH2F(N4+"_hpt_rdphiA", "pt vs rdphi @ME", 50, 5., 205, 200, 0.0, 0.1);
  }

 H2DRecHit4(TString name_, TFile* file) {
    name=name_;

    heta_dphiA = (TH2F *) file->Get(name+"_heta_dphiA");
    heta_rdphiA = (TH2F *) file->Get(name+"_heta_rdphiA");
    heta_detaA = (TH2F *) file->Get(name+"_heta_detaA");
    heta_rdetaA = (TH2F *) file->Get(name+"_heta_rdetaA");

    heta_dphiPtA = (TH2F *) file->Get(name+"_heta_dphiPtA");
    heta_rdphiPtA = (TH2F *) file->Get(name+"_heta_rdphiPtA");

    heta_Pt    = (TH2F *) file->Get(name+"_heta_Pt");
    heta_dphir = (TH2F *) file->Get(name+"_heta_dphir");
    hpt_dphiA  = (TH2F *) file->Get(name+"_hpt_dphiA");
    hpt_rdphiA = (TH2F *) file->Get(name+"_hpt_rdphiA");
  } 

  /// Destructor
  virtual ~H2DRecHit4() {

    delete heta_dphiA;
    delete heta_rdphiA;
    delete heta_detaA;
    delete heta_rdetaA;

    delete heta_dphiPtA;
    delete heta_rdphiPtA;

    delete heta_dphir;
    delete heta_Pt;
    delete hpt_dphiA;
    delete hpt_rdphiA;
  }

  void Fill8(double ptxdphi, double dphi, double deta, double eta, double pt)
  {
     heta_dphiA->Fill(eta,dphi);
     heta_detaA->Fill(eta,deta);
     heta_dphiPtA->Fill(eta,ptxdphi);
     hpt_dphiA->Fill(pt,dphi);
  }
  void Fill8a(double ptxrdphi, double rdphi, double rdeta, double reta, double pt, double dphir)
  {
     heta_rdphiA->Fill(reta,rdphi);
     heta_rdetaA->Fill(reta,rdeta);
     heta_rdphiPtA->Fill(reta,ptxrdphi);
     hpt_rdphiA->Fill(pt,rdphi);
     heta_Pt->Fill(reta,pt);
     heta_dphir->Fill(reta,dphir);
  }

  void Write() {
       heta_dphiA->Write();
       heta_rdphiA->Write();
       heta_detaA->Write();
       heta_rdetaA->Write();

       heta_dphiPtA->Write();
       heta_rdphiPtA->Write();

       hpt_dphiA->Write();
       hpt_rdphiA->Write();
       heta_Pt->Write();
       heta_dphir->Write();
  }

  TH2F *heta_dphiA;
  TH2F *heta_rdphiA;
  TH2F *heta_detaA;
  TH2F *heta_rdetaA;

  TH2F *heta_dphiPtA;
  TH2F *heta_rdphiPtA;

  TH2F *hpt_dphiA;
  TH2F *hpt_rdphiA;
  TH2F *heta_Pt;
  TH2F *heta_dphir;

  TString name;
};

class H2DRecHit5 {
public:

 H2DRecHit5(std::string name_) {
    TString N5 = name_.c_str();
    name=N5;

    heta_dphiA  = new TH2F(N5+"_heta_dphiA", "eta vs dphi @MB", 110, 0., 1.1, 200, -0.03, 0.03);
    heta_rdphiA = new TH2F(N5+"_heta_rdphiA", "eta vs rdphi @MB", 110, 0., 1.1, 200, -0.03, 0.03);
    heta_detaA  = new TH2F(N5+"_heta_detaA", "eta vs deta  @MB", 110, 0., 1.1, 101, -0.01, 0.01);
    heta_rdetaA  = new TH2F(N5+"_heta_rdetaA", "eta vs rdeta  @MB", 110, 0., 1.1, 101, -0.01, 0.01);

    heta_dphiPtA  = new TH2F(N5+"_heta_dphiPtA", "eta vs dphi*Pt @MB", 110, 0., 1.1, 200, -1.0, 1.0);
    heta_rdphiPtA = new TH2F(N5+"_heta_rdphiPtA", "eta vs rdphi*Pt @MB", 110, 0., 1.1, 200, -1.0, 1.0);

    heta_Pt    = new TH2F(N5+"_heta_Pt", "eta vs Pt @MB", 110, 0., 1.1, 50, 5.0, 205.0);
    heta_dphir = new TH2F(N5+"_heta_dphir", "eta vs dphi resolution @MB", 110, 0., 1.1, 101, -0.0005, 0.0005);
    hpt_dphiA  = new TH2F(N5+"_hpt_dphiA", "pt vs dphi @MB", 50, 5., 205.0, 200, 0.0, 0.1);
    hpt_rdphiA = new TH2F(N5+"_hpt_rdphiA", "pt vs rdphi @MB", 50, 5., 205.0, 200, 0.0, 0.1);
  }

 H2DRecHit5(TString name_, TFile* file) {
    name=name_;

    heta_dphiA = (TH2F *) file->Get(name+"_heta_dphiA");
    heta_rdphiA = (TH2F *) file->Get(name+"_heta_rdphiA");
    heta_detaA = (TH2F *) file->Get(name+"_heta_detaA");
    heta_rdetaA = (TH2F *) file->Get(name+"_heta_rdetaA");

    heta_dphiPtA = (TH2F *) file->Get(name+"_heta_dphiPtA");
    heta_rdphiPtA = (TH2F *) file->Get(name+"_heta_rdphiPtA");

    heta_Pt    = (TH2F *) file->Get(name+"_heta_Pt");
    heta_dphir = (TH2F *) file->Get(name+"_heta_dphir");
    hpt_dphiA  = (TH2F *) file->Get(name+"_hpt_dphiA");
    hpt_rdphiA = (TH2F *) file->Get(name+"_hpt_rdphiA");
  } 

  /// Destructor
  virtual ~H2DRecHit5() {

    delete heta_dphiA;
    delete heta_rdphiA;
    delete heta_detaA;
    delete heta_rdetaA;

    delete heta_dphiPtA;
    delete heta_rdphiPtA;

    delete heta_Pt;
    delete heta_dphir;
    delete hpt_dphiA;
    delete hpt_rdphiA;
  }

  void Fill9(double ptxdphi,double dphi, double deta, double eta, double pt)
  {
     heta_dphiA->Fill(eta,dphi);
     heta_detaA->Fill(eta,deta);
     heta_dphiPtA->Fill(eta,ptxdphi);
     hpt_dphiA->Fill(pt,dphi);
  }
  void Fill9a(double ptxrdphi, double rdphi, double rdeta, double reta, double pt, double dphir)
  {
     heta_rdphiA->Fill(reta,rdphi);
     heta_rdetaA->Fill(reta,rdeta);
     heta_rdphiPtA->Fill(reta,ptxrdphi);
     hpt_rdphiA->Fill(pt,rdphi);
     heta_Pt->Fill(reta,pt);
     heta_dphir->Fill(reta,dphir);
  }

  void Write() {

       heta_dphiA->Write();
       heta_rdphiA->Write();
       heta_detaA->Write();
       heta_rdetaA->Write();

       heta_dphiPtA->Write();
       heta_rdphiPtA->Write();

       hpt_dphiA->Write();
       hpt_rdphiA->Write();
       heta_Pt->Write();
       heta_dphir->Write();
  }

  TH2F *heta_dphiA;
  TH2F *heta_rdphiA;
  TH2F *heta_detaA;
  TH2F *heta_rdetaA;

  TH2F *heta_dphiPtA;
  TH2F *heta_rdphiPtA;

  TH2F *hpt_dphiA;
  TH2F *hpt_rdphiA;
  TH2F *heta_Pt;
  TH2F *heta_dphir;

  TString name;
};

class H2DRecHit6 {
public:

 H2DRecHit6(std::string name_) {
    TString N6 = name_.c_str();
    name=N6;

    heta_dphi1  = new TH2F(N6+"_heta_dphi1", "eta vs dphit @ME", 150, 1., 2.5, 200, -0.03, 0.03);
    heta_rdphi1 = new TH2F(N6+"_heta_rdphi1", "eta vs rdphi @ME", 150, 1., 2.5, 200, -0.03, 0.03);

    heta_phir  = new TH2F(N6+"_heta_phir", "eta vs phi resolution  @ME", 150, 1., 2.5, 101, -0.01, 0.01);
    heta_etar  = new TH2F(N6+"_heta_etar", "eta vs eta resolution  @ME", 150, 1., 2.5, 101, -0.01, 0.01);

    heta_dphiPt1 = new TH2F(N6+"_heta_dphiPt1", "eta vs dphi*Pt 1", 150, 1., 2.5, 400, -3., 3.);
    heta_rdphiPt1 = new TH2F(N6+"_heta_rdphiPt1", "eta vs rdphi*Pt 1", 150, 1., 2.5, 400, -3., 3.);

    heta_Pt1   = new TH2F(N6+"_heta_Pt1", "eta vs Pt @ME", 150, 1., 2.5, 50, 5.0, 205.0);
    hpt_dphiA1 = new TH2F(N6+"_hpt_dphiA1", "pt vs dphi @ME", 50, 5., 205, 200, 0.0, 0.1);
    hpt_rdphiA1 = new TH2F(N6+"_hpt_rdphiA1", "pt vs rdphi @ME", 50, 5., 205, 200, 0.0, 0.1);
  }

 H2DRecHit6(TString name_, TFile* file) {
    name=name_;

    heta_dphi1  = (TH2F *) file->Get(name+"_heta_dphi1");
    heta_rdphi1 = (TH2F *) file->Get(name+"_heta_rdphi1");

    heta_phir = (TH2F *) file->Get(name+"_heta_phir");
    heta_etar = (TH2F *) file->Get(name+"_heta_etar");

    heta_Pt1 = (TH2F *) file->Get(name+"_heta_Pt1");
    heta_dphiPt1 = (TH2F *) file->Get(name+"_heta_dphiPt1");
    heta_rdphiPt1 = (TH2F *) file->Get(name+"_heta_rdphiPt1");

    hpt_dphiA1 = (TH2F *) file->Get(name+"_hpt_dphiA");
    hpt_rdphiA1 = (TH2F *) file->Get(name+"_hpt_rdphiA");
  } 

  /// Destructor
  virtual ~H2DRecHit6() {

    delete  heta_dphi1;
    delete  heta_rdphi1;

    delete  heta_phir;
    delete  heta_etar;

    delete heta_dphiPt1;
    delete heta_rdphiPt1;

    delete heta_Pt1;
    delete hpt_dphiA1;
    delete hpt_rdphiA1;
  }

  void Fill8b(double ptxdphi1, double dphi, double eta, double pt) 
  {
     heta_dphiPt1->Fill(eta,ptxdphi1);
     heta_dphi1->Fill(eta,dphi);
     hpt_dphiA1->Fill(pt,dphi);
  }
  void Fill8c(double ptxrdphi1, double rdphi, double phir, double etar, double reta, double pt) 
  {
     heta_rdphiPt1->Fill(reta,ptxrdphi1);
     heta_rdphi1->Fill(reta,rdphi);

     heta_phir->Fill(reta,phir);
     heta_etar->Fill(reta,etar);
     hpt_rdphiA1->Fill(pt,rdphi);
     heta_Pt1->Fill(reta,pt);
  }

  void Write() {
 
       heta_dphi1->Write();
       heta_rdphi1->Write();

       heta_phir->Write();
       heta_etar->Write();

       heta_dphiPt1->Write();
       heta_rdphiPt1->Write();

       hpt_dphiA1->Write();
       hpt_rdphiA1->Write();
       heta_Pt1->Write();
  }
 
  TH2F *heta_dphi1;
  TH2F *heta_rdphi1;
 
  TH2F *heta_phir;
  TH2F *heta_etar;

  TH2F *heta_dphiPt1;
  TH2F *heta_rdphiPt1;

  TH2F *hpt_dphiA1;
  TH2F *hpt_rdphiA1;
  TH2F *heta_Pt1;

  TString name;
};

class H2DRecHit7 {
public:

 H2DRecHit7(std::string name_) {
    TString N7 = name_.c_str();
    name=N7;

    heta_dphi1  = new TH2F(N7+"_heta_dphi1", "eta vs dphit @MB", 110, 0., 1.1, 200, -0.03, 0.03);
    heta_rdphi1 = new TH2F(N7+"_heta_rdphi1", "eta vs rdphi @MB", 110, 0., 1.1, 200, -0.03, 0.03);

    heta_phir  = new TH2F(N7+"_heta_phir", "eta vs phi resolution  @MB", 110, 0., 1.1, 101, -0.01, 0.01);
    heta_etar  = new TH2F(N7+"_heta_etar", "eta vs eta resolution  @MB", 110, 0., 1.1, 101, -0.01, 0.01);

    heta_dphiPt1  = new TH2F(N7+"_heta_dphiPt1", "eta vs dphi*Pt @MB", 110, 0., 1.1, 400, -3., 3.);
    heta_rdphiPt1 = new TH2F(N7+"_heta_rdphiPt1", "eta vs rdphi*Pt @MB", 110, 0., 1.1, 400, -3., 3.);

    heta_Pt1   = new TH2F(N7+"_heta_Pt1", "eta vs Pt @MB", 150, 1., 2.5, 50, 5.0, 205.0);
    hpt_dphiA1 = new TH2F(N7+"_hpt_dphiA1", "pt vs dphi @ME", 50, 5., 205, 200, 0.0, 0.1);
    hpt_rdphiA1 = new TH2F(N7+"_hpt_rdphiA1", "pt vs rdphi @ME", 50, 5., 205, 200, 0.0, 0.1);
  }

 H2DRecHit7(TString name_, TFile* file) {
    name=name_;

    heta_dphi1  = (TH2F *) file->Get(name+"_heta_dphi1");
    heta_rdphi1 = (TH2F *) file->Get(name+"_heta_rdphi1");

    heta_phir = (TH2F *) file->Get(name+"_heta_phir");
    heta_etar = (TH2F *) file->Get(name+"_heta_etar");

    heta_dphiPt1 = (TH2F *) file->Get(name+"_heta_dphiPt1");
    heta_rdphiPt1 = (TH2F *) file->Get(name+"_heta_rdphiPt1");

    heta_Pt1 = (TH2F *) file->Get(name+"_heta_Pt1");
    hpt_dphiA1 = (TH2F *) file->Get(name+"_hpt_dphiA");
    hpt_rdphiA1 = (TH2F *) file->Get(name+"_hpt_rdphiA");
  } 

  /// Destructor
  virtual ~H2DRecHit7() {

    delete  heta_dphi1;
    delete  heta_rdphi1;

    delete  heta_phir;
    delete  heta_etar;

    delete heta_dphiPt1;
    delete heta_rdphiPt1;

    delete heta_Pt1;
    delete hpt_dphiA1;
    delete hpt_rdphiA1;
  }

  void Fill9b(double ptxdphi1, double dphi, double eta, double pt)
  {
     heta_dphiPt1->Fill(eta,ptxdphi1);
     heta_dphi1->Fill(eta,dphi);
     hpt_dphiA1->Fill(pt,dphi);
  }
  void Fill9c(double ptxrdphi1, double rdphi, double phir, double etar, double reta, double pt) 
  {
     heta_rdphiPt1->Fill(reta,ptxrdphi1);
     heta_rdphi1->Fill(reta,rdphi);

     heta_phir->Fill(reta,phir);
     heta_etar->Fill(reta,etar);
     hpt_rdphiA1->Fill(pt,rdphi);
     heta_Pt1->Fill(reta,pt);
  }

  void Write() {

       heta_dphi1->Write();
       heta_rdphi1->Write();

       heta_phir->Write();
       heta_etar->Write();

       heta_dphiPt1->Write();
       heta_rdphiPt1->Write();

       hpt_dphiA1->Write();
       hpt_rdphiA1->Write();
       heta_Pt1->Write();
  }

  TH2F *heta_dphi1;
  TH2F *heta_rdphi1;
 
  TH2F *heta_phir;
  TH2F *heta_etar;

  TH2F *heta_dphiPt1;
  TH2F *heta_rdphiPt1;

  TH2F *hpt_dphiA1;
  TH2F *hpt_rdphiA1;
  TH2F *heta_Pt1;

  TString name;
};

class H2DRecHit8 {
public:

 H2DRecHit8(std::string name_) {
    TString N8 = name_.c_str();
    name=N8;

    heta_dphiRatio = new TH2F(N8+"_heta_dphiRatio", "eta vs dphi Ratio", 150, 1., 2.5, 500, -0.5, 1.5);
    heta_rdphiRatio = new TH2F(N8+"_heta_rdphiRatio", "eta vs rdphi Ratio", 150, 1., 2.5, 500, -0.5, 1.5);
    hpt_dphiRatio = new TH2F(N8+"_hpt_dphiRatio", "pt vs dphi Ratio", 2000, 0., 200.0, 500, -0.5, 1.5);
    hpt_rdphiRatio = new TH2F(N8+"_hpt_rdphiRatio", "pt vs rdphi Ratio", 200, 0., 200.0, 500, -0.5, 1.5);

  }

 H2DRecHit8(TString name_, TFile* file) {
    name=name_;

    heta_dphiRatio = (TH2F *) file->Get(name+"_heta_dphiRatio");
    heta_rdphiRatio = (TH2F *) file->Get(name+"_heta_rdphiRatio");
    hpt_dphiRatio = (TH2F *) file->Get(name+"_hpt_dphiRatio");
    hpt_rdphiRatio = (TH2F *) file->Get(name+"_hpt_rdphiRatio");

  } 

  /// Destructor
  virtual ~H2DRecHit8() {

    delete heta_dphiRatio;
    delete heta_rdphiRatio;
    delete hpt_dphiRatio;
    delete hpt_rdphiRatio;
  }

  void Fill10(double dphi_rt, double rdphi_rt, double eta, double reta, double pt)
  {
     heta_dphiRatio->Fill(eta,dphi_rt);
     heta_rdphiRatio->Fill(reta,rdphi_rt);
     hpt_dphiRatio->Fill(pt,dphi_rt);
     hpt_rdphiRatio->Fill(pt,rdphi_rt);
  }

  void Write() {

       heta_dphiRatio->Write();
       heta_rdphiRatio->Write();
       hpt_dphiRatio->Write();
       hpt_rdphiRatio->Write();
  }

  TH2F *heta_dphiRatio;
  TH2F *heta_rdphiRatio;
  TH2F *hpt_dphiRatio;
  TH2F *hpt_rdphiRatio;

  TString name;
};

class H2DRecHit9 {
public:

 H2DRecHit9(std::string name_) {
    TString N9 = name_.c_str();
    name=N9;

    heta_dphiRatio = new TH2F(N9+"_heta_dphiRatio", "eta vs dphi Ratio", 110, 0., 1.1, 500, -0.5, 1.5);
    heta_rdphiRatio = new TH2F(N9+"_heta_rdphiRatio", "eta vs rdphi Ratio", 110, 0., 1.1, 500, -0.5, 1.5);
    hpt_dphiRatio = new TH2F(N9+"_hpt_dphiRatio", "pt vs dphi Ratio", 2000, 0., 200.0, 500, -0.5, 1.5);
    hpt_rdphiRatio = new TH2F(N9+"_hpt_rdphiRatio", "pt vs rdphi Ratio", 200, 0., 200.0, 500, -0.5, 1.5);

  }

 H2DRecHit9(TString name_, TFile* file) {
    name=name_;

    heta_dphiRatio = (TH2F *) file->Get(name+"_heta_dphiRatio");
    heta_rdphiRatio = (TH2F *) file->Get(name+"_heta_rdphiRatio");
    hpt_dphiRatio = (TH2F *) file->Get(name+"_hpt_dphiRatio");
    hpt_rdphiRatio = (TH2F *) file->Get(name+"_hpt_rdphiRatio");

  } 

  /// Destructor
  virtual ~H2DRecHit9() {

    delete heta_dphiRatio;
    delete heta_rdphiRatio;
    delete hpt_dphiRatio;
    delete hpt_rdphiRatio;
  }

  void Fill11(double dphi_rt, double rdphi_rt, double eta, double reta, double pt)
  {
     heta_dphiRatio->Fill(eta,dphi_rt);
     heta_rdphiRatio->Fill(reta,rdphi_rt);
     hpt_dphiRatio->Fill(pt,dphi_rt);
     hpt_rdphiRatio->Fill(pt,rdphi_rt);
  }

  void Write() {

       heta_dphiRatio->Write();
       heta_rdphiRatio->Write();
       hpt_dphiRatio->Write();
       hpt_rdphiRatio->Write();
  }

  TH2F *heta_dphiRatio;
  TH2F *heta_rdphiRatio;
  TH2F *hpt_dphiRatio;
  TH2F *hpt_rdphiRatio;

  TString name;
};

class H2DRecHit10 {
public:

 H2DRecHit10(std::string name_) {
    TString N10 = name_.c_str();
    name=N10;

    heta_dphiA  = new TH2F(N10+"_heta_dphiA", "eta vs dphi @MB", 100, 0.5, 1.5, 200, -0.08, 0.08);
    heta_rdphiA = new TH2F(N10+"_heta_rdphiA", "eta vs rdphi @MB", 100, 0.5, 1.5, 200, -0.08, 0.08);
    heta_detaA  = new TH2F(N10+"_heta_detaA", "eta vs deta  @MB", 100, 0.5, 1.5, 101, -0.01, 0.01);
    heta_rdetaA  = new TH2F(N10+"_heta_rdetaA", "eta vs rdeta  @MB", 100, 0.5, 1.5, 101, -0.01, 0.01);

    heta_dphiPtA  = new TH2F(N10+"_heta_dphiPtA", "eta vs dphi*Pt @MB", 100, 0.5, 1.5, 200, -1.0, 1.0);
    heta_rdphiPtA = new TH2F(N10+"_heta_rdphiPtA", "eta vs rdphi*Pt @MB", 100, 0.5, 1.5, 200, -1.0, 1.0);

    heta_Pt    = new TH2F(N10+"_heta_Pt", "eta vs Pt @MB", 100, 0.5, 1.5, 50, 5.0, 205.0);
    heta_dphir = new TH2F(N10+"_heta_dphir", "eta vs dphi resolution @MB", 100, 0.5, 1.5, 101, -0.0005, 0.0005);
    hpt_dphiA  = new TH2F(N10+"_hpt_dphiA", "pt vs dphi @MB", 50, 5., 205.0, 200, 0.0, 0.1);
    hpt_rdphiA = new TH2F(N10+"_hpt_rdphiA", "pt vs rdphi @MB", 50, 5., 205.0, 200, 0.0, 0.1);
  }

 H2DRecHit10(TString name_, TFile* file) {
    name=name_;

    heta_dphiA = (TH2F *) file->Get(name+"_heta_dphiA");
    heta_rdphiA = (TH2F *) file->Get(name+"_heta_rdphiA");
    heta_detaA = (TH2F *) file->Get(name+"_heta_detaA");
    heta_rdetaA = (TH2F *) file->Get(name+"_heta_rdetaA");

    heta_dphiPtA = (TH2F *) file->Get(name+"_heta_dphiPtA");
    heta_rdphiPtA = (TH2F *) file->Get(name+"_heta_rdphiPtA");

    heta_Pt    = (TH2F *) file->Get(name+"_heta_Pt");
    heta_dphir = (TH2F *) file->Get(name+"_heta_dphir");
    hpt_dphiA  = (TH2F *) file->Get(name+"_hpt_dphiA");
    hpt_rdphiA = (TH2F *) file->Get(name+"_hpt_rdphiA");
  } 

  /// Destructor
  virtual ~H2DRecHit10() {

    delete heta_dphiA;
    delete heta_rdphiA;
    delete heta_detaA;
    delete heta_rdetaA;

    delete heta_dphiPtA;
    delete heta_rdphiPtA;

    delete heta_Pt;
    delete heta_dphir;
    delete hpt_dphiA;
    delete hpt_rdphiA;
  }

  void Fill12(double ptxdphi,double dphi, double deta, double eta, double pt)
  {
     heta_dphiA->Fill(eta,dphi);
     heta_detaA->Fill(eta,deta);
     heta_dphiPtA->Fill(eta,ptxdphi);
     hpt_dphiA->Fill(pt,dphi);
  }
  void Fill12a(double ptxrdphi, double rdphi, double rdeta, double reta, double pt, double dphir)
  {
     heta_rdphiA->Fill(reta,rdphi);
     heta_rdetaA->Fill(reta,rdeta);
     heta_rdphiPtA->Fill(reta,ptxrdphi);
     hpt_rdphiA->Fill(pt,rdphi);
     heta_Pt->Fill(reta,pt);
     heta_dphir->Fill(reta,dphir);
  }

  void Write() {

       heta_dphiA->Write();
       heta_rdphiA->Write();
       heta_detaA->Write();
       heta_rdetaA->Write();

       heta_dphiPtA->Write();
       heta_rdphiPtA->Write();

       hpt_dphiA->Write();
       hpt_rdphiA->Write();
       heta_Pt->Write();
       heta_dphir->Write();
  }

  TH2F *heta_dphiA;
  TH2F *heta_rdphiA;
  TH2F *heta_detaA;
  TH2F *heta_rdetaA;

  TH2F *heta_dphiPtA;
  TH2F *heta_rdphiPtA;

  TH2F *hpt_dphiA;
  TH2F *hpt_rdphiA;
  TH2F *heta_Pt;
  TH2F *heta_dphir;

  TString name;
};

#endif
