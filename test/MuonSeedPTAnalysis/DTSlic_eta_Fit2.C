void DTSlic_eta_Fit2(int s1, int s2, double f1, double f2){


  /* Get the profile histogram for Pt*dphi vs eta for non-station 1 case
   * Slice each pt bin and do the gaussian fit
   * plot the mean and the sigma against the middle of the Y
   *
   * Author:  Shih-Chuan Kao  --- UCR
   */


// Number of sigmas around mean to fit gaussian.  It uses 2 iterations 
// i.e. range is = [mu - nsigmas * sigma, mu + nsigmas * sigma]

float  nsigmas = 1.25;  

TFile *file = TFile::Open("seed_5-200DTCSC_1.root");
FILE *dfile = fopen("SMB_Data.txt","a");

TString suffixps = ".jpg";

char dphi_case[19];
sprintf(dphi_case,"SMB_%d%d_heta_dphiPt1",s1,s2);
//char dphi_case[20];
//sprintf(dphi_case,"SMB_%d%d_heta_rdphiPt1",s1,s2);
TString dphi_type = dphi_case;

char dphi_case1[17];
sprintf(dphi_case1,"SMB_%d%d_heta_dphi1",s1,s2);
//char dphi_case1[18];
//sprintf(dphi_case1,"SMB_%d%d_heta_rdphi1",s1,s2);
TString dphi_type1 = dphi_case1;
char deta_case1[16];
sprintf(deta_case1,"SMB_%d%d_heta_etar",s1,s2);
TString deta_type1 = deta_case1;


 TString plot01 = dphi_type+"_1"+suffixps;  
 TString plot02 = dphi_type1+"_2"+suffixps;
 TString plot05 = dphi_type+"_5"+suffixps;
 TString plot06 = dphi_type+"_6"+suffixps;


// ********************************************************************
// Pointers to histograms
// ********************************************************************

    heta_dphiPt  = (TH2F *) file->Get(dphi_type);
    heta_dphi    = (TH2F *) file->Get(dphi_type1);
    heta_deta    = (TH2F *) file->Get(deta_type1);

// *****************************************************************
// main program -- for pt : 5 ~ 200 GeV
// *****************************************************************

 // define the fitting function
 Double_t fitf(Double_t *x, Double_t *par)
{
   
          Double_t fitval = (par[0]*x[0])+par[1]
; 
          return fitval;

 }



 // *****************************************
 // ***** 1 hdeta vs. Pt*dphi   Low eta *****
 // *****************************************

 double r=0.01;

 double prf1[55]={0.0};
 double xbin1[55]={0.0};
 double prfErrY1[55]={0.0};
 double prfErrX1[55]={0.0};
 for (int i=0; i<55; i++) {
     prfErrX1[i]=0.01;
 }

 int j=0;
 for (int i=1; i<56; i++) {
     j=(2*i)-1;

     heta_dphiPt->ProjectionY("heta_dphiPt_pjy",j,j+1);

     //double mean1 = heta_dphiPt_pjy->GetMean();
     //double rms1 = heta_dphiPt_pjy->GetRMS();
     //double L1 = mean1 - rms1;
     //double H1 = mean1 + rms1;

     heta_dphiPt_pjy->Fit("gaus","N0");
     double par0 = gaus->GetParameter(0);
     double par1 = gaus->GetParameter(1);
     double par2 = gaus->GetParameter(2);
     cout << "================="<<i+1<<" First  fit ================" << endl;
     cout << "Parameters are: " << "P0: " <<par0<< " P1: " <<par1<< " P2: " <<par2<< endl;
     double L2 = par1 -nsigmas * par2;
     double H2 = par1 + nsigmas * par2;

     //gaus->SetParLimits(1,0.,3.);

     heta_dphiPt_pjy->Fit("gaus","N0R","",L2,H2);
     par0 = gaus->GetParameter(0);
     par1 = gaus->GetParameter(1);
     par2 = gaus->GetParameter(2);
     cout << "************ Second fit ***********" << endl;
     cout << "Parameters are: " << "P0: " <<par0<< " P1: " <<par1<< " P2: " <<par2<< endl;

     double nu = heta_dphiPt_pjy->GetEntries();
     if (nu < 20.) {
        par2 =0.0;
        par1 =0.0;
     }

     prf1[i-1]=par1;
     prfErrY1[i-1]=par2;
     xbin1[i-1]=r;

     r=r+0.02;

 }

 gStyle->SetOptStat(kTRUE);
 gStyle->SetOptFit(0111);  
 c1 = new TCanvas("c1");
 c1->SetFillColor(10);
 c1->SetGrid(); 
 c1->GetFrame()->SetFillColor(21);
 c1->GetFrame()->SetBorderSize(12);
 c1->cd();

 eta_dphiPt_prf1 = new TGraphErrors(55,xbin1,prf1,prfErrX1,prfErrY1);
 eta_dphiPt_prf1->SetTitle(dphi_type);
 eta_dphiPt_prf1->SetMarkerColor(4);
 eta_dphiPt_prf1->SetMarkerStyle(21);

 TF1 *func = new TF1("fitf",fitf,f1,f2,2);
 eta_dphiPt_prf1->Fit("fitf","R","",f1,f2);

 double para0 = func->GetParameter(0);
 double para1 = func->GetParameter(1);
 double perr0 = func->GetParError(0);
 double perr1 = func->GetParError(1);

 //eta_dphiPt_prf1->GetXaxis()->SetTitle("Pt*dphi in different eta ranges");
 //eta_dphiPt_prf1->GetYaxis()->SetTitle(" ");
 eta_dphiPt_prf1->Draw("AP");
 eta_dphiPt_prf1->Print();

 c1->Update();
 c1->Print(plot01);

 // *****************************
 // ****** 2 dphi Fitting *******
 // *****************************
 
 double prf1a[55]={0.0};
 double xbin1a[55]={0.0};
 double prfErrY1a[55]={0.0};
 double prfErrX1a[55]={0.0};


 for (int i=0; i<55; i++) {
     prfErrX1a[i]=0.01;
 }


 double pt_rt1[55]={0.0};
 r=0.01;
 j=0;
 for (int i=1; i<56; i++) {
     j=(2*i)-1;

     heta_dphi->ProjectionY("heta_dphi_pjy",j,j+1);
     heta_deta->ProjectionY("heta_deta_pjy",j,j+1);

     heta_dphi_pjy->Fit("landau","N0");
     double pr1 = landau->GetParameter(1);
     double pr2 = landau->GetParameter(2);
     //double L2 = pr1 - nsigmas * pr2;
     //double H2 = pr1 + nsigmas * pr2;
     //heta_dphi_pjy->Fit("gaus","N0R","",L2,H2);
     //pr1 = gaus->GetParameter(1);
     //pr2 = gaus->GetParameter(2);

     heta_deta_pjy->Fit("gaus","N0");
     double pa1 = landau->GetParameter(1);
     double pa2 = landau->GetParameter(2);
     double L3 = pa1 - nsigmas * pa2;
     double H3 = pa1 + nsigmas * pa2;
     heta_deta_pjy->Fit("gaus","N0R","",L3,H3);
     pa1 = gaus->GetParameter(1);
     pa2 = gaus->GetParameter(2);

   
     //gaus->SetParLimits(1,0.,1.);
       

     double nu1 = heta_dphi_pjy->GetEntries();
     //if ((pr1 <= 0.)||(nu1 < 20.)) {
     //    pr1 =0.00001;
     if (nu1 < 20.) { 
           pr1 = 0.0;
           pr2 =0.0;

     }

     prf1a[i-1]=pr1;
     prfErrY1a[i-1]=pr2;
     xbin1a[i-1]=r;


     // **************************************
     // ***** Pt uncertainty calculation *****
     // **************************************

     if ((r > f2)||(r < f1)) {
        pt_rt1[i-1]=0.0;
     }
     else {
         //double sigma_etaf = (perr1*perr1) + (perr0*perr0*xbin1a[i-1]*xbin1a[i-1]) + (para0*para0*pa2*pa2);
         double sigma_etaf = (para0*para0*pa2*pa2);
         double etaf = (para1+ (para0*xbin1a[i-1]))*(para1+ (para0*xbin1a[i-1]));
        if ((pr1!=0.0)&&(etaf!=0.0)){
           double eta_rt = sigma_etaf / etaf ;
           double phi_rt = (pr2/pr1)*(pr2/pr1);
           pt_rt1[i-1] = sqrt(eta_rt + phi_rt);
        }
        if (pt_rt1[i-1] > 1.0) {
           pt_rt1[i-1] = 0.0;
        }
     }
     if (pt_rt1[i-1] != 0.0) { 
        fprintf (dfile, "%d%d %f %f %f %f\n",s1,s2,r,para0,para1,pt_rt1[i-1]);
     }
     r=r+0.02;

 }



 gStyle->SetOptStat(kTRUE);
 gStyle->SetOptFit(0111);  
 c2 = new TCanvas("c2");
 c2->SetFillColor(10);
 c2->SetGrid(); 
 c2->GetFrame()->SetFillColor(21);
 c2->GetFrame()->SetBorderSize(12);
 c2->cd();

 eta_dphi_prf1 = new TGraphErrors(55,xbin1a,prf1a,prfErrX1a,prfErrY1a);
 eta_dphi_prf1->SetTitle("eta_dphi_profile  eta 0.0 ~ 1.1");
 eta_dphi_prf1->SetMarkerColor(4);
 eta_dphi_prf1->SetMarkerStyle(21);

 //eta_dphi_prf1->GetXaxis()->SetTitle("Pt*dphi in different eta ranges");
 //eta_dphi_prf1->GetYaxis()->SetTitle(" ");
 eta_dphi_prf1->Draw("AP");
 eta_dphi_prf1->Print();

 c2->Update();
 c2->Print(plot02);

 gStyle->SetOptStat(kTRUE);
 gStyle->SetOptFit(0111);  
 c6 = new TCanvas("c6");
 c6->SetFillColor(10);
 c6->SetGrid();
 c6->GetFrame()->SetFillColor(21);
 c6->GetFrame()->SetBorderSize(12);
 c6->cd();

 pt_accuracy = new TGraph(55,xbin1a,pt_rt1);
 pt_accuracy->SetTitle("pT accuracy  eta 0.0 ~ 1.1");
 pt_accuracy->SetMarkerColor(4);
 pt_accuracy->SetMarkerStyle(21);
 pt_accuracy->Draw("AP");
 pt_accuracy->Print();

 c6->Update();
 c6->Print(plot06);


// ********************************************************************
// Draw the origin scalar plots
// ********************************************************************

 gStyle->SetOptStat(kTRUE);
 TCanvas *c5 = new TCanvas("c5","");
 c5->SetFillColor(10);
 c5->SetFillColor(10);
 heta_dphiPt->SetTitle(dphi_type);
 heta_dphiPt->Draw();
 heta_dphiPt->GetXaxis()->SetTitle("eta  ");
 heta_dphiPt->GetYaxis()->SetTitle("Pt x dphi  ");
 c5->Update();
 c5->Print(plot05);


 fclose(dfile);
 file->Close();
// gROOT->ProcessLine(".q");

}
