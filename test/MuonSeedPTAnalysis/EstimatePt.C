void EstimatePt(){

  float  nsigmas = 1.5;

  TFile *file = TFile::Open("castor:/castor/cern.ch/user/d/dfortin/pt_estimate_5-200.root");
  
  TString suffixps = ".jpg";

  /******************************************************************
      Define the fitting function
  *******************************************************************/
  
  Double_t fitf(Double_t *x, Double_t *par) {
    
    Double_t fitval = (par[0]*x[0])+par[1]; 
    return fitval;
    
  }

  /******************************************************************
      Define Pointers to histograms
  *******************************************************************/
    
  TString name[16];
  name[0]  = "ME11_1";
  name[1]  = "ME11_2";
  name[2]  = "ME11_3";
  name[3]  = "ME11_4";
  name[4]  = "ME12_2";
  name[5]  = "ME12_3";
  name[6]  = "ME12_4";
  name[7]  = "ME2_3";
  name[8]  = "ME2_4";
  name[9]  = "ME3_4";
  name[10] = "MB1_2";
  name[11] = "MB1_3";
  name[12] = "MB1_4";
  name[13] = "MB2_3";
  name[14] = "MB2_4";
  name[15] = "MB3_4";


  /**********************************************
      Project  Pt dphi vs Eta in 120 bins of eta         
  ***********************************************/

  for ( int i = 13; i < 14; ++i) {
    
    TString nplot = name[i]+suffixps;
    hDphiPtvsEta = (TH2F *) file->Get(name[i]+"_hDPhiPt");

    double mean[100];
    double width[100];
    double errwidth[100];
    double eta[100];
    double err[100];

    // CSC: 60 bins from 0.9 to 2.4
    // DT:  55 bins from 0.0 to 1.1
    int nBins;
    float start;
    float end;
    float stepSize;
    // Dealing with CSC
    if ( i < 10 ) {
      nBins = 60;
      start = 0.9;
      end = 2.4;
      // Dealing with DT  
    } else {
      nBins = 55;
      start = 0.;
      end = 1.1;    
    }
    stepSize = (end - start) / nBins;
    for (int j=0; j < nBins; ++j) {
      err[j] = stepSize/2.;
      eta[j] = stepSize * j + stepSize/2. + start;
    }
  
    // Project histogram on pT * Dphi in slices of eta to fit mean pT * Dphi.
    for (int j=1; j<nBins+1; ++j) {
      hDphiPtvsEta->ProjectionY("hTemp",j,j+1);

      double nEntries = hTemp->GetEntries()-0.5;
      std::cout << "number of entries is " << nEntries << std::endl;

      // Check that have enough entries for fit
      if (nEntries < 50.) {
	mean[j] = 0.;
	width[j]= 0.;
	errwidth[j]= 0.;
      } else {
	// Perform 3 Gaussian fits
	// Initial on full range to get sigma
	// Following 2 fits using +/- x sigma range
	hTemp->Fit("gaus");
	TF1 *myfunc = hTemp->GetFunction("gaus");
	double par0 = myfunc->GetParameter(0);
	double par1 = myfunc->GetParameter(1);
	double par2 = myfunc->GetParameter(2);
	float Lo = par1 - nsigmas * par2;
	float Hi = par1 + nsigmas * par2;
	hTemp->Fit("gaus","R","",Lo,Hi);
	TF1 *myfunc = hTemp->GetFunction("gaus");
	par0 = myfunc->GetParameter(0);
	par1 = myfunc->GetParameter(1);
	par2 = myfunc->GetParameter(2);
	Lo = par1 - nsigmas * par2;
	Hi = par1 + nsigmas * par2;
	hTemp->Fit("gaus","R","",Lo,Hi);
	TF1 *myfunc = hTemp->GetFunction("gaus");
	par0 = myfunc->GetParameter(0);
	par1 = myfunc->GetParameter(1);
	par2 = myfunc->GetParameter(2);
	
	mean[j]    = par1;
	width[j]   = par2;
	errwidth[j]= par2/sqrt(nEntries); 
	
      }
    }
    

    int nGoodPoints = 0;
    float mean2[100];	 
    float width2[100];	 
    float errwidth2[100];
    float eta2[100];	 
    float err2[100];
    int k = 0;
    for (int j = 0; j < nBins; ++j) {
      if (width[j] > 0. && eta[j] > 0.41 && eta[j] < 0.61) {
	mean2[k]    = mean[j];	 
	width2[k]   = width[j];	 
	errwidth2[k]= errwidth[j];
	eta2[k]     = eta[j];	 
	err2[k]     = err[j];
	k++;
      }
    }
    int nGoodPoints = k;

    mean2[k]    = 0.;	 
    width2[k]   = 0.;	 
    errwidth2[k]= 0.;
    eta2[k]     = 2.5;
    if (i > 9) eta2[k]= 1.2;
    err2[k]     = 0.;

    std::cout << "Number of Good points is " << nGoodPoints << std::endl;
    if (nGoodPoints < 3) continue;

    gStyle->SetOptStat(kTRUE);
    gStyle->SetOptFit(0111);  
    c1 = new TCanvas("c1");
    c1->SetFillColor(10);
    c1->SetGrid(); 
    c1->GetFrame()->SetFillColor(21);
    c1->GetFrame()->SetBorderSize(12);
    c1->cd();
    
    h2 = new TGraphErrors(nGoodPoints+1,eta2,mean2,err2,width2);
    h2->SetTitle(name[i]);
    h2->SetMarkerColor(4);
    h2->SetMarkerStyle(21);
    
    TF1 *func = new TF1("fitf",fitf,eta2[0],eta2[nGoodPoints-1],2);
    h2->Fit("fitf","R","",eta2[0],eta2[nGoodPoints-1]);
    // h2->Fit("fitf","R","",1.2,2.4);
    
    double para0 = func->GetParameter(0);
    double para1 = func->GetParameter(1);
    double perr0 = func->GetParError(0);
    double perr1 = func->GetParError(1);
    
    h2->Draw("AP");
    h2->Print();
    
    c1->Update();
    c1->Print(nplot);
  
  }

  //  gROOT->ProcessLine(".q");

}
