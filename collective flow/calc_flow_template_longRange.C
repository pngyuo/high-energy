const double PI=3.1415926;
TH1D *hGlobalPeriph;
double gYIntegral = 0;

void divideHist(TH2D *hin){
  for(int i=1; i<=hin->GetNbinsX(); ++i){
    for(int j=1; j<=hin->GetNbinsY(); ++j){
      double dphi = hin->GetXaxis()->GetBinCenter(i);
      double deta = hin->GetYaxis()->GetBinCenter(j);
      double ratio = (1-1e-4)*(4.75-fabs(deta))/4.75 +1e-4;
      double val = hin->GetBinContent(i,j);
      cout<<dphi<<" "<<deta<<" "<<ratio<<" "<<val<<endl;
      val = val/ratio;
      if(fabs(deta)>4.75) val=0;
      hin->SetBinContent(i, j, val);
    }
  }
}

void removeBkg(TH1D *hin){
	//double bkg = hin->GetBinContent(hin->FindBin(3.1416/2.));
	double bkg = hin->GetMinimum();
	for(int i=1; i<=hin->GetNbinsX(); i++){
		hin->SetBinContent(i, hin->GetBinContent(i)-bkg);
	}
}

TH1D *sumRidge(TH2D *hin, TString name){
  TH1D *hout = new TH1D(name, "", hin->GetNbinsX(), hin->GetXaxis()->GetXmin(), hin->GetXaxis()->GetXmax()); 
//	TH1D *htmp1 = (TH1D*)hin->ProjectionX("htmp1", 1, hin->GetYaxis()->FindBin(-2));
//	TH1D *htmp2 = (TH1D*)hin->ProjectionX("htmp2", hin->GetYaxis()->FindBin(2), hin->GetNbinsY());

	TH1D *htmp1 = (TH1D*)hin->ProjectionX("htmp1", hin->GetYaxis()->FindBin(-6.0), hin->GetYaxis()->FindBin(-3.0));
	TH1D *htmp2 = (TH1D*)hin->ProjectionX("htmp2", hin->GetYaxis()->FindBin(3.0), hin->GetYaxis()->FindBin(6.0));


  hout->Add(htmp1, htmp2);
  //get yield per deta bin, so divide by the |delta_eta| bin number, (-4~-2)+(2~4)=4
//  double nbins = hin->GetYaxis()->FindBin(4.0)-hin->GetYaxis()->FindBin(2.0);
//  nbins=nbins*2;
//  cout<<nbins<<endl;
//  hout->Scale(1./nbins);
  return hout;
}

TH1D *getY(TH1D *hSame, TH1D *hMix, TString name, double Ntrig){
  TH1D *hout = new TH1D(name, "", hSame->GetNbinsX(), hSame->GetXaxis()->GetXmin(), hSame->GetXaxis()->GetXmax()); 
  //get 1-D correlation function
  hout->Divide(hSame,hMix);
  double Btot=hMix->Integral();
  hout->Scale(1./Ntrig/2/3.14159*Btot);
  return hout;
}

TH1D *getYFromHist(TH2D *hSame, TH2D *hMix, TString name, double Ntrig){
  TH1D *hSameRidge = sumRidge(hSame, "same"+name);
  TH1D *hMixRidge = sumRidge(hMix, "mix"+name);
  TH1D *hY = getY(hSameRidge, hMixRidge, "Y"+name, Ntrig);
  return hY;
}



/*
double getJetYield(TH2D *hin){

  //divide by the bin width
  //hin->Scale(1, "width");

  TH1D *htmp_ridge = getRidge(hin, "htmp_ridge");
  TH1D *h_jet = (TH1D*)hin->ProjectionX("h_jet", hin->GetYaxis()->FindBin(-1.0),hin->GetYaxis()->FindBin(1.0));
  //cout<<hin->GetYaxis()->GetBinCenter(hin->GetYaxis()->FindBin(-1))<<" "<<hin->GetYaxis()->GetBinCenter(hin->GetYaxis()->FindBin(-1.2))<<endl;
//  TCanvas *ctmp = new TCanvas();
//  removeBkg(h_jet);
//  removeBkg(htmp_ridge);
  h_jet->Scale(1./2.0);//dEta range=2.0: (-1.0,1.0)
  htmp_ridge->Scale(1./4.);//dEta range=4: (-4,-2)+(2,4)
  htmp_ridge->SetMarkerColor(kRed);
//  h_jet->DrawClone();
//  htmp_ridge->DrawClone("same");

  h_jet->Add(htmp_ridge, -1);
//  removeBkg(h_jet);
  //TCanvas *ctmp = new TCanvas();
  //h_jet->DrawClone();

  double err_val=-999;
  //double YJet=h_jet->IntegralAndError(h_jet->GetXaxis()->FindBin(-1.2), h_jet->GetXaxis()->FindBin(1.2), err_val);
  double YJet=h_jet->Integral(h_jet->GetXaxis()->FindBin(-1.2), h_jet->GetXaxis()->FindBin(1.2));
  //scale the delta_phi width (-1.2,1.2)
  YJet/=2.4;

  //cout<<"YJet_val="<<YJet<<" err="<<err_val<<endl;

  delete htmp_ridge;
  delete h_jet;

  return YJet;
}
*/

double funTemplate(double *x, double *par){
  double dPhi = x[0];
  double F = par[0];
  double v2 = par[1];
  double v3 = par[2];
  double v4 = par[3];
  //double G = par[4];
  TF1 *funRidge = new TF1("funRidge", "[0]*(1+2*[1]*cos(2*x)+2*[2]*cos(3*x)+2*[3]*cos(4*x))", -PI*0.5, PI*1.5);
  //TF1 *funRidge = new TF1("funRidge", "[0]*(1+2*[1]*cos(2*x))", -PI*0.5, PI*1.5);
  double YIntegPeriph = hGlobalPeriph->Integral("width");
  double G = (gYIntegral-F*YIntegPeriph)/2./PI; //dDeltaPhi P.S. factor  considered
  funRidge->SetParameters(G, v2, v3, v4);
  //funRidge->SetParameters(G, v2);
  return funRidge->Eval(dPhi)+F*hGlobalPeriph->Interpolate(dPhi);
}


void SavePlot(TH1D *hYHigh, TH1D *hYLow, TString name){

  TCanvas c0;
  hYHigh->Draw();
  double par[4];
  for(int i=0; i<4; i++){
    par[i]=hYHigh->GetFunction("fitFun")->GetParameter(i);
  }
  hYHigh->GetFunction("fitFun")->GetParameters(par);
  double chi2 = hYHigh->GetFunction("fitFun")->GetChisquare();
  int ndf = hYHigh->GetFunction("fitFun")->GetNDF();
  TLatex *text = new TLatex();
  text->SetTextSize(0.04);
  text->DrawLatexNDC(0.25, 0.85, Form("F=%f",par[0]));
  text->DrawLatexNDC(0.25, 0.80, Form("v2,2=%f",par[1]));
  text->DrawLatexNDC(0.25, 0.75, Form("v3,3=%f",par[2]));
  text->DrawLatexNDC(0.25, 0.70, Form("v4,4=%f",par[3]));
  text->DrawLatexNDC(0.25, 0.65, Form("chi2=%.2f,ndf=%d",chi2,ndf));

  TF1 *funPeriph = new TF1("funPeriph", funTemplate, -0.5*PI, 1.5*PI, 4);
  funPeriph->SetParameters(par);
  funPeriph->SetParameter(1,0);
  funPeriph->SetParameter(2,0);
  funPeriph->SetParameter(3,0);
  funPeriph->SetLineColor(kRed);
  funPeriph->Draw("same");

  TF1 *funV2 = new TF1("funV2", "[0]*(2*[1]*cos(2*x))+[2]", -PI*0.5, PI*1.5);
  double F = par[0];
  double YIntegPeriph = hGlobalPeriph->Integral("width");
  double G = (gYIntegral-F*YIntegPeriph)/2./PI; //dDeltaPhi P.S. factor  considered
  funV2->SetParameter(0,G);
  funV2->SetParameter(1,par[1]);
  funV2->SetParameter(2,funPeriph->Eval(0));
  funV2->SetLineColor(kBlue);
  funV2->Draw("same");

  c0.SaveAs("test_figure/dphi_template_"+name+".png");

}



double getV22(TH1D *hYHigh, TH1D *hYLow, double &err, double &v33, double &v33err){
  gYIntegral=hYHigh->Integral("width");
  hGlobalPeriph = hYLow;
  TF1 *fitFun = new TF1("fitFun", funTemplate, -0.5*PI, 1.5*PI, 4);
  fitFun->SetParameters(1,1,1,1);
  hYHigh->Fit("fitFun");
  err = fitFun->GetParError(1);
  v33 = fitFun->GetParameter(2);
  v33err = fitFun->GetParError(2);

  SavePlot(hYHigh, hYLow, hYHigh->GetName());
  return fitFun->GetParameter(1);
}


void calc_flow_template_longRange(){
//	TFile *file_1 = new TFile("hist_Low0to20/hist_pp_13TeV_HIForced_minbias_0mb_geoNC3_MPT_flow.root");
//	TFile *file_2 = new TFile("hist_Low0to20/hist_pp_13TeV_HIForced_Trig200_0mb_geoNC3_MPT_flow.root");

	TFile *file_1 = new TFile("hist_root/hist_pp_13TeV_HIForced_minbias_0.15mb_geoNC3_MPT_flow_longRange.root");
	TFile *file_2 = new TFile("hist_root/hist_pp_13TeV_HIForced_minbias_0.15mb_geoNC3_MPT_flow_longRange.root");



	TH2D *hsame_1 = (TH2D*)file_1->Get("hDEtaDPhiSameEvent");
	TH2D *hmixed_1 = (TH2D*)file_1->Get("hDEtaDPhiMixEvent");

	TH2D *hsameHigh_1 = (TH2D*)file_2->Get("hDEtaDPhiSameEventHighMid");
	TH2D *hmixedHigh_1 = (TH2D*)file_2->Get("hDEtaDPhiMixEventHighMid");
  TH1D *hTrigPtHigh_1 = (TH1D*)file_2->Get("hTrigPtHigh");


	TH2D *hsameLow_1 = (TH2D*)file_1->Get("hDEtaDPhiSameEventLowMid");
	TH2D *hmixedLow_1 = (TH2D*)file_1->Get("hDEtaDPhiMixEventLowMid");
  TH1D *hTrigPtLow_1 = (TH1D*)file_1->Get("hTrigPtLow");


  TH1D *hSameLowRidge = sumRidge(hsameLow_1, "sameLow");
  TH1D *hMixLowRidge = sumRidge(hmixedLow_1, "mixLow");
  TH1D *hYLow = getY(hSameLowRidge, hMixLowRidge, "YLow", hTrigPtLow_1->GetEntries());
 
  TH1D *hSameHighRidge = sumRidge(hsameHigh_1, "sameHigh");
  TH1D *hMixHighRidge = sumRidge(hmixedHigh_1, "mixHigh");
  TH1D *hYHigh = getY(hSameHighRidge, hMixHighRidge, "YHigh", hTrigPtHigh_1->GetEntries());

  TCanvas *c_Y = new TCanvas();
  hYHigh->Draw();
  hYLow->Draw("same");

  double v33Ref = 0;
  double v33Ref_err = 0;
  double v22Ref_err = 0;
  //0.3-3 Ref particle v22
  double v22Ref = getV22(hYHigh, hYLow, v22Ref_err, v33Ref, v33Ref_err);
  cout<<v22Ref<<" "<<v22Ref_err<<endl;

//	TFile *infile_Pid_Low = new TFile("hist_Low0to20/hist_pp_13TeV_HIForced_minbias_0mb_geoNC3_MPT_flow.root");
//	TFile *infile_Pid_High = new TFile("hist_Low0to20/hist_pp_13TeV_HIForced_Trig200_0mb_geoNC3_MPT_flow.root");
	TFile *infile_Pid_Low = new TFile("hist_root/hist_pp_13TeV_HIForced_minbias_0.15mb_geoNC3_MPT_flow_longRange.root");
	TFile *infile_Pid_High = new TFile("hist_root/hist_pp_13TeV_HIForced_minbias_0.15mb_geoNC3_MPT_flow_longRange.root");



  const int N_ptbin=10;
  double pt_min=0;
  double pt_bin=0.5;
  double ptbin_val[N_ptbin]={0};
  double ptbin_edge[N_ptbin+1]={0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2.0, 2.5, 3.5, 4.5};
  double ptbin_width[N_ptbin]={0};
  double v2_pt[N_ptbin]={0};
  double v2_pt_err[N_ptbin]={0};
  double v3_pt[N_ptbin]={0};
  double v3_pt_err[N_ptbin]={0};
  double v22_High_pt[N_ptbin]={0};
  double v22_High_pt_err[N_ptbin]={0};
  TH2D *hsameHigh_ptbin[N_ptbin];
  TH2D *hmixedHigh_ptbin[N_ptbin];
  TH1D *hTrigPtHigh_ptbin[N_ptbin];
  TH2D *hsameLow_ptbin[N_ptbin];
  TH2D *hmixedLow_ptbin[N_ptbin];
  TH1D *hTrigPtLow_ptbin[N_ptbin];

  for(int i=0; i<N_ptbin; i++){
    ptbin_val[i]=(ptbin_edge[i]+ptbin_edge[i+1])*0.5;
    hsameHigh_ptbin[i]=(TH2D*)infile_Pid_High->Get(Form("hDEtaDPhiSameEventHigh_ptbin%d",i));
    hmixedHigh_ptbin[i]=(TH2D*)infile_Pid_High->Get(Form("hDEtaDPhiMixEventHigh_ptbin%d",i));
    hTrigPtHigh_ptbin[i]=(TH1D*)infile_Pid_High->Get(Form("hTrigPtHigh_ptbin%d",i));
    ptbin_val[i]=hTrigPtHigh_ptbin[i]->GetMean();

    TH1D *hYHigh_ptbin = getYFromHist(hsameHigh_ptbin[i], hmixedHigh_ptbin[i], Form("High_ptbin%d",i), hTrigPtHigh_ptbin[i]->GetEntries());

    hsameLow_ptbin[i]=(TH2D*)infile_Pid_Low->Get(Form("hDEtaDPhiSameEventLow_ptbin%d",i));
    hmixedLow_ptbin[i]=(TH2D*)infile_Pid_Low->Get(Form("hDEtaDPhiMixEventLow_ptbin%d",i));
    hTrigPtLow_ptbin[i]=(TH1D*)infile_Pid_Low->Get(Form("hTrigPtLow_ptbin%d",i));
    TH1D *hYLow_ptbin = getYFromHist(hsameLow_ptbin[i], hmixedLow_ptbin[i], Form("Low_ptbin%d",i), hTrigPtLow_ptbin[i]->GetEntries());


    double v33 = 0;
    double v33_err = 0;
    double v22_err = 0;
    double v22 = getV22(hYHigh_ptbin, hYLow_ptbin, v22_err, v33, v33_err);

    v2_pt[i] = v22/sqrt(v22Ref);
    v2_pt_err[i] = v2_pt[i]*sqrt( pow(v22_err/v22,2.)+ pow(0.5*v22Ref_err/v22Ref,2.) );

    v3_pt[i] = v33/sqrt(v33Ref);
    v3_pt_err[i] = v3_pt[i]*sqrt( pow(v33_err/v33,2.)+ pow(0.5*v33Ref_err/v33Ref,2.) );


    v22_High_pt[i]=v22;
    v22_High_pt_err[i]=v22_err;
  }

  TGraphErrors *gr_v2_pt = new TGraphErrors(N_ptbin, ptbin_val, v2_pt, ptbin_width, v2_pt_err);
  TGraphErrors *gr_v22_High_pt = new TGraphErrors(N_ptbin, ptbin_val, v22_High_pt, ptbin_width, v22_High_pt_err);
  TGraphErrors *gr_v3_pt = new TGraphErrors(N_ptbin, ptbin_val, v3_pt, ptbin_width, v3_pt_err);

  TH2D *hframe = new TH2D("hframe", ";p_{T}", 100, 0, 4, 100, -0.1, 0.3);
  TCanvas *c_v2 = new TCanvas();
  hframe->Draw();
  gr_v2_pt->SetMarkerColor(kRed);
  gr_v2_pt->Draw("same p");
  //gr_v22_High_pt->Draw("same p");

	TH2D *hsameLow_Pid = (TH2D*)infile_Pid_Low->Get("hDEtaDPhiSameEventLowMid");
	TH2D *hmixedLow_Pid = (TH2D*)infile_Pid_Low->Get("hDEtaDPhiMixEventLowMid");
  TH1D *hTrigPtLow_Pid = (TH1D*)infile_Pid_Low->Get("hTrigPtLow");
  TH1D *hYLow_Pid = getYFromHist(hsameLow_Pid, hmixedLow_Pid, "YLow_Pid", hTrigPtLow_Pid->GetEntries());


  const int N_cent=8;
  double Nsel_min=10;
  double Nsel_bin=10;
  TH2D *hsame_cent[N_cent];
  TH2D *hmixed_cent[N_cent];
  TH1D *hTrigPt_cent[N_cent];
  TH2D *hsame_cent_Pid[N_cent];
  TH2D *hmixed_cent_Pid[N_cent];
  TH1D *hTrigPt_cent_Pid[N_cent];
  double Nsel_val[N_cent]={0};
  double Nsel_bin_width[N_cent]={0};
  double v22Pid_Nch[N_cent]={0};
  double v22Pid_Nch_err[N_cent]={0};
  double v33Pid_Nch[N_cent]={0};
  double v33Pid_Nch_err[N_cent]={0};
  double v22Ref_Nch[N_cent]={0};
  double v22Ref_Nch_err[N_cent]={0};
  double v33Ref_Nch[N_cent]={0};
  double v33Ref_Nch_err[N_cent]={0};
  double v2_Nch[N_cent]={0};
  double v2_Nch_err[N_cent]={0};
  double v3_Nch[N_cent]={0};
  double v3_Nch_err[N_cent]={0};
  for(int icent=0; icent<N_cent; icent++){
    Nsel_val[icent]=Nsel_min+icent*Nsel_bin+0.5*Nsel_bin;
    Nsel_bin_width[icent]=0.5*Nsel_bin;
    if(Nsel_val[icent]<85){
      hsame_cent[icent]=(TH2D*)file_1->Get(Form("hDEtaDPhiSameEvent_Cent%d",icent));
      hmixed_cent[icent]=(TH2D*)file_1->Get(Form("hDEtaDPhiMixEvent_Cent%d",icent));
      hTrigPt_cent[icent]=(TH1D*)file_1->Get(Form("hTrigPt_Cent%d",icent));

      hsame_cent_Pid[icent]=(TH2D*)infile_Pid_Low->Get(Form("hDEtaDPhiSameEvent_Cent%d",icent));
      hmixed_cent_Pid[icent]=(TH2D*)infile_Pid_Low->Get(Form("hDEtaDPhiMixEvent_Cent%d",icent));
      hTrigPt_cent_Pid[icent]=(TH1D*)infile_Pid_Low->Get(Form("hTrigPt_Cent%d",icent));

    }
    else{
      hsame_cent[icent]=(TH2D*)file_2->Get(Form("hDEtaDPhiSameEvent_Cent%d",icent));
      hmixed_cent[icent]=(TH2D*)file_2->Get(Form("hDEtaDPhiMixEvent_Cent%d",icent));
      hTrigPt_cent[icent]=(TH1D*)file_2->Get(Form("hTrigPt_Cent%d",icent));

      hsame_cent_Pid[icent]=(TH2D*)infile_Pid_High->Get(Form("hDEtaDPhiSameEvent_Cent%d",icent));
      hmixed_cent_Pid[icent]=(TH2D*)infile_Pid_High->Get(Form("hDEtaDPhiMixEvent_Cent%d",icent));
      hTrigPt_cent_Pid[icent]=(TH1D*)infile_Pid_High->Get(Form("hTrigPt_Cent%d",icent));

    }
    TH1D *hYHigh_centbin = getYFromHist(hsame_cent[icent], hmixed_cent[icent], Form("High_cent%d",icent), hTrigPt_cent[icent]->GetEntries());

    double v33Ref_cent = 0;
    double v33Ref_cent_err = 0;
    double v22Ref_cent_err = 0;
    double v22Ref_cent = getV22(hYHigh_centbin, hYLow, v22Ref_cent_err, v33Ref_cent, v33Ref_cent_err);

    TH1D *hYHigh_centbin_Pid = getYFromHist(hsame_cent_Pid[icent], hmixed_cent_Pid[icent], Form("High_cent_Pid%d",icent), hTrigPt_cent_Pid[icent]->GetEntries());

    double v33Pid_cent = 0;
    double v33Pid_cent_err = 0;
    double v22Pid_cent_err = 0;
    double v22Pid_cent = getV22(hYHigh_centbin_Pid, hYLow_Pid, v22Pid_cent_err, v33Pid_cent, v33Pid_cent_err);

    v22Ref_Nch[icent]=v22Ref_cent;
    v22Ref_Nch_err[icent]=v22Ref_cent_err;
    v33Ref_Nch[icent]=v33Ref_cent;
    v33Ref_Nch_err[icent]=v33Ref_cent_err;

    v22Pid_Nch[icent]=v22Pid_cent;
    v22Pid_Nch_err[icent]=v22Pid_cent_err;
    v33Pid_Nch[icent]=v33Pid_cent;
    v33Pid_Nch_err[icent]=v33Pid_cent_err;

    v2_Nch[icent] = v22Pid_cent/sqrt(v22Ref_cent);
    v2_Nch_err[icent] = v2_Nch[icent]*sqrt( pow(v22Pid_cent_err/v22Pid_cent,2.)+ pow(0.5*v22Ref_cent_err/v22Ref_cent,2.) );

    v3_Nch[icent] = v33Pid_cent/sqrt(v33Ref_cent);
    v3_Nch_err[icent] = v3_Nch[icent]*sqrt( pow(v33Pid_cent_err/v33Pid_cent,2.)+ pow(0.5*v33Ref_cent_err/v33Ref_cent,2.) );
  }
 
  TGraphErrors *gr_v2_Nch = new TGraphErrors(N_cent, Nsel_val, v2_Nch, Nsel_bin_width, v2_Nch_err);
  TGraphErrors *gr_v3_Nch = new TGraphErrors(N_cent, Nsel_val, v3_Nch, Nsel_bin_width, v3_Nch_err);

  TGraphErrors *gr_v22_Pid_Nch = new TGraphErrors(N_cent, Nsel_val, v22Pid_Nch, Nsel_bin_width, v22Pid_Nch_err);
  TGraphErrors *gr_v33_Pid_Nch = new TGraphErrors(N_cent, Nsel_val, v33Pid_Nch, Nsel_bin_width, v33Pid_Nch_err);

  TGraphErrors *gr_v22_Ref_Nch = new TGraphErrors(N_cent, Nsel_val, v22Ref_Nch, Nsel_bin_width, v22Ref_Nch_err);
  TGraphErrors *gr_v33_Ref_Nch = new TGraphErrors(N_cent, Nsel_val, v33Ref_Nch, Nsel_bin_width, v33Ref_Nch_err);


  TCanvas *c_v2_Nch = new TCanvas();
  gr_v2_Nch->Draw("ape");

//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_Rope_charge.root", "recreate");
//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_Rope_PiCh.root", "recreate");
//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_Rope_KCh.root", "recreate");
//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_Rope_P.root", "recreate");
//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_Rope_Lam.root", "recreate");

//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_py8mpt_0.17mb_ntmax2_charge.root", "recreate");
//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_py8mpt_0.17mb_ntmax2_PiCh.root", "recreate");
//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_py8mpt_0.17mb_ntmax2_KCh.root", "recreate");
//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_py8mpt_0.17mb_ntmax2_P.root", "recreate");

//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_py8mpt_0.17mb_charge.root", "recreate");
//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_py8mpt_0.17mb_PiCh.root", "recreate");
//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_py8mpt_0.17mb_KCh.root", "recreate");
//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_py8mpt_0.17mb_P.root", "recreate");
//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_py8mpt_0.17mb_Lam.root", "recreate");

//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_py8mpt_0.1mb_charge.root", "recreate");
//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_py8mpt_0.15mb_charge.root", "recreate");
//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_py8mpt_0.15mb_PiCh.root", "recreate");
//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_py8mpt_0.15mb_KCh.root", "recreate");
//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_py8mpt_0.15mb_P.root", "recreate");
//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_std_charge.root", "recreate");

//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_py8mpt_0.15mb_p8hResc_charge.root", "recreate");
//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_py8mpt_0.15mb_p8hResc_PiCh.root", "recreate");
//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_py8mpt_0.15mb_p8hResc_KCh.root", "recreate");
//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_py8mpt_0.15mb_p8hResc_P.root", "recreate");

//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_py8mpt_0.15mb_p8decayOnly_charge.root", "recreate");
//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_py8mpt_0.15mb_p8decayOnly_PiCh.root", "recreate");
//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_py8mpt_0.15mb_p8decayOnly_KCh.root", "recreate");
//  TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_py8mpt_0.15mb_p8decayOnly_P.root", "recreate");

  //TFile *outfile = new TFile("flow_figures/Low0to20/pp_13TeV_v2_template_py8mpt_0mb_charge.root", "recreate");
  TFile *outfile = new TFile("longRange_flow.root", "recreate");

  gr_v2_pt->Write("gr_v2_pt");
  gr_v2_Nch->Write("gr_v2_Nch");
  gr_v3_pt->Write("gr_v3_pt");
  gr_v3_Nch->Write("gr_v3_Nch");

  gr_v22_Ref_Nch->Write("gr_v22_Ref_Nch");
  gr_v33_Ref_Nch->Write("gr_v33_Ref_Nch");
  gr_v22_Pid_Nch->Write("gr_v22_Pid_Nch");
  gr_v33_Pid_Nch->Write("gr_v33_Pid_Nch");

}
