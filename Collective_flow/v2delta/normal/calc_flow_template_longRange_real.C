const double PI=3.1415926;
TH1D *hGlobalPeriph;
double gYIntegral = 0;

TH1D *sumRidge(TH2D *hin, TString name){
  TH1D *hout = new TH1D(name, "", hin->GetNbinsX(), hin->GetXaxis()->GetXmin(), hin->GetXaxis()->GetXmax()); 
	TH1D *htmp1 = (TH1D*)hin->ProjectionX("htmp1", hin->GetYaxis()->FindBin(-6.0), hin->GetYaxis()->FindBin(-5.0));
	TH1D *htmp2 = (TH1D*)hin->ProjectionX("htmp2", hin->GetYaxis()->FindBin(5.0), hin->GetYaxis()->FindBin(6.0));
  hout->Add(htmp1, htmp2);
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

double funTemplate(double *x, double *par){
  double dPhi = x[0];
  double F = par[0];
  double v2 = par[1];
  double v3 = par[2];
  //double G = par[4];
  static TF1 *funRidge = nullptr;
if (!funRidge) {
  funRidge = new TF1("funRidge", "[0]*(1+2*[1]*cos(2*x)+2*[2]*cos(3*x))", -PI*0.5, PI*1.5);
}
  //TF1 *funRidge = new TF1("funRidge", "[0]*(1+2*[1]*cos(2*x))", -PI*0.5, PI*1.5);
  double YIntegPeriph = hGlobalPeriph->Integral("width");
  double G = (gYIntegral-F*YIntegPeriph)/2./PI; //dDeltaPhi P.S. factor  considered
  funRidge->SetParameters(G, v2, v3);
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
  TF1 *fitFun = new TF1("fitFun", funTemplate, -0.5*PI, 1.5*PI, 3);
  fitFun->SetParameters(1, 0.1,-0.1);
  hYHigh->Fit("fitFun");
  err = fitFun->GetParError(1);
  v33 = fitFun->GetParameter(2);
  v33err = fitFun->GetParError(2);
  SavePlot(hYHigh, hYLow, hYHigh->GetName());
  return fitFun->GetParameter(1);
}

void SaveSumRidgePlots(TH1D *hSumRidge, TString name) {
  TCanvas c;
  hSumRidge->Draw();
  c.SaveAs(Form("sumridge_plots/%s.png", name.Data()));
}

void processFile(TString filename, TString label, int color, TGraphErrors** gr_v22_Nch, TGraphErrors** gr_v33_Nch, TGraphErrors** gr_v2_pt, TGraphErrors** gr_v3_pt) {
	TFile *file_1 = new TFile(filename);
	TFile *file_2 = new TFile(filename);

  //FMD12FMD3
  TH2D *hsameHigh_1FMD12FMD3 = (TH2D*)file_2->Get("hDEtaDPhiSameEventHighMidFMD12FMD3");
	TH2D *hmixedHigh_1FMD12FMD3 = (TH2D*)file_2->Get("hDEtaDPhiMixEventHighMidFMD12FMD3");
  TH1D *hTrigPtHigh_1FMD12FMD3 = (TH1D*)file_2->Get("hTrigPtHighFMD12FMD3");

	TH2D *hsameLow_1FMD12FMD3 = (TH2D*)file_1->Get("hDEtaDPhiSameEventLowMidFMD12FMD3");
	TH2D *hmixedLow_1FMD12FMD3 = (TH2D*)file_1->Get("hDEtaDPhiMixEventLowMidFMD12FMD3");
  TH1D *hTrigPtLow_1FMD12FMD3 = (TH1D*)file_1->Get("hTrigPtLowFMD12FMD3");

  TH1D *hSameLowRidgeFMD12FMD3 = sumRidge(hsameLow_1FMD12FMD3, "sameLowFMD12FMD3");
  TH1D *hMixLowRidgeFMD12FMD3 = sumRidge(hmixedLow_1FMD12FMD3, "mixLowFMD12FMD3");
  TH1D *hYLowFMD12FMD3 = getY(hSameLowRidgeFMD12FMD3, hMixLowRidgeFMD12FMD3, "YLowFMD12FMD3", hTrigPtLow_1FMD12FMD3->GetEntries());
  hYLowFMD12FMD3->Rebin(2);

  TH1D *hSameHighRidgeFMD12FMD3 = sumRidge(hsameHigh_1FMD12FMD3, "sameHighFMD12FMD3");
  TH1D *hMixHighRidgeFMD12FMD3 = sumRidge(hmixedHigh_1FMD12FMD3, "mixHighFMD12FMD3");
  TH1D *hYHighFMD12FMD3 = getY(hSameHighRidgeFMD12FMD3, hMixHighRidgeFMD12FMD3, "YHighFMD12FMD3", hTrigPtHigh_1FMD12FMD3->GetEntries());

  //FMD12FMD3
  double v33FMD12FMD3 = 0;
  double v33_errFMD12FMD3 = 0;
  double v22_errFMD12FMD3 = 0;
  //0.3-3  particle v22
  double v22FMD12FMD3 = getV22(hYHighFMD12FMD3, hYLowFMD12FMD3, v22_errFMD12FMD3, v33FMD12FMD3, v33_errFMD12FMD3);
  cout<<v22FMD12FMD3<<" "<<v22_errFMD12FMD3<<endl;

  //FMD12FMD3 - NCh dependence
  const int N_cent=8;
  double Nsel_min=0;
  double Nsel_bin=10;
  TH2D *hsame_centFMD12FMD3[N_cent];
  TH2D *hmixed_centFMD12FMD3[N_cent];
  TH1D *hTrigPt_centFMD12FMD3[N_cent];
  double Nsel_valFMD12FMD3[N_cent]={0};
  double Nsel_bin_widthFMD12FMD3[N_cent]={0};
  double v22_NchFMD12FMD3[N_cent]={0};
  double v22_Nch_errFMD12FMD3[N_cent]={0};
  double v33_NchFMD12FMD3[N_cent]={0};
  double v33_Nch_errFMD12FMD3[N_cent]={0};

  for(int icent=1; icent<N_cent; icent++){
    //FMD12FMD3
    Nsel_valFMD12FMD3[icent]=Nsel_min+icent*Nsel_bin+0.5*Nsel_bin;
    Nsel_bin_widthFMD12FMD3[icent]=0.5*Nsel_bin;
    if(Nsel_valFMD12FMD3[icent]<85){
      hsame_centFMD12FMD3[icent]=(TH2D*)file_1->Get(Form("hDEtaDPhiSameEvent_CentFMD12FMD3%d",icent));
      hmixed_centFMD12FMD3[icent]=(TH2D*)file_1->Get(Form("hDEtaDPhiMixEvent_CentFMD12FMD3%d",icent));
      hTrigPt_centFMD12FMD3[icent]=(TH1D*)file_1->Get(Form("hTrigPt_CentFMD12FMD3%d",icent));
    }
    else{
      hsame_centFMD12FMD3[icent]=(TH2D*)file_2->Get(Form("hDEtaDPhiSameEvent_CentFMD12FMD3%d",icent));
      hmixed_centFMD12FMD3[icent]=(TH2D*)file_2->Get(Form("hDEtaDPhiMixEvent_CentFMD12FMD3%d",icent));
      hTrigPt_centFMD12FMD3[icent]=(TH1D*)file_2->Get(Form("hTrigPt_CentFMD12FMD3%d",icent));
    }
    TH1D *hYHigh_centbinFMD12FMD3 = getYFromHist(hsame_centFMD12FMD3[icent], hmixed_centFMD12FMD3[icent], Form("High_centFMD12FMD3%d",icent), hTrigPt_centFMD12FMD3[icent]->GetEntries());
    hYHigh_centbinFMD12FMD3->Rebin(2);
    TH1D *hSumRidgeFMD12FMD3 = sumRidge(hsame_centFMD12FMD3[icent], Form("sumRidgeFMD12FMD3_cent%d", icent));
    SaveSumRidgePlots(hSumRidgeFMD12FMD3, Form("sumRidgeFMD12FMD3_cent%d", icent));
    TH1D *hmixSumRidgeFMD12FMD3 = sumRidge(hmixed_centFMD12FMD3[icent], Form("summixRidgeFMD12FMD3_cent%d", icent));
    SaveSumRidgePlots(hmixSumRidgeFMD12FMD3, Form("summixRidgeFMD12FMD3_cent%d", icent));
 
    double v33_centFMD12FMD3 = 0;
    double v33_cent_errFMD12FMD3 = 0;
    double v22_cent_errFMD12FMD3 = 0;
    double v22_centFMD12FMD3 = getV22(hYHigh_centbinFMD12FMD3, hYLowFMD12FMD3, v22_cent_errFMD12FMD3, v33_centFMD12FMD3, v33_cent_errFMD12FMD3);

    v22_NchFMD12FMD3[icent]=v22_centFMD12FMD3;
    v22_Nch_errFMD12FMD3[icent]=v22_cent_errFMD12FMD3;
    v33_NchFMD12FMD3[icent]=v33_centFMD12FMD3;
    v33_Nch_errFMD12FMD3[icent]=v33_cent_errFMD12FMD3;
  }
  
  *gr_v22_Nch = new TGraphErrors(N_cent, Nsel_valFMD12FMD3, v22_NchFMD12FMD3, Nsel_bin_widthFMD12FMD3, v22_Nch_errFMD12FMD3);
  *gr_v33_Nch = new TGraphErrors(N_cent, Nsel_valFMD12FMD3, v33_NchFMD12FMD3, Nsel_bin_widthFMD12FMD3, v33_Nch_errFMD12FMD3);
  
  (*gr_v22_Nch)->SetMarkerColor(color);
  (*gr_v22_Nch)->SetLineColor(color);
  (*gr_v33_Nch)->SetMarkerColor(color);
  (*gr_v33_Nch)->SetLineColor(color);

  // pT dependence
  const int N_ptbin=10;
  double ptbin_valFMD12FMD3[N_ptbin]={0};
  double ptbin_edgeFMD12FMD3[N_ptbin+1]={0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2.0, 2.5, 3.5, 4.5};
  double ptbin_widthFMD12FMD3[N_ptbin]={0};
  double v22_High_ptFMD12FMD3[N_ptbin]={0};
  double v22_High_pt_errFMD12FMD3[N_ptbin]={0};
  double v33_High_ptFMD12FMD3[N_ptbin]={0};
  double v33_High_pt_errFMD12FMD3[N_ptbin]={0};
  double v2_pt[N_ptbin]={0};
  double v2_pt_err[N_ptbin]={0};
  double v3_pt[N_ptbin]={0};
  double v3_pt_err[N_ptbin]={0};
  
  TH2D *hsameHigh_ptbinFMD12FMD3[N_ptbin];
  TH2D *hmixedHigh_ptbinFMD12FMD3[N_ptbin];
  TH1D *hTrigPtHigh_ptbinFMD12FMD3[N_ptbin];
  TH2D *hsameLow_ptbinFMD12FMD3[N_ptbin];
  TH2D *hmixedLow_ptbinFMD12FMD3[N_ptbin];
  TH1D *hTrigPtLow_ptbinFMD12FMD3[N_ptbin];
  
  for(int i=0; i<N_ptbin; i++){
    //FMD12FMD3
    ptbin_valFMD12FMD3[i]=(ptbin_edgeFMD12FMD3[i]+ptbin_edgeFMD12FMD3[i+1])*0.5;
    ptbin_widthFMD12FMD3[i]=(ptbin_edgeFMD12FMD3[i+1]-ptbin_edgeFMD12FMD3[i])*0.5;
    
    hsameHigh_ptbinFMD12FMD3[i]=(TH2D*)file_1->Get(Form("hDEtaDPhiSameEventHigh_ptbinFMD12FMD3%d",i));
    hmixedHigh_ptbinFMD12FMD3[i]=(TH2D*)file_1->Get(Form("hDEtaDPhiMixEventHigh_ptbinFMD12FMD3%d",i));
    hTrigPtHigh_ptbinFMD12FMD3[i]=(TH1D*)file_1->Get(Form("hTrigPtHigh_ptbinFMD12FMD3%d",i));
    
    if(hTrigPtHigh_ptbinFMD12FMD3[i] && hTrigPtHigh_ptbinFMD12FMD3[i]->GetEntries() > 0) {
      ptbin_valFMD12FMD3[i]=hTrigPtHigh_ptbinFMD12FMD3[i]->GetMean();
    }

    TH1D *hYHigh_ptbinFMD12FMD3 = getYFromHist(hsameHigh_ptbinFMD12FMD3[i], hmixedHigh_ptbinFMD12FMD3[i], Form("High_ptbinFMD12FMD3%d",i), hTrigPtHigh_ptbinFMD12FMD3[i]->GetEntries());

    hsameLow_ptbinFMD12FMD3[i]=(TH2D*)file_1->Get(Form("hDEtaDPhiSameEventLow_ptbinFMD12FMD3%d",i));
    hmixedLow_ptbinFMD12FMD3[i]=(TH2D*)file_1->Get(Form("hDEtaDPhiMixEventLow_ptbinFMD12FMD3%d",i));
    hTrigPtLow_ptbinFMD12FMD3[i]=(TH1D*)file_1->Get(Form("hTrigPtLow_ptbinFMD12FMD3%d",i));
    TH1D *hYLow_ptbinFMD12FMD3 = getYFromHist(hsameLow_ptbinFMD12FMD3[i], hmixedLow_ptbinFMD12FMD3[i], Form("Low_ptbinFMD12FMD3%d",i), hTrigPtLow_ptbinFMD12FMD3[i]->GetEntries());

    double v33FMD12FMD3 = 0;
    double v33_errFMD12FMD3 = 0;
    double v22_errFMD12FMD3 = 0;
    double v22FMD12FMD3 = getV22(hYHigh_ptbinFMD12FMD3, hYLow_ptbinFMD12FMD3, v22_errFMD12FMD3, v33FMD12FMD3, v33_errFMD12FMD3);

    v22_High_ptFMD12FMD3[i]=v22FMD12FMD3;
    v22_High_pt_errFMD12FMD3[i]=v22_errFMD12FMD3;
    v33_High_ptFMD12FMD3[i]=v33FMD12FMD3;
    v33_High_pt_errFMD12FMD3[i]=v33_errFMD12FMD3;

    // Simplified calculation - using the flow coefficients directly
    v2_pt[i] = v22_High_ptFMD12FMD3[i];
    v2_pt_err[i] = v22_High_pt_errFMD12FMD3[i];
    v3_pt[i] = v33_High_ptFMD12FMD3[i];
    v3_pt_err[i] = v33_High_pt_errFMD12FMD3[i];

    cout<<v2_pt[i]<<" "<<v2_pt_err[i]<<endl;
    cout<<v3_pt[i]<<" "<<v3_pt_err[i]<<endl;
  }

  *gr_v2_pt = new TGraphErrors(N_ptbin, ptbin_valFMD12FMD3, v2_pt, ptbin_widthFMD12FMD3, v2_pt_err);
  *gr_v3_pt = new TGraphErrors(N_ptbin, ptbin_valFMD12FMD3, v3_pt, ptbin_widthFMD12FMD3, v3_pt_err);
  
  (*gr_v2_pt)->SetMarkerColor(color);
  (*gr_v2_pt)->SetLineColor(color);
  (*gr_v3_pt)->SetMarkerColor(color);
  (*gr_v3_pt)->SetLineColor(color);
  
  file_1->Close();
  file_2->Close();
}

void calc_flow_template_longRange_real(){
  gStyle->SetOptStat(kFALSE);
  // Create output directories
  gSystem->mkdir("test_figure", kTRUE);
  gSystem->mkdir("test_figure/low_templates", kTRUE);
  gSystem->mkdir("sumridge_plots", kTRUE);
  
  // Define the three files and their properties
  TString filenames[4] = {
    "hist_ampt_normal_1.5mb_longRange_yuhao.root",
    "hist_ampt_normal_0.15mb_longRange_yuhao.root", 
    "hist_ampt_normal_1.5mb_a_0.8_b_0.4_longRange_yuhao.root",
    "hist_ampt_normal_0.15mb_a_0.8_b_0.4_longRange_yuhao.root"
  };
  
  TString labels[4] = {"normal_1.5mb", "normal_0.15mb", "normal_1.5mb_a0.8_b0.4","normal_0.15mb_a0.8_b0.4"};
  int colors[4] = {kBlack,kRed, kBlue, kGreen};
  int markers[4] = {20, 21, 22,23};
  
  // Arrays to store graphs for all files
  TGraphErrors *gr_v22_Nch[4];
  TGraphErrors *gr_v33_Nch[4];
  TGraphErrors *gr_v2_pt[4];
  TGraphErrors *gr_v3_pt[4];
  
  // Process all three files
  for(int ifile = 0; ifile < 4; ifile++) {
    cout << "\n=== Processing file " << ifile+1 << ": " << filenames[ifile] << " ===" << endl;
    processFile(filenames[ifile], labels[ifile], colors[ifile], 
                &gr_v22_Nch[ifile], &gr_v33_Nch[ifile], 
                &gr_v2_pt[ifile], &gr_v3_pt[ifile]);
    
    // Set marker styles
    gr_v22_Nch[ifile]->SetMarkerStyle(markers[ifile]);
    gr_v33_Nch[ifile]->SetMarkerStyle(markers[ifile]);
    gr_v2_pt[ifile]->SetMarkerStyle(markers[ifile]);
    gr_v3_pt[ifile]->SetMarkerStyle(markers[ifile]);
  }
  
  // Create output file and save individual graphs
  TFile *outfile = new TFile("longRange_flow_comparison.root", "recreate");
  
  for(int ifile = 0; ifile < 4; ifile++) {
    gr_v22_Nch[ifile]->Write(Form("gr_v22_Nch_%s", labels[ifile].Data()));
    gr_v33_Nch[ifile]->Write(Form("gr_v33_Nch_%s", labels[ifile].Data()));
    gr_v2_pt[ifile]->Write(Form("gr_v2_pt_%s", labels[ifile].Data()));
    gr_v3_pt[ifile]->Write(Form("gr_v3_pt_%s", labels[ifile].Data()));
  }
  
  // Create comparison plots
  
  // 1. v22 vs NCh comparison
  TCanvas *c_v22_Nch = new TCanvas("c_v22_Nch", "v22 vs Nch Comparison", 800, 600);
  
  TH2D *hframe_Nch = new TH2D("hframe_Nch", ";N_{ch};v_{2}{2}", 100, 0, 80, 100, 0, 0.0072);
  hframe_Nch->Draw();
  
  TLegend *leg_Nch = new TLegend(0.15, 0.4, 0.4, 0.7);
  leg_Nch->SetTextSize(0.030);
  leg_Nch->SetBorderSize(0);
  for(int ifile = 0; ifile < 4; ifile++) {
    gr_v22_Nch[ifile]->Draw("same p");
    leg_Nch->AddEntry(gr_v22_Nch[ifile], labels[ifile], "p");
  }
  leg_Nch->Draw();
  c_v22_Nch->SaveAs("v22_Nch_comparison.png");
  
  // 2. v33 vs NCh comparison
  //TCanvas *c_v33_Nch = new TCanvas("c_v33_Nch", "v33 vs Nch Comparison", 800, 600);
  
  //TH2D *hframe_v33_Nch = new TH2D("hframe_v33_Nch", ";N_{ch};v_{3}{3}", 100, 0, 80, 100, 0, 0.0032);
  //hframe_v33_Nch->Draw();
  
  // TLegend *leg_v33_Nch = new TLegend(0.1, 0.6, 0.3, 0.9);
  // for(int ifile = 0; ifile < 4; ifile++) {
  //   gr_v33_Nch[ifile]->Draw("same p");
  //   leg_v33_Nch->AddEntry(gr_v33_Nch[ifile], labels[ifile], "p");
  // }
  // leg_v33_Nch->Draw();
  //c_v33_Nch->SaveAs("v33_Nch_comparison.png");
  
  // 3. v2 vs pT comparison
  TCanvas *c_v2_pt = new TCanvas("c_v2_pt", "v2 vs pT Comparison", 800, 600);
  
  TH2D *hframe_pt = new TH2D("hframe_pt", ";p_{T} (GeV/c);v_{2}", 100, 0, 4, 100, 0, 0.0152);
  hframe_pt->Draw();
  TLegend *leg_pt = new TLegend(0.1, 0.6, 0.25, 0.9);
  leg_pt->SetTextSize(0.030);
  leg_pt->SetBorderSize(0);
  for(int ifile = 0; ifile < 4; ifile++) {
    gr_v2_pt[ifile]->Draw("same p");
    leg_pt->AddEntry(gr_v2_pt[ifile], labels[ifile], "p");
  }
  leg_pt->Draw();
  c_v2_pt->SaveAs("v2_pt_comparison.png");
  
  // 4. v3 vs pT comparison
  // TCanvas *c_v3_pt = new TCanvas("c_v3_pt", "v3 vs pT Comparison", 800, 600);
  // c_v3_pt->SetGridx();
  // c_v3_pt->SetGridy();
  
  //TH2D *hframe_v3_pt = new TH2D("hframe_v3_pt", ";p_{T} (GeV/c);v_{3}", 100, 0, 4, 100, 0, 0.0032);
  //hframe_v3_pt->Draw();
  
  // TLegend *leg_v3_pt = new TLegend(0.6, 0.7, 0.85, 0.85);
  // for(int ifile = 0; ifile < 4; ifile++) {
  //   gr_v3_pt[ifile]->Draw("same p");
  //   leg_v3_pt->AddEntry(gr_v3_pt[ifile], labels[ifile], "p");
  // }
  // leg_v3_pt->Draw();
  //c_v3_pt->SaveAs("v3_pt_comparison.png");
  
  // Save combined canvas with all four plots
  TCanvas *c_combined = new TCanvas("c_combined", "Flow Coefficients Comparison", 1600, 1200);
  c_combined->Divide(2, 2);
  c_combined->cd(1);
  hframe_Nch->Draw();
  for(int ifile = 0; ifile < 4; ifile++) {
    gr_v22_Nch[ifile]->Draw("same p");
  }
  leg_Nch->Draw();
  
  c_combined->cd(2);
  // hframe_v33_Nch->Draw();
  // for(int ifile = 0; ifile < 4; ifile++) {
  //   gr_v33_Nch[ifile]->Draw("same p");
  // }
  // leg_v33_Nch->Draw();
  
  c_combined->cd(3);
  hframe_pt->Draw();
  for(int ifile = 0; ifile < 4; ifile++) {
    gr_v2_pt[ifile]->Draw("same p");
  }
  leg_pt->Draw();
  
  c_combined->cd(4);
  //hframe_v3_pt->Draw();
  // for(int ifile = 0; ifile < 4; ifile++) {
  //   gr_v3_pt[ifile]->Draw("same p");
  // }
  // leg_v3_pt->Draw();
  
  c_combined->SaveAs("flow_comparison_all.png");
  
  outfile->Close();
  
  cout << "\n=== Analysis completed ===" << endl;
  cout << "Output files created:" << endl;
  cout << "- longRange_flow_comparison.root (contains all TGraphErrors)" << endl;
  cout << "- v22_Nch_comparison.png" << endl;
  cout << "- v33_Nch_comparison.png" << endl;
  cout << "- v2_pt_comparison.png" << endl;
  cout << "- v3_pt_comparison.png" << endl;
  cout << "- flow_comparison_all.png (combined plot)" << endl;
}