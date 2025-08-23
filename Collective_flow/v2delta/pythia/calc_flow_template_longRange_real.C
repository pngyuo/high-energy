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
  fitFun->SetParameters(2, 0.03,0.01);
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

void calc_flow_template_longRange_real(){
	TFile *file_1 = new TFile("hist_outputallFSI.root");
	TFile *file_2 = new TFile("hist_outputallFSI.root");

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
  {
    TCanvas cLow("cLowFMD12FMD3", "Low Multiplicity Template - FMD12FMD3", 800, 600);
    hYLowFMD12FMD3->SetTitle("Low Multiplicity Template - FMD12FMD3; #Delta#phi; Y(#Delta#phi)");
    hYLowFMD12FMD3->SetLineColor(kGreen);
    hYLowFMD12FMD3->SetMarkerColor(kGreen);
    hYLowFMD12FMD3->SetMarkerStyle(22);
    hYLowFMD12FMD3->Draw("P");
    cLow.SaveAs("test_figure/low_templates/low_template_FMD12FMD3.png");
}
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

  //FMD12FMD3
  const int N_cent=8;
  double Nsel_min=0;
  double Nsel_bin=10;
  TH2D *hsame_centFMD12FMD3[N_cent];
  TH2D *hmixed_centFMD12FMD3[N_cent];
  TH1D *hTrigPt_centFMD12FMD3[N_cent];
  TH2D *hsame_cent_PidFMD12FMD3[N_cent];
  TH2D *hmixed_cent_PidFMD12FMD3[N_cent];
  TH1D *hTrigPt_cent_PidFMD12FMD3[N_cent];
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
  TGraphErrors *gr_v22_NchFMD12FMD3 = new TGraphErrors(N_cent, Nsel_valFMD12FMD3, v22_NchFMD12FMD3, Nsel_bin_widthFMD12FMD3, v22_Nch_errFMD12FMD3);
  TGraphErrors *gr_v33_NchFMD12FMD3 = new TGraphErrors(N_cent, Nsel_valFMD12FMD3, v33_NchFMD12FMD3, Nsel_bin_widthFMD12FMD3, v33_Nch_errFMD12FMD3);

  TFile *outfile = new TFile("longRange_flow.root", "recreate");

  gr_v22_NchFMD12FMD3->Write("gr_v22_NchFMD12FMD3");
  gr_v33_NchFMD12FMD3->Write("gr_v33_NchFMD12FMD3");
  
  // 保存 gr_v22_NchFMD12FMD3 为图片
  TCanvas *c_v22_NchFMD12FMD3 = new TCanvas("c_v22_NchFMD12FMD3", "v22 vs Nch FMD12FMD3", 800, 600);
  gr_v22_NchFMD12FMD3->Draw("ape");
  c_v22_NchFMD12FMD3->SaveAs("v22_NchFMD12FMD3.png");
  
  // 保存 gr_v33_NchFMD12FMD3 为图片
  // TCanvas *c_v33_NchFMD12FMD3 = new TCanvas("c_v33_NchFMD12FMD3", "v33 vs Nch FMD12FMD3", 800, 600);
  // gr_v33_NchFMD12FMD3->Draw("ape");
  // c_v33_NchFMD12FMD3->SaveAs("v33_NchFMD12FMD3.png");
}