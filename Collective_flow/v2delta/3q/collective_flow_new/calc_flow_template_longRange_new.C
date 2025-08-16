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

double funTemplate(double *x, double *par){
  double dPhi = x[0];
  double F = par[0];
  double v2 = par[1];
  double v3 = par[2];
  //double G = par[4];
  TF1 *funRidge = new TF1("funRidge", "[0]*(1+2*[1]*cos(2*x)+2*[2]*cos(3*x))", -PI*1, PI*2);
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
  fitFun->SetParameters(0.8,0.003,0.002);
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

void calc_flow_template_longRange_new(){
	TFile *file_1 = new TFile("hist_outputallFSI.root");
	TFile *file_2 = new TFile("hist_outputallFSI.root");
  //TPCFMD12
	TH2D *hsameHigh_1TPCFMD12 = (TH2D*)file_2->Get("hDEtaDPhiSameEventHighMidTPCFMD12");
	TH2D *hmixedHigh_1TPCFMD12 = (TH2D*)file_2->Get("hDEtaDPhiMixEventHighMidTPCFMD12");
  TH1D *hTrigPtHigh_1TPCFMD12 = (TH1D*)file_2->Get("hTrigPtHighTPCFMD12");

	TH2D *hsameLow_1TPCFMD12 = (TH2D*)file_1->Get("hDEtaDPhiSameEventLowMidTPCFMD12");
	TH2D *hmixedLow_1TPCFMD12 = (TH2D*)file_1->Get("hDEtaDPhiMixEventLowMidTPCFMD12");
  TH1D *hTrigPtLow_1TPCFMD12 = (TH1D*)file_1->Get("hTrigPtLowTPCFMD12");

  TH1D *hSameLowRidgeTPCFMD12 = sumRidge(hsameLow_1TPCFMD12, "sameLowTPCFMD12");
  TH1D *hMixLowRidgeTPCFMD12 = sumRidge(hmixedLow_1TPCFMD12, "mixLowTPCFMD12");
  TH1D *hYLowTPCFMD12 = getY(hSameLowRidgeTPCFMD12, hMixLowRidgeTPCFMD12, "YLowTPCFMD12", hTrigPtLow_1TPCFMD12->GetEntries());
 
  TH1D *hSameHighRidgeTPCFMD12 = sumRidge(hsameHigh_1TPCFMD12, "sameHighTPCFMD12");
  TH1D *hMixHighRidgeTPCFMD12 = sumRidge(hmixedHigh_1TPCFMD12, "mixHighTPCFMD12");
  TH1D *hYHighTPCFMD12 = getY(hSameHighRidgeTPCFMD12, hMixHighRidgeTPCFMD12, "YHighTPCFMD12", hTrigPtHigh_1TPCFMD12->GetEntries());

  //TPCFMD3
	TH2D *hsameHigh_1TPCFMD3 = (TH2D*)file_2->Get("hDEtaDPhiSameEventHighMidTPCFMD3");
	TH2D *hmixedHigh_1TPCFMD3 = (TH2D*)file_2->Get("hDEtaDPhiMixEventHighMidTPCFMD3");
  TH1D *hTrigPtHigh_1TPCFMD3 = (TH1D*)file_2->Get("hTrigPtHighTPCFMD3");

	TH2D *hsameLow_1TPCFMD3 = (TH2D*)file_1->Get("hDEtaDPhiSameEventLowMidTPCFMD3");
	TH2D *hmixedLow_1TPCFMD3 = (TH2D*)file_1->Get("hDEtaDPhiMixEventLowMidTPCFMD3");
  TH1D *hTrigPtLow_1TPCFMD3 = (TH1D*)file_1->Get("hTrigPtLowTPCFMD3");

  TH1D *hSameLowRidgeTPCFMD3 = sumRidge(hsameLow_1TPCFMD3, "sameLowTPCFMD3");
  TH1D *hMixLowRidgeTPCFMD3 = sumRidge(hmixedLow_1TPCFMD3, "mixLowTPCFMD3");
  TH1D *hYLowTPCFMD3 = getY(hSameLowRidgeTPCFMD3, hMixLowRidgeTPCFMD3, "YLow", hTrigPtLow_1TPCFMD3->GetEntries());
 
  TH1D *hSameHighRidgeTPCFMD3 = sumRidge(hsameHigh_1TPCFMD3, "sameHighTPCFMD3");
  TH1D *hMixHighRidgeTPCFMD3 = sumRidge(hmixedHigh_1TPCFMD3, "mixHighTPCFMD3");
  TH1D *hYHighTPCFMD3 = getY(hSameHighRidgeTPCFMD3, hMixHighRidgeTPCFMD3, "YHighTPCFMD3", hTrigPtHigh_1TPCFMD3->GetEntries());
  TCanvas c0;
  hYHighTPCFMD3->Draw();
  c0.SaveAs("test_figure/dphi_tem"".png");
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
 
  TH1D *hSameHighRidgeFMD12FMD3 = sumRidge(hsameHigh_1FMD12FMD3, "sameHighFMD12FMD3");
  TH1D *hMixHighRidgeFMD12FMD3 = sumRidge(hmixedHigh_1FMD12FMD3, "mixHighFMD12FMD3");
  TH1D *hYHighFMD12FMD3 = getY(hSameHighRidgeFMD12FMD3, hMixHighRidgeFMD12FMD3, "YHighFMD12FMD3", hTrigPtHigh_1FMD12FMD3->GetEntries());

  TCanvas *c_Y = new TCanvas();
  //TPCFMD12
  double v33TPCFMD12 = 0;
  double v33_errTPCFMD12 = 0;
  double v22_errTPCFMD12 = 0;
  //0.3-3  particle v22
  double v22TPCFMD12 = getV22(hYHighTPCFMD12, hYLowTPCFMD12, v22_errTPCFMD12, v33TPCFMD12, v33_errTPCFMD12);
  cout<<v22TPCFMD12<<" "<<v22_errTPCFMD12<<endl;
  //TPCFMD3
  double v33TPCFMD3 = 0;
  double v33_errTPCFMD3 = 0;
  double v22_errTPCFMD3 = 0;
  //0.3-3  particle v22
  double v22TPCFMD3 = getV22(hYHighTPCFMD3, hYLowTPCFMD3, v22_errTPCFMD3, v33TPCFMD3, v33_errTPCFMD3);
  cout<<v22TPCFMD3<<" "<<v22_errTPCFMD3<<endl;

  //FMD12FMD3
  double v33FMD12FMD3 = 0;
  double v33_errFMD12FMD3 = 0;
  double v22_errFMD12FMD3 = 0;
  //0.3-3  particle v22
  double v22FMD12FMD3 = getV22(hYHighFMD12FMD3, hYLowFMD12FMD3, v22_errFMD12FMD3, v33FMD12FMD3, v33_errFMD12FMD3);
  cout<<v22FMD12FMD3<<" "<<v22_errFMD12FMD3<<endl;

  if (v22FMD12FMD3 != 0) {
    double v2_3x2PC = (v22TPCFMD12 * v22TPCFMD3) / v22FMD12FMD3;
    double v2_err = v2_3x2PC* sqrt(
        pow(v22_errTPCFMD12 / v22TPCFMD12, 2) +
        pow(v22_errTPCFMD3 / v22TPCFMD3, 2) +
        pow(v22_errFMD12FMD3 / v22FMD12FMD3, 2)
    );
    // 输出结果
    cout << "v2(3x2PC) = " << v2_3x2PC << " ± " << v2_err << endl;
} else {
    cerr << "Error: Division by zero in 3x2PC formula." << endl;
}

	TFile *infile_Pid_Low = new TFile("hist_outputallFSI.root");
	TFile *infile_Pid_High = new TFile("hist_outputallFSI.root");
  const int N_ptbin=10;
  double v2_pt[N_ptbin]={0};
  double v2_pt_err[N_ptbin]={0};
  double v3_pt[N_ptbin]={0};
  double v3_pt_err[N_ptbin]={0};
  //TPCFMD12
  double ptbin_valTPCFMD12[N_ptbin]={0};
  double ptbin_edgeTPCFMD12[N_ptbin+1]={0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2.0, 2.5, 3.5, 4.5};
  double ptbin_widthTPCFMD12[N_ptbin]={0};
  double v22_High_ptTPCFMD12[N_ptbin]={0};
  double v22_High_pt_errTPCFMD12[N_ptbin]={0};
  double v33_High_ptTPCFMD12[N_ptbin]={0};
  double v33_High_pt_errTPCFMD12[N_ptbin]={0}; 
  TH2D *hsameHigh_ptbinTPCFMD12[N_ptbin];
  TH2D *hmixedHigh_ptbinTPCFMD12[N_ptbin];
  TH1D *hTrigPtHigh_ptbinTPCFMD12[N_ptbin];
  TH2D *hsameLow_ptbinTPCFMD12[N_ptbin];
  TH2D *hmixedLow_ptbinTPCFMD12[N_ptbin];
  TH1D *hTrigPtLow_ptbinTPCFMD12[N_ptbin];
  //TPCFMD3
  double ptbin_valTPCFMD3[N_ptbin]={0};
  double ptbin_edgeTPCFMD3[N_ptbin+1]={0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2.0, 2.5, 3.5, 4.5};
  double ptbin_widthTPCFMD3[N_ptbin]={0};
  double v22_High_ptTPCFMD3[N_ptbin]={0};
  double v22_High_pt_errTPCFMD3[N_ptbin]={0};
  double v33_High_ptTPCFMD3[N_ptbin]={0};
  double v33_High_pt_errTPCFMD3[N_ptbin]={0};
  TH2D *hsameHigh_ptbinTPCFMD3[N_ptbin];
  TH2D *hmixedHigh_ptbinTPCFMD3[N_ptbin];
  TH1D *hTrigPtHigh_ptbinTPCFMD3[N_ptbin];
  TH2D *hsameLow_ptbinTPCFMD3[N_ptbin];
  TH2D *hmixedLow_ptbinTPCFMD3[N_ptbin];
  TH1D *hTrigPtLow_ptbinTPCFMD3[N_ptbin];
  //FMD12FMD3
  double ptbin_valFMD12FMD3[N_ptbin]={0};
  double ptbin_edgeFMD12FMD3[N_ptbin+1]={0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2.0, 2.5, 3.5, 4.5};
  double ptbin_widthFMD12FMD3[N_ptbin]={0};
  double v22_High_ptFMD12FMD3[N_ptbin]={0};
  double v22_High_pt_errFMD12FMD3[N_ptbin]={0};
  double v33_High_ptFMD12FMD3[N_ptbin]={0};
  double v33_High_pt_errFMD12FMD3[N_ptbin]={0};
  TH2D *hsameHigh_ptbinFMD12FMD3[N_ptbin];
  TH2D *hmixedHigh_ptbinFMD12FMD3[N_ptbin];
  TH1D *hTrigPtHigh_ptbinFMD12FMD3[N_ptbin];
  TH2D *hsameLow_ptbinFMD12FMD3[N_ptbin];
  TH2D *hmixedLow_ptbinFMD12FMD3[N_ptbin];
  TH1D *hTrigPtLow_ptbinFMD12FMD3[N_ptbin];
  
  for(int i=0; i<N_ptbin; i++){
    //TPCFMD12
    ptbin_valTPCFMD12[i]=(ptbin_edgeTPCFMD12[i]+ptbin_edgeTPCFMD12[i+1])*0.5;
    hsameHigh_ptbinTPCFMD12[i]=(TH2D*)infile_Pid_High->Get(Form("hDEtaDPhiSameEventHigh_ptbinTPCFMD12%d",i));
    hmixedHigh_ptbinTPCFMD12[i]=(TH2D*)infile_Pid_High->Get(Form("hDEtaDPhiMixEventHigh_ptbinTPCFMD12%d",i));
    hTrigPtHigh_ptbinTPCFMD12[i]=(TH1D*)infile_Pid_High->Get(Form("hTrigPtHigh_ptbinTPCFMD12%d",i));
    ptbin_valTPCFMD12[i]=hTrigPtHigh_ptbinTPCFMD12[i]->GetMean();

    TH1D *hYHigh_ptbinTPCFMD12 = getYFromHist(hsameHigh_ptbinTPCFMD12[i], hmixedHigh_ptbinTPCFMD12[i], Form("High_ptbinTPCFMD12%d",i), hTrigPtHigh_ptbinTPCFMD12[i]->GetEntries());

    hsameLow_ptbinTPCFMD12[i]=(TH2D*)infile_Pid_Low->Get(Form("hDEtaDPhiSameEventLow_ptbinTPCFMD12%d",i));
    hmixedLow_ptbinTPCFMD12[i]=(TH2D*)infile_Pid_Low->Get(Form("hDEtaDPhiMixEventLow_ptbinTPCFMD12%d",i));
    hTrigPtLow_ptbinTPCFMD12[i]=(TH1D*)infile_Pid_Low->Get(Form("hTrigPtLow_ptbinTPCFMD12%d",i));
    TH1D *hYLow_ptbinTPCFMD12 = getYFromHist(hsameLow_ptbinTPCFMD12[i], hmixedLow_ptbinTPCFMD12[i], Form("Low_ptbinTPCFMD12%d",i), hTrigPtLow_ptbinTPCFMD12[i]->GetEntries());

    double v33TPCFMD12 = 0;
    double v33_errTPCFMD12 = 0;
    double v22_errTPCFMD12 = 0;
    double v22TPCFMD12 = getV22(hYHigh_ptbinTPCFMD12, hYLow_ptbinTPCFMD12, v22_errTPCFMD12, v33TPCFMD12, v33_errTPCFMD12);

    v22_High_ptTPCFMD12[i]=v22TPCFMD12;
    v22_High_pt_errTPCFMD12[i]=v22_errTPCFMD12;
    v33_High_ptTPCFMD12[i]=v33TPCFMD12;
    v33_High_pt_errTPCFMD12[i]=v33_errTPCFMD12;

    //TPCFMD3
    ptbin_valTPCFMD3[i]=(ptbin_edgeTPCFMD3[i]+ptbin_edgeTPCFMD3[i+1])*0.5;
    hsameHigh_ptbinTPCFMD3[i]=(TH2D*)infile_Pid_High->Get(Form("hDEtaDPhiSameEventHigh_ptbinTPCFMD3%d",i));
    hmixedHigh_ptbinTPCFMD3[i]=(TH2D*)infile_Pid_High->Get(Form("hDEtaDPhiMixEventHigh_ptbinTPCFMD3%d",i));
    hTrigPtHigh_ptbinTPCFMD3[i]=(TH1D*)infile_Pid_High->Get(Form("hTrigPtHigh_ptbinTPCFMD3%d",i));
    ptbin_valTPCFMD3[i]=hTrigPtHigh_ptbinTPCFMD3[i]->GetMean();

    TH1D *hYHigh_ptbinTPCFMD3 = getYFromHist(hsameHigh_ptbinTPCFMD3[i], hmixedHigh_ptbinTPCFMD3[i], Form("High_ptbinTPCFMD3%d",i), hTrigPtHigh_ptbinTPCFMD3[i]->GetEntries());

    hsameLow_ptbinTPCFMD3[i]=(TH2D*)infile_Pid_Low->Get(Form("hDEtaDPhiSameEventLow_ptbinTPCFMD3%d",i));
    hmixedLow_ptbinTPCFMD3[i]=(TH2D*)infile_Pid_Low->Get(Form("hDEtaDPhiMixEventLow_ptbinTPCFMD3%d",i));
    hTrigPtLow_ptbinTPCFMD3[i]=(TH1D*)infile_Pid_Low->Get(Form("hTrigPtLow_ptbinTPCFMD3%d",i));
    TH1D *hYLow_ptbinTPCFMD3 = getYFromHist(hsameLow_ptbinTPCFMD3[i], hmixedLow_ptbinTPCFMD3[i], Form("Low_ptbin%d",i), hTrigPtLow_ptbinTPCFMD3[i]->GetEntries());

    double v33TPCFMD3 = 0;
    double v33_errTPCFMD3 = 0;
    double v22_errTPCFMD3 = 0;
    double v22TPCFMD3 = getV22(hYHigh_ptbinTPCFMD3, hYLow_ptbinTPCFMD3, v22_errTPCFMD3, v33TPCFMD3, v33_errTPCFMD3);

    v22_High_ptTPCFMD3[i]=v22TPCFMD3;
    v22_High_pt_errTPCFMD3[i]=v22_errTPCFMD3;
    v33_High_ptTPCFMD3[i]=v33TPCFMD3;
    v33_High_pt_errTPCFMD3[i]=v33_errTPCFMD3;

    //FMD12FMD3
    ptbin_valFMD12FMD3[i]=(ptbin_edgeFMD12FMD3[i]+ptbin_edgeFMD12FMD3[i+1])*0.5;
    hsameHigh_ptbinFMD12FMD3[i]=(TH2D*)infile_Pid_High->Get(Form("hDEtaDPhiSameEventHigh_ptbinFMD12FMD3%d",i));
    hmixedHigh_ptbinFMD12FMD3[i]=(TH2D*)infile_Pid_High->Get(Form("hDEtaDPhiMixEventHigh_ptbinFMD12FMD3%d",i));
    hTrigPtHigh_ptbinFMD12FMD3[i]=(TH1D*)infile_Pid_High->Get(Form("hTrigPtHigh_ptbinFMD12FMD3%d",i));
    ptbin_valFMD12FMD3[i]=hTrigPtHigh_ptbinFMD12FMD3[i]->GetMean();

    TH1D *hYHigh_ptbinFMD12FMD3 = getYFromHist(hsameHigh_ptbinFMD12FMD3[i], hmixedHigh_ptbinFMD12FMD3[i], Form("High_ptbinFMD12FMD3%d",i), hTrigPtHigh_ptbinFMD12FMD3[i]->GetEntries());

    hsameLow_ptbinFMD12FMD3[i]=(TH2D*)infile_Pid_Low->Get(Form("hDEtaDPhiSameEventLow_ptbinFMD12FMD3%d",i));
    hmixedLow_ptbinFMD12FMD3[i]=(TH2D*)infile_Pid_Low->Get(Form("hDEtaDPhiMixEventLow_ptbinFMD12FMD3%d",i));
    hTrigPtLow_ptbinFMD12FMD3[i]=(TH1D*)infile_Pid_Low->Get(Form("hTrigPtLow_ptbinFMD12FMD3%d",i));
    TH1D *hYLow_ptbinFMD12FMD3 = getYFromHist(hsameLow_ptbinFMD12FMD3[i], hmixedLow_ptbinFMD12FMD3[i], Form("Low_ptbinFMD12FMD3%d",i), hTrigPtLow_ptbinFMD12FMD3[i]->GetEntries());

    double v33FMD12FMD3 = 0;
    double v33_errFMD12FMD3 = 0;
    double v22_errFMD12FMD3 = 0;
    double v22FMD12FMD3 = getV22(hYHigh_ptbinFMD12FMD3, hYLow_ptbinFMD12FMD3, v22_errFMD12FMD3, v33FMD12FMD3, v33_errFMD12FMD3);

    v22_High_ptFMD12FMD3[i]=v22FMD12FMD3;
    v22_High_pt_errFMD12FMD3[i]=v22_errFMD12FMD3;
    v33_High_ptFMD12FMD3[i]=v33FMD12FMD3;
    v33_High_pt_errFMD12FMD3[i]=v33_errFMD12FMD3;

    if (v22_High_ptFMD12FMD3[i] != 0) {
      v2_pt[i] = (v22_High_ptTPCFMD12[i] * v22_High_ptTPCFMD3[i]) / v22_High_ptFMD12FMD3[i];
      v2_pt_err[i] = v2_pt[i] * sqrt(
          pow(v22_High_pt_errTPCFMD12[i] / v22_High_ptTPCFMD12[i], 2) +
          pow(v22_High_pt_errTPCFMD3[i] / v22_High_ptTPCFMD3[i], 2) +
          pow(v22_High_pt_errFMD12FMD3[i] / v22_High_ptFMD12FMD3[i], 2)
      );
  } else {
      v2_pt[i] = 0;
      v2_pt_err[i] = 0;
  }
  
  if (v33_High_ptFMD12FMD3[i] != 0) {
      v3_pt[i] = (v33_High_ptTPCFMD12[i] * v33_High_ptTPCFMD3[i]) / v33_High_ptFMD12FMD3[i];
      v3_pt_err[i] = v3_pt[i] * sqrt(
          pow(v33_High_pt_errTPCFMD12[i] / v33_High_ptTPCFMD12[i], 2) +
          pow(v33_High_pt_errTPCFMD3[i] / v33_High_ptTPCFMD3[i], 2) +
          pow(v33_High_pt_errFMD12FMD3[i] / v33_High_ptFMD12FMD3[i], 2)
      );
  } else {
      v3_pt[i] = 0;
      v3_pt_err[i] = 0;
  }
  cout<<v2_pt[i]<<" "<<v2_pt_err[i]<<endl;
  cout<<v3_pt[i]<<" "<<v3_pt_err[i]<<endl;
  }

  TGraphErrors *gr_v2_pt = new TGraphErrors(N_ptbin, ptbin_valFMD12FMD3, v2_pt, ptbin_widthFMD12FMD3, v2_pt_err);
  TGraphErrors *gr_v3_pt = new TGraphErrors(N_ptbin, ptbin_valFMD12FMD3, v3_pt, ptbin_widthFMD12FMD3, v3_pt_err);

  TH2D *hframe = new TH2D("hframe", ";p_{T}", 100, 0, 4, 100, -0.1, 0.3);
  TCanvas *c_v2 = new TCanvas();
  hframe->Draw();
  gr_v2_pt->SetMarkerColor(kRed);
  gr_v2_pt->Draw("same p");
  //gr_v22_High_pt->Draw("same p");

  //TPCFMD12
  const int N_cent=8;
  double Nsel_min=0;
  double Nsel_bin=10;
  TH2D *hsame_centTPCFMD12[N_cent];
  TH2D *hmixed_centTPCFMD12[N_cent];
  TH1D *hTrigPt_centTPCFMD12[N_cent];
  TH2D *hsame_cent_PidTPCFMD12[N_cent];
  TH2D *hmixed_cent_PidTPCFMD12[N_cent];
  TH1D *hTrigPt_cent_PidTPCFMD12[N_cent];
  double Nsel_valTPCFMD12[N_cent]={0};
  double Nsel_bin_widthTPCFMD12[N_cent]={0};
  double v22_NchTPCFMD12[N_cent]={0};
  double v22_Nch_errTPCFMD12[N_cent]={0};
  double v33_NchTPCFMD12[N_cent]={0};
  double v33_Nch_errTPCFMD12[N_cent]={0};
  //TPCFMD3
  TH2D *hsame_centTPCFMD3[N_cent];
  TH2D *hmixed_centTPCFMD3[N_cent];
  TH1D *hTrigPt_centTPCFMD3[N_cent];
  TH2D *hsame_cent_PidTPCFMD3[N_cent];
  TH2D *hmixed_cent_PidTPCFMD3[N_cent];
  TH1D *hTrigPt_cent_PidTPCFMD3[N_cent];
  double Nsel_valTPCFMD3[N_cent]={0};
  double Nsel_bin_widthTPCFMD3[N_cent]={0};
  double v22_NchTPCFMD3[N_cent]={0};
  double v22_Nch_errTPCFMD3[N_cent]={0};
  double v33_NchTPCFMD3[N_cent]={0};
  double v33_Nch_errTPCFMD3[N_cent]={0};
  //FMD12FMD3
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

  double v2_Nch[N_cent]={0};
  double v2_Nch_err[N_cent]={0};
  double v3_Nch[N_cent]={0};
  double v3_Nch_err[N_cent]={0};

  for(int icent=1; icent<N_cent; icent++){
    //TPCFMD12
    Nsel_valTPCFMD12[icent]=Nsel_min+icent*Nsel_bin+0.5*Nsel_bin;
    Nsel_bin_widthTPCFMD12[icent]=0.5*Nsel_bin;
    if(Nsel_valTPCFMD12[icent]<85){
      hsame_centTPCFMD12[icent]=(TH2D*)file_1->Get(Form("hDEtaDPhiSameEvent_CentTPCFMD12%d",icent));
      hmixed_centTPCFMD12[icent]=(TH2D*)file_1->Get(Form("hDEtaDPhiMixEvent_CentTPCFMD12%d",icent));
      hTrigPt_centTPCFMD12[icent]=(TH1D*)file_1->Get(Form("hTrigPt_CentTPCFMD12%d",icent));
    }
    else{
      hsame_centTPCFMD12[icent]=(TH2D*)file_2->Get(Form("hDEtaDPhiSameEvent_CentTPCFMD12%d",icent));
      hmixed_centTPCFMD12[icent]=(TH2D*)file_2->Get(Form("hDEtaDPhiMixEvent_CentTPCFMD12%d",icent));
      hTrigPt_centTPCFMD12[icent]=(TH1D*)file_2->Get(Form("hTrigPt_CentTPCFMD12%d",icent));
    }
    TH1D *hYHigh_centbinTPCFMD12 = getYFromHist(hsame_centTPCFMD12[icent], hmixed_centTPCFMD12[icent], Form("High_centTPCFMD12%d",icent), hTrigPt_centTPCFMD12[icent]->GetEntries());
    TH1D *hSumRidgeTPCFMD12 = sumRidge(hsame_centTPCFMD12[icent], Form("sumRidgeTPCFMD12_cent%d", icent));
    SaveSumRidgePlots(hSumRidgeTPCFMD12, Form("sumRidgeTPCFMD12_cent%d", icent));
    TH1D *hmixSumRidgeTPCFMD12 = sumRidge(hmixed_centTPCFMD12[icent], Form("summixRidgeTPCFMD12_cent%d", icent));
    SaveSumRidgePlots(hmixSumRidgeTPCFMD12, Form("summixRidgeTPCFMD12_cent%d", icent));
 
    double v33_centTPCFMD12 = 0;
    double v33_cent_errTPCFMD12 = 0;
    double v22_cent_errTPCFMD12 = 0;
    double v22_centTPCFMD12 = getV22(hYHigh_centbinTPCFMD12, hYLowTPCFMD12, v22_cent_errTPCFMD12, v33_centTPCFMD12, v33_cent_errTPCFMD12);

    v22_NchTPCFMD12[icent]=v22_centTPCFMD12;
    v22_Nch_errTPCFMD12[icent]=v22_cent_errTPCFMD12;
    v33_NchTPCFMD12[icent]=v33_centTPCFMD12;
    v33_Nch_errTPCFMD12[icent]=v33_cent_errTPCFMD12;

    //TPCFMD3
    Nsel_valTPCFMD3[icent]=Nsel_min+icent*Nsel_bin+0.5*Nsel_bin;
    Nsel_bin_widthTPCFMD3[icent]=0.5*Nsel_bin;
    if(Nsel_valTPCFMD3[icent]<85){
      hsame_centTPCFMD3[icent]=(TH2D*)file_1->Get(Form("hDEtaDPhiSameEvent_CentTPCFMD3%d",icent));
      hmixed_centTPCFMD3[icent]=(TH2D*)file_1->Get(Form("hDEtaDPhiMixEvent_CentTPCFMD3%d",icent));
      hTrigPt_centTPCFMD3[icent]=(TH1D*)file_1->Get(Form("hTrigPt_CentTPCFMD3%d",icent));
    }
    else{
      hsame_centTPCFMD3[icent]=(TH2D*)file_2->Get(Form("hDEtaDPhiSameEvent_CentTPCFMD3%d",icent));
      hmixed_centTPCFMD3[icent]=(TH2D*)file_2->Get(Form("hDEtaDPhiMixEvent_CentTPCFMD3%d",icent));
      hTrigPt_centTPCFMD3[icent]=(TH1D*)file_2->Get(Form("hTrigPt_CentTPCFMD3%d",icent));
    }
    TH1D *hYHigh_centbinTPCFMD3 = getYFromHist(hsame_centTPCFMD3[icent], hmixed_centTPCFMD3[icent], Form("High_centTPCFMD3%d",icent), hTrigPt_centTPCFMD3[icent]->GetEntries());
    TH1D *hSumRidgeTPCFMD3 = sumRidge(hsame_centTPCFMD3[icent], Form("sumRidgeTPCFMD3_cent%d", icent));
    SaveSumRidgePlots(hSumRidgeTPCFMD3, Form("sumRidgeTPCFMD3_cent%d", icent));
    TH1D *hmixSumRidgeTPCFMD3 = sumRidge(hmixed_centTPCFMD3[icent], Form("summixRidgeTPCFMD3_cent%d", icent));
    SaveSumRidgePlots(hmixSumRidgeTPCFMD3, Form("summixRidgeTPCFMD3_cent%d", icent));
    
    double v33_centTPCFMD3 = 0;
    double v33_cent_errTPCFMD3 = 0;
    double v22_cent_errTPCFMD3 = 0;
    double v22_centTPCFMD3 = getV22(hYHigh_centbinTPCFMD3, hYLowTPCFMD3, v22_cent_errTPCFMD3, v33_centTPCFMD3, v33_cent_errTPCFMD3);

    v22_NchTPCFMD3[icent]=v22_centTPCFMD3;
    v22_Nch_errTPCFMD3[icent]=v22_cent_errTPCFMD3;
    v33_NchTPCFMD3[icent]=v33_centTPCFMD3;
    v33_Nch_errTPCFMD3[icent]=v33_cent_errTPCFMD3;

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

    if (v22_NchFMD12FMD3[icent] != 0) {
      v2_Nch[icent] = (v22_NchTPCFMD12[icent] * v22_NchTPCFMD3[icent]) / v22_NchFMD12FMD3[icent];
      v2_Nch_err[icent] = v2_Nch[icent] * sqrt(
          pow(v22_Nch_errTPCFMD12[icent] / v22_NchTPCFMD12[icent], 2) +
          pow(v22_Nch_errTPCFMD3[icent]/ v22_NchTPCFMD3[icent], 2) +
          pow(v22_Nch_errFMD12FMD3[icent] / v22_NchFMD12FMD3[icent], 2)
      );
  } else {
      v2_Nch[icent] = 0;
      v2_Nch_err[icent] = 0;
  }
  
  if (v33_NchFMD12FMD3[icent] != 0) {
      v3_Nch[icent] = (v33_NchTPCFMD12[icent] * v33_NchTPCFMD3[icent]) / v33_NchFMD12FMD3[icent];
      v3_Nch_err[icent] = v3_Nch[icent] * sqrt(
          pow(v33_Nch_errTPCFMD12[icent] / v33_NchTPCFMD12[icent], 2) +
          pow(v33_Nch_errTPCFMD3[icent]/ v33_NchTPCFMD3[icent], 2) +
          pow(v33_Nch_errFMD12FMD3[icent] / v33_NchFMD12FMD3[icent], 2)
      );
  } else {
      v3_Nch[icent] = 0;
      v3_Nch_err[icent] = 0;
  }
  cout<<v2_Nch[icent]<<" "<<v2_Nch_err[icent]<<endl;
  cout<<v3_Nch[icent]<<" "<<v3_Nch_err[icent]<<endl;
  }
  TGraphErrors *gr_v2_Nch = new TGraphErrors(N_cent, Nsel_valFMD12FMD3, v2_Nch, Nsel_bin_widthFMD12FMD3, v2_Nch_err);
  TGraphErrors *gr_v3_Nch = new TGraphErrors(N_cent, Nsel_valFMD12FMD3, v3_Nch, Nsel_bin_widthFMD12FMD3, v3_Nch_err);

  TGraphErrors *gr_v22_NchTPCFMD12 = new TGraphErrors(N_cent, Nsel_valTPCFMD12, v22_NchTPCFMD12, Nsel_bin_widthTPCFMD12, v22_Nch_errTPCFMD12);
  TGraphErrors *gr_v33_NchTPCFMD12 = new TGraphErrors(N_cent, Nsel_valTPCFMD12, v33_NchTPCFMD12, Nsel_bin_widthTPCFMD12, v33_Nch_errTPCFMD12);
  TGraphErrors *gr_v22_NchTPCFMD3 = new TGraphErrors(N_cent, Nsel_valTPCFMD3, v22_NchTPCFMD3, Nsel_bin_widthTPCFMD3, v22_Nch_errTPCFMD3);
  TGraphErrors *gr_v33_NchTPCFMD3 = new TGraphErrors(N_cent, Nsel_valTPCFMD3, v33_NchTPCFMD3, Nsel_bin_widthTPCFMD3, v33_Nch_errTPCFMD3);
  TGraphErrors *gr_v22_NchFMD12FMD3 = new TGraphErrors(N_cent, Nsel_valFMD12FMD3, v22_NchFMD12FMD3, Nsel_bin_widthFMD12FMD3, v22_Nch_errFMD12FMD3);
  TGraphErrors *gr_v33_NchFMD12FMD3 = new TGraphErrors(N_cent, Nsel_valFMD12FMD3, v33_NchFMD12FMD3, Nsel_bin_widthFMD12FMD3, v33_Nch_errFMD12FMD3);


  // TCanvas *c_v2_Nch = new TCanvas();
  // gr_v2_Nch->Draw("ape");

  TFile *outfile = new TFile("longRange_flow.root", "recreate");

  gr_v2_pt->Write("gr_v2_pt");
  gr_v2_Nch->Write("gr_v2_Nch");
  gr_v3_pt->Write("gr_v3_pt");
  gr_v3_Nch->Write("gr_v3_Nch");

  gr_v22_NchTPCFMD12->Write("gr_v22_NchTPCFMD12");
  gr_v33_NchTPCFMD12->Write("gr_v33_NchTPCFMD12");
  gr_v22_NchTPCFMD3->Write("gr_v22_NchTPCFMD3");
  gr_v33_NchTPCFMD3->Write("gr_v33_NchTPCFMD3");
  gr_v22_NchFMD12FMD3->Write("gr_v22_NchFMD12FMD3");
  gr_v33_NchFMD12FMD3->Write("gr_v33_NchFMD12FMD3");

  // 保存 gr_v2_pt 为图片
  TCanvas *c_v2_pt = new TCanvas("c_v2_pt", "v2 vs pT", 800, 600);
  gr_v2_pt->Draw("ape");
  c_v2_pt->SaveAs("v2_pt.png");
  
  // 保存 gr_v2_Nch 为图片
  TCanvas *c_v2_Nch = new TCanvas("c_v2_Nch", "v2 vs Nch", 800, 600);
  gr_v2_Nch->Draw("ape");
  c_v2_Nch->SaveAs("v2_Nch.png");
  
  // 保存 gr_v3_pt 为图片
  TCanvas *c_v3_pt = new TCanvas("c_v3_pt", "v3 vs pT", 800, 600);
  gr_v3_pt->Draw("ape");
  c_v3_pt->SaveAs("v3_pt.png");
  
  // 保存 gr_v3_Nch 为图片
  TCanvas *c_v3_Nch = new TCanvas("c_v3_Nch", "v3 vs Nch", 800, 600);
  gr_v3_Nch->Draw("ape");
  c_v3_Nch->SaveAs("v3_Nch.png");
  
  // 保存 gr_v22_NchTPCFMD12 为图片
  TCanvas *c_v22_NchTPCFMD12 = new TCanvas("c_v22_NchTPCFMD12", "v22 vs Nch TPCFMD12", 800, 600);
  gr_v22_NchTPCFMD12->Draw("ape");
  c_v22_NchTPCFMD12->SaveAs("v22_NchTPCFMD12.png");
  
  // 保存 gr_v33_NchTPCFMD12 为图片
  TCanvas *c_v33_NchTPCFMD12 = new TCanvas("c_v33_NchTPCFMD12", "v33 vs Nch TPCFMD12", 800, 600);
  gr_v33_NchTPCFMD12->Draw("ape");
  c_v33_NchTPCFMD12->SaveAs("v33_NchTPCFMD12.png");
  
  // 保存 gr_v22_NchTPCFMD3 为图片
  TCanvas *c_v22_NchTPCFMD3 = new TCanvas("c_v22_NchTPCFMD3", "v22 vs Nch TPCFMD3", 800, 600);
  gr_v22_NchTPCFMD3->Draw("ape");
  c_v22_NchTPCFMD3->SaveAs("v22_NchTPCFMD3.png");
  
  // 保存 gr_v33_NchTPCFMD3 为图片
  TCanvas *c_v33_NchTPCFMD3 = new TCanvas("c_v33_NchTPCFMD3", "v33 vs Nch TPCFMD3", 800, 600);
  gr_v33_NchTPCFMD3->Draw("ape");
  c_v33_NchTPCFMD3->SaveAs("v33_NchTPCFMD3.png");
  
  // 保存 gr_v22_NchFMD12FMD3 为图片
  TCanvas *c_v22_NchFMD12FMD3 = new TCanvas("c_v22_NchFMD12FMD3", "v22 vs Nch FMD12FMD3", 800, 600);
  gr_v22_NchFMD12FMD3->Draw("ape");
  c_v22_NchFMD12FMD3->SaveAs("v22_NchFMD12FMD3.png");
  
  // 保存 gr_v33_NchFMD12FMD3 为图片
  TCanvas *c_v33_NchFMD12FMD3 = new TCanvas("c_v33_NchFMD12FMD3", "v33 vs Nch FMD12FMD3", 800, 600);
  gr_v33_NchFMD12FMD3->Draw("ape");
  c_v33_NchFMD12FMD3->SaveAs("v33_NchFMD12FMD3.png");
}
