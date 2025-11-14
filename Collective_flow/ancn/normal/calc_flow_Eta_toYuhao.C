const double PI=3.1415926;
TH1D *hGlobalPeriph;
double gYIntegral = 0;

struct CentralityResults {
    TGraphErrors* gr_v2_eta[9];
    TGraphErrors* gr_v22_High_eta[9];
    TGraphErrors* gr_v3_eta[9];
    TGraphErrors* gr_v33_High_eta[9];
    double fitPar1[9];  // Stores par[1] from fitEta for each centrality
    double fitPar2[9];  // Stores par[2] from fitEta for each centrality
    double fitParErr1[9];  // Stores par[1] from fitEta for each centrality
    double fitParErr2[9];  // Stores par[2] from fitEta for each centrality
};

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
    TH1D *htmp1 = (TH1D*)hin->ProjectionX("htmp1", hin->GetYaxis()->FindBin(-7.0), hin->GetYaxis()->FindBin(-2.0));
    TH1D *htmp2 = (TH1D*)hin->ProjectionX("htmp2", hin->GetYaxis()->FindBin(2.0), hin->GetYaxis()->FindBin(7.0)); 
    hout->Add(htmp1, htmp2);
    delete htmp1;
    delete htmp2;
    return hout;
}

TH1D *getY(TH1D *hSame, TH1D *hMix, TString name, double Ntrig) {
    TH1D *hout = new TH1D(name, "", hSame->GetNbinsX(), hSame->GetXaxis()->GetXmin(), hSame->GetXaxis()->GetXmax());
    hout->Divide(hSame, hMix);
    double Btot = hMix->Integral();
    hout->Scale(1. / Ntrig / 2 / PI * Btot);
    return hout;
}

TH1D *getYFromHist(TH2D *hSame, TH2D *hMix, TString name, double Ntrig) {
    TH1D *hSameRidge = sumRidge(hSame, "same" + name);
    TH1D *hMixRidge = sumRidge(hMix, "mix" + name);
    TH1D *hY = getY(hSameRidge, hMixRidge, "Y" + name, Ntrig);
    delete hSameRidge;
    delete hMixRidge;
    return hY;
}



TF1 *gFunRidge = new TF1("gFunRidge", "[0]*(1+2*[1]*cos(2*x)+2*[2]*cos(3*x)+2*[3]*cos(4*x))", -PI*0.5, PI*1.5);

double funTemplate(double *x, double *par) {
    double dPhi = x[0];
    double G = par[0];
    double a1 = par[1];
    double a2 = par[2];
    double a3 = par[3];
    double a4 = par[4];
    TF1 *funRidge = new TF1("funRidge", "[0]*(1 + 2*[1]*cos(x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x) + 2*[4]*cos(4*x))", -0.5*PI, 1.5*PI);
    funRidge->SetParameters(G, a1, a2, a3, a4);
    double result = funRidge->Eval(dPhi);
    delete funRidge;
    return result;
}

void SavePlot(TH1D *hYHigh, TH1D *hYLow, TString name){

  TCanvas c0;
  double ymin=hYHigh->GetMinimum();
  double ymax=hYHigh->GetMaximum();
  double ydiff=abs(ymax-ymin);
  if(ydiff>1e-3) hYHigh->GetYaxis()->SetRangeUser(hYHigh->GetMinimum()-0.5*ydiff,hYHigh->GetMaximum()+1.0*ydiff);
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
  text->DrawLatexNDC(0.25, 0.6, Form("G=%f",G));
  funV2->SetParameter(0,G);
  funV2->SetParameter(1,par[1]);
  funV2->SetParameter(2,funPeriph->Eval(0));
  funV2->SetLineColor(kBlue);
  funV2->Draw("same");

  c0.SaveAs("test_figure/dphi_template_"+name+".png");

}


void SavePlot_nosub(TH1D *hYHigh, TString name){

  TCanvas c0;
  double ymin=hYHigh->GetMinimum();
  double ymax=hYHigh->GetMaximum();
  double ydiff=abs(ymax-ymin);
  if(ydiff>1e-3) hYHigh->GetYaxis()->SetRangeUser(hYHigh->GetMinimum()-0.5*ydiff,hYHigh->GetMaximum()+1.0*ydiff);
  hYHigh->Draw();
  double par[4];
  for(int i=0; i<4; i++){
    par[i]=hYHigh->GetFunction("fitFun_nosub")->GetParameter(i);
  }
  hYHigh->GetFunction("fitFun_nosub")->GetParameters(par);
  double chi2 = hYHigh->GetFunction("fitFun_nosub")->GetChisquare();
  int ndf = hYHigh->GetFunction("fitFun_nosub")->GetNDF();
  TLatex *text = new TLatex();
  text->SetTextSize(0.04);
  text->DrawLatexNDC(0.25, 0.85, Form("F=%f",par[0]));
  text->DrawLatexNDC(0.25, 0.80, Form("v2,2=%f",par[1]));
  text->DrawLatexNDC(0.25, 0.75, Form("v3,3=%f",par[2]));
  text->DrawLatexNDC(0.25, 0.70, Form("v1,1=%f",par[3]));
  text->DrawLatexNDC(0.25, 0.65, Form("chi2=%.2f,ndf=%d",chi2,ndf));

  TF1 *funV2 = new TF1("funV2", "[0]*(1+2*[1]*cos(2*x))", -PI*0.5, PI*1.5);
  double F = par[0];
  funV2->SetParameter(0,F);
  funV2->SetParameter(1,par[1]);
  funV2->SetLineColor(kBlue);
  funV2->Draw("same");

  TF1 *funV1 = new TF1("funV1", "[0]*(1+2*[1]*cos(x))", -PI*0.5, PI*1.5);
  funV1->SetParameter(0,F);
  funV1->SetParameter(1,par[3]);
  funV1->SetLineColor(kRed);
  funV1->Draw("same");


  c0.SaveAs("test_figure/dphi_nosub_"+name+".png");

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
double getV22_nosub(TH1D *hYHighIn, double &err, double &v33, double &v33err){
  TString hName=TString(hYHighIn->GetName())+"_nosub";
  TH1D *hYHigh = (TH1D*)hYHighIn->Clone(hName.Data());
  TF1 *fitFun_nosub = new TF1("fitFun_nosub", "[0]*(1 + 2*[1]*cos(x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x) + 2*[4]*cos(4*x))", -PI*0.5, PI*1.5);
  fitFun_nosub->SetParameters(100, -0.1, 0.01, 0.1, 0.1);
  hYHigh->Fit("fitFun_nosub");
  err = fitFun_nosub->GetParError(2);
  v33 = fitFun_nosub->GetParameter(3);
  v33err = fitFun_nosub->GetParError(3);

  SavePlot_nosub(hYHigh, hYHigh->GetName());
  return fitFun_nosub->GetParameter(2);
}

double getAn(TH1D *hY, double &a1_out, double &a2_out, double &a3_out, double &a4_out, 
             double &G_out, double &chi2_out, int &ndf_out, double &a2err_out, double &a3err_out, // 添加a3err_out
             TString name, TString multiplicity, TString fileLabel) {
    gYIntegral = hY->Integral("width");
    TF1 *fitFun = new TF1("fitFun", funTemplate, -0.5 * PI, 1.5 * PI, 5);
    fitFun->SetParameters(100, -0.1, 0.01, 0.1, 0.1);
    hY->Fit("fitFun", "Q");
    
    a1_out = fitFun->GetParameter(1);
    a2_out = fitFun->GetParameter(2);
    a3_out = fitFun->GetParameter(3);
    a4_out = fitFun->GetParameter(4);
    a2err_out = fitFun->GetParError(2);
    a3err_out = fitFun->GetParError(3); // 获取v3的误差
    G_out = gYIntegral / (2.0 * PI);
    chi2_out = fitFun->GetChisquare();
    ndf_out = fitFun->GetNDF();
    
    // SavePlot(hY, name, multiplicity, fileLabel); // 确保SavePlot函数存在且参数匹配
    double result = fitFun->GetParameter(2);
    delete fitFun;
    return result;
}


/*
double getV22_nosub(TH1D *hYHighIn, TH1D *hYLow, double &err, double &v33, double &v33err){
  double V2=0;
  double V3=0;
  double count=0;
  for(int i=1; i<=hYHighIn->GetNbinsX(); i++){
    double binCont=hYHighIn->GetBinContent(i);
    double phi=hYHighIn->GetBinCenter(i);
    V2+=cos(2*phi)*binCont;
    V3+=cos(3*phi)*binCont;
    count+=binCont;
  }
  if(count>0){
    V2/=count;
    V3/=count;
  }

  v33=V3;

  return V2;
}
*/


TH2D* getProjectEta(TH3D* hinput, int etabin_low, int etabin_high, TString hName){
  hinput->GetZaxis()->SetRange(etabin_low, etabin_high);
  TString option = hName+"_yxoe";
  //yx means y vs x, o keeps original axis but only bins in SetRange, e recalculates the error
  //h3d->Project3D("name_yxoe") will generate a hist with name: h3dname_name_yx
  TH2D *hout = (TH2D*)hinput->Project3D(option);
  return hout;
}

CentralityResults calcFlowEtaCentrality(TFile* file_1, TFile* file_2, int centrality_bins = 9) {
    CentralityResults results;
    const int N_etabin=10;
    double eta_min=-2.5;
    double eta_max=2.5;
    double eta_bin=(eta_max-eta_min)/N_etabin;
    TH1D *hNChMid = (TH1D *)file_1->Get("hNChMid");
    const int N_cent = 8;
    double Nsel_cut[N_cent + 1] = {10, 30, 40, 50, 60, 80, 100, 120, 150};
    double Nsel_val[N_cent] = {10};
    for(int cent = 0; cent < centrality_bins; cent++) {
        // Initialize arrays for each centrality bin
                Nsel_val[cent] = (Nsel_cut[cent] + Nsel_cut[cent + 1]) / 2.0;
                double centMin = Nsel_cut[cent];
        double centMax = Nsel_cut[cent + 1];
        double etabin_val[N_etabin]={0};
        double etabin_edge[N_etabin+1]={0};
        for(int i=0; i<N_etabin+1; i++){
            etabin_edge[i]=eta_min+i*eta_bin;
        }
        double etabin_width[N_etabin]={0};
        double v2_eta[N_etabin]={0};
        double v2_eta_err[N_etabin]={0};
        double v3_eta[N_etabin]={0};
        double v3_eta_err[N_etabin]={0};
        double v22_High_eta[N_etabin]={0};
        double v22_High_eta_err[N_etabin]={0};
        double v33_High_eta[N_etabin]={0};
        double v33_High_eta_err[N_etabin]={0};

        // Get histograms for current centrality bin
        TH3D* hsameHighEta = (TH3D*)file_2->Get(Form("hDEtaDPhiTrigEtaSameEvent_Cent%d", cent));
        TH3D* hmixHighEta = (TH3D*)file_2->Get(Form("hDEtaDPhiTrigEtaMixEvent_Cent%d", cent));
        TH2D* hTrigPtHighEta = (TH2D*)file_2->Get(Form("hTrigPtEta_Cent%d", cent));

        TH3D* hsameLowEta = (TH3D*)file_1->Get("hDEtaDPhiTrigEtaSameEventLowMid");
        TH3D* hmixLowEta = (TH3D*)file_1->Get("hDEtaDPhiTrigEtaMixEventLowMid");
        TH2D* hTrigPtLowEta = (TH2D*)file_1->Get("hTrigPtEtaLow");
        double totalChargedParticleslow_cent = 0;
        double nEventslow_cent = 0;
        for (int bin = 1; bin <= hNChMid->GetNbinsX(); bin++) {
            double binCenter = hNChMid->GetXaxis()->GetBinCenter(bin);
            double binContent = hNChMid->GetBinContent(bin);
            if (binCenter >= 10 && binCenter < 30) {
                totalChargedParticleslow_cent += binCenter * binContent;
                nEventslow_cent += binContent;
            }
        }
        double averageNTlow_cent = totalChargedParticleslow_cent / nEventslow_cent;
        
        // Calculate averageNThigh_cent
        double totalChargedParticleshigh_cent = 0;
        double nEventshigh_cent = 0;
        for (int bin = 1; bin <= hNChMid->GetNbinsX(); bin++) {
            double binCenter = hNChMid->GetXaxis()->GetBinCenter(bin);
            double binContent = hNChMid->GetBinContent(bin);
            if (binCenter >= centMin && binCenter < centMax) {
                totalChargedParticleshigh_cent += binCenter * binContent;
                nEventshigh_cent += binContent;
            }
        }
        double averageNThigh_cent = totalChargedParticleshigh_cent / nEventshigh_cent;
        double k = averageNTlow_cent / averageNThigh_cent;
        // Process each eta bin
        for(int i=0; i<N_etabin; i++){
            etabin_val[i]=(etabin_edge[i]+etabin_edge[i+1])*0.5;
            int etabin_low=hTrigPtHighEta->GetYaxis()->FindBin(etabin_edge[i]);
            int etabin_high=hTrigPtHighEta->GetYaxis()->FindBin(etabin_edge[i+1])-1;

            hTrigPtHighEta->GetYaxis()->SetRange(etabin_low, etabin_high);
            hTrigPtLowEta->GetYaxis()->SetRange(etabin_low, etabin_high);

            TH2D* hsameHigh_etabin = getProjectEta(hsameHighEta,etabin_low, etabin_high,Form("_cent%d_etabin%d",cent,i));
            TH2D* hmixedHigh_etabin = getProjectEta(hmixHighEta,etabin_low, etabin_high,Form("_cent%d_etabin%d",cent,i));
            TH1D* hTrigPtHigh_etabin = hTrigPtHighEta->ProjectionX(Form("hTrigPtHigh_cent%d_etabin%d",cent,i), etabin_low, etabin_high);

            TH1D *hYHigh_etabin = getYFromHist(hsameHigh_etabin, hmixedHigh_etabin, Form("High_cent%d_etabin%d",cent,i), hTrigPtHigh_etabin->Integral());

            TH2D* hsameLow_etabin = getProjectEta(hsameLowEta,etabin_low, etabin_high,Form("_cent%d_etabin%d",cent,i));
            TH2D* hmixedLow_etabin = getProjectEta(hmixLowEta,etabin_low, etabin_high,Form("_cent%d_etabin%d",cent,i));
            TH1D* hTrigPtLow_etabin = hTrigPtLowEta->ProjectionX(Form("hTrigPtLow_cent%d_etabin%d",cent,i), etabin_low, etabin_high);

            TH1D *hYLow_etabin = getYFromHist(hsameLow_etabin, hmixedLow_etabin, Form("Low_cent%d_etabin%d",cent,i), hTrigPtLow_etabin->Integral());

            double v33 = 0;
            double v33_err = 0;
            double v22_err = 0;
            double v22high = getV22_nosub(hYHigh_etabin,  v22_err, v33, v33_err);
            double v22low = getV22_nosub(hYLow_etabin,  v22_err, v33, v33_err);
            double v22 = v22high - v22low * k;
       
            v2_eta[i] = sqrt(v22);
            v2_eta_err[i] = v2_eta[i]*sqrt( pow(v22_err/v22,2.)+ pow(0.5*v22_err/v22,2.) );

            v3_eta[i] = sqrt(v33);
            v3_eta_err[i] = v3_eta[i]*sqrt( pow(v33_err/v33,2.)+ pow(0.5*v33_err/v33,2.) );

            v22_High_eta[i]=v22;
            v22_High_eta_err[i]=v22_err;
            v33_High_eta[i]=v33;
            v33_High_eta_err[i]=v33_err;
        }

        // Create graphs for current centrality bin
        results.gr_v2_eta[cent] = new TGraphErrors(N_etabin, etabin_val, v2_eta, etabin_width, v2_eta_err);
        results.gr_v22_High_eta[cent] = new TGraphErrors(N_etabin, etabin_val, v22_High_eta, etabin_width, v22_High_eta_err);
        results.gr_v3_eta[cent] = new TGraphErrors(N_etabin, etabin_val, v3_eta, etabin_width, v3_eta_err);
        results.gr_v33_High_eta[cent] = new TGraphErrors(N_etabin, etabin_val, v33_High_eta, etabin_width, v22_High_eta_err);

        // Perform fits and store parameters
        TF1 *fitEta = new TF1("fitEta", "[0]*(1+[1]*x+[2]*x*x)", -3, 3);
        results.gr_v22_High_eta[cent]->Fit("fitEta");
        results.fitPar1[cent] = fitEta->GetParameter(1);
        results.fitParErr1[cent] = fitEta->GetParError(1);     
        results.gr_v33_High_eta[cent]->Fit("fitEta");         
        results.fitPar2[cent] = fitEta->GetParameter(1);
        results.fitParErr2[cent] = fitEta->GetParError(1);
    }

    return results;
}

void calc_flow_Eta_toYuhao(){

	TFile *file_1 = new TFile("hist_ampt_normal_0.15mb_a_0.8_b_0.4_Decorr_yuhao.root");
	TFile *file_2 = new TFile("hist_ampt_normal_0.15mb_a_0.8_b_0.4_Decorr_yuhao.root");



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

//	TFile *infile_Pid_Low = new TFile("hist_root/hist_ampt_3q_1.5mb_Decorr.root");
//	TFile *infile_Pid_High = new TFile("hist_root/hist_ampt_3q_1.5mb_Decorr.root");
	TFile *infile_Pid_Low = new TFile("hist_ampt_normal_1.5mb_a_0.8_b_0.4_Decorr_yuhao.root");
	TFile *infile_Pid_High = new TFile("hist_ampt_normal_1.5mb_a_0.8_b_0.4_Decorr_yuhao.root");



  CentralityResults results = calcFlowEtaCentrality(file_1, file_2);
  double Nchsel[8]={35, 45, 55, 70, 90, 110, 135, 175};
  double Nchsel_err[8]={0};
  double F22_Nch[8]={0};
  double F33_Nch[8]={0};
  double F22Err_Nch[8]={0};
  double F33Err_Nch[8]={0}; 
  for(int igr=0; igr<8; igr++){
    F22_Nch[igr]=results.fitPar1[igr+1];
    F33_Nch[igr]=results.fitPar2[igr+1];
    F22Err_Nch[igr]=results.fitParErr1[igr+1];
    F33Err_Nch[igr]=results.fitParErr2[igr+1];
  }
	TGraphErrors *gr_F22_Nch = new TGraphErrors(8, Nchsel, F22_Nch, Nchsel_err, F22Err_Nch);
	TGraphErrors *gr_F33_Nch = new TGraphErrors(8, Nchsel, F33_Nch, Nchsel_err, F33Err_Nch);


  const int N_etabin=30;
  double eta_min=-3;
  double eta_max=3;
  double eta_bin=(eta_max-eta_min)/N_etabin;
  double etabin_val[N_etabin]={0};
  double etabin_edge[N_etabin+1]={0};
  for(int i=0; i<N_etabin+1; i++){
    etabin_edge[i]=eta_min+i*eta_bin;
  }
  double etabin_width[N_etabin]={0};
  double v2_eta[N_etabin]={0};
  double v2_eta_err[N_etabin]={0};
  double v3_eta[N_etabin]={0};
  double v3_eta_err[N_etabin]={0};
  double v22_High_eta[N_etabin]={0};
  double v22_High_eta_err[N_etabin]={0};
  double v33_High_eta[N_etabin]={0};
  double v33_High_eta_err[N_etabin]={0};

  double v22_nosub_High_eta[N_etabin]={0};
  double v22_nosub_High_eta_err[N_etabin]={0};
  double v33_nosub_High_eta[N_etabin]={0};
  double v33_nosub_High_eta_err[N_etabin]={0};

  TH2D *hsameHigh_etabin[N_etabin];
  TH2D *hmixedHigh_etabin[N_etabin];
  TH1D *hTrigPtHigh_etabin[N_etabin];
  TH2D *hsameLow_etabin[N_etabin];
  TH2D *hmixedLow_etabin[N_etabin];
  TH1D *hTrigPtLow_etabin[N_etabin];

//  TH3D* hsameHighEta = (TH3D*)file_2->Get("hDEtaDPhiTrigEtaSameEvent_Cent4");
//  TH3D* hmixHighEta = (TH3D*)file_2->Get("hDEtaDPhiTrigEtaMixEvent_Cent4");
//  TH2D* hTrigPtHighEta = (TH2D*)file_2->Get("hTrigPtEta_Cent4");
  TH3D* hsameHighEta = (TH3D*)file_2->Get("hDEtaDPhiTrigEtaSameEventHighMid");
  TH3D* hmixHighEta = (TH3D*)file_2->Get("hDEtaDPhiTrigEtaMixEventHighMid");
  TH2D* hTrigPtHighEta = (TH2D*)file_2->Get("hTrigPtEtaHigh");

  TH3D* hsameLowEta = (TH3D*)file_1->Get("hDEtaDPhiTrigEtaSameEventLowMid");
  TH3D* hmixLowEta = (TH3D*)file_1->Get("hDEtaDPhiTrigEtaMixEventLowMid");
  TH2D* hTrigPtLowEta = (TH2D*)file_1->Get("hTrigPtEtaLow");


  for(int i=0; i<N_etabin; i++){
    etabin_val[i]=(etabin_edge[i]+etabin_edge[i+1])*0.5;
    int etabin_low=hTrigPtHighEta->GetYaxis()->FindBin(etabin_edge[i]);
    int etabin_high=hTrigPtHighEta->GetYaxis()->FindBin(etabin_edge[i+1])-1;
    cout<<etabin_val[i]<<" "<<etabin_low<<" "<<etabin_high<<" "<<hTrigPtHighEta->GetNbinsY()<<endl;
    /*
    hsameHighEta->GetZaxis()->SetRange(etabin_low, etabin_high);
    hmixHighEta->GetZaxis()->SetRange(etabin_low, etabin_high);
    hsameLowEta->GetZaxis()->SetRange(etabin_low, etabin_high);
    hmixLowEta->GetZaxis()->SetRange(etabin_low, etabin_high);
    */
    hTrigPtHighEta->GetYaxis()->SetRange(etabin_low, etabin_high);
    hTrigPtLowEta->GetYaxis()->SetRange(etabin_low, etabin_high);

    hsameHigh_etabin[i]=getProjectEta(hsameHighEta,etabin_low, etabin_high,Form("_etabin%d",i));
    hmixedHigh_etabin[i]=getProjectEta(hmixHighEta,etabin_low, etabin_high,Form("_etabin%d",i));
    hTrigPtHigh_etabin[i]=hTrigPtHighEta->ProjectionX(Form("hTrigPtHigh_etabin%d",i), etabin_low, etabin_high);
    //etabin_val[i]=hTrigPtHigh_etabin[i]->GetMean();

    TH1D *hYHigh_etabin = getYFromHist(hsameHigh_etabin[i], hmixedHigh_etabin[i], Form("High_etabin%d",i),10);

    hsameLow_etabin[i]=getProjectEta(hsameLowEta,etabin_low, etabin_high,Form("_etabin%d",i));
    hmixedLow_etabin[i]=getProjectEta(hmixLowEta,etabin_low, etabin_high,Form("_etabin%d",i));
    hTrigPtLow_etabin[i]=hTrigPtLowEta->ProjectionX(Form("hTrigPtLow_etabin%d",i), etabin_low, etabin_high);

    TH1D *hYLow_etabin = getYFromHist(hsameLow_etabin[i], hmixedLow_etabin[i], Form("Low_etabin%d",i), 10);


    double v33 = 0;
    double v33_err = 0;
    double v22_err = 0;
    double v22 = getV22(hYHigh_etabin, hYLow_etabin, v22_err, v33, v33_err);

    double v22_nosub_err = 0;
    double v33_nosub = 0;
    double v33_nosub_err = 0;
    double v22_nosub = getV22_nosub(hYHigh_etabin, v22_nosub_err, v33_nosub, v33_nosub_err);

    double v22_low_nosub_err = 0;
    double v33_low_nosub = 0;
    double v33_low_nosub_err = 0;
    double v22_low_nosub = getV22_nosub(hYLow_etabin, v22_low_nosub_err, v33_low_nosub, v33_low_nosub_err);

    /*
    //v2_eta[i] = v22/sqrt(v22Ref);
    v2_eta[i] = sqrt(v22);
    v2_eta_err[i] = v2_eta[i]*sqrt( pow(v22_err/v22,2.)+ pow(0.5*v22_err/v22,2.) );

    //v3_eta[i] = v33/sqrt(v33Ref);
    v3_eta[i] = sqrt(v33);
    //v3_eta_err[i] = v3_eta[i]*sqrt( pow(v33_err/v33,2.)+ pow(0.5*v33Ref_err/v33Ref,2.) );
    v3_eta_err[i] = v3_eta[i]*sqrt( pow(v33_err/v33,2.)+ pow(0.5*v33_err/v33,2.) );
    */

    //this k scaling needs to be implemented!
    //double k=hTrigPtLow_etabin[i]->Integral()/hTrigPtHigh_etabin[i]->Integral();
    double k=2./7.;
    cout<<"k="<<k<<endl;
    
    if(k<1e-6) {
      cout<<"strange k factor"<<endl;
      exit(0);
    }
    

    
    v22_High_eta[i]=v22;
    v22_High_eta_err[i]=v22_err;
    v33_High_eta[i]=v33;
    v33_High_eta_err[i]=v33_err;
    

    v2_eta[i]=v22_low_nosub;
    v2_eta_err[i]=v22_low_nosub_err;

    /*
    v22_High_eta[i]=v22_nosub-k*v22_low_nosub;
    v22_High_eta_err[i]=sqrt(v22_nosub_err*v22_nosub_err+k*k*v22_low_nosub_err*v22_low_nosub_err);
    v33_High_eta[i]=v33_nosub-k*v33_low_nosub;
    v33_High_eta_err[i]=sqrt(v33_nosub_err*v33_nosub_err+k*k*v33_low_nosub_err*v33_low_nosub_err);
    */

    v22_nosub_High_eta[i]=v22_nosub;
    v22_nosub_High_eta_err[i]=v22_nosub_err;
    v33_nosub_High_eta[i]=v33_nosub;
    v33_nosub_High_eta_err[i]=v33_nosub_err;
  }

  TGraphErrors *gr_v2_eta = new TGraphErrors(N_etabin, etabin_val, v2_eta, etabin_width, v2_eta_err);
  TGraphErrors *gr_v22_High_eta = new TGraphErrors(N_etabin, etabin_val, v22_High_eta, etabin_width, v22_High_eta_err);
  TGraphErrors *gr_v3_eta = new TGraphErrors(N_etabin, etabin_val, v3_eta, etabin_width, v3_eta_err);
  TGraphErrors *gr_v33_High_eta = new TGraphErrors(N_etabin, etabin_val, v33_High_eta, etabin_width, v33_High_eta_err);

  TGraphErrors *gr_v22_nosub_High_eta = new TGraphErrors(N_etabin, etabin_val, v22_nosub_High_eta, etabin_width, v22_nosub_High_eta_err);
  TGraphErrors *gr_v33_nosub_High_eta = new TGraphErrors(N_etabin, etabin_val, v33_nosub_High_eta, etabin_width, v33_nosub_High_eta_err);

  TF1 *fitEta = new TF1("fitEta", "[0]*(1+[1]*x+[2]*x*x)", -3, 3);
  gr_v22_High_eta->Fit("fitEta");
  cout<<"F2="<<fitEta->GetParameter(1)<<endl;
  gr_v33_High_eta->Fit("fitEta");
  cout<<"F3="<<fitEta->GetParameter(1)<<endl;



  TFile *outfile = new TFile("Decorr_ampt_normal_a_0.8_b_0.4_1.5mb.root", "recreate");

  gr_v33_nosub_High_eta->Write("gr_v33_nosub_High_eta");
  gr_v2_eta->Write("gr_v22_low_eta");
  gr_v22_High_eta->Write("gr_v22_High_eta");
  gr_v33_High_eta->Write("gr_v33_High_eta");
  gr_v22_nosub_High_eta->Write("gr_v22_nosub_High_eta");
  gr_v33_nosub_High_eta->Write("gr_v33_nosub_High_eta");
 
  for(int igr=0; igr<9; igr++){
    results.gr_v22_High_eta[igr]->Write(Form("gr_v22_High_eta_Cent%d",igr));
    results.gr_v33_High_eta[igr]->Write(Form("gr_v33_High_eta_Cent%d",igr));
  }
  gr_F22_Nch->Write("gr_F22_Nch");
  gr_F33_Nch->Write("gr_F33_Nch");
  // 在文件写入之前添加画图代码
TCanvas *c1 = new TCanvas("c1", "v2 and v3 vs eta", 1200, 800);
c1->Divide(2,2);

// 绘制 gr_v22_High_eta
c1->cd(1);
gr_v22_High_eta->SetTitle("v_{2}{2} vs #eta; #eta; v_{2}{2}");
gr_v22_High_eta->SetMarkerStyle(20);
gr_v22_High_eta->SetMarkerColor(kBlue);
gr_v22_High_eta->Draw("AP");

// 绘制 gr_v33_High_eta
c1->cd(2);
gr_v33_High_eta->SetTitle("v_{3}{2} vs #eta; #eta; v_{3}{2}");
gr_v33_High_eta->SetMarkerStyle(21);
gr_v33_High_eta->SetMarkerColor(kRed);
gr_v33_High_eta->Draw("AP");

// 绘制 gr_v22_nosub_High_eta
c1->cd(3);
gr_v22_nosub_High_eta->SetTitle("v_{2}{2} (no subtraction) vs #eta; #eta; v_{2}{2}");
gr_v22_nosub_High_eta->SetMarkerStyle(22);
gr_v22_nosub_High_eta->GetXaxis()->SetRangeUser(-3, 3);
gr_v22_nosub_High_eta->SetMarkerColor(kGreen+2);
gr_v22_nosub_High_eta->Draw("AP");

// 绘制 gr_v33_nosub_High_eta
c1->cd(4);
gr_v33_nosub_High_eta->SetTitle("v_{3}{2} (no subtraction) vs #eta; #eta; v_{3}{2}");
gr_v33_nosub_High_eta->SetMarkerStyle(23);
gr_v33_nosub_High_eta->SetMarkerColor(kMagenta+2);
gr_v33_nosub_High_eta->Draw("AP");

c1->SaveAs("flow_vn_vs_eta.png");

// 另外创建一个画布分别显示每个图形，更大更清晰
TCanvas *c2 = new TCanvas("c2", "v22_High_eta", 800, 600);
gr_v22_High_eta->Draw("AP");
c2->SaveAs("v22_High_eta.png");

TCanvas *c3 = new TCanvas("c3", "v33_High_eta", 800, 600);
gr_v33_High_eta->Draw("AP");
c3->SaveAs("v33_High_eta.png");

TCanvas *c4 = new TCanvas("c4", "v22_nosub_High_eta", 800, 600);
gr_v22_nosub_High_eta->Draw("AP");
c4->SaveAs("v22_nosub_High_eta.png");

TCanvas *c5 = new TCanvas("c5", "v33_nosub_High_eta", 800, 600);
gr_v33_nosub_High_eta->Draw("AP");
c5->SaveAs("v33_nosub_High_eta.png");

// 创建对比图：有减背景 vs 无减背景
TCanvas *c6 = new TCanvas("c6", "v2 comparison", 800, 600);
gr_v22_High_eta->SetTitle("v_{2}{2} comparison; #eta; v_{2}{2}");
gr_v22_High_eta->Draw("AP");
gr_v22_nosub_High_eta->SetMarkerStyle(24);
gr_v22_nosub_High_eta->SetMarkerColor(kRed);
gr_v22_nosub_High_eta->Draw("P same");

TLegend *leg1 = new TLegend(0.7, 0.7, 0.9, 0.9);
leg1->AddEntry(gr_v22_High_eta, "With background subtraction", "p");
leg1->AddEntry(gr_v22_nosub_High_eta, "No background subtraction", "p");
leg1->Draw();

c6->SaveAs("v2_comparison.png");

TCanvas *c7 = new TCanvas("c7", "v3 comparison", 800, 600);
gr_v33_High_eta->SetTitle("v_{3}{2} comparison; #eta; v_{3}{2}");
gr_v33_High_eta->Draw("AP");
gr_v33_nosub_High_eta->SetMarkerStyle(24);
gr_v33_nosub_High_eta->SetMarkerColor(kRed);
gr_v33_nosub_High_eta->Draw("P same");

TLegend *leg2 = new TLegend(0.7, 0.7, 0.9, 0.9);
leg2->AddEntry(gr_v33_High_eta, "With background subtraction", "p");
leg2->AddEntry(gr_v33_nosub_High_eta, "No background subtraction", "p");
leg2->Draw();

c7->SaveAs("v3_comparison.png");

// 清理画布
delete c1;
delete c2;
delete c3;
delete c4;
delete c5;
delete c6;
delete c7;
  // 添加F2和Nch关系图
  TCanvas *c8 = new TCanvas("c8", "F2 vs Nch", 800, 600);
  gr_F22_Nch->SetTitle("F_{2} vs N_{ch}; N_{ch}; F_{2}");
  gr_F22_Nch->SetMarkerStyle(20);
  gr_F22_Nch->SetMarkerColor(kBlue);
  gr_F22_Nch->SetLineColor(kBlue);
  gr_F22_Nch->Draw("AP");
  
  // 添加拟合线（可选）
  TF1 *fitF2 = new TF1("fitF2", "[0]+[1]*x", 0, 200);
  gr_F22_Nch->Fit("fitF2", "Q");
  fitF2->SetLineColor(kRed);
  fitF2->Draw("same");
  
  // 添加图例
  TLegend *leg3 = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg3->AddEntry(gr_F22_Nch, "F_{2} data", "p");
  leg3->AddEntry(fitF2, Form("Fit: %.4f + %.4f*Nch", fitF2->GetParameter(0), fitF2->GetParameter(1)), "l");
  leg3->Draw();
  
  c8->SaveAs("F2_vs_Nch.png");
  
  // 绘制F3和Nch关系图
  TCanvas *c9 = new TCanvas("c9", "F3 vs Nch", 800, 600);
  gr_F33_Nch->SetTitle("F_{3} vs N_{ch}; N_{ch}; F_{3}");
  gr_F33_Nch->SetMarkerStyle(21);
  gr_F33_Nch->SetMarkerColor(kRed);
  gr_F33_Nch->SetLineColor(kRed);
  gr_F33_Nch->Draw("AP");
  
  // 添加拟合线（可选）
  TF1 *fitF3 = new TF1("fitF3", "[0]+[1]*x", 0, 200);
  gr_F33_Nch->Fit("fitF3", "Q");
  fitF3->SetLineColor(kBlue);
  fitF3->Draw("same");
  
  // 添加图例
  TLegend *leg4 = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg4->AddEntry(gr_F33_Nch, "F_{3} data", "p");
  leg4->AddEntry(fitF3, Form("Fit: %.4f + %.4f*Nch", fitF3->GetParameter(0), fitF3->GetParameter(1)), "l");
  leg4->Draw();
  
  c9->SaveAs("F3_vs_Nch.png");
  
  // 绘制F2和F3在同一张图上进行比较
  TCanvas *c10 = new TCanvas("c10", "F2 and F3 vs Nch", 800, 600);
  
  // 先画F2
  gr_F22_Nch->SetTitle("F_{2} and F_{3} vs N_{ch}; N_{ch}; F_{2}, F_{3}");
  gr_F22_Nch->Draw("AP");
  
  // 再画F3
  gr_F33_Nch->Draw("P same");
  
  // 添加图例
  TLegend *leg5 = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg5->AddEntry(gr_F22_Nch, "F_{2}", "p");
  leg5->AddEntry(gr_F33_Nch, "F_{3}", "p");
  leg5->Draw();
  
  c10->SaveAs("F2_F3_vs_Nch_comparison.png");
  
  // 将新创建的画布也写入输出文件
  c8->Write("F2_vs_Nch");
  c9->Write("F3_vs_Nch");
  c10->Write("F2_F3_vs_Nch_comparison");
  
  // 清理新创建的对象
  delete c8;
  delete c9;
  delete c10;
  delete leg3;
  delete leg4;
  delete leg5;
  delete fitF2;
  delete fitF3;
//============================  画图代码结束  ============================//
  // gr_v22_Ref_Nch->Write("gr_v22_Ref_Nch");
  // gr_v33_Ref_Nch->Write("gr_v33_Ref_Nch");
  // gr_v22_Pid_Nch->Write("gr_v22_Pid_Nch");
  // gr_v33_Pid_Nch->Write("gr_v33_Pid_Nch");

}
