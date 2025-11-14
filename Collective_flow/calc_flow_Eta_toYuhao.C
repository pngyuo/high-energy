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
//	TH1D *htmp1 = (TH1D*)hin->ProjectionX("htmp1", 1, hin->GetYaxis()->FindBin(-2));
//	TH1D *htmp2 = (TH1D*)hin->ProjectionX("htmp2", hin->GetYaxis()->FindBin(2), hin->GetNbinsY());

//	TH1D *htmp1 = (TH1D*)hin->ProjectionX("htmp1", hin->GetYaxis()->FindBin(-12.0), hin->GetYaxis()->FindBin(-2.5));
//	TH1D *htmp2 = (TH1D*)hin->ProjectionX("htmp2", hin->GetYaxis()->FindBin(2.5), hin->GetYaxis()->FindBin(12.0));


  //hout->Add(htmp1, htmp2);
  hout = (TH1D*)hin->ProjectionX();
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

//follow method in Eq.5-7 arXiv:1609.06213
TH1D *getYFromHist(TH2D *hSame, TH2D *hMix, TString name, double Ntrig){
  /*
  TH1D *hSameRidge = sumRidge(hSame, "same"+name);
  TH1D *hMixRidge = sumRidge(hMix, "mix"+name);
  TH1D *hY = getY(hSameRidge, hMixRidge, "Y"+name, Ntrig);
  */

  hMix->Scale(1./hMix->Integral());
  hSame->Divide(hMix);
  TCanvas c0;
  hSame->Draw("surf1");
  c0.SaveAs("test_figure/2D_"+TString(hSame->GetName())+".png");
  TH1D *hout = new TH1D(name, "", hSame->GetNbinsX(), hSame->GetXaxis()->GetXmin(), hSame->GetXaxis()->GetXmax()); 
  hout = (TH1D*)hSame->ProjectionX("_pjx");
  hout->SetName("Y"+name);
  hout->Scale(1./Ntrig/2/3.14159);
  return hout;
}



TF1 *gFunRidge = new TF1("gFunRidge", "[0]*(1+2*[1]*cos(2*x)+2*[2]*cos(3*x)+2*[3]*cos(4*x))", -PI*0.5, PI*1.5);

double funTemplate(double *x, double *par){
  double dPhi = x[0];
  double F = par[0];
  double v2 = par[1];
  double v3 = par[2];
  double v4 = par[3];
  //double G = par[4];
  //TF1 *gFunRidge = new TF1("gFunRidge", "[0]*(1+2*[1]*cos(2*x))", -PI*0.5, PI*1.5);
  double YIntegPeriph = hGlobalPeriph->Integral("width");
  //double YIntegPeriph = 0;
  double G = (gYIntegral-F*YIntegPeriph)/2./PI; //dDeltaPhi P.S. factor  considered
  //double G = F; 
  gFunRidge->SetParameters(G, v2, v3, v4);
  //gFunRidge->SetParameters(G, v2);
  return gFunRidge->Eval(dPhi)+F*hGlobalPeriph->Interpolate(dPhi);
  //return gFunRidge->Eval(dPhi);
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
  TF1 *fitFun_nosub = new TF1("fitFun_nosub", "[0]*(1+2*[1]*cos(2*x)+2*[2]*cos(3*x)+2*[3]*cos(x))", -PI*0.5, PI*1.5);
  fitFun_nosub->SetParameters(1,1,1,1);
  hYHigh->Fit("fitFun_nosub");
  err = fitFun_nosub->GetParError(1);
  v33 = fitFun_nosub->GetParameter(2);
  v33err = fitFun_nosub->GetParError(2);

  SavePlot_nosub(hYHigh, hYHigh->GetName());
  return fitFun_nosub->GetParameter(1);
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

    for(int cent = 0; cent < centrality_bins; cent++) {
        // Initialize arrays for each centrality bin
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
            double v22 = getV22(hYHigh_etabin, hYLow_etabin, v22_err, v33, v33_err);

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

void calc_flow_Eta(){

	TFile *file_1 = new TFile("hist_ampt_normal_1.5mb_a_0.8_b_0.4_Decorr_yuhao.root");
	TFile *file_2 = new TFile("hist_ampt_normal_1.5mb_a_0.8_b_0.4_Decorr_yuhao.root");



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


  //perform Centrality dependent analysis in one high Multiplicity bin e.g. Nch>60
  const int N_etabin=10;
  double eta_min=-2.5;
  double eta_max=2.5;
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

    TH1D *hYHigh_etabin = getYFromHist(hsameHigh_etabin[i], hmixedHigh_etabin[i], Form("High_etabin%d",i), hTrigPtHigh_etabin[i]->Integral());

    hsameLow_etabin[i]=getProjectEta(hsameLowEta,etabin_low, etabin_high,Form("_etabin%d",i));
    hmixedLow_etabin[i]=getProjectEta(hmixLowEta,etabin_low, etabin_high,Form("_etabin%d",i));
    hTrigPtLow_etabin[i]=hTrigPtLowEta->ProjectionX(Form("hTrigPtLow_etabin%d",i), etabin_low, etabin_high);

    TH1D *hYLow_etabin = getYFromHist(hsameLow_etabin[i], hmixedLow_etabin[i], Form("Low_etabin%d",i), hTrigPtLow_etabin[i]->Integral());


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
  

  gr_v22_Ref_Nch->Write("gr_v22_Ref_Nch");
  gr_v33_Ref_Nch->Write("gr_v33_Ref_Nch");
  gr_v22_Pid_Nch->Write("gr_v22_Pid_Nch");
  gr_v33_Pid_Nch->Write("gr_v33_Pid_Nch");

}
