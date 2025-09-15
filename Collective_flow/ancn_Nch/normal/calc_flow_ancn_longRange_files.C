#include <iostream>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TLegend.h>
#include <TDirectory.h>

const double PI = 3.1415926;
double gYIntegral = 0;

void divideHist(TH2D *hin) {
    for (int i = 1; i <= hin->GetNbinsX(); ++i) {
        for (int j = 1; j <= hin->GetNbinsY(); ++j) {
            double dphi = hin->GetXaxis()->GetBinCenter(i);
            double deta = hin->GetYaxis()->GetBinCenter(j);
            double ratio = (1 - 1e-4) * (4.75 - fabs(deta)) / 4.75 + 1e-4;
            double val = hin->GetBinContent(i, j);
            val = val / ratio;
            if (fabs(deta) > 4.75) {
                val = 0;
            }
            hin->SetBinContent(i, j, val);
        }
    }
}

void removeBkg(TH1D *hin) {
    double bkg = hin->GetMinimum();
    for (int i = 1; i <= hin->GetNbinsX(); ++i) {
        hin->SetBinContent(i, hin->GetBinContent(i) - bkg);
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
  
void SavePlot(TH1D *hY, TString name, TString multiplicity, TString fileLabel) {
    TCanvas c0;
    hY->Draw();
    double par[5];
    for (int i = 0; i < 5; i++) {
        par[i] = hY->GetFunction("fitFun")->GetParameter(i);
    }
    double chi2 = hY->GetFunction("fitFun")->GetChisquare();
    int ndf = hY->GetFunction("fitFun")->GetNDF();
    TLatex *text = new TLatex();
    text->SetTextSize(0.04);
    text->DrawLatexNDC(0.25, 0.85, Form("G=%f", par[0]));
    text->DrawLatexNDC(0.25, 0.80, Form("a1=%f", par[1]));
    text->DrawLatexNDC(0.25, 0.75, Form("a2=%f", par[2]));
    text->DrawLatexNDC(0.25, 0.70, Form("a3=%f", par[3]));
    text->DrawLatexNDC(0.25, 0.65, Form("a4=%f", par[4]));
    text->DrawLatexNDC(0.25, 0.60, Form("chi2=%.2f, ndf=%d", chi2, ndf));
    c0.SaveAs(Form("test_figure/dphi_template_%s_%s_%s.png", multiplicity.Data(), name.Data(), fileLabel.Data()));
    delete text;
}

double getAn(TH1D *hY, double &a1, double &a2, double &a3, double &a4, double &G, double &chi2, int &ndf, double &a2err, TString name, TString multiplicity, TString fileLabel) {
    gYIntegral = hY->Integral("width");
    TF1 *fitFun = new TF1("fitFun", funTemplate, -0.5 * PI, 1.5 * PI, 5);
    fitFun->SetParameters(100, -0.1, 0.01, 0.1, 0.1);
    hY->Fit("fitFun", "Q");
    a1 = fitFun->GetParameter(1);
    a2 = fitFun->GetParameter(2);
    a3 = fitFun->GetParameter(3);
    a4 = fitFun->GetParameter(4);
    a2err = fitFun->GetParError(2);
    G = gYIntegral / (2.0 * PI);
    chi2 = fitFun->GetChisquare();
    ndf = fitFun->GetNDF();
    SavePlot(hY, name, multiplicity, fileLabel);
    double result = fitFun->GetParameter(2);
    delete fitFun;
    return result;
}

// 处理单个文件的函数
void processFile(TString filename,
                 TString label,
                 int color,
                 TGraphErrors **gr_a2_Nch,
                 TGraphErrors **gr_a3_Nch
                 )
{
    TFile *file_1 = new TFile(filename);
    TFile *file_2 = new TFile(filename); 

    //-------------------------------------------------------
    // 3. N_ch 依赖
    //-------------------------------------------------------
    const int N_cent = 11;
    double Nsel_min      = 0;
    double Nsel_bin      = 10;
    double Nsel_val[N_cent]      = {0};
    double Nsel_bin_width[N_cent]= {0};
    double a2_Nch[N_cent]        = {0};
    double a2_Nch_err[N_cent]    = {0};
    double a3_Nch[N_cent]        = {0};
    double a3_Nch_err[N_cent]    = {0};

    TH2D *hsame_cent[N_cent];
    TH2D *hmixed_cent[N_cent];
    TH1D *hTrigPt_cent[N_cent];

    for (int icent = 0; icent < N_cent; ++icent) {
        Nsel_val[icent]      = Nsel_min + icent * Nsel_bin + 0.5 * Nsel_bin;
        Nsel_bin_width[icent]= 0.5 * Nsel_bin;

        if (Nsel_val[icent] < 85) {
            hsame_cent[icent] = (TH2D *)file_1->Get(Form("hDEtaDPhiSameEvent_Centall%d", icent));
            hmixed_cent[icent]= (TH2D *)file_1->Get(Form("hDEtaDPhiMixEvent_Centall%d", icent));
            hTrigPt_cent[icent]= (TH1D *)file_1->Get(Form("hTrigPt_Cent%d", icent));
        } else {
            hsame_cent[icent] = (TH2D *)file_2->Get(Form("hDEtaDPhiSameEvent_Centall%d", icent));
            hmixed_cent[icent]= (TH2D *)file_2->Get(Form("hDEtaDPhiMixEvent_Centall%d", icent));
            hTrigPt_cent[icent]= (TH1D *)file_2->Get(Form("hTrigPt_Cent%d", icent));
        }

        TH1D *hYHigh_centbin = getYFromHist(hsame_cent[icent], hmixed_cent[icent],
                                            Form("High_cent%d", icent),
                                           100);
        hYHigh_centbin->Rebin(2);

        // 这些函数如果不需要可注释掉
        TH1D *hSumRidge    = sumRidge(hsame_cent[icent],
                                      Form("sumRidge_cent%d", icent));
        TH1D *hmixSumRidge = sumRidge(hmixed_cent[icent],
                                      Form("summixRidge_cent%d", icent));

        double a1_cent = 0, a2_cent = 0, a3_cent = 0, a4_cent = 0;
        double G_cent  = 0, chi2_cent = 0;
        int    ndf_cent = 0;
        double a2_cent_err = 0, a3_cent_err = 0;

        a2_cent = getAn(hYHigh_centbin,
                        a1_cent, a2_cent, a3_cent, a4_cent,
                        G_cent, chi2_cent, ndf_cent, a2_cent_err,
                        Form("cent%d", icent), label, Form("cent%d", icent));

        a2_Nch[icent]     = a2_cent;
        a2_Nch_err[icent] = a2_cent_err;
        a3_Nch[icent]     = a3_cent;
        a3_Nch_err[icent] = a3_cent_err;
    }

    //-------------------------------------------------------
    // 4. 保存结果图
    //-------------------------------------------------------
    *gr_a2_Nch = new TGraphErrors(N_cent,
                                  Nsel_val, a2_Nch,
                                  Nsel_bin_width, a2_Nch_err);
    *gr_a3_Nch = new TGraphErrors(N_cent,
                                  Nsel_val, a3_Nch,
                                  Nsel_bin_width, a3_Nch_err);

    (*gr_a2_Nch)->SetMarkerColor(color);
    (*gr_a2_Nch)->SetLineColor(color);
    (*gr_a3_Nch)->SetMarkerColor(color);
    (*gr_a3_Nch)->SetLineColor(color);

    file_1->Close();
    file_2->Close();
}

void calc_flow_ancn_longRange_files(){
  gStyle->SetOptStat(kFALSE);
  // Create output directories
  gSystem->mkdir("test_figure", kTRUE);
  gSystem->mkdir("test_figure/low_templates", kTRUE);
  gSystem->mkdir("sumridge_plots", kTRUE);
  
  // Define the three files and their properties
  TString filenames[4] = {
    "hist_outputallFSI_Nch.root",
    "hist_ampt_normal_0.15mb_longRange_yuhao.root", 
    "hist_ampt_normal_1.5mb_a_0.8_b_0.4_longRange_yuhao.root",
    "hist_ampt_normal_0.15mb_a_0.8_b_0.4_longRange_yuhao.root"
  };
  
  TString labels[4] = {"normal_1.5mb", "normal_0.15mb", "normal_1.5mb_a0.8_b0.4","normal_0.15mb_a0.8_b0.4"};
  int colors[4] = {kBlack,kRed, kBlue, kGreen};
  int markers[4] = {20, 21, 22,23};
  
  // Arrays to store graphs for all files
  TGraphErrors *gr_a2_Nch[4];
  TGraphErrors *gr_a3_Nch[4];
  
  // Process all three files
  for(int ifile = 0; ifile < 1; ifile++) {
    cout << "\n=== Processing file " << ifile+1 << ": " << filenames[ifile] << " ===" << endl;
    processFile(filenames[ifile], labels[ifile], colors[ifile], 
                &gr_a2_Nch[ifile], &gr_a3_Nch[ifile]);
    
    // Set marker styles
    gr_a2_Nch[ifile]->SetMarkerStyle(markers[ifile]);
    gr_a3_Nch[ifile]->SetMarkerStyle(markers[ifile]);
  }
  
  // Create output file and save individual graphs
  TFile *outfile = new TFile("longRange_flow_comparison.root", "recreate");
  
  for(int ifile = 0; ifile < 1; ifile++) {
    gr_a2_Nch[ifile]->Write(Form("gr_a2_Nch_%s", labels[ifile].Data()));
    gr_a3_Nch[ifile]->Write(Form("gr_a3_Nch_%s", labels[ifile].Data()));
  }
  
  // Create comparison plots
  
  // 1. a2 vs NCh comparison
  TCanvas *c_a2_Nch = new TCanvas("c_a2_Nch", "a2 vs Nch Comparison", 800, 600);
  
  TH2D *hframe_Nch = new TH2D("hframe_Nch", ";N_{ch};a_{n}", 100, 0, 80, 100, 0, 0.0072);
  hframe_Nch->Draw();
  
  TLegend *leg_Nch = new TLegend(0.15, 0.2, 0.4, 0.5);
  leg_Nch->SetTextSize(0.030);
  leg_Nch->SetBorderSize(0);
  for(int ifile = 0; ifile < 1; ifile++) {
    gr_a2_Nch[ifile]->Draw("same p");
    leg_Nch->AddEntry(gr_a2_Nch[ifile], labels[ifile], "p");
  }
  leg_Nch->Draw();
  c_a2_Nch->SaveAs("a2_Nch_comparison.png");
  
  // 2. a3 vs NCh comparison
  //TCanvas *c_a3_Nch = new TCanvas("c_a3_Nch", "a3 vs Nch Comparison", 800, 600);
  
  //TH2D *hframe_a3_Nch = new TH2D("hframe_a3_Nch", ";N_{ch};v_{3}{3}", 100, 0, 80, 100, 0, 0.0032);
  //hframe_a3_Nch->Draw();
  
  // TLegend *leg_a3_Nch = new TLegend(0.1, 0.6, 0.3, 0.9);
  // for(int ifile = 0; ifile < 1; ifile++) {
  //   gr_a3_Nch[ifile]->Draw("same p");
  //   leg_a3_Nch->AddEntry(gr_a3_Nch[ifile], labels[ifile], "p");
  // }
  // leg_a3_Nch->Draw();
  //c_a3_Nch->SaveAs("a3_Nch_comparison.png");
  
  // Save combined canvas with all four plots
  TCanvas *c_combined = new TCanvas("c_combined", "Flow Coefficients Comparison", 1600, 1200);
  c_combined->Divide(2, 1);
  c_combined->cd(1);
  hframe_Nch->Draw();
  for(int ifile = 0; ifile < 1; ifile++) {
    gr_a2_Nch[ifile]->Draw("same p");
  }
  leg_Nch->Draw();
  
  c_combined->cd(2);
  // hframe_a3_Nch->Draw();
  // for(int ifile = 0; ifile < 1; ifile++) {
  //   gr_a3_Nch[ifile]->Draw("same p");
  // }
  // leg_a3_Nch->Draw();
  
  c_combined->SaveAs("flow_comparison_all.png");
  
  outfile->Close();
}