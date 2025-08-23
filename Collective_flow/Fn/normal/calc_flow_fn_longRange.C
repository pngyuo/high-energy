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
    static TF1 *funRidge = nullptr;
    if (!funRidge) {
        funRidge = new TF1("funRidge", "[0]*(1 + 2*[1]*cos(x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x) + 2*[4]*cos(4*x))", -2.5 , 2.5 );
        funRidge->SetParameters(par[0], par[1], par[2], par[3], par[4]);
    } else {
        funRidge->SetParameters(par[0], par[1], par[2], par[3], par[4]);
    }
    return funRidge->Eval(x[0]);
}
  
void fitAndSave(TGraphErrors* graph, const TString& name, double& F_n, double& F_n_err) {
    TCanvas c;
    graph->Draw("AP");
    
    TF1 *fitFunc = new TF1("fitFunc", "[0]*(1 + [1]*x + [2]*x*x)", 
                          graph->GetXaxis()->GetXmin(), 
                          graph->GetXaxis()->GetXmax());
    
    fitFunc->SetParameters(0.01, 0.003, 0.003);
    
    TVirtualFitter::SetMaxIterations(10000);
    int fitStatus = graph->Fit("fitFunc", "QRS");
    
    F_n = fitFunc->GetParameter(1);
    F_n_err = fitFunc->GetParError(1);
    double p0 = fitFunc->GetParameter(0);
    double p0_err = fitFunc->GetParError(0);
    double p2 = fitFunc->GetParameter(2);
    double p2_err = fitFunc->GetParError(2);
    double chi2 = fitFunc->GetChisquare();
    int ndf = fitFunc->GetNDF();
    
    TLatex text;
    text.SetTextSize(0.03);
    text.SetTextAlign(12);
    text.DrawLatexNDC(0.25, 0.85, Form("Fit status: %d", fitStatus));
    text.DrawLatexNDC(0.25, 0.80, Form("F_{n} = %.4f #pm %.4f", F_n, F_n_err));
    text.DrawLatexNDC(0.25, 0.75, Form("p0 = %.4f #pm %.4f", p0, p0_err));
    text.DrawLatexNDC(0.25, 0.70, Form("p2 = %.4f #pm %.4f", p2, p2_err));
    text.DrawLatexNDC(0.25, 0.65, Form("#chi^{2}/NDF = %.2f/%d", chi2, ndf));
    text.DrawLatexNDC(0.25, 0.60, Form("Entries: %d", graph->GetN()));
    
    TGraph *rawPoints = new TGraph(*graph);
    rawPoints->SetMarkerStyle(3);
    rawPoints->SetMarkerColor(kGray);
    rawPoints->Draw("P SAME");
    
    c.SaveAs(Form("fit_results/%s.png", name.Data()));
    delete fitFunc;
    delete rawPoints;
}

void SavePlot(TH1D *hY, TString name) {
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
}

double getA2(TH1D *hY, double &a1, double &a2, double &a3, double &a4, double &G, double &chi2, int &ndf, double &a2err) {
    gYIntegral = hY->Integral("width");
    TF1 *fitFun = new TF1("fitFun", funTemplate, -0.5 * PI, 1.5 * PI, 5); 
    fitFun->SetParameters(0.1,0.1, 0.1, 0.1, 0.1); 
    hY->Fit("fitFun");
    a1 = fitFun->GetParameter(1);
    a2 = fitFun->GetParameter(2);
    a3 = fitFun->GetParameter(3);
    a4 = fitFun->GetParameter(4);
    a2err=fitFun->GetParError(2);
    G = gYIntegral / (2.0 * PI);
    chi2 = fitFun->GetChisquare();
    ndf = fitFun->GetNDF();
    delete fitFun;
    return a2;
}

double getA3(TH1D *hY, double &a1, double &a2, double &a3, double &a4, double &G, double &chi2, int &ndf, double &a3err) {
    gYIntegral = hY->Integral("width");
    TF1 *fitFun = new TF1("fitFun", funTemplate, -0.5 * PI, 1.5 * PI, 5); 
    fitFun->SetParameters(0.1,0.1, 0.1, 0.1, 0.1); 
    hY->Fit("fitFun");
    a1 = fitFun->GetParameter(1);
    a2 = fitFun->GetParameter(2);
    a3 = fitFun->GetParameter(3);
    a4 = fitFun->GetParameter(4);
    a3err=fitFun->GetParError(3);  // Fixed: was GetParError(2), should be GetParError(3)
    G = gYIntegral / (2.0 * PI);
    chi2 = fitFun->GetChisquare();
    ndf = fitFun->GetNDF();
    delete fitFun;
    return a3;
}

void drawProjection(TH2D *h2D, TString name) {
    TCanvas c;
    h2D->Draw("colz");
    c.SaveAs("projections/" + name + ".png");
}

void SaveSumRidgePlots(TH1D *hSumRidge, TString name) {
    TCanvas c;
    hSumRidge->Draw();
    c.SaveAs(Form("sumridge_plots/%s.png", name.Data()));
}

// Function to process a single file
void processFile(TFile* file_1, TFile* file_2, TString fileLabel, 
                double cent_an_f2[], double cent_cn_f2[], 
                double err_an_f2[], double err_cn_f2[],
                double cent_an_f3[], double cent_cn_f3[], 
                double err_an_f3[], double err_cn_f3[]) {
    
    TH1D *hYHigh = nullptr;
    TH1D *hYLow = nullptr;

    double G_high, G_low;
    double a1_high, a2_high, a3_high, a4_high;
    double a1_low, a2_low, a3_low, a4_low;
    double chi2_high, chi2_low;
    int ndf_high, ndf_low;

    TH1D *hTrigPtLow = (TH1D *)file_1->Get("hTrigPtLow");
    TH1D *hTrigPtHigh = (TH1D *)file_1->Get("hTrigPtHigh");
    TH1D *hNChMid = (TH1D *)file_1->Get("hNChMid");

    int nEtaBins = 50;
    double etaMin = -2.5;
    double etaMax = 2.5;
    double etaValues[nEtaBins];
    double c2_Eta[50];
    double c2_Eta_err[50];
    double a2_Eta_high[50];
    double a2_Eta_err[50];
    double a3_Eta_high[50];
    double a3_Eta_err[50];
    double c3_Eta[50];
    double c3_Eta_err[50];
    double etaBinWidths[nEtaBins];

    const int N_cent = 8;
    double Nsel_cut[N_cent + 1] = {10, 30, 40, 50, 60, 80, 100, 120, 150};
    double Nsel_val[N_cent] = {10};

    TH3D *hsame_cent[N_cent];
    TH3D *hmixed_cent[N_cent];
    TH1D *hTrigPt_cent[N_cent];
    TH3D *hsame_centlow = (TH3D *)file_1->Get("hDEtaDPhiTrigEtaSameEvent_Cent0");
    TH3D *hmixed_centlow = (TH3D *)file_1->Get("hDEtaDPhiTrigEtaMixEvent_Cent0");

    for (int icent = 1; icent < N_cent; icent++) {
        Nsel_val[icent] = (Nsel_cut[icent] + Nsel_cut[icent + 1]) / 2.0;
        double centMin = Nsel_cut[icent];
        double centMax = Nsel_cut[icent + 1];
        
        if (Nsel_val[icent] < 85) {
            if (icent == 0) {
                hTrigPt_cent[icent] = (TH1D *)file_1->Get(Form("hTrigPt_Cent%d", icent));
            } else {
                hsame_cent[icent] = (TH3D *)file_1->Get(Form("hDEtaDPhiTrigEtaSameEvent_Cent%d", icent));
                hmixed_cent[icent] = (TH3D *)file_1->Get(Form("hDEtaDPhiTrigEtaMixEvent_Cent%d", icent));
                hTrigPt_cent[icent] = (TH1D *)file_1->Get(Form("hTrigPt_Cent%d", icent));
            }
        } else {
            if (icent == 0) {
                hTrigPt_cent[icent] = (TH1D *)file_2->Get(Form("hTrigPt_Cent%d", icent));
            } else {
                hsame_cent[icent] = (TH3D *)file_2->Get(Form("hDEtaDPhiTrigEtaSameEvent_Cent%d", icent));
                hmixed_cent[icent] = (TH3D *)file_2->Get(Form("hDEtaDPhiTrigEtaMixEvent_Cent%d", icent));
                hTrigPt_cent[icent] = (TH1D *)file_2->Get(Form("hTrigPt_Cent%d", icent));
            }
        }
        
        // Calculate averageNTlow_cent
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

        for (int iEta = 0; iEta < nEtaBins; iEta++) {
            double eta = etaMin + iEta * (etaMax - etaMin) / nEtaBins;
            etaBinWidths[iEta] = (etaMax - etaMin) / nEtaBins/2;
            etaValues[iEta] = eta;
            int zBin = iEta + 6;
            
            hsame_centlow->GetZaxis()->SetRange(zBin, zBin);
            hmixed_centlow->GetZaxis()->SetRange(zBin, zBin);
            hsame_cent[icent]->GetZaxis()->SetRange(zBin, zBin);
            hmixed_cent[icent]->GetZaxis()->SetRange(zBin, zBin);

            TH2D *hsame_eta_high = (TH2D *)hsame_cent[icent]->Project3D("yx");
            TH2D *hmix_eta_high = (TH2D *)hmixed_cent[icent]->Project3D("yx");
            TH2D *hsame_eta_low = (TH2D *)hsame_centlow->Project3D("yx");
            TH2D *hmix_eta_low = (TH2D *)hmixed_centlow->Project3D("yx");

            TH1D *hsame_eta_high1D = sumRidge(hsame_eta_high, Form("hsame_eta_high_cent%d_eta%d_%s", icent, iEta, fileLabel.Data()));
            TH1D *hmix_eta_high1D = sumRidge(hmix_eta_high, Form("hmix_eta_high_cent%d_eta%d_%s", icent, iEta, fileLabel.Data()));
            TH1D *hsame_eta_low1D = sumRidge(hsame_eta_low, Form("hsame_eta_low_cent%d_eta%d_%s", icent, iEta, fileLabel.Data()));
            TH1D *hmix_eta_low1D = sumRidge(hmix_eta_low, Form("hmix_eta_low_eta%d_%s", iEta, fileLabel.Data()));

            if (hsame_eta_high->GetEntries() > 0 && hmix_eta_high->GetEntries() > 0) {
                hYHigh = getYFromHist(hsame_eta_high, hmix_eta_high, Form("YHigh_eta%d_%s", iEta, fileLabel.Data()), hTrigPtHigh->GetEntries());
                hYLow = getYFromHist(hsame_eta_low, hmix_eta_low, Form("YLow_eta%d_%s", iEta, fileLabel.Data()), hTrigPtLow->GetEntries());
                
                double a2higherr = 0, a2lowerr = 0;
                double a2_high = getA2(hYHigh, a1_high, a2_high, a3_high, a4_high, G_high, chi2_high, ndf_high, a2higherr);
                double a2_low = getA2(hYLow, a1_low, a2_low, a3_low, a4_low, G_low, chi2_low, ndf_low, a2lowerr);
                a2_Eta_high[iEta] = a2_high;
                a2_Eta_err[iEta] = a2higherr;
                double c2 = a2_high - a2_low * k;
                c2_Eta[iEta] = c2;
                double c2_err = sqrt(pow(hYHigh->GetFunction("fitFun")->GetParError(2), 2) + pow(hYLow->GetFunction("fitFun")->GetParError(2) * k, 2));
                c2_Eta_err[iEta] = c2_err;

                double a3higherr = 0, a3lowerr = 0;
                double a3_high = getA3(hYHigh, a1_high, a2_high, a3_high, a4_high, G_high, chi2_high, ndf_high, a3higherr);
                double a3_low = getA3(hYLow, a1_low, a2_low, a3_low, a4_low, G_low, chi2_low, ndf_low, a3lowerr);
                a3_Eta_high[iEta] = a3_high;
                a3_Eta_err[iEta] = a3higherr;
                double c3 = a3_high - a3_low * k;
                c3_Eta[iEta] = c3;
                double c3_err = sqrt(pow(hYHigh->GetFunction("fitFun")->GetParError(3), 2) + pow(hYLow->GetFunction("fitFun")->GetParError(3) * k, 2));
                c3_Eta_err[iEta] = c3_err;
                
                delete hYHigh;
                delete hYLow;
            }
            delete hsame_eta_high;
            delete hmix_eta_high;
            delete hsame_eta_low;
            delete hmix_eta_low;
            delete hsame_eta_high1D;
            delete hmix_eta_high1D;
            delete hsame_eta_low1D;
            delete hmix_eta_low1D;
        }

        // Fit a2 and a3
        TGraphErrors *gr_an_f2 = new TGraphErrors(nEtaBins, etaValues, a2_Eta_high, etaBinWidths, a2_Eta_err);
        TGraphErrors *gr_an_f3 = new TGraphErrors(nEtaBins, etaValues, a3_Eta_high, etaBinWidths, a3_Eta_err);
        double Fn_cent_an_f2, Fn_err_an_f2, Fn_cent_an_f3, Fn_err_an_f3;
        fitAndSave(gr_an_f2, Form("an_f2_cent%d_%s", icent, fileLabel.Data()), Fn_cent_an_f2, Fn_err_an_f2);
        fitAndSave(gr_an_f3, Form("an_f3_cent%d_%s", icent, fileLabel.Data()), Fn_cent_an_f3, Fn_err_an_f3);
        cent_an_f2[icent] = Fn_cent_an_f2;
        err_an_f2[icent] = Fn_err_an_f2;
        cent_an_f3[icent] = Fn_cent_an_f3;
        err_an_f3[icent] = Fn_err_an_f3;

        // Fit c2 and c3
        TGraphErrors *gr_cn_f2 = new TGraphErrors(nEtaBins, etaValues, c2_Eta, etaBinWidths, c2_Eta_err);
        TGraphErrors *gr_cn_f3 = new TGraphErrors(nEtaBins, etaValues, c3_Eta, etaBinWidths, c3_Eta_err);
        double Fn_cent_cn_f2, Fn_err_cn_f2, Fn_cent_cn_f3, Fn_err_cn_f3;
        fitAndSave(gr_cn_f2, Form("cn_f2_cent%d_%s", icent, fileLabel.Data()), Fn_cent_cn_f2, Fn_err_cn_f2);
        fitAndSave(gr_cn_f3, Form("cn_f3_cent%d_%s", icent, fileLabel.Data()), Fn_cent_cn_f3, Fn_err_cn_f3);
        cent_cn_f2[icent] = Fn_cent_cn_f2;
        err_cn_f2[icent] = Fn_err_cn_f2;
        cent_cn_f3[icent] = Fn_cent_cn_f3;
        err_cn_f3[icent] = Fn_err_cn_f3;

        delete gr_an_f2; delete gr_an_f3; delete gr_cn_f2; delete gr_cn_f3;
    }
}

void calc_flow_fn_longRange() {
    // Define the three files to process
    TString fileNames[3] = {
        "hist_ampt_normal_1.5mb_Decorr_yuhao.root",
        "hist_ampt_normal_0.15mb_Decorr_yuhao.root", 
        "hist_ampt_normal_1.5mb_a_0.8_b_0.4_Decorr_yuhao.root"
    };
    
    TString fileLabels[3] = {"1p5mb", "0p15mb", "1p5mb_a0p8_b0p4"};
    
    const int N_cent = 8;
    double Nsel_cut[N_cent + 1] = {10, 30, 40, 50, 60, 80, 100, 120, 150};
    double Nsel_val[N_cent] = {0};
    double Nsel_bin_width[N_cent] = {0};
    
    // Calculate centrality values and bin widths
    for (int icent = 1; icent < N_cent; icent++) {
        Nsel_val[icent] = (Nsel_cut[icent] + Nsel_cut[icent + 1]) / 2.0;
        Nsel_bin_width[icent] = (Nsel_cut[icent + 1] - Nsel_cut[icent]) / 2.0;
    }
    
    // Arrays to store results for all three files
    double cent_an_f2[3][N_cent], cent_cn_f2[3][N_cent];
    double err_an_f2[3][N_cent], err_cn_f2[3][N_cent];
    double cent_an_f3[3][N_cent], cent_cn_f3[3][N_cent];
    double err_an_f3[3][N_cent], err_cn_f3[3][N_cent];
    
    // Process each file
    for (int iFile = 0; iFile < 3; iFile++) {
        std::cout << "Processing file: " << fileNames[iFile] << std::endl;
        
        TFile *file_1 = new TFile(fileNames[iFile]);
        TFile *file_2 = new TFile(fileNames[iFile]); // Using same file for both
        
        processFile(file_1, file_2, fileLabels[iFile], 
                   cent_an_f2[iFile], cent_cn_f2[iFile], 
                   err_an_f2[iFile], err_cn_f2[iFile],
                   cent_an_f3[iFile], cent_cn_f3[iFile], 
                   err_an_f3[iFile], err_cn_f3[iFile]);
        
        file_1->Close();
        file_2->Close();
    }
    
    // Create combined plot
    TCanvas *cCombined = new TCanvas("cCombined", "Fn vs Nch - All Files", 1200, 800);
    
    // Define colors and markers for different files
    Color_t colors[3] = {kRed, kBlue, kGreen+2};
    Style_t markers_raw[3] = {20, 21, 22};
    Style_t markers_sub[3] = {24, 25, 26};
    
    TGraphErrors *gr_f2_raw[3], *gr_f2_sub[3];
    
    bool firstDraw = true;
    for (int iFile = 2; iFile < 3; iFile++) {
        gr_f2_raw[iFile] = new TGraphErrors(N_cent-1, &Nsel_val[1], &cent_an_f2[iFile][1], &Nsel_bin_width[1], &err_an_f2[iFile][1]);
        gr_f2_sub[iFile] = new TGraphErrors(N_cent-1, &Nsel_val[1], &cent_cn_f2[iFile][1], &Nsel_bin_width[1], &err_cn_f2[iFile][1]);
        
        // Set styles for raw data
        gr_f2_raw[iFile]->SetMarkerStyle(markers_raw[iFile]);
        gr_f2_raw[iFile]->SetMarkerColor(colors[iFile]);
        gr_f2_raw[iFile]->SetLineColor(colors[iFile]);
        
        // Set styles for subtracted data
        gr_f2_sub[iFile]->SetMarkerStyle(markers_sub[iFile]);
        gr_f2_sub[iFile]->SetMarkerColor(colors[iFile]);
        gr_f2_sub[iFile]->SetLineColor(colors[iFile]);
        
        if (firstDraw) {
            gr_f2_raw[iFile]->Draw("AP");
            gr_f2_raw[iFile]->GetXaxis()->SetTitle("N_{ch}");
            gr_f2_raw[iFile]->SetTitle("Long-range flow coefficients - All Files");
            gr_f2_raw[iFile]->GetXaxis()->SetRangeUser(0, 150);
            gr_f2_raw[iFile]->GetYaxis()->SetRangeUser(-0.05, 0.1);
            firstDraw = false;
        } else {
            gr_f2_raw[iFile]->Draw("P SAME");
        }
        gr_f2_sub[iFile]->Draw("P SAME");
    }
    
    // Create legend
    TLegend *leg = new TLegend(0.15, 0.15, 0.45, 0.35);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    
    for (int iFile = 2; iFile < 3; iFile++) {
        TString legendLabel;
        if (iFile == 0) legendLabel = "1.5mb";
        else if (iFile == 1) legendLabel = "0.15mb";
        else legendLabel = "1.5mb a=0.8 b=0.4";
        
        leg->AddEntry(gr_f2_raw[iFile], Form("f2^{raw} %s", legendLabel.Data()), "p");
        leg->AddEntry(gr_f2_sub[iFile], Form("f2^{sub} %s", legendLabel.Data()), "p");
    }
    leg->Draw();
    
    cCombined->SaveAs("Fn_1.5_a0.8_b0.4.png");
    
    // Save results to ROOT file
    TFile *outfile = new TFile("longRange_flow_combined.root", "recreate");
    
    for (int iFile = 0; iFile < 3; iFile++) {
        gr_f2_raw[iFile]->Write(Form("gr_f2_raw_%s", fileLabels[iFile].Data()));
        gr_f2_sub[iFile]->Write(Form("gr_f2_sub_%s", fileLabels[iFile].Data()));
    }
    
    outfile->Close();
    
    // Clean up
    for (int iFile = 0; iFile < 3; iFile++) {
        delete gr_f2_raw[iFile];
        delete gr_f2_sub[iFile];
    }
    
    delete cCombined;
}