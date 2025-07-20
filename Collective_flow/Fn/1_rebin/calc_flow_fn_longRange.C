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
    TH1D *htmp1 = (TH1D*)hin->ProjectionX("htmp1", hin->GetYaxis()->FindBin(-11.0), hin->GetYaxis()->FindBin(-1.0));
    TH1D *htmp2 = (TH1D*)hin->ProjectionX("htmp2", hin->GetYaxis()->FindBin(1.0), hin->GetYaxis()->FindBin(11.0)); 
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
        funRidge = new TF1("funRidge", "[0]*(1 + 2*[1]*cos(x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x) + 2*[4]*cos(4*x))", -0.5*PI, 1.5*PI);
        funRidge->SetParameters(par[0], par[1], par[2], par[3], par[4]);
    } else {
        funRidge->SetParameters(par[0], par[1], par[2], par[3], par[4]);
    }
    return funRidge->Eval(x[0]);
}
  
void fitAndSave(TGraphErrors* graph, const TString& name, double& F_n, double& F_n_err) {
    TCanvas c;
    graph->Draw("AP");
    
    TF1 *fitFunc = new TF1("fitFunc", "[0]*(1 + [1]*x + [2]*x*x)", -2.5,2.5);
    fitFunc->SetParameters(0.1, 0.1, -0.1);
    graph->Fit("fitFunc", "Q");
    
    F_n = fitFunc->GetParameter(1);
    F_n_err = fitFunc->GetParError(1);
    
    TLatex text;
    text.SetTextSize(0.04);
    text.DrawLatexNDC(0.25, 0.85, Form("F_{n} = %.4f #pm %.4f", F_n, F_n_err));
    
    c.SaveAs(Form("fit_results/%s.png", name.Data()));
    delete fitFunc;
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
    a3err=fitFun->GetParError(3);
    G = gYIntegral / (2.0 * PI);
    chi2 = fitFun->GetChisquare();
    ndf = fitFun->GetNDF();
    delete fitFun;
    return a3;
}

void calc_flow_fn_longRange() {
    TFile *file_1 = new TFile("hist_outputallFSI_new.root");
    TFile *file_2 = new TFile("hist_outputallFSI_new.root");

    TH1D *hYHigh = nullptr;
    TH1D *hYLow = nullptr;

    double G_high, G_low;
    double a1_high, a2_high, a3_high, a4_high;
    double a1_low, a2_low, a3_low, a4_low;
    double chi2_high, chi2_low;
    int ndf_high, ndf_low;
    double totalChargedParticleslow, totalChargedParticleshigh;
    double nEventslow, nEventshigh;
    double averageNThigh, averageNTlow;

    TH1D *hTrigPtLow = (TH1D *)file_1->Get("hTrigPtLow");
    TH1D *hTrigPtHigh = (TH1D *)file_1->Get("hTrigPtHigh");
    TH1D *hNChMid = (TH1D *)file_1->Get("hNChMid");

    int nEtaBins = 25;
    double etaMin = -2.5;
    double etaMax = 2.5;
    double etaValues[nEtaBins];
    double c2_Eta[25];
    double c2_Eta_err[25];
    double a2_Eta_high[25];
    double a2_Eta_err[25];
    double a3_Eta_high[25];
    double a3_Eta_err[25];
    double c3_Eta[25];
    double c3_Eta_err[25];
    double etaBinWidths[nEtaBins];

    const int N_cent = 8;
    double Nsel_cut[N_cent + 1] = {10, 30, 40, 50, 60, 80, 100, 120, 150};
    double Nsel_val[N_cent] = {0};
    double Nsel_bin_width[N_cent] = {0};

    TH3D *hsame_cent[N_cent];
    TH3D *hmixed_cent[N_cent];
    TH1D *hTrigPt_cent[N_cent];
    double cent_an_f2[N_cent] = {0}, cent_an_f3[N_cent] = {0};
    double cent_cn_f2[N_cent] = {0}, cent_cn_f3[N_cent] = {0};
    double err_an_f2[N_cent] = {0}, err_an_f3[N_cent] = {0};
    double err_cn_f2[N_cent] = {0}, err_cn_f3[N_cent] = {0};
    TH3D *hsame_centlow = (TH3D *)file_1->Get("hDEtaDPhiTrigEtaSameEvent_Cent0");
    TH3D *hmixed_centlow = (TH3D *)file_1->Get("hDEtaDPhiTrigEtaMixEvent_Cent0");

    for (int icent = 1; icent < N_cent; icent++) {
        Nsel_val[icent] = (Nsel_cut[icent] + Nsel_cut[icent + 1]) / 2.0;
        Nsel_bin_width[icent] = (Nsel_cut[icent + 1] - Nsel_cut[icent]) / 2.0;
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

            if (hsame_eta_high->GetEntries() > 0 && hmix_eta_high->GetEntries() > 0) {
                hYHigh = getYFromHist(hsame_eta_high, hmix_eta_high, Form("YHigh_eta%d", iEta), hTrigPtHigh->GetEntries());
                hYLow = getYFromHist(hsame_eta_low, hmix_eta_low, Form("YLow_eta%d", iEta), hTrigPtLow->GetEntries());
                double a2higherr = 0;
                double a2lowerr = 0;
                double a2_high = getA2(hYHigh, a1_high, a2_high, a3_high, a4_high, G_high, chi2_high, ndf_high, a2higherr);
                double a2_low = getA2(hYLow, a1_low, a2_low, a3_low, a4_low, G_low, chi2_low, ndf_low, a2lowerr);
                a2_Eta_high[iEta] = a2_high;
                a2_Eta_err[iEta] = a2higherr;
                double c2 = a2_high - a2_low * k;
                c2_Eta[iEta] = c2;
                double c2_err = sqrt(pow(a2higherr, 2) + pow(a2lowerr * k, 2));
                c2_Eta_err[iEta] = c2_err;

                double a3higherr = 0;
                double a3lowerr = 0;
                double a3_high = getA3(hYHigh, a1_high, a2_high, a3_high, a4_high, G_high, chi2_high, ndf_high, a3higherr);
                double a3_low = getA3(hYLow, a1_low, a2_low, a3_low, a4_low, G_low, chi2_low, ndf_low, a3lowerr);
                a3_Eta_high[iEta] = a3_high;
                a3_Eta_err[iEta] = a3higherr;
                double c3 = a3_high - a3_low * k;
                c3_Eta[iEta] = c3;
                double c3_err = sqrt(pow(a3higherr, 2) + pow(a3lowerr * k, 2));
                c3_Eta_err[iEta] = c3_err;
                delete hYHigh;
                delete hYLow;
            }
            delete hsame_eta_high;
            delete hmix_eta_high;
            delete hsame_eta_low;
            delete hmix_eta_low;
        }

        // ==================== 添加 rebin 操作 ====================
        const int rebinFactor = 1; // 每5个点合并为1个点
        const int new_nBins = nEtaBins / rebinFactor;
        
        double etaValues_rebin[new_nBins];
        double x_err_rebin[new_nBins];
        double a2_Eta_high_rebin[new_nBins], a2_Eta_err_rebin[new_nBins];
        double a3_Eta_high_rebin[new_nBins], a3_Eta_err_rebin[new_nBins];
        double c2_Eta_rebin[new_nBins], c2_Eta_err_rebin[new_nBins];
        double c3_Eta_rebin[new_nBins], c3_Eta_err_rebin[new_nBins];
        
        for (int new_i = 0; new_i < new_nBins; new_i++) {
            int startBin = new_i * rebinFactor;
            double sum_a2 = 0, sum_a3 = 0, sum_c2 = 0, sum_c3 = 0;
            double sum_err2_a2 = 0, sum_err2_a3 = 0, sum_err2_c2 = 0, sum_err2_c3 = 0;
            double sum_x = 0;
            
            for (int j = 0; j < rebinFactor; j++) {
                int old_i = startBin + j;
                if (old_i >= nEtaBins) continue;
                
                sum_a2 += a2_Eta_high[old_i];
                sum_err2_a2 += pow(a2_Eta_err[old_i], 2);
                sum_a3 += a3_Eta_high[old_i];
                sum_err2_a3 += pow(a3_Eta_err[old_i], 2);
                sum_c2 += c2_Eta[old_i];
                sum_err2_c2 += pow(c2_Eta_err[old_i], 2);
                sum_c3 += c3_Eta[old_i];
                sum_err2_c3 += pow(c3_Eta_err[old_i], 2);
                sum_x += etaValues[old_i];
            }
            
            etaValues_rebin[new_i] = sum_x / rebinFactor;
            a2_Eta_high_rebin[new_i] = sum_a2 / rebinFactor;
            a2_Eta_err_rebin[new_i] = sqrt(sum_err2_a2) / rebinFactor;
            a3_Eta_high_rebin[new_i] = sum_a3 / rebinFactor;
            a3_Eta_err_rebin[new_i] = sqrt(sum_err2_a3) / rebinFactor;
            c2_Eta_rebin[new_i] = sum_c2 / rebinFactor;
            c2_Eta_err_rebin[new_i] = sqrt(sum_err2_c2) / rebinFactor;
            c3_Eta_rebin[new_i] = sum_c3 / rebinFactor;
            c3_Eta_err_rebin[new_i] = sqrt(sum_err2_c3) / rebinFactor;
            
            x_err_rebin[new_i] = rebinFactor * (etaMax - etaMin) / nEtaBins / 2.0;
        }
        
        // ==================== 使用 rebin 后的数据 ====================
        TGraphErrors *gr_an_f2 = new TGraphErrors(
            new_nBins, etaValues_rebin, a2_Eta_high_rebin, 
            x_err_rebin, a2_Eta_err_rebin
        );
        TGraphErrors *gr_an_f3 = new TGraphErrors(
            new_nBins, etaValues_rebin, a3_Eta_high_rebin, 
            x_err_rebin, a3_Eta_err_rebin
        );
        TGraphErrors *gr_cn_f2 = new TGraphErrors(
            new_nBins, etaValues_rebin, c2_Eta_rebin, 
            x_err_rebin, c2_Eta_err_rebin
        );
        TGraphErrors *gr_cn_f3 = new TGraphErrors(
            new_nBins, etaValues_rebin, c3_Eta_rebin, 
            x_err_rebin, c3_Eta_err_rebin
        );
        
        double Fn_cent_an_f2, Fn_err_an_f2, Fn_cent_an_f3, Fn_err_an_f3;
        fitAndSave(gr_an_f2, Form("an_f2_cent%d_rebin", icent), Fn_cent_an_f2, Fn_err_an_f2);
        fitAndSave(gr_an_f3, Form("an_f3_cent%d_rebin", icent), Fn_cent_an_f3, Fn_err_an_f3);
        cent_an_f2[icent] = Fn_cent_an_f2;
        err_an_f2[icent] = Fn_err_an_f2;
        cent_an_f3[icent] = Fn_cent_an_f3;
        err_an_f3[icent] = Fn_err_an_f3;

        double Fn_cent_cn_f2, Fn_err_cn_f2, Fn_cent_cn_f3, Fn_err_cn_f3;
        fitAndSave(gr_cn_f2, Form("cn_f2_cent%d_rebin", icent), Fn_cent_cn_f2, Fn_err_cn_f2);
        fitAndSave(gr_cn_f3, Form("cn_f3_cent%d_rebin", icent), Fn_cent_cn_f3, Fn_err_cn_f3);
        cent_cn_f2[icent] = Fn_cent_cn_f2;
        err_cn_f2[icent] = Fn_err_cn_f2;
        cent_cn_f3[icent] = Fn_cent_cn_f3;
        err_cn_f3[icent] = Fn_err_cn_f3;
        
        delete gr_an_f2; 
        delete gr_an_f3;
        delete gr_cn_f2; 
        delete gr_cn_f3;
    }

    TCanvas *cCombined = new TCanvas("cCombined", "Fn vs Nch", 800, 600);
    cCombined->SetLogy();
    TGraphErrors *gr_f2_raw = new TGraphErrors(N_cent, Nsel_val, cent_an_f2, Nsel_bin_width, err_an_f2);
    TGraphErrors *gr_f2_sub = new TGraphErrors(N_cent, Nsel_val, cent_cn_f2, Nsel_bin_width, err_cn_f2); 
    TGraphErrors *gr_f3_raw = new TGraphErrors(N_cent-1, &Nsel_val[1], &cent_an_f3[1], 0, &err_an_f3[1]);
    TGraphErrors *gr_f3_sub = new TGraphErrors(N_cent-1, &Nsel_val[1], &cent_cn_f3[1], 0, &err_cn_f3[1]);

    gr_f2_raw->SetMarkerStyle(20); gr_f2_raw->SetMarkerColor(kRed);
    //gr_f3_raw->SetMarkerStyle(21); gr_f3_raw->SetMarkerColor(kBlue);
    gr_f2_sub->SetMarkerStyle(22); gr_f2_sub->SetMarkerColor(kGreen+2);
    //gr_f3_sub->SetMarkerStyle(23); gr_f3_sub->SetMarkerColor(kMagenta);

    gr_f2_raw->Draw("AP");
    //gr_f3_raw->Draw("P SAME");
    gr_f2_sub->Draw("P SAME");
    //gr_f3_sub->Draw("P SAME");

    gr_f2_raw->GetXaxis()->SetTitle("N_{ch}");
    gr_f2_raw->GetYaxis()->SetTitle("F_{n}");
    gr_f2_raw->SetTitle("Long-range flow coefficients");
    gr_f2_raw->GetYaxis()->SetRangeUser(0.002, 2);

    TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(gr_f2_raw, "f2^{raw}", "p");
    //leg->AddEntry(gr_f3_raw, "f3^{raw}", "p");
    leg->AddEntry(gr_f2_sub, "f2^{sub}", "p");
    //leg->AddEntry(gr_f3_sub, "f3^{sub}", "p");
    leg->Draw();

    cCombined->SaveAs("Fn_rebin.png");

    TFile *outfile = new TFile("longRange_flow.root", "recreate");
    gr_f2_raw->Write("gr_f2_raw");
    //gr_f3_raw->Write("gr_f3_raw");
    gr_f2_sub->Write("gr_f2_sub");
    //gr_f3_sub->Write("gr_f3_sub");
    outfile->Close();
}