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
  
void SavePlot(TH1D *hY, TString name, TString multiplicity) {
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
    c0.SaveAs(Form("test_figure/dphi_template_%s_%s.png", multiplicity.Data(), name.Data()));
    delete text;
}

double getAn(TH1D *hY, double &a1, double &a2, double &a3, double &a4, double &G, double &chi2, int &ndf, double &a2err, TString name, TString multiplicity) {
    gYIntegral = hY->Integral("width");
    TF1 *fitFun = new TF1("fitFun", funTemplate, -0.5 * PI, 1.5 * PI, 5);
    fitFun->SetParameters(0.1, 0.1, 0.1, 0.1, 0.1);
    hY->Fit("fitFun", "Q");
    a1 = fitFun->GetParameter(1);
    a2 = fitFun->GetParameter(2);
    a3 = fitFun->GetParameter(3);
    a4 = fitFun->GetParameter(4);
    a2err = fitFun->GetParError(2);
    G = gYIntegral / (2.0 * PI);
    chi2 = fitFun->GetChisquare();
    ndf = fitFun->GetNDF();
    SavePlot(hY, name, multiplicity);
    double result = fitFun->GetParameter(2);
    delete fitFun;
    return result;
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

void calc_flow_ancn_longRange() {
    TFile *file_1 = new TFile("hist_ampt_3q_1.5mb_Decorr_yuhao.root");
    TFile *file_2 = new TFile("hist_ampt_3q_1.5mb_Decorr_yuhao.root");

    // 创建输出文件，提前打开以便保存每个eta bin的结果
    TFile *outfile = new TFile("longRange_flow.root", "RECREATE");
    
    TH1D *hYHigh = nullptr;
    TH1D *hYLow = nullptr;

    double G_high, G_low;
    double a1_high, a2_high, a3_high, a4_high;
    double a1_low, a2_low, a3_low, a4_low;
    double chi2_high, chi2_low;
    int ndf_high, ndf_low;
    double totalChargedParticleslow = 0, totalChargedParticleshigh = 0;
    double nEventslow = 0, nEventshigh = 0;
    double averageNThigh, averageNTlow;

    TH1D *hTrigPtLow = (TH1D *)file_1->Get("hTrigPtLow");
    TH1D *hTrigPtHigh = (TH1D *)file_1->Get("hTrigPtHigh");

    TH1D *chargedParticlesCount = (TH1D*)file_1->Get("hNChMid");
    double totalEvents = chargedParticlesCount->GetEntries();
    std::cout << "Total events in h2D: " << totalEvents << std::endl;
    for (int bin = 1; bin <= chargedParticlesCount->GetNbinsX(); bin++) {
        double binCenter = chargedParticlesCount->GetXaxis()->GetBinCenter(bin);
        double binContent = chargedParticlesCount->GetBinContent(bin);
        
        if(binCenter < 30 && binCenter > 10){
            totalChargedParticleslow += binCenter * binContent;
            nEventslow += binContent;
        }
        else if(binCenter > 60){
            totalChargedParticleshigh += binCenter * binContent;
            nEventshigh += binContent;
        }
    }
    averageNTlow = nEventslow > 0 ? totalChargedParticleslow / nEventslow : 0;
    std::cout << "averageNTlow: " << averageNTlow << std::endl;
    averageNThigh = nEventshigh > 0 ? totalChargedParticleshigh / nEventshigh : 0;
    std::cout << "averageNThigh: " << averageNThigh << std::endl;
    double k = (nEventshigh > 0 && nEventslow > 0) ? averageNTlow / averageNThigh : 0;
    std::cout << "k = " << k << std::endl;

    TH3D *hDEtaDPhiTrigEtaSameEventHighMid = (TH3D *)file_1->Get("hDEtaDPhiTrigEtaSameEventHighMid");
    TH3D *hDEtaDPhiTrigEtaMixEventHighMid = (TH3D *)file_1->Get("hDEtaDPhiTrigEtaMixEventHighMid");
    TH3D *hDEtaDPhiTrigEtaSameEventLowMid = (TH3D *)file_1->Get("hDEtaDPhiTrigEtaSameEventLowMid");
    TH3D *hDEtaDPhiTrigEtaMixEventLowMid = (TH3D *)file_1->Get("hDEtaDPhiTrigEtaMixEventLowMid");

    int nEtaBins = 50; 
    double etaMin = -2.5;
    double etaMax = 2.5;
    double etaValues[nEtaBins];
    double c2_Eta[50];
    double c2_Eta_err[50];
    double a2_Eta_high[50];
    double a2_Eta_err_high[50];
    double a2_Eta_low[50];
    double a2_Eta_err_low[50];
    double etaBinHalfWidths[nEtaBins];  // 改为存储半宽

    for (int iEta = 0; iEta < nEtaBins; iEta++) {
        double eta = etaMin + (iEta + 0.5) * (etaMax - etaMin) / nEtaBins;
        double binWidth = (etaMax - etaMin) / nEtaBins;
        etaBinHalfWidths[iEta] = binWidth / 2.0;  // 计算半宽
        etaValues[iEta] = eta;
        int zBin = iEta + 6;
        
        hDEtaDPhiTrigEtaSameEventHighMid->GetZaxis()->SetRange(zBin, zBin);
        hDEtaDPhiTrigEtaMixEventHighMid->GetZaxis()->SetRange(zBin, zBin);
        hDEtaDPhiTrigEtaSameEventLowMid->GetZaxis()->SetRange(zBin, zBin);
        hDEtaDPhiTrigEtaMixEventLowMid->GetZaxis()->SetRange(zBin, zBin);
        
        TH2D *hsame_eta_high = (TH2D*)hDEtaDPhiTrigEtaSameEventHighMid->Project3D("yx");
        TH2D *hmix_eta_high = (TH2D*)hDEtaDPhiTrigEtaMixEventHighMid->Project3D("yx");
        TH2D *hsame_eta_low = (TH2D*)hDEtaDPhiTrigEtaSameEventLowMid->Project3D("yx");
        TH2D *hmix_eta_low = (TH2D*)hDEtaDPhiTrigEtaMixEventLowMid->Project3D("yx");

        drawProjection(hsame_eta_high, Form("hsame_eta_high_eta%d", iEta));
        drawProjection(hmix_eta_high, Form("hmix_eta_high_eta%d", iEta));
        drawProjection(hsame_eta_low, Form("hsame_eta_low_eta%d", iEta));
        drawProjection(hmix_eta_low, Form("hmix_eta_low_eta%d", iEta));
        
        TH1D *hsame_eta_high1D = sumRidge(hsame_eta_high, Form("hsame_eta_high_eta%d", iEta));
        SaveSumRidgePlots(hsame_eta_high1D, Form("hsame_eta_high_eta%d", iEta));

        TH1D *hmix_eta_high1D = sumRidge(hmix_eta_high, Form("hmix_eta_high_eta%d", iEta));
        SaveSumRidgePlots(hmix_eta_high1D, Form("hmix_eta_high_eta%d", iEta));

        TH1D *hsame_eta_low1D = sumRidge(hsame_eta_low, Form("hsame_eta_low_eta%d", iEta));
        SaveSumRidgePlots(hsame_eta_low1D, Form("hsame_eta_low_eta%d", iEta));

        TH1D *hmix_eta_low1D = sumRidge(hmix_eta_low, Form("hmix_eta_low_eta%d", iEta));
        SaveSumRidgePlots(hmix_eta_low1D, Form("hmix_eta_low_eta%d", iEta));

        if (hsame_eta_high->GetEntries() > 0 && hmix_eta_high->GetEntries() > 0) {
            hYHigh = getYFromHist(hsame_eta_high, hmix_eta_high, Form("YHigh_eta%d", iEta), hTrigPtHigh->GetEntries());
            hYLow = getYFromHist(hsame_eta_low, hmix_eta_low, Form("YLow_eta%d", iEta), hTrigPtLow->GetEntries());
            
            double a2higherr = 0, a2lowerr = 0;
            double a2_high = getAn(hYHigh, a1_high, a2_high, a3_high, a4_high, G_high, chi2_high, ndf_high, a2higherr, Form("eta%d", iEta), "high");
            double a2_low = getAn(hYLow, a1_low, a2_low, a3_low, a4_low, G_low, chi2_low, ndf_low, a2lowerr, Form("eta%d", iEta), "low");
            
            std::cout << "a2_high: " << a2_high << std::endl;
            std::cout << "a2_low: " << a2_low << std::endl;
            
            a2_Eta_high[iEta] = a2_high;
            a2_Eta_err_high[iEta] = a2higherr;
            a2_Eta_low[iEta] = a2_low;
            a2_Eta_err_low[iEta] = a2lowerr;
            
            double c2 = a2_high - a2_low * k;
            std::cout << "c2: " << c2 << std::endl;
            c2_Eta[iEta] = c2;
            
            double c2_err = sqrt(a2higherr*a2higherr + (a2lowerr*k)*(a2lowerr*k));
            c2_Eta_err[iEta] = c2_err;
        }

        // 保存每个eta bin的直方图到输出文件
        if (hYHigh) {
            outfile->cd();
            hYHigh->Write(Form("hYHigh_eta%d", iEta));
            delete hYHigh;
            hYHigh = nullptr;
        }
        
        if (hYLow) {
            outfile->cd();
            hYLow->Write(Form("hYLow_eta%d", iEta));
            delete hYLow;
            hYLow = nullptr;
        }

        delete hsame_eta_high;
        delete hmix_eta_high;
        delete hsame_eta_low;
        delete hmix_eta_low;
    }

    int rebinFactor = 2;
    int newNEtaBins = nEtaBins / rebinFactor;
    double newEtaValues[newNEtaBins], newC2_Eta[newNEtaBins], newC2_Eta_err[newNEtaBins];
    double newA2_Eta_high[newNEtaBins], newA2_Eta_err_high[newNEtaBins];
    double newA2_Eta_low[newNEtaBins], newA2_Eta_err_low[newNEtaBins];
    double newEtaBinHalfWidths[newNEtaBins];  // 存储合并后的半宽
    
    double binWidth = (etaMax - etaMin) / nEtaBins;
    double rebinnedBinWidth = binWidth * rebinFactor;
    double rebinnedHalfWidth = rebinnedBinWidth / 2.0;

    for (int i = 0; i < newNEtaBins; i++) {
        int start = i * rebinFactor;
        double sumC2 = 0, sumC2Err2 = 0;
        double sumA2_high = 0, sumA2_high_Err2 = 0;
        double sumA2_low = 0, sumA2_low_Err2 = 0;
        double sumEta = 0;
        
        for (int j = 0; j < rebinFactor; j++) {
            int idx = start + j;
            if (idx >= nEtaBins) break;
            
            sumC2 += c2_Eta[idx];
            sumC2Err2 += c2_Eta_err[idx] * c2_Eta_err[idx];
            
            sumA2_high += a2_Eta_high[idx];
            sumA2_high_Err2 += a2_Eta_err_high[idx] * a2_Eta_err_high[idx];
            
            sumA2_low += a2_Eta_low[idx];
            sumA2_low_Err2 += a2_Eta_err_low[idx] * a2_Eta_err_low[idx];
            
            sumEta += etaValues[idx];
        }
        
        newC2_Eta[i] = sumC2 / rebinFactor;
        newC2_Eta_err[i] = sqrt(sumC2Err2) / rebinFactor;
        
        newA2_Eta_high[i] = sumA2_high / rebinFactor;
        newA2_Eta_err_high[i] = sqrt(sumA2_high_Err2) / rebinFactor;
        
        newA2_Eta_low[i] = sumA2_low / rebinFactor;
        newA2_Eta_err_low[i] = sqrt(sumA2_low_Err2) / rebinFactor;
        
        newEtaValues[i] = sumEta / rebinFactor;
        newEtaBinHalfWidths[i] = rebinnedHalfWidth;  // 使用合并后的半宽
    }

    // 创建图形时使用半宽作为横坐标误差
    TGraphErrors *gr_c2_rebinned = new TGraphErrors(newNEtaBins, newEtaValues, newC2_Eta, newEtaBinHalfWidths, newC2_Eta_err);
    TGraphErrors *gr_a2_high_rebinned = new TGraphErrors(newNEtaBins, newEtaValues, newA2_Eta_high, newEtaBinHalfWidths, newA2_Eta_err_high);
    TGraphErrors *gr_a2_low_rebinned = new TGraphErrors(newNEtaBins, newEtaValues, newA2_Eta_low, newEtaBinHalfWidths, newA2_Eta_err_low);
    
    TCanvas *c = new TCanvas("c", "c", 800, 600);
    c->cd();
    
    gr_c2_rebinned->SetTitle("Long-range flow coefficients vs #eta");
    gr_c2_rebinned->GetXaxis()->SetTitle("#eta");
    gr_c2_rebinned->GetYaxis()->SetTitle("Coefficient value");
    gr_c2_rebinned->GetXaxis()->SetRangeUser(etaMin, etaMax);
    gr_c2_rebinned->GetYaxis()->SetRangeUser(0.000, 0.008);
    
    gr_c2_rebinned->SetMarkerStyle(20);
    gr_c2_rebinned->SetMarkerColor(kRed);
    gr_c2_rebinned->SetLineColor(kRed);
    gr_c2_rebinned->Draw("AP");
    
    gr_a2_high_rebinned->SetMarkerStyle(21);
    gr_a2_high_rebinned->SetMarkerColor(kBlue);
    gr_a2_high_rebinned->SetLineColor(kBlue);
    gr_a2_high_rebinned->Draw("P SAME");
    
    gr_a2_low_rebinned->SetMarkerStyle(22);
    gr_a2_low_rebinned->SetMarkerColor(kGreen+2);
    gr_a2_low_rebinned->SetLineColor(kGreen+2);
    gr_a2_low_rebinned->Draw("P SAME");
    
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(gr_c2_rebinned, "c_{2}", "p");
    leg->AddEntry(gr_a2_high_rebinned, "a_{2}^{high}", "p");
    leg->AddEntry(gr_a2_low_rebinned, "a_{2}^{low}", "p");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();
    
    c->SaveAs("c2_and_a2_vs_eta_rebinned.png");

    // 保存最终结果到输出文件
    outfile->cd();
    gr_c2_rebinned->Write("gr_c2_rebinned");
    gr_a2_high_rebinned->Write("gr_a2_high_rebinned");
    gr_a2_low_rebinned->Write("gr_a2_low_rebinned");
    
    c->Write("c2_and_a2_vs_eta_rebinned");
    
    outfile->Close();
    
    // 清理内存
    delete c;
    delete leg;
    delete file_1;
    delete file_2;
}