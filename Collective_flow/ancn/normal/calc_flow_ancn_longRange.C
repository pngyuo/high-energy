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
    return funRidge->Eval(dPhi);
}
  
void SavePlot(TH1D *hY, TString name) {
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
    c0.SaveAs("test_figure/dphi_template_" + name + ".png");
}

double getAn(TH1D *hY, double &a1, double &a2, double &a3, double &a4, double &G, double &chi2, int &ndf, double &a2err) {
    gYIntegral = hY->Integral("width");
    TF1 *fitFun = new TF1("fitFun", funTemplate, -0.5 * PI, 1.5 * PI, 4); // 参数数量改为4
    fitFun->SetParameters(0.1,0.1, 0.1, 0.1, 0.1); // 初始化a1-a4
    hY->Fit("fitFun");
    a1 = fitFun->GetParameter(1);
    a2 = fitFun->GetParameter(2);
    a3 = fitFun->GetParameter(3);
    a4 = fitFun->GetParameter(4);
    a2err=fitFun->GetParError(2);
    G = gYIntegral / (2.0 * PI); // 计算最终的G值
    chi2 = fitFun->GetChisquare();
    ndf = fitFun->GetNDF();
    SavePlot(hY, hY->GetName());
    return fitFun->GetParameter(2);
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
    TFile *file_1 = new TFile("hist_ampt_normal_1.5mb_Decorr_yuhao.root");
    TFile *file_2 = new TFile("hist_ampt_normal_1.5mb_Decorr_yuhao.root");

    TH1D *hYHigh = nullptr;
    TH1D *hYLow = nullptr;

    double G_high, G_low;
    double a1_high, a2_high, a3_high, a4_high;
    double a1_low, a2_low, a3_low, a4_low;
    double chi2_high, chi2_low;
    int ndf_high, ndf_low;
    double totalChargedParticleslow,totalChargedParticleshigh;
    double nEventslow,nEventshigh;
    double averageNThigh,averageNTlow;

    TH1D *hTrigPtLow = (TH1D *)file_1->Get("hTrigPtLow");
    TH1D *hTrigPtHigh = (TH1D *)file_1->Get("hTrigPtHigh");

    TH1D *chargedParticlesCount = (TH1D*)file_1->Get("hNChMid");
    double totalEvents = chargedParticlesCount->GetEntries();
    std::cout << "Total events in h2D: " << totalEvents << std::endl;
    double totalChargedParticles = 0; // 加权总和
    double nEvents = 0;  // 总的事件数
    for (int bin = 1; bin <= 2000; bin++) {
        double binCenter = chargedParticlesCount->GetXaxis()->GetBinCenter(bin); // 获取当前 bin 中心的 NT 值
        if(binCenter<30&&binCenter>10){
        double binContent = chargedParticlesCount->GetBinContent(bin);  // 获取当前 bin 的分布数量
        totalChargedParticleslow += binCenter * binContent;  // 加权累加
        nEventslow += binContent;  // 累加事件数
        }
        else if(binCenter>60){
            double binContent = chargedParticlesCount->GetBinContent(bin);  // 获取当前 bin 的分布数量
            totalChargedParticleshigh += binCenter * binContent;  // 加权累加
            nEventshigh += binContent;  // 累加事件数
            }
    }
    averageNTlow = totalChargedParticleslow / nEventslow;  // 计算加权平均NT
    cout<<"averageNTlow: "<<averageNTlow<<endl;
    averageNThigh = totalChargedParticleshigh / nEventshigh;  // 计算加权平均NT
    cout<<"averageNThigh: "<<averageNThigh<<endl;
    // double M_low = hTrigPtLow->GetMean();
    // double M = hTrigPtHigh->GetMean();
    double k = averageNTlow / averageNThigh;
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
    double a2_Eta_err[50];
    double etaBinWidths[nEtaBins];

    for (int iEta = 0; iEta < nEtaBins; iEta++) {
        double eta = etaMin + iEta * (etaMax - etaMin) / nEtaBins;
        etaBinWidths[iEta] = (etaMax - etaMin) / nEtaBins;
        etaValues[iEta] = eta;
        // 计算对应的Z bin（从6开始）
        int zBin = iEta + 6; // TH3D的Z bin索引从1开始，有效数据在6到55
        hDEtaDPhiTrigEtaSameEventHighMid->GetZaxis()->SetRange(zBin, zBin);
        hDEtaDPhiTrigEtaMixEventHighMid->GetZaxis()->SetRange(zBin, zBin);
        hDEtaDPhiTrigEtaSameEventLowMid->GetZaxis()->SetRange(zBin, zBin);
        hDEtaDPhiTrigEtaMixEventLowMid->GetZaxis()->SetRange(zBin, zBin);
        TH2D *hsame_eta_high = (TH2D*)hDEtaDPhiTrigEtaSameEventHighMid->Project3D("yx");
        TH2D *hmix_eta_high = (TH2D*)hDEtaDPhiTrigEtaMixEventHighMid->Project3D("yx");
        TH2D *hsame_eta_low = (TH2D*)hDEtaDPhiTrigEtaSameEventLowMid->Project3D("yx");
        TH2D *hmix_eta_low = (TH2D*)hDEtaDPhiTrigEtaMixEventLowMid->Project3D("yx");

        // 绘制投影图
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
            double a2higherr=0;
            double a2lowerr=0;
            double a2_high = getAn(hYHigh, a1_high, a2_high, a3_high, a4_high, G_high, chi2_high, ndf_high,a2higherr);
            double a2_low = getAn(hYLow, a1_low, a2_low, a3_low, a4_low, G_low, chi2_low, ndf_low,a2lowerr);
            cout<<"a2_high"<<a2_high<<endl;
            cout<<"a2_low"<<a2_low<<endl;
            a2_Eta_high[iEta] = a2_high;
            a2_Eta_err[iEta] = a2higherr;
            double c2 = a2_high - a2_low * k;
            cout<<"c2"<<c2<<endl;
            c2_Eta[iEta] = c2;
            double c2_err = sqrt(pow(hYHigh->GetFunction("fitFun")->GetParError(1), 2) + pow(hYLow->GetFunction("fitFun")->GetParError(2) * k, 2));
            c2_Eta_err[iEta] = c2_err;
        }

        delete hsame_eta_high;
        delete hmix_eta_high;
        delete hsame_eta_low;
        delete hmix_eta_low;
    }

    TGraphErrors *gr_c2 = new TGraphErrors(nEtaBins, etaValues, c2_Eta, etaBinWidths, c2_Eta_err); // 如果你有误差值，可以传递误差数组
    TGraphErrors *gr_a2_high = new TGraphErrors(nEtaBins, etaValues, a2_Eta_high, etaBinWidths, a2_Eta_err); // 高多重性下的 a2 和 eta 的关系
    
    int rebinFactor = 2; // 合并的bin数
    int newNEtaBins = nEtaBins / rebinFactor;
    double newEtaValues[newNEtaBins], newC2_Eta[newNEtaBins], newC2_Eta_err[newNEtaBins];
    double newA2_Eta_high[newNEtaBins], newA2_Eta_err[newNEtaBins], newEtaBinWidths[newNEtaBins];
    double deltaEta = (etaMax - etaMin) / nEtaBins;

    for (int i = 0; i < newNEtaBins; i++) {
        int start = i * rebinFactor;
        double sumC2 = 0, sumC2Err2 = 0, sumA2 = 0, sumA2Err2 = 0, sumEta = 0;
        for (int j = 0; j < rebinFactor; j++) {
            int idx = start + j;
            if (idx >= nEtaBins) break;
            sumC2 += c2_Eta[idx];
            sumC2Err2 += pow(c2_Eta_err[idx], 2);
            sumA2 += a2_Eta_high[idx];
            sumA2Err2 += pow(a2_Eta_err[idx], 2);
            sumEta += etaValues[idx];
        }
        newC2_Eta[i] = sumC2 / rebinFactor;
        newC2_Eta_err[i] = sqrt(sumC2Err2) / rebinFactor;
        newA2_Eta_high[i] = sumA2 / rebinFactor;
        newA2_Eta_err[i] = sqrt(sumA2Err2) / rebinFactor;
        newEtaValues[i] = sumEta / rebinFactor;
        newEtaBinWidths[i] = deltaEta * rebinFactor;
    }

    // 使用合并后的数据创建新图形
    TGraphErrors *gr_c2_rebinned = new TGraphErrors(newNEtaBins, newEtaValues, newC2_Eta, newEtaBinWidths, newC2_Eta_err);
    TGraphErrors *gr_a2_high_rebinned = new TGraphErrors(newNEtaBins, newEtaValues, newA2_Eta_high, newEtaBinWidths, newA2_Eta_err);
    TCanvas *c = new TCanvas("c", "c");
    c->cd();
    gr_c2_rebinned->SetTitle("c_2 and a_2 vs Eta (Rebinned)");
    gr_c2_rebinned->GetXaxis()->SetTitle("Eta");
    gr_c2_rebinned->GetYaxis()->SetTitle("c_2 / a_2");
    gr_c2_rebinned->GetXaxis()->SetRangeUser(etaMin, etaMax);
    gr_c2_rebinned->GetYaxis()->SetRangeUser(0.006, 0.012); // 调整范围
    gr_c2_rebinned->SetMarkerColor(kRed);
    gr_c2_rebinned->SetLineColor(kRed);
    gr_c2_rebinned->Draw("AP");

    gr_a2_high_rebinned->SetMarkerColor(kBlue);
    gr_a2_high_rebinned->SetLineColor(kBlue);
    gr_a2_high_rebinned->Draw("P SAME");

    c->SaveAs("c2_and_a2_vs_eta_rebinned.png");

    TFile *outfile = new TFile("longRange_flow.root", "recreate");
    if (hYHigh) hYHigh->Write("hYHigh");
    if (hYLow) hYLow->Write("hYLow");
    outfile->Close();
}