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

double getAn(TH1D *hY, double &a1_out, double &a2_out, double &a3_out, double &a4_out, 
             double &G_out, double &chi2_out, int &ndf_out, double &a2err_out, 
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
    G_out = gYIntegral / (2.0 * PI);
    chi2_out = fitFun->GetChisquare();
    ndf_out = fitFun->GetNDF();
    
    SavePlot(hY, name, multiplicity, fileLabel);
    double result = fitFun->GetParameter(2);
    delete fitFun;
    return result;
}

// 处理单个文件的函数
void processFile(const char* fileName, const char* fileLabel, TFile* outfile, 
                 TGraphErrors** gr_c2, TGraphErrors** gr_a2_high, TGraphErrors** gr_a2_low) {
    
    std::cout << "Processing file: " << fileName << std::endl;
    
    TFile *file = new TFile(fileName);
    if (!file || file->IsZombie()) {
        std::cout << "Error opening file: " << fileName << std::endl;
        return;
    }

    TH1D *hYHigh = nullptr;
    TH1D *hYLow = nullptr;

    // 为每个eta bin分别定义变量，避免冲突
    double G_high_temp, G_low_temp;
    double a1_high_temp, a2_high_temp, a3_high_temp, a4_high_temp;
    double a1_low_temp, a2_low_temp, a3_low_temp, a4_low_temp;
    double chi2_high_temp, chi2_low_temp;
    int ndf_high_temp, ndf_low_temp;
    
    double totalChargedParticleslow = 0, totalChargedParticleshigh = 0;
    double nEventslow = 0, nEventshigh = 0;
    double averageNThigh, averageNTlow;

    TH1D *hTrigPtLow = (TH1D *)file->Get("hTrigPtLow");
    TH1D *hTrigPtHigh = (TH1D *)file->Get("hTrigPtHigh");
    
    TH1D *chargedParticlesCount = (TH1D*)file->Get("hNChMid");
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

    TH3D *hDEtaDPhiTrigEtaSameEventHighMid = (TH3D *)file->Get("hDEtaDPhiTrigEtaSameEventHighMid");
    TH3D *hDEtaDPhiTrigEtaMixEventHighMid = (TH3D *)file->Get("hDEtaDPhiTrigEtaMixEventHighMid");
    TH3D *hDEtaDPhiTrigEtaSameEventLowMid = (TH3D *)file->Get("hDEtaDPhiTrigEtaSameEventLowMid");
    TH3D *hDEtaDPhiTrigEtaMixEventLowMid = (TH3D *)file->Get("hDEtaDPhiTrigEtaMixEventLowMid");

    int nEtaBins = 25; 
    double etaMin = -2.5;
    double etaMax = 2.5;
    double etaValues[nEtaBins];
    double c2_Eta[25];
    double c2_Eta_err[25];
    double a2_Eta_high[25];
    double a2_Eta_err_high[25];
    double a2_Eta_low[25];
    double a2_Eta_err_low[25];
    double etaBinHalfWidths[nEtaBins];

    for (int iEta = 0; iEta < nEtaBins; iEta++) {
        double eta = etaMin + (iEta + 0.5) * (etaMax - etaMin) / nEtaBins;
        double binWidth = (etaMax - etaMin) / nEtaBins;
        etaBinHalfWidths[iEta] = binWidth / 2.0;
        etaValues[iEta] = eta;
        int zBinmin = 2*iEta + 6;
        int zBinmax = 2*iEta + 7;
        
        std::cout << "Processing eta bin " << iEta << ", eta = " << eta << std::endl;
        double zMin = hDEtaDPhiTrigEtaSameEventHighMid->GetZaxis()->GetBinLowEdge(zBinmin);
        double zMax = hDEtaDPhiTrigEtaSameEventHighMid->GetZaxis()->GetBinUpEdge(zBinmax);
        std::cout << "  Trigger eta range: " << zMin << " to " << zMax << std::endl;
        hDEtaDPhiTrigEtaSameEventHighMid->GetZaxis()->SetRange(zBinmin, zBinmax);
        hDEtaDPhiTrigEtaMixEventHighMid->GetZaxis()->SetRange(zBinmin, zBinmax);
        hDEtaDPhiTrigEtaSameEventLowMid->GetZaxis()->SetRange(zBinmin, zBinmax);
        hDEtaDPhiTrigEtaMixEventLowMid->GetZaxis()->SetRange(zBinmin, zBinmax);

        TH2D *hsame_eta_high = (TH2D*)hDEtaDPhiTrigEtaSameEventHighMid->Project3D("yx");
        TH2D *hmix_eta_high = (TH2D*)hDEtaDPhiTrigEtaMixEventHighMid->Project3D("yx");
        TH2D *hsame_eta_low = (TH2D*)hDEtaDPhiTrigEtaSameEventLowMid->Project3D("yx");
        TH2D *hmix_eta_low = (TH2D*)hDEtaDPhiTrigEtaMixEventLowMid->Project3D("yx");

        TH2D *hTrigPtEtaLow = (TH2D *)file->Get("hTrigPtEtaLow");
        TH2D *hTrigPtEtaHigh = (TH2D *)file->Get("hTrigPtEtaHigh");
        int bin_eta = hTrigPtEtaHigh->GetYaxis()->FindBin(eta);
        double Ntrig_high_bin = hTrigPtEtaHigh->Integral(5, hTrigPtEtaHigh->GetXaxis()->GetNbins(), bin_eta, bin_eta);
        double Ntrig_low_bin = hTrigPtEtaLow->Integral(5, hTrigPtEtaLow->GetXaxis()->GetNbins(), bin_eta, bin_eta);

        if (hsame_eta_high->GetEntries() > 0 && hmix_eta_high->GetEntries() > 0) {
            // 清理之前的直方图
            if (hYHigh) delete hYHigh;
            if (hYLow) delete hYLow;
            
            hYHigh = getYFromHist(hsame_eta_high, hmix_eta_high, Form("YHigh_eta%d_%s", iEta, fileLabel), 10);
            hYLow = getYFromHist(hsame_eta_low, hmix_eta_low, Form("YLow_eta%d_%s", iEta, fileLabel), 10);
            
            // 使用临时变量接收getAn的返回值
            double a2higherr_temp = 0, a2lowerr_temp = 0;
            
            // 修正：使用不同的变量名避免冲突
            double a2_result_high = getAn(hYHigh, a1_high_temp, a2_high_temp, a3_high_temp, a4_high_temp, 
                                         G_high_temp, chi2_high_temp, ndf_high_temp, a2higherr_temp, 
                                         Form("eta%d", iEta), "high", fileLabel);
            
            double a2_result_low = getAn(hYLow, a1_low_temp, a2_low_temp, a3_low_temp, a4_low_temp, 
                                        G_low_temp, chi2_low_temp, ndf_low_temp, a2lowerr_temp, 
                                        Form("eta%d", iEta), "low", fileLabel);
            
            std::cout << "File " << fileLabel << " - Eta bin " << iEta << " (eta=" << eta << ")" << std::endl;
            std::cout << "  a2_high: " << a2_result_high << std::endl;
            std::cout << "  a2_low: " << a2_result_low << std::endl;
            
            // 将结果存储到数组中
            a2_Eta_high[iEta] = a2_result_high;
            a2_Eta_err_high[iEta] = a2higherr_temp;
            a2_Eta_low[iEta] = a2_result_low;
            a2_Eta_err_low[iEta] = a2lowerr_temp;
            
            double c2 = a2_result_high - a2_result_low * k;
            std::cout << "  c2: " << c2 << std::endl;
            c2_Eta[iEta] = c2;
            
            double c2_err = sqrt(a2higherr_temp*a2higherr_temp + (a2lowerr_temp*k)*(a2lowerr_temp*k));
            c2_Eta_err[iEta] = c2_err;
        }

                // 保存每个eta bin的直方图到输出文件
        if (hYHigh) {
            outfile->cd();
            hYHigh->Write(Form("hYHigh_eta%d_%s", iEta, fileLabel));
            delete hYHigh;
            hYHigh = nullptr;
        }
        
        if (hYLow) {
            outfile->cd();
            hYLow->Write(Form("hYLow_eta%d_%s", iEta, fileLabel));
            delete hYLow;
            hYLow = nullptr;
        }

        delete hsame_eta_high;
        delete hmix_eta_high;
        delete hsame_eta_low;
        delete hmix_eta_low;
    }

    // 创建图形时使用半宽作为横坐标误差
    *gr_c2 = new TGraphErrors(nEtaBins, etaValues, c2_Eta, etaBinHalfWidths, c2_Eta_err);
    *gr_a2_high = new TGraphErrors(nEtaBins, etaValues, a2_Eta_high, etaBinHalfWidths, a2_Eta_err_high);
    *gr_a2_low = new TGraphErrors(nEtaBins, etaValues, a2_Eta_low, etaBinHalfWidths, a2_Eta_err_low);

    // 保存单个文件的结果到输出文件
    outfile->cd();
    (*gr_c2)->Write(Form("gr_c2_%s", fileLabel));
    (*gr_a2_high)->Write(Form("gr_a2_high_%s", fileLabel));
    (*gr_a2_low)->Write(Form("gr_a2_low_%s", fileLabel));

    // 清理内存
    if (hYHigh) delete hYHigh;
    if (hYLow) delete hYLow;
    
    file->Close();
    delete file;
}

void calc_flow_ancn_longRange_files() {
    // 文件列表
    const char* fileNames[4] = {
        "hist_ampt_3q_1.5mb_Decorr_yuhao.root",
        "hist_ampt_3q_0.15mb_Decorr_yuhao.root", 
        "hist_ampt_3q_1.5mb_a_0.8_b_0.4_Decorr_yuhao.root",
        "hist_ampt_normal_0.15mb_a_0.8_b_0.4_Decorr_yuhao.root"
    };
    
    const char* fileLabels[4] = {"3q_1.5mb", "3q_0.15mb", "normal_1.5mb_a08_b04","normal_0.15mb_a08_b04"};
    const char* legendLabels[4] = {"3q 1.5mb", "3q 0.15mb", "normal 1.5mb a=0.8 b=0.4", "normal 0.15mb a=0.8 b=0.4"};
    
    // 创建输出文件
    TFile *outfile = new TFile("longRange_flow_multifile.root", "RECREATE");
    
    // 存储三个文件的结果
    TGraphErrors *gr_c2[4];
    TGraphErrors *gr_a2_high[4];
    TGraphErrors *gr_a2_low[4];
    
    // 处理三个文件
    for (int iFile = 0; iFile < 2; iFile++) {
        processFile(fileNames[iFile], fileLabels[iFile], outfile, 
                   &gr_c2[iFile], &gr_a2_high[iFile], &gr_a2_low[iFile]);
    }
    
    // 创建综合图表
    // TCanvas *c = new TCanvas("c", "c", 1600, 800);
    // c->Divide(4, 1); // 分为4个子画布
    
    // // 颜色和标记样式设置
    int colors[4] = {kBlack,kRed, kBlue, kGreen+2};
    int markerStyles[4] = {20, 21, 22,23};
    
    // // 子画布1: c2对比
    // c->cd(1);
    // gPad->SetLeftMargin(0.15);
    // gPad->SetBottomMargin(0.15);
    
     bool firstDraw = true;
    // for (int iFile = 0; iFile < 2; iFile++) {
    //     if (gr_c2[iFile]) {
    //         gr_c2[iFile]->SetTitle("c_{2} vs #eta");
    //         gr_c2[iFile]->GetXaxis()->SetTitle("#eta");
    //         gr_c2[iFile]->GetYaxis()->SetTitle("c_{2}");
    //         gr_c2[iFile]->GetXaxis()->SetRangeUser(-2.5, 2.5);
    //         gr_c2[iFile]->GetYaxis()->SetRangeUser(0.0025, 0.01);
            
    //         gr_c2[iFile]->SetMarkerStyle(markerStyles[iFile]);
    //         gr_c2[iFile]->SetMarkerColor(colors[iFile]);
    //         gr_c2[iFile]->SetLineColor(colors[iFile]);
            
    //         if (firstDraw) {
    //             gr_c2[iFile]->Draw("AP");
    //             firstDraw = false;
    //         } else {
    //             gr_c2[iFile]->Draw("P SAME");
    //         }
    //     }
    // }
    
    // TLegend *leg1 = new TLegend(0.2, 0.7, 0.8, 0.9);
    // for (int iFile = 0; iFile < 2; iFile++) {
    //     if (gr_c2[iFile]) {
    //         leg1->AddEntry(gr_c2[iFile], legendLabels[iFile], "p");
    //     }
    // }
    // leg1->SetBorderSize(0);
    // leg1->SetFillStyle(0);
    // leg1->Draw();
    
    // // 子画布2: a2_high对比
    // c->cd(2);
    // gPad->SetLeftMargin(0.15);
    // gPad->SetBottomMargin(0.15);
    
    // firstDraw = true;
    // for (int iFile = 0; iFile < 2; iFile++) {
    //     if (gr_a2_high[iFile]) {
    //         gr_a2_high[iFile]->SetTitle("a_{2}^{high} vs #eta");
    //         gr_a2_high[iFile]->GetXaxis()->SetTitle("#eta");
    //         gr_a2_high[iFile]->GetYaxis()->SetTitle("a_{2}^{high}");
    //         gr_a2_high[iFile]->GetXaxis()->SetRangeUser(-2.5, 2.5);
            
    //         gr_a2_high[iFile]->SetMarkerStyle(markerStyles[iFile]);
    //         gr_a2_high[iFile]->SetMarkerColor(colors[iFile]);
    //         gr_a2_high[iFile]->SetLineColor(colors[iFile]);
            
    //         if (firstDraw) {
    //             gr_a2_high[iFile]->Draw("AP");
    //             firstDraw = false;
    //         } else {
    //             gr_a2_high[iFile]->Draw("P SAME");
    //         }
    //     }
    // }
    
    // TLegend *leg2 = new TLegend(0.2, 0.7, 0.8, 0.9);
    // for (int iFile = 0; iFile < 2; iFile++) {
    //     if (gr_a2_high[iFile]) {
    //         leg2->AddEntry(gr_a2_high[iFile], legendLabels[iFile], "p");
    //     }
    // }
    // leg2->SetBorderSize(0);
    // leg2->SetFillStyle(0);
    // leg2->Draw();
    
    // // 子画布3: a2_low对比
    // c->cd(3);
    // gPad->SetLeftMargin(0.15);
    // gPad->SetBottomMargin(0.15);
    
    // firstDraw = true;
    // for (int iFile = 0; iFile < 2; iFile++) {
    //     if (gr_a2_low[iFile]) {
    //         gr_a2_low[iFile]->SetTitle("a_{2}^{low} vs #eta");
    //         gr_a2_low[iFile]->GetXaxis()->SetTitle("#eta");
    //         gr_a2_low[iFile]->GetYaxis()->SetTitle("a_{2}^{low}");
    //         gr_a2_low[iFile]->GetXaxis()->SetRangeUser(-2.5, 2.5);
            
    //         gr_a2_low[iFile]->SetMarkerStyle(markerStyles[iFile]);
    //         gr_a2_low[iFile]->SetMarkerColor(colors[iFile]);
    //         gr_a2_low[iFile]->SetLineColor(colors[iFile]);
            
    //         if (firstDraw) {
    //             gr_a2_low[iFile]->Draw("AP");
    //             firstDraw = false;
    //         } else {
    //             gr_a2_low[iFile]->Draw("P SAME");
    //         }
    //     }
    // }
    
    // TLegend *leg3 = new TLegend(0.2, 0.7, 0.8, 0.9);
    // for (int iFile = 0; iFile < 2; iFile++) {
    //     if (gr_a2_low[iFile]) {
    //         leg3->AddEntry(gr_a2_low[iFile], legendLabels[iFile], "p");
    //     }
    // }
    // leg3->SetBorderSize(0);
    // leg3->SetFillStyle(0);
    // leg3->Draw();
    
    // // 子画布4: 所有系数在一起显示（9条线）
    // c->cd(4);
    // gPad->SetLeftMargin(0.15);
    // gPad->SetBottomMargin(0.15);
    
    // // 线型设置：实线、虚线、点线
    // int lineStyles[3] = {1, 2, 3}; // 实线、虚线、点线
    // int fillStyles[3] = {1001, 0, 0}; // 实心、空心、空心
    
    // firstDraw = true;
    // double yMin = 1e10, yMax = -1e10;
    
    // // 先找到y轴范围
    // for (int iFile = 0; iFile < 2; iFile++) {
    //     if (gr_c2[iFile]) {
    //         for (int i = 0; i < gr_c2[iFile]->GetN(); i++) {
    //             double x, y;
    //             gr_c2[iFile]->GetPoint(i, x, y);
    //             if (y < yMin) yMin = y;
    //             if (y > yMax) yMax = y;
    //         }
    //     }
    //     if (gr_a2_high[iFile]) {
    //         for (int i = 0; i < gr_a2_high[iFile]->GetN(); i++) {
    //             double x, y;
    //             gr_a2_high[iFile]->GetPoint(i, x, y);
    //             if (y < yMin) yMin = y;
    //             if (y > yMax) yMax = y;
    //         }
    //     }
    //     if (gr_a2_low[iFile]) {
    //         for (int i = 0; i < gr_a2_low[iFile]->GetN(); i++) {
    //             double x, y;
    //             gr_a2_low[iFile]->GetPoint(i, x, y);
    //             if (y < yMin) yMin = y;
    //             if (y > yMax) yMax = y;
    //         }
    //     }
    // }
    
    // // 画c2
    // for (int iFile = 0; iFile < 2; iFile++) {
    //     if (gr_c2[iFile]) {
    //         gr_c2[iFile]->SetTitle("Flow Coefficients vs #eta");
    //         gr_c2[iFile]->GetXaxis()->SetTitle("#eta");
    //         gr_c2[iFile]->GetYaxis()->SetTitle("Coefficient value");
    //         gr_c2[iFile]->GetXaxis()->SetRangeUser(-2.5, 2.5);
    //         gr_c2[iFile]->GetYaxis()->SetRangeUser(yMin*0.8, yMax*1.2);
            
    //         gr_c2[iFile]->SetMarkerStyle(markerStyles[iFile]);
    //         gr_c2[iFile]->SetMarkerColor(colors[iFile]);
    //         gr_c2[iFile]->SetLineColor(colors[iFile]);
    //         gr_c2[iFile]->SetLineStyle(lineStyles[iFile]);
    //         gr_c2[iFile]->SetFillStyle(fillStyles[iFile]);
            
    //         if (firstDraw) {
    //             gr_c2[iFile]->Draw("AP");
    //             firstDraw = false;
    //         } else {
    //             gr_c2[iFile]->Draw("PL SAME");
    //         }
    //     }
    // }
    
    // // 画a2_high
    // for (int iFile = 0; iFile < 2; iFile++) {
    //     if (gr_a2_high[iFile]) {
    //         // 使用更深的颜色变体
    //         int darkColors[3] = {kRed+2, kBlue+2, kGreen+4};
            
    //         gr_a2_high[iFile]->SetMarkerStyle(markerStyles[iFile]);
    //         gr_a2_high[iFile]->SetMarkerColor(darkColors[iFile]);
    //         gr_a2_high[iFile]->SetLineColor(darkColors[iFile]);
    //         gr_a2_high[iFile]->SetLineStyle(lineStyles[iFile]);
    //         gr_a2_high[iFile]->SetFillStyle(fillStyles[iFile]);
            
    //         gr_a2_high[iFile]->Draw("P SAME");
    //     }
    // }
    
    // // 画a2_low
    // for (int iFile = 0; iFile < 2; iFile++) {
    //     if (gr_a2_low[iFile]) {
    //         // 使用更浅的颜色变体
    //         int lightColors[3] = {kRed-7, kBlue-7, kGreen-7};
            
    //         gr_a2_low[iFile]->SetMarkerStyle(markerStyles[iFile]);
    //         gr_a2_low[iFile]->SetMarkerColor(lightColors[iFile]);
    //         gr_a2_low[iFile]->SetLineColor(lightColors[iFile]);
    //         gr_a2_low[iFile]->SetLineStyle(lineStyles[iFile]);
    //         gr_a2_low[iFile]->SetFillStyle(fillStyles[iFile]);
            
    //         gr_a2_low[iFile]->Draw("P SAME");
    //     }
    // }
    
    // // 创建综合图例
    // TLegend *leg4 = new TLegend(0.15, 0.75, 0.85, 0.95);
    // leg4->SetNColumns(3); // 3列显示
    
    // for (int iFile = 0; iFile < 2; iFile++) {
    //     if (gr_c2[iFile]) {
    //         leg4->AddEntry(gr_c2[iFile], Form("c_{2} %s", legendLabels[iFile]), "p");
    //     }
    // }
    // for (int iFile = 0; iFile < 2; iFile++) {
    //     if (gr_a2_high[iFile]) {
    //         int darkColors[3] = {kRed+2, kBlue+2, kGreen+4};
    //         gr_a2_high[iFile]->SetLineColor(darkColors[iFile]); // 确保图例颜色正确
    //         leg4->AddEntry(gr_a2_high[iFile], Form("a_{2}^{high} %s", legendLabels[iFile]), "p");
    //     }
    // }
    // for (int iFile = 0; iFile < 2; iFile++) {
    //     if (gr_a2_low[iFile]) {
    //         int lightColors[3] = {kRed-7, kBlue-7, kGreen-7};
    //         gr_a2_low[iFile]->SetLineColor(lightColors[iFile]); // 确保图例颜色正确
    //         leg4->AddEntry(gr_a2_low[iFile], Form("a_{2}^{low} %s", legendLabels[iFile]), "p");
    //     }
    // }
    
    // leg4->SetBorderSize(0);
    // leg4->SetFillStyle(0);
    // leg4->SetTextSize(0.025); // 调小字体以适应更多条目
    // leg4->Draw();
    
    // c->SaveAs("c2_and_a2_vs_eta_multifile_comparison.png");
    
    // 保存每个子画布为单独的文件
// 子画布1 - c2对比
TCanvas *c1_individual = new TCanvas("c1_individual", "c2 vs eta", 800, 600);
c1_individual->cd();
gPad->SetLeftMargin(0.15);
gPad->SetBottomMargin(0.15);

firstDraw = true;
for (int iFile = 0; iFile < 2; iFile++) {
    if (gr_c2[iFile]) {
        gr_c2[iFile]->SetTitle("c_{2} vs #eta");
        gr_c2[iFile]->GetXaxis()->SetTitle("#eta");
        gr_c2[iFile]->GetYaxis()->SetTitle("c_{2}");
        gr_c2[iFile]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        gr_c2[iFile]->GetYaxis()->SetRangeUser(0.0025, 0.01);
        
        // 添加颜色和标记样式设置
        gr_c2[iFile]->SetMarkerStyle(markerStyles[iFile]);
        gr_c2[iFile]->SetMarkerColor(colors[iFile]);
        gr_c2[iFile]->SetLineColor(colors[iFile]);
        
        if (firstDraw) {
            gr_c2[iFile]->Draw("AP");
            firstDraw = false;
        } else {
            gr_c2[iFile]->Draw("P SAME");
        }
    }
}

TLegend *leg1_copy = new TLegend(0.2, 0.7, 0.8, 0.9);
for (int iFile = 0; iFile < 2; iFile++) {
    if (gr_c2[iFile]) {
        leg1_copy->AddEntry(gr_c2[iFile], legendLabels[iFile], "p");
    }
}
leg1_copy->SetBorderSize(0);
leg1_copy->SetFillStyle(0);
leg1_copy->Draw();
c1_individual->SaveAs("subplot1_c2_vs_eta.png");

// 子画布2 - a2_high对比
TCanvas *c2_individual = new TCanvas("c2_individual", "a2_high vs eta", 800, 600);
c2_individual->cd();
gPad->SetLeftMargin(0.15);
gPad->SetBottomMargin(0.15);

firstDraw = true;
for (int iFile = 0; iFile < 2; iFile++) {
    if (gr_a2_high[iFile]) {
        gr_a2_high[iFile]->SetTitle("a_{2}^{high} vs #eta");
        gr_a2_high[iFile]->GetXaxis()->SetTitle("#eta");
        gr_a2_high[iFile]->GetYaxis()->SetTitle("a_{2}^{high}");
        gr_a2_high[iFile]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        gr_a2_high[iFile]->GetYaxis()->SetRangeUser(0.0015, 0.012);
        
        // 添加颜色和标记样式设置
        gr_a2_high[iFile]->SetMarkerStyle(markerStyles[iFile]);
        gr_a2_high[iFile]->SetMarkerColor(colors[iFile]);
        gr_a2_high[iFile]->SetLineColor(colors[iFile]);
        
        if (firstDraw) {
            gr_a2_high[iFile]->Draw("AP");
            firstDraw = false;
        } else {
            gr_a2_high[iFile]->Draw("P SAME");
        }
    }
}

TLegend *leg2_copy = new TLegend(0.2, 0.7, 0.8, 0.9);
for (int iFile = 0; iFile < 2; iFile++) {
    if (gr_a2_high[iFile]) {
        leg2_copy->AddEntry(gr_a2_high[iFile], legendLabels[iFile], "p");
    }
}
leg2_copy->SetBorderSize(0);
leg2_copy->SetFillStyle(0);
leg2_copy->Draw();
c2_individual->SaveAs("subplot2_a2_high_vs_eta.png");

// 子画布3 - a2_low对比
TCanvas *c3_individual = new TCanvas("c3_individual", "a2_low vs eta", 800, 600);
c3_individual->cd();
gPad->SetLeftMargin(0.15);
gPad->SetBottomMargin(0.15);

firstDraw = true;
for (int iFile = 0; iFile < 2; iFile++) {
    if (gr_a2_low[iFile]) {
        gr_a2_low[iFile]->SetTitle("a_{2}^{low} vs #eta");
        gr_a2_low[iFile]->GetXaxis()->SetTitle("#eta");
        gr_a2_low[iFile]->GetYaxis()->SetTitle("a_{2}^{low}");
        gr_a2_low[iFile]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        gr_a2_low[iFile]->GetYaxis()->SetRangeUser(0.0015, 0.012);
        
        // 添加颜色和标记样式设置
        gr_a2_low[iFile]->SetMarkerStyle(markerStyles[iFile]);
        gr_a2_low[iFile]->SetMarkerColor(colors[iFile]);
        gr_a2_low[iFile]->SetLineColor(colors[iFile]);
        
        if (firstDraw) {
            gr_a2_low[iFile]->Draw("AP");
            firstDraw = false;
        } else {
            gr_a2_low[iFile]->Draw("P SAME");
        }
    }
}

TLegend *leg3_copy = new TLegend(0.2, 0.7, 0.8, 0.9);
for (int iFile = 0; iFile < 2; iFile++) {
    if (gr_a2_low[iFile]) {
        leg3_copy->AddEntry(gr_a2_low[iFile], legendLabels[iFile], "p");
    }
}
leg3_copy->SetBorderSize(0);
leg3_copy->SetFillStyle(0);
leg3_copy->Draw();
c3_individual->SaveAs("subplot3_a2_low_vs_eta.png");
    
    // 子画布4 - 所有系数在一起显示（9条线）
    TCanvas *c4_individual = new TCanvas("c4_individual", "All Flow Coefficients vs eta", 800, 600);
    c4_individual->cd();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    
    // 线型设置：实线、虚线、点线
    int lineStyles_copy[3] = {1, 2, 3}; // 实线、虚线、点线
    int fillStyles_copy[3] = {1001, 0, 0}; // 实心、空心、空心
    
    firstDraw = true;
    double yMin_copy = 1e10, yMax_copy = -1e10;
    
    // 先找到y轴范围
    for (int iFile = 0; iFile < 2; iFile++) {
        if (gr_c2[iFile]) {
            for (int i = 0; i < gr_c2[iFile]->GetN(); i++) {
                double x, y;
                gr_c2[iFile]->GetPoint(i, x, y);
                if (y < yMin_copy) yMin_copy = y;
                if (y > yMax_copy) yMax_copy = y;
            }
        }
        if (gr_a2_high[iFile]) {
            for (int i = 0; i < gr_a2_high[iFile]->GetN(); i++) {
                double x, y;
                gr_a2_high[iFile]->GetPoint(i, x, y);
                if (y < yMin_copy) yMin_copy = y;
                if (y > yMax_copy) yMax_copy = y;
            }
        }
        if (gr_a2_low[iFile]) {
            for (int i = 0; i < gr_a2_low[iFile]->GetN(); i++) {
                double x, y;
                gr_a2_low[iFile]->GetPoint(i, x, y);
                if (y < yMin_copy) yMin_copy = y;
                if (y > yMax_copy) yMax_copy = y;
            }
        }
    }
    
    // 画c2
    for (int iFile = 0; iFile < 2; iFile++) {
        if (gr_c2[iFile]) {
            gr_c2[iFile]->SetTitle("Flow Coefficients vs #eta");
            gr_c2[iFile]->GetXaxis()->SetTitle("#eta");
            gr_c2[iFile]->GetYaxis()->SetTitle("Coefficient value");
            gr_c2[iFile]->GetXaxis()->SetRangeUser(-2.5, 2.5);
            gr_c2[iFile]->GetYaxis()->SetRangeUser(0.005, 0.013);
            gr_c2[iFile]->GetYaxis()->SetRangeUser(yMin_copy*0.8, yMax_copy*1.2);
            
            gr_c2[iFile]->SetMarkerStyle(markerStyles[iFile]);
            gr_c2[iFile]->SetMarkerColor(colors[iFile]);
            gr_c2[iFile]->SetLineColor(colors[iFile]);
            gr_c2[iFile]->SetLineStyle(lineStyles_copy[iFile]);
            gr_c2[iFile]->SetFillStyle(fillStyles_copy[iFile]);
            
            if (firstDraw) {
                gr_c2[iFile]->Draw("AP");
                firstDraw = false;
            } else {
                gr_c2[iFile]->Draw("P SAME");
            }
        }
    }
    
    // 画a2_high
    for (int iFile = 0; iFile < 2; iFile++) {
        if (gr_a2_high[iFile]) {
            // 使用更深的颜色变体
            int darkColors_copy[3] = {kRed+2, kBlue+2, kGreen+4};
            
            gr_a2_high[iFile]->SetMarkerStyle(markerStyles[iFile]);
            gr_a2_high[iFile]->SetMarkerColor(darkColors_copy[iFile]);
            gr_a2_high[iFile]->SetLineColor(darkColors_copy[iFile]);
            gr_a2_high[iFile]->SetLineStyle(lineStyles_copy[iFile]);
            gr_a2_high[iFile]->SetFillStyle(fillStyles_copy[iFile]);
            
            gr_a2_high[iFile]->Draw("P SAME");
        }
    }
    
    // 画a2_low
    for (int iFile = 0; iFile < 2; iFile++) {
        if (gr_a2_low[iFile]) {
            // 使用更浅的颜色变体
            int lightColors_copy[3] = {kRed-7, kBlue-7, kGreen-7};
            
            gr_a2_low[iFile]->SetMarkerStyle(markerStyles[iFile]);
            gr_a2_low[iFile]->SetMarkerColor(lightColors_copy[iFile]);
            gr_a2_low[iFile]->SetLineColor(lightColors_copy[iFile]);
            gr_a2_low[iFile]->SetLineStyle(lineStyles_copy[iFile]);
            gr_a2_low[iFile]->SetFillStyle(fillStyles_copy[iFile]);
            
            gr_a2_low[iFile]->Draw("P SAME");
        }
    }
    
    // 创建综合图例
    TLegend *leg4_copy = new TLegend(0.15, 0.75, 0.85, 0.95);
    leg4_copy->SetNColumns(3); // 3列显示
    
    for (int iFile = 0; iFile < 2; iFile++) {
        if (gr_c2[iFile]) {
            leg4_copy->AddEntry(gr_c2[iFile], Form("c_{2} %s", legendLabels[iFile]), "p");
        }
    }
    for (int iFile = 0; iFile < 2; iFile++) {
        if (gr_a2_high[iFile]) {
            int darkColors_copy[3] = {kRed+2, kBlue+2, kGreen+4};
            gr_a2_high[iFile]->SetLineColor(darkColors_copy[iFile]); // 确保图例颜色正确
            leg4_copy->AddEntry(gr_a2_high[iFile], Form("a_{2}^{high} %s", legendLabels[iFile]), "p");
        }
    }
    for (int iFile = 0; iFile < 2; iFile++) {
        if (gr_a2_low[iFile]) {
            int lightColors_copy[3] = {kRed-7, kBlue-7, kGreen-7};
            gr_a2_low[iFile]->SetLineColor(lightColors_copy[iFile]); // 确保图例颜色正确
            leg4_copy->AddEntry(gr_a2_low[iFile], Form("a_{2}^{low} %s", legendLabels[iFile]), "p");
        }
    }
    
    leg4_copy->SetBorderSize(0);
    leg4_copy->SetFillStyle(0);
    leg4_copy->SetTextSize(0.025); // 调小字体以适应更多条目
    leg4_copy->Draw();
    c4_individual->SaveAs("subplot4_all_coefficients_vs_eta.png");
    
    // 创建只显示c2的单独图
    TCanvas *c_c2_only = new TCanvas("c_c2_only", "c2 Comparison", 800, 600);
    c_c2_only->cd();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    
    firstDraw = true;
    for (int iFile = 0; iFile < 2; iFile++) {
        if (gr_c2[iFile]) {
            gr_c2[iFile]->SetTitle("Long-range flow coefficient c_{2} vs #eta");
            gr_c2[iFile]->GetXaxis()->SetTitle("#eta");
            gr_c2[iFile]->GetYaxis()->SetTitle("c_{2}");
            gr_c2[iFile]->GetXaxis()->SetRangeUser(-2.5, 2.5);
            gr_c2[iFile]->GetYaxis()->SetRangeUser(0.0025, 0.01);
            
            if (firstDraw) {
                gr_c2[iFile]->Draw("AP");
                firstDraw = false;
            } else {
                gr_c2[iFile]->Draw("P SAME");
            }
        }
    }
    
    TLegend *leg_c2_only = new TLegend(0.2, 0.7, 0.8, 0.9);
    for (int iFile = 0; iFile < 2; iFile++) {
        if (gr_c2[iFile]) {
            leg_c2_only->AddEntry(gr_c2[iFile], legendLabels[iFile], "p");
        }
    }
    leg_c2_only->SetBorderSize(0);
    leg_c2_only->SetFillStyle(0);
    leg_c2_only->Draw();
    
    c_c2_only->SaveAs("c2_vs_eta_comparison_only.png");
    
    // 保存最终结果到输出文件
    outfile->cd();
    //c->Write("multifile_comparison");
    c_c2_only->Write("c2_comparison_only");
    c1_individual->Write("subplot1_c2");
    c2_individual->Write("subplot2_a2_high"); 
    c3_individual->Write("subplot3_a2_low");
    c4_individual->Write("subplot4_all_coefficients");
    
    outfile->Close();
    
    // 清理内存
    //delete c;
    delete c_c2_only;
    delete c1_individual;
    delete c2_individual;
    delete c3_individual;
    delete c4_individual;
    // delete leg1;
    // delete leg2; 
    // delete leg3;
    // delete leg4;
    delete leg_c2_only;
    delete leg1_copy;
    delete leg2_copy;
    delete leg3_copy;
    delete leg4_copy;
    
    std::cout << "Analysis completed successfully!" << std::endl;
}