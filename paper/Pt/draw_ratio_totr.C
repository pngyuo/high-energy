#include <iostream>
#include <vector>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TPaveText.h>

TGraphErrors* drawConnectedPoints(TH1D* hist, Color_t color, Style_t markerStyle) {
    int n = hist->GetNbinsX();
    double* x = new double[n];
    double* y = new double[n];
    double* ey = new double[n]; // y误差数组

    int count = 0; // 用于记录满足条件的点的数量
    for (int i = 1; i <= n; ++i) {
        double binCenter = hist->GetBinCenter(i);
        if (binCenter < 3) { // 只保留 x < 3 的点
            x[count] = binCenter;
            y[count] = hist->GetBinContent(i);
            ey[count] = hist->GetBinError(i);
            count++;
        }
    }

    TGraphErrors* graph = new TGraphErrors(count, x, y, 0, ey); // 使用实际的点数 count
    graph->SetLineColor(color);
    graph->SetLineWidth(2);
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerColor(color);

    delete[] x; // 释放动态分配的内存
    delete[] y;
    delete[] ey;

    return graph;  // 返回 TGraphErrors 对象
}

void addText(Double_t x, Double_t y, const char* text)
  {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);  // 将坐标转换为百分比
    latex->SetTextSize(0.04);  // 设置字体大小
    latex->Draw();
  }
  void addText1(Double_t x, Double_t y, const char* text)
  {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);  // 将坐标转换为百分比
    latex->SetTextSize(0.05);  // 设置字体大小
    latex->Draw();
  }
void addText2(Double_t x, Double_t y, const char* text)
  {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);  // 将坐标转换为百分比
    latex->SetTextSize(0.06);  // 设置字体大小
    latex->Draw();
  }
void addText3(Double_t x, Double_t y, const char* text)
  {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);  // 将坐标转换为百分比
    latex->SetTextSize(0.08);  // 设置字体大小
    latex->Draw();
  }
void addText4(Double_t x, Double_t y, const char* text)
  {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);  // 将坐标转换为百分比
    latex->SetTextSize(0.10);  // 设置字体大小
    latex->Draw();
  }

  void addText5(Double_t x, Double_t y, const char* text)
  {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);  // 将坐标转换为百分比
    latex->SetTextSize(0.062);  // 设置字体大小
    latex->Draw();
  }

  std::vector<TH1D*> getCentPtProjections(const TString& filename, const TString& histname) {
    TFile* infile = new TFile(filename, "READ");
    TH3D* hist3D = (TH3D*)infile->Get(histname);
    TH1D *chargedParticlesCount = (TH1D*)infile->Get("nTHist");
    
    double totalChargedParticles = 0; // 加权总和
    double nEvents = 0;  // 总的事件数
    for (int bin = 1; bin <= chargedParticlesCount->GetNbinsX(); bin++) {
        double binCenter = chargedParticlesCount->GetXaxis()->GetBinCenter(bin);
        double binContent = chargedParticlesCount->GetBinContent(bin);
        totalChargedParticles += binCenter * binContent;
        nEvents += binContent;
    }
    double averageNT = totalChargedParticles / nEvents;
    std::cout << "averageNT: " << averageNT << std::endl;
    
    std::vector<TH1D*> projections;
    double rtRanges[2] = {0, 2.5};  // RT 的下限
    double rtUpperRanges[2] = {0.5, 5};  // RT 的上限
  
    for (int i = 0; i < 2; i++) {
        double ntRangeLow = rtRanges[i] * averageNT;
        double ntRangeUp = rtUpperRanges[i] * averageNT;

        int binLow = hist3D->GetZaxis()->FindBin(ntRangeLow);
        int binUp = hist3D->GetZaxis()->FindBin(ntRangeUp);
        binUp -= 1;

        // 计算该RT范围内的事件数
        double totalEventsInRange = 0;
        for (int binIndex = binLow; binIndex <= binUp; ++binIndex) {
            totalEventsInRange += chargedParticlesCount->GetBinContent(binIndex);
        }
        
        // 输出事件数信息
        std::cout << "RT Range [" << rtRanges[i] << ", " << rtUpperRanges[i] 
                  << "] has " << totalEventsInRange << " events" << std::endl;

        TH1D* histProjX = hist3D->ProjectionX(Form("histProjX_%i", i + 1), 1, hist3D->GetNbinsY(), binLow, binUp);

        histProjX->Scale(1.0 / totalEventsInRange, "width");
        projections.push_back(histProjX);
    }
    return projections;
}

TH1D* getCentPtProjectionsAll(const TString& filename, const TString& histname) {
    TFile* infile = new TFile(filename, "READ");
    TH3D* hist3D = (TH3D*)infile->Get(histname);
    TH1D *chargedParticlesCount = (TH1D*)infile->Get("nTHist");
    
    double totalChargedParticles = 0;
    double nEvents = 0;
    for (int bin = 1; bin <= chargedParticlesCount->GetNbinsX(); bin++) {
        double binCenter = chargedParticlesCount->GetXaxis()->GetBinCenter(bin);
        double binContent = chargedParticlesCount->GetBinContent(bin);
        totalChargedParticles += binCenter * binContent;
        nEvents += binContent;
    }
    double averageNT = totalChargedParticles / nEvents;
    std::cout << "averageNT: " << averageNT << std::endl;

    // 计算总事件数
    double totalEvents = 0;
    for (int binIndex = 1; binIndex <= chargedParticlesCount->GetNbinsX(); ++binIndex) {
        totalEvents += chargedParticlesCount->GetBinContent(binIndex);
    }
    
    // 输出总事件数
    std::cout << "Total events: " << totalEvents << std::endl;

    TH1D* histProjX = hist3D->ProjectionX("histProjX_all", 1, hist3D->GetNbinsY());
    
    histProjX->Scale(1.0 / totalEvents, "width");
    return histProjX;
}

TH1D* createRatioHistogram(TH1D* numerator, TH1D* denominator, const char* name) {
    TH1D* ratio = (TH1D*)numerator->Clone(name);
    ratio->Divide(denominator);
    return ratio;
}

void draw_ratio_totr() {
    TCanvas *c1 = new TCanvas("c1", "c1",400, 420);
    c1->SetTicks(1, 1);
    c1->SetMargin(0.17, 0.02, 0.18, 0.1);
    gStyle->SetOptStat(kFALSE);
    
    // 通用样式设置
    std::vector<Color_t> colors = {kBlack, kRed, kBlue}; // 三种曲线颜色
    std::vector<Style_t> markerStyles = {20, 21, 22};      // 三种标记样
  
        c1->SetTicks(1,1);
        TH2D *axisFrame1 = new TH2D("axisFrame1", "", 100,-0.1, 2.99, 100, 0.02, 3);
        axisFrame1->GetYaxis()->SetTickLength(0.02);
        axisFrame1->GetXaxis()->SetTickLength(0.02);
        axisFrame1->GetXaxis()->SetLabelSize(0.070);
        axisFrame1->GetYaxis()->SetLabelSize(0.070);
        axisFrame1->Draw("axis");

        TLatex *label = new TLatex();
        label->SetTextSize(0.09);
        label->SetTextAngle(90);
        label->SetTextAlign(22);
        label->DrawLatexNDC(0.045, 0.5, "To/Tr");

        TString filename = "hist_outputallFSI_liang.root";
        TString histname1 = "hPiCh_dPhi0";
        TString histname2 = "hPiCh_dPhi1";
        // 绘制所有RT区间
        TH1D* hist_all_toward = getCentPtProjectionsAll(filename, histname1);
        TH1D* hist_all_transverse = getCentPtProjectionsAll(filename, histname2);
        TH1D* hist_all = createRatioHistogram(hist_all_toward, hist_all_transverse, "hist_all_ratio");
        TGraphErrors* graph_all = drawConnectedPoints(hist_all, colors[0], markerStyles[0]);
        graph_all->Draw("LP");
        
        std::vector<TH1D*> projections_toward = getCentPtProjections(filename, histname1); 
        std::vector<TH1D*> projections_transverse = getCentPtProjections(filename, histname2);         
        std::vector<TH1D*> ratio_hists; // 存储比率直方图
        
for (int i = 0; i < projections_toward.size(); i++) {
    TH1D* ratio_hist = createRatioHistogram(projections_toward[i], projections_transverse[i], 
                                           Form("ratio_%d", i));

    // 对第二条线（i == 0，即红色）乘以 0.1
    if (i == 0) {
        for (int bin = 1; bin <= ratio_hist->GetNbinsX(); ++bin) {
            ratio_hist->SetBinContent(bin, ratio_hist->GetBinContent(bin) * 0.1);
            ratio_hist->SetBinError(bin, ratio_hist->GetBinError(bin) * 0.1);
        }
    }

    ratio_hists.push_back(ratio_hist);
    TGraphErrors* graph_rt = drawConnectedPoints(ratio_hist, colors[i+1], markerStyles[i+1]);
    graph_rt->Draw("LP same");
}
        
        axisFrame1->Draw("axis same");

        // 创建图例
        TLegend *legend = new TLegend(0.20, 0.55, 0.70, 0.75);
        legend->SetNColumns(1);
        legend->SetFillStyle(0);
        legend->SetBorderSize(0);
        legend->SetMargin(0.20);
        legend->SetTextSize(0.056);

TLatex *text2 = new TLatex(0.59, 0.06, "p_{T} (GeV/c)");
text2->SetNDC();
text2->SetTextSize(0.08);
text2->SetTextAlign(22);
text2->Draw();

        addText5(0.25, 0.75, "AMPT");
        legend->AddEntry(graph_all, "R_{T}>0", "lp");
        legend->AddEntry(drawConnectedPoints(ratio_hists[0], colors[1], markerStyles[1]), "0<R_{T}<0.5", "lp");
        legend->AddEntry(drawConnectedPoints(ratio_hists[1], colors[2], markerStyles[2]), "2.5<R_{T}<5", "lp");
        legend->Draw();
        addText4(0.50, 0.81, "#pi^{+}+#pi^{-}");
        addText1(0.92, 0.20, "(a)");

    c1->SaveAs("pi_totrratio.png");
    c1->SaveAs("pi_totrratio.pdf");
}