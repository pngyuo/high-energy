#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TMath.h>

TGraphErrors* drawConnectedPoints(TH1D* hist, Color_t color, Style_t markerStyle) {
  int n = hist->GetNbinsX();
  double* x = new double[n];
  double* y = new double[n];
  double* ey = new double[n]; // 误差数组

  for (int i = 1; i <= n; ++i) {
    x[i-1] = hist->GetBinCenter(i);
    y[i-1] = hist->GetBinContent(i);
    ey[i-1] = hist->GetBinError(i); // 获取每个bin的误差
  }
  TGraphErrors* graph = new TGraphErrors(n, x, y, 0,ey);
  graph->SetLineColor(color);
  graph->SetLineWidth(2);
  graph->SetMarkerStyle(markerStyle);
  graph->SetMarkerColor(color);
  return graph;  // 返回 TGraphErrors 对象
}

void addText(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);  // 将坐标转换为百分比
    latex->SetTextSize(0.05);  // 设置字体大小
    latex->Draw();
}
void addText2(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);  // 将坐标转换为百分比
    latex->SetTextSize(0.08);  // 设置字体大小
    latex->Draw();
}
void addText3(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);  // 将坐标转换为百分比
    latex->SetTextSize(0.06);  // 设置字体大小
    latex->Draw();
}
void addText4(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);  // 将坐标转换为百分比
    latex->SetTextSize(0.15);  // 设置字体大小
    latex->Draw();
}

std::vector<TH1D*> getCentPtProjections(const TString& filename, const TString& histname) {
    TFile* infile = new TFile(filename, "READ");
    TH3D* hist3D = (TH3D*)infile->Get(histname);

    TH1D *chargedParticlesCount = (TH1D*)infile->Get("chargedParticlesHist");
    
    double totalChargedParticles = 0; // 加权总和
    double nEvents = 0;  // 总的事件数
    for (int bin = 1; bin <= chargedParticlesCount->GetNbinsX(); bin++) {
        double binCenter = chargedParticlesCount->GetXaxis()->GetBinCenter(bin); // 获取当前 bin 中心的 NT 值
        double binContent = chargedParticlesCount->GetBinContent(bin);  // 获取当前 bin 的分布数量
        totalChargedParticles += binCenter * binContent;  // 加权累加
        nEvents += binContent;  // 累加事件数
    }
    double averageNT = totalChargedParticles / nEvents;  // 计算加权平均NT
    std::cout << "averageNT: " << averageNT << std::endl;
    std::vector<TH1D*> projections;
   double rtRanges[2] = {0,1.5};  // RT 的范围
    double rtUpperRanges[2] = {1,5};  // RT 上限
  
   for (int i = 0; i < 2; i++) {
        // 将 RT 转换为 NT 的范围
        double ntRangeLow = rtRanges[i] * averageNT;  // RT 转为 NT 的下限
        double ntRangeUp = rtUpperRanges[i] * averageNT;  // RT 转为 NT 的上限

        // 打印调试信息
        cout << "ntRangeLow: " << ntRangeLow << endl;
        cout << "ntRangeUp: " << ntRangeUp << endl;

        // 获取 Z 轴上对应的 bin 范围
        int binLow = hist3D->GetZaxis()->FindBin(ntRangeLow);  // NT 下限对应的 bin
        int binUp = hist3D->GetZaxis()->FindBin(ntRangeUp);  // NT 上限对应的 bin
        binUp -=1;
        // 打印 bin 范围调试信息
        cout << "RT Range[" << rtRanges[i] << ", " << rtUpperRanges[i] << "] mapped to NT Bins[" 
             << binLow << ", " << binUp << "]" << endl;

        // 根据 bin 范围投影 X 轴上的直方图
        TH1D* histProjX = hist3D->ProjectionX(Form("histProjX_%i", i + 1), 1, hist3D->GetNbinsY(), binLow, binUp);

        double totalEventsInRange = 0;
        for (int binIndex = binLow; binIndex <= binUp; ++binIndex) {
            double eventsWithParticleCount = chargedParticlesCount->GetBinContent(binIndex);
            totalEventsInRange += eventsWithParticleCount;
        }

        histProjX->Scale(1.0 / totalEventsInRange, "width");
        projections.push_back(histProjX);
        }
    return projections;
}

void drawGraphs(const std::vector<TH1D*>& projections1,const std::vector<TH1D*>& projections2, const std::vector<Color_t>& colors, const std::vector<Style_t>& markerStyles) {
  for (size_t i = 0; i < projections1.size(); ++i) {
    TH1D* hist1 = projections1[i];
    TH1D* hist2 = projections2[i];
    Color_t color = colors[i];
    Style_t markerStyle = markerStyles[i];

     TH1D* ratioHist = (TH1D*)hist1->Clone();
    ratioHist->Divide(hist2); // 除以hist4

    // 绘制相减和相除的结果
    TGraph* graphRatio = drawConnectedPoints(ratioHist, color, markerStyle);

    if(i == 0) {
      graphRatio->Draw("LP");
    } else {
      graphRatio->Draw("LP same");
    }
  }
}

void drawK() {
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600); // 创建Canvas
  c1->SetTicks(1,1); // 关键修改：开启X/Y轴的双边刻度
   
  TH2D *axisFrame = new TH2D("axisFrame", "; ; ", 100, -0.1, 5.1, 100, 0, 0.6);
  axisFrame->GetYaxis()->SetTickLength(0.02); // 调整刻度长度
  axisFrame->GetXaxis()->SetTickLength(0.02);
  axisFrame->Draw("axis");
  TPaveText *text = new TPaveText(0.425, 0, 0.575, 0.06, "NDC");
  text->AddText("p_{T}(GeV/c)");
  text->SetFillColor(0);
  text->SetTextAlign(22); // 文本居中对齐
  text->Draw("same");

  TLatex *label = new TLatex();
  label->SetTextSize(0.06);
  label->SetTextAngle(90);
  label->SetTextAlign(22); // 设置文本居中对齐
  label->DrawLatexNDC(0.03, 0.5, "p/#pi"); // 使用NDC坐标

  TString filenames= "hist_outputallFSI.root";
  TString histname1 = "hProton_dPhi0";
  TString histname2 = "hPiCh_dPhi0";

  std::vector<Color_t> colors = {kBlack,kRed, kBlue, kMagenta};
  std::vector<Style_t> markerStyles = {20, 21, 22, 23};
  gStyle->SetOptStat(kFALSE);

     std::vector<TH1D*> projection1 = getCentPtProjections(filenames, histname1); // 获取投影
     std::vector<TH1D*> projection2 = getCentPtProjections(filenames, histname2); // 获取投影
    drawGraphs(projection1, projection2,colors, markerStyles);

  TLegend *legend = new TLegend(0.15, 0.60, 0.32,0.71);
    legend->SetNColumns(1);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);

  TGraphErrors *marker1 = new TGraphErrors();
  marker1->SetMarkerStyle(20);
  marker1->SetMarkerColor(kBlack);
  marker1->SetMarkerSize(1);
  marker1->SetLineWidth(2);
  marker1->SetLineColor(kBlack);
  legend->AddEntry(marker1, "0<R_{T}<1", "lp");

  TGraphErrors *marker2 = new TGraphErrors();
  marker2->SetMarkerStyle(21);
  marker2->SetMarkerColor(kRed);
  marker2->SetMarkerSize(1);
  marker2->SetLineWidth(2);
  marker2->SetLineColor(kRed);
  legend->AddEntry(marker2, "1.5<R_{T}<5", "lp");

  TGraphErrors *marker3 = new TGraphErrors();
  marker3->SetMarkerStyle(21);
  marker3->SetMarkerColor(kBlue);
  marker3->SetMarkerSize(1);
  marker3->SetLineWidth(2);
  marker3->SetLineColor(kBlue);
  //legend->AddEntry(marker3, "1.5<R_{T}<5", "lp");

  TGraphErrors *marker4 = new TGraphErrors();
  marker4->SetMarkerStyle(22);
  marker4->SetMarkerColor(kMagenta);
  marker4->SetMarkerSize(1);
  marker4->SetLineWidth(2);
  marker4->SetLineColor(kMagenta);
  //legend->AddEntry(marker4, "nohFSI", "lp");
  legend->Draw();

  addText(0.15, 0.73, "allFSI");
  addText2(0.4, 0.8, "Toward");

  c1->SaveAs("P_PI_To.png");
}