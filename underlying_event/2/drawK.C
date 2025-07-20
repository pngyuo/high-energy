#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TFile.h>
TGraph* drawConnectedPoints(TH1D* hist, Color_t color, Style_t markerStyle) {
  int n = hist->GetNbinsX();
  double* x = new double[n];
  double* y = new double[n];
  for (int i = 1; i <= n; ++i) {
    x[i-1] = hist->GetBinCenter(i);
    y[i-1] = hist->GetBinContent(i);
  }
  TGraph* graph = new TGraph(n, x, y);
  graph->SetLineColor(color);
  graph->SetLineWidth(2);
  graph->SetMarkerStyle(markerStyle);
  graph->SetMarkerColor(color);
  return graph;  // 返回 TGraph 对象
}

void addText(Double_t x, Double_t y, const char* text)
  {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);  // 将坐标转换为百分比
    latex->SetTextSize(0.04);  // 设置字体大小
    latex->Draw();
  }
void addText2(Double_t x, Double_t y, const char* text)
  {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);  // 将坐标转换为百分比
    latex->SetTextSize(0.08);  // 设置字体大小
    latex->Draw();
  }
void addText3(Double_t x, Double_t y, const char* text)
  {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);  // 将坐标转换为百分比
    latex->SetTextSize(0.06);  // 设置字体大小
    latex->Draw();
  }
void addText4(Double_t x, Double_t y, const char* text)
  {
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
    cout<<"averageNT: "<<averageNT<<endl;
    std::vector<TH1D*> projections;

    // 根据整个Z轴投影X轴上的直方图
    TH1D* histProjX = hist3D->ProjectionX("histProjX", 1, hist3D->GetNbinsY());

    double totalEventsInRange = 0;
    for (int binIndex = 1; binIndex <= chargedParticlesCount->GetNbinsX(); ++binIndex) {
        double eventsWithParticleCount = chargedParticlesCount->GetBinContent(binIndex);
        totalEventsInRange += eventsWithParticleCount;
    }
    
    histProjX->Scale(1.0 / totalEventsInRange, "width");
    projections.push_back(histProjX);
    return projections;
}

void drawGraphs(const std::vector<TH1D*>& projections, Color_t color, Style_t markerStyle) {


  for (size_t i = 0; i < projections.size(); ++i) {
    TH1D* hist = projections[i];
    TGraph* graph = drawConnectedPoints(hist, color, markerStyle);
     if(i == 0) {
      graph->Draw("LP"); // L: 折线, P: 数据点
    } else {
      graph->Draw("LP same"); // 将新的图形绘制在同一画布上
    }
}
   
}

void drawK(){
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600); // 创建Canvas
  TH2D *axisFrame = new TH2D("axisFrame", "; ; ", 100, -0.1, 5.1, 100, 1e-3, 1);
  axisFrame->Draw("axis");  // 只绘制坐标轴，不绘制直方图
  TPaveText *text = new TPaveText(0.425, 0, 0.575, 0.06, "NDC");
  text->AddText("p_{T}(GeV/c)");
  text->SetFillColor(0);
  text->SetTextAlign(22); // 文本居中对齐
  text->Draw("same");

TLatex *label = new TLatex();
label->SetTextSize(0.06);
label->SetTextAngle(90);
label->SetTextAlign(22); // 设置文本居中对齐
label->DrawLatexNDC(0.03, 0.5, "dN/(dp_{T}dy)(GeV/c)^{-1}"); // 使用NDC坐标


  TString filenames="hist_outputCR.root";
  TString histname1 = "hKshort_dPhi0"; 
  TString histname2 = "hKshort_dPhi1"; 

  std::vector<Color_t> colors = {kBlack, kRed, kBlue, kMagenta};
  // 定义标记样式数组
  std::vector<Style_t> markerStyles = {20, 21, 22, 23};
    gStyle->SetOptStat(kFALSE);


        std::vector<TH1D*> projection1 = getCentPtProjections(filenames, histname1); // 获取投影
        drawGraphs(projection1, colors[0], markerStyles[0]);
         std::vector<TH1D*> projection2 = getCentPtProjections(filenames, histname2); // 获取投影
        drawGraphs(projection2, colors[1], markerStyles[1]);
        

    TLegend *legend = new TLegend(0.65, 0.18, 0.85,0.3);
    legend->SetNColumns(1);
    legend->SetFillStyle(0);
    legend->SetBorderSize(1);

    TGraphErrors *marker1 = new TGraphErrors();
    marker1->SetMarkerStyle(20);
    marker1->SetMarkerColor(kBlack);
    marker1->SetMarkerSize(1);
    marker1->SetLineWidth(2); 
    marker1->SetLineColor(kBlack);
    legend->AddEntry(marker1, "Toward", "lp");

    TGraphErrors *marker2 = new TGraphErrors();
    marker2->SetMarkerStyle(21);
    marker2->SetMarkerColor(kRed);
    marker2->SetMarkerSize(1);
    marker2->SetLineWidth(2); 
    marker2->SetLineColor(kRed);
    legend->AddEntry(marker2, "Transverse", "lp");

    TGraphErrors *marker3 = new TGraphErrors();
    marker3->SetMarkerStyle(21);
    marker3->SetMarkerColor(kBlue);
    marker3->SetMarkerSize(1);
    marker3->SetLineWidth(2); 
    marker3->SetLineColor(kBlue);


    TGraphErrors *marker4 = new TGraphErrors();
    marker4->SetMarkerStyle(22);
    marker4->SetMarkerColor(kMagenta);
    marker4->SetMarkerSize(1);
    marker4->SetLineWidth(2); 
    marker4->SetLineColor(kMagenta);
    legend->Draw(); 
  addText2(0.4, 0.8, "CR");
  addText4(0.33, 0.25, "K_{S}");
  c1->SetLogy(); // 设置y轴为对数刻度
  c1->SaveAs("Ks_CRPt.png");

}