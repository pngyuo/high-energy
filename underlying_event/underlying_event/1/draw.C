#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TH3D.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <vector>
#include <TMath.h>

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

void addText(Double_t x, Double_t y, const char* text) {
  TLatex* latex = new TLatex(x, y, text);
  latex->SetNDC(kTRUE);  // 将坐标转换为百分比
  latex->SetTextSize(0.04);  // 设置字体大小
  latex->Draw();
}

std::vector<TH1D*> getCentPtProjections(const TString& filename, const TString& histname) {
    TFile* infile = new TFile(filename, "READ");
    TH3D* hist3D = (TH3D*)infile->Get(histname);

    TH1D *chargedParticlesCount = (TH1D*)infile->Get("chargedParticlesHist");
    chargedParticlesCount->Draw();
    double totalChargedParticles = 0; // 加权总和
    double nEvents = 0;  // 总的事件数
    for (int bin = 1; bin <= chargedParticlesCount->GetNbinsX(); bin++) {
        double binCenter = chargedParticlesCount->GetXaxis()->GetBinCenter(bin); // 获取当前 bin 中心的 NT 值
        double binContent = chargedParticlesCount->GetBinContent(bin);  // 获取当前 bin 的分布数量
        totalChargedParticles += binCenter * binContent;  // 加权累加
        nEvents += binContent;  // 累加事件数
    }
    double averageNT = totalChargedParticles / nEvents;  // 计算加权平均NT
    std::cout<<"averageNT: "<<averageNT<<std::endl;
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

void drawGraphs(TVirtualPad *pad, TH1D* hist1, TH1D* hist2, Color_t color, Style_t markerStyle) {
  pad->cd(); // 切换到当前子画布

  // 绘制坐标轴
  TH2D *axisFrame = new TH2D("axisFrame", "; ; ", 100, -0.1, 5.1, 100, 1e-5, 1e2);
  axisFrame->Draw("axis");  // 只绘制坐标轴，不绘制直方图

  // 计算两个直方图的差值
  TH1D* diffHist = (TH1D*)hist1->Clone();
  diffHist->SetName("diffHist");
  for (int i = 1; i <= diffHist->GetNbinsX(); ++i) {
    double content1 = hist1->GetBinContent(i);
    double content2 = hist2->GetBinContent(i);
    diffHist->SetBinContent(i, content1 - content2);
    diffHist->SetBinError(i, TMath::Sqrt(TMath::Power(hist1->GetBinError(i), 2) + TMath::Power(hist2->GetBinError(i), 2)));
  }

  TGraph* graph = drawConnectedPoints(diffHist, color, markerStyle);
  graph->Draw("LP"); // L: 折线, P: 数据点

  // 重新绘制坐标轴以确保它们在最上面
  axisFrame->Draw("axis same");
    TLatex *label = new TLatex(0.1, 0.9, "#Delta p_{T} (Toward - Transverse)");
  label->SetTextSize(0.06);
  label->Draw();
  TPaveText *text = new TPaveText(0.4, 0.01, 0.6, 0.1, "NDC");
  text->AddText("p_{T}(GeV/c)");
  text->SetFillColor(0);
  text->SetTextAlign(22); // 文本居中对齐
  text->Draw("same");
    TLatex *label1 = new TLatex(0.01, 0.6, "dN/(dp_{T}dy)(GeV/c)^{-1}");
  label1->SetTextSize(0.06);
  label1->SetTextAngle(90);
  label1->SetTextAlign(23);
  label1->Draw();
   TFile *outputFile = new TFile("pi-jian-cr.root", "RECREATE");
  graph->Write(); // 将Canvas写入文件
  outputFile->Close(); // 关闭文件
}

void draw(){
  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 500); // 创建Canvas

  TString filename = "hist_outputCR.root"; 
  TString histnameToward = "hPiCh_dPhi0"; 
  TString histnameTransverse = "hPiCh_dPhi1"; 

  // 定义颜色和标记样式
  Color_t color = kBlack;
  Style_t markerStyle = 20;

  // 获取Toward和Transverse区域的投影
  std::vector<TH1D*> projectionsToward = getCentPtProjections(filename, histnameToward);
  std::vector<TH1D*> projectionsTransverse = getCentPtProjections(filename, histnameTransverse);

  // 绘制Toward区域减去Transverse区域的pt谱
  c1->cd();
  drawGraphs(gPad, projectionsToward[0], projectionsTransverse[0], color, markerStyle);

  // 添加标题和标签

  c1->SetLogy();
  c1->SaveAs("INEL proton_Monash_Pt_Diff.png");
  c1->SaveAs("INEL proton_Monash_Pt_Diff.pdf");
}