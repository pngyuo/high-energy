#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TFile.h>

TGraphErrors* drawConnectedPoints(TH1D* hist, Color_t color, Style_t markerStyle) {
  int n = hist->GetNbinsX();
  double* x = new double[n];
  double* y = new double[n];
  double* ey = new double[n]; // 存储误差值

  for (int i = 1; i <= n; ++i) {
    x[i-1] = hist->GetBinCenter(i);
    y[i-1] = hist->GetBinContent(i);
    ey[i-1] = hist->GetBinError(i); // 获取每个箱子的误差
    cout<<"ey[i-1]: "<<ey[i-1]<<endl;
  }

  TGraphErrors* graph = new TGraphErrors(n, x, y, 0,ey);
  graph->SetLineColor(color);
  graph->SetLineWidth(2);
  graph->SetMarkerStyle(markerStyle);
  graph->SetMarkerColor(color);
  return graph;  // 返回 TGraphErrors 对象
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

  double rtRanges[4] = {0, 0.5, 2};  // RT 的范围
    double rtUpperRanges[4] = {0.5, 2,5};  // RT 上限
  
   for (int i = 0; i < 4; i++) {
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

void drawGraphs(TVirtualPad *pad, const std::vector<TH1D*>& projections1, std::vector<TH1D*> projections2,const std::vector<Color_t>& colors, const std::vector<Style_t>& markerStyles) {
  pad->cd(); // 切换到当前子画布


  for (size_t i = 0; i < projections1.size(); ++i) {
    TH1D* hist1 = projections1[i];
    TH1D* hist2 = projections2[i];
    Color_t color = colors[i];
    Style_t markerStyle = markerStyles[i];
      TH1D* diffHist = (TH1D*)hist1->Clone();
  diffHist->SetName("diffHist");
  for (int i = 1; i <= diffHist->GetNbinsX(); ++i) {
    double content1 = hist1->GetBinContent(i);
    double content2 = hist2->GetBinContent(i);
    diffHist->SetBinContent(i, content1 - content2);
  }
  TGraph* graph = drawConnectedPoints(diffHist, color, markerStyle);
    if(i == 0) {
      graph->Draw("LP"); // L: 折线, P: 数据点
    } else {
      graph->Draw("LP same"); // 将新的图形绘制在同一画布上
    }
  }

}

void drawPiLog(){
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600); // 创建Canvas
  TH2D *axisFrame = new TH2D("axisFrame", "; ; ", 100, -0.1, 5.1, 100, 1e-2, 1e1);
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

  TString filename = "hist_outputnohFSI.root"; 
  TString histname1 = "hPiCh_dPhi0"; 
  TString histname2 = "hPiCh_dPhi1"; 

  // 定义颜色数组
  std::vector<Color_t> colors = {kBlack, kRed, kBlue, kMagenta};
  // 定义标记样式数组
  std::vector<Style_t> markerStyles = {20, 20, 21, 22};
   gStyle->SetOptStat(kFALSE);
  std::vector<TH1D*> projections1 = getCentPtProjections(filename, histname1);
  std::vector<TH1D*> projections2 = getCentPtProjections(filename, histname2);
  drawGraphs(gPad,projections1, projections2,colors, markerStyles);
 
    c1->cd(1);
    TLegend *legend = new TLegend(0.65, 0.70, 0.90, 0.81);
    legend->SetNColumns(1);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);

    TGraphErrors *marker1 = new TGraphErrors();
    marker1->SetMarkerStyle(20);
    marker1->SetMarkerColor(kBlack);
    marker1->SetMarkerSize(1);
    marker1->SetLineWidth(2); 
    marker1->SetLineColor(kBlack);
    legend->AddEntry(marker1, "0<R_{T}<0.5", "lp");

    TGraphErrors *marker2 = new TGraphErrors();
    marker2->SetMarkerStyle(20);
    marker2->SetMarkerColor(kRed);
    marker2->SetMarkerSize(1);
    marker2->SetLineWidth(2); 
    marker2->SetLineColor(kRed);
    legend->AddEntry(marker2, "0.5<R_{T}<2", "lp");

    TGraphErrors *marker3 = new TGraphErrors();
    marker3->SetMarkerStyle(21);
    marker3->SetMarkerColor(kBlue);
    marker3->SetMarkerSize(1);
    marker3->SetLineWidth(2); 
    marker3->SetLineColor(kBlue);
    legend->AddEntry(marker3, "2<R_{T}<5", "lp");

    TGraphErrors *marker4 = new TGraphErrors();
    marker4->SetMarkerStyle(22);
    marker4->SetMarkerColor(kMagenta);
    marker4->SetMarkerSize(1);
    marker4->SetLineWidth(2); 
    marker4->SetLineColor(kMagenta);
    //legend->AddEntry(marker4, "2<R_{T}<5", "lp");
    legend->Draw(); 
  addText2(0.4, 0.8, "nohFSI");
  addText4(0.33, 0.20, "#pi^{+}+#pi^{-}");
   c1->SetLogy(); // 设置y轴为对数刻度
  c1->SaveAs("pi_nohFSI_ Pt.png");
  //c1->SaveAs("proton_CR_Pt.pdf");
}
