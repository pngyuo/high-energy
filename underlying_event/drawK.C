#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>

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
    chargedParticlesCount->Draw();
    double totalChargedParticles = 0; // 加权总和
    double nEvents = 0;  // 总的事件数
    for (int bin = 1; bin <= 36; bin++) {
        double binCenter = chargedParticlesCount->GetXaxis()->GetBinCenter(bin); // 获取当前 bin 中心的 NT 值
        double binContent = chargedParticlesCount->GetBinContent(bin);  // 获取当前 bin 的分布数量
        totalChargedParticles += binCenter * binContent;  // 加权累加
        nEvents += binContent;  // 累加事件数
    }
    double averageNT = totalChargedParticles / nEvents;  // 计算加权平均NT
    cout<<"averageNT: "<<averageNT<<endl;
    std::vector<TH1D*> projections;

    // 计算新的范围
    double rtRanges[4] = {0, 0.5, 1.5, 2.5};  // RT 的范围
    double rtUpperRanges[4] = {0.5, 1.5, 2.5, 30 / averageNT};  // RT 上限

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
        TH1D* histProjX = hist3D->ProjectionX(Form("histProjX_%i", i+1), 1, hist3D->GetNbinsY(), binLow, binUp);

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

void drawGraphs(TVirtualPad *pad, const std::vector<TH1D*>& projections, const std::vector<Color_t>& colors, const std::vector<Style_t>& markerStyles) {
  pad->cd(); // 切换到当前子画布

  // 绘制坐标轴
  TH2D *axisFrame = new TH2D("axisFrame", "; ; ", 100, -0.1, 5.1, 100, 1e-5, 1e2);
  axisFrame->Draw("axis");  // 只绘制坐标轴，不绘制直方图

  for (size_t i = 0; i < projections.size(); ++i) {
    TH1D* hist = projections[i];
    Color_t color = colors[i];
    Style_t markerStyle = markerStyles[i];
    TGraph* graph = drawConnectedPoints(hist, color, markerStyle);

    if(i == 0) {
      graph->Draw("LP"); // L: 折线, P: 数据点
    } else {
      graph->Draw("LP same"); // 将新的图形绘制在同一画布上
    }
  }
  
  // 重新绘制坐标轴以确保它们在最上面
  axisFrame->Draw("axis same");
}

void drawK(){
  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 500); // 创建Canvas

  TPaveText *text = new TPaveText(0.4, 0.01, 0.6, 0.1, "NDC");
  text->AddText("p_{T}(GeV/c)");
  text->SetFillColor(0);
  text->SetTextAlign(22); // 文本居中对齐
  text->Draw("same");

  TLatex *label = new TLatex(0.01, 0.6, "dN/(dp_{T}dy)(GeV/c)^{-1}");
  label->SetTextSize(0.06);
  label->SetTextAngle(90);
  label->SetTextAlign(23);
  label->Draw();

  c1->Divide(3, 1);
  c1->cd(1)->SetPad(0.04, 0.1, 0.39, 1); 
  c1->cd(1)->SetRightMargin(0);
  //c1->cd(1)->SetLeftMargin(0);
  c1->cd(2)->SetPad(0.39, 0.1, 0.69, 1); 
  c1->cd(2)->SetRightMargin(0);
  c1->cd(2)->SetLeftMargin(0);
  c1->cd(3)->SetPad(0.69, 0.1, 0.99, 1);
  c1->cd(3)->SetRightMargin(0);
  c1->cd(3)->SetLeftMargin(0);

  TString filename = "hist_outputMonash.root"; 
  TString histname1 = "hKCh_dPhi0"; 
  TString histname2 = "hKCh_dPhi2"; 
  TString histname3 = "hKCh_dPhi1"; 

  // 定义颜色数组
  std::vector<Color_t> colors = {kBlack, kRed, kBlue, kMagenta};
  // 定义标记样式数组
  std::vector<Style_t> markerStyles = {20, 21, 22, 23};

  c1->cd(1);
  std::vector<TH1D*> projections1 = getCentPtProjections(filename, histname1);

  // 绘制图表
  drawGraphs(gPad,projections1, colors, markerStyles);

  
  TLegend *legend = new TLegend(0.45, 0.18, 0.85, 0.39);
  legend->SetNColumns(1);
  legend->SetFillStyle(0);

  TGraphErrors *marker1 = new TGraphErrors();
  marker1->SetMarkerStyle(20); // 黑色实心圆标记
  marker1->SetMarkerColor(kBlack);
  marker1->SetMarkerSize(1);
  marker1->SetLineWidth(2); 
  marker1->SetLineColor(kBlack);
  legend->AddEntry(marker1, " ", "lp");



  TGraphErrors *marker3 = new TGraphErrors();
  marker3->SetMarkerStyle(22); // 蓝色实心三角形
  marker3->SetMarkerColor(kBlue);
  marker3->SetMarkerSize(1);
  marker3->SetLineWidth(2); 
  marker3->SetLineColor(kBlue);
  legend->AddEntry(marker3, " ", "lp");

  TGraphErrors *marker4 = new TGraphErrors();
  marker4->SetMarkerStyle(23); // 
  marker4->SetMarkerColor(kMagenta);
  marker4->SetMarkerSize(1);
  marker4->SetLineWidth(2); 
  marker4->SetLineColor(kMagenta);
  legend->AddEntry(marker4, " ", "lp");

  legend->Draw("same");

  addText3(0.2, 0.4, "Monash");
  addText(0.2, 0.35, "0 #leq R_{T} < 0.5");
  addText(0.2, 0.3, "0.5 #leq R_{T} < 1.5");
  addText(0.2, 0.25, "1.5 #leq R_{T} < 2.5");
  addText(0.2, 0.2, "2.5 #leq R_{T} < 5");
  addText2(0.4, 0.8, "Toward");

  c1->cd(2);
  std::vector<TH1D*> projections2 = getCentPtProjections(filename, histname2);
  // 绘制图表
  drawGraphs(gPad,projections2, colors, markerStyles);

  addText4(0.33, 0.25, "#Kappa^{+}+#Kappa^{-}");
  addText2(0.4, 0.8, "Away");
  addText2(0.65, 0.8, "PYTHIA8");

  c1->cd(3);
  std::vector<TH1D*> projections3 = getCentPtProjections(filename, histname3);
  // 绘制图表
  drawGraphs(gPad,projections3, colors, markerStyles);
  addText2(0.4, 0.8, "Transverse");
  c1->cd(1);
  gPad->SetLogy();
  c1->cd(2);
  gPad->SetLogy();
  c1->cd(3);
  gPad->SetLogy();

  c1->SaveAs("INEL K_Monash_Pt.png");
  c1->SaveAs("INEL K_Monash_Pt.pdf");
}


