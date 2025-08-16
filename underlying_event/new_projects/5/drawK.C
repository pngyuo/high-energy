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
  int validPoints = 0; // 用于记录有效点的数量

  for (int i = 1; i <= n; ++i) {
    double binCenter = hist->GetBinCenter(i);
    if (binCenter > 0.25) { // 只保留横轴大于0.25的点
      x[validPoints] = binCenter;
      y[validPoints] = hist->GetBinContent(i);
      ey[validPoints] = hist->GetBinError(i); // 获取每个箱子的误差
      validPoints++;
    }
  }

  TGraphErrors* graph = new TGraphErrors(validPoints, x, y, 0, ey); // 只使用有效点
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

    TH1D *chargedParticlesCount = (TH1D*)infile->Get("nTHist");

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

double rtRanges[2] = {0,0.5};  // RT 的范围
    double rtUpperRanges[2] = {0.5,1.5};  // RT 上限

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

std::vector<TH1D*> getCentPtProjections1(const TString& filename, const TString& histname) {
    TFile* infile = new TFile(filename, "READ");
    TH3D* hist3D = (TH3D*)infile->Get(histname);

    TH1D *chargedParticlesCount = (TH1D*)infile->Get("nTHist");

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

double rtRanges[2] = {0,0.5};  // RT 的范围
    double rtUpperRanges[2] = {0.5,1.5};  // RT 上限
  
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

void drawGraphs(TVirtualPad *pad, const std::vector<TH1D*>& projections1, std::vector<TH1D*> projections2, std::vector<TH1D*> projections3, std::vector<TH1D*> projections4, const std::vector<Color_t>& colors, const std::vector<Style_t>& markerStyles) {
  pad->cd(); // 切换到当前子画布
  for (size_t i = 0; i < projections1.size(); ++i) {
    TH1D* hist1 = projections1[i];
    TH1D* hist2 = projections2[i];
    TH1D* hist3 = projections3[i];
    TH1D* hist4 = projections4[i];
     Color_t color = colors[i];
    Style_t markerStyle = markerStyles[i];
    // 相减操作：hist1 - hist2
    TH1D* diffHist1 = (TH1D*)hist1->Clone();
    diffHist1->Add(hist2, -1); // 减去hist2
  // 相减操作：hist3 - hist4
    TH1D* diffHist2 = (TH1D*)hist3->Clone();
    diffHist2->Add(hist4, -1); // 减去hist2
    // 相除操作：hist3 / hist4
    TH1D* ratioHist = (TH1D*)diffHist1->Clone();
    ratioHist->Divide(diffHist2); // 除以hist4

    // 绘制相减和相除的结果
    TGraph* graphRatio = drawConnectedPoints(ratioHist, color, markerStyle);

    if(i == 0) {
      graphRatio->Draw("LP");
    } else {
      graphRatio->Draw("LP same");
    }
  }
}

void drawGraphs1(TVirtualPad *pad, const std::vector<TH1D*>& projections1,  std::vector<TH1D*> projections3, const std::vector<Color_t>& colors1, const std::vector<Style_t>& markerStyles1) {
  pad->cd(); // 切换到当前子画布
  for (size_t i = 0; i < projections1.size(); ++i) {
    TH1D* hist1 = projections1[i];
    TH1D* hist3 = projections3[i];
     Color_t color = colors1[i];
    Style_t markerStyle = markerStyles1[i];
 
    TH1D* ratioHist = (TH1D*)hist1->Clone();
    ratioHist->Divide(hist3); // 除以hist3

    // 绘制相减和相除的结果
    TGraph* Ratio = drawConnectedPoints(ratioHist, color, markerStyle);

    if(i == 0) {
      Ratio->Draw("LP");
    } else {
      Ratio->Draw("LP same");
    }
  }
}

void drawK(){
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600); // 创建Canvas
  c1->SetTicks(1,1); // 关键修改：开启X/Y轴的双边刻度
  TH2D *axisFrame = new TH2D("axisFrame", "; ; ", 100, -0.1, 5.1, 100, 0, 0.5);
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
label->DrawLatexNDC(0.03, 0.5, "#kappa/#pi");// 使用NDC坐标

  TString filename[1] = {"hist_outputallFSI_liang.root"}; 
  TString histname1 = "hKCh_dPhi0"; 
  TString histname2 = "hKCh_dPhi1"; 
  TString histname3 = "hPiCh_dPhi0"; 
  TString histname4 = "hPiCh_dPhi1"; 

  // 定义颜色数组
  std::vector<Color_t> colors = {kBlack, kRed, kBlue, kMagenta};
  // 定义标记样式数组
  std::vector<Style_t> markerStyles = {24, 24, 21, 22};
    // 定义颜色数组
  std::vector<Color_t> colors1 = {kBlack, kRed, kBlue, kMagenta};
  // 定义标记样式数组
  std::vector<Style_t> markerStyles1 = {20, 20, 21, 22};
   gStyle->SetOptStat(kFALSE);

  std::vector<TH1D*> projections1 = getCentPtProjections(filename[0], histname1);
  std::vector<TH1D*> projections2 = getCentPtProjections1(filename[0], histname2);
   std::vector<TH1D*> projections3 = getCentPtProjections(filename[0], histname3);
  std::vector<TH1D*> projections4 = getCentPtProjections1(filename[0], histname4);
  drawGraphs(gPad,projections1, projections2,projections3, projections4,colors, markerStyles);
drawGraphs1(gPad,projections1, projections3,colors1, markerStyles1);

  TLegend *legend = new TLegend(0.65, 0.70, 0.82,0.81);
    legend->SetNColumns(1);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);

    TGraphErrors *marker1 = new TGraphErrors();
    marker1->SetMarkerStyle(20);
    marker1->SetMarkerColor(kBlack);
    marker1->SetMarkerSize(1);
    marker1->SetLineWidth(0); 
    marker1->SetLineColor(kBlack);
    legend->AddEntry(marker1, "Toward", "lp");

    TGraphErrors *marker2 = new TGraphErrors();
    marker2->SetMarkerStyle(24);
    marker2->SetMarkerColor(kBlack);
    marker2->SetMarkerSize(1);
    marker2->SetLineWidth(0); 
    marker2->SetLineColor(kRed);
    legend->AddEntry(marker2, "InJet", "lp");

    TGraphErrors *marker3 = new TGraphErrors();
    marker3->SetMarkerStyle(20);
    marker3->SetMarkerColor(kBlack);
    marker3->SetMarkerSize(1);
    marker3->SetLineWidth(0); 
    marker3->SetLineColor(kBlue);
    //legend->AddEntry(marker3, "nohFSI", "lp");

    TGraphErrors *marker4 = new TGraphErrors();
    marker4->SetMarkerStyle(20);
    marker4->SetMarkerColor(kRed);
    marker4->SetMarkerSize(1);
    marker4->SetLineWidth(0); 
    marker4->SetLineColor(kMagenta);
    //legend->AddEntry(marker4, "2.5<R_{T}<5", "lp");
    legend->Draw(); 

  TLegend *legend1 = new TLegend(0.15, 0.70, 0.32,0.81);
    legend1->SetNColumns(1);
    legend1->SetFillStyle(0);
    legend1->SetBorderSize(0);
   legend1->AddEntry(marker3, "0<R_{T}<0.5", "lp");
   legend1->AddEntry(marker4, "0.5<R_{T}<1.5", "lp");
      legend1->Draw(); 
  //addText2(0.4, 0.8, "Monash");
  //addText4(0.33, 0.20, "#Kappa^{+}+#Kappa^{-}");
   //c1->SetLogy(); // 设置y轴为对数刻度
  c1->SaveAs("2K.png");
  //c1->SaveAs("proton_CR_Pt.pdf");
}
