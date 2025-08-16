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
    cout << "averageNT: " << averageNT << endl;
    std::vector<TH1D*> projections;

    // 根据 bin 范围投影 X 轴上的直方图
    TH1D* histProjX = hist3D->ProjectionX(Form("histProjX"), 1, hist3D->GetNbinsY());

    double totalEventsInRange = 0;
    for (int binIndex = 1; binIndex <= chargedParticlesCount->GetNbinsX(); ++binIndex) {
        double eventsWithParticleCount = chargedParticlesCount->GetBinContent(binIndex);
        totalEventsInRange += eventsWithParticleCount;
    }

    histProjX->Scale(1.0 / totalEventsInRange, "width");
    histProjX->Rebin(2);  // 这里将每个2个bin合并为一个，你可以根据需要调整这个参数
    projections.push_back(histProjX);
    
    return projections;
}

std::vector<TH1D*> getCentPtProjections1(const TString& filename, const TString& histname) {
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

        // 根据 bin 范围投影 X 轴上的直方图
        TH1D* histProjX = hist3D->ProjectionX(Form("histProjX"), 1, hist3D->GetNbinsY());

        double totalEventsInRange = 0;
         for (int binIndex = 1; binIndex <= chargedParticlesCount->GetNbinsX(); ++binIndex) {
        double eventsWithParticleCount = chargedParticlesCount->GetBinContent(binIndex);
        totalEventsInRange += eventsWithParticleCount;
    }

        histProjX->Scale(2.0 / totalEventsInRange, "width");
        projections.push_back(histProjX);
        
    return projections;
}

void drawGraphs(TVirtualPad *pad, const std::vector<TH1D*>& projections1, std::vector<TH1D*> projections2, std::vector<TH1D*> projections3, std::vector<TH1D*> projections4, Color_t color, Style_t markerStyle) {
  pad->cd(); // 切换到当前子画布
  for (size_t i = 0; i < projections1.size(); ++i) {
    TH1D* hist1 = projections1[i];
    TH1D* hist2 = projections2[i];
    TH1D* hist3 = projections3[i];
    TH1D* hist4 = projections4[i];
 
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

void drawKLog(){
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600); // 创建Canvas
  TH2D *axisFrame = new TH2D("axisFrame", "; ; ", 100, -0.1, 5.1, 100, 0, 0.4);
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
label->DrawLatexNDC(0.03, 0.5, "(#Kappa^{+}+#Kappa^{-})/(#pi^{+}+#pi^{-})");// 使用NDC坐标

  TString filename[3] = {"hist_outputallFSI.root","hist_outputnoFSI.root","hist_outputnohFSI.root"}; 
  TString histname1 = "hKCh_dPhi0"; 
  TString histname2 = "hKCh_dPhi1"; 
    TString histname3 = "hPiCh_dPhi0"; 
  TString histname4 = "hPiCh_dPhi1"; 

  // 定义颜色数组
  std::vector<Color_t> colors = {kBlack, kRed, kBlue, kMagenta};
  // 定义标记样式数组
  std::vector<Style_t> markerStyles = {20, 20, 21, 22};
   gStyle->SetOptStat(kFALSE);
   for(int i=0;i<3;++i){
        Color_t color = colors[i];
    Style_t markerStyle = markerStyles[i];
  std::vector<TH1D*> projections1 = getCentPtProjections(filename[i], histname1);
  std::vector<TH1D*> projections2 = getCentPtProjections1(filename[i], histname2);
   std::vector<TH1D*> projections3 = getCentPtProjections(filename[i], histname3);
  std::vector<TH1D*> projections4 = getCentPtProjections1(filename[i], histname4);
  drawGraphs(gPad,projections1, projections2,projections3, projections4,color, markerStyle);

   }
    TLegend *legend = new TLegend(0.65, 0.18, 0.85,0.3);
    legend->SetNColumns(1);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);

    TGraphErrors *marker1 = new TGraphErrors();
    marker1->SetMarkerStyle(20);
    marker1->SetMarkerColor(kBlack);
    marker1->SetMarkerSize(1);
    marker1->SetLineWidth(2); 
    marker1->SetLineColor(kBlack);
    legend->AddEntry(marker1, "allFSI", "lp");

    TGraphErrors *marker2 = new TGraphErrors();
    marker2->SetMarkerStyle(20);
    marker2->SetMarkerColor(kRed);
    marker2->SetMarkerSize(1);
    marker2->SetLineWidth(2); 
    marker2->SetLineColor(kRed);
    legend->AddEntry(marker2, "noFSI", "lp");

    TGraphErrors *marker3 = new TGraphErrors();
    marker3->SetMarkerStyle(21);
    marker3->SetMarkerColor(kBlue);
    marker3->SetMarkerSize(1);
    marker3->SetLineWidth(2); 
    marker3->SetLineColor(kBlue);
    legend->AddEntry(marker3, "nohFSI", "lp");

    TGraphErrors *marker4 = new TGraphErrors();
    marker4->SetMarkerStyle(22);
    marker4->SetMarkerColor(kMagenta);
    marker4->SetMarkerSize(1);
    marker4->SetLineWidth(2); 
    marker4->SetLineColor(kMagenta);
    //legend->AddEntry(marker4, "2.5<R_{T}<5", "lp");
    legend->Draw(); 
  //addText2(0.4, 0.8, "Monash");
  //addText4(0.33, 0.20, "#Kappa^{+}+#Kappa^{-}");
   //c1->SetLogy(); // 设置y轴为对数刻度
  c1->SaveAs("K_Monash.png");
  //c1->SaveAs("proton_CR_Pt.pdf");
}
