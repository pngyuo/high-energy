#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TFile.h>

TGraph* drawConnectedPoints(TH1D* hist, Color_t color, Style_t markerStyle,double& averageNT) {
  int n = hist->GetNbinsX();
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> ey; // y误差

  for (int i = 1; i <= n; ++i) {
    double binContent = hist->GetBinContent(i);
    if (binContent != 0) { // 只考虑非零数据点
      double binCenter = hist->GetBinCenter(i);
      double binError = hist->GetBinError(i); // 获取bin的误差
      cout<<"binError: "<<binError<<endl;
      if ( binContent>0&&binCenter < 2) { 
        x.push_back(binCenter);
        y.push_back(binContent);
        ey.push_back(binError);
      }
    }
  }

  int nPoints = x.size();
  if (nPoints == 0) return nullptr; // 如果没有非零点，则返回nullptr

  double* xArr = new double[nPoints];
  double* yArr = new double[nPoints];
  double* eyArr = new double[nPoints];
  for (int i = 0; i < nPoints; ++i) {
    xArr[i] = x[i];
    yArr[i] = y[i];
    eyArr[i] = ey[i];
  }
   TGraphErrors* graph = new TGraphErrors(nPoints, xArr, yArr, 0, eyArr);
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

std::vector<TH1D*> getCentPtProjections(const TString& filename, const TString& histname,double& averageNT) {
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
   averageNT = totalChargedParticles / nEvents;  // 计算加权平均NT
    cout<<"averageNT: "<<averageNT<<endl;
    std::vector<TH1D*> projections;

    // 投影到Z轴，并限制X和Y轴的范围
    TH1D* histProjZ = hist3D->ProjectionZ("histProjZ", 1, hist3D->GetNbinsX(), 1, hist3D->GetNbinsY());

    // Rebin the histogram if needed, for example, rebin by a factor of 2
    int rebinFactor = 3;
    TH1D* histProjZRebin = (TH1D*)histProjZ->Rebin(rebinFactor);

    // 获取原始直方图的bins数量和X轴范围
    int nbins = histProjZRebin->GetNbinsX();
    double xmin = histProjZRebin->GetXaxis()->GetBinLowEdge(1) / averageNT;
    double xmax = histProjZRebin->GetXaxis()->GetBinUpEdge(nbins) / averageNT;

    // 创建一个新的TH1D对象，其横轴范围缩小了averageNT倍
    TH1D *scaledHistProjZ = new TH1D("scaledHistProjZ", "Scaled Projection Z", nbins, xmin, xmax);

    // 填充新的TH1D对象
    for (int i = 1; i <= nbins; ++i) {
        double binContent = histProjZRebin->GetBinContent(i);
        scaledHistProjZ->SetBinContent(i, binContent);
        scaledHistProjZ->SetBinError(i, histProjZRebin->GetBinError(i));
    }

    // 设置新的TH1D对象的X轴和Y轴标签
    scaledHistProjZ->GetXaxis()->SetTitle(histProjZRebin->GetXaxis()->GetTitle());
    scaledHistProjZ->GetYaxis()->SetTitle(histProjZRebin->GetYaxis()->GetTitle());

    projections.push_back(scaledHistProjZ);
    return projections;
}

void drawGraphs(TVirtualPad *pad, const std::vector<TH1D*>& projections1, std::vector<TH1D*> projections2, std::vector<TH1D*> projections3, std::vector<TH1D*> projections4, Color_t color, Style_t markerStyle,double& averageNT) {
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
    TGraph* graphRatio = drawConnectedPoints(ratioHist, color, markerStyle,averageNT);

    if(i == 0) {
      graphRatio->Draw("LP");
    } else {
      graphRatio->Draw("LP same");
    }
  }
}

void drawGraphs1(TVirtualPad *pad, const std::vector<TH1D*>& projections1, std::vector<TH1D*> projections2, std::vector<TH1D*> projections3, std::vector<TH1D*> projections4, Color_t color, Style_t markerStyle,double& averageNT) {
  pad->cd(); // 切换到当前子画布
  for (size_t i = 0; i < projections1.size(); ++i) {
    TH1D* hist1 = projections1[i];
    TH1D* hist2 = projections2[i];
    TH1D* hist3 = projections3[i];
    TH1D* hist4 = projections4[i];
    TH1D* ratioHist1 = (TH1D*)hist1->Clone();
    ratioHist1->Divide(hist3); // 除以hist3
    // 相减操作：hist1 - hist2
    TH1D* diffHist1 = (TH1D*)hist1->Clone();
    diffHist1->Add(hist2, -1); // 减去hist2
  // 相减操作：hist3 - hist4
    TH1D* diffHist2 = (TH1D*)hist3->Clone();
    diffHist2->Add(hist4, -1); // 减去hist2
    // 相除操作：hist3 / hist4
    TH1D* ratioHist = (TH1D*)diffHist1->Clone();
    ratioHist->Divide(diffHist2); // 除以hist4

 TH1D* ratio = (TH1D*)ratioHist->Clone();
    ratio->Divide(ratioHist1); // 除以hist4
    // 绘制相减和相除的结果
    TGraph* graphRatio = drawConnectedPoints(ratio, color, markerStyle,averageNT);

    if(i == 0) {
      graphRatio->Draw("LP");
    } else {
      graphRatio->Draw("LP same");
    }
  }
}

void drawK(){
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600); // 创建Canvas
    c1->Divide(1, 2,0,0); // 将画布分为上下两部分


  c1->cd(1); // 切换到上半部分
  gPad->SetTicks(1,1); // 开启X/Y轴的双边刻度
  gPad->SetTopMargin(0.05); // 减小底部边距
  gPad->SetRightMargin(0.05); // 减小右侧边距
  TH2D *axisFrame = new TH2D("axisFrame", "; ; ", 100, 0,2, 100, 0,0.5);
  axisFrame->GetYaxis()->SetTickLength(0.02); // 调整刻度长度
  axisFrame->GetXaxis()->SetTickLength(0.02);
  axisFrame->Draw("axis");  // 只绘制坐标轴，不绘制直方图


TLatex *label = new TLatex();
label->SetTextSize(0.06);
label->SetTextAngle(90);
label->SetTextAlign(22); // 设置文本居中对齐
label->DrawLatexNDC(0.03, 0.5, "#Kappa/#pi"); // 使用NDC坐标

  TString filenames[3] = {"hist_outputallFSI_liang.root", "hist_outputnoFSI_liang.root", "hist_outputnohFSI_liang.root"};
  TString histname1 = "hKCh_dPhi0"; 
  TString histname2 = "hKCh_dPhi1"; 
  TString histname3 = "hPiCh_dPhi0"; 
  TString histname4 = "hPiCh_dPhi1"; 

  // 定义颜色数组
  std::vector<Color_t> colors = {kBlack, kRed, kBlue, kMagenta};
  // 定义标记样式数组
  std::vector<Style_t> markerStyles = {20, 20, 21, 22};
     double averageNT;
   gStyle->SetOptStat(kFALSE);
   for(int i=0;i<3;++i){
    Color_t color = colors[i];
    Style_t markerStyle = markerStyles[i];
  std::vector<TH1D*> projections1 = getCentPtProjections(filenames[i], histname1,averageNT);
  std::vector<TH1D*> projections2 = getCentPtProjections(filenames[i], histname2,averageNT);
  std::vector<TH1D*> projections3 = getCentPtProjections(filenames[i], histname3,averageNT);
  std::vector<TH1D*> projections4 = getCentPtProjections(filenames[i], histname4,averageNT);
  drawGraphs(gPad,projections1,projections2,projections3, projections4,color, markerStyle, averageNT);
   }

    TLegend *legend = new TLegend(0.15, 0.18, 0.35,0.3);
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
    //legend->AddEntry(marker4, "nohFSI", "lp");
    legend->Draw(); 


     c1->cd(2); // 切换到下半部分
  gPad->SetTicks(1,1); // 开启X/Y轴的双边刻度
  gPad->SetRightMargin(0.05); // 减小右侧边距
   TH2D *axisFrame1 = new TH2D("axisFrame", "; ; ", 100, 0,2, 100, 0,5);
    axisFrame1->GetYaxis()->SetTickLength(0.02); // 调整刻度长度
  axisFrame1->GetXaxis()->SetTickLength(0.02);
  axisFrame1->Draw("axis");  // 只绘制坐标轴，不绘制直方图
   for(int i=0;i<3;++i){
    Color_t color = colors[i];
    Style_t markerStyle = markerStyles[i];
  std::vector<TH1D*> projections1 = getCentPtProjections(filenames[i], histname1,averageNT);
  std::vector<TH1D*> projections2 = getCentPtProjections(filenames[i], histname2,averageNT);
  std::vector<TH1D*> projections3 = getCentPtProjections(filenames[i], histname3,averageNT);
  std::vector<TH1D*> projections4 = getCentPtProjections(filenames[i], histname4,averageNT);
  drawGraphs1(gPad,projections1,projections2,projections3, projections4,color, markerStyle, averageNT);
   }
  //addText2(0.4, 0.8, "allFSI");
  //addText4(0.33, 0.20, "#Kappa^{+}+#Kappa^{-}");
   //c1->SetLogy(); // 设置y轴为对数刻度
     TPaveText *text = new TPaveText(0.425, 0, 0.575, 0.06, "NDC");
  text->AddText("R_{T}");
  text->SetFillColor(0);
  text->SetTextAlign(22); // 文本居中对齐
  text->Draw("same");
     TLatex *label2 = new TLatex();
  label2->SetTextSize(0.06);
  label2->SetTextAngle(90);
  label2->SetTextAlign(22); // 设置文本居中对齐
  label2->DrawLatexNDC(0.03, 0.5, "Ratio to Toward"); //
  // 在纵坐标等于1处增加一条虚线
  TLine *line2 = new TLine(-0.1, 1, 5.1, 1);
  line2->SetLineStyle(2); // 设置虚线样式
  line2->SetLineColor(kBlack);
  line2->Draw();
  c1->SaveAs("K__allFSI_ Pt.png");
  //c1->SaveAs("proton_CR_Pt.pdf");
}
