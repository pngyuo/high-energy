#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TFile.h>

TGraphErrors* drawConnectedPoints(TH1D* hist, Color_t color, Style_t markerStyle) {
    int n = hist->GetNbinsX();
    double* x = new double[n];
    double* y = new double[n];
    double* ey = new double[n]; // y误差数组

    int count = 0; // 用于记录满足条件的点的数量
    for (int i = 1; i <= n; ++i) {
        double binCenter = hist->GetBinCenter(i);
        if (binCenter < 3&&binCenter >0.25) { // 只保留 x < 3 的点
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

void drawGraphs1(TVirtualPad *pad, const std::vector<TH1D*>& projections1, std::vector<TH1D*> projections2, std::vector<TH1D*> projections3, std::vector<TH1D*> projections4, const std::vector<Color_t>& colors, const std::vector<Style_t>& markerStyles) {
  pad->cd(); // 切换到当前子画布
  for (size_t i = 0; i < projections1.size(); ++i) {
    TH1D* hist1 = projections1[i];
    TH1D* hist2 = projections2[i];
    TH1D* hist3 = projections3[i];
    TH1D* hist4 = projections4[i];
     Color_t color = colors[i];
    Style_t markerStyle = markerStyles[i];
    TH1D* Hist1 = (TH1D*)hist1->Clone();
    Hist1->Divide(hist3); // 除以hist3求toward

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
    ratio->Divide(Hist1); // 除以hist4

    // 绘制相减和相除的结果
    TGraph* graphRatio = drawConnectedPoints(ratio, color, markerStyle);

    if(i == 0) {
      graphRatio->Draw("LP");
    } else {
      graphRatio->Draw("LP same");
    }
  }
}

void drawKLog(){
      TCanvas *c1 = new TCanvas("c1", "c1", 600, 800); // 调整画布尺寸为1200x800
    c1->SetTicks(1, 1); // 开启X/Y轴的双边刻度
    c1->SetMargin(0.17, 0.17, 0.18, 0.1); // 调整边距
  c1->Divide(1, 2,0,0); // 将画布分为上下两部分


  c1->cd(1); // 切换到上半部分
  gPad->SetTicks(1,1); // 开启X/Y轴的双边刻度
  gPad->SetTopMargin(0.05); // 减小底部边距
  gPad->SetRightMargin(0.05); // 减小右侧边距
  TH2D *axisFrame = new TH2D("axisFrame", "; ; ", 100, 0, 3, 100, 0.02, 0.3);
  axisFrame->GetYaxis()->SetTickLength(0.02); // 调整刻度长度
  axisFrame->GetXaxis()->SetTickLength(0.02);
  axisFrame->GetXaxis()->SetLabelSize(0.06); // 将X轴标签字体大小调整为0.045
  axisFrame->GetYaxis()->SetLabelSize(0.06); // 将Y轴标签字体大小调整为0.045
  axisFrame->Draw("axis");

TLatex *label = new TLatex();
label->SetTextSize(0.08);
label->SetTextAngle(90);
label->SetTextAlign(22); // 设置文本居中对齐
label->DrawLatexNDC(0.05, 0.5, "#Kappa/#pi");// 使用NDC坐标
addText2(0.51, 0.85, "InJet");
  TString filename[1] = {"hist_outputallFSI_liang.root"}; 
  TString histname1 = "hKCh_dPhi0"; 
  TString histname2 = "hKCh_dPhi1"; 
  TString histname3 = "hPiCh_dPhi0"; 
  TString histname4 = "hPiCh_dPhi1"; 

  // 定义颜色数组
  std::vector<Color_t> colors = {kBlack, kRed, kBlue, kMagenta};
  // 定义标记样式数组
  std::vector<Style_t> markerStyles = {20, 20, 21, 22};
    // 定义颜色数组
   gStyle->SetOptStat(kFALSE);

  std::vector<TH1D*> projections1 = getCentPtProjections(filename[0], histname1);
  std::vector<TH1D*> projections2 = getCentPtProjections(filename[0], histname2);
  std::vector<TH1D*> projections3 = getCentPtProjections(filename[0], histname3);
  std::vector<TH1D*> projections4 = getCentPtProjections(filename[0], histname4);
  drawGraphs(gPad,projections1, projections2,projections3, projections4,colors, markerStyles);
  //drawGraphs1(gPad,projections1, projections3,colors1, markerStyles1);

  TLegend *legend = new TLegend(0.20, 0.67, 0.48,0.85);
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
  legend->AddEntry(marker2, "0.5<R_{T}<1.5", "lp");
  legend->Draw(); 

  // 下半部分绘制除以Toward区域数据的结果
  c1->cd(2); // 切换到下半部分
  gPad->SetTicks(1,1); // 开启X/Y轴的双边刻度
  
  gPad->SetRightMargin(0.05); // 减小右侧边距
  TH2D *axisFrame2 = new TH2D("axisFrame2", "; ; ", 100, 0,3, 100, 0, 1.3); // Y轴范围调整为0到2
  axisFrame2->GetYaxis()->SetTickLength(0.02); // 调整刻度长度
  axisFrame2->GetXaxis()->SetTickLength(0.02);
  axisFrame2->GetXaxis()->SetLabelSize(0.06); // 将X轴标签字体大小调整为0.045
  axisFrame2->GetYaxis()->SetLabelSize(0.06); // 将Y轴标签字体大小调整为0.045
  axisFrame2->Draw("axis");
  TPaveText *text2 = new TPaveText(0.54, 0.05, 0.575, 0.06, "NDC");
  text2->AddText("p_{T}(GeV/c)");
  text2->SetFillColor(0);
   text2->SetTextSize(0.08); // 增大字体大小
  text2->SetTextAlign(22); // 文本居中对齐
  text2->Draw("same");

  TLatex *label2 = new TLatex();
  label2->SetTextSize(0.08);
  label2->SetTextAngle(90);
  label2->SetTextAlign(22); // 设置文本居中对齐
  label2->DrawLatexNDC(0.05, 0.5, "Ratio to Toward"); //

 drawGraphs1(gPad,projections1, projections2,projections3, projections4,colors, markerStyles);

  TLegend *legend2 = new TLegend(0.65, 0.70, 0.82, 0.81);
  legend2->SetNColumns(1);
  legend2->SetFillStyle(0);
  legend2->SetBorderSize(0);
  //legend2->AddEntry(marker3, "0<R_{T}<0.5", "lp");
  //legend2->AddEntry(marker4, "0.5<R_{T}<1.5", "lp");
  legend2->Draw();

  // 在纵坐标等于1处增加一条虚线
  TLine *line2 = new TLine(0, 1, 3, 1);
  line2->SetLineStyle(2); // 设置虚线样式
  line2->SetLineColor(kBlack);
  line2->Draw();
  
  c1->SaveAs("1K.png");
  //c1->SaveAs("proton_CR_Pt.pdf");
}
