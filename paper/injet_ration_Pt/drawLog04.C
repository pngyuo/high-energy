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
        if (binCenter < 5&&binCenter >0.25) { // 只保留 x < 3 的点
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
    latex->SetTextSize(0.08);  // 设置字体大小
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

double rtRanges[2] = {0,0.7};  // RT 的范围
    double rtUpperRanges[2] = {0.3,1.1};  // RT 上限

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

// ... 前面的代码保持不变 ...

void drawLog04(){
  TCanvas *c1 = new TCanvas("c1", "c1", 1500, 900);
  c1->SetTicks(1, 1);
  c1->SetMargin(0.17, 0.10, 0.2, 0.2);
  c1->Divide(3,2,0,0);

  // ================ 第一行：K/π ================
  
  // 第一子图 (K/π, no FSI)
  {
      c1->cd(1);
      gPad->SetTicks(1,1);
      TH2D *axisFrame1 = new TH2D("axisFrame1", "; ; ", 100, -0.1, 5.1, 100,  0.01, 0.38);
      axisFrame1->GetYaxis()->SetTickLength(0.02);
      axisFrame1->GetXaxis()->SetTickLength(0.02);
      axisFrame1->GetXaxis()->SetLabelSize(0.075);
      axisFrame1->GetYaxis()->SetLabelSize(0.075);
      axisFrame1->Draw("axis");

      TLatex *label1 = new TLatex();
      label1->SetTextSize(0.09);
      label1->SetTextAngle(90);
      label1->SetTextAlign(22);
      label1->DrawLatexNDC(0.05, 0.44, "#Kappa/#pi");
      addText3(0.79, 0.90, "noFSI");
      addText3(0.24, 0.90, "In-Jet");
      TString filename1 = "hist_outputnoFSI_mb.root"; 
      TString histname1_1 = "hKCh_dPhi0"; 
      TString histname2_1 = "hKCh_dPhi1"; 
      TString histname3_1 = "hPiCh_dPhi0"; 
      TString histname4_1 = "hPiCh_dPhi1"; 

      std::vector<Color_t> colors1 = {kBlack, kRed, kBlue, kMagenta};
      std::vector<Style_t> markerStyles1 = {20, 20, 21, 22};
      gStyle->SetOptStat(kFALSE);

      std::vector<TH1D*> projections1_1 = getCentPtProjections(filename1, histname1_1);
      std::vector<TH1D*> projections2_1 = getCentPtProjections(filename1, histname2_1);
      std::vector<TH1D*> projections3_1 = getCentPtProjections(filename1, histname3_1);
      std::vector<TH1D*> projections4_1 = getCentPtProjections(filename1, histname4_1);

      Int_t rebinFactor1 = 2;
      for (auto& hist : projections1_1) hist->Rebin(rebinFactor1);
      for (auto& hist : projections2_1) hist->Rebin(rebinFactor1);
      for (auto& hist : projections3_1) hist->Rebin(rebinFactor1);
      for (auto& hist : projections4_1) hist->Rebin(rebinFactor1);

      drawGraphs(gPad, projections1_1, projections2_1, projections3_1, projections4_1, colors1, markerStyles1);

  }

  // 第二子图 (K/π, no hadronic FSI)
  {
      c1->cd(2);
      gPad->SetTicks(1,1);
      TH2D *axisFrame2 = new TH2D("axisFrame2", "; ; ", 100, -0.1, 5.1, 100,  0.01, 0.38);
      axisFrame2->GetYaxis()->SetTickLength(0.02);
      axisFrame2->GetXaxis()->SetTickLength(0.02);
      axisFrame2->GetXaxis()->SetLabelSize(0.06);
      axisFrame2->GetYaxis()->SetLabelSize(0.06);
      axisFrame2->Draw("axis");

      addText3(0.71, 0.90, "nohFSI");

      TString filename2 = "hist_outputnohFSI_hadd.root"; 
      TString histname1_2 = "hKCh_dPhi0"; 
      TString histname2_2 = "hKCh_dPhi1"; 
      TString histname3_2 = "hPiCh_dPhi0"; 
      TString histname4_2 = "hPiCh_dPhi1"; 

      std::vector<Color_t> colors2 = {kBlack, kRed, kBlue, kMagenta};
      std::vector<Style_t> markerStyles2 = {20, 20, 21, 22};
      gStyle->SetOptStat(kFALSE);

      std::vector<TH1D*> projections1_2 = getCentPtProjections(filename2, histname1_2);
      std::vector<TH1D*> projections2_2 = getCentPtProjections(filename2, histname2_2);
      std::vector<TH1D*> projections3_2 = getCentPtProjections(filename2, histname3_2);
      std::vector<TH1D*> projections4_2 = getCentPtProjections(filename2, histname4_2);

      Int_t rebinFactor2 = 2;
      for (auto& hist : projections1_2) hist->Rebin(rebinFactor2);
      for (auto& hist : projections2_2) hist->Rebin(rebinFactor2);
      for (auto& hist : projections3_2) hist->Rebin(rebinFactor2);
      for (auto& hist : projections4_2) hist->Rebin(rebinFactor2);

      drawGraphs(gPad, projections1_2, projections2_2, projections3_2, projections4_2, colors2, markerStyles2);

      TLegend *legend2 = new TLegend(0.02, 0.65, 0.39,0.90);
      legend2->SetNColumns(1);
      legend2->SetFillStyle(0);
      legend2->SetBorderSize(0);

      TGraphErrors *marker1_2 = new TGraphErrors();
      marker1_2->SetMarkerStyle(20);
      marker1_2->SetMarkerColor(kBlack);
      marker1_2->SetMarkerSize(1);
      marker1_2->SetLineWidth(2); 
      marker1_2->SetLineColor(kBlack);
      legend2->AddEntry(marker1_2, "0<R_{T}<0.3", "lp");

      TGraphErrors *marker2_2 = new TGraphErrors();
      marker2_2->SetMarkerStyle(20);
      marker2_2->SetMarkerColor(kRed);
      marker2_2->SetMarkerSize(1);
      marker2_2->SetLineWidth(2); 
      marker2_2->SetLineColor(kRed);
      legend2->AddEntry(marker2_2, "0.7<R_{T}<1.1", "lp");
      legend2->Draw(); 
  }

  // 第三子图 (K/π, all FSI)
  {
      c1->cd(3);
      gPad->SetTicks(1,1);
      gPad->SetRightMargin(0.01);
      TH2D *axisFrame3 = new TH2D("axisFrame3", "; ; ", 100, -0.1, 5.1, 100,  0.01, 0.38);
      axisFrame3->GetYaxis()->SetTickLength(0.02);
      axisFrame3->GetXaxis()->SetTickLength(0.02);
      axisFrame3->GetXaxis()->SetLabelSize(0.06);
      axisFrame3->GetYaxis()->SetLabelSize(0.06);
      axisFrame3->Draw("axis");
      addText3(0.74, 0.90, "allFSI");

      TString filename3 = "hist_outputallFSI_more.root"; 
      TString histname1_3 = "hKCh_dPhi0"; 
      TString histname2_3 = "hKCh_dPhi1"; 
      TString histname3_3 = "hPiCh_dPhi0"; 
      TString histname4_3 = "hPiCh_dPhi1"; 

      std::vector<Color_t> colors3 = {kBlack, kRed, kBlue, kMagenta};
      std::vector<Style_t> markerStyles3 = {20, 20, 21, 22};
      gStyle->SetOptStat(kFALSE);

      std::vector<TH1D*> projections1_3 = getCentPtProjections(filename3, histname1_3);
      std::vector<TH1D*> projections2_3 = getCentPtProjections(filename3, histname2_3);
      std::vector<TH1D*> projections3_3 = getCentPtProjections(filename3, histname3_3);
      std::vector<TH1D*> projections4_3 = getCentPtProjections(filename3, histname4_3);

      Int_t rebinFactor3 = 2;
      for (auto& hist : projections1_3) hist->Rebin(rebinFactor3);
      for (auto& hist : projections2_3) hist->Rebin(rebinFactor3);
      for (auto& hist : projections3_3) hist->Rebin(rebinFactor3);
      for (auto& hist : projections4_3) hist->Rebin(rebinFactor3);

      drawGraphs(gPad, projections1_3, projections2_3, projections3_3, projections4_3, colors3, markerStyles3);
  }

  // ================ 第二行：p/π ================
  
  // 第四子图 (p/π, no FSI)
  {
      c1->cd(4);
      gPad->SetTicks(1,1);
      TH2D *axisFrame4 = new TH2D("axisFrame4", "; ; ", 100, -0.1, 5.1, 100,  0, 0.33);
      axisFrame4->GetYaxis()->SetTickLength(0.02);
      axisFrame4->GetXaxis()->SetTickLength(0.02);
      axisFrame4->GetXaxis()->SetLabelSize(0.06);
      axisFrame4->GetYaxis()->SetLabelSize(0.06);
      axisFrame4->Draw("axis");

      TLatex *label4 = new TLatex();
      label4->SetTextSize(0.08);
      label4->SetTextAngle(90);
      label4->SetTextAlign(22);
      label4->DrawLatexNDC(0.05, 0.58, "p/#pi");
      //addText2(0.52, 0.85, "In-Jet");
      TString filename4 = "hist_outputnoFSI_mb.root"; 
      TString histname1_4 = "hProton_dPhi0"; 
      TString histname2_4 = "hProton_dPhi1"; 
      TString histname3_4 = "hPiCh_dPhi0"; 
      TString histname4_4 = "hPiCh_dPhi1"; 

      std::vector<Color_t> colors4 = {kBlack, kRed, kBlue, kMagenta};
      std::vector<Style_t> markerStyles4 = {20, 20, 21, 22};
      gStyle->SetOptStat(kFALSE);

      std::vector<TH1D*> projections1_4 = getCentPtProjections(filename4, histname1_4);
      std::vector<TH1D*> projections2_4 = getCentPtProjections(filename4, histname2_4);
      std::vector<TH1D*> projections3_4 = getCentPtProjections(filename4, histname3_4);
      std::vector<TH1D*> projections4_4 = getCentPtProjections(filename4, histname4_4);

      Int_t rebinFactor4 = 2;
      for (auto& hist : projections1_4) hist->Rebin(rebinFactor4);
      for (auto& hist : projections2_4) hist->Rebin(rebinFactor4);
      for (auto& hist : projections3_4) hist->Rebin(rebinFactor4);
      for (auto& hist : projections4_4) hist->Rebin(rebinFactor4);

      drawGraphs(gPad, projections1_4, projections2_4, projections3_4, projections4_4, colors4, markerStyles4);
  }

  // 第五子图 (p/π, no hadronic FSI)
  {
      c1->cd(5);
      gPad->SetTicks(1,1);
      TH2D *axisFrame5 = new TH2D("axisFrame5", "; ; ", 100, -0.1, 5.1, 100,  0, 0.33);
      axisFrame5->GetYaxis()->SetTickLength(0.02);
      axisFrame5->GetXaxis()->SetTickLength(0.02);
      axisFrame5->GetXaxis()->SetLabelSize(0.06);
      axisFrame5->GetYaxis()->SetLabelSize(0.06);
      axisFrame5->Draw("axis");

      TString filename5 = "hist_outputnohFSI_hadd.root"; 
      TString histname1_5 = "hProton_dPhi0"; 
      TString histname2_5 = "hProton_dPhi1"; 
      TString histname3_5 = "hPiCh_dPhi0"; 
      TString histname4_5 = "hPiCh_dPhi1"; 

      TPaveText *text2 = new TPaveText(0.42, 0.05, 0.595, 0.06, "NDC");
      text2->AddText("p_{T}(GeV/c)");
      text2->SetFillColor(0);
       text2->SetTextSize(0.10); // 增大字体大小
      text2->SetTextAlign(22); // 文本居中对齐
      text2->Draw("same");

      std::vector<Color_t> colors5 = {kBlack, kRed, kBlue, kMagenta};
      std::vector<Style_t> markerStyles5 = {20, 20, 21, 22};
      gStyle->SetOptStat(kFALSE);

      std::vector<TH1D*> projections1_5 = getCentPtProjections(filename5, histname1_5);
      std::vector<TH1D*> projections2_5 = getCentPtProjections(filename5, histname2_5);
      std::vector<TH1D*> projections3_5 = getCentPtProjections(filename5, histname3_5);
      std::vector<TH1D*> projections4_5 = getCentPtProjections(filename5, histname4_5);

      Int_t rebinFactor5 = 2;
      for (auto& hist : projections1_5) hist->Rebin(rebinFactor5);
      for (auto& hist : projections2_5) hist->Rebin(rebinFactor5);
      for (auto& hist : projections3_5) hist->Rebin(rebinFactor5);
      for (auto& hist : projections4_5) hist->Rebin(rebinFactor5);

      drawGraphs(gPad, projections1_5, projections2_5, projections3_5, projections4_5, colors5, markerStyles5);
  }

  // 第六子图 (p/π, all FSI)
  {
      c1->cd(6);
      gPad->SetTicks(1,1);
      gPad->SetRightMargin(0.01);
      TH2D *axisFrame6 = new TH2D("axisFrame6", "; ; ", 100, -0.1, 5.1, 100,  0, 0.33);
      axisFrame6->GetYaxis()->SetTickLength(0.02);
      axisFrame6->GetXaxis()->SetTickLength(0.02);
      axisFrame6->GetXaxis()->SetLabelSize(0.06);
      axisFrame6->GetYaxis()->SetLabelSize(0.06);
      axisFrame6->Draw("axis");

      TString filename6 = "hist_outputallFSI_more.root"; 
      TString histname1_6 = "hProton_dPhi0"; 
      TString histname2_6 = "hProton_dPhi1"; 
      TString histname3_6 = "hPiCh_dPhi0"; 
      TString histname4_6 = "hPiCh_dPhi1"; 

      std::vector<Color_t> colors6 = {kBlack, kRed, kBlue, kMagenta};
      std::vector<Style_t> markerStyles6 = {20, 20, 21, 22};
      gStyle->SetOptStat(kFALSE);

      std::vector<TH1D*> projections1_6 = getCentPtProjections(filename6, histname1_6);
      std::vector<TH1D*> projections2_6 = getCentPtProjections(filename6, histname2_6);
      std::vector<TH1D*> projections3_6 = getCentPtProjections(filename6, histname3_6);
      std::vector<TH1D*> projections4_6 = getCentPtProjections(filename6, histname4_6);

      Int_t rebinFactor6 = 2;
      for (auto& hist : projections1_6) hist->Rebin(rebinFactor6);
      for (auto& hist : projections2_6) hist->Rebin(rebinFactor6);
      for (auto& hist : projections3_6) hist->Rebin(rebinFactor6);
      for (auto& hist : projections4_6) hist->Rebin(rebinFactor6);

      drawGraphs(gPad, projections1_6, projections2_6, projections3_6, projections4_6, colors6, markerStyles6);
  }

  c1->SaveAs("K_and_P_hadd04015_rebin.png");
  c1->SaveAs("K_and_P_hadd04015_rebin.pdf");
}