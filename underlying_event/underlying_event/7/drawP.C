#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>

TGraphErrors* drawConnectedPoints(TH1D* hist, Color_t color, Style_t markerStyle, int step) {
  int n = hist->GetNbinsX();
  std::vector<double> x, y, ey;  // x, y, 和误差 ey

  for (int i = 1; i <= n; ++i) {
    double binContent = hist->GetBinContent(i);
    if (binContent != 0) { // 只考虑非零数据点
      double binCenter = hist->GetBinCenter(i);
      double binError = hist->GetBinError(i);  // 获取每个bin的误差
      if (i % step == 0) { // 每step个点取一个点
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

double averageNT = 0.0;  // 全局变量，初始值为 0.0

 TH1D* getCentPtProjections(const TString& filename, const TString& histname) {
    TFile* infile = new TFile(filename, "READ");
    TH3D* hist3D = (TH3D*)infile->Get(histname);

    TH1D *chargedParticlesCount = (TH1D*)infile->Get("chargedParticlesHist");

    double totalChargedParticles = 0; // 加权总和
    double nEvents = 0;  // 总的事件数
    for (int bin = 1; bin <= 36; bin++) {
        double binCenter = chargedParticlesCount->GetXaxis()->GetBinCenter(bin); // 获取当前 bin 中心的 NT 值
        double binContent = chargedParticlesCount->GetBinContent(bin);  // 获取当前 bin 的分布数量
        totalChargedParticles += binCenter * binContent;  // 加权累加
        nEvents += binContent;  // 累加事件数
    }
    averageNT = totalChargedParticles / nEvents;  // 计算加权平均NT
    //cout<<"averageNT: "<<averageNT<<endl;

     // 根据 bin 范围投影 X 轴上的直方图
        TH1D* histProjX = hist3D->ProjectionX(Form("histProjX"), 1, hist3D->GetNbinsY());

        double totalEventsInRange = 0;
       for (int binIndex = 1; binIndex <= chargedParticlesCount->GetNbinsX(); ++binIndex) {
        double eventsWithParticleCount = chargedParticlesCount->GetBinContent(binIndex);
        totalEventsInRange += eventsWithParticleCount;
    }

        histProjX->Scale(1.0 / totalEventsInRange, "width");
        
        
    return histProjX;
}

TGraph*drawGraphs( TH1D* projections, TH1D* projections1,TVirtualPad *pad, Color_t color, Style_t markerStyle) {
    pad->cd(); // 切换到当前子画布
      projections->Divide(projections1);
      TGraph* graph = drawConnectedPoints(projections, color, markerStyle,1);
    return graph;
}


void drawP() {
    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 500); // 创建Canvas

    TPaveText *text = new TPaveText(0.4, 0.01, 0.6, 0.1, "NDC");
    text->AddText("P_{T}");
    text->SetFillColor(0);
    text->SetTextAlign(22); // 文本居中对齐
    text->Draw("same");
  TLatex *label = new TLatex();
  label->SetTextSize(0.06);
  label->SetTextAngle(90);
  label->SetTextAlign(22); // 设置文本居中对齐
  label->DrawLatexNDC(0.03, 0.5,"(#Kappa^{+}+#Kappa^{-})/(#pi^{+}+#pi^{-})"); 
  //label->DrawLatexNDC(0.03, 0.5,"(#Kappa^{+}+#Kappa^{-})/(#pi^{+}+#pi^{-})"); 
  
    c1->Divide(3, 1);

    c1->cd(1)->SetPad(0.04, 0.1, 0.39, 1); 
    c1->cd(1)->SetRightMargin(0);
    c1->cd(2)->SetPad(0.39, 0.1, 0.69, 1); 
    c1->cd(2)->SetRightMargin(0);
    c1->cd(2)->SetLeftMargin(0);
    c1->cd(3)->SetPad(0.69, 0.1, 0.99, 1);
    c1->cd(3)->SetRightMargin(0);
    c1->cd(3)->SetLeftMargin(0);

    TString filenames[3] = {"hist_outputMonash.root", "hist_outputRope.root", "hist_outputCR.root"};
    TString histname[3] = {"hKCh_dPhi0", "hKCh_dPhi2", "hKCh_dPhi1"};
    TString histname1[3] = {"hPiCh_dPhi0", "hPiCh_dPhi2", "hPiCh_dPhi1"};
    TString filenames1[3]= {"K_Pi_To.root", "K_Pi_Aw.root", "K_Pi_Tr.root"};

    std::vector<Color_t> colors = {kRed, kBlue, kMagenta};
    std::vector<Style_t> markerStyles = {20, 21, 22};

    gStyle->SetOptStat(kFALSE);
    for (int j = 0; j < 3; j++) {
    c1->cd(j+1);

    TH2D *axisFrame = new TH2D("axisFrame", "; ; ", 100, -0.1, 5.1, 100, 0, 0.6);
    axisFrame->Draw("axis");  // 只绘制坐标轴，不绘制直方图

    for (int i = 0; i < 3; i++) {
 cout << "i: " << i << endl;
        Color_t color = colors[i];
        Style_t markerStyle = markerStyles[i];
        TH1D* projections = getCentPtProjections(filenames[i], histname[j]); // 获取投影
        TH1D* projections1 = getCentPtProjections(filenames[i], histname1[j]); // 获取投影
        TGraph* graph1 = drawGraphs(projections, projections1,gPad, color, markerStyle); // 分析并绘制每个投影
        TFile* infile = new TFile(filenames1[j], "READ");
        TGraphErrors *graph2 = (TGraphErrors*)infile->Get("Pt");
        graph2->SetMarkerStyle(24);
        graph2->SetMarkerSize(1);
        graph2->SetLineColor(kBlack);
        graph2->SetMarkerColor(graph2->GetLineColor());
        graph2->Draw("P"); 
        graph1->Draw("LP same"); // 绘制到同一画布上
    }

    axisFrame->Draw("axis same");

    if (j == 0) {
        addText2(0.4, 0.8, "Toward");
    } else if (j == 1) {
        addText2(0.4, 0.8, "Away");
        //addText4(0.33, 0.25, "#Kappa^{+}+#Kappa^{-}");
    } else if (j == 2) {
        addText2(0.4, 0.8, "Transverse");
    }
}

    c1->cd(1);
    TLegend *legend = new TLegend(0.65, 0.15, 1.1, 0.36);
    legend->SetNColumns(1);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);

    TGraphErrors *marker1 = new TGraphErrors();
    marker1->SetMarkerStyle(24);
    marker1->SetMarkerColor(kBlack);
    marker1->SetMarkerSize(1);
    marker1->SetLineWidth(2); 
    marker1->SetLineColor(kBlack);
    legend->AddEntry(marker1, "ALICE", "lp");

    TGraphErrors *marker2 = new TGraphErrors();
    marker2->SetMarkerStyle(20);
    marker2->SetMarkerColor(kRed);
    marker2->SetMarkerSize(1);
    marker2->SetLineWidth(2); 
    marker2->SetLineColor(kRed);
    legend->AddEntry(marker2, "Monash", "lp");

    TGraphErrors *marker3 = new TGraphErrors();
    marker3->SetMarkerStyle(21);
    marker3->SetMarkerColor(kBlue);
    marker3->SetMarkerSize(1);
    marker3->SetLineWidth(2); 
    marker3->SetLineColor(kBlue);
    legend->AddEntry(marker3, "Rope", "lp");

    TGraphErrors *marker4 = new TGraphErrors();
    marker4->SetMarkerStyle(22);
    marker4->SetMarkerColor(kMagenta);
    marker4->SetMarkerSize(1);
    marker4->SetLineWidth(2); 
    marker4->SetLineColor(kMagenta);
    legend->AddEntry(marker4, "CR", "lp");
    legend->Draw(); 

    c1->SaveAs("k_pi.png");
    //c1->SaveAs("p_pi_P.pdf");
}