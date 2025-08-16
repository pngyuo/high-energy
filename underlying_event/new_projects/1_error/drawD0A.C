#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>


TGraph* drawConnectedPoints(TH1D* hist, Color_t color, Style_t markerStyle, int step, double averageNT) {
  int n = hist->GetNbinsX();
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> ey; // 存储误差值
  
  for (int i = 1; i <= n; ++i) {
    double binContent = hist->GetBinContent(i);
    if (binContent != 0) { // 只考虑非零数据点
      double binCenter = hist->GetBinCenter(i);
      if (i % step == 0 && binCenter < 30.0/averageNT) { // 每step个点取一个点
        x.push_back(binCenter);
        y.push_back(binContent);
        ey.push_back(hist->GetBinError(i)); // 获取误差值
      }
    }
  }
  
  int nPoints = x.size();
  if (nPoints == 0) return nullptr; // 如果没有非零点，则返回nullptr
  
  double* xArr = new double[nPoints];
  double* yArr = new double[nPoints];
  double* eyArr = new double[nPoints]; // 错误数组
  for (int i = 0; i < nPoints; ++i) {
    xArr[i] = x[i];
    yArr[i] = y[i];
    eyArr[i] = ey[i];
  }
  
  TGraphErrors* graph = new TGraphErrors(nPoints, xArr, yArr, 0, eyArr); // 使用误差数组
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

TH1D* getCentPtProjections(const TString& filename, const TString& histname, double& averageNT) {
    TFile* infile = new TFile(filename, "READ");
    TProfile *profile= (TProfile*)infile->Get(histname);

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
    std::cout << "averageNT: " << averageNT << std::endl;

    // 获取TProfile的bins数量和X轴范围
    int nbins = profile->GetNbinsX();
    double xmin = profile->GetXaxis()->GetBinLowEdge(1) / averageNT;
    double xmax = profile->GetXaxis()->GetBinUpEdge(nbins) / averageNT;
    TH1D *hist = new TH1D("hist", "Extracted Profile Data", nbins, xmin, xmax);

    // 填充TH1D对象
    for (int i = 1; i <= nbins; ++i) {
        double binContent = profile->GetBinContent(i);
        double binError = profile->GetBinError(i);
        hist->SetBinContent(i, binContent);
        hist->SetBinError(i, binError);
    }

    // 设置TH1D对象的X轴和Y轴标签
    hist->GetXaxis()->SetTitle(profile->GetXaxis()->GetTitle());
    hist->GetYaxis()->SetTitle(profile->GetYaxis()->GetTitle());
    return hist;
}

TGraph* drawGraphs(TVirtualPad* pad, TH1D* projections, Color_t color, Style_t markerStyle, int step,double averageNT) {
  TGraph* graph = drawConnectedPoints(projections, color, markerStyle, step,averageNT);
  if (graph) {
    graph->Draw("LP"); // L: 折线, P: 数据点
  }
  return graph;
}

void drawD0A(){
  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 500); // 创建Canvas

  TPaveText *text = new TPaveText(0.4, 0.01, 0.6, 0.1, "NDC");
  text->AddText("R_{T}");
  text->SetFillColor(0);
  text->SetTextAlign(22); // 文本居中对齐
  text->Draw("same");

  TLatex *label = new TLatex(0.01, 0.6, "<p_{T}>(GeV/c)");
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

    TString filenames[3] = {"hist_outputallFSI.root", "hist_outputnoFSI.root", "hist_outputnohFSI.root"};
    TString histname[3] = {"pD0AvgPt_NT0_ptcut", "pD0AvgPt_NT2_ptcut", "pD0AvgPt_NT1_ptcut"};
    //TString filenames1[3]= {"pi_toward.root", "pi_away.root", "pi_transverse.root"};

    std::vector<Color_t> colors = {kRed, kBlue, kMagenta};
    std::vector<Style_t> markerStyles = {20, 21, 22};
 
    gStyle->SetOptStat(kFALSE);
    double averageNT;
    for (int j = 0; j < 3; j++) {
        c1->cd(j+1);

        TH2D *axisFrame = new TH2D("axisFrame", "; ; ", 100, -0.1, 5.1, 100, 0.5, 3.9);
        axisFrame->Draw("axis");  // 只绘制坐标轴，不绘制直方图

       for (int i = 0; i < 3; i++) {
  Color_t color = colors[i];
  Style_t markerStyle = markerStyles[i];
  TH1D* projections = getCentPtProjections(filenames[i], histname[j], averageNT); // 获取投影
  TGraph* graph = drawGraphs(gPad, projections, color, markerStyle, 2,averageNT); 
       //TFile* infile = new TFile(filenames1[j], "READ");
        //TGraphErrors *graph1 = (TGraphErrors*)infile->Get("Pt");
        //graph1->SetMarkerStyle(24);
        //graph1->SetMarkerSize(1);
        //graph1->SetLineColor(kBlack);
        //graph1->SetMarkerColor(graph1->GetLineColor());
       //graph1->Draw("P"); 
        graph->Draw("LP SAME"); // 绘制到同一画布上
}

    axisFrame->Draw("axis same");

    if (j == 0) {
        addText2(0.4, 0.8, "Toward");
    } else if (j == 1) {
        addText2(0.4, 0.8, "Away");
        addText4(0.33, 0.25, "D^{0}");
    } else if (j == 2) {
        addText2(0.4, 0.8, "Transverse");
    }
}

    c1->cd(1);
    TLegend *legend = new TLegend(0.65, 0.18, 1.1, 0.39);
    legend->SetNColumns(1);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0); // 设置边框大小为0，去掉边框


    TGraphErrors *marker1 = new TGraphErrors();
    marker1->SetMarkerStyle(24);
    marker1->SetMarkerColor(kBlack);
    marker1->SetMarkerSize(1);
    marker1->SetLineWidth(2); 
    marker1->SetLineColor(kBlack);
    //legend->AddEntry(marker1, "Monash", "lp");

    TGraphErrors *marker2 = new TGraphErrors();
    marker2->SetMarkerStyle(20);
    marker2->SetMarkerColor(kRed);
    marker2->SetMarkerSize(1);
    marker2->SetLineWidth(2); 
    marker2->SetLineColor(kRed);
    legend->AddEntry(marker2, "allFSI", "lp");

    TGraphErrors *marker3 = new TGraphErrors();
    marker3->SetMarkerStyle(21);
    marker3->SetMarkerColor(kBlue);
    marker3->SetMarkerSize(1);
    marker3->SetLineWidth(2); 
    marker3->SetLineColor(kBlue);
    legend->AddEntry(marker3, "noFSI", "lp");

    TGraphErrors *marker4 = new TGraphErrors();
    marker4->SetMarkerStyle(22);
    marker4->SetMarkerColor(kMagenta);
    marker4->SetMarkerSize(1);
    marker4->SetLineWidth(2); 
    marker4->SetLineColor(kMagenta);
    legend->AddEntry(marker4, "nohFSI", "lp");
    legend->Draw(); 

    c1->SaveAs("D0.A.png");
    c1->SaveAs("D0.A.pdf");
}


