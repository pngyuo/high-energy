#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>

TGraph* drawConnectedPoints(TH1D* hist, Color_t color, Style_t markerStyle,double& averageNT) {
  int n = hist->GetNbinsX();
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> ey; // y误差

  for (int i = 1; i <= n; ++i) {
    double binContent = hist->GetBinContent(i);
    double binCenter = hist->GetBinCenter(i);
    if (binContent != 0&&binCenter<30.0/averageNT) { // 只考虑非零数据点
      double binCenter = hist->GetBinCenter(i);
      double binError = hist->GetBinError(i); // 获取bin的误差
      cout<<"binError: "<<binError<<endl;
      if (binCenter < 40.0/averageNT) { 
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
TH1D* getCentPtProjections(const TString& filename, const TString& histname, double& averageNT, int rebinFactor) {
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
    cout<<"averageNT: "<<averageNT<<endl;

    // 投影到Z轴，并限制X和Y轴的范围
    TH1D* histProjZ = hist3D->ProjectionZ("histProjZ", 1, hist3D->GetNbinsX(), 1, hist3D->GetNbinsY());

    // Rebin the histogram
    histProjZ->Rebin(rebinFactor);

    // 获取原始直方图的bins数量和X轴范围
    int nbins = histProjZ->GetNbinsX();
    double xmin = histProjZ->GetXaxis()->GetBinLowEdge(1) / averageNT;
    double xmax = histProjZ->GetXaxis()->GetBinUpEdge(nbins) / averageNT;

    // 创建一个新的TH1D对象，其横轴范围缩小了averageNT倍
    TH1D *scaledHistProjZ = new TH1D("scaledHistProjZ", "Scaled Projection Z", nbins, xmin, xmax);

    // 填充新的TH1D对象
    for (int i = 1; i <= nbins; ++i) {
        double binContent = histProjZ->GetBinContent(i);
        scaledHistProjZ->SetBinContent(i, binContent);
        scaledHistProjZ->SetBinError(i, histProjZ->GetBinError(i));
    }

    // 设置新的TH1D对象的X轴和Y轴标签
    scaledHistProjZ->GetXaxis()->SetTitle(histProjZ->GetXaxis()->GetTitle());
    scaledHistProjZ->GetYaxis()->SetTitle(histProjZ->GetYaxis()->GetTitle());

    return scaledHistProjZ; // 返回新的缩放后的直方图对象
}

TGraph*drawGraphs( TH1D* projections, TH1D* projections1,TVirtualPad *pad, Color_t color, Style_t markerStyle,double& averageNT) {
    pad->cd(); // 切换到当前子画布
      projections->Divide(projections1);
      TGraph* graph = drawConnectedPoints(projections, color, markerStyle,averageNT);
    return graph;
}


void drawK1() {
    TCanvas *c1 = new TCanvas("c1", "c1", 800,600); // 创建Canvas

    TString filenames[3] = {"hist_outputallFSI.root", "hist_outputnoFSI.root", "hist_outputnohFSI.root"};
    TString histname = "hKCh_dPhi0";
    TString histname1= "hPiCh_dPhi0";
    TString filenames1= "k_pi_to.root";

    std::vector<Color_t> colors = {kRed, kBlue, kMagenta};
    std::vector<Style_t> markerStyles = {20, 21, 22};
    int rebinFactor = 3;
    double averageNT;
    gStyle->SetOptStat(kFALSE);
    
    TH2D *axisFrame = new TH2D("axisFrame", "; ; ", 100, -0.1, 5, 100, 0.06,0.25);
    axisFrame->Draw("axis");  // 只绘制坐标轴，不绘制直方图
for(int i=0;i<3;++i){
            Color_t color = colors[i];
            Style_t markerStyle = markerStyles[i];
            TH1D* projections = getCentPtProjections(filenames[i], histname, averageNT, rebinFactor); // 获取投影并rebin
            TH1D* projections1 = getCentPtProjections(filenames[i], histname1, averageNT, rebinFactor); // 获取投影并rebin
            TGraph* graph1 = drawGraphs(projections, projections1, gPad, color, markerStyle, averageNT); // 分析并绘制每个投影        TFile* infile = new TFile(filenames1[j], "READ");
                TFile* infile = new TFile(filenames1, "READ");
        TGraphErrors *graph2 = (TGraphErrors*)infile->Get("Pt");
        graph2->SetMarkerStyle(24);
        graph2->SetMarkerSize(1);
        graph2->SetLineColor(kBlack);
        graph2->SetMarkerColor(graph2->GetLineColor());
        graph2->Draw("P"); 
        graph1->Draw("LP same"); // 绘制到同一画布上


}
        addText2(0.4, 0.8, "Toward");
    TLegend *legend = new TLegend(0.15, 0.62, 0.6, 0.83);
    legend->SetNColumns(1);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
 
 TPaveText *text = new TPaveText(0.4, 0.01, 0.6, 0.06, "NDC");
    text->AddText("R_{T}");
    text->SetFillColor(0);
    text->SetTextAlign(22); // 文本居中对齐
    text->Draw("same");

TLatex* latex = new TLatex(0.01, 0.6, "#Kappa/#pi"); // x, y 为画布的相对坐标（NDC）
latex->SetNDC(kTRUE); // 使用相对坐标
latex->SetTextSize(0.05); // 设置字体大小
latex->SetTextAngle(90); // 
latex->SetTextAlign(23); 
latex->Draw();
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

    c1->SaveAs("1k_to.png");
    //c1->SaveAs("p_pi.pdf");
}