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
      if (binCenter < 70.0/averageNT) { 
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

double averageNT = 0.0;  // 全局变量，初始值为 0.0
TH1D* getCentPtProjections(const TString& filename, const TString& histname, double& averageNT, int rebinFactor1, int rebinFactor2) {
  TFile* infile = new TFile(filename, "READ");
  TH3D* hist3D = (TH3D*)infile->Get(histname);

  TH1D *chargedParticlesCount = (TH1D*)infile->Get("nTHist");

  double totalChargedParticles = 0; // 加权总和
  double nEvents = 0;  // 总的事件数
  for (int bin = 1; bin <= 36; bin++) {
      double binCenter = chargedParticlesCount->GetXaxis()->GetBinCenter(bin); // 获取当前 bin 中心的 NT 值
      double binContent = chargedParticlesCount->GetBinContent(bin);  // 获取当前 bin 的分布数量
      totalChargedParticles += binCenter * binContent;  // 加权累加
      nEvents += binContent;  // 累加事件数
  }
  averageNT = totalChargedParticles / nEvents;  // 计算加权平均NT
  cout << "averageNT: " << averageNT << endl;
  int binEnd = hist3D->GetXaxis()->FindBin(3.0);
  // 投影到Z轴，并限制X和Y轴的范围
  TH1D* histProjZ = hist3D->ProjectionZ("histProjZ", 1, binEnd, 1, hist3D->GetNbinsY());

  // 获取原始直方图的bins数量和X轴范围
  int nbins = histProjZ->GetNbinsX();
  double xmin = histProjZ->GetXaxis()->GetBinLowEdge(1) / averageNT;
  double xmax = histProjZ->GetXaxis()->GetBinUpEdge(nbins) / averageNT;

  // 定义新的bin边界
  std::vector<double> newBins;
  double currentX = xmin;
  while (currentX < 2.5) {
      newBins.push_back(currentX);
      currentX += rebinFactor1 * (histProjZ->GetXaxis()->GetBinWidth(1) / averageNT);
  }
  while (currentX < 5.0) {
      newBins.push_back(currentX);
      currentX += rebinFactor2 * (histProjZ->GetXaxis()->GetBinWidth(1) / averageNT);
  }
  newBins.push_back(xmax); // 确保包含最后一个bin

  // 创建一个新的TH1D对象，使用自定义的bin边界
  TH1D *scaledHistProjZ = new TH1D("scaledHistProjZ", "Scaled Projection Z", newBins.size() - 1, &newBins[0]);

  // 填充新的TH1D对象
  for (int i = 1; i <= nbins; ++i) {
      double binContent = histProjZ->GetBinContent(i);
      double binCenter = histProjZ->GetBinCenter(i) / averageNT;
      int newBin = scaledHistProjZ->FindBin(binCenter);
      scaledHistProjZ->SetBinContent(newBin, scaledHistProjZ->GetBinContent(newBin) + binContent);
      scaledHistProjZ->SetBinError(newBin, sqrt(scaledHistProjZ->GetBinContent(newBin)));
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
  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800); // 调整画布尺寸为1200x800
  c1->SetTicks(1, 1); // 开启X/Y轴的双边刻度
  c1->SetMargin(0.17, 0.17, 0.18, 0.1); // 调整边距
  TH2D *axisFrame = new TH2D("axisFrame", "; ; ", 100, 0, 5, 100, 0.11,0.21);
axisFrame->GetYaxis()->SetTickLength(0.02); // 调整刻度长度
axisFrame->GetXaxis()->SetTickLength(0.02);
axisFrame->GetXaxis()->SetLabelSize(0.045); // 将X轴标签字体大小调整为0.045
axisFrame->GetYaxis()->SetLabelSize(0.045); // 将Y轴标签字体大小调整为0.045
axisFrame->Draw("axis");  // 只绘制坐标轴，不绘制直方图
    TString filenames[3] = {"hist_outputallFSI_liang.root", "hist_outputnoFSI_liang.root", "hist_outputnohFSI_liang.root"};
    TString histname = "hKCh_dPhi1";
    TString histname1= "hPiCh_dPhi1";
    TString filenames1= "k_pi_tr.root";

    std::vector<Color_t> colors = {kRed, kBlue, kMagenta};
    std::vector<Style_t> markerStyles = {20, 21, 22};
    int rebinFactor1 = 2;
    int rebinFactor2 = 10;
    double averageNT;
    gStyle->SetOptStat(kFALSE);
    
for(int i=0;i<3;++i){
            Color_t color = colors[i];
            Style_t markerStyle = markerStyles[i];
            TH1D* projections = getCentPtProjections(filenames[i], histname, averageNT, rebinFactor1,rebinFactor2); // 获取投影并rebin
            TH1D* projections1 = getCentPtProjections(filenames[i], histname1, averageNT,rebinFactor1,rebinFactor2); // 获取投影并rebin
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
        addText2(0.43, 0.83, "Toward");
        TLegend *legend = new TLegend(0.19, 0.70, 0.90, 0.87);
legend->SetNColumns(1);
legend->SetFillStyle(0);
legend->SetBorderSize(0);
legend->SetMargin(0.10);
 
  TPaveText *text = new TPaveText(0.425, 0.05, 0.575, 0.06, "NDC");
  text->AddText("R_{T}");
  text->SetFillColor(0);
    text->SetTextSize(0.08); // 增大字体大小
  text->SetTextAlign(22); // 文本居中对齐
  text->Draw("same");


  TLatex *label = new TLatex();
  label->SetTextSize(0.08);
  label->SetTextAngle(90);
  label->SetTextAlign(22); // 设置文本居中对齐
  label->DrawLatexNDC(0.05, 0.5, "#Kappa/#pi"); // 使用NDC坐标
    TGraphErrors *marker1 = new TGraphErrors();
    marker1->SetMarkerStyle(24);
    marker1->SetMarkerColor(kBlack);
    marker1->SetMarkerSize(1);
    marker1->SetLineWidth(2); 
    marker1->SetLineColor(kBlack);
    legend->AddEntry(marker1, "ALICE", "p");


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

    c1->SaveAs("k_to.png");
    //c1->SaveAs("p_pi.pdf");
}