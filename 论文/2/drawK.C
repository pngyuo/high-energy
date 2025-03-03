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
    latex->SetTextSize(0.10);  // 设置字体大小
    latex->Draw();
  }

 TH1D*  getCentPtProjections(const TString& filename, const TString& histname) {
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

    // 根据整个Z轴投影X轴上的直方图
    TH1D* histProjX = hist3D->ProjectionX("histProjX", 1, hist3D->GetNbinsY());

    double totalEventsInRange = 0;
    for (int binIndex = 1; binIndex <= chargedParticlesCount->GetNbinsX(); ++binIndex) {
        double eventsWithParticleCount = chargedParticlesCount->GetBinContent(binIndex);
        totalEventsInRange += eventsWithParticleCount;
    }
    
    histProjX->Scale(1.0 / totalEventsInRange, "width");
    return histProjX;
}

TGraph* drawGraphs(TH1D* projections, Color_t color, Style_t markerStyle) {



    TGraph* graph = drawConnectedPoints(projections, color, markerStyle);
 
   return graph;
}

void drawK(){
      TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800); // 调整画布尺寸为1200x800
    c1->SetTicks(1, 1); // 开启X/Y轴的双边刻度
    c1->SetMargin(0.17, 0.17, 0.18, 0.1); // 调整边距
  TH2D *axisFrame = new TH2D("axisFrame", "; ; ", 100, 0, 3, 100, 1e-2, 1);
  axisFrame->GetYaxis()->SetTickLength(0.02); // 调整刻度长度
  axisFrame->GetXaxis()->SetTickLength(0.02);
  axisFrame->GetXaxis()->SetLabelSize(0.045); // 将X轴标签字体大小调整为0.045
  axisFrame->GetYaxis()->SetLabelSize(0.045); // 将Y轴标签字体大小调整为0.045
  axisFrame->Draw("axis");
  TPaveText *text = new TPaveText(0.425, 0.05, 0.575, 0.06, "NDC");
  text->AddText("p_{T}(GeV/c)");
  text->SetFillColor(0);
  text->SetTextSize(0.08); // 增大字体大小
  text->SetTextAlign(22); // 文本居中对齐
  text->Draw("same");

TLatex *label = new TLatex();
label->SetTextSize(0.08);
label->SetTextAngle(90);
label->SetTextAlign(22); // 设置文本居中对齐
label->DrawLatexNDC(0.05, 0.5, "dN/(dp_{T}dy)(GeV/c)^{-1}"); // 使用NDC坐标


TString filenames[1] = {"hist_outputallFSI_liang.root"};
    TString histname1 ="hKCh_dPhi0";
    TString filenames1= "K_To.root";

    std::vector<Color_t> colors = {kRed, kBlue, kMagenta};
    std::vector<Style_t> markerStyles = {20, 21, 22};

    gStyle->SetOptStat(kFALSE);
  

    for (int i = 0; i < 1; i++) {
 cout << "i: " << i << endl;
        Color_t color = colors[i];
        Style_t markerStyle = markerStyles[i];
       TH1D*projections = getCentPtProjections(filenames[i], histname1); // 获取投影
        TGraph* graph1=drawGraphs(projections, color, markerStyle); // 分析并绘制每个投影
        if(i == 0) {
      graph1->Draw("LP"); // L: 折线, P: 数据点
    } else {
      graph1->Draw("LP same"); // 将新的图形绘制在同一画布上
    }
    }
       TFile* infile = new TFile(filenames1, "READ");
        TGraphErrors *graph2 = (TGraphErrors*)infile->Get("Pt");
        graph2->SetMarkerStyle(24);
        graph2->SetMarkerSize(1);
        graph2->SetLineColor(kBlack);
        graph2->SetMarkerColor(graph2->GetLineColor());
        graph2->Draw("P same"); 
    axisFrame->Draw("axis same");
        
      TLegend *legend = new TLegend(0.18, 0.21, 0.56,0.35);
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
    legend->AddEntry(marker2, "allFSI", "lp");

    TGraphErrors *marker3 = new TGraphErrors();
    marker3->SetMarkerStyle(21);
    marker3->SetMarkerColor(kBlue);
    marker3->SetMarkerSize(1);
    marker3->SetLineWidth(2); 
    marker3->SetLineColor(kBlue);
    //legend->AddEntry(marker3, "noFSI", "lp");

    TGraphErrors *marker4 = new TGraphErrors();
    marker4->SetMarkerStyle(22);
    marker4->SetMarkerColor(kMagenta);
    marker4->SetMarkerSize(1);
    marker4->SetLineWidth(2); 
    marker4->SetLineColor(kMagenta);
    //legend->AddEntry(marker4, "nohFSI", "lp");
    legend->Draw(); 

  addText2(0.45, 0.85, "Toward");
  addText4(0.40, 0.24, "#Kappa^{+}+#Kappa^{-}");
  c1->SetLogy(); // 设置y轴为对数刻度
  c1->SaveAs("Ks_1To.png");

}