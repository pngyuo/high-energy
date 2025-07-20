#include <iostream>
#include <vector>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TPaveText.h>

TGraphErrors* drawConnectedPoints(TH1D* hist, Color_t color, Style_t markerStyle) {
    int n = hist->GetNbinsX();
    double* x = new double[n];
    double* y = new double[n];
    double* ey = new double[n]; // y误差数组

    int count = 0; // 用于记录满足条件的点的数量
    for (int i = 1; i <= n; ++i) {
        double binCenter = hist->GetBinCenter(i);
        if (binCenter < 3) { // 只保留 x < 3 的点
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
  void addText1(Double_t x, Double_t y, const char* text)
  {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);  // 将坐标转换为百分比
    latex->SetTextSize(0.05);  // 设置字体大小
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
    latex->SetTextSize(0.10);  // 设置字体大小
    latex->Draw();
  }

  std::vector<TH1D*> getCentPtProjections(const TString& filename, const TString& histname) {
    TFile* infile = new TFile(filename, "READ");
    TH3D* hist3D = (TH3D*)infile->Get(histname);
    TH1D *chargedParticlesCount = (TH1D*)infile->Get("nTHist");
    
    double totalChargedParticles = 0; // 加权总和
    double nEvents = 0;  // 总的事件数
    for (int bin = 1; bin <= chargedParticlesCount->GetNbinsX(); bin++) {
        double binCenter = chargedParticlesCount->GetXaxis()->GetBinCenter(bin);
        double binContent = chargedParticlesCount->GetBinContent(bin);
        totalChargedParticles += binCenter * binContent;
        nEvents += binContent;
    }
    double averageNT = totalChargedParticles / nEvents;
    std::cout << "averageNT: " << averageNT << std::endl;
    
    std::vector<TH1D*> projections;
    double rtRanges[2] = {0, 2.5};  // RT 的下限
    double rtUpperRanges[2] = {0.5, 5};  // RT 的上限
  
    for (int i = 0; i < 2; i++) {
        double ntRangeLow = rtRanges[i] * averageNT;
        double ntRangeUp = rtUpperRanges[i] * averageNT;

        int binLow = hist3D->GetZaxis()->FindBin(ntRangeLow);
        int binUp = hist3D->GetZaxis()->FindBin(ntRangeUp);
        binUp -= 1;

        // 计算该RT范围内的事件数
        double totalEventsInRange = 0;
        for (int binIndex = binLow; binIndex <= binUp; ++binIndex) {
            totalEventsInRange += chargedParticlesCount->GetBinContent(binIndex);
        }
        
        // 输出事件数信息
        std::cout << "RT Range [" << rtRanges[i] << ", " << rtUpperRanges[i] 
                  << "] has " << totalEventsInRange << " events" << std::endl;

        TH1D* histProjX = hist3D->ProjectionX(Form("histProjX_%i", i + 1), 1, hist3D->GetNbinsY(), binLow, binUp);

        histProjX->Scale(1.0 / totalEventsInRange, "width");
        projections.push_back(histProjX);
    }
    return projections;
}

TH1D* getCentPtProjectionsAll(const TString& filename, const TString& histname) {
    TFile* infile = new TFile(filename, "READ");
    TH3D* hist3D = (TH3D*)infile->Get(histname);
    TH1D *chargedParticlesCount = (TH1D*)infile->Get("nTHist");
    
    double totalChargedParticles = 0;
    double nEvents = 0;
    for (int bin = 1; bin <= chargedParticlesCount->GetNbinsX(); bin++) {
        double binCenter = chargedParticlesCount->GetXaxis()->GetBinCenter(bin);
        double binContent = chargedParticlesCount->GetBinContent(bin);
        totalChargedParticles += binCenter * binContent;
        nEvents += binContent;
    }
    double averageNT = totalChargedParticles / nEvents;
    std::cout << "averageNT: " << averageNT << std::endl;

    // 计算总事件数
    double totalEvents = 0;
    for (int binIndex = 1; binIndex <= chargedParticlesCount->GetNbinsX(); ++binIndex) {
        totalEvents += chargedParticlesCount->GetBinContent(binIndex);
    }
    
    // 输出总事件数
    std::cout << "Total events: " << totalEvents << std::endl;

    TH1D* histProjX = hist3D->ProjectionX("histProjX_all", 1, hist3D->GetNbinsY());
    
    histProjX->Scale(1.0 / totalEvents, "width");
    return histProjX;
}

void draw() {
    TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1200);
    c1->SetTicks(1, 1);
    c1->SetMargin(0.17, 0.02, 0.18, 0.1);
    c1->Divide(2, 3, 0, 0);
    gStyle->SetOptStat(kFALSE);
    
    // 通用样式设置
    std::vector<Color_t> colors = {kBlack, kRed, kBlue}; // 三种曲线颜色
    std::vector<Style_t> markerStyles = {20, 21, 22};      // 三种标记样式
    
    // ============= Pad 1 =============
    c1->cd(1);
    {
        gPad->SetTicks(1,1);
        gPad->SetLogy();
        TH2D *axisFrame1 = new TH2D("axisFrame1", "", 100,-0.1, 2.99, 100, 0.002, 16);
        axisFrame1->GetYaxis()->SetTickLength(0.02);
        axisFrame1->GetXaxis()->SetTickLength(0.02);
        axisFrame1->GetXaxis()->SetLabelSize(0.070);
        axisFrame1->GetYaxis()->SetLabelSize(0.070);
        axisFrame1->Draw("axis");

        TString filename = "hist_outputallFSI_liang.root";
        TString histname = "hPiCh_dPhi0";
        TString datafile = "Pi_To.root";

        // 绘制所有RT区间
        TH1D* hist_all = getCentPtProjectionsAll(filename, histname);
        TGraphErrors* graph_all = drawConnectedPoints(hist_all, colors[0], markerStyles[0]);
        graph_all->Draw("LP");
        
        // 绘制分RT区间
        std::vector<TH1D*> projections = getCentPtProjections(filename, histname);
        for (int i = 0; i < projections.size(); i++) {
            TGraphErrors* graph_rt = drawConnectedPoints(projections[i], colors[i+1], markerStyles[i+1]);
            graph_rt->Draw("LP same");
        }

        // 绘制实验数据
        TFile* data_infile = new TFile(datafile, "READ");
        TGraphErrors *data_graph = (TGraphErrors*)data_infile->Get("Pt");
        data_graph->SetMarkerStyle(24);
        data_graph->SetMarkerSize(1.3);
        data_graph->SetLineColor(kBlack);
        data_graph->SetMarkerColor(kBlack);
        data_graph->Draw("P same");
        
        axisFrame1->Draw("axis same");

        // 创建图例
        TLegend *legend = new TLegend(0.23, 0.06, 0.90, 0.40);
        legend->SetNColumns(1);
        legend->SetFillStyle(0);
        legend->SetBorderSize(0);
        legend->SetMargin(0.20);
        legend->SetTextSize(0.076);

        addText3(0.28, 0.42, "AMPT");
        legend->AddEntry(graph_all, "R_{T}>0", "lp");
        legend->AddEntry(drawConnectedPoints(projections[0], colors[1], markerStyles[1]), "0<R_{T}<0.5", "lp");
        legend->AddEntry(drawConnectedPoints(projections[1], colors[2], markerStyles[2]), "2.5<R_{T}<5", "lp");
        legend->Draw();
        addText3(0.52, 0.90, "Toward");
        addText4(0.80, 0.10, "#pi^{+}+#pi^{-}");
        addText2(0.92, 0.03, "(a)");
    }

    // ============= Pad 2 =============
    c1->cd(2);
    {
        gPad->SetTicks(1,1);
        gPad->SetLogy();
        gPad->SetRightMargin(0.01);
        TH2D *axisFrame2 = new TH2D("axisFrame2", "", 100, -0.1, 2.99, 100, 0.002, 16);
        axisFrame2->GetYaxis()->SetTickLength(0.02);
        axisFrame2->GetXaxis()->SetTickLength(0.02);
        axisFrame2->GetXaxis()->SetLabelSize(0.070);
        axisFrame2->GetYaxis()->SetLabelSize(0.070);
        axisFrame2->Draw("axis");

        TString filename = "hist_outputallFSI_liang.root";
        TString histname = "hPiCh_dPhi1";
        TString datafile = "Pi_Tr.root";

        TH1D* hist_all = getCentPtProjectionsAll(filename, histname);
        TGraphErrors* graph_all = drawConnectedPoints(hist_all, colors[0], markerStyles[0]);
        graph_all->Draw("LP");
        
        std::vector<TH1D*> projections = getCentPtProjections(filename, histname);
        for (int i = 0; i < projections.size(); i++) {
            TGraphErrors* graph_rt = drawConnectedPoints(projections[i], colors[i+1], markerStyles[i+1]);
            graph_rt->Draw("LP same");
        }

        TFile* data_infile = new TFile(datafile, "READ");
        TGraphErrors *data_graph = (TGraphErrors*)data_infile->Get("Pt");
        data_graph->SetMarkerStyle(24);
        data_graph->SetMarkerSize(1.3);
        data_graph->SetLineColor(kBlack);
        data_graph->SetMarkerColor(kBlack);
        data_graph->Draw("P same");
        
        axisFrame2->Draw("axis same");

        TLegend *legend2 = new TLegend(0.05, 0.20, 0.72, 0.50);
        legend2->SetNColumns(1);
        legend2->SetFillStyle(0);
        legend2->SetBorderSize(0);
        legend2->SetMargin(0.20);
        legend2->SetTextSize(0.076);
        addText3(0.09, 0.42, "ALICE");
        legend2->AddEntry(data_graph, "R_{T}>0", "p");
        legend2->Draw();
        addText3(0.365, 0.90, "Transverse");
        addText2(0.90, 0.03, "(d)");
    }

    // ============= Pad 3 =============
    c1->cd(3);
    {
        gPad->SetTicks(1,1);
        gPad->SetLogy();
        TH2D *axisFrame3 = new TH2D("axisFrame3", "", 100, -0.1, 2.99, 100, 0.0007, 1.6);
        axisFrame3->GetYaxis()->SetTickLength(0.02);
        axisFrame3->GetXaxis()->SetTickLength(0.02);
        axisFrame3->GetXaxis()->SetLabelSize(0.070);
        axisFrame3->GetYaxis()->SetLabelSize(0.070);
        axisFrame3->Draw("axis");

        TLatex *label3 = new TLatex();
        label3->SetTextSize(0.10);
        label3->SetTextAngle(90);
        label3->SetTextAlign(22);
        label3->DrawLatexNDC(0.05, 0.5, "dN/(dp_{T}dy)(GeV/c)^{-1}");

        TString filename = "hist_outputallFSI_liang.root";
        TString histname = "hKCh_dPhi0";
        TString datafile = "K_To.root";

        TH1D* hist_all = getCentPtProjectionsAll(filename, histname);
        TGraphErrors* graph_all = drawConnectedPoints(hist_all, colors[0], markerStyles[0]);
        graph_all->Draw("LP");
        
        std::vector<TH1D*> projections = getCentPtProjections(filename, histname);
        for (int i = 0; i < projections.size(); i++) {
            TGraphErrors* graph_rt = drawConnectedPoints(projections[i], colors[i+1], markerStyles[i+1]);
            graph_rt->Draw("LP same");
        }

        TFile* data_infile = new TFile(datafile, "READ");
        TGraphErrors *data_graph = (TGraphErrors*)data_infile->Get("Pt");
        data_graph->SetMarkerStyle(24);
        data_graph->SetMarkerSize(1.3);
        data_graph->SetLineColor(kBlack);
        data_graph->SetMarkerColor(kBlack);
        data_graph->Draw("P same");
        
        axisFrame3->Draw("axis same");

        addText4(0.80, 0.11, "K^{+}+K^{-}");
        addText2(0.92, 0.03, "(b)");
    }

    // ============= Pad 4 =============
    c1->cd(4);
    {
        gPad->SetTicks(1,1);
        gPad->SetLogy();
        gPad->SetRightMargin(0.01);
        TH2D *axisFrame4 = new TH2D("axisFrame4", "", 100, -0.1, 2.99, 100, 0.0007, 1.6);
        axisFrame4->GetYaxis()->SetTickLength(0.02);
        axisFrame4->GetXaxis()->SetTickLength(0.02);
        axisFrame4->GetXaxis()->SetLabelSize(0.070);
        axisFrame4->GetYaxis()->SetLabelSize(0.070);
        axisFrame4->Draw("axis");

        TString filename = "hist_outputallFSI_liang.root";
        TString histname = "hKCh_dPhi1";
        TString datafile = "K_Tr.root";

        TH1D* hist_all = getCentPtProjectionsAll(filename, histname);
        TGraphErrors* graph_all = drawConnectedPoints(hist_all, colors[0], markerStyles[0]);
        graph_all->Draw("LP");
        
        std::vector<TH1D*> projections = getCentPtProjections(filename, histname);
        for (int i = 0; i < projections.size(); i++) {
            TGraphErrors* graph_rt = drawConnectedPoints(projections[i], colors[i+1], markerStyles[i+1]);
            graph_rt->Draw("LP same");
        }

        TFile* data_infile = new TFile(datafile, "READ");
        TGraphErrors *data_graph = (TGraphErrors*)data_infile->Get("Pt");
        data_graph->SetMarkerStyle(24);
        data_graph->SetMarkerSize(1.3);
        data_graph->SetLineColor(kBlack);
        data_graph->SetMarkerColor(kBlack);
        data_graph->Draw("P same");
        
        axisFrame4->Draw("axis same");
        addText2(0.90, 0.03, "(e)");
    }

    // ============= Pad 5 =============
    c1->cd(5);
    {
        gPad->SetLogy();
        TH2D *axisFrame5 = new TH2D("axisFrame5", "", 100, -0.1, 2.99, 100, 0.0003, 0.6);
        axisFrame5->GetYaxis()->SetTickLength(0.02);
        axisFrame5->GetXaxis()->SetTickLength(0.02);
        axisFrame5->GetXaxis()->SetLabelSize(0.060);
        axisFrame5->GetYaxis()->SetLabelSize(0.060);
        axisFrame5->Draw("axis");

        TString filename = "hist_outputallFSI_liang.root";
        TString histname = "hProton_dPhi0";
        TString datafile = "P_To.root";
        TPaveText *text6 = new TPaveText(0.535, 0.05, 0.685, 0.06, "NDC");
        text6->AddText("p_{T}(GeV/c)");
        text6->SetFillColor(0);
        text6->SetTextSize(0.08);
        text6->SetTextAlign(22);
        text6->Draw("same");
        TH1D* hist_all = getCentPtProjectionsAll(filename, histname);
        TGraphErrors* graph_all = drawConnectedPoints(hist_all, colors[0], markerStyles[0]);
        graph_all->Draw("LP");
        
        std::vector<TH1D*> projections = getCentPtProjections(filename, histname);
        for (int i = 0; i < projections.size(); i++) {
            TGraphErrors* graph_rt = drawConnectedPoints(projections[i], colors[i+1], markerStyles[i+1]);
            graph_rt->Draw("LP same");
        }

        TFile* data_infile = new TFile(datafile, "READ");
        TGraphErrors *data_graph = (TGraphErrors*)data_infile->Get("Pt");
        data_graph->SetMarkerStyle(24);
        data_graph->SetMarkerSize(1.3);
        data_graph->SetLineColor(kBlack);
        data_graph->SetMarkerColor(kBlack);
        data_graph->Draw("P same");
        
        axisFrame5->Draw("axis same");

        addText4(0.80, 0.28, "p+#bar{p}");
        addText1(0.92, 0.22, "(c)");
    }

    // ============= Pad 6 =============
    c1->cd(6);
    {
        gPad->SetLogy();
        gPad->SetRightMargin(0.01);
        TH2D *axisFrame6 = new TH2D("axisFrame6", "", 100, -0.1, 2.99, 100,  0.0003, 0.6);
        axisFrame6->GetYaxis()->SetTickLength(0.02);
        axisFrame6->GetXaxis()->SetTickLength(0.02);
        axisFrame6->GetXaxis()->SetLabelSize(0.060);
        axisFrame6->GetYaxis()->SetLabelSize(0.060);
        axisFrame6->Draw("axis");


        TPaveText *text6 = new TPaveText(0.435, 0.05, 0.585, 0.06, "NDC");
        text6->AddText("p_{T}(GeV/c)");
        text6->SetFillColor(0);
        text6->SetTextSize(0.08);
        text6->SetTextAlign(22);
        text6->Draw("same");

        TString filename = "hist_outputallFSI_liang.root";
        TString histname = "hProton_dPhi1";
        TString datafile = "P_Tr.root";

        TH1D* hist_all = getCentPtProjectionsAll(filename, histname);
        TGraphErrors* graph_all = drawConnectedPoints(hist_all, colors[0], markerStyles[0]);
        graph_all->Draw("LP");
        
        std::vector<TH1D*> projections = getCentPtProjections(filename, histname);
        for (int i = 0; i < projections.size(); i++) {
            TGraphErrors* graph_rt = drawConnectedPoints(projections[i], colors[i+1], markerStyles[i+1]);
            graph_rt->Draw("LP same");
        }

        TFile* data_infile = new TFile(datafile, "READ");
        TGraphErrors *data_graph = (TGraphErrors*)data_infile->Get("Pt");
        data_graph->SetMarkerStyle(24);
        data_graph->SetMarkerSize(1.3);
        data_graph->SetLineColor(kBlack);
        data_graph->SetMarkerColor(kBlack);
        data_graph->Draw("P same");
        
        axisFrame6->Draw("axis same");
        addText1(0.90, 0.22, "(f)");
    }

    c1->SaveAs("pi_k_p.png");
    c1->SaveAs("pi_k_p.pdf");
}