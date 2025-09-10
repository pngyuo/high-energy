#include <iostream>
#include <vector>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TPaveText.h>

// 改进的分区重绑定函数：正确处理误差
TH1* rebinHistogramThreshold(TH1* hist, double threshold, int rebinFactorBelow) {
    if (!hist) return nullptr;

    TAxis* axis = hist->GetXaxis();
    int nBins = axis->GetNbins();
    
    // 找到阈值点所在的bin
    int thresholdBin = axis->FindBin(threshold);
    while (thresholdBin > 1 && axis->GetBinCenter(thresholdBin) >= threshold) {
        thresholdBin--;
    }

    // 分区处理
    std::vector<double> newBins;
    newBins.push_back(axis->GetBinLowEdge(1)); // 起始边界

    // 处理pT<2.5区域：应用rebin
    for (int i = 1; i <= thresholdBin; i += rebinFactorBelow) {
        int endBin = std::min(i + rebinFactorBelow - 1, thresholdBin);
        newBins.push_back(axis->GetBinUpEdge(endBin));
    }

    // 处理pT≥2.5区域：保持原bin
    for (int i = thresholdBin + 1; i <= nBins; i++) {
        newBins.push_back(axis->GetBinUpEdge(i));
    }

    // 创建新直方图
    TH1D* hnew = new TH1D(
        Form("%s_rebinned_threshold", hist->GetName()), 
        hist->GetTitle(),
        newBins.size() - 1,
        &newBins[0]
    );
    hnew->SetDirectory(0);
    hnew->Sumw2(); // 启用误差计算

    // 正确填充数据并计算误差
    for (int bin = 1; bin <= hnew->GetNbinsX(); bin++) {
        double lowEdge = hnew->GetBinLowEdge(bin);
        double upEdge = hnew->GetBinLowEdge(bin+1);
        
        int firstBinOld = axis->FindBin(lowEdge + 1e-5); // 避免边界问题
        int lastBinOld = axis->FindBin(upEdge - 1e-5);
        
        double content = 0;
        double errorSq = 0;
        
        // 累加原始bin的内容和误差平方
        for (int j = firstBinOld; j <= lastBinOld; j++) {
            content += hist->GetBinContent(j);
            double err = hist->GetBinError(j);
            errorSq += err * err;
        }
        
        hnew->SetBinContent(bin, content);
        hnew->SetBinError(bin, TMath::Sqrt(errorSq));
    }

    return hnew;
}

// 创建连接点的TGraphErrors (只保留x<3的点)
TGraphErrors* drawConnectedPoints(TH1* hist, Color_t color, Style_t markerStyle) {
    int n = hist->GetNbinsX();
    double* x = new double[n];
    double* y = new double[n];
    double* ey = new double[n]; // y误差数组

    int count = 0; // 记录有效点数量
    for (int i = 1; i <= n; ++i) {
        double binCenter = hist->GetBinCenter(i);
        if (binCenter < 3.0) { // 只保留x<3的点
            x[count] = binCenter;
            y[count] = hist->GetBinContent(i);
            ey[count] = hist->GetBinError(i);
            count++;
        }
    }

    TGraphErrors* graph = new TGraphErrors(count, x, y, 0, ey);
    graph->SetLineColor(color);
    graph->SetLineWidth(2);
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerColor(color);

    delete[] x;
    delete[] y;
    delete[] ey;

    return graph;
}

// 文本添加函数
void addText(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);
    latex->SetTextSize(0.05);
    latex->Draw();
}
void addText2(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);
    latex->SetTextSize(0.06);
    latex->Draw();
}
void addText3(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);
    latex->SetTextSize(0.08);
    latex->Draw();
}
void addText4(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);
    latex->SetTextSize(0.10);
    latex->Draw();
}

// 绘制图形：先分区处理直方图，再计算比例
void drawGraphs(const std::vector<TH1*>& projections1, 
                const std::vector<TH1*>& projections2, 
                const std::vector<Color_t>& colors, 
                const std::vector<Style_t>& markerStyles, 
                int rebinFactor = 1) {
    for (size_t i = 0; i < projections1.size(); ++i) {
        TH1* hist1 = rebinHistogramThreshold(projections1[i], 2.5, rebinFactor);
        TH1* hist2 = rebinHistogramThreshold(projections2[i], 2.5, rebinFactor);
        Color_t color = colors[i];
        Style_t markerStyle = markerStyles[i];

        TH1D* ratioHist = (TH1D*)hist1->Clone();
        ratioHist->Divide((TH1D*)hist2); // 计算比例

        TGraphErrors* graphRatio = drawConnectedPoints(ratioHist, color, markerStyle);

        if (i == 0) {
            graphRatio->Draw("LP");
        } else {
            graphRatio->Draw("LP same");
        }

        delete hist1;
        delete hist2;
        delete ratioHist;
    }
}

// 主绘图函数
void drawhadd_rebin_data_transverse() {
    TCanvas *c1 = new TCanvas("c1", "c1", 1290, 800);
    c1->SetTicks(1, 1);
    c1->SetMargin(0.17, 0.10, 0.2, 0.2);
    c1->Divide(3, 2, 0, 0);
    gStyle->SetOptStat(kFALSE);

    // 定义可重用变量
    std::vector<Color_t> colors = {kBlack, kRed, kBlue, kMagenta};
    std::vector<Style_t> markerStyles = {20, 21, 22, 23};
    TFile *file = TFile::Open("output_histograms_transverse.root", "READ");
    int rebinFactor = 2; // 对pT<2.5区域应用rebin因子

    // 第一列：K/π (noFSI)
    {
        c1->cd(1);
        gPad->SetTicks(1, 1);
        TH2D *axisFrame1 = new TH2D("axisFrame1", "; ; ", 100, -0.1, 2.99, 100, 0.01, 0.53);
        axisFrame1->GetYaxis()->SetTickLength(0.02);
        axisFrame1->GetXaxis()->SetTickLength(0.02);
        axisFrame1->GetXaxis()->SetLabelSize(0.065);
        axisFrame1->GetYaxis()->SetLabelSize(0.065);
        axisFrame1->Draw("axis");
        
        TLatex *label = new TLatex();
        label->SetTextSize(0.09);
        label->SetTextAngle(90);
        label->SetTextAlign(22);
        label->DrawLatexNDC(0.05, 0.5, "#Kappa/#pi");
        
        TH1F *hist = (TH1F*)file->Get("histProjX_1_KCh_noFSI");
        TH1F *hist1 = (TH1F*)file->Get("histProjX_2_KCh_noFSI");
        TH1F *hist2 = (TH1F*)file->Get("histProjX_1_PiCh_noFSI");
        TH1F *hist3 = (TH1F*)file->Get("histProjX_2_PiCh_noFSI");
        
        std::vector<TH1*> projections = {hist, hist1};
        std::vector<TH1*> projections1 = {hist2, hist3};
        
        drawGraphs(projections, projections1, colors, markerStyles, rebinFactor);
        
        TString filenames= "k_pi_ratio_transverse_05.root";
        TFile* infile1 = new TFile(filenames, "READ");
        TGraphErrors *graph1 = (TGraphErrors*)infile1->Get("Pt");
        graph1->SetMarkerStyle(24);
        graph1->SetMarkerSize(1.3);
        graph1->SetLineColor(kBlack);
        graph1->SetMarkerColor(graph1->GetLineColor());
        graph1->Draw("P"); 

        TString filenames1= "k_pi_ratio_transverse_255.root";
        TFile* infile2 = new TFile(filenames1, "READ");
        TGraphErrors *graph2 = (TGraphErrors*)infile2->Get("Pt");
        graph2->SetMarkerStyle(24);
        graph2->SetMarkerSize(1.3);
        graph2->SetLineColor(kRed);
        graph2->SetMarkerColor(graph2->GetLineColor());
        graph2->Draw("P same"); 

        TLegend *legend2 = new TLegend(0.6, 0.08, 1.25, 0.35);
        legend2->SetNColumns(1);
        legend2->SetFillStyle(0);
        legend2->SetBorderSize(0);
        legend2->SetMargin(0.20);
        legend2->SetTextSize(0.078);

        TGraphErrors *marker1_2 = new TGraphErrors();
        marker1_2->SetMarkerStyle(20);
        marker1_2->SetMarkerColor(kBlack);
        marker1_2->SetLineColor(kBlack);
        marker1_2->SetLineWidth(2);
        legend2->AddEntry(marker1_2, "0<R_{T}<0.5", "lp");

        TGraphErrors *marker2_2 = new TGraphErrors();
        marker2_2->SetMarkerStyle(20);
        marker2_2->SetMarkerColor(kRed);
        marker2_2->SetLineColor(kRed);
        marker2_2->SetLineWidth(2);
        legend2->AddEntry(marker2_2, "2.5<R_{T}<5", "lp");
        legend2->Draw();
        addText3(0.65, 0.35, "AMPT");
        addText4(0.28, 0.90, "Transverse");
        addText4(0.75, 0.90, "noFSI");
        addText2(0.92, 0.03, "(a)");
    }

    // 第二列：K/π (nohFSI)
    {
        c1->cd(2);
        gPad->SetTicks(1, 1);
        TH2D *axisFrame2 = new TH2D("axisFrame2", "; ; ", 100, -0.1, 2.99, 100, 0.01, 0.53);
        axisFrame2->GetYaxis()->SetTickLength(0.02);
        axisFrame2->GetXaxis()->SetTickLength(0.02);
        axisFrame2->GetXaxis()->SetLabelSize(0.065);
        axisFrame2->GetYaxis()->SetLabelSize(0.065);
        axisFrame2->Draw("axis");

        TH1F *hist = (TH1F*)file->Get("histProjX_1_KCh_nohFSI");
        TH1F *hist1 = (TH1F*)file->Get("histProjX_2_KCh_nohFSI");
        TH1F *hist2 = (TH1F*)file->Get("histProjX_1_PiCh_nohFSI");
        TH1F *hist3 = (TH1F*)file->Get("histProjX_2_PiCh_nohFSI");
        
        std::vector<TH1*> projections = {hist, hist1};
        std::vector<TH1*> projections1 = {hist2, hist3};
        
        drawGraphs(projections, projections1, colors, markerStyles, rebinFactor);

        TString filenames= "k_pi_ratio_transverse_05.root";
        TFile* infile1 = new TFile(filenames, "READ");
        TGraphErrors *graph1 = (TGraphErrors*)infile1->Get("Pt");
        graph1->SetMarkerStyle(24);
        graph1->SetMarkerSize(1.3);
        graph1->SetLineColor(kBlack);
        graph1->SetMarkerColor(graph1->GetLineColor());
        graph1->Draw("P"); 

        TString filenames1= "k_pi_ratio_transverse_255.root";
        TFile* infile2 = new TFile(filenames1, "READ");
        TGraphErrors *graph2 = (TGraphErrors*)infile2->Get("Pt");
        graph2->SetMarkerStyle(24);
        graph2->SetMarkerSize(1.3);
        graph2->SetLineColor(kRed);
        graph2->SetMarkerColor(graph2->GetLineColor());
        graph2->Draw("P same"); 

        addText4(0.75, 0.90, "pFSI");
        addText2(0.90, 0.03, "(c)");
    }

    // 第三列：K/π (allFSI)
    {
        c1->cd(3);
        gPad->SetTicks(1, 1);
        gPad->SetRightMargin(0.01);
        TH2D *axisFrame3 = new TH2D("axisFrame3", "; ; ", 100, -0.1, 2.99, 100,0.01, 0.53);
        axisFrame3->GetYaxis()->SetTickLength(0.02);
        axisFrame3->GetXaxis()->SetTickLength(0.02);
        axisFrame3->GetXaxis()->SetLabelSize(0.065);
        axisFrame3->GetYaxis()->SetLabelSize(0.065);
        axisFrame3->Draw("axis");

        addText4(0.71, 0.90, "allFSI");

        TH1F *hist = (TH1F*)file->Get("histProjX_1_KCh_allFSI");
        TH1F *hist1 = (TH1F*)file->Get("histProjX_2_KCh_allFSI");
        TH1F *hist2 = (TH1F*)file->Get("histProjX_1_PiCh_allFSI");
        TH1F *hist3 = (TH1F*)file->Get("histProjX_2_PiCh_allFSI");
        
        std::vector<TH1*> projections = {hist, hist1};
        std::vector<TH1*> projections1 = {hist2, hist3};
        
        drawGraphs(projections, projections1, colors, markerStyles, rebinFactor);

        TString filenames= "k_pi_ratio_transverse_05.root";
        TFile* infile1 = new TFile(filenames, "READ");
        TGraphErrors *graph1 = (TGraphErrors*)infile1->Get("Pt");
        graph1->SetMarkerStyle(24);
        graph1->SetMarkerSize(1.3);
        graph1->SetLineColor(kBlack);
        graph1->SetMarkerColor(graph1->GetLineColor());
        graph1->Draw("P"); 

        TString filenames1= "k_pi_ratio_transverse_255.root";
        TFile* infile2 = new TFile(filenames1, "READ");
        TGraphErrors *graph2 = (TGraphErrors*)infile2->Get("Pt");
        graph2->SetMarkerStyle(24);
        graph2->SetMarkerSize(1.3);
        graph2->SetLineColor(kRed);
        graph2->SetMarkerColor(graph2->GetLineColor());
        graph2->Draw("P same"); 

        addText3(0.59, 0.35, "ALICE");

        TLegend *legend2 = new TLegend(0.55, 0.08, 1.2, 0.35);
        legend2->SetNColumns(1);
        legend2->SetFillStyle(0);
        legend2->SetBorderSize(0);
        legend2->SetMargin(0.20);
        legend2->SetTextSize(0.078);

        TGraphErrors *marker1_2 = new TGraphErrors();
        marker1_2->SetMarkerStyle(24);
        marker1_2->SetMarkerColor(kBlack);
        marker1_2->SetMarkerSize(1.3);
        marker1_2->SetLineWidth(2); 
        marker1_2->SetLineColor(kBlack);
        legend2->AddEntry(marker1_2, "0<R_{T}<0.5", "p");

        TGraphErrors *marker2_2 = new TGraphErrors();
        marker2_2->SetMarkerStyle(24);
        marker2_2->SetMarkerColor(kRed);
        marker2_2->SetMarkerSize(1.3);
        marker2_2->SetLineWidth(2); 
        marker2_2->SetLineColor(kRed);
        legend2->AddEntry(marker2_2, "2.5<R_{T}<5", "p");
        legend2->Draw();
        addText2(0.90, 0.03, "(e)");
    }

    // 第四列：p/π (noFSI)
    {
        c1->cd(4);
        gPad->SetTicks(1, 1);
        TH2D *axisFrame4 = new TH2D("axisFrame4", "; ; ", 100, -0.1, 2.99, 100, 0, 0.33);
        axisFrame4->GetYaxis()->SetTickLength(0.02);
        axisFrame4->GetXaxis()->SetTickLength(0.02);
        axisFrame4->GetXaxis()->SetLabelSize(0.055);
        axisFrame4->GetYaxis()->SetLabelSize(0.055);
        axisFrame4->Draw("axis");

        TLatex *label3 = new TLatex();
        label3->SetTextSize(0.08);
        label3->SetTextAngle(90);
        label3->SetTextAlign(22);
        label3->DrawLatexNDC(0.05, 0.5, "p/#pi");

        TH1F *hist4 = (TH1F*)file->Get("histProjX_1_Proton_noFSI");
        TH1F *hist5 = (TH1F*)file->Get("histProjX_2_Proton_noFSI");
        TH1F *hist6 = (TH1F*)file->Get("histProjX_1_PiCh_noFSI");
        TH1F *hist7 = (TH1F*)file->Get("histProjX_2_PiCh_noFSI");
        
        std::vector<TH1*> projections = {hist4, hist5};
        std::vector<TH1*> projections1 = {hist6, hist7};
        
        drawGraphs(projections, projections1, colors, markerStyles, rebinFactor);

        TString filenames= "p_pi_ratio_transverse_05.root";
        TFile* infile1 = new TFile(filenames, "READ");
        TGraphErrors *graph1 = (TGraphErrors*)infile1->Get("Pt");
        graph1->SetMarkerStyle(24);
        graph1->SetMarkerSize(1.3);
        graph1->SetLineColor(kBlack);
        graph1->SetMarkerColor(graph1->GetLineColor());
        graph1->Draw("P"); 

        TString filenames1= "p_pi_ratio_transverse_255.root";
        TFile* infile2 = new TFile(filenames1, "READ");
        TGraphErrors *graph2 = (TGraphErrors*)infile2->Get("Pt");
        graph2->SetMarkerStyle(24);
        graph2->SetMarkerSize(1.3);
        graph2->SetLineColor(kRed);
        graph2->SetMarkerColor(graph2->GetLineColor());
        graph2->Draw("P same"); 
        addText(0.92, 0.23, "(b)");
    }

    // 第五列：p/π (nohFSI)
    {
        c1->cd(5);
        gPad->SetTicks(1, 1);
        TH2D *axisFrame5 = new TH2D("axisFrame5", "; ; ", 100, -0.1, 2.99, 100, 0, 0.33);
        axisFrame5->GetYaxis()->SetTickLength(0.02);
        axisFrame5->GetXaxis()->SetTickLength(0.02);
        axisFrame5->GetXaxis()->SetLabelSize(0.055);
        axisFrame5->GetYaxis()->SetLabelSize(0.055);
        axisFrame5->Draw("axis");

        TPaveText *text2 = new TPaveText(0.45, 0.05, 0.6, 0.06, "NDC");
        text2->AddText("p_{T}(GeV/c)");
        text2->SetFillColor(0);
        text2->SetTextSize(0.08);
        text2->SetTextAlign(22);
        text2->Draw("same");

        TH1F *hist8 = (TH1F*)file->Get("histProjX_1_Proton_nohFSI");
        TH1F *hist9 = (TH1F*)file->Get("histProjX_2_Proton_nohFSI");
        TH1F *hist10 = (TH1F*)file->Get("histProjX_1_PiCh_nohFSI");
        TH1F *hist11 = (TH1F*)file->Get("histProjX_2_PiCh_nohFSI");
        
        std::vector<TH1*> projections = {hist8, hist9};
        std::vector<TH1*> projections1 = {hist10, hist11};
        
        drawGraphs(projections, projections1, colors, markerStyles, rebinFactor);

        TString filenames= "p_pi_ratio_transverse_05.root";
        TFile* infile1 = new TFile(filenames, "READ");
        TGraphErrors *graph1 = (TGraphErrors*)infile1->Get("Pt");
        graph1->SetMarkerStyle(24);
        graph1->SetMarkerSize(1.3);
        graph1->SetLineColor(kBlack);
        graph1->SetMarkerColor(graph1->GetLineColor());
        graph1->Draw("P"); 

        TString filenames1= "p_pi_ratio_transverse_255.root";
        TFile* infile2 = new TFile(filenames1, "READ");
        TGraphErrors *graph2 = (TGraphErrors*)infile2->Get("Pt");
        graph2->SetMarkerStyle(24);
        graph2->SetMarkerSize(1.3);
        graph2->SetLineColor(kRed);
        graph2->SetMarkerColor(graph2->GetLineColor());
        graph2->Draw("P same"); 
        addText(0.90, 0.23, "(d)");
    }

    // 第六列：p/π (allFSI)
    {
        c1->cd(6);
        gPad->SetTicks(1, 1);
        gPad->SetRightMargin(0.01);
        TH2D *axisFrame6 = new TH2D("axisFrame6", "; ; ", 100, -0.1, 2.99, 100, 0, 0.33);
        axisFrame6->GetYaxis()->SetTickLength(0.02);
        axisFrame6->GetXaxis()->SetTickLength(0.02);
        axisFrame6->GetXaxis()->SetLabelSize(0.055);
        axisFrame6->GetYaxis()->SetLabelSize(0.055);
        axisFrame6->Draw("axis");

        TH1F *hist12 = (TH1F*)file->Get("histProjX_1_Proton_allFSI");
        TH1F *hist13 = (TH1F*)file->Get("histProjX_2_Proton_allFSI");
        TH1F *hist14 = (TH1F*)file->Get("histProjX_1_PiCh_allFSI");
        TH1F *hist15 = (TH1F*)file->Get("histProjX_2_PiCh_allFSI");
        
        std::vector<TH1*> projections = {hist12, hist13};
        std::vector<TH1*> projections1 = {hist14, hist15};
        
        drawGraphs(projections, projections1, colors, markerStyles, rebinFactor);

        TString filenames= "p_pi_ratio_transverse_05.root";
        TFile* infile1 = new TFile(filenames, "READ");
        TGraphErrors *graph1 = (TGraphErrors*)infile1->Get("Pt");
        graph1->SetMarkerStyle(24);
        graph1->SetMarkerSize(1.3);
        graph1->SetLineColor(kBlack);
        graph1->SetMarkerColor(graph1->GetLineColor());
        graph1->Draw("P"); 

        TString filenames1= "p_pi_ratio_transverse_255.root";
        TFile* infile2 = new TFile(filenames1, "READ");
        TGraphErrors *graph2 = (TGraphErrors*)infile2->Get("Pt");
        graph2->SetMarkerStyle(24);
        graph2->SetMarkerSize(1.3);
        graph2->SetLineColor(kRed);
        graph2->SetMarkerColor(graph2->GetLineColor());
        graph2->Draw("P same"); 
        addText(0.90, 0.23, "(f)");
    }

    c1->SaveAs("pi_k_p_transverse_hadd_rebin_new_data.png");
    c1->SaveAs("pi_k_p_transverse_hadd_rebin_new_data.pdf");
}