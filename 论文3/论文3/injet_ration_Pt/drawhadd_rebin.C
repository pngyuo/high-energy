#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TMath.h>

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

TGraphErrors* drawConnectedPoints(TH1* hist, Color_t color, Style_t markerStyle) {
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

void addText(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);  // 将坐标转换为百分比
    latex->SetTextSize(0.05);  // 设置字体大小
    latex->Draw();
}
void addText2(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);  // 将坐标转换为百分比
    latex->SetTextSize(0.06);  // 设置字体大小
    latex->Draw();
}
void addText3(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);  // 将坐标转换为百分比
    latex->SetTextSize(0.08);  // 设置字体大小
    latex->Draw();
}
void addText5(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);  // 将坐标转换为百分比
    latex->SetTextSize(0.10);  // 设置字体大小
    latex->Draw();
}
void addText4(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);  // 将坐标转换为百分比
    latex->SetTextSize(0.15);  // 设置字体大小
    latex->Draw();
}

void drawGraphs(TVirtualPad *pad, const std::vector<TH1*>& projections1, 
                std::vector<TH1*> projections2, 
                std::vector<TH1*> projections3, 
                std::vector<TH1*> projections4, 
                const std::vector<Color_t>& colors, 
                const std::vector<Style_t>& markerStyles,
                int rebinFactor = 2) {
    pad->cd(); // 切换到当前子画布
    for (size_t i = 0; i < projections1.size(); ++i) {
        TH1* hist1 = rebinHistogramThreshold(projections1[i], 2.5, rebinFactor);
        TH1* hist2 = rebinHistogramThreshold(projections2[i], 2.5, rebinFactor);
        TH1* hist3 = rebinHistogramThreshold(projections3[i], 2.5, rebinFactor);
        TH1* hist4 = rebinHistogramThreshold(projections4[i], 2.5, rebinFactor);
        Color_t color = colors[i];
        Style_t markerStyle = markerStyles[i];
        // 相减操作：hist1 - hist2
        TH1* diffHist1 = (TH1*)hist1->Clone();
        diffHist1->Add(hist2, -1); // 减去hist2
        // 相减操作：hist3 - hist4
        TH1* diffHist2 = (TH1*)hist3->Clone();
        diffHist2->Add(hist4, -1); // 减去hist2
        // 相除操作：hist3 / hist4
        TH1* ratioHist = (TH1*)diffHist1->Clone();
        ratioHist->Divide(diffHist2); // 除以hist4
    
        // 绘制相减和相除的结果
        TGraphErrors* graphRatio = drawConnectedPoints(ratioHist, color, markerStyle);
    
        if(i == 0) {
            graphRatio->Draw("LP");
        } else {
            graphRatio->Draw("LP same");
        }
        
        // 删除临时创建的直方图
        delete diffHist1;
        delete diffHist2;
        delete ratioHist;
    }
}

void drawhadd_rebin() {
    TCanvas *c1 = new TCanvas("c1", "c1", 1290, 800);
    c1->SetTicks(1, 1);
    c1->SetMargin(0.17, 0.10, 0.2, 0.2);
    c1->Divide(3, 2, 0, 0);
    gStyle->SetOptStat(kFALSE);

    // 定义重分箱因子
    const int rebinFactor = 2;

    // 存储重分箱后的直方图指针，用于后续删除
    std::vector<TH1*> rebinnedHistograms;

    // Define reusable variables
    std::vector<Color_t> colors = {kBlack, kRed, kBlue, kMagenta};
    std::vector<Style_t> markerStyles = {20, 21, 22, 23};
    TFile *file = TFile::Open("output_histograms_toward_RT0to11.root", "READ");
    TFile *file1 = TFile::Open("output_histograms_transverse_RT0to11.root", "READ");

    // 第一子图 (K/π, no FSI)
    {
        c1->cd(1);
        gPad->SetTicks(1, 1);
        TH2D *axisFrame1 = new TH2D("axisFrame1", "; ; ", 100, -0.1, 2.99, 100, 0.01, 0.38);
        axisFrame1->GetYaxis()->SetTickLength(0.02);
        axisFrame1->GetXaxis()->SetTickLength(0.02);
        axisFrame1->GetXaxis()->SetLabelSize(0.065);
        axisFrame1->GetYaxis()->SetLabelSize(0.065);
        axisFrame1->Draw("axis");

        TLatex *label1 = new TLatex();
        label1->SetTextSize(0.09);
        label1->SetTextAngle(90);
        label1->SetTextAlign(22);
        label1->DrawLatexNDC(0.05, 0.44, "#Kappa/#pi");
        addText5(0.79, 0.90, "noFSI");
        addText5(0.24, 0.90, "In-Jet");

        std::vector<Color_t> colors1 = {kBlack, kRed, kBlue, kMagenta};
        std::vector<Style_t> markerStyles1 = {20, 20, 21, 22};
        gStyle->SetOptStat(kFALSE);
        
        // 获取并重分箱直方图
        TH1F *hist = (TH1F*)file->Get("histProjX_1_KCh_noFSI");
        TH1F *hist1 = (TH1F*)file->Get("histProjX_2_KCh_noFSI");
        TH1F *hist_1 = (TH1F*)file1->Get("histProjX_1_KCh_noFSI");
        TH1F *hist1_1 = (TH1F*)file1->Get("histProjX_2_KCh_noFSI");
        TH1F *hist2 = (TH1F*)file->Get("histProjX_1_PiCh_noFSI");
        TH1F *hist3 = (TH1F*)file->Get("histProjX_2_PiCh_noFSI");
        TH1F *hist2_1 = (TH1F*)file1->Get("histProjX_1_PiCh_noFSI");
        TH1F *hist3_1 = (TH1F*)file1->Get("histProjX_2_PiCh_noFSI");
    
        std::vector<TH1*> projections;
        projections.push_back(hist);
        projections.push_back(hist1);
        std::vector<TH1*> projections1;
        projections1.push_back(hist_1);
        projections1.push_back(hist1_1);
        std::vector<TH1*> projections2;
        projections2.push_back(hist2);
        projections2.push_back(hist3);
        std::vector<TH1*> projections3;
        projections3.push_back(hist2_1);
        projections3.push_back(hist3_1);

        drawGraphs(gPad, projections, projections1, projections2, projections3, colors1, markerStyles1,rebinFactor);
        addText2(0.92, 0.03, "(a)");
    }

    // 第二子图 (K/π, no hadronic FSI)
    {
        c1->cd(2);
        gPad->SetTicks(1, 1);
        TH2D *axisFrame2 = new TH2D("axisFrame2", "; ; ", 100, -0.1, 2.99, 100, 0.01, 0.38);
        axisFrame2->GetYaxis()->SetTickLength(0.02);
        axisFrame2->GetXaxis()->SetTickLength(0.02);
        axisFrame2->GetXaxis()->SetLabelSize(0.065);
        axisFrame2->GetYaxis()->SetLabelSize(0.065);
        axisFrame2->Draw("axis");

        addText5(0.76, 0.90, "pFSI");

        std::vector<Color_t> colors2 = {kBlack, kRed, kBlue, kMagenta};
        std::vector<Style_t> markerStyles2 = {20, 20, 21, 22};
        gStyle->SetOptStat(kFALSE);

        // 获取并重分箱直方图
        TH1F *hist = (TH1F*)file->Get("histProjX_1_KCh_nohFSI");
        TH1F *hist1 = (TH1F*)file->Get("histProjX_2_KCh_nohFSI");
        TH1F *hist_1 = (TH1F*)file1->Get("histProjX_1_KCh_nohFSI");
        TH1F *hist1_1 = (TH1F*)file1->Get("histProjX_2_KCh_nohFSI");
        TH1F *hist2 = (TH1F*)file->Get("histProjX_1_PiCh_nohFSI");
        TH1F *hist3 = (TH1F*)file->Get("histProjX_2_PiCh_nohFSI");
        TH1F *hist2_1 = (TH1F*)file1->Get("histProjX_1_PiCh_nohFSI");
        TH1F *hist3_1 = (TH1F*)file1->Get("histProjX_2_PiCh_nohFSI");
        
        std::vector<TH1*> projections;
        projections.push_back(hist);
        projections.push_back(hist1);
        std::vector<TH1*> projections1;
        projections1.push_back(hist_1);
        projections1.push_back(hist1_1);
        std::vector<TH1*> projections2;
        projections2.push_back(hist2);
        projections2.push_back(hist3);
        std::vector<TH1*> projections3;
        projections3.push_back(hist2_1);
        projections3.push_back(hist3_1);

        drawGraphs(gPad, projections, projections1, projections2, projections3, colors2, markerStyles2,rebinFactor);
        
        TLegend *legend2 = new TLegend(0.02, 0.62, 0.75, 0.89);
        legend2->SetNColumns(1);
        legend2->SetFillStyle(0);
        legend2->SetBorderSize(0);
        legend2->SetMargin(0.20);
        legend2->SetTextSize(0.078);
  
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
        addText2(0.90, 0.03, "(c)");
    }

    // 第三子图 (K/π, all FSI)
    {
        c1->cd(3);
        gPad->SetTicks(1, 1);
        gPad->SetRightMargin(0.01);
        TH2D *axisFrame3 = new TH2D("axisFrame3", "; ; ", 100, -0.1, 2.99, 100, 0.01, 0.38);
        axisFrame3->GetYaxis()->SetTickLength(0.02);
        axisFrame3->GetXaxis()->SetTickLength(0.02);
        axisFrame3->GetXaxis()->SetLabelSize(0.065);
        axisFrame3->GetYaxis()->SetLabelSize(0.065);
        axisFrame3->Draw("axis");
        addText5(0.74, 0.90, "allFSI");

        std::vector<Color_t> colors3 = {kBlack, kRed, kBlue, kMagenta};
        std::vector<Style_t> markerStyles3 = {20, 20, 21, 22};
        gStyle->SetOptStat(kFALSE);

        // 获取并重分箱直方图
        TH1F *hist = (TH1F*)file->Get("histProjX_1_KCh_allFSI");
        TH1F *hist1 = (TH1F*)file->Get("histProjX_2_KCh_allFSI");
        TH1F *hist_1 = (TH1F*)file1->Get("histProjX_1_KCh_allFSI");
        TH1F *hist1_1 = (TH1F*)file1->Get("histProjX_2_KCh_allFSI");
        TH1F *hist2 = (TH1F*)file->Get("histProjX_1_PiCh_allFSI");
        TH1F *hist3 = (TH1F*)file->Get("histProjX_2_PiCh_allFSI");
        TH1F *hist2_1 = (TH1F*)file1->Get("histProjX_1_PiCh_allFSI");
        TH1F *hist3_1 = (TH1F*)file1->Get("histProjX_2_PiCh_allFSI");
        
        std::vector<TH1*> projections;
        projections.push_back(hist);
        projections.push_back(hist1);
        std::vector<TH1*> projections1;
        projections1.push_back(hist_1);
        projections1.push_back(hist1_1);
        std::vector<TH1*> projections2;
        projections2.push_back(hist2);
        projections2.push_back(hist3);
        std::vector<TH1*> projections3;
        projections3.push_back(hist2_1);
        projections3.push_back(hist3_1);

        drawGraphs(gPad, projections, projections1, projections2, projections3, colors3, markerStyles3,rebinFactor);
        addText2(0.90, 0.03, "(e)");
    }

    // ================ 第二行：p/π ================
    
    // 第四子图 (p/π, no FSI)
    {
        c1->cd(4);
        gPad->SetTicks(1, 1);
        TH2D *axisFrame4 = new TH2D("axisFrame4", "; ; ", 100, -0.1, 2.99, 100, 0, 0.19);
        axisFrame4->GetYaxis()->SetTickLength(0.02);
        axisFrame4->GetXaxis()->SetTickLength(0.02);
        axisFrame4->GetXaxis()->SetLabelSize(0.055);
        axisFrame4->GetYaxis()->SetLabelSize(0.055);
        axisFrame4->Draw("axis");

        TLatex *label4 = new TLatex();
        label4->SetTextSize(0.08);
        label4->SetTextAngle(90);
        label4->SetTextAlign(22);
        label4->DrawLatexNDC(0.05, 0.58, "p/#pi");
        //addText2(0.52, 0.85, "In-Jet");

        std::vector<Color_t> colors4 = {kBlack, kRed, kBlue, kMagenta};
        std::vector<Style_t> markerStyles4 = {20, 20, 21, 22};
        gStyle->SetOptStat(kFALSE);
        
        // 获取并重分箱直方图
        TH1F *hist = (TH1F*)file->Get("histProjX_1_Proton_noFSI");
        TH1F *hist1 = (TH1F*)file->Get("histProjX_2_Proton_noFSI");
        TH1F *hist_1 = (TH1F*)file1->Get("histProjX_1_Proton_noFSI");
        TH1F *hist1_1 = (TH1F*)file1->Get("histProjX_2_Proton_noFSI");
        TH1F *hist2 = (TH1F*)file->Get("histProjX_1_PiCh_noFSI");
        TH1F *hist3 = (TH1F*)file->Get("histProjX_2_PiCh_noFSI");
        TH1F *hist2_1 = (TH1F*)file1->Get("histProjX_1_PiCh_noFSI");
        TH1F *hist3_1 = (TH1F*)file1->Get("histProjX_2_PiCh_noFSI");
        
        std::vector<TH1*> projections;
        projections.push_back(hist);
        projections.push_back(hist1);
        std::vector<TH1*> projections1;
        projections1.push_back(hist_1);
        projections1.push_back(hist1_1);
        std::vector<TH1*> projections2;
        projections2.push_back(hist2);
        projections2.push_back(hist3);
        std::vector<TH1*> projections3;
        projections3.push_back(hist2_1);
        projections3.push_back(hist3_1);
        drawGraphs(gPad, projections, projections1, projections2, projections3, colors4, markerStyles4,rebinFactor);
        addText(0.92, 0.23, "(b)");
    }

    // 第五子图 (p/π, no hadronic FSI)
    {
        c1->cd(5);
        gPad->SetTicks(1, 1);
        TH2D *axisFrame5 = new TH2D("axisFrame5", "; ; ", 100, -0.1, 2.99, 100, 0, 0.19);
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

        std::vector<Color_t> colors5 = {kBlack, kRed, kBlue, kMagenta};
        std::vector<Style_t> markerStyles5 = {20, 20, 21, 22};
        gStyle->SetOptStat(kFALSE);
        // 修正：使用质子(Proton)而非K介子(KCh)
        
        // 获取并重分箱直方图
        TH1F *hist = (TH1F*)file->Get("histProjX_1_Proton_nohFSI");
        TH1F *hist1 = (TH1F*)file->Get("histProjX_2_Proton_nohFSI");
        TH1F *hist_1 = (TH1F*)file1->Get("histProjX_1_Proton_nohFSI");
        TH1F *hist1_1 = (TH1F*)file1->Get("histProjX_2_Proton_nohFSI");
        TH1F *hist2 = (TH1F*)file->Get("histProjX_1_PiCh_nohFSI");
        TH1F *hist3 = (TH1F*)file->Get("histProjX_2_PiCh_nohFSI");
        TH1F *hist2_1 = (TH1F*)file1->Get("histProjX_1_PiCh_nohFSI");
        TH1F *hist3_1 = (TH1F*)file1->Get("histProjX_2_PiCh_nohFSI");
        
        
        std::vector<TH1*> projections;
        projections.push_back(hist);
        projections.push_back(hist1);
        std::vector<TH1*> projections1;
        projections1.push_back(hist_1);
        projections1.push_back(hist1_1);
        std::vector<TH1*> projections2;
        projections2.push_back(hist2);
        projections2.push_back(hist3);
        std::vector<TH1*> projections3;
        projections3.push_back(hist2_1);
        projections3.push_back(hist3_1);

        drawGraphs(gPad, projections, projections1, projections2, projections3, colors5, markerStyles5,rebinFactor);
        addText(0.90, 0.23, "(d)");
    }

    // 第六子图 (p/π, all FSI)
    {
        c1->cd(6);
        gPad->SetTicks(1, 1);
        gPad->SetRightMargin(0.01);
        TH2D *axisFrame6 = new TH2D("axisFrame6", "; ; ", 100, -0.1, 2.99, 100, 0, 0.19);
        axisFrame6->GetYaxis()->SetTickLength(0.02);
        axisFrame6->GetXaxis()->SetTickLength(0.02);
        axisFrame6->GetXaxis()->SetLabelSize(0.055);
        axisFrame6->GetYaxis()->SetLabelSize(0.055);
        axisFrame6->Draw("axis");

        std::vector<Color_t> colors6 = {kBlack, kRed, kBlue, kMagenta};
        std::vector<Style_t> markerStyles6 = {20, 20, 21, 22};
        gStyle->SetOptStat(kFALSE);

        // 获取并重分箱直方图
        TH1F *hist = (TH1F*)file->Get("histProjX_1_Proton_allFSI");
        TH1F *hist1 = (TH1F*)file->Get("histProjX_2_Proton_allFSI");
        TH1F *hist_1 = (TH1F*)file1->Get("histProjX_1_Proton_allFSI");
        TH1F *hist1_1 = (TH1F*)file1->Get("histProjX_2_Proton_allFSI");
        TH1F *hist2 = (TH1F*)file->Get("histProjX_1_PiCh_allFSI");
        TH1F *hist3 = (TH1F*)file->Get("histProjX_2_PiCh_allFSI");
        TH1F *hist2_1 = (TH1F*)file1->Get("histProjX_1_PiCh_allFSI");
        TH1F *hist3_1 = (TH1F*)file1->Get("histProjX_2_PiCh_allFSI");
        
        std::vector<TH1*> projections;
        projections.push_back(hist);
        projections.push_back(hist1);
        std::vector<TH1*> projections1;
        projections1.push_back(hist_1);
        projections1.push_back(hist1_1);
        std::vector<TH1*> projections2;
        projections2.push_back(hist2);
        projections2.push_back(hist3);
        std::vector<TH1*> projections3;
        projections3.push_back(hist2_1);
        projections3.push_back(hist3_1);

        drawGraphs(gPad, projections, projections1, projections2, projections3, colors6, markerStyles6,rebinFactor);
        addText(0.90, 0.23, "(f)");
    }

    c1->SaveAs("K_and_P_hadd_rebin_new.png");
    c1->SaveAs("K_and_P_hadd_rebin_new.pdf");
    
    // 删除所有重分箱后的直方图
    for (TH1* h : rebinnedHistograms) {
        delete h;
    }
    rebinnedHistograms.clear();
    
    // 关闭文件
    file->Close();
    file1->Close();
}