#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TMath.h>
#include <vector>
#include <map>

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

TGraphErrors* drawConnectedPoints(TH1D* hist, Color_t color, Style_t markerStyle) {
    int n = hist->GetNbinsX();
    double* x = new double[n];
    double* y = new double[n];
    double* ey = new double[n]; // y误差数组

    int count = 0; // 用于记录满足条件的点的数量
    for (int i = 1; i <= n; ++i) {
        double binCenter = hist->GetBinCenter(i);
        if (binCenter<3) { // 只保留 x < 3 的点
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
void addText6(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);
    latex->SetTextSize(0.09);
    latex->Draw();
}
void addText5(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);
    latex->SetTextSize(0.07);
    latex->Draw();
}

std::vector<TH1D*> getCentPtProjections(const TString& filename, const TString& histname, const TString& suffix) {
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
    std::cout << "averageNT: " << averageNT << std::endl;
    std::vector<TH1D*> projections;
    double rtRanges[2] = {0,2.5};  // RT 的范围
    double rtUpperRanges[2] = {0.5,5};  // RT 上限
  
    for (int i = 0; i < 2; i++) {
        double ntRangeLow = rtRanges[i] * averageNT;  // RT 转为 NT 的下限
        double ntRangeUp = rtUpperRanges[i] * averageNT;  // RT 转为 NT 的上限

        std::cout << "ntRangeLow: " << ntRangeLow << std::endl;
        std::cout << "ntRangeUp: " << ntRangeUp << std::endl;

        int binLow = hist3D->GetZaxis()->FindBin(ntRangeLow);  // NT 下限对应的 bin
        int binUp = hist3D->GetZaxis()->FindBin(ntRangeUp);  // NT 上限对应的 bin;
        binUp-=1;

        std::cout << "RT Range[" << rtRanges[i] << ", " << rtUpperRanges[i] << "] mapped to NT Bins[" 
             << binLow << ", " << binUp << "]" << std::endl;

        TH1D* histProjX = hist3D->ProjectionX(Form("histProjX_%i_%s", i + 1, suffix.Data()), 1, hist3D->GetNbinsY(), binLow, binUp);
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

void drawGraphs(const std::vector<TH1D*>& projections1,const std::vector<TH1D*>& projections2, const std::vector<Color_t>& colors, const std::vector<Style_t>& markerStyles) {
  for (size_t i = 0; i < projections1.size(); ++i) {
    int rebinFactor=2;
    TH1* hist1 = rebinHistogramThreshold(projections1[i], 2.5, rebinFactor);
    TH1* hist2 = rebinHistogramThreshold(projections2[i], 2.5, rebinFactor);
    Color_t color = colors[i];
    Style_t markerStyle = markerStyles[i];

     TH1D* ratioHist = (TH1D*)hist1->Clone();
    ratioHist->Divide(hist2); // 除以hist4

    // 绘制相减和相除的结果
    TGraphErrors* graphRatio = drawConnectedPoints(ratioHist, color, markerStyle);

    if(i == 0) {
      graphRatio->Draw("LP");
    } else {
      graphRatio->Draw("LP same");
    }
  }
}


// 主绘图函数
void draw_tr() {
    TCanvas *c1 = new TCanvas("c1", "c1", 1290, 800);
    c1->SetTicks(1, 1);
    c1->SetMargin(0.17, 0.10, 0.2, 0.2);
    c1->Divide(3, 2, 0, 0);
    gStyle->SetOptStat(kFALSE);

    // 定义可重用变量
    std::vector<Color_t> colors = {kBlack, kRed, kBlue, kMagenta};
    std::vector<Style_t> markerStyles = {20, 21, 22, 23};

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
        
        TString filename = "hist_outputnoFSI_hadd.root";
        TString histname1 = "hKCh_dPhi1";
        TString histname2 = "hPiCh_dPhi1";
        gStyle->SetOptStat(kFALSE);

        std::vector<TH1D*> projections = getCentPtProjections(filename, histname1, "KCh_noFSI");
        std::vector<TH1D*> projections1 = getCentPtProjections(filename, histname2, "PiCh_noFSI");
        drawGraphs(projections, projections1, colors, markerStyles);
        
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

        addText6(0.23, 0.90, "Transverse");
        addText6(0.62, 0.90, "0mb w/o ART");
        addText5(0.92, 0.03, "(a)");
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
        TString filename = "hist_outputnohFSI_hadd.root";
        TString histname1 = "hKCh_dPhi1";
        TString histname2 = "hPiCh_dPhi1";
        gStyle->SetOptStat(kFALSE);
        
        std::vector<TH1D*> projections = getCentPtProjections(filename, histname1, "KCh_nohFSI");
        std::vector<TH1D*> projections1 = getCentPtProjections(filename, histname2, "PiCh_nohFSI");
        drawGraphs(projections, projections1, colors, markerStyles);

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

        addText6(0.45, 0.90, "0.15mb w/o ART");
        addText5(0.90, 0.03, "(c)");
    }

    // 第三列：K/π (allFSI)
    {
        c1->cd(3);
        gPad->SetTicks(1, 1);
        gPad->SetRightMargin(0.01);
        TH2D *axisFrame3 = new TH2D("axisFrame3", "; ; ", 100, -0.1, 2.99, 100, 0.01, 0.53);
        axisFrame3->GetYaxis()->SetTickLength(0.02);
        axisFrame3->GetXaxis()->SetTickLength(0.02);
        axisFrame3->GetXaxis()->SetLabelSize(0.065);
        axisFrame3->GetYaxis()->SetLabelSize(0.065);
        axisFrame3->Draw("axis");

        addText6(0.48, 0.90, "0.15mb w/ ART");

        TString filename = "hist_outputallFSI_more.root";
        TString histname1 = "hKCh_dPhi1";
        TString histname2 = "hPiCh_dPhi1";
        
        gStyle->SetOptStat(kFALSE);
        
        std::vector<TH1D*> projections = getCentPtProjections(filename, histname1, "KCh_allFSI");
        std::vector<TH1D*> projections1 = getCentPtProjections(filename, histname2, "PiCh_allFSI");
        drawGraphs(projections, projections1, colors, markerStyles);


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

        addText5(0.90, 0.03, "(e)");
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

        TString filename = "hist_outputnoFSI_hadd.root";
        TString histname1 = "hProton_dPhi1";
        TString histname2 = "hPiCh_dPhi1";
        
        std::vector<TH1D*> projections = getCentPtProjections(filename, histname1, "Proton_noFSI");
        std::vector<TH1D*> projections1 = getCentPtProjections(filename, histname2, "PiCh_noFSI");
        drawGraphs(projections, projections1, colors, markerStyles);

        TString filenames= "p_pi_ratio_transverse_05.root";
        TFile* infile1 = new TFile(filenames, "READ");
        TGraphErrors *graph1 = (TGraphErrors*)infile1->Get("Pt");
        graph1->SetMarkerStyle(24);
        graph1->SetMarkerSize(1.3);
        graph1->SetLineColor(kBlack);
        graph1->SetMarkerColor(graph1->GetLineColor());
        graph1->Draw("P"); 

        TLegend *legend2 = new TLegend(0.21, 0.68, 0.86, 0.89);
        legend2->SetNColumns(1);
        legend2->SetFillStyle(0);
        legend2->SetBorderSize(0);
        legend2->SetMargin(0.20);
        legend2->SetTextSize(0.068);

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
        addText5(0.26, 0.89, "AMPT");

        TString filenames1= "p_pi_ratio_transverse_255.root";
        TFile* infile2 = new TFile(filenames1, "READ");
        TGraphErrors *graph2 = (TGraphErrors*)infile2->Get("Pt");
        graph2->SetMarkerStyle(24);
        graph2->SetMarkerSize(1.3);
        graph2->SetLineColor(kRed);
        graph2->SetMarkerColor(graph2->GetLineColor());
        graph2->Draw("P same"); 
        addText2(0.92, 0.23, "(b)");
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

TLatex *text2 = new TLatex(0.53, 0.06, "p_{T} (GeV/c)");
text2->SetNDC();
text2->SetTextSize(0.08);
text2->SetTextAlign(22);
text2->Draw();

        TString filename = "hist_outputnohFSI_hadd.root";
        TString histname1 = "hProton_dPhi1";
        TString histname2 = "hPiCh_dPhi1";
        
        std::vector<TH1D*> projections = getCentPtProjections(filename, histname1, "Proton_nohFSI");
        std::vector<TH1D*> projections1 = getCentPtProjections(filename, histname2, "PiCh_nohFSI");
        drawGraphs(projections, projections1, colors, markerStyles);

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
        addText2(0.90, 0.23, "(d)");

        addText5(0.09, 0.89, "ALICE");

        TLegend *legend2 = new TLegend(0.05, 0.68, 0.7, 0.89);
        legend2->SetNColumns(1);
        legend2->SetFillStyle(0);
        legend2->SetBorderSize(0);
        legend2->SetMargin(0.20);
        legend2->SetTextSize(0.068);

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

        TString filename = "hist_outputallFSI_more.root";
        TString histname1 = "hProton_dPhi1";
        TString histname2 = "hPiCh_dPhi1";
        
        std::vector<TH1D*> projections = getCentPtProjections(filename, histname1, "Proton_allFSI");
        std::vector<TH1D*> projections1 = getCentPtProjections(filename, histname2, "PiCh_allFSI");
        drawGraphs(projections, projections1, colors, markerStyles);

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
        addText2(0.90, 0.23, "(f)");
    }

    c1->SaveAs("pi_k_p_transverse.png");
    c1->SaveAs("pi_k_p_transverse.pdf");
}