#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>

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

std::vector<TH1D*> getCentPtProjections(const TString& filename, const TString& histname) {
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
    std::vector<TH1D*> projections;

    // 计算新的范围
   double rtRanges[8] = {0, 0.5, 1.0,1.5, 2.0,2.5,3.0,3.5};  // RT 的范围
   double rtUpperRanges[8] = {0.5, 1.0,1.5, 2.0,2.5, 3.0,3.5,30 / averageNT};  // RT 上限
   std::vector<double> avgPt; // 存储每个 RT 范围的平均 p_T
   std::vector<double> avgPtErr;
   for (int i = 0; i < 8; i++) {
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
        double meanPt = histProjX->GetMean();
        double meanPtErr = histProjX->GetMeanError();
        std::cout << "Mean pT for RT range [" << rtRanges[i] << ", " << rtUpperRanges[i] << "]: " << meanPt << std::endl;
        std::cout << "Mean pT Error for RT range [" << rtRanges[i] << ", " << rtUpperRanges[i] << "]: " << meanPtErr << std::endl;
        avgPt.push_back(meanPt);
        avgPtErr.push_back(meanPtErr);
    }
    return projections;
}

TGraphErrors* drawAveragePt(const std::vector<double>& avgPt, const std::vector<double>& avgPtErr, Color_t color, Style_t markerStyle) {
    float rtRanges[] = {0.25f, 0.75f, 1.25f, 1.75f, 2.25f, 2.75f, 3.25f, static_cast<float>(1.75 + 15.0 / averageNT)};

    int n = avgPt.size();
    float* x = new float[n];
    float* y = new float[n];
    float* ey = new float[n]; // 错误数组

    for (int i = 0; i < n; ++i) {
        x[i] = rtRanges[i];
        y[i] = static_cast<float>(avgPt[i]);
        ey[i] = static_cast<float>(avgPtErr[i]); // 设置误差值
    }

    TGraphErrors* gr = new TGraphErrors(n, x, y, 0, ey); // 使用误差数组
    gr->SetTitle("");
    gr->SetMarkerStyle(markerStyle);
    gr->SetMarkerColor(color);
    gr->SetLineColor(color);
    gr->SetLineWidth(2);

    return gr;
}

TGraphErrors* analyzeAndDraw(const std::vector<TH1D*>& projections, TVirtualPad* pad, Color_t color, Style_t markerStyle) {
    pad->cd();
    std::vector<double> avgPt;
    std::vector<double> avgPtErr;

    for (const auto& h : projections) {
        double meanPt = h->GetMean();
        double meanPtErr = h->GetMeanError();
        //std::cout << "Mean pT: " << meanPt << std::endl;
        //std::cout << "Mean pT Error: " << meanPtErr << std::endl;
        avgPt.push_back(meanPt);
        avgPtErr.push_back(meanPtErr);
    }

    if (!avgPt.empty()) {
        return drawAveragePt(avgPt, avgPtErr, color, markerStyle);
    }

    return nullptr;
}


void drawK() {
    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 500); // 创建Canvas

    TPaveText *text = new TPaveText(0.4, 0.01, 0.6, 0.1, "NDC");
    text->AddText("R_{T}");
    text->SetFillColor(0);
    text->SetTextAlign(22); // 文本居中对齐
    text->Draw("same");

    TLatex *label = new TLatex(0.01, 0.6, " <p_{T}>(GeV/c)");
    label->SetTextSize(0.06);
    label->SetTextAngle(90);
    label->SetTextAlign(23);
    label->Draw();

    c1->Divide(3, 1);

    c1->cd(1)->SetPad(0.04, 0.1, 0.39, 1); 
    c1->cd(1)->SetRightMargin(0);
    c1->cd(2)->SetPad(0.39, 0.1, 0.69, 1); 
    c1->cd(2)->SetRightMargin(0);
    c1->cd(2)->SetLeftMargin(0);
    c1->cd(3)->SetPad(0.69, 0.1, 0.99, 1);
    c1->cd(3)->SetRightMargin(0);
    c1->cd(3)->SetLeftMargin(0);

    TString filenames[3] = {"hist_outputallFSI.root", "hist_outputnoFSI.root", "hist_outputnohFSI.root"};
    TString histname[3] = {"hPiCh_dPhi0", "hPiCh_dPhi2", "hPiCh_dPhi1"};
    TString filenames1[3]= {"pi_toward.root", "pi_away.root", "pi_transverse.root"};

    std::vector<Color_t> colors = {kRed, kBlue, kMagenta};
    std::vector<Style_t> markerStyles = {20, 21, 22};

    gStyle->SetOptStat(kFALSE);
    for (int j = 0; j < 3; j++) {
    c1->cd(j+1);

    TH2D *axisFrame = new TH2D("axisFrame", "; ; ", 100, -0.1, 5.1, 100, 0.4, 0.9);
    axisFrame->Draw("axis");  // 只绘制坐标轴，不绘制直方图

    for (int i = 0; i < 3; i++) {
        Color_t color = colors[i];
        Style_t markerStyle = markerStyles[i];
        std::vector<TH1D*> projections = getCentPtProjections(filenames[i], histname[j]); // 获取投影
        TGraphErrors* graph = analyzeAndDraw(projections, gPad, color, markerStyle); // 分析并绘制每个投影
        TFile* infile = new TFile(filenames1[j], "READ");
        TGraphErrors *graph1 = (TGraphErrors*)infile->Get("Pt");
        graph1->SetMarkerStyle(24);
        graph1->SetMarkerSize(1);
        graph1->SetLineColor(kBlack);
        graph1->SetMarkerColor(graph1->GetLineColor());
        graph1->Draw("P"); 
        graph->Draw("LP same"); // 绘制到同一画布上
    }

    axisFrame->Draw("axis same");

    if (j == 0) {
        addText2(0.4, 0.8, "Toward");
    } else if (j == 1) {
        addText2(0.4, 0.8, "Away");
        addText4(0.33, 0.25, "#pi^{+}+#pi^{-}");
    } else if (j == 2) {
        addText2(0.4, 0.8, "Transverse");
    }
}

    c1->cd(1);
    TLegend *legend = new TLegend(0.65, 0.18, 1.1, 0.39);
    legend->SetNColumns(1);
    legend->SetFillStyle(0);

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

    c1->SaveAs("pi_RT_A.png");
    c1->SaveAs("Pi_RT_A.pdf");
}