#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TH2D.h>
#include <TFile.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TDirectory.h>

void addText(Double_t x, Double_t y, const char* text, Double_t size) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);
    latex->SetTextSize(size);
    latex->Draw();
}

void draw() {
    TCanvas *c1 = new TCanvas("c1", "c1", 1500, 1200);
    c1->SetTicks(1, 1);
    c1->SetMargin(0.17, 0.10, 0.2, 0.2);
    c1->Divide(3, 3, 0, 0); // 不增加分割线间距
    gStyle->SetOptStat(kFALSE);

    // 打开三个粒子的结果文件
    TFile *piFile = new TFile("pi_results.root", "READ");
    TFile *kFile = new TFile("k_results.root", "READ");
    TFile *protonFile = new TFile("p_results.root", "READ");
    
    std::vector<Color_t> colors = {kRed, kBlue, kMagenta};
    std::vector<Style_t> markerStyles = {20, 21, 22};
    
    // 第一子图: Toward区域 - π介子
    c1->cd(1);
    {
        gPad->SetTicks(1,1);
        TH2D *axisFrame1 = new TH2D("axisFrame1", "", 100, -0.1, 5.1, 100, 0.28, 2.19);
        axisFrame1->GetYaxis()->SetTickLength(0.02);
        axisFrame1->GetXaxis()->SetTickLength(0.02);
        axisFrame1->GetXaxis()->SetLabelSize(0.085);
        axisFrame1->GetYaxis()->SetLabelSize(0.085);
        axisFrame1->Draw("axis");

        // 绘制π介子Toward区域数据
        piFile->cd("Toward");
        TGraphErrors *piToward0 = (TGraphErrors*)gDirectory->Get("Pi_Toward_0");
        TGraphErrors *piToward1 = (TGraphErrors*)gDirectory->Get("Pi_Toward_1");
        TGraphErrors *piToward2 = (TGraphErrors*)gDirectory->Get("Pi_Toward_2");
        TGraphErrors *piToward3 = (TGraphErrors*)gDirectory->Get("Pi_Toward_4");
        
        if (piToward0) piToward0->Draw("LP");
        if (piToward1) piToward1->Draw("P same");
        if (piToward2) piToward2->Draw("LP same");
        if (piToward3) piToward3->Draw("LP same");
        
        addText(0.49, 0.9, "Toward", 0.10);
    }

    // 第二子图: Transverse区域 - π介子
    c1->cd(2);
    {
        gPad->SetTicks(1,1);
        TH2D *axisFrame2 = new TH2D("axisFrame2", "", 100, -0.1, 5.1, 50,0.28, 2.19);
        axisFrame2->GetYaxis()->SetTickLength(0.02);
        axisFrame2->GetXaxis()->SetTickLength(0.02);
        axisFrame2->GetXaxis()->SetLabelSize(0.085);
        axisFrame2->GetYaxis()->SetLabelSize(0.085);
        axisFrame2->Draw("axis");

        // 绘制π介子Transverse区域数据
        piFile->cd("Transverse");
        TGraphErrors *piTransverse0 = (TGraphErrors*)gDirectory->Get("Pi_Transverse_0");
        TGraphErrors *piTransverse1 = (TGraphErrors*)gDirectory->Get("Pi_Transverse_1");
        TGraphErrors *piTransverse2 = (TGraphErrors*)gDirectory->Get("Pi_Transverse_2");
        TGraphErrors *piTransverse3 = (TGraphErrors*)gDirectory->Get("Pi_Transverse_4");
        
        if (piTransverse0) piTransverse0->Draw("LP");
        if (piTransverse1) piTransverse1->Draw("P same");
        if (piTransverse2) piTransverse2->Draw("LP same");
        if (piTransverse3) piTransverse3->Draw("LP same");
        
        addText(0.33, 0.9, "Transverse", 0.10);
        addText(0.45, 0.05, "#pi^{+} + #pi^{-}", 0.15);
    }

    // 第三子图: In-Jet区域 - π介子
    c1->cd(3);
    {
        gPad->SetTicks(1,1);
        gPad->SetRightMargin(0.01);
        TH2D *axisFrame3 = new TH2D("axisFrame3", "", 100, -0.1,1.3, 50,0.28,  2.19);
        axisFrame3->GetYaxis()->SetTickLength(0.02);
        axisFrame3->GetXaxis()->SetTickLength(0.02);
        axisFrame3->GetXaxis()->SetLabelSize(0.085);
        axisFrame3->GetYaxis()->SetLabelSize(0.085);
        axisFrame3->Draw("axis");

        // 绘制π介子In-Jet区域数据
        piFile->cd("InJet");
        TGraphErrors *piInJet0 = (TGraphErrors*)gDirectory->Get("Pi_InJet_0");
        TGraphErrors *piInJet1 = (TGraphErrors*)gDirectory->Get("Pi_InJet_1");
        TGraphErrors *piInJet2 = (TGraphErrors*)gDirectory->Get("Pi_InJet_2");
        TGraphErrors *piInJet3 = (TGraphErrors*)gDirectory->Get("Pi_InJet_4");
        
        if (piInJet0) piInJet0->Draw("LP");
        if (piInJet1) piInJet1->Draw("LP same");
        if (piInJet2) piInJet2->Draw("LP same");
        if (piInJet3) piInJet3->Draw("LP same");
        
        addText(0.42, 0.9, "In-Jet", 0.10);
    }

    // 第四子图: Toward区域 - K介子
    c1->cd(4);
    {
        gPad->SetTicks(1,1);
        TH2D *axisFrame4 = new TH2D("axisFrame4", "", 100, -0.1, 5.1, 50,0.28,  2.19);
        axisFrame4->GetYaxis()->SetTickLength(0.02);
        axisFrame4->GetXaxis()->SetTickLength(0.02);
        axisFrame4->GetXaxis()->SetLabelSize(0.085);
        axisFrame4->GetYaxis()->SetLabelSize(0.085);
        axisFrame4->Draw("axis");

        // 绘制K介子Toward区域数据
        kFile->cd("Toward");
        TGraphErrors *kToward0 = (TGraphErrors*)gDirectory->Get("K_Toward_0");
        // TGraphErrors *kToward1 = (TGraphErrors*)gDirectory->Get("K_Toward_1");
        TGraphErrors *kToward2 = (TGraphErrors*)gDirectory->Get("K_Toward_2");
        TGraphErrors *kToward3 = (TGraphErrors*)gDirectory->Get("K_Toward_4");
        
        if (kToward0) kToward0->Draw("LP");
        // if (kToward1) kToward1->Draw("P same");
        if (kToward2) kToward2->Draw("LP same");
        if (kToward3) kToward3->Draw("LP same");
        
        TLatex* label1 = new TLatex();
        label1->SetTextSize(0.1);
        label1->SetTextAngle(90);
        label1->SetTextAlign(22);
        label1->DrawLatexNDC(0.05, 0.54, "<p_{T}>(GeV/c)");
    }

    // 第五子图: Transverse区域 - K介子
    c1->cd(5);
    {
        gPad->SetTicks(1,1);
        TH2D *axisFrame5 = new TH2D("axisFrame5", "", 100, -0.1, 5.1, 50,0.28,  2.19);
        axisFrame5->GetYaxis()->SetTickLength(0.02);
        axisFrame5->GetXaxis()->SetTickLength(0.02);
        axisFrame5->GetXaxis()->SetLabelSize(0.085);
        axisFrame5->GetYaxis()->SetLabelSize(0.085);
        axisFrame5->Draw("axis");

        // 绘制K介子Transverse区域数据
        kFile->cd("Transverse");
        TGraphErrors *kTransverse0 = (TGraphErrors*)gDirectory->Get("K_Transverse_0");
        // TGraphErrors *kTransverse1 = (TGraphErrors*)gDirectory->Get("K_Transverse_1");
        TGraphErrors *kTransverse2 = (TGraphErrors*)gDirectory->Get("K_Transverse_2");
        TGraphErrors *kTransverse3 = (TGraphErrors*)gDirectory->Get("K_Transverse_4");
        
        if (kTransverse0) kTransverse0->Draw("LP");
        // if (kTransverse1) kTransverse1->Draw("P same");
        if (kTransverse2) kTransverse2->Draw("LP same");
        if (kTransverse3) kTransverse3->Draw("LP same");
        
        addText(0.45, 0.08, "K^{+} + K^{-}", 0.15);
    }

    // 第六子图: In-Jet区域 - K介子
    c1->cd(6);
    {
        gPad->SetTicks(1,1);
        gPad->SetRightMargin(0.01);
        TH2D *axisFrame6 = new TH2D("axisFrame6", "", 100, -0.1, 1.3, 50, 0.28,  2.19);
        axisFrame6->GetYaxis()->SetTickLength(0.02);
        axisFrame6->GetXaxis()->SetTickLength(0.02);
        axisFrame6->GetXaxis()->SetLabelSize(0.085);
        axisFrame6->GetYaxis()->SetLabelSize(0.085);
        axisFrame6->Draw("axis");

        // 绘制K介子In-Jet区域数据
        kFile->cd("InJet");
        TGraphErrors *kInJet0 = (TGraphErrors*)gDirectory->Get("K_InJet_0");
        TGraphErrors *kInJet1 = (TGraphErrors*)gDirectory->Get("K_InJet_1");
        TGraphErrors *kInJet2 = (TGraphErrors*)gDirectory->Get("K_InJet_2");
        TGraphErrors *kInJet3 = (TGraphErrors*)gDirectory->Get("K_InJet_4");
        
        if (kInJet0) kInJet0->Draw("LP");
        if (kInJet1) kInJet1->Draw("LP same");
        if (kInJet2) kInJet2->Draw("LP same");
        if (kInJet3) kInJet3->Draw("LP same");
    }

    // 第七子图: Toward区域 - 质子
    c1->cd(7);
    {
        gPad->SetTicks(1,1);
        TH2D *axisFrame7 = new TH2D("axisFrame7", "", 100, -0.1, 5.1,50,0.28,  2.19);
        axisFrame7->GetYaxis()->SetTickLength(0.02);
        axisFrame7->GetXaxis()->SetTickLength(0.02);
        axisFrame7->GetXaxis()->SetLabelSize(0.075);
        axisFrame7->GetYaxis()->SetLabelSize(0.075);
        axisFrame7->Draw("axis");

        // 绘制质子Toward区域数据
        protonFile->cd("Toward");
        TGraphErrors *protonToward0 = (TGraphErrors*)gDirectory->Get("Proton_Toward_0");
        // TGraphErrors *protonToward1 = (TGraphErrors*)gDirectory->Get("Proton_Toward_1");
        TGraphErrors *protonToward2 = (TGraphErrors*)gDirectory->Get("Proton_Toward_2");
        TGraphErrors *protonToward3 = (TGraphErrors*)gDirectory->Get("Proton_Toward_4");
        
        if (protonToward0) protonToward0->Draw("LP");
        // if (protonToward1) protonToward1->Draw("P same");
        if (protonToward2) protonToward2->Draw("LP same");
        if (protonToward3) protonToward3->Draw("LP same");
    }

    // 第八子图: Transverse区域 - 质子
    c1->cd(8);
    {
        gPad->SetTicks(1,1);
        TH2D *axisFrame8 = new TH2D("axisFrame8", "", 100, -0.1, 5.1,50,0.28, 2.19);
        axisFrame8->GetYaxis()->SetTickLength(0.02);
        axisFrame8->GetXaxis()->SetTickLength(0.02);
        axisFrame8->GetXaxis()->SetLabelSize(0.075);
        axisFrame8->GetYaxis()->SetLabelSize(0.075);
        axisFrame8->Draw("axis");

        // 绘制质子Transverse区域数据
        protonFile->cd("Transverse");
        TGraphErrors *protonTransverse0 = (TGraphErrors*)gDirectory->Get("Proton_Transverse_0");
        // TGraphErrors *protonTransverse1 = (TGraphErrors*)gDirectory->Get("Proton_Transverse_1");
        TGraphErrors *protonTransverse2 = (TGraphErrors*)gDirectory->Get("Proton_Transverse_2");
        TGraphErrors *protonTransverse3 = (TGraphErrors*)gDirectory->Get("Proton_Transverse_4");
        
        if (protonTransverse0) protonTransverse0->Draw("LP");
        // if (protonTransverse1) protonTransverse1->Draw("P same");
        if (protonTransverse2) protonTransverse2->Draw("LP same");
        if (protonTransverse3) protonTransverse3->Draw("LP same");
        
TLatex *text2 = new TLatex(0.53, 0.06, "R_{T}");
text2->SetNDC();
text2->SetTextSize(0.12);
text2->SetTextAlign(22);
text2->Draw();
        
        addText(0.45, 0.27, "p + #bar{p}", 0.15);
    }

    // 第九子图: In-Jet区域 - 质子
    c1->cd(9);
    {
        gPad->SetTicks(1,1);
        gPad->SetRightMargin(0.01);
        TH2D *axisFrame9 = new TH2D("axisFrame9", "", 100, -0.1, 1.3, 50, 0.28,  2.19);
        axisFrame9->GetYaxis()->SetTickLength(0.02);
        axisFrame9->GetXaxis()->SetTickLength(0.02);
        axisFrame9->GetXaxis()->SetLabelSize(0.075);
        axisFrame9->GetYaxis()->SetLabelSize(0.075);
        axisFrame9->Draw("axis");

        // 绘制质子In-Jet区域数据
        protonFile->cd("InJet");
        TGraphErrors *protonInJet0 = (TGraphErrors*)gDirectory->Get("Proton_InJet_0");
        TGraphErrors *protonInJet1 = (TGraphErrors*)gDirectory->Get("Proton_InJet_1");
        TGraphErrors *protonInJet2 = (TGraphErrors*)gDirectory->Get("Proton_InJet_2");
        TGraphErrors *protonInJet3 = (TGraphErrors*)gDirectory->Get("Proton_InJet_4");
        
        if (protonInJet0) protonInJet0->Draw("LP");
        if (protonInJet1) protonInJet1->Draw("LP same");
        if (protonInJet2) protonInJet2->Draw("LP same");
        if (protonInJet3) protonInJet3->Draw("LP same");
    }

    // 添加图例到第一个pad
    c1->cd(1);
    TLegend *legend = new TLegend(0.2, 0.52, 0.70, 0.82);
    legend->SetNColumns(1);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.08); 

    // TGraphErrors *marker1 = new TGraphErrors();
    // marker1->SetMarkerStyle(24);
    // marker1->SetMarkerColor(kBlack);
    // marker1->SetMarkerSize(1);
    // marker1->SetLineWidth(2); 
    // marker1->SetLineColor(kBlack);
    // legend->AddEntry(marker1, "ALICE", "lp");

    TGraphErrors *marker2 = new TGraphErrors();
    marker2->SetMarkerStyle(20);
    marker2->SetMarkerColor(kRed);
    marker2->SetMarkerSize(1);
    marker2->SetLineWidth(2); 
    marker2->SetLineColor(kRed);
    legend->AddEntry(marker2, "nohFSI 0.10mb", "lp");

    TGraphErrors *marker3 = new TGraphErrors();
    marker3->SetMarkerStyle(21);
    marker3->SetMarkerColor(kBlue);
    marker3->SetMarkerSize(1);
    marker3->SetLineWidth(2); 
    marker3->SetLineColor(kBlue);
    legend->AddEntry(marker3, "nohFSI 0.15mb", "lp");

    TGraphErrors *marker4 = new TGraphErrors();
    marker4->SetMarkerStyle(22);
    marker4->SetMarkerColor(kMagenta);
    marker4->SetMarkerSize(1);
    marker4->SetLineWidth(2); 
    marker4->SetLineColor(kMagenta);
    legend->AddEntry(marker4, "nohFSI 0.20mb", "lp");
    legend->Draw(); 

    c1->Update();
    c1->SaveAs("pi_k_p_05_mb.png");
    c1->SaveAs("pi_k_p_05_mb.pdf");

    // 清理资源
    delete piFile;
    delete kFile;
    delete protonFile;
}