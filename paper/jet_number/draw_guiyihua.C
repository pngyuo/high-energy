#include <iostream>
#include <vector>
#include <string>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TString.h>
#include <TPaveText.h>

void draw_guiyihua(){ // 函数名已修改
    // ========== 1. Setup Global Plotting Style ==========
    gStyle->SetOptStat(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(kWhite);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameBorderSize(1);
    gStyle->SetFrameFillColor(0);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetAxisColor(1, "XYZ");
    gStyle->SetStripDecimals(kTRUE);
    gStyle->SetTickLength(0.03, "XYZ");
    gStyle->SetNdivisions(510, "XYZ");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetLegendFont(42);
    gStyle->SetMarkerSize(1.2);
    gStyle->SetLineWidth(2);

    // ========== 2. Open File and Get Objects ==========
    TFile *f = TFile::Open("hist_outputallFSI_more01jet.root");
    if (!f || f->IsZombie()){
        std::cout << "Error: Could not open file 'hist_outputallFSI_liangjet.root'!\n";
        return;
    }

    TH1D *hNchTR_hist = (TH1D*)f->Get("nTHist");
    if (!hNchTR_hist){
        std::cout << "Error: Histogram 'nTHist' not found!\n";
        f->Close();
        return;
    }
    double avgNTR = hNchTR_hist->GetMean();
    std::cout << "Calculated average charged particle multiplicity in Transverse Region: <N_ch^TR> = " << avgNTR << std::endl;

    TTree *t = (TTree*)f->Get("tTRjet");
    if (!t){
        std::cout << "Error: TTree 'tTRjet' not found!\n";
        f->Close();
        return;
    }

    // ========== 3. Define Bins and Initialize Arrays ==========
    const int nRtBins = 20;
    double rtEdges[nRtBins + 1];
    for (int i = 0; i <= nRtBins; ++i) rtEdges[i] = 0.0 + i * 0.6;

    double rtCenter[nRtBins];
    double avgNjTR[nRtBins] = {0}, avgMultTR[nRtBins] = {0};
    double avgNjTow[nRtBins] = {0}, avgMultTow[nRtBins] = {0};

    for (int i = 0; i < nRtBins; ++i){
        rtCenter[i] = (rtEdges[i] + rtEdges[i+1]) / 2.0;
    }

    // ========== 4. Setup TTree Branches ==========
    int NchTR, NjetTR, NchTow, NjetTow;
    std::vector<double> *tJetPt = nullptr, *tJetTowPt = nullptr;
    std::vector<int> *tJetConstN = nullptr, *tJetTowConstN = nullptr;

    t->SetBranchAddress("NchTR", &NchTR);
    t->SetBranchAddress("NjetTR", &NjetTR);
    t->SetBranchAddress("jetPt", &tJetPt);
    t->SetBranchAddress("tJetConstN", &tJetConstN);

    t->SetBranchAddress("NchTow", &NchTow);
    t->SetBranchAddress("NjetTow", &NjetTow);
    t->SetBranchAddress("jetTowPt", &tJetTowPt);
    t->SetBranchAddress("tJetTowConstN", &tJetTowConstN);

    // ========== 5. Create Histograms for Distributions ==========
    TH1D *hJetNumDistTR[nRtBins], *hJetMultDistTR[nRtBins];
    TH1D *hJetNumDistTow[nRtBins], *hJetMultDistTow[nRtBins];
    TH1D *hJetPtDistTR[nRtBins], *hJetPtDistTow[nRtBins];

    for (int i = 0; i < nRtBins; ++i) {
        hJetNumDistTR[i]  = new TH1D(Form("hJetNumDistTR_%d", i), "", 20, 0, 20);
        hJetMultDistTR[i] = new TH1D(Form("hJetMultDistTR_%d", i), "", 50, 0, 50);
        hJetNumDistTow[i]  = new TH1D(Form("hJetNumDistTow_%d", i), "", 20, 0, 20);
        hJetMultDistTow[i] = new TH1D(Form("hJetMultDistTow_%d", i), "", 50, 0, 50);
        
        // 确保 pT 谱的 bin 设置与原始代码中的一致
        hJetPtDistTR[i] = new TH1D(Form("hJetPtDistTR_%d", i), ";p_{T} (GeV/c);Normalized to Unity", 100, 0, 5);
        hJetPtDistTow[i] = new TH1D(Form("hJetPtDistTow_%d", i), ";p_{T} (GeV/c);Normalized to Unity", 100, 0, 5);
    }

    // ========== 6. Loop Over TTree and Fill Data (不变) ==========
    std::vector<double> sumNjTR(nRtBins, 0.0), sumTotConstTR(nRtBins, 0.0), sumEvt(nRtBins, 0.0);
    std::vector<double> sumNjTow(nRtBins, 0.0), sumTotConstTow(nRtBins, 0.0);

    Long64_t nEntries = t->GetEntries();
    std::cout << "Processing " << nEntries << " events..." << std::endl;
    for (Long64_t i=0; i<nEntries; ++i) {
        t->GetEntry(i);
        if (i % 10000 == 0) std::cout << "  Event: " << i << "/" << nEntries << "\r" << std::flush;

        double rt = (avgNTR > 0) ? NchTR / avgNTR : 0;
        
        for (int j = 0; j < nRtBins; ++j){
            if (rt >= rtEdges[j] && rt < rtEdges[j+1]){
                sumEvt[j]++;
                
                // --- Process TR Region ---
                sumNjTR[j] += NjetTR;
                hJetNumDistTR[j]->Fill(NjetTR);
                if (tJetConstN) {
                    for (int nConst : *tJetConstN) {
                        sumTotConstTR[j] += nConst;
                        hJetMultDistTR[j]->Fill(nConst);
                    }
                }
                if (tJetPt) {
                    for (double pt : *tJetPt) {
                        hJetPtDistTR[j]->Fill(pt);
                    }
                }
                
                // --- Process Toward Region ---
                sumNjTow[j] += NjetTow;
                hJetNumDistTow[j]->Fill(NjetTow);
                if (tJetTowConstN){
                    for (int nConst : *tJetTowConstN) {
                        sumTotConstTow[j] += nConst;
                        hJetMultDistTow[j]->Fill(nConst);
                    }
                }
                if (tJetTowPt) {
                    for (double pt : *tJetTowPt) {
                        hJetPtDistTow[j]->Fill(pt);
                    }
                }
                break;
            }
        }
    }
    std::cout << "\nEvent processing complete." << std::endl;

    // ========== 7. Calculate Final Averages (不变) ==========
    for (int i = 0; i < nRtBins; ++i){
        if (sumEvt[i] > 0){
            avgNjTR[i]   = sumNjTR[i] / sumEvt[i];
            avgMultTR[i] = (sumNjTR[i] > 0) ? sumTotConstTR[i] / sumNjTR[i] : 0;
            avgNjTow[i]   = sumNjTow[i] / sumEvt[i];
            avgMultTow[i] = (sumNjTow[i] > 0) ? sumTotConstTow[i] / sumNjTow[i] : 0;
        }
    }

    // ========== 8. Normalize and Create Main Comparison Plot (不变) ==========
    double maxNjTR   = TMath::MaxElement(nRtBins, avgNjTR);
    double maxMultTR = TMath::MaxElement(nRtBins, avgMultTR);
    double maxNjTow   = TMath::MaxElement(nRtBins, avgNjTow);
    double maxMultTow = TMath::MaxElement(nRtBins, avgMultTow);
    
    double avgNjTRNorm[nRtBins], avgMultTRNorm[nRtBins], avgNjTowNorm[nRtBins], avgMultTowNorm[nRtBins];
    for(int i = 0; i < nRtBins; ++i) {
        avgNjTRNorm[i]   = (maxNjTR > 0) ? avgNjTR[i] / maxNjTR : 0;
        avgMultTRNorm[i] = (maxMultTR > 0) ? avgMultTR[i] / maxMultTR : 0;
        avgNjTowNorm[i]   = (maxNjTow > 0) ? avgNjTow[i] / maxNjTow : 0;
        avgMultTowNorm[i] = (maxMultTow > 0) ? avgMultTow[i] / maxMultTow : 0;
    }
    
    TCanvas *c1 = new TCanvas("c1", "Normalized Jet Properties vs Rt", 900, 700);
    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Normalized Jet Properties vs. R_{T}");

    TGraph *gNjTRNorm   = new TGraph(nRtBins, rtCenter, avgNjTRNorm);
    gNjTRNorm->SetMarkerStyle(kFullCircle);
    gNjTRNorm->SetMarkerColor(kRed+1);
    gNjTRNorm->SetLineColor(kRed+1);
    gNjTRNorm->SetTitle("N_{jet}^{TR}");
    
    TGraph *gMultTRNorm = new TGraph(nRtBins, rtCenter, avgMultTRNorm);
    gMultTRNorm->SetMarkerStyle(kFullSquare);
    gMultTRNorm->SetMarkerColor(kAzure+2);
    gMultTRNorm->SetLineColor(kAzure+2);
    gMultTRNorm->SetTitle("n_{const}^{TR}");

    TGraph *gNjTowNorm   = new TGraph(nRtBins, rtCenter, avgNjTowNorm);
    gNjTowNorm->SetMarkerStyle(kOpenCircle);
    gNjTowNorm->SetMarkerColor(kRed+1);
    gNjTowNorm->SetLineColor(kRed+1);
    gNjTowNorm->SetLineStyle(kDashed);
    gNjTowNorm->SetTitle("N_{jet}^{Toward}");

    TGraph *gMultTowNorm = new TGraph(nRtBins, rtCenter, avgMultTowNorm);
    gMultTowNorm->SetMarkerStyle(kOpenSquare);
    gMultTowNorm->SetMarkerColor(kAzure+2);
    gMultTowNorm->SetLineColor(kAzure+2);
    gMultTowNorm->SetLineStyle(kDashed);
    gMultTowNorm->SetTitle("n_{const}^{Toward}");

    mg->Add(gNjTRNorm, "PL");
    mg->Add(gMultTRNorm, "PL");
    mg->Add(gNjTowNorm, "PL");
    mg->Add(gMultTowNorm, "PL");
    mg->Draw("A");
    mg->GetXaxis()->SetTitle("R_{T} = N_{ch}^{TR} / #LTN_{ch}^{TR}#GT");
    mg->GetYaxis()->SetTitle("Normalized Yield");
    mg->GetXaxis()->SetTitleSize(0.045);
    mg->GetXaxis()->SetRangeUser(0, 5);
    mg->GetYaxis()->SetTitleSize(0.045);
    mg->GetYaxis()->SetRangeUser(0, 1.4);

    TLegend *legend = new TLegend(0.55, 0.70, 0.90, 0.90);
    legend->AddEntry(gNjTRNorm, "#LTN_{jet}#GT_{TR} (Solid)", "pl");
    legend->AddEntry(gMultTRNorm, "#LTn_{const}#GT_{TR} (Solid)", "pl");
    legend->AddEntry(gNjTowNorm, "#LTN_{jet}#GT_{Toward} (Dashed)", "pl");
    legend->AddEntry(gMultTowNorm, "#LTn_{const}#GT_{Toward} (Dashed)", "pl");
    legend->Draw();
    
    c1->SaveAs("normalized_comparison_with_toward_more01.png");

    // ========== 9. Plot Overlaid Distributions for Representative Bins ==========
    int selectedBins[] = {0, 2, 4, 6};
    
    // --- Jet Number Distributions ---
    TCanvas *c2 = new TCanvas("c2", "Jet Number Distributions", 1200, 900);
    c2->Divide(2, 2);
    for (int i = 0; i < 4; ++i) {
        c2->cd(i + 1);
        gPad->SetLeftMargin(0.15);
        int binIdx = selectedBins[i];
        TH1D* hTR = hJetNumDistTR[binIdx];
        TH1D* hTow = hJetNumDistTow[binIdx];

        if (hTR->GetEntries() > 0) hTR->Scale(1.0 / hTR->Integral());
        if (hTow->GetEntries() > 0) hTow->Scale(1.0 / hTow->Integral());

        hTR->SetLineColor(kAzure+2);
        hTR->SetFillColorAlpha(kAzure-9, 0.5);
        hTow->SetLineColor(kRed+1);
        hTow->SetFillStyle(0);

        double maxY = TMath::Max(hTR->GetMaximum(), hTow->GetMaximum()) * 1.2;
        hTR->SetMaximum(maxY);
        hTR->GetYaxis()->SetTitle("Normalized to Unity");
        hTR->GetXaxis()->SetTitle("N_{jet}");
        hTR->Draw("HIST");
        hTow->Draw("HIST SAME");

        TPaveText *pt = new TPaveText(0.4, 0.8, 0.9, 0.9, "NDC");
        pt->SetFillColor(0);
        pt->SetBorderSize(0);
        pt->AddText(Form("%.1f #leq R_{T} < %.1f", rtEdges[binIdx], rtEdges[binIdx+1]));
        pt->Draw();
        if (i == 0) {
            TLegend* leg_dist = new TLegend(0.5, 0.6, 0.9, 0.75);
            leg_dist->AddEntry(hTR, "Transverse", "lf");
            leg_dist->AddEntry(hTow, "Toward", "l");
            leg_dist->Draw();
        }
    }
    c2->SaveAs("jet_number_distributions_overlay_more01.png");

    // --- Jet Multiplicity Distributions ---
    TCanvas *c3 = new TCanvas("c3", "Jet Constituent Multiplicity Distributions", 1200, 900);
    c3->Divide(2, 2);
    for (int i = 0; i < 4; ++i) {
        c3->cd(i + 1);
        gPad->SetLeftMargin(0.15);
        int binIdx = selectedBins[i];
        TH1D* hTR = hJetMultDistTR[binIdx];
        TH1D* hTow = hJetMultDistTow[binIdx];

        if (hTR->GetEntries() > 0) hTR->Scale(1.0 / hTR->Integral());
        if (hTow->GetEntries() > 0) hTow->Scale(1.0 / hTow->Integral());

        hTR->SetLineColor(kAzure+2);
        hTR->SetFillColorAlpha(kAzure-9, 0.5);
        hTow->SetLineColor(kRed+1);
        hTow->SetFillStyle(0);
        
        double maxY = TMath::Max(hTR->GetMaximum(), hTow->GetMaximum()) * 1.2;
        hTR->SetMaximum(maxY);
        hTR->GetYaxis()->SetTitle("Normalized to Unity");
        hTR->GetXaxis()->SetTitle("n_{const}");
        hTR->GetXaxis()->SetRangeUser(0, 20);
        hTR->Draw("HIST");
        hTow->Draw("HIST SAME");

        TPaveText *pt = new TPaveText(0.4, 0.8, 0.9, 0.9, "NDC");
        pt->SetFillColor(0);
        pt->SetBorderSize(0);
        pt->AddText(Form("%.1f #leq R_{T} < %.1f", rtEdges[binIdx], rtEdges[binIdx+1]));
        pt->Draw();
        if (i == 0) {
            TLegend* leg_dist = new TLegend(0.5, 0.6, 0.9, 0.75);
            leg_dist->AddEntry(hTR, "Transverse", "lf");
            leg_dist->AddEntry(hTow, "Toward", "l");
            leg_dist->Draw();
        }
    }
    c3->SaveAs("jet_multiplicity_distributions_overlay_more01.png");

 // ========== 10. Plot Jet pT Distributions with Log Y-axis ==========
    TCanvas *c4 = new TCanvas("c4", "Jet pT Distributions", 1200, 900);
    c4->Divide(2, 2);
    const double pt_xmin = 0.0;
    const double pt_xmax = 5.0;

    for (int i = 0; i < 4; ++i) {
        c4->cd(i + 1);
        gPad->SetLogy(1);
        gPad->SetLeftMargin(0.15);
        int binIdx = selectedBins[i];
        TH1D* hTR = hJetPtDistTR[binIdx];
        TH1D* hTow = hJetPtDistTow[binIdx];
        
        double N_events = sumEvt[binIdx];

        if (hTR->GetEntries() > 0 && N_events > 0) {
            hTR->Scale(1.0 / N_events , "width");
        }
        if (hTow->GetEntries() > 0 && N_events > 0) {
            hTow->Scale(1.0 / N_events , "width");
        }

        hTR->SetLineColor(kAzure+2);
        // hTR->SetFillColorAlpha(kAzure-9, 0.5); // <-- 移除填充
        hTow->SetLineColor(kRed+1);
        hTow->SetFillStyle(0);
        
        double maxY = TMath::Max(hTR->GetMaximum(), hTow->GetMaximum()) * 5.0; 
        if(maxY == 0) maxY = 1;
        hTR->SetMaximum(maxY);
        hTR->SetMinimum(1e-5); 
        
        hTR->GetYaxis()->SetTitle("d^{2}N_{jet} / (N_{evt} dp_{T} d#eta)");
        hTR->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        
        hTR->Draw("L");
        hTow->Draw("L SAME");
        hTR->GetXaxis()->SetRangeUser(pt_xmin, pt_xmax);

        TPaveText *pt = new TPaveText(0.4, 0.8, 0.9, 0.9, "NDC");
        pt->SetFillColor(0);
        pt->SetBorderSize(0);
        pt->AddText(Form("%.1f #leq R_{T} < %.1f", rtEdges[binIdx], rtEdges[binIdx+1]));
        pt->Draw();

        if (i == 0) {
            TLegend* leg_dist_pt = new TLegend(0.18, 0.18, 0.5, 0.33); 
            leg_dist_pt->AddEntry(hTR, "Transverse", "l"); 
            leg_dist_pt->AddEntry(hTow, "Toward", "l");  
            leg_dist_pt->Draw();
        }
    }
    c4->SaveAs("jet_pt_distributions_overlay_logy_normalized_to_event.png"); // 更新文件名
    f->Close();
}