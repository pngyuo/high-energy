#include <iostream>
#include <vector>
#include <string>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TString.h>
#include <TPaveText.h>
#include <TLatex.h>

// Helper function to convert a TH1D to a TGraphErrors.
// Horizontal errors are set to 0.
TGraphErrors* convertToGraph(TH1D* hist, Color_t color, Style_t markerStyle) {
    if (!hist) return nullptr;
    int n = hist->GetNbinsX();
    std::vector<double> x, y, ex, ey;

    for (int i = 1; i <= n; ++i) {
        if (hist->GetBinContent(i) != 0 || hist->GetBinError(i) != 0) {
            x.push_back(hist->GetBinCenter(i));
            y.push_back(hist->GetBinContent(i));
            ex.push_back(0); // --- Set horizontal error to zero ---
            ey.push_back(hist->GetBinError(i));
        }
    }

    if (x.empty()) return new TGraphErrors();

    TGraphErrors* graph = new TGraphErrors(x.size(), &x[0], &y[0], &ex[0], &ey[0]);
    graph->SetLineColor(color);
    graph->SetLineWidth(2);
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerColor(color);
    graph->SetMarkerSize(1.0);
    return graph;
}


void draw_guiyihua(){
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
    gStyle->SetMarkerSize(1.0);
    gStyle->SetLineWidth(2);

    // ========== 2. Open File and Get Objects ==========
    TFile *f = TFile::Open("hist_outputallFSI_more01jet.root");
    if (!f || f->IsZombie()){
        std::cout << "Error: Could not open file 'hist_outputallFSI_more01jet.root'!\n";
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
        
        hJetPtDistTR[i] = new TH1D(Form("hJetPtDistTR_%d", i), ";p_{T} (GeV/c);Normalized to Unity", 100, 0, 5);
        hJetPtDistTow[i] = new TH1D(Form("hJetPtDistTow_%d", i), ";p_{T} (GeV/c);Normalized to Unity", 100, 0, 5);
    }

    // ========== 6. Loop Over TTree and Fill Data ==========
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
                
                sumNjTR[j] += NjetTR;
                hJetNumDistTR[j]->Fill(NjetTR);
                if (tJetConstN) { for (int nConst : *tJetConstN) { sumTotConstTR[j] += nConst; hJetMultDistTR[j]->Fill(nConst); } }
                if (tJetPt) { for (double pt : *tJetPt) { hJetPtDistTR[j]->Fill(pt); } }
                
                sumNjTow[j] += NjetTow;
                hJetNumDistTow[j]->Fill(NjetTow);
                if (tJetTowConstN){ for (int nConst : *tJetTowConstN) { sumTotConstTow[j] += nConst; hJetMultDistTow[j]->Fill(nConst); } }
                if (tJetTowPt) { for (double pt : *tJetTowPt) { hJetPtDistTow[j]->Fill(pt); } }
                break;
            }
        }
    }
    std::cout << "\nEvent processing complete." << std::endl;

    // ========== 7. Calculate Final Averages ==========
    for (int i = 0; i < nRtBins; ++i){
        if (sumEvt[i] > 0){
            avgNjTR[i]   = sumNjTR[i] / sumEvt[i];
            avgMultTR[i] = (sumNjTR[i] > 0) ? sumTotConstTR[i] / sumNjTR[i] : 0;
            avgNjTow[i]   = sumNjTow[i] / sumEvt[i];
            avgMultTow[i] = (sumNjTow[i] > 0) ? sumTotConstTow[i] / sumNjTow[i] : 0;
        }
    }

    // ========== 8. Normalize and Create Main Comparison Plot ==========
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
    c1->SetTicks(1,1);
    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Normalized Jet Properties vs. R_{T};R_{T} = N_{ch}^{TR} / #LTN_{ch}^{TR}#GT;Normalized Yield");

    TGraph *gNjTRNorm   = new TGraph(nRtBins, rtCenter, avgNjTRNorm);
    gNjTRNorm->SetMarkerStyle(kFullCircle); gNjTRNorm->SetMarkerColor(kRed+1); gNjTRNorm->SetLineColor(kRed+1);
    
    TGraph *gMultTRNorm = new TGraph(nRtBins, rtCenter, avgMultTRNorm);
    gMultTRNorm->SetMarkerStyle(kFullSquare); gMultTRNorm->SetMarkerColor(kAzure+2); gMultTRNorm->SetLineColor(kAzure+2);

    TGraph *gNjTowNorm   = new TGraph(nRtBins, rtCenter, avgNjTowNorm);
    gNjTowNorm->SetMarkerStyle(kOpenCircle); gNjTowNorm->SetMarkerColor(kRed+1); gNjTowNorm->SetLineColor(kRed+1); gNjTowNorm->SetLineStyle(kDashed);

    TGraph *gMultTowNorm = new TGraph(nRtBins, rtCenter, avgMultTowNorm);
    gMultTowNorm->SetMarkerStyle(kOpenSquare); gMultTowNorm->SetMarkerColor(kAzure+2); gMultTowNorm->SetLineColor(kAzure+2); gMultTowNorm->SetLineStyle(kDashed);

    mg->Add(gNjTRNorm, "PL"); mg->Add(gMultTRNorm, "PL"); mg->Add(gNjTowNorm, "PL"); mg->Add(gMultTowNorm, "PL");
    mg->Draw("A");
    mg->GetXaxis()->SetTitleSize(0.045); mg->GetXaxis()->SetRangeUser(0, 5);
    mg->GetYaxis()->SetTitleSize(0.045); mg->GetYaxis()->SetRangeUser(0, 1.4);

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
    TCanvas *c2 = new TCanvas("c2", "Jet Number Distributions", 1000, 800);
    c2->Divide(2, 2, 0, 0); // --- Set spacing to zero ---
    c2->SetTicks(1, 1);
    for (int i = 0; i < 4; ++i) {
        c2->cd(i + 1);

        bool isLeft = (i % 2 == 0);
        bool isBottom = (i >= 2);
        gPad->SetRightMargin(!isLeft ? 0.1 : 0.0);
        gPad->SetTopMargin(!isBottom ? 0.05 : 0.0);
        gPad->SetLeftMargin(isLeft ? 0.2 : 0.0);
        gPad->SetBottomMargin(isBottom ? 0.2 : 0.0);
        gPad->SetTicks(1,1);

        int binIdx = selectedBins[i];
             
        TH1D* hTR = hJetNumDistTR[binIdx];
        TH1D* hTow = hJetNumDistTow[binIdx];
        if (hTR->GetEntries() > 0) hTR->Scale(1.0 / hTR->Integral(), "width");
        if (hTow->GetEntries() > 0) hTow->Scale(1.0 / hTow->Integral(), "width");

        double yMax = (i < 2) ? 0.77 : 0.37;
        double yMin = (i < 2) ? 0.0001 : 0;
        double XMin = (i % 2 == 0) ? 0 : 0.01;
        TH2D* frame = new TH2D(Form("frame_c2_%d", i), "", 10, XMin, 19.7, 10, yMin, yMax);
        frame->GetXaxis()->SetTitle("N_{jet}");
        if (isLeft) frame->GetYaxis()->SetTitle("Normalized to Unity");
        
        // Hide labels and titles for inner axes
        if (!isBottom) { frame->GetXaxis()->SetLabelSize(0); frame->GetXaxis()->SetTitleSize(0); }
        if (!isLeft)   { frame->GetYaxis()->SetLabelSize(0); frame->GetYaxis()->SetTitleSize(0); }
        
        // Adjust fonts and offsets for visibility
        frame->GetXaxis()->SetLabelSize(0.07); 
        frame->GetYaxis()->SetLabelSize(0.07);
        frame->GetXaxis()->SetTitleSize(0.08); 
        if (isLeft) {
            if (i == 0) { 
                frame->GetYaxis()->SetTitleSize(0.09); 
            } else {  
                frame->GetYaxis()->SetTitleSize(0.08); 
            }
        }
        frame->GetXaxis()->SetTitleOffset(1.1); 
        frame->GetYaxis()->SetTitleOffset(1.1); 
        frame->Draw("axis");

        TGraphErrors* gTR  = convertToGraph(hJetNumDistTR[binIdx], kAzure+2, kFullCircle);
        TGraphErrors* gTow = convertToGraph(hJetNumDistTow[binIdx], kRed+1, kFullCircle);
        gTR->Draw("LP SAME");
        gTow->Draw("LP SAME");

        TPaveText *pt = nullptr; 
        switch (i) {
            case 0: // 子画布 1 (左上)
                pt = new TPaveText(0.45, 0.78, 1.05, 0.93, "NDC");
                pt->SetTextSize(0.09);
                break;
            case 1: // 子画布 2 (右上)
                pt = new TPaveText(0.35, 0.78, 0.95, 0.93, "NDC");
                pt->SetTextSize(0.09);
                break;
            case 2: // 子画布 3 (左下)
                pt = new TPaveText(0.45, 0.82, 1.05, 0.97, "NDC");
                pt->SetTextSize(0.08);
                break;
            case 3: // 子画布 4 (右下)
                pt = new TPaveText(0.35, 0.82, 0.95, 0.97, "NDC"); // 移动到左下角
                pt->SetTextSize(0.08);
                break;
        }

        if (pt) {
            pt->SetFillStyle(0); 
            pt->SetBorderSize(0);
            pt->SetTextAlign(12); // 文本左对齐、居中
            pt->AddText(Form("%.1f #leq R_{T} < %.1f", rtEdges[binIdx], rtEdges[binIdx+1]));
            pt->Draw();
        }

        if (i == 0) {
            TLegend* leg_dist = new TLegend(0.5, 0.6, 0.9, 0.75);
            leg_dist->SetTextSize(0.08); 
            leg_dist->AddEntry(gTR, "Transverse", "lp");
            leg_dist->AddEntry(gTow, "Toward", "lp");
            leg_dist->Draw();
        }
    }
    c2->SaveAs("jet_number_distributions.png");

    // --- Jet Multiplicity Distributions ---
    TCanvas *c3 = new TCanvas("c3", "Jet Constituent Multiplicity Distributions", 1000, 800);
    c3->Divide(2, 2, 0, 0); // --- Set spacing to zero ---
    c3->SetTicks(1, 1);
    for (int i = 0; i < 4; ++i) {
        c3->cd(i + 1);

        bool isLeft = (i % 2 == 0);
        bool isBottom = (i >= 2);

        gPad->SetRightMargin(!isLeft ? 0.1 : 0.0);
        gPad->SetTopMargin(!isBottom ? 0.05 : 0.0);
        gPad->SetLeftMargin(isLeft ? 0.2 : 0.0);
        gPad->SetBottomMargin(isBottom ? 0.2 : 0.0);
        gPad->SetTicks(1,1);
        
        int binIdx = selectedBins[i];

        TH1D* hTR = hJetMultDistTR[binIdx];
        TH1D* hTow = hJetMultDistTow[binIdx];
        if (hTR->GetEntries() > 0) hTR->Scale(1.0 / hTR->Integral(), "width");
        if (hTow->GetEntries() > 0) hTow->Scale(1.0 / hTow->Integral(), "width");

        double yMax = (i < 2) ? 0.77 : 0.37;
        double yMin = (i < 2) ? 0.0001 : 0;
        double XMin = (i % 2 == 0) ? 0 : 0.01;
        TH2D* frame = new TH2D(Form("frame_c3_%d", i), "", 10, XMin, 19.7, 10, yMin, yMax);
        frame->GetXaxis()->SetTitle("n_{const}");
        if (isLeft) frame->GetYaxis()->SetTitle("Normalized to Unity");

        if (!isBottom) { frame->GetXaxis()->SetLabelSize(0); frame->GetXaxis()->SetTitleSize(0); }
        if (!isLeft)   { frame->GetYaxis()->SetLabelSize(0); frame->GetYaxis()->SetTitleSize(0); }
        
        frame->GetXaxis()->SetLabelSize(0.07);
        frame->GetYaxis()->SetLabelSize(0.07);
        frame->GetXaxis()->SetTitleSize(0.08);
        if (isLeft) {
            if (i == 0) { // 第一个子画板 (左上)
                frame->GetYaxis()->SetTitleSize(0.09); // 增大
            } else {      // 第三个子画板 (左下, i=2)
                frame->GetYaxis()->SetTitleSize(0.08); // 减小
            }
        }
        frame->GetXaxis()->SetTitleOffset(1.1);
        frame->GetYaxis()->SetTitleOffset(1.1);
        frame->Draw("axis");

        TGraphErrors* gTR  = convertToGraph(hJetMultDistTR[binIdx], kAzure+2, kFullCircle);
        TGraphErrors* gTow = convertToGraph(hJetMultDistTow[binIdx], kRed+1, kFullCircle);
        gTR->Draw("LP SAME");
        gTow->Draw("LP SAME");
        TPaveText *pt = nullptr; 
        switch (i) {
            case 0: // 子画布 1 (左上)
                pt = new TPaveText(0.45, 0.78, 1.05, 0.93, "NDC");
                pt->SetTextSize(0.09);
                break;
            case 1: // 子画布 2 (右上)
                pt = new TPaveText(0.35, 0.78, 0.95, 0.93, "NDC");
                pt->SetTextSize(0.09);
                break;
            case 2: // 子画布 3 (左下)
                pt = new TPaveText(0.45, 0.82, 1.05, 0.97, "NDC");
                pt->SetTextSize(0.08);
                break;
            case 3: // 子画布 4 (右下)
                pt = new TPaveText(0.35, 0.82, 0.95, 0.97, "NDC"); // 移动到左下角
                pt->SetTextSize(0.08);
                break;
        }

        if (pt) {
            pt->SetFillStyle(0); 
            pt->SetBorderSize(0);
            pt->SetTextAlign(12); // 文本左对齐、居中
            pt->AddText(Form("%.1f #leq R_{T} < %.1f", rtEdges[binIdx], rtEdges[binIdx+1]));
            pt->Draw();
        }

        if (i == 0) {
            TLegend* leg_dist = new TLegend(0.5, 0.6, 0.9, 0.75);
            leg_dist->SetTextSize(0.08); // 您可以调整这个值
            leg_dist->AddEntry(gTR, "Transverse", "lp");
            leg_dist->AddEntry(gTow, "Toward", "lp");
            leg_dist->Draw();
        }
    }
    c3->SaveAs("jet_multiplicity_distributions.png");

    // ========== 10. Plot Jet pT Distributions with Log Y-axis ==========
    TCanvas *c4 = new TCanvas("c4", "Jet pT Distributions", 1000, 800);
    c4->Divide(2, 2, 0, 0); // --- Set spacing to zero ---
    c4->SetTicks(1, 1);
    for (int i = 0; i < 4; ++i) {
        c4->cd(i + 1);

        bool isLeft = (i % 2 == 0);
        bool isBottom = (i >= 2);

        gPad->SetRightMargin(!isLeft ? 0.1 : 0.0);
        gPad->SetTopMargin(!isBottom ? 0.05 : 0.0);
        gPad->SetLeftMargin(isLeft ? 0.2 : 0.0);
        gPad->SetBottomMargin(isBottom ? 0.2 : 0.0);
        gPad->SetTicks(1,1);
        gPad->SetLogy(1);

        int binIdx = selectedBins[i];
        
        TH1D* hTR = hJetPtDistTR[binIdx];
        TH1D* hTow = hJetPtDistTow[binIdx];
        
        hTR->Rebin(4);
        hTow->Rebin(4);

        double N_events = sumEvt[binIdx];
        if (hTR->GetEntries() > 0 && N_events > 0) hTR->Scale(1.0 / N_events, "width");
        if (hTow->GetEntries() > 0 && N_events > 0) hTow->Scale(1.0 / N_events, "width");
        
        double yMin = (i < 2) ? 0.002 : 1e-1;
        double yMax = 9.7;
        double XMin = (i % 2 == 0) ? 0 : 0.01;
        TH2D* frame = new TH2D(Form("frame_c4_%d", i), "", 10, XMin, 5.4, 10, yMin, yMax);
        frame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        if (isLeft) frame->GetYaxis()->SetTitle("dN/(dp_{T}dy) (GeV/c)^{-1}");

        if (!isBottom) { frame->GetXaxis()->SetLabelSize(0); frame->GetXaxis()->SetTitleSize(0); }
        if (!isLeft)   { frame->GetYaxis()->SetLabelSize(0); frame->GetYaxis()->SetTitleSize(0); }
        
        frame->GetXaxis()->SetLabelSize(0.07);
        frame->GetYaxis()->SetLabelSize(0.07);
        frame->GetXaxis()->SetTitleSize(0.08);
        if (isLeft) {
            if (i == 0) { 
                frame->GetYaxis()->SetTitleSize(0.09); 
            } else { 
                frame->GetYaxis()->SetTitleSize(0.08);
            }
        }
        frame->GetXaxis()->SetTitleOffset(1.1);
        frame->GetYaxis()->SetTitleOffset(1.1);
        frame->Draw("axis");

        TGraphErrors* gTR  = convertToGraph(hTR, kAzure+2, kFullCircle);
        TGraphErrors* gTow = convertToGraph(hTow, kRed+1, kFullCircle);
        gTR->Draw("LP SAME");
        gTow->Draw("LP SAME");

        TPaveText *pt = nullptr; 
        switch (i) {
            case 0: // 子画布 1 (左上)
                pt = new TPaveText(0.45, 0.78, 1.05, 0.93, "NDC");
                pt->SetTextSize(0.09);
                break;
            case 1: // 子画布 2 (右上)
                pt = new TPaveText(0.35, 0.78, 0.95, 0.93, "NDC");
                pt->SetTextSize(0.09);
                break;
            case 2: // 子画布 3 (左下)
                pt = new TPaveText(0.45, 0.82, 1.05, 0.97, "NDC");
                pt->SetTextSize(0.08);
                break;
            case 3: // 子画布 4 (右下)
                pt = new TPaveText(0.35, 0.82, 0.95, 0.97, "NDC"); // 移动到左下角
                pt->SetTextSize(0.08);
                break;
        }

        if (pt) {
            pt->SetFillStyle(0); 
            pt->SetBorderSize(0);
            pt->SetTextAlign(12); // 文本左对齐、居中
            pt->AddText(Form("%.1f #leq R_{T} < %.1f", rtEdges[binIdx], rtEdges[binIdx+1]));
            pt->Draw();
        }

        if (i == 0) {
            TLegend* leg_dist_pt = new TLegend(0.25, 0.10, 0.55, 0.30);
            leg_dist_pt->SetTextSize(0.08);
            leg_dist_pt->AddEntry(gTR, "Transverse", "lp"); 
            leg_dist_pt->AddEntry(gTow, "Toward", "lp");  
            leg_dist_pt->Draw();
        }
    }
    c4->SaveAs("jet_pt_distributions.png");
    f->Close();
}