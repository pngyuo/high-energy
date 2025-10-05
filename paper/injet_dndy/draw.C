#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>

TGraph* drawConnectedPoints(TH1D* hist, Color_t color, Style_t markerStyle, double& averageNT) {
    int n = hist->GetNbinsX();
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> ey;

    for (int i = 1; i <= n; ++i) {
        double binContent = hist->GetBinContent(i);
        double binCenter = hist->GetBinCenter(i);
        if (binContent != 0 && binCenter < 30.0 / averageNT) {
            double binError = hist->GetBinError(i);
            if (binCenter < 70.0 / averageNT) {
                x.push_back(binCenter);
                y.push_back(binContent);
                ey.push_back(binError);
            }
        }
    }

    int nPoints = x.size();
    if (nPoints == 0) return nullptr;

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
    return graph;
}

void addText(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);
    latex->SetTextSize(0.04);
    latex->Draw();
}

void addText1(Double_t x, Double_t y, const char* text) {
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
void addText5(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);
    latex->SetTextSize(0.10);
    latex->Draw();
}
void addText4(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);
    latex->SetTextSize(0.15);
    latex->Draw();
}

void addText6(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);
    latex->SetTextSize(0.07);
    latex->Draw();
}

void addText7(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);
    latex->SetTextSize(0.14);
    latex->Draw();
}

void addText8(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);
    latex->SetTextSize(0.09);
    latex->Draw();
}

double averageNT = 0.0;

TH1D* getCentPtProjections(const TString& filename, const TString& histname, double& averageNT, int rebinFactor1, int rebinFactor2) {
    TFile* infile = new TFile(filename, "READ");
    TH3D* hist3D = (TH3D*)infile->Get(histname);
    TH1D* chargedParticlesCount = (TH1D*)infile->Get("nTHist");

    double totalChargedParticles = 0;
    double nEvents = 0;
    for (int bin = 1; bin <= 36; bin++) {
        double binCenter = chargedParticlesCount->GetXaxis()->GetBinCenter(bin);
        double binContent = chargedParticlesCount->GetBinContent(bin);
        totalChargedParticles += binCenter * binContent;
        nEvents += binContent;
    }
    averageNT = totalChargedParticles / nEvents;
    std::cout << "averageNT: " << averageNT << std::endl;

    TH1D* histProjZ = hist3D->ProjectionZ("histProjZ", 1, hist3D->GetNbinsX(), 1, hist3D->GetNbinsY());
    int nbins = histProjZ->GetNbinsX();
    double xmin = histProjZ->GetXaxis()->GetBinLowEdge(1);
    double xmax = histProjZ->GetXaxis()->GetBinUpEdge(nbins);

    std::vector<double> newBins;
    double currentX = xmin / averageNT;
    double binWidth = histProjZ->GetXaxis()->GetBinWidth(1) / averageNT;

    while (currentX < 2.5) {
        newBins.push_back(currentX);
        currentX += rebinFactor1 * binWidth;
    }
    while (currentX < 5.0) {
        newBins.push_back(currentX);
        currentX += rebinFactor2 * binWidth;
    }
    newBins.push_back(xmax / averageNT);

    TH1D* scaledHistProjZ = new TH1D("scaledHistProjZ", "Scaled Projection Z", newBins.size() - 1, &newBins[0]);
    for (int i = 1; i <= scaledHistProjZ->GetNbinsX(); ++i) {
        scaledHistProjZ->SetBinContent(i, 0.0);
        scaledHistProjZ->SetBinError(i, 0.0);
    }

    for (int i = 1; i <= nbins; ++i) {
        double binContent = histProjZ->GetBinContent(i);
        double binError = histProjZ->GetBinError(i);
        double binCenter = histProjZ->GetBinCenter(i);
        double rtValue = binCenter / averageNT;
        int ntBin = chargedParticlesCount->FindBin(binCenter);
        double eventsInBin = chargedParticlesCount->GetBinContent(ntBin);

        if (eventsInBin > 0) {
            double normContent = binContent / eventsInBin;
            double normError = binError / eventsInBin;
            int newBin = scaledHistProjZ->FindBin(rtValue);
            scaledHistProjZ->SetBinContent(newBin, normContent);
            scaledHistProjZ->SetBinError(newBin, normError);
        }
    }

    return scaledHistProjZ;
}

void drawGraphs(TH1D* projectionsToward, TH1D* projectionsTransverse, TVirtualPad* pad, double& averageNT) {
    pad->cd();
    std::vector<Color_t> colors = {kBlack, kRed, kBlue, kMagenta};
    std::vector<Style_t> markerStyles = {20, 20, 22};

    TGraph* graphToward = drawConnectedPoints(projectionsToward, colors[0], markerStyles[0], averageNT);
    if (graphToward) {
        graphToward->SetLineWidth(2);
        graphToward->Draw("LP");
    }

    TGraph* graphTransverse = drawConnectedPoints(projectionsTransverse, colors[1], markerStyles[1], averageNT);
    if (graphTransverse) {
        graphTransverse->SetLineWidth(2);
        graphTransverse->Draw("LP same");
    }
}

void draw() {
    TCanvas* c1 = new TCanvas("c1", "c1", 1500, 1200);
    c1->SetTicks(1, 1);
    c1->SetMargin(0.17, 0.10, 0.2, 0.2);
    c1->Divide(3, 3, 0, 0);
    gStyle->SetOptStat(kFALSE);

    // 第一行：Pi±
    {
        // 第一个子画布：noFSI
        c1->cd(1);
        gPad->SetLogy();
        gPad->SetTicks(1, 1);
        TH2D* axisFrame1 = new TH2D("axisFrame1", "; ; ", 100, -0.1, 4.99, 100, 0.6, 40);
        axisFrame1->GetYaxis()->SetTickLength(0.02);
        axisFrame1->GetXaxis()->SetTickLength(0.02);
        axisFrame1->GetXaxis()->SetLabelSize(0.085);
        axisFrame1->GetYaxis()->SetLabelSize(0.085);
        axisFrame1->Draw("axis");

        TString filename = "hist_outputnohFSI_01mb.root";
        TString histnameToward = "hPiCh_dPhi0";
        TString histnameTransverse = "hPiCh_dPhi1";
        int rebinFactor1 = 2;
        int rebinFactor2 = 2;
        double averageNT;

        TLegend* legend = new TLegend(0.60, 0.27, 1.13, 0.52);
        legend->SetNColumns(1);
        legend->SetFillStyle(0);
        legend->SetBorderSize(0);
        legend->SetTextSize(0.082);

        TGraphErrors* marker1 = new TGraphErrors();
        marker1->SetMarkerStyle(20);
        marker1->SetMarkerColor(kBlack);
        marker1->SetMarkerSize(1);
        marker1->SetLineWidth(2);
        marker1->SetLineColor(kBlack);
        legend->AddEntry(marker1, "Toward", "lp");

        TGraphErrors* marker2 = new TGraphErrors();
        marker2->SetMarkerStyle(20);
        marker2->SetMarkerColor(kRed);
        marker2->SetMarkerSize(1);
        marker2->SetLineWidth(2);
        marker2->SetLineColor(kRed);
        legend->AddEntry(marker2, "Transverse", "lp");
        legend->Draw();

        TH1D* projectionsToward = getCentPtProjections(filename, histnameToward, averageNT, rebinFactor1, rebinFactor2);
        TH1D* projectionsTransverse = getCentPtProjections(filename, histnameTransverse, averageNT, rebinFactor1, rebinFactor2);
        drawGraphs(projectionsToward, projectionsTransverse, gPad, averageNT);
        addText5(0.52, 0.90, "0.10mb w/o ART");
        addText6(0.92, 0.03, "(a)");
        addText4(0.23, 0.85, "#pi^{+}+#pi^{-}");
    }

    {
        // 第二个子画布：nohFSI
        c1->cd(2);
        gPad->SetLogy();
        gPad->SetTicks(1, 1);
        TH2D* axisFrame2 = new TH2D("axisFrame2", "; ; ", 100, -0.1, 4.99, 100, 0.6, 40);
        axisFrame2->GetYaxis()->SetTickLength(0.02);
        axisFrame2->GetXaxis()->SetTickLength(0.02);
        axisFrame2->GetXaxis()->SetLabelSize(0.085);
        axisFrame2->GetYaxis()->SetLabelSize(0.085);
        axisFrame2->Draw("axis");

        TString filename = "hist_outputnohFSI_liang.root";
        TString histnameToward = "hPiCh_dPhi0";
        TString histnameTransverse = "hPiCh_dPhi1";
        int rebinFactor1 = 2;
        int rebinFactor2 = 2;
        double averageNT;

        TH1D* projectionsToward = getCentPtProjections(filename, histnameToward, averageNT, rebinFactor1, rebinFactor2);
        TH1D* projectionsTransverse = getCentPtProjections(filename, histnameTransverse, averageNT, rebinFactor1, rebinFactor2);
        drawGraphs(projectionsToward, projectionsTransverse, gPad, averageNT);
        addText5(0.42, 0.90, "0.15mb w/o ART");

        addText6(0.90, 0.03, "(d)");
    }

    {
        // 第三个子画布：allFSI
        c1->cd(3);
        gPad->SetLogy();
        gPad->SetRightMargin(0.01);
        gPad->SetTicks(1, 1);
        TH2D* axisFrame3 = new TH2D("axisFrame3", "; ; ", 100, -0.1, 4.99, 100, 0.6, 40);
        axisFrame3->GetYaxis()->SetTickLength(0.02);
        axisFrame3->GetXaxis()->SetTickLength(0.02);
        axisFrame3->GetXaxis()->SetLabelSize(0.085);
        axisFrame3->GetYaxis()->SetLabelSize(0.085);
        axisFrame3->Draw("axis");

        TString filename = "hist_outputnohFSI_02mb.root";
        TString histnameToward = "hPiCh_dPhi0";
        TString histnameTransverse = "hPiCh_dPhi1";
        int rebinFactor1 = 2;
        int rebinFactor2 = 2;
        double averageNT;

        TH1D* projectionsToward = getCentPtProjections(filename, histnameToward, averageNT, rebinFactor1, rebinFactor2);
        TH1D* projectionsTransverse = getCentPtProjections(filename, histnameTransverse, averageNT, rebinFactor1, rebinFactor2);
        drawGraphs(projectionsToward, projectionsTransverse, gPad, averageNT);
        addText5(0.42, 0.90, "0.20mb w/o ART");
        addText6(0.90, 0.03, "(j)");
    }

    // 第二行：p±
    {
        // 第四个子画布：noFSI
        c1->cd(4);
        gPad->SetLogy();
        gPad->SetTicks(1, 1);
        TH2D* axisFrame4 = new TH2D("axisFrame4", "; ; ", 100, -0.1, 4.99, 100, 0.06, 12);
        axisFrame4->GetYaxis()->SetTickLength(0.02);
        axisFrame4->GetXaxis()->SetTickLength(0.02);
        axisFrame4->GetXaxis()->SetLabelSize(0.085);
        axisFrame4->GetYaxis()->SetLabelSize(0.085);
        axisFrame4->Draw("axis");

        TString filename = "hist_outputnohFSI_01mb.root";
        TString histnameToward = "hKCh_dPhi0";
        TString histnameTransverse = "hKCh_dPhi1";
        int rebinFactor1 = 2;
        int rebinFactor2 = 2;
        double averageNT;

        TLatex* label1 = new TLatex();
        label1->SetTextSize(0.12);
        label1->SetTextAngle(90);
        label1->SetTextAlign(22);
        label1->DrawLatexNDC(0.05, 0.54, "dN/dy");

        TH1D* projectionsToward = getCentPtProjections(filename, histnameToward, averageNT, rebinFactor1, rebinFactor2);
        TH1D* projectionsTransverse = getCentPtProjections(filename, histnameTransverse, averageNT, rebinFactor1, rebinFactor2);
        drawGraphs(projectionsToward, projectionsTransverse, gPad, averageNT);
        addText6(0.92, 0.03, "(b)");
        addText4(0.23, 0.85, "#Kappa^{+}+#Kappa^{-}");

    }

    {
        // 第五个子画布：nohFSI
        c1->cd(5);
        gPad->SetLogy();
        gPad->SetTicks(1, 1);
        TH2D* axisFrame5 = new TH2D("axisFrame5", "; ; ", 100, -0.1, 4.99, 100, 0.06, 12);
        axisFrame5->GetYaxis()->SetTickLength(0.02);
        axisFrame5->GetXaxis()->SetTickLength(0.02);
        axisFrame5->GetXaxis()->SetLabelSize(0.085);
        axisFrame5->GetYaxis()->SetLabelSize(0.085);
        axisFrame5->Draw("axis");

        TString filename = "hist_outputnohFSI_liang.root";
        TString histnameToward = "hKCh_dPhi0";
        TString histnameTransverse = "hKCh_dPhi1";
        int rebinFactor1 = 2;
        int rebinFactor2 = 2;
        double averageNT;

        TH1D* projectionsToward = getCentPtProjections(filename, histnameToward, averageNT, rebinFactor1, rebinFactor2);
        TH1D* projectionsTransverse = getCentPtProjections(filename, histnameTransverse, averageNT, rebinFactor1, rebinFactor2);
        drawGraphs(projectionsToward, projectionsTransverse, gPad, averageNT);
        addText6(0.90, 0.03, "(e)");
    }

    {
        // 第六个子画布：allFSI
        c1->cd(6);
        gPad->SetLogy();
        gPad->SetRightMargin(0.01);
        gPad->SetTicks(1, 1);
        TH2D* axisFrame6 = new TH2D("axisFrame6", "; ; ", 100, -0.1, 4.99, 100, 0.06, 12);
        axisFrame6->GetYaxis()->SetTickLength(0.02);
        axisFrame6->GetXaxis()->SetTickLength(0.02);
        axisFrame6->GetXaxis()->SetLabelSize(0.085);
        axisFrame6->GetYaxis()->SetLabelSize(0.085);
        axisFrame6->Draw("axis");

        TString filename = "hist_outputnohFSI_02mb.root";
        TString histnameToward = "hKCh_dPhi0";
        TString histnameTransverse = "hKCh_dPhi1";
        int rebinFactor1 = 2;
        int rebinFactor2 = 2;
        double averageNT;

        TH1D* projectionsToward = getCentPtProjections(filename, histnameToward, averageNT, rebinFactor1, rebinFactor2);
        TH1D* projectionsTransverse = getCentPtProjections(filename, histnameTransverse, averageNT, rebinFactor1, rebinFactor2);
        drawGraphs(projectionsToward, projectionsTransverse, gPad, averageNT);
        addText6(0.90, 0.03, "(h)");
    }

    // 第三行：π±
    {
        // 第七个子画布：noFSI
        c1->cd(7);
        gPad->SetLogy();
        gPad->SetTicks(1, 1);
        TH2D* axisFrame7 = new TH2D("axisFrame7", "; ; ", 100, -0.1, 4.99, 100, 0.006, 8);
        axisFrame7->GetYaxis()->SetTickLength(0.02);
        axisFrame7->GetXaxis()->SetTickLength(0.02);
        axisFrame7->GetXaxis()->SetLabelSize(0.075);
        axisFrame7->GetYaxis()->SetLabelSize(0.075);
        axisFrame7->Draw("axis");

        TString filename = "hist_outputnohFSI_01mb.root";
        TString histnameToward = "hProton_dPhi0";
        TString histnameTransverse = "hProton_dPhi1";
        int rebinFactor1 = 2;
        int rebinFactor2 = 2;
        double averageNT;

        TH1D* projectionsToward = getCentPtProjections(filename, histnameToward, averageNT, rebinFactor1, rebinFactor2);
        TH1D* projectionsTransverse = getCentPtProjections(filename, histnameTransverse, averageNT, rebinFactor1, rebinFactor2);
        drawGraphs(projectionsToward, projectionsTransverse, gPad, averageNT);
        addText2(0.92, 0.23, "(c)");
        addText7(0.23, 0.87, "p+#bar{p}");
    }

    {
        // 第八个子画布：nohFSI
        c1->cd(8);
        gPad->SetLogy();
        gPad->SetTicks(1, 1);
        TH2D* axisFrame8 = new TH2D("axisFrame8", "; ; ", 100, -0.1, 4.99, 100, 0.006, 8);
        axisFrame8->GetYaxis()->SetTickLength(0.02);
        axisFrame8->GetXaxis()->SetTickLength(0.02);
        axisFrame8->GetXaxis()->SetLabelSize(0.075);
        axisFrame8->GetYaxis()->SetLabelSize(0.075);
        axisFrame8->Draw("axis");

        TString filename = "hist_outputnohFSI_liang.root";
        TString histnameToward = "hProton_dPhi0";
        TString histnameTransverse = "hProton_dPhi1";
        int rebinFactor1 = 2;
        int rebinFactor2 = 2;
        double averageNT;

        TH1D* projectionsToward = getCentPtProjections(filename, histnameToward, averageNT, rebinFactor1, rebinFactor2);
        TH1D* projectionsTransverse = getCentPtProjections(filename, histnameTransverse, averageNT, rebinFactor1, rebinFactor2);
        drawGraphs(projectionsToward, projectionsTransverse, gPad, averageNT);

TLatex *text2 = new TLatex(0.53, 0.06, "R_{T}");
text2->SetNDC();
text2->SetTextSize(0.12);
text2->SetTextAlign(22);
text2->Draw();
        addText2(0.90, 0.23, "(f)");
    }

    {
        // 第九个子画布：allFSI
        c1->cd(9);
        gPad->SetLogy();
        gPad->SetRightMargin(0.01);
        gPad->SetTicks(1, 1);
        TH2D* axisFrame9 = new TH2D("axisFrame9", "; ; ", 100, -0.1, 4.99, 100, 0.006, 8);
        axisFrame9->GetYaxis()->SetTickLength(0.02);
        axisFrame9->GetXaxis()->SetTickLength(0.02);
        axisFrame9->GetXaxis()->SetLabelSize(0.075);
        axisFrame9->GetYaxis()->SetLabelSize(0.075);
        axisFrame9->Draw("axis");

        TString filename = "hist_outputnohFSI_02mb.root";
        TString histnameToward = "hProton_dPhi0";
        TString histnameTransverse = "hProton_dPhi1";
        int rebinFactor1 = 2;
        int rebinFactor2 = 2;
        double averageNT;

        TH1D* projectionsToward = getCentPtProjections(filename, histnameToward, averageNT, rebinFactor1, rebinFactor2);
        TH1D* projectionsTransverse = getCentPtProjections(filename, histnameTransverse, averageNT, rebinFactor1, rebinFactor2);
        drawGraphs(projectionsToward, projectionsTransverse, gPad, averageNT);
        addText2(0.90, 0.23, "(i)");
    }

    c1->SaveAs("pi_k_p_Yield_RT_mb.png");
    c1->SaveAs("pi_k_p_Yield_RT_mb.pdf");
}