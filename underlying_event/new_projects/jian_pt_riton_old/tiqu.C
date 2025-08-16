#include <iostream>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>

// Subtract two TGraphAsymmErrors
TGraphAsymmErrors* SubtractGraphs(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2) {
    int nPoints = g1->GetN();
    TGraphAsymmErrors* result = new TGraphAsymmErrors(nPoints);

    for (int i = 0; i < nPoints; ++i) {
        double x = g1->GetX()[i];
        double y1 = g1->GetY()[i];
        double y2 = g2->GetY()[i];
        double dy1L = g1->GetEYlow()[i];
        double dy1H = g1->GetEYhigh()[i];
        double dy2L = g2->GetEYlow()[i];
        double dy2H = g2->GetEYhigh()[i];

        double y = y1 - y2;
        double dyL = TMath::Sqrt(dy1L * dy1L + dy2L * dy2L);
        double dyH = TMath::Sqrt(dy1H * dy1H + dy2H * dy2H);

        result->SetPoint(i, x, y);
        result->SetPointError(i, 0, 0, dyL, dyH); // Set symmetric errors for simplicity
    }

    return result;
}

// Divide two TGraphAsymmErrors at the same x points
TGraphAsymmErrors* DivideGraphsAtSameX(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2, double tolerance = 1e-6) {
    int nPoints1 = g1->GetN();
    int nPoints2 = g2->GetN();
    TGraphAsymmErrors* result = new TGraphAsymmErrors();

    for (int i = 0; i < nPoints1; ++i) {
        double x1 = g1->GetX()[i];
        double y1 = g1->GetY()[i];
        double dy1L = g1->GetEYlow()[i];
        double dy1H = g1->GetEYhigh()[i];

        for (int j = 0; j < nPoints2; ++j) {
            double x2 = g2->GetX()[j];
            double y2 = g2->GetY()[j];
            double dy2L = g2->GetEYlow()[j];
            double dy2H = g2->GetEYhigh()[j];

            if (TMath::Abs(x1 - x2) < tolerance) { // Check if x values are within the tolerance
                if (y2 != 0) {
                    double y = y1 / y2;
                    double dyL = y * TMath::Sqrt((dy1L * dy1L) / (y1 * y1) + (dy2L * dy2L) / (y2 * y2));
                    double dyH = y * TMath::Sqrt((dy1H * dy1H) / (y1 * y1) + (dy2H * dy2H) / (y2 * y2));

                    result->SetPoint(result->GetN(), x1, y);
                    result->SetPointError(result->GetN() - 1, 0, 0, dyL, dyH);
                } else {
                    std::cerr << "Warning: Denominator is zero at point (" << x1 << ", " << y1 << "). Skipping this point." << std::endl;
                }
                break; // Exit the inner loop once a match is found
            }
        }
    }

    return result;
}

void tiqu() {
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600); // 创建Canvas

    TFile *input1 = new TFile("HEPData-ins2626034-v1-Table_33.root"); // 修改一：改为下载的名字
    TDirectoryFile *file1 = (TDirectoryFile*)input1->Get("Table 33"); // 修改三：提取的数据是 table 什么就改为什么
    TGraphAsymmErrors *hist1 = (TGraphAsymmErrors*)file1->Get("Graph1D_y1"); // 修改四

    TFile *input2 = new TFile("HEPData-ins2626034-v1-Table_43.root"); // 修改一：改为下载的名字
    TDirectoryFile *file2 = (TDirectoryFile*)input2->Get("Table 43"); // 修改三：提取的数据是 table 什么就改为什么
    TGraphAsymmErrors *hist2 = (TGraphAsymmErrors*)file2->Get("Graph1D_y1"); // 修改四

    TFile *input3 = new TFile("HEPData-ins2626034-v1-Table_3.root"); // 修改一：改为下载的名字
    TDirectoryFile *file3 = (TDirectoryFile*)input3->Get("Table 3"); // 修改三：提取的数据是 table 什么就改为什么
    TGraphAsymmErrors *hist3 = (TGraphAsymmErrors*)file3->Get("Graph1D_y1"); // 修改四

    TFile *input4 = new TFile("HEPData-ins2626034-v1-Table_13.root"); // 修改一：改为下载的名字
    TDirectoryFile *file4 = (TDirectoryFile*)input4->Get("Table 13"); // 修改三：提取的数据是 table 什么就改为什么
    TGraphAsymmErrors *hist4 = (TGraphAsymmErrors*)file4->Get("Graph1D_y1"); // 修改四

    // Perform the subtraction operation: hist1 - hist2
    TGraphAsymmErrors* diffGraph1 = SubtractGraphs(hist1, hist2);

    // Perform the subtraction operation: hist3 - hist4
    TGraphAsymmErrors* diffGraph2 = SubtractGraphs(hist3, hist4);

    // Perform the division operation only at points with the same x value
    TGraphAsymmErrors* divGraph = DivideGraphsAtSameX(diffGraph1, diffGraph2, 1e-6);

    // Draw the result
    divGraph->Draw();

    // Save the result to a ROOT file
    TFile *outputFile = new TFile("p.pi.root", "RECREATE");
    divGraph->Write("Pt");
    outputFile->Close();

    // Close input files
    input1->Close();
    input2->Close();
    input3->Close();
    input4->Close();
}