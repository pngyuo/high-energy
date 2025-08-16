#include <iostream>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TFile.h>

void jian() {
    TString filenames1 = "Pi_To.root";
    TString filenames2 = "Pi_Tr.root";

    TFile* infile1 = new TFile(filenames1, "READ");
    TGraphErrors *graph1 = (TGraphErrors*)infile1->Get("Pt");

    TFile* infile2 = new TFile(filenames2, "READ");
    TGraphErrors *graph2 = (TGraphErrors*)infile2->Get("Pt");

    if (graph1->GetN() != graph2->GetN()) {
        std::cerr << "Error: Graphs do not have the same number of points." << std::endl;
        return;
    }

    int n = graph1->GetN();
    double* x = new double[n];
    double* y = new double[n];
    double* ey = new double[n];

    for (int i = 0; i < n; ++i) {
        x[i] = graph1->GetX()[i];
        y[i] = graph1->GetY()[i] - graph2->GetY()[i];
        ey[i] = std::sqrt(std::pow(graph1->GetEY()[i], 2) + std::pow(graph2->GetEY()[i], 2));
    }

    TGraphErrors* diffGraph = new TGraphErrors(n, x, y, 0, ey);
    diffGraph->SetName("Pt");
    diffGraph->SetTitle("Difference between Graph 1 and Graph 2");

    TCanvas* c = new TCanvas();
    diffGraph->Draw("AP");

    // 指定输出文件名
    TFile* outfile = new TFile("Pi.root", "RECREATE");
    diffGraph->Write();
    outfile->Close();

    delete [] x;
    delete [] y;
    delete [] ey;

    infile1->Close();
    delete infile1;
    infile2->Close();
    delete infile2;
}