#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TProfile.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>

// Custom step+expo function for negative trace
Double_t step_expo_function_neg(Double_t *x, Double_t *par)
{
    Double_t xx = x[0];
    // range for constant part
    Double_t constant_part_range_left = 0;
    Double_t constant_part_range_right = 100;
    // range for exponential decay part
    Double_t expo_part_range_left = 120;
    Double_t expo_part_range_right = 400;

    if (xx >= constant_part_range_left && xx <= constant_part_range_right) // range for constant part
        return par[0];
    else if (xx > expo_part_range_left && xx <= expo_part_range_right) // range for exponential decay part
        return par[0] + (par[1] * (1.0 - TMath::Exp(-(xx - expo_part_range_left) / par[2]))) - par[1];
    else
        // return par[0] - par[1]/2;
        // return TMath::QuietNaN();
        return 0;
}

// Custom step+expo function for positive trace
Double_t step_expo_function_pos(Double_t *x, Double_t *par)
{
    Double_t xx = x[0];
    // range for constant part
    Double_t constant_part_range_left = 0;
    Double_t constant_part_range_right = 105;
    // range for exponential decay part
    Double_t expo_part_range_left = 110;
    Double_t expo_part_range_right = 400;

    if (xx >= constant_part_range_left && xx <= constant_part_range_right) // range for constant part
        return par[0];
    else if (xx > expo_part_range_left && xx <= expo_part_range_right) // range for exponential decay part
        return par[0] + (par[1] * (TMath::Exp(-(xx - expo_part_range_left) / par[2])));
    else
        // return par[0] - par[1]/2;
        // return TMath::QuietNaN();
        return 0;
}

void fit_trace()
{
    //     Open ROOT file
    TFile *file = TFile::Open("/u/ddas/c4/shiyan_test/test_0003_tree.root");
    if (!file || file->IsZombie())
    {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    //     Set up TTreeReader
    TTreeReader reader("evt", file);
    TTreeReaderArray<LisaCalItem> LisaItem(reader, "LisaCalData");
    TTreeReaderValue<EventHeader> Header(reader, "EventHeader.");

    //     Graph to hold (x, y) points
    TGraph *graph = new TGraph();
    int pointIndex = 0;

    //     Process only the first event
    if (reader.Next())
    {
        std::cout << "Processing Event no: " << Header->GetEventno() << std::endl;

        for (const auto &item : LisaItem)
        {
            if (item.Get_layer_id() == 9 && item.Get_xposition() == 2 && item.Get_yposition() == 1)
            {
                std::vector<int16_t> trace_x = item.Get_trace_x();
                std::vector<float> trace_febex = item.Get_trace_febex();

                // std::vector<short> trace_x = item.Get_trace_x();
                // std::vector<short> trace_febex = item.Get_trace_febex();

                for (UShort_t i = 0; i < trace_febex.size(); ++i)
                // for (size_t i = 0; i < trace_febex.size(); ++i)
                {
                    graph->SetPoint(pointIndex++, trace_x[i], trace_febex[i]);
                }
            }
        }
    }

    //     Draw graph and fit
    TCanvas *c1 = new TCanvas("c1", "Trace Fit", 800, 600);
    graph->SetTitle("Trace Data;Time (x);Amplitude (y)");
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.6);
    graph->Draw("A*");

    TF1 *fitFunc = new TF1("fitFunc", step_expo_function_pos, 0, 400, 3);
    fitFunc->SetParNames("Const", "Expo_Amp", "Tau");
    fitFunc->SetParameters(8150, 300, 1000); // Initial guesses
    fitFunc->SetParLimits(0, 8100, 8200);    // Const
    fitFunc->SetParLimits(1, 200, 4000);     // Expo_Amp
    fitFunc->SetParLimits(2, 10, 2000);      // Tau
    fitFunc->SetLineColor(kRed);
    fitFunc->SetLineWidth(2);
    graph->SetMinimum(8000); // Lower bound of Y-axis
    graph->SetMaximum(9000); // Upper bound of Y-axis
    graph->GetYaxis()->SetTitle("Amplitude");
    graph->Fit(fitFunc, "Rrob=0.8", "", 0, 400);
    c1->Update();
}

void fit_trace_profile()
{
    //     Open file
    TFile *file = TFile::Open("/u/ddas/c4/shiyan_test/run_0003_0001_tree.root");
    if (!file || file->IsZombie())
    {
        std::cerr << "Failed to open file!" << std::endl;
        return;
    }

    //     Set up TTreeReader
    TTreeReader reader("evt", file);
    TTreeReaderArray<LisaCalItem> LisaItem(reader, "LisaCalData");
    TTreeReaderValue<EventHeader> Header(reader, "EventHeader.");

    //     Create a TProfile: avg y for each x bin
    TProfile *prof = new TProfile("prof", "TProfile of Layer 9; Sample point (x) [each point 10ns]; Mean Trace ADC Channel (y)", 400, 0, 400);

    //     Read only the first event
    if (reader.Next())
    {
        std::cout << "Processing event no: " << Header->GetEventno() << std::endl;

        for (const auto &item : LisaItem)
        {
            if (item.Get_layer_id() != 9)
                continue;

            auto trace_x = item.Get_trace_x();
            auto trace_febex = item.Get_trace_febex();

            for (size_t i = 0; i < trace_febex.size(); ++i)
            {
                prof->Fill(trace_x[i], trace_febex[i]);
            }
        }
    }

    //     Draw and fit the TProfile
    TCanvas *c1 = new TCanvas("c1", "TProfile Fit", 800, 600);
    prof->SetMarkerStyle(21);
    prof->SetMarkerColor(kBlue);
    prof->Draw();

    // TF1 *fitFunc = new TF1("fitFunc", step_expo_function_pos, 0, 400, 3);
    TF1 *fitFunc = new TF1("fitFunc", step_expo_function_neg, 0, 400, 3);
    fitFunc->SetParNames("Const", "Expo_Amp", "Tau");
    fitFunc->SetParameters(8150, 300, 1000); // Initial guesses
    fitFunc->SetParLimits(0, 8100, 8200);    // Const
    fitFunc->SetParLimits(1, 200, 4000);     // Expo_Amp
    fitFunc->SetParLimits(2, 10, 2000);      // Tau
    // fitFunc->SetLineColor(kRed);
    // fitFunc->SetLineWidth(2);
    prof->GetXaxis()->SetRangeUser(0, 500); // Full X range
    // prof->SetMinimum(8000);                 // Lower bound of Y-axis
    // prof->SetMaximum(9000);                 // Upper bound of Y-axis
    prof->GetYaxis()->SetTitle("Amplitude (ADC counts)");
    prof->Fit(fitFunc, "Rrob=0.8", "", 0, 400);
    TF1 *fittedFunc = prof->GetFunction("fitFunc");
    fittedFunc->SetLineColor(kRed);
    fittedFunc->SetLineWidth(2);
    fittedFunc->Draw("same");
    c1->Update();
}

// void fit_trace2()
// {
//     Open file
//     TFile *file = TFile::Open("/u/ddas/c4/shiyan_test/test_0003_tree.root");
//     if (!file || file->IsZombie()) {
//         std::cerr << "Failed to open file!" << std::endl;
//         return;
//     }
//     Get the tree
// 	TTree *tree = (TTree*)file->Get("evt");  // <-- change this if needed
// 	if (!tree) {
// 		printf("Tree not found.\n");
// 		return;
// 	}
//
//     TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
//
//     Fill the histogram with a cut
// 	tree->Draw("LisaCalData.trace_febex:LisaCalData.trace_x>>h2()","LisaCalData.layer_id==2&&LisaCalData.xposition==2&&LisaCalData.yposition==2","",1);
//
// 	TH2D *hh = (TH2D*)gROOT->FindObject("h2")->Clone();
//
//     TF1* fexp = new TF1("fexp", step_expo_function, 0, 400, 3);
// 	fexp->SetParNames("Const", "Expo_Amp", "Tau");
// 	fexp->SetParameters(8150, 300, 1000); // initial guesses
// 	fexp->SetParLimits(0, 8100, 8200);
// 	fexp->SetParLimits(1, 200, 400);
// 	fexp->SetParLimits(2, 1000, 2000);
// 	fexp->SetLineColor(kBlue);
// 	hh->Fit(fexp, "Rrob=0.5", "", 0, 1000);
// 	fexp->Draw("same");
//
// }
