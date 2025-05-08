/* :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 * @author: Debajyoti Das
 * Last Updated: 2025-05-07
 *
 * This script is used to calculate the resolution of the MWD energy
 * for a given layer and x-y position of the diamond grid
 *
 *
 * ::: Before running this script, make sure to initialise C4 :::
 * >> cd c4/build/
 * >> make -j
 * >> . ./config.sh
 * >> cd ../c4Root/macros/lisa/trace_analysis/
 *
 *
 * ::: Run this script using the following command :::
 * >> root -l -q -b "mwd.C(params)"
 ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <sys/time.h>
#include <signal.h>
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <iostream>
#include <cmath>
// #include <chrono>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
// #include "TSpectrum.h"

using namespace std;

// ++++++++++++++++++++++++++++++++++++++++++++++++++++
//     :::::  Define Different MWD Functions  :::::
// ++++++++++++++++++++++++++++++++++++++++++++++++++++

std::vector<float> calcBaselineCorrectedTrace(std::vector<float> &trace_febex)
{
    double sum = 0.0;
    int count = 0;
    double average_baseline = 0.0;
    for (int i = 0; i < 10; i++)
    {
        sum += trace_febex.at(i);
        count++;
    }
    average_baseline = sum / count;

    for (int i = 0; i < trace_febex.size(); i++)
    {
        trace_febex.at(i) = (trace_febex.at(i) - average_baseline);
        // trace_febex.at(i) = (trace_febex.at(i) - average_baseline) / 8.0; // dividing by 8 makes ADC channel to mV conversion. 8000 ADC = 1000 mV
    }
    return trace_febex;
}

std::vector<double> calcMWDTrace(std::vector<float> &trace_febex, int k0, int kend, int L, int M, double tau)
{
    std::vector<double> trace_MWD; // Output vector for MWD values
    // Loop through the trace
    for (int k = k0; k < kend; ++k)
    {
        double DM = 0.0;
        for (int j = k - L; j <= k - 1; ++j)
        {
            if (j < 1)
                continue;      // Skip if out-of-bounds
            double sum0 = 0.0; // Initialize sum for baseline subtraction
            for (int i = j - M; i <= j - 1; ++i)
            {
                if (i < 0)
                    continue;              // Skip if out-of-bounds
                sum0 += trace_febex.at(i); // Sum over the moving window
            }
            DM += trace_febex.at(j) - trace_febex.at(j - M) + sum0 / tau;
        }
        // Calculate MWD value and index
        double mwd_value = DM / L; // Average over the rising time
        trace_MWD.push_back(mwd_value);
    }
    return trace_MWD;
}

// std::vector<double> calcMWDTrace_noShift(std::vector<float> &trace_febex, int k0, int kend, int L, int M, double tau)
// {
//     std::vector<double> trace_MWD; // Output vector for MWD values
//     // Loop through the trace
//     for (int k = 0; k < kend; ++k)
//     {
//         double DM = 0.0;
//         if (k >= k0)
//         {
//             for (int j = k - L; j <= k - 1; ++j)
//             {
//                 if (j < 1)
//                     continue;      // Skip if out-of-bounds
//                 double sum0 = 0.0; // Initialize sum for baseline subtraction
//                 for (int i = j - M; i <= j - 1; ++i)
//                 {
//                     if (i < 0)
//                         continue;              // Skip if out-of-bounds
//                     sum0 += trace_febex.at(i); // Sum over the moving window
//                 }
//                 DM += trace_febex.at(j) - trace_febex.at(j - M) + sum0 / tau;
//             }
//         }

//         // Calculate MWD value and index
//         double mwd_value = DM / L; // Average over the rising time
//         trace_MWD.push_back(mwd_value);
//     }
//     return trace_MWD;
// }

std::vector<double> calcMWDTrace_noShift(std::vector<float> &trace_febex, int k0, int kend, int L, int M, double tau)
{
    std::vector<double> trace_MWD; // Output vector for MWD values
    // Loop through the trace
    for (int i = 0; i < kend; ++i)
    {
        double DM = 0.0;
        if (i >= k0)
        {
            for (int j = i - L; j <= i - 1; ++j)
            {
                if (j < 1)
                    continue;      // Skip if out-of-bounds
                double sum0 = 0.0; // Initialize sum for baseline subtraction
                for (int k = j - M; k <= j - 1; ++k)
                {
                    if (k < 0)
                        continue;              // Skip if out-of-bounds
                    sum0 += trace_febex.at(k); // Sum over the moving window
                }
                DM += trace_febex.at(j) - trace_febex.at(j - M) + sum0 / tau;
            }
        }

        // Calculate MWD value and index
        double mwd_value = DM / L; // Average over the rising time
        trace_MWD.push_back(mwd_value);
    }
    return trace_MWD;
}

double calcMWDEnergy(std::vector<double> trace_MWD, int amp_start_idx, int amp_stop_idx,
                     int baseline_start_idx, int baseline_stop_idx, int k0)
{
    // amp_start_idx -= k0; // for shifting the MWD trace to x=0, but for Shiyan Pulsar it's not needed
    // amp_stop_idx -= k0;
    // baseline_start_idx -= k0;
    // baseline_stop_idx -= k0;

    // cout << "Adjusted Indices: amp_start_idx = " << amp_start_idx
    //      << ", amp_stop_idx = " << amp_stop_idx
    //      << ", baseline_start_idx = " << baseline_start_idx
    //      << ", baseline_stop_idx = " << baseline_stop_idx
    //      << ", k0 = " << k0 << std::endl;

    // Boundary check
    if (amp_start_idx < 0 || amp_stop_idx > trace_MWD.size() ||
        baseline_start_idx < 0 || baseline_stop_idx > trace_MWD.size())
    {
        std::cerr << "Error: Index out of bounds in calcMWDEnergy." << std::endl;
        return 0.0;
    }

    double baseline_sum = 0.0, baseline_avg = 0.0, energy_sum = 0.0;
    int amp_count = 0, baseline_count = 0;

    for (int i = amp_start_idx; i < amp_stop_idx; ++i)
    {
        energy_sum += trace_MWD.at(i);
        amp_count++;
    }
    double energy_avg = energy_sum / amp_count;

    for (int i = baseline_start_idx; i < baseline_stop_idx; ++i)
    {
        baseline_sum += trace_MWD.at(i);
        baseline_count++;
    }
    baseline_avg = baseline_sum / baseline_count;

    return abs(energy_avg - baseline_avg);
}

// +++++++++++++++++++++++++++++++++++++
//     :::::  Main function  :::::
// +++++++++++++++++++++++++++++++++++++

double mwd_shiyan(int channelID, double smoothing_L, double MWD_length, double MWD_trace_start, double MWD_trace_stop,
                  double MWD_amp_start, double MWD_amp_stop, double MWD_baseline_start, double MWD_baseline_stop,
                  double sampling, double decay_time, double FIT_RANGE_PAR)
{
    string INPUT_FILE = "/u/ddas/c4/shiyan_test/test_0003_tree.root";
    string HISTOGRAM_FILE_PATH = "/u/ddas/c4/c4Root/macros/lisa/trace_analysis/mwd_histos/";
    // string HISTOGRAM_FILE_PATH = "/u/ddas/c4/c4Root/macros/lisa/trace_analysis/mwd_histos2/";
    // TFile *file = new TFile("/u/ddas/Lustre/gamma/ddas/LISA/c4_output/run_0075_0001_c4MWD.root");
    TString inputFile = INPUT_FILE;
    TFile *file = new TFile(inputFile);
    TTree *tree = (TTree *)file->Get("evt");
    Int_t entries = tree->GetEntries();
    // Int_t entries = 200;

    TTreeReader reader("evt", file);
    TTreeReaderArray<LisaCalItem> LisaItem(reader, "LisaCalData");
    TTreeReaderValue<EventHeader> Header(reader, "EventHeader.");

    // Define a unique histogram name based on parameters
    TString lisa_trace_histName = Form("lisa_Trace_L_%d_M_%d_k0_%d_kend_%d",
                                       int(smoothing_L), int(MWD_length), int(MWD_trace_start), int(MWD_trace_stop));

    TString lisa_MWD_histName = Form("lisa_MWD_L_%d_M_%d_k0_%d_kend_%d",
                                     int(smoothing_L), int(MWD_length), int(MWD_trace_start), int(MWD_trace_stop));

    TString lisa_Energy_histName = Form("lisa_Energy_L_%d_M_%d_k0_%d_kend_%d",
                                        int(smoothing_L), int(MWD_length), int(MWD_trace_start), int(MWD_trace_stop));

    // TString lisa_Energy_histTitle = Form("lisa_Energy >> L:%d M:%d k0:%d kend:%d a:%d b:%d c:%d d:%d",
    //                                      int(smoothing_L), int(MWD_length), int(MWD_trace_start), int(MWD_trace_stop) int(MWD_amp_start), int(MWD_amp_stop), int(MWD_baseline_start), int(MWD_baseline_stop));

    TH2F *lisa_Trace = new TH2F("lisa_Trace", "lisa_Trace", 500, 0, 500, 1000, 0, -1000); // change y-range depending on positive or negative polarity of Trace

    TH2F *lisa_MWD = new TH2F(lisa_MWD_histName, "lisa_MWD", 500, 0, 500, 1000, 0, -1000); // change y-range depending on positive or negative polarity of Trace

    TH1F *lisa_Energy = new TH1F(lisa_Energy_histName, "lisa_Energy", 1000, 10, 1010);

    int k = 0;
    while (reader.Next())
    {
        // cout << "Event No.: " << k << endl;
        k++;
        uint evtno = Header->GetEventno();

        double energy_MWD = 0.0;

        for (auto const &LisaItem : LisaItem)
        {
            std::vector<float> trace_febex = LisaItem.Get_trace_febex();
            std::vector<int16_t> trace_x = LisaItem.Get_trace_x();
            // int channel_id = LisaItem.Get_channel_id_traces();

            // cout << "Channel ID: " << channel_id << endl;

            if (LisaItem.Get_layer_id() == 2 && LisaItem.Get_xposition() == 2 && LisaItem.Get_yposition() == 2)
            {
                // Convert params to sample points ----------------------------
                int L = static_cast<int>(smoothing_L / sampling);       // Smoothing time in samples
                int M = static_cast<int>(MWD_length / sampling);        // MWD length in samples
                float tau = decay_time / sampling;                      // Decay constant in samples
                int k0 = static_cast<int>(MWD_trace_start / sampling);  // Start of MWD in samples
                int kend = static_cast<int>(MWD_trace_stop / sampling); // Stop of MWD in samples
                if (kend > trace_febex.size())
                {
                    kend = trace_febex.size(); // If kend out of bound, replace it with trace_febex limit
                }

                // cout << "L:" << L << " M: " << M << " tau: " << tau << " k0: " << k0 << " kend: " << kend << endl;

                //    ::: Convert MWD Energy Parameters to sample points :::
                int amp_start_idx = static_cast<int>(MWD_amp_start / sampling);
                int amp_stop_idx = static_cast<int>(MWD_amp_stop / sampling);
                int baseline_start_idx = static_cast<int>(MWD_baseline_start / sampling);
                int baseline_stop_idx = static_cast<int>(MWD_baseline_stop / sampling);

                // std::cout << "Adjusted Indices: amp_start_idx = " << amp_start_idx - k0
                //           << ", amp_stop_idx = " << amp_stop_idx - k0
                //           << ", baseline_start_idx = " << baseline_start_idx - k0
                //           << ", baseline_stop_idx = " << baseline_stop_idx - k0 << std::endl;

                // ------------------------------------------------------------

                // // ::: Baseline correction :::
                std::vector<float> trace_febex_baselineCorrected = calcBaselineCorrectedTrace(trace_febex);
                // // cout << "Trace Size: " << trace_febex_corrected.size() << endl;

                // //  ::: Do MWD Optimisation :::
                // std::vector<double> trace_MWD = calcMWDTrace(trace_febex_baselineCorrected, k0, kend, L, M, tau); // for run_0075_0001_c4MWD.root it's already baseline corrected. No need to do it again
                std::vector<double> trace_MWD = calcMWDTrace_noShift(trace_febex_baselineCorrected, k0, kend, L, M, tau);
                // std::vector<double> trace_MWD = calcMWDTrace(trace_febex, k0, kend, L, M, tau);
                // // cout << "Trace MWD Size: " << trace_MWD.size() << endl;

                energy_MWD = calcMWDEnergy(trace_MWD, amp_start_idx, amp_stop_idx, baseline_start_idx, baseline_stop_idx, k0);
                // // cout << "MWD Energy: " << energy_MWD << endl;

                //  ::: Fill Different Histograms :::

                // +++ Fill the RAW/LISA Trace (2D Histogram) +++
                for (size_t i = 0; i < trace_febex.size(); i++)
                {
                    // lisa_Trace->Fill(i, trace_febex.at(i));
                    lisa_Trace->Fill(i, trace_febex_baselineCorrected.at(i));
                }

                // +++ Fill the MWD Trace (2D Histogram) +++
                for (int i = 0; i < trace_MWD.size(); i++)
                {
                    // cout << trace_MWD.at(i) << endl;
                    lisa_MWD->Fill(i, trace_MWD.at(i));
                }

                // +++ Fill the Energy Histogram +++
                lisa_Energy->Fill(energy_MWD);
            }
        }
    }

    // +++++++++++++++++++++++++++++++++++++++++++
    //    ::: Calculate Resolution and FWHM :::
    // +++++++++++++++++++++++++++++++++++++++++++

    int maxBin = lisa_Energy->GetMaximumBin();
    double maxEnergy = lisa_Energy->GetBinContent(maxBin);
    double maxEnergyPos = lisa_Energy->GetBinCenter(maxBin);

    // Rough Gaussian Fit
    TF1 *roughGausFit = new TF1("roughGausFit", "gaus", maxEnergyPos - 5, maxEnergyPos + 5);
    roughGausFit->SetParameter(0, maxEnergy);
    roughGausFit->SetParameter("Mean", maxEnergyPos);
    lisa_Energy->Fit("roughGausFit", "RQ");
    double roughMean = roughGausFit->GetParameter("Mean");
    double roughSigma = roughGausFit->GetParameter("Sigma");

    // First fine Gaussian Fit
    TF1 *fineGausFit = new TF1("fineGausFit", "gaus", maxEnergyPos - 2 * roughSigma, maxEnergyPos + 2 * roughSigma);
    fineGausFit->SetParameter(0, maxEnergy);
    fineGausFit->SetParameter("Mean", roughMean);
    fineGausFit->SetParameter("Sigma", roughSigma);
    lisa_Energy->Fit("fineGausFit", "RQ");
    double fineMean = fineGausFit->GetParameter("Mean");
    double fineSigma = fineGausFit->GetParameter("Sigma");

    double FWHM = fineSigma * 2.355;
    double Resolution = (FWHM / fineMean) * 100;

    // Second fine Gaussian Fit
    TF1 *fineGausFit2 = new TF1("fineGausFit2", "gaus", fineMean - FIT_RANGE_PAR * fineSigma, fineMean + FIT_RANGE_PAR * fineSigma);
    fineGausFit2->SetParameter(0, maxEnergy);
    fineGausFit2->SetParameter("Mean", fineMean);
    fineGausFit2->SetParameter("Sigma", fineSigma);
    lisa_Energy->Fit("fineGausFit2", "RQ");
    double fineMean2 = fineGausFit2->GetParameter("Mean");
    double fineSigma2 = fineGausFit2->GetParameter("Sigma");

    double FWHM_2 = fineSigma2 * 2.355;
    double Resolution_2 = (FWHM_2 / fineMean2) * 100;

    // cout << "Max Energy Position: " << maxEnergyPos << " Mean Position: " << fineMean2 << " Resolution : " << Resolution_2 << endl;

    TString histFileName = Form("%shist_%d_%d_%d_%d.root", HISTOGRAM_FILE_PATH.c_str(), int(smoothing_L), int(MWD_length), int(MWD_trace_start), int(MWD_trace_stop));
    TFile f(histFileName, "RECREATE");
    lisa_Energy->SetTitle(Form("lisa_Energy >> Res:%.6f L:%d M:%d k0:%d kend:%d a:%d b:%d c:%d d:%d",
                               Resolution_2, int(smoothing_L), int(MWD_length), int(MWD_trace_start), int(MWD_trace_stop),
                               int(MWD_amp_start), int(MWD_amp_stop), int(MWD_baseline_start), int(MWD_baseline_stop)));
    lisa_Energy->Write();
    f.Close();

    std::cout << "LISA_ENERGY:" << lisa_Energy_histName << std::endl;
    std::cout << "LISA_ENERGY_FILE: " << histFileName << std::endl;

    TString lisa2DHistFileName = Form("%slisaTrace_%d_%d_%d_%d.root", HISTOGRAM_FILE_PATH.c_str(), int(smoothing_L), int(MWD_length), int(MWD_trace_start), int(MWD_trace_stop));
    TFile f2(lisa2DHistFileName, "RECREATE");
    lisa_Trace->Write();
    lisa_MWD->Write();
    f2.Close();

    // lisa_Trace->Draw("colz");
    // lisa_MWD->Draw("colz");
    // lisa_Energy->Draw();

    // std::cout << "LISA_TRACE:" << lisa_trace_histName << std::endl;
    std::cout << "LISA_TRACE:" << "lisa_Trace" << std::endl;
    std::cout << "LISA_MWD:" << lisa_MWD_histName << std::endl;
    std::cout << "LISA_TRACE_FILE: " << lisa2DHistFileName << std::endl;

    return Resolution_2;
}
