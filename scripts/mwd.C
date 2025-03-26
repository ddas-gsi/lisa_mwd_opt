/* :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 * @author: Debajyoti Das
 * Last Updated: 2025-03-24
 *
 * This script is used to calculate the resolution of the MWD energy
 * for a given channel ID.
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

std::vector<double> calcBaselineCorrectedTrace(std::vector<double> &trace_febex)
{
    double sum = 0.0;
    int count = 0;
    double average_baseline = 0.0;
    for (int i = 20; i < 100; i++)
    {
        sum += trace_febex.at(i);
        count++;
    }
    average_baseline = sum / count;

    for (int i = 0; i < trace_febex.size(); i++)
    {
        trace_febex.at(i) = (trace_febex.at(i) - average_baseline) / 8;
    }
    return trace_febex;
}

std::vector<double> calcMWDTrace(std::vector<short> &trace_febex, int k0, int kend, int L, int M, double tau)
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

double calcMWDEnergy(std::vector<double> trace_MWD, int amp_start_idx, int amp_stop_idx,
                     int baseline_start_idx, int baseline_stop_idx, int k0)
{
    amp_start_idx -= k0;
    amp_stop_idx -= k0;
    baseline_start_idx -= k0;
    baseline_stop_idx -= k0;

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

double mwd(int channelID, double smoothing_L, double MWD_length, double MWD_trace_start, double MWD_trace_stop,
           double MWD_amp_start, double MWD_amp_stop, double MWD_baseline_start, double MWD_baseline_stop,
           double sampling, double decay_time)
{
    TFile *file = new TFile("/u/ddas/Lustre/gamma/ddas/LISA/c4_output/run_0075_0001_c4MWD.root");
    TTree *tree = (TTree *)file->Get("evt");
    Int_t entries = tree->GetEntries();
    // Int_t entries = 200;

    TTreeReader reader("evt", file);
    TTreeReaderArray<LisaAnaItem> LisaItem(reader, "LisaAnaData");
    TTreeReaderValue<EventHeader> Header(reader, "EventHeader.");

    // Define a unique histogram name based on parameters
    TString lisa_Energy_histName = Form("lisa_Energy_ch%d_sL%d_MWD%d_start%d_stop%d",
                                        channelID, int(smoothing_L), int(MWD_length), int(MWD_trace_start), int(MWD_trace_stop));

    TH2F *lisa_Trace = new TH2F("lisa_Trace", "lisa_Trace", 2000, 0, 2000, 4000, -500, 100);
    TH2F *lisa_MWD = new TH2F("lisa_MWD", "lisa_MWD", 2000, 0, 1000, 4000, -500, 100);
    // TH1F *lisa_Energy = new TH1F("lisa_Energy", "lisa_Energy", 1000, 10, 500);
    TH1F *lisa_Energy = new TH1F(lisa_Energy_histName, "lisa_Energy", 1000, 10, 500);
    // TH1F *lisa_Energy_cloned = new TH1F("lisa_Energy_cloned", "lisa_Energy", 1000, -10, 500);

    int k = 0;
    while (reader.Next())
    {
        // cout << "Event No.: " << k << endl;
        k++;
        uint evtno = Header->GetEventno();

        double energy_MWD = 0.0;

        for (auto const &LisaItem : LisaItem)
        {
            std::vector<short> trace_febex = LisaItem.Get_trace_febex();
            std::vector<short> trace_x = LisaItem.Get_trace_x();
            int channel_id = LisaItem.Get_channel_id_traces();

            // cout << "Channel ID: " << channel_id << endl;

            if (channel_id == channelID)
            // if (channel_id == 5)
            {
                // Convert params to sample points ----------------------------
                int L = static_cast<int>(smoothing_L / sampling);       // Smoothing time in samples
                int M = static_cast<int>(MWD_length / sampling);        // MWD length in samples
                float tau = decay_time / sampling;                      // Decay constant in samples
                int k0 = static_cast<int>(MWD_trace_start / sampling);  // Start of MWD in samples
                int kend = static_cast<int>(MWD_trace_stop / sampling); // Stop of MWD in samples
                if (kend > trace_febex.size())
                    kend = trace_febex.size(); // If kend out of bound, replace it with trace_febex limit

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

                //  ::: Do MWD Optimisation :::
                std::vector<double> trace_MWD = calcMWDTrace(trace_febex, k0, kend, L, M, tau);
                // cout << "Trace MWD Size: " << trace_MWD.size() << endl;

                energy_MWD = calcMWDEnergy(trace_MWD, amp_start_idx, amp_stop_idx, baseline_start_idx, baseline_stop_idx, k0);
                // cout << "MWD Energy: " << energy_MWD << endl;

                //  ::: Fill Different Histograms :::

                // +++ Fill the RAW/LISA Trace (2D Histogram) +++
                for (UShort_t i = 0; i < trace_febex.size(); i++)
                {
                    // int x = trace_x.at(i);
                    // int y = trace_febex.at(i);
                    lisa_Trace->Fill(i, trace_febex.at(i));
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

        // // for finite number of events ---------------
        // if (k % 10000 == 0)
        // {
        //     cout << setw(5) << setiosflags(ios::fixed) << setprecision(1) << (100. * k) / entries << " % done\r" << flush;
        // }
        // if (k > entries)
        // {
        //     break;
        // }
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
    TF1 *fineGausFit2 = new TF1("fineGausFit2", "gaus", fineMean - 2 * fineSigma, fineMean + 2 * fineSigma);
    fineGausFit2->SetParameter(0, maxEnergy);
    fineGausFit2->SetParameter("Mean", fineMean);
    fineGausFit2->SetParameter("Sigma", fineSigma);
    lisa_Energy->Fit("fineGausFit2", "RQ");
    double fineMean2 = fineGausFit2->GetParameter("Mean");
    double fineSigma2 = fineGausFit2->GetParameter("Sigma");

    double FWHM_2 = fineSigma2 * 2.355;
    double Resolution_2 = (FWHM_2 / fineMean2) * 100;

    // cout << "Max Energy Position: " << maxEnergyPos << " Mean Position: " << fineMean2 << " Resolution : " << Resolution_2 << endl;

    TString histFileName = Form("./histos/hist_%d_%d_%d_%d.root", int(smoothing_L), int(MWD_length), int(MWD_trace_start), int(MWD_trace_stop));
    TFile f(histFileName, "RECREATE");
    lisa_Energy->Write();
    f.Close();

    std::cout << "HISTOGRAM:" << lisa_Energy_histName << std::endl;
    std::cout << "HISTFILE: " << histFileName << std::endl;

    return Resolution_2;

    // +++++++++++++++++++++++++++++++++++++
    // Do TSpectrum stuff here +++++++
    // +++++++++++++++++++++++++++++++++++++
    // TH1F *lisa_Energy_cloned = (TH1F *)lisa_Energy->Clone("clone_for_search");
    // TSpectrum *s = new TSpectrum(1);
    // Int_t nfound = s->Search(lisa_Energy_cloned, 2, "", 0.008);
    // printf("Found %d candidate peaks to fit\n", nfound);

    // // Determine the xy-coordinates of the peak positions
    // Double_t *xpeaks, *ypeaks;
    // xpeaks = s->GetPositionX();
    // ypeaks = s->GetPositionY();
    // for (int p = 0; p < nfound; p++)
    // {
    //     cout << "xPeakPosition: " << xpeaks[p] << " yPeakPosition: " << ypeaks[p] << endl;
    //     Int_t bin = lisa_Energy_cloned->GetXaxis()->FindBin(xpeaks[p]);
    //     Double_t yp = lisa_Energy_cloned->GetBinContent(bin);
    //     cout << "binNbr: " << bin << " yPos: " << yp << endl;
    // }
    // vector<double> Amp, mean, sigma;
    // double par[30];

    // for (int p = 0; p < nfound; p++)
    // {
    //     TF1 *roughFitFunc = new TF1("roughFitFunc", "gaus", xpeaks[p] - 5, xpeaks[p] + 5);
    //     roughFitFunc->SetParameter(0, ypeaks[p]);
    //     roughFitFunc->SetParameter("Mean", xpeaks[p]);
    //     lisa_Energy_cloned->Fit("roughFitFunc", "R");
    //     double roughMean = roughFitFunc->GetParameter("Mean");
    //     double roughSigma = roughFitFunc->GetParameter("Sigma");

    //     TF1 *fineTunedFitFunc = new TF1("fineTunedFitFunc", "gaus", xpeaks[p] - 2 * roughSigma, xpeaks[p] + 2 * roughSigma);
    //     fineTunedFitFunc->SetParameter(0, ypeaks[p]);
    //     fineTunedFitFunc->SetParameter("Mean", roughMean);
    //     fineTunedFitFunc->SetParameter("Sigma", roughSigma);
    //     lisa_Energy_cloned->Fit("fineTunedFitFunc", "R");
    //     Amp.push_back(fineTunedFitFunc->GetParameter(0));
    //     mean.push_back(fineTunedFitFunc->GetParameter("Mean"));
    //     sigma.push_back(fineTunedFitFunc->GetParameter("Sigma"));
    //     fineTunedFitFunc->GetParameters(&par[p * 3]);
    // }

    // vector<double> FWHM, Resolution;
    // for (int p = 0; p < nfound; p++)
    // {
    //     printf("xPeakPosition: %f,\tMean: %f,\tSigma: %f,\tAmplitude: %f\n", xpeaks[p], mean[p], sigma[p], Amp[p]);
    //     FWHM.push_back(sigma[p] * 2.355);
    //     Resolution.push_back((FWHM[p] / mean[p]) * 100);

    //     printf("xPeakPosition: %f,\t FWHM: %f,\t Resolution: %f %%\n", xpeaks[p], FWHM[p], Resolution[p]);
    // }

    // // Record the end time
    // // auto scriptEndTime = std::chrono::high_resolution_clock::now();

    // // Calculate the duration
    // // std::chrono::duration<double> elapsed = scriptEndTime - scriptStartTime;

    // // Output the runtime
    // std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;

    // //  ::: Draw Histograms :::
    // c1->cd(1);
    // lisa_Trace->Draw("colz");
    // c1->cd(2);
    // lisa_MWD->Draw("colz");
    // c1->cd(3);
    // lisa_Energy->Draw();
    // c1->Update();

    // c2->cd();
    // lisa_Energy_cloned->Draw();
    // c2->Update();
}
