# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# @author: Debajyoti Das                  
# Last Updated: 2025-03-24                
#                                         
# 
# This script runs MWD calculations in multiple parallel processes
# to speed up the brute-force search of the best parameter set 
# for optimal energy resolution of LISA Diamond detectors
# 
# The lowest resolution is 1.3761951, with parameters: (200.0, 700.0, 2000.0, 4500.0, 2825.0, 3075.0, 2125.0, 2375.0)
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


import ROOT
import os
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from itertools import product
from datetime import datetime
from multiprocessing import Pool, Manager
import subprocess
import re
import load_mwd 
import threading
import time




MWD_SETTINGS = load_mwd.load_mwd_settings("./set_mwd.txt")

INPUT_FILE = str(MWD_SETTINGS.get("INPUT_FILE", "Not Found"))
# print(INPUT_FILE)
HISTOGRAM_FILE_PATH = str(MWD_SETTINGS.get("HISTOGRAM_FILE_PATH", "Not Found"))
# print(HISTOGRAM_FILE_PATH)
channelID = MWD_SETTINGS.get("channelID", "Not Found")
Sampling = MWD_SETTINGS.get("Sampling", "Not Found")
Decaytime_ch = MWD_SETTINGS.get("Decaytime_ch", "Not Found")
FIT_RANGE_PAR = MWD_SETTINGS.get("FIT_RANGE_PAR", "Not Found")

smoothing_L_Axis = np.linspace(MWD_SETTINGS.get("Smoothing_L_Range.start", "Not Found"), 
                               MWD_SETTINGS.get("Smoothing_L_Range.stop", "Not Found"), 
                               MWD_SETTINGS.get("Smoothing_L_Range.split", "Not Found"))

MWD_length_Axis = np.linspace(MWD_SETTINGS.get("MWD_Length_Range.start", "Not Found"), 
                              MWD_SETTINGS.get("MWD_Length_Range.stop", "Not Found"), 
                              MWD_SETTINGS.get("MWD_Length_Range.split", "Not Found"))

MWD_trace_start_Axis = np.linspace(MWD_SETTINGS.get("MWD_Trace_Start_Range.start", "Not Found"), 
                                   MWD_SETTINGS.get("MWD_Trace_Start_Range.stop", "Not Found"), 
                                   MWD_SETTINGS.get("MWD_Trace_Start_Range.split", "Not Found"))

MWD_trace_stop_Axis = np.linspace(MWD_SETTINGS.get("MWD_Trace_Stop_Range.start", "Not Found"), 
                                  MWD_SETTINGS.get("MWD_Trace_Stop_Range.stop", "Not Found"), 
                                  MWD_SETTINGS.get("MWD_Trace_Stop_Range.split", "Not Found"))


CPU_Usage = MWD_SETTINGS.get("CPU_Usage", "Not Found")
    
use_cores = int(os.cpu_count()*CPU_Usage) 

def ensure_directory_exists(dir_path):
    """Check if a directory exists, and create it if it doesn't."""
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
        print(f"Created directory: {dir_path}")
    else:
        print(f"Directory already exists: {dir_path}")


def monitor_progress(filtered_param_combinations):
        while True:
            histo_count = int(os.popen("ls -lhrt mwd_histos/ | wc -l").read().strip())
            params_computed = int(histo_count/2)
            # print("\n")
            print(f"Parameters remaining: {len(filtered_param_combinations) - params_computed} ...", end='\r', flush=True)
            time.sleep(5)  # Adjust the interval as needed


def worker_function(params):
    smoothing_L, MWD_length, MWD_trace_start, MWD_trace_stop, MWD_amp_start, MWD_amp_stop, MWD_baseline_start, MWD_baseline_stop = params
    root_command = ["root", "-l", "-q", "-b", f"mwd.C({channelID},{smoothing_L},{MWD_length},{MWD_trace_start},{MWD_trace_stop},{MWD_amp_start},{MWD_amp_stop},{MWD_baseline_start},{MWD_baseline_stop},{Sampling},{Decaytime_ch},{FIT_RANGE_PAR})"]
    print(f"MWD running for parameters {params}...")
    start_time = datetime.now()
    result = subprocess.run(root_command, capture_output=True, text=True)

    print("STDOUT:", result.stdout)  # for debugging
    print("STDERR:", result.stderr)  # for debugging

    match = re.search(r"\(double\)\s*([\d.eE+-]+)", result.stdout)
    resolution = None
    if match:
        resolution = float(match.group(1))
        print(f"Parameters {params} -> Resolution = {resolution}")
    else:
        print(f"Error: Resolution not found in output for parameters {params}")

    stop_time = datetime.now()
    print('Duration:', stop_time - start_time)

    lisaEnergyMatch = re.search(r"LISA_ENERGY:\s*(\S+)", result.stdout)
    lisaEnergyHistName = lisaEnergyMatch.group(1) if lisaEnergyMatch else None

    lisaEnergyFileMatch = re.search(r"LISA_ENERGY_FILE:\s*(\S+)", result.stdout)
    lisaEnergyHistFile = lisaEnergyFileMatch.group(1) if lisaEnergyFileMatch else None

    lisaTraceMatch = re.search(r"LISA_TRACE:\s*(\S+)", result.stdout)
    lisaTraceHistName = lisaTraceMatch.group(1) if lisaTraceMatch else None

    lisaMWDMatch = re.search(r"LISA_MWD:\s*(\S+)", result.stdout)
    lisaMWDHistName = lisaMWDMatch.group(1) if lisaMWDMatch else None

    lisaTraceFileMatch = re.search(r"LISA_TRACE_FILE:\s*(\S+)", result.stdout)
    lisaTraceHistFile = lisaTraceFileMatch.group(1) if lisaTraceFileMatch else None

    # # Update progress counter
    # with progress_counter.get_lock():
    #     progress_counter.value -= 1
    #     print(f"Remaining parameter sets: {progress_counter.value}", end='\r', flush=True)

    return resolution, params, lisaEnergyHistName, lisaEnergyHistFile, lisaTraceHistName, lisaMWDHistName, lisaTraceHistFile


def mwd_multicore():

    start_time = datetime.now()
    ensure_directory_exists(HISTOGRAM_FILE_PATH)

    param_combinations = list(product(smoothing_L_Axis, MWD_length_Axis, MWD_trace_start_Axis, MWD_trace_stop_Axis))

    print(f"Total sets of parameters {len(param_combinations)} ...")

    filtered_param_combinations = []

    for params in param_combinations:
        smoothing_L, MWD_length, MWD_trace_start, MWD_trace_stop = params

        if (smoothing_L<MWD_length) and ((smoothing_L+MWD_length)<(MWD_trace_stop-MWD_trace_start)) and (MWD_trace_start<MWD_trace_stop) and (MWD_trace_start>=0):
            # Calculate the MWD Energy Parameters
            MWD_amp_start = (MWD_trace_start + (3 / 2) * MWD_length - (1 / 2) * smoothing_L) - ((MWD_length - smoothing_L) / 4); # Trapez_amp_calc_window_0[ns]
            MWD_amp_stop = (MWD_trace_start + (3 / 2) * MWD_length - (1 / 2) * smoothing_L) + ((MWD_length - smoothing_L) / 4);  # Trapez_amp_calc_window_1[ns]
            MWD_baseline_start = MWD_trace_start + ((MWD_length - smoothing_L) / 2) - ((MWD_length - smoothing_L) / 4);          # Trapez_baseline_window_0[ns]
            MWD_baseline_stop = MWD_trace_start + ((MWD_length - smoothing_L) / 2) + ((MWD_length - smoothing_L) / 4);           # Trapez_baseline_window_1[ns]

            if (MWD_amp_start>(MWD_trace_start+MWD_length)) and (MWD_amp_stop<(MWD_trace_start+2*MWD_length-smoothing_L)) and (MWD_baseline_start>MWD_trace_start) and (MWD_baseline_stop<(MWD_trace_start+MWD_length-smoothing_L)):
                new_params = params + (MWD_amp_start, MWD_amp_stop, MWD_baseline_start, MWD_baseline_stop)
                filtered_param_combinations.append(new_params)

    print(f"Running {len(filtered_param_combinations)} sets of parameter on {use_cores} cores in parallel...")

    

    # Start monitoring in a separate thread
    monitor_thread = threading.Thread(target=monitor_progress, daemon=True, args=(filtered_param_combinations,))
    monitor_thread.start()

    with Pool(processes=use_cores) as pool:
        results = pool.map(worker_function, filtered_param_combinations)

    # global progress_counter
    # with Manager() as manager:
    #     progress_counter = manager.Value('i', len(filtered_param_combinations))  # Initialize progress counter
    #     with Pool(processes=use_cores) as pool:
    #         results = pool.map(worker_function, filtered_param_combinations)

    Res_dict = {}
    for resolution, params, _, _, _, _, _ in results:
        if resolution is not None:
            if resolution > 0:
                Res_dict[resolution] = params
    
    sorted_res = sorted(Res_dict.items())

    output_filename = "resolution_results.txt"
    with open(output_filename, "w") as file:
        file.write("Resolution\tParameters\n")
        for res, params in sorted_res:
            file.write(f"{res:.6f}\t{params}\n")

    if sorted_res:
        opt_res, params_for_opt_res = sorted_res[0]
        print(f"The lowest resolution is {opt_res}, with parameters: {params_for_opt_res}")
        print(f"All resolutions and parameters saved to {output_filename}")
    else:
        print("No valid resolutions found.")
                
    # Create output ROOT files
    energy_hist_file = ROOT.TFile("merged_lisa_Energy_histos.root", "RECREATE")
    trace_hist_file = ROOT.TFile("merged_lisa_Trace_histos.root", "RECREATE")

    hList = ROOT.TList()

    positive_res_results = [res for res in results if res[0] is not None and res[0] > 0]
    non_positive_res_results = [res for res in results if res[0] is not None and res[0] <= 0]

    sorted_results = sorted(positive_res_results) + sorted(non_positive_res_results)

    # Flag to track if lisaTraceHist has been written
    trace_hist_written = False

    # Process histograms in the required order
    for resolution, params, lisaEnergyHistName, lisaEnergyHistFile, lisaTraceHistName, lisaMWDHistName, lisaTraceHistFile in sorted_results:
        if lisaEnergyHistFile:
            temp_Energy_file = ROOT.TFile(lisaEnergyHistFile, "READ")
            energyHist = temp_Energy_file.Get(lisaEnergyHistName)  
            if energyHist:
                energy_hist_file.cd()
                energyHist.Write()
                hList.Add(energyHist)
            temp_Energy_file.Close()

        if lisaTraceHistFile:
            temp_Trace_file = ROOT.TFile(lisaTraceHistFile, "READ")
        
            if not trace_hist_written:  # Only write lisaTraceHist once
                lisaTraceHist = temp_Trace_file.Get(lisaTraceHistName)
                if lisaTraceHist:
                    trace_hist_file.cd()
                    lisaTraceHist.Write()
                    trace_hist_written = True  # Mark as written

            lisaMWDHist = temp_Trace_file.Get(lisaMWDHistName)
            if lisaMWDHist:
                trace_hist_file.cd()
                lisaMWDHist.Write()
        
            temp_Trace_file.Close()

                
    # if hList:
    #     output_file.cd()
    #     hList.Write("Energy_histos", ROOT.TObject.kSingleKey)

    energy_hist_file.Close()
    print("All LISA Energy histograms saved to merged_lisa_Energy_histos.root")

    trace_hist_file.Close()
    print("All LISA Trace histograms saved to merged_lisa_Trace_histos.root")
    
    stop_time = datetime.now()
    print('Total Script RunTime:', stop_time - start_time)

  

if __name__ == '__main__':
    mwd_multicore()


