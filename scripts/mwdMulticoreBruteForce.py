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
from multiprocessing import Pool
import subprocess
import re

# Parameter ranges
smoothing_L_range = np.linspace(200, 600, 3)              # Smoothing_L[ns]          
MWD_length_range = np.linspace(500, 1500, 2)             # Trapez_moving_window_length[ns]   
MWD_trace_start_range = np.linspace(1800, 2000, 3)        # Trapez_sample_window_0[ns]
MWD_trace_stop_range = np.linspace(4000, 4500, 3)         # Trapez_sample_window_1[ns]

# Define CPU usage for parallel processing
cpu_usage = 0.5                                           # Use 50% of available CPU cores

# ::: Constant Values ::: 
# Channel ID of Diamond
channelID = 5

sampling = 10
decay_time = 30500

use_cores = int(os.cpu_count()*cpu_usage) 


def worker_function(params):
    smoothing_L, MWD_length, MWD_trace_start, MWD_trace_stop, MWD_amp_start, MWD_amp_stop, MWD_baseline_start, MWD_baseline_stop = params
    root_command = ["root", "-l", "-q", "-b", f"mwd.C({channelID},{smoothing_L},{MWD_length},{MWD_trace_start},{MWD_trace_stop}, {MWD_amp_start}, {MWD_amp_stop}, {MWD_baseline_start}, {MWD_baseline_stop}, {sampling}, {decay_time})"]
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

    # histMatch = re.search(r"HISTOGRAM:\s*(\S+)", result.stdout)
    # histName = histMatch.group(1) if histMatch else None

    # histFileMatch = re.search(r"HISTFILE:\s*(\S+)", result.stdout)
    # histFile = histFileMatch.group(1) if histFileMatch else None

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

    # return resolution, params, histName, histFile
    return resolution, params, lisaEnergyHistName, lisaEnergyHistFile, lisaTraceHistName, lisaMWDHistName, lisaTraceHistFile
    
    
    


def mwd_multicore():

    start_time = datetime.now()

    param_combinations = list(product(smoothing_L_range, MWD_length_range, MWD_trace_start_range, MWD_trace_stop_range))

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

    with Pool(processes=use_cores) as pool:
        results = pool.map(worker_function, filtered_param_combinations)

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

    # Separate results into two lists: 
    # - `positive_res_results` for histograms with resolution > 0
    # - `non_positive_res_results` for histograms with resolution <= 0
    positive_res_results = [res for res in results if res[0] is not None and res[0] > 0]
    non_positive_res_results = [res for res in results if res[0] is not None and res[0] <= 0]

    # Merge the two lists, ensuring that positive resolutions are written first
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






# =========================================================================================

# import ROOT
# import os
# import numpy as np
# import scipy as sc
# import matplotlib.pyplot as plt
# from itertools import product
# from datetime import datetime
# from multiprocessing import Pool
# from functools import partial
# import subprocess
# import re
# from datetime import datetime


# def chi2(data, model):
#     return np.sum((data - model)**2)

# # Parameter ranges
# smoothing_L_range = np.linspace(200, 600, 2)              # Smoothing_L[ns]          # np.linspace(200, 600, 10, dtype=int) 
# MWD_length_range = np.linspace(500, 1500, 1)              # Trapez_moving_window_length[ns]   
# MWD_trace_start_range = np.linspace(2000, 6200, 1)        # Trapez_sample_window_0[ns]
# MWD_trace_stop_range = np.linspace(4200, 4500, 1)         # Trapez_sample_window_1[ns]

# param_combinations = list(product(smoothing_L_range, MWD_length_range, MWD_trace_start_range, MWD_trace_stop_range))

# Res_dict = {}

# # Run the script for all parameter sets
# for params in param_combinations:

#     # Channel ID
#     channelID = 5

#     smoothing_L, MWD_length, MWD_trace_start, MWD_trace_stop = params
#     root_command = ["root", "-l", "-q", "-b", "-n", f"mwd.C({channelID},{smoothing_L},{MWD_length},{MWD_trace_start},{MWD_trace_stop})"]
#     print(f"mwd script is running for parameters {params}...")
    
#     start_time = datetime.now()
#     # Run the command without `check=True` to prevent errors from stopping execution
#     result = subprocess.run(root_command, capture_output=True, text=True)
#     # print("STDOUT:", result.stdout)  # From stdout, we can check for which parameters we get this resolution... Important in case of parallel processing
#     # print("STDERR:", result.stderr)  # Just for debugging

#     # Extract Resolution from output
#     match = re.search(r"\(double\)\s*([\d.eE+-]+)", result.stdout)
#     if match:
#         Resolution = float(match.group(1))
#         print(f"Parameters {params} give Resolution = {Resolution}")
#         Res_dict[Resolution] = params
#     else:
#         print(f"Error: Resolution not found in output for parameters {params}.")

#     stop_time = datetime.now()
#     print('Duration:', stop_time - start_time)


# # Find the parameters corresponding to the lowest resolution
# opt_res = min(Res_dict.keys())  # Find the lowest resolution (key)
# params_for_opt_res = Res_dict[opt_res]  # Get the parameters for that resolution

# print(f"The lowest resolution is {opt_res}, with parameters: {params_for_opt_res}")

# =========================================================================================

# start_time = datetime.now()
# # Channel ID
# channelID = 10
# # Parameters
# smoothing_L = 400
# MWD_length = 1000
# MWD_trace_start = 2000  
# MWD_trace_stop = 4200

# # parameter sets
# param_sets = [
#     (2,400,1000,2000,4200),
#     (5,400,1000,2000,4200),
#     (10,400,1000,2000,4200)
#     ]
# processes = []
# for params in param_sets:
#     channelID, smoothing_L, MWD_length, MWD_trace_start, MWD_trace_stop = params
#     root_command = ["root", "-l", "-q", "-b", f"mwd.C({channelID},{smoothing_L},{MWD_length},{MWD_trace_start},{MWD_trace_stop})"]
#     print("mwd script is running ...")
    
#     # Start subprocess in parallel
#     process = subprocess.Popen(root_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
#     processes.append((process,params))
    

# # Collect results
# for process, params in processes:
#     stdout, stderr = process.communicate()
#     match = re.search(r"\(double\)\s*([\d.eE+-]+)", stdout)
#     if match:
#         Resolution = float(match.group(1))
#         print(f"Parameters {params} give Resolution = {Resolution}")
#     else:
#         print(f"Error: Resolution not found in output for parameters {params}.")

# stop_time = datetime.now()
# print('Duration:', stop_time - start_time)

# =========================================================================================