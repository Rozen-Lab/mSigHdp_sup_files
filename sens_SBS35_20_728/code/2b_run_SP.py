"""
Run SigProfilerExtractor, using parameters specified in Supp Note 1
of Mishu's SigProfilerExtractor pre-print.

Some parameters may be different from the default.

SigProfilerExtractor accepts ICAMS-formatted catalogs in CSV format.

However, we still converted ICAMS-formatted spectra CSV catalogs 
into SigPro formatted TSV catalogs.

The input and output directories is required to be provided 
as arguments from stdin.
"""
# Please run this script from the top-level directory, mSigHdp_paper_sup_files_x1

###############################################################################
#%% Cell 1: load prerequisites
###############################################################################
from pathlib import Path         # Make directory recursively
import gc                        # Garbage collector
import os
import pandas as pd              # DataFrame I/O
import pickle                    # Serialze object and save
import random                    # Required to designate random seed
import sys
import tracemalloc               # RAM profiler



###############################################################################
# From Cell 2 and onwards:
# Must add the statement if __name__=="__main__".
# Otherwise SigProfilerExtractor will raise a RunTimeError
# see GitHub repo issue #88:
# https://github.com/AlexandrovLab/SigProfilerExtractor/issues/88
###############################################################################

###############################################################################
#%% Cell 2: Set global variables and import SigProfilerExtractor
###############################################################################

if __name__ == "__main__":

    ###########################################################################
    #  2.1: Set global variables
    ###########################################################################
    # Read old working directory
    old_wd = os.getcwd()
    # os.path.abspath extends a relative path to an absolute one.
    home_for_data = os.path.abspath("./sens_SBS35_20_728/input")
    home_for_run = os.path.abspath("./sens_SBS35_20_728/raw_results")
    # Seed numbers
    seed_numbers = (145879, 200437, 310111, 528401, 1076753)
    seed_numbers = [(sn % 10000000) for sn in seed_numbers]
    # Naming the datasets for cycling
    dataset_names = ("Realistic",)

    ###########################################################################
    #  2.2: Set NON-DEFAULT argument values
    ###########################################################################
    # 96 - signatures and spectra in SBS context
    context_type = "96"
    # Range for number of signatures to extract.
    # We set to 10..40 as there are 23 ground-truth signatures.
    min_K = 10
    max_K = 40
    # Argument values from Supp Note 1 of BioRxiv Preprint by Islam et al.
    min_nmf_iterations = 1000
    max_nmf_iterations = 200000
    nmf_test_conv = 1000
    nmf_tolerance = 1e-08
    # number of CPUs to use. Free to change.
    numcpus = 60

    ###########################################################################
    #  2.3: Import SigProfilerExtractor
    ###########################################################################
    from SigProfilerExtractor import sigpro as sig



###############################################################################
#%% Cell 3: Fetch argument values for each run,
#   and run SigProfilerExtractor using these values.
###############################################################################

if __name__ == "__main__":
    for dataset_name in dataset_names:
        for seed_number in seed_numbers:
            ###################################################################
            #  3.1: Set argument values specific for each run
            ###################################################################
            input_dir = home_for_data+"/"+dataset_name
            input_catalog = input_dir+"/ground.truth.syn.catalog.tsv"
            output_dir = "".join([home_for_run,"/SigProfilerExtractor.results/",
                                  dataset_name,"/seed.",str(seed_number),"/"])
            path=Path(output_dir)
            if path.exists()==False:
                path.mkdir(parents=True, exist_ok=True)
            # Skip the current job if it has already finished.
            path=Path(output_dir+"/profiling_info.pickle")
            if path.exists()==True:
                continue
            # Set seed
            random.seed(seed_number)
            # Generate tab separated seed file required by SigPro.
            seeds_df = pd.DataFrame([seed_number], columns=["Seed"])
            seeds_df.to_csv(os.path.abspath(output_dir+"/Seeds.txt"), sep="\t")
            seeds = os.path.abspath(output_dir+"/Seeds.txt")

            ###################################################################
            #  3.2: Enable profiling before running SigProfilerExtractor
            ###################################################################
            # Messages
            print("\n\n#####################")
            print("Start running catalog "+str(input_catalog)+"\n")
            print("Seed number:"+str(seed_number)+"\n")
            print("#####################\n")
            # Garbage collection before RAM profiling
            gc.collect()
            tracemalloc.reset_peak()
            tracemalloc.start()
            # Start timing
            #
            # os.times(): returns an object with 5 attributes:
            #
            # user: user time
            # system: system time
            # children_user: user time of all child processes
            # children_system: system time of all child processes
            # elapsed: elapsed real time since a fixed point in the past
            #
            times_start = os.times()
            
            ###################################################################
            #  3.3: Run sigProfilerExtractor -
            #  extract signatures and infer exposures
            #
            #  Argument values provided as constants are default values;
            #  whereas non-default values are listed as variables 
            #  in sections 2.2 and 3.1.
            ###################################################################
            sig.sigProfilerExtractor(
                "matrix",
                output_dir,
                input_data=input_catalog, 
                reference_genome="GRCh37",
                cosmic_version="3.2",
                context_type=context_type,
                exome=False,
                minimum_signatures=min_K,
                maximum_signatures=max_K,
                nmf_replicates=100,
                resample=True,
                seeds=seeds,
                matrix_normalization="gmm",
                nmf_init="random",
                precision="single",
                min_nmf_iterations=min_nmf_iterations,
                max_nmf_iterations=max_nmf_iterations,
                nmf_test_conv=nmf_test_conv,
                nmf_tolerance=nmf_tolerance,
                cpu=numcpus)
            
            ###################################################################
            #  3.4: Stop profiling, and calculate time and RAM usage.
            ###################################################################
            times_end=os.times()
            # Diff in os.times() in seconds
            times_diff=[times_end[ii]-times_start[ii] 
                        for ii in range(0, len(times_start))]
            # CPU time consumed by SigProfilerExtractor in seconds
            cpu_time=sum(times_diff[0:4])
            # Elapsed time consumed by SigProfilerExtractor in seconds
            elapsed_time=times_diff[4]    
            # Current and Peak memory in Bytes
            current_mem, peak_mem = tracemalloc.get_traced_memory()  
            tracemalloc.stop()
            # Storing Profiling data into a dict.
            # Store all info objects to a pickle file.
            profiling_info = {
                "cpu_time": cpu_time,
                "elapsed_time": elapsed_time,
                "times_diff": times_diff,
                "current_mem": current_mem,
                "peak_mem": peak_mem}
            fi=open(os.path.abspath(output_dir+"/profiling_info.pickle"), "wb")
            pickle.dump(obj=profiling_info, file=fi)
            fi.close()
