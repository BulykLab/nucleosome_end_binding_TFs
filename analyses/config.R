# Update variables to modify paths

SELEX_data_dir = "../data/SELEX_data/"
SELEX_data_report_tsv = "../data/filereport_read_run_PRJEB22684_tsv.txt"
SELEX_tf_list_file = "../data/lig147_tfs_test.txt"
SELEX_library_file = "../data/lig147_libs_test.txt"
SELEX_ligand = "lig147"
SELEX_cache_dir = "../data/cache/"

FW_primer = "CCCTACACGACGCTCTTCCGATCT"
RV_primer = "AGATCGGAAGAGCACACGTCTG"

FW_primer_rc = reverse(chartr("ATGC","TACG",FW_primer))
RV_primer_rc = reverse(chartr("ATGC","TACG",RV_primer))

info_gain_path <- "../data/info_gain.txt"
SELEX_end_binding_perc = "../data/end_binding_percent.csv"
SELEX_end_binder_list = "../data/end_binding.csv"
SELEX_non_end_binder_list = "../data/non_end_binding.csv"
  
  
SELEX_cyc_result_csv = "../data/SELEX_test_cyc_results.csv"
SELEX_internal_cyc_result_csv = "../data/SELEX_test_internal_cyc_results.csv"
SELEX_cyc_slope_avg = "../data/SELEX_cyc_slope_avg.csv"
SELEX_cyc_ht_thres_tf_csv = "../data/ht_thres_TF.csv"

pear_path = "/home/xil493/capstone/pear-0.9.11-linux-x86_64/bin/pear" ##### Update path to pear
