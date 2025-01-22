# Update variables to modify paths

data_dir = "/home/katrina_liu/harvard/nucleosome_end_binding_TFs/data/"
SELEX_data_dir = paste0(data_dir, "SELEX_data/")
SELEX_data_report_tsv = paste0(data_dir, "filereport_read_run_PRJEB22684_tsv.txt") # Need to have
SELEX_tf_list_file = paste0(data_dir, "lig147_tfs_test.txt") # Need to have
SELEX_library_file = paste0(data_dir, "lig147_libs_test.txt") # Need to have
SELEX_ligand = "lig147"
SELEX_cache_dir = paste0(data_dir, "cache/")

tri_mer_bend_file = paste0(data_dir, "rod_trinuc.txt")

FW_primer = "CCCTACACGACGCTCTTCCGATCT"
RV_primer = "AGATCGGAAGAGCACACGTCTG"

FW_primer_rc = reverse(chartr("ATGC","TACG",FW_primer))
RV_primer_rc = reverse(chartr("ATGC","TACG",RV_primer))

info_gain_path <- paste0(data_dir, "info_gain.txt")
SELEX_end_binding_perc = paste0(data_dir, "end_binding_percent.csv")
SELEX_end_binder_list = paste0(data_dir, "end_binding.csv")
SELEX_non_end_binder_list = paste0(data_dir, "non_end_binding.csv")
  
  
SELEX_cyc_result_csv = paste0(data_dir, "SELEX_test_cyc_results.csv")
SELEX_internal_cyc_result_csv = paste0(data_dir, "SELEX_test_internal_cyc_results.csv")
SELEX_cyc_slope_avg = paste0(data_dir, "SELEX_cyc_slope_avg.csv")
SELEX_cyc_ht_thres_tf_csv = paste0(data_dir, "ht_thres_TF.csv")

left_region=35:74
right_region=76:115
nuc="A"
SELEX_fft_left = paste0(data_dir, "SELEX_fft_A_35_74.csv")
SELEX_fft_right = paste0(data_dir, "SELEX_fft_A_76_115.csv")


pear_path = "/home/katrina_liu/harvard/capstone/pear-0.9.11-linux-x86_64/bin/pear" ##### Update path to pear
