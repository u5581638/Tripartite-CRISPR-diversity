# program to facilitate spacer mapping parallelisation

import pilercr_pos_extractor_annotation_full_spacers
import crispr_detect_tabulation
import crisprdetect_inversion
import crispr_mapping_parser_draft_v7_functions_standalone_upgraded
import csv
import subprocess
import crt_reconciliation7
import pilercr_reconciliation5
import fcntl
import os
import shutil

def rejoiner (in_str):
	i = 0
	ret_str = "" 
	while (i < len(in_str) - 1):
		ret_str += in_str[i] + "/"
		i += 1
	return ret_str	


def parallel_spacer_partition (output_dir, partition_dir, db_directory_path, b_basename, crispr_detect=0, crispr_orientation=0, partition_cores=1):
	dir_dir = rejoiner(output_dir.split("/"))
	sir_dir = partition_dir + '_dir/'
	tmp_dir = partition_dir + '_tmp/'
	subprocess.run(["/g/data/va71/crispr_pipeline_annotation/pilercr1.06/pilercr", "-in", partition_dir, "-out", partition_dir + "_crisprs" + ".lst", "-quiet"]) # run all spacers
	
	pilercr_real_table = pilercr_pos_extractor_annotation_full_spacers.all_info_extractor(partition_dir + "_crisprs" + ".lst") # need to add a line for the final repeat!!
	if (crispr_detect == 1):
		os.mkdir(tmp_dir)
		subprocess.run(["/g/data/va71/crispr_pipeline_annotation/CrisprDetect/CRISPRDetect_2.2-master/crispr_detection_running_annotation.sh",partition_dir,partition_dir + "_crispr_detect_results.txt", str(partition_cores), tmp_dir])
		subprocess.run(["rm -r " + tmp_dir ], shell=True)
		crispr_detect_tabulation.tabulate(partition_dir + "_crispr_detect_results.txt") # need to add a line for the final repeat!!!
		crisprdetect_inversion.invert(partition_dir + "_crispr_detect_results.txt" + "_crisprdetect_results.csv") # need to add a line of the final repeat!!
	if (crispr_orientation == 1):
		os.mkdir(sir_dir)
		shutil.copyfile(partition_dir, sir_dir + b_basename)
		subprocess.run(["python3 " + "/g/data/va71/crispr_pipeline_annotation/annotation_database/annotation_upgrades/CRISPRleader1.0.3/CRISPRleader.py " + "f_c_o " + sir_dir + b_basename + " partial" + " b"],shell=True)
		with open(sir_dir + "Output/CRISPRs_Repeats_Spacers.out") as csvfile:
			spacers_table = list(csv.reader(csvfile))
		with open(sir_dir + "Output/CRISPR_Info.txt") as csvfile2:
			info_table = list(csv.reader(csvfile2))	
		crispr_mapping_parser_draft_v7_functions_standalone_upgraded.spacer_crt_table_generation( spacers_table,info_table, partition_dir ) # need to add a line for the final repeat.
	#	crispr_mapping_parser_draft_v7_functions_standalone_upgraded.spacer_crt_table_generation( output_dir + "Output/CRISPRs_Repeats_Spacers.out",output_dir + "Output/CRISPR_info.txt", partition_dir )
		subprocess.run(["rm -r " + sir_dir],shell=True)
	if (crispr_detect == 1 and crispr_orientation == 1):
		detect_dir = partition_dir + "_crispr_detect_results.txt" + "_crisprdetect_results.csv" + "_inverted.csv"
		crt_dir = partition_dir + "_orientation_out.csv"
		pilercr_dir = partition_dir + "_crisprs" + ".lst" + "_full_real_arr_positions.csv"
		crt_reconciliation7.crt_reconcile(detect_dir,crt_dir, partition_dir + "_fused_crdetect_crt.csv") # need to modify to account for reconciling the final repeat!!
		pilercr_reconciliation5.pilercr_reconcile(partition_dir + "_fused_crdetect_crt.csv", pilercr_dir, partition_dir + "_reconciled_all.csv") # need to modify to account for reconciling the final repeat!!!


	if (crispr_detect == 1 and crispr_orientation == 1):
	# add code regarding catting reconciliated table
	# need to change this to use python append function and fcntl instead!!
			
		output_path_handle = db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_reconciled_full_arr_positions.csv"
		append_path_handle = partition_dir + "_reconciled_all.csv"
	#	subprocess.run(["cat " + partition_dir + "_reconciled_all.csv" + " >> " + db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_reconciled_full_arr_positions.csv"], shell=True)
	else:	
		output_path = db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_full_real_arr_positions.csv"
		append_path_handle = partition_dir + "_crisprs" + ".lst" + "_full_real_arr_positions.csv"
		#subprocess.run(["cat " + partition_dir + "_crisprs" + ".lst" + "_full_real_arr_positions.csv" + " >> " + db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_full_real_arr_positions.csv"], shell=True)
	
	out_path = open(output_path_handle, "a")
	append_path = open(append_path_handle, "r")
	fcntl.flock(out_path.fileno(), fcntl.LOCK_EX)
	my_data = append_path.read()
	out_path.write(my_data) # may need to add an end of line seperator
	append_path.close()
	out_path.close()
	# rm statements	
#	if (crispr_detect == 1):
#		subprocess.run(["rm -r " + partition_dir + "_crispr_detect_results.txt" ], shell=True)
#		subprocess.run(["rm -r " + partition_dir + "_crispr_detect_results.txt" + "_crisprdetect_results.csv"], shell=True)
#	if (crispr_orientation == 1 ):
#		subprocess.run(["rm -r " + partition_dir + "_orientation_out.csv" ], shell=True)
#	if (crispr_detect == 1 and crispr_orientation == 1):
#		subprocess.run(["rm -r " + partition_dir + "_crisprs" + ".lst" + "_full_real_arr_positions.csv"], shell=True)
#		subprocess.run(["rm -r " + partition_dir + "_fused_crdetect_crt.csv"], shell=True)
#		subprocess.run(["rm -r " + partition_dir + "_crispr_detect_results.txt" + "_crisprdetect_results.csv" + "_inverted.csv"], shell=True)
	return 0
'''	if (tony_mapping_switch == 1):
		pilercr_dr_real_table = pilercr_pos_extractor_annotation_full_spacers.all_info_extractor_tony_DR_spacers(partition_dir + "_crisprs" + ".lst")
		if (crispr_detect == 1):
			subprocess.run(["/g/data/va71/crispr_pipeline_annotation/CrisprDetect/CRISPRDetect_2.2-master/crispr_detection_running_annotation.sh",partition_dir,partition_dir + "_dr_crispr_detect_results.txt", str(cores)])
			crispr_detect_tabulation.tabulate(partition_dir + "_dr_crispr_detect_results.txt")
			crisprdetect_inversion.invert(partition_dir + "_dr_crispr_detect_results.txt" + "_crisprdetect_results.csv")
		if (crispr_orientation == 1):
			subprocess.run(["/g/data/va71/crispr_pipeline_annotation/annotation_database/annotation_upgrades/CRISPRleader1.0.3/CRISPRleader.py", "f_c_o", partition_dir, "partial", "b"])
			crt_orientation_table = crispr_mapping_parser_draft_v7_functions_working_standalone_upgraded.spacer_crt_table_generation(output_dir + "Output/CRISPR_info.txt", output_dir + "Output/CRISPRs_Repeats_Spacers.out")
			subprocess.run(["rm -r " + output_dir + "Output/"],shell=True)
		if (crispr_detect == 1 and crispr_orientation == 1):
			detect_dir = partition_dir + "_dr_crispr_detect_results.txt" + "_crisprdetect_results.csv" + "inverted.csv"
			crt_dir = partition_dir + "_orientation_out.csv"
			pilercr_dir = partition_dir + "_crisprs" + ".lst" + "_full_real_arr_positions.csv"
			crt_reconciliation.crt_reconcile(detect_dir,crt_dir, partition_dir + "_fused_crdetect_crt.csv")
			pilercr_reconciliation.pilercr_reconcile(partition_dir + "_fused_crdetect_crt.csv", pilercr_dir, partition_dir + "_dr_reconciled_all.csv")	
		#	crt_reconciliation.reconcile(pilercr_real_table, crt_orientation_table) # should add CRISPRdetect/direction predictions as well!!

	if (crispr_detect == 1 and crispr_orientation == 1):

		output_path_handle = db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_drs_reconciled_full_arr_positions.csv"
		append_path_handle = partition_dir + "_reconciled_all.csv"
		#	subprocess.run(["cat " + partition_dir + "_reconciled_all.csv" + " >> " + db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_drs_reconciled_full_arr_positions.csv"], shell=True)
	else:
		output_path_handle = db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_full_real_arr_positions_tony.csv"
		append_path_handle = partition_dir + "_crisprs" + ".lst" + "_full_real_arr_positions_tony.csv"
	#	subprocess.run(["cat " + partition_dir + "_crisprs" + ".lst" + "_full_real_arr_positions_tony.csv" + " >> " + db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_full_real_arr_positions_tony.csv"], shell=True)
	out_path = open(output_path_handle, "a")
	append_path = open(append_path_handle, "r")
	fcntl.flock(out_path.fileno(), fcntl.LOCK_EX)
	my_data = append_path.read()
	output_path.write(my_data) # may need to add an end of line seperator
	append_path.close()
	out_path.close()
	

	if (crispr_detect == 1):
		subprocess.run(["rm -r " + partition_dir + "_dr_crispr_detect_results.txt" ], shell=True)
		subprocess.run([partition_dir + "_dr_crispr_detect_results.txt" + "_crisprdetect_results.csv"], shell=True)
	if (crispr_orientation == 1 ):
		subprocess.run(["rm -r " + partition_dir + "_orientation_out.csv" ], shell=True)
	if (crispr_detect == 1 and crispr_orientation == 1):
		subprocess.run(["rm -r " + partition_dir + "_crisprs" + ".lst" + "_full_real_arr_positions.csv"], shell=True)
		subprocess.run(["rm -r " + partition_dir + "_fused_crdetect_crt.csv"], shell=True)
		subprocess.run(["rm -r " + partition_dir + "_dr_reconciled_all.csv" ], shell=True)
		subprocess.run(["rm -r " + partition_dir + "_crispr_detect_results.txt" + "_crisprdetect_results.csv" + "inverted.csv"], shell=True)
	'''	
# need to change the return value to downstream input variables for annotation!!
