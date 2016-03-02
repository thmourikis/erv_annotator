#!/usr/bin/python

#######################################################################
## Author: 	Thanos Mourikis
## Date: 	2014
## Purpose:	Python module to run LTRharvest on masked extracted hits
#######################################################################


import sys
import subprocess
import os
import glob
import shutil
import DB

def prepare(job_name, master_path, db_name):

	global job_path
	job_path = master_path + '/scan_results/' + db_name + '/' + job_name


	global fl
	fl = job_path + '/repeatmasker/extracted_hits_masked.fasta'
	global ltrharvest_directory
	ltrharvest_directory = job_path + '/ltrharvest'


def run_ltrharvest(gt_path, ltrharvest_sim, mintsd, threads, trnas_path, hmms_path):

	''' Run LTRharvest and LTRdigest'''

	print "\n\n- Running LTRharvest on masked extracted hits"

	## Local copy distributed with ERVannotator
	gt_exe = gt_path + '/gt'

	## Run Suffixerator (see genometools manual)
	suffix_args = [gt_exe, 'suffixerator', '-db', fl, '-indexname', fl, '-tis', '-suf', '-lcp', '-des', '-ssp', '-sds', '-dna']	
	suffix = subprocess.Popen(suffix_args,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	Popen_stdout,Popen_stderr=suffix.communicate()
	if len(Popen_stderr)!=0 or suffix.returncode!=0:
		sys.stderr.write(Popen_stderr)
		exit('Suffixerator error')
	
	## Run LTRharvest (see genometools manual)
	ltrharvest_args = [gt_exe, 'ltrharvest', '-index', fl, '-similar', ltrharvest_sim, '-mintsd', mintsd, '-v', '-out', fl+'_ltrharvest', '-outinner', fl+'_ltrharvest_inner', '-gff3', fl+'_ltrharvest.gff3']
	ltrharvest = subprocess.Popen(ltrharvest_args,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	Popen_stdout,Popen_stderr=ltrharvest.communicate()
	if len(Popen_stderr)!=0 or ltrharvest.returncode!=0:
		sys.stderr.write(Popen_stderr)
		exit('LTRharvest error')
	with open(fl+'_ltrharvest.log', 'w') as ltrharvest_log:
		for line in Popen_stdout:
			ltrharvest_log.write(line)
	ltrharvest.wait()

	## Run GFF sorting (see genometools manual)
	gff_sort_args = [gt_exe, 'gff3', '-sort', fl+'_ltrharvest.gff3']	
	gff_sort = subprocess.Popen(gff_sort_args,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	Popen_stdout,Popen_stderr=gff_sort.communicate()	
	if len(Popen_stderr)!=0 or gff_sort.returncode!=0:
		sys.stderr.write(Popen_stderr)
	with open(fl+'_sorted_ltrharvest.gff3', 'w') as sorted_gff:
		for line in Popen_stdout:
			sorted_gff.write(line)
	gff_sort.wait()

	print "\n\nDone!"

	print "\n\n- Running LTRdigest on masked extracted hits"

	## Run LTRdigest (see genometools manual)
	## Can't run with subprocess because wildcard in hmm specification is not recognised
	ltrdigest = gt_exe + ' -j ' + threads + ' ltrdigest' + ' -trnas ' + trnas_path + ' -hmms ' + hmms_path \
			+ ' -outfileprefix ' + 'ltrdigest ' + fl+'_sorted_ltrharvest.gff3 ' + fl + ' > ltrs_after_ltrdigest.gff3'

	os.chdir(ltrharvest_directory)
	os.system(ltrdigest)

	## Clean up the directories
	ltrharvest_files = job_path + '/repeatmasker/extracted_hits_masked.fasta.*'
	for f in glob.glob(ltrharvest_files):
		shutil.move(f, ltrharvest_directory)

	ltrharvest_files = job_path + '/repeatmasker/extracted_hits_masked.fasta_*'
	for f in glob.glob(ltrharvest_files):
		shutil.move(f, ltrharvest_directory)


	print "\n\nDone!"


def domain_check_per_gi(list):
    count = 0
    for i in range(len(list)):
        values = list[i].split("\t")
        if values[-1] and not values[-1].isspace():
            count += 1
    return count


def ltrharvest_parsing(db_name):

	''' Parse LTRdigest output for filtering. Add to mysql table the appropriate values'''

	print "\n\n- Filtering LTRdigest output"

	ltrdigest_csv = ltrharvest_directory + '/ltrdigest_tabout.csv'
	ltrdigest_handle = open(ltrdigest_csv, "rU")
	out_file = ltrharvest_directory + '/ltrdigest_tabout_filtered.csv'
	out_handle = open(out_file, "w") 

	filtering_dict = {}
	erv_dict = {}


	## Filter elements with domain hits
	for (line_no, line) in enumerate(ltrdigest_handle):
		if line_no == 0:
			out_handle.write(line)
		elif line_no > 0 and not line.startswith('\s'):
			value_list = line.split("\t")
			if value_list[-1] and not value_list[-1].isspace():
				out_handle.write(line)

	out_handle.close()

	## Reopen it and pass values to mysql database
	## A bit of ugly code but probably clear
	filtered_handle = open(out_file, 'rU') ## Our new input


	for line_no,line in enumerate(filtered_handle):
		if line_no > 0:
			values_dict = {}
			values = line.split("\t")
			values_dict['Element_start'] = int(values[0])
			values_dict['Element_end'] = int(values[1])
			## IMPORTANT because LTRdigest takes 20 chars of sequence name
			## and fill gaps with '_' which interferes with our naming scheme
			if values[3].count('_') > 2:
				sequence_id = '_'.join(values[3].split('_')[:3])
			else:
				sequence_id = values[3]
			values_dict['lLTR_start'] = int(values[4])
			values_dict['lLTR_end'] = int(values[5])
			
			values_dict['rLTR_start'] = int(values[7])
			values_dict['rLTR_end'] = int(values[8])
			
			values_dict['lTSD_start'] = int(values[10])
			values_dict['lTSD_end'] = int(values[11])			
			values_dict['rTSD_start'] = int(values[13])
			values_dict['rTSD_end'] = int(values[14])
			

			## Correct values, add appropriate keys to dict and pass it to mysql
			## Return a list with Extracted hits entries
			res = DB.select_pattern(db_name, "Extracted_hits", "Sequence", sequence_id)
        	
        	## Check point for returned values
        	## Here, unlike Annotation module,
        	## always MUST return just one value
			if len(res) > 1:
				print "ERVannotator: ERROR detected!"
				print "Returned more than one entries from Extracted_hits table!"
				exit()

        	## Add tblastn start to all values to correct them
			for key in values_dict.keys():				
				values_dict[key] += res[0][3]

			try:
				values_dict['PPT_start'] = int(values[16]) + res[0][3]				
			except ValueError:
				values_dict['PPT_start'] = 'NA'

			try:
				values_dict['PPT_end'] = int(values[17]) + res[0][3]				
			except ValueError:
				values_dict['PPT_end'] = 'NA'

			try:
				values_dict['PBS_start'] = int(values[21]) + res[0][3]				
			except ValueError:
				values_dict['PBS_start'] = 'NA'

			try:
				values_dict['PBS_end'] = int(values[22]) + res[0][3]				
			except ValueError:
				values_dict['PBS_end'] = 'NA'

			## Add missing elements
			values_dict['Element_length'] = int(values[2])
			values_dict['lLTR_length'] = int(values[6])
			values_dict['rLTR_length'] = int(values[9])
			values_dict['Sequence'] = res[0][1]
			values_dict['lTSD_motif'] = values[12]
			values_dict['rTSD_motif'] = values[15]
			values_dict['PPT_motif'] = values[18]
			values_dict['PPT_strand'] = values[19]
			values_dict['PBS_strand'] = values[23]
			values_dict['Domains'] = values[29]
			
			DB.add_to_erv_table(db_name, values_dict)

	print "\n\n Done!"












