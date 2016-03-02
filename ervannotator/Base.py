#!/usr/bin/python

#######################################################################
## Author: 	Thanos Mourikis
## Date: 	2014
## Purpose:	Python module for basic functions
#######################################################################

import os

## Create job folder -- Possible extension to let flush additional jobs to an existing one

def prepare(job_name, master_path, db_name):

	job_path = master_path + '/scan_results/' + db_name + '/' + job_name

	if not os.path.exists(job_path):
		os.makedirs(job_path)
		#log_fn = job_path + '/run.log' ## This has to be replaced with proper logger object
		#open(log_fn, 'a').close()
		blast_db = job_path + '/blast_db'
		os.makedirs(blast_db)
		blast_results = job_path + '/blast_results'
		os .makedirs(blast_results)
		repeatmasker = job_path + '/repeatmasker'
		os .makedirs(repeatmasker)
		ltrharvest = job_path + '/ltrharvest'
		os .makedirs(ltrharvest)
		annotation = job_path + '/annotation'
		os.makedirs(annotation)
	else:
		print "\n\n ERROR: Job folder exists. ERVannotator cannot overwrite result files. Please rename previous results or change genome file name."
		exit(1)


