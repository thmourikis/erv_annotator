#!/usr/bin/python

#######################################################################
## Author: 	Thanos Mourikis
## Date: 	2014
## Purpose:	Python module to run RepeatMasker on extracted hits
#######################################################################


from blastdb_fetch import blastdb_fetch
import subprocess
from StringIO import StringIO
from Bio import SeqIO
import os
import re
import glob
import shutil
import DB

def prepare(job_name, master_path, db_name):

	global job_path
	job_path = master_path + '/scan_results/' + db_name + '/' + job_name


def run_repeatmasker(rm_path, threads, blast_path):
	print "\n\n- Running RepeatMasker on extracted hits"
	rm_directory = job_path + '/repeatmasker'
	RepeatMasker = rm_path + ' -species mammal' + ' -pa ' + threads + ' ' + job_path + '/blast_results/extracted_hits.fasta'
	os.system(RepeatMasker)
	## Keep files organised
	## RM files stored in the directory of input, let's move them
	rm_files = job_path + '/blast_results/extracted_hits.fasta.*'
	for f in glob.glob(rm_files):
		shutil.move(f, rm_directory)

	## Create BLASTdb extracted hits for the next step
	makeblastdb_extracted_hits = blast_path + '/makeblastdb -in ' + job_path + '/blast_results/extracted_hits.fasta'  + ' -dbtype nucl ' \
				+ '-title "Extracted hits BLAST database" ' + '-parse_seqids ' \
				+ '-out ' + job_path + '/blast_db/Extracted_hits_blastdb'
	os.system(makeblastdb_extracted_hits)


def mask_reps_not_LTRs(blast_path, db_name):

	'''Masks all reps but LTRs and prepares extracted hits for LTR de novo identification '''

	print "\n\n- Masking known repeats except from LTRs"

	rm_file = job_path + '/repeatmasker/extracted_hits.fasta.out'
	rm_handle = open(rm_file, "rU")
	output_file = job_path + '/repeatmasker/extracted_hits_masked.fasta'
	out_handle = open(output_file, 'w')
	log_file = job_path + '/repeatmasker/extracted_hits_masked.log'
	log_handle = open(log_file, 'w')
	blast_db = job_path + '/blast_db/Extracted_hits_blastdb'
	gis = {}

	for (line_no, line) in enumerate(rm_handle):
		if line_no > 2:
			rm = line.split()
			rm_id = rm[4].split('|')
			#if not rm[10].startswith('LTR'):
			#not re.search("ERV", rm[10]) or
			if not re.search("LTR", rm[10]): ## VERY IMPORTANT
				if rm_id[1] in gis.keys():
					#The first time a gi is met cannot have '*'
					gis[rm_id[1]].append(str(line_no)+':'+str(rm[5])+':'+str(rm[6])+':'+str(rm[10])+':'+str(rm[-1]))
					if rm[-1] is '*':
						log_handle.write(rm[5] + ' - ' + rm[6] + ' => ' + rm[10] + rm[-1] + "\n")
					else:
						log_handle.write(rm[5] + ' - ' + rm[6] + ' => ' + rm[10] + "\n")
				else:
					gis[rm_id[1]] = [str(line_no)+':'+str(rm[5])+':'+str(rm[6])+':'+str(rm[10])]
					log_handle.write("Regions removed from: " + rm[4] + "\n")
					if rm[-1] is '*':
						log_handle.write(rm[5] + ' - ' + rm[6] + ' => ' + rm[10] + rm[-1] + "\n")
					else:
						log_handle.write(rm[5] + ' - ' + rm[6] + ' => ' + rm[10] + "\n")

	for key in gis:
		record = blastdb_fetch(blast_path, blast_db, key)
		sequence = ''
		for i in range(len(gis[key])):
			rm_dict = {}
			entry = gis[key][i].split(':')
			start = entry[1]
			stop = entry[2]
			rm_dict['Sequence'] = "ref|" + key + "|"
			rm_dict['Rep_start'] = start
			rm_dict['Rep_end'] = stop
			rm_dict['Rep_type'] = entry[3]

			## Maybe re-establish the DB connection
			## At some point if the RM run is too long it gets lost
			## Insert a try statement to re-connect to the MySQL DB
			DB.add_to_rm_table(db_name, rm_dict)
			length = int(stop) - int(start)
			if i == 0:
				sequence = list(record[:int(start)-1] + (length + 1) * 'N' + record[int(stop):])
			else:
				for y in range(int(start)-1, int(stop)):
					sequence[y] = 'N'

		## Turn the fasta string from blastdbcmd into a temporary fake file handle.
		## and parse it with Biopython SeqIO, to then return a SeqRecord.
		fasta_handle = StringIO('>' + record.description + "\n" + "".join(sequence))
		temp_record = SeqIO.read(fasta_handle, 'fasta')
		out_handle.write(temp_record.format("fasta"))


	#Now retrieve the sequences that are not in the RepeatMasker output (if any)
	p = subprocess.Popen(['blastdbcmd', '-entry', 'all', '-db', blast_db,  '-outfmt', '%a'],
						stdout=subprocess.PIPE, stderr=subprocess.PIPE)

	out, err = p.communicate()
	gis_total = out.split()
	set_difference = list(set(gis_total).difference(gis.keys()))

	if len(set_difference) > 0:
		for i in set_difference:
			record = blastdb_fetch(blast_path, blast_db, i)
			out_handle.write(record.format('fasta'))
		log_handle.write(str(len(set_difference)) + """ gis were not in RepeatMasker's output""")
	else:
		log_handle.write("All the gis were in RepeatMasker's output")



	print "\n\nDone!"



