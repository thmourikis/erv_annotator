#!/usr/bin/python

#######################################################################
## Author: 	Thanos Mourikis
## Date: 	2014
## Purpose:	Python module for genome manipulation
#######################################################################

from Bio import SeqIO
from StringIO import StringIO
import os
import sys
import DB


def prepare(job_name, master_path, db_name):

	global job_path
	job_path = master_path + '/scan_results/' + db_name + '/' + job_name

## Read fasta file and fill in the target table
def insert_to_basic_tables(genome_path, organism, probe_file_path,
						reference_file_path, db_name, job_name, probe_type,
						reference_type, min_hsp_distance, e_value, ext_buffer, 
						ltrharvest_sim, mintsd):

	## Target_table
	genome_handle = open(genome_path, "r")
	genome_records = SeqIO.parse(genome_handle, "fasta")

	for genome_record in genome_records:
		data_targets = {}
		data_targets['Record_ID'] = genome_record.id
		data_targets['Organism'] = organism
		data_targets['Record_desc'] = " ".join(genome_record.description.split()[1:])
		data_targets['Length'] = len(genome_record)
		DB.add_to_target_table(db_name, data_targets)

	## Probes_table
	probe_handle = open(probe_file_path, "r")
	probe_records = SeqIO.parse(probe_handle, "fasta")

	for probe_record in probe_records:
		probe_dict = {}
		probe_dict['Probe_name'] = probe_record.id
		probe_dict['Length'] = len(probe_record)
		DB.add_to_probes_table(db_name, probe_dict)

	## References_table
	reference_handle = open(reference_file_path, "r")
	reference_records = SeqIO.parse(reference_handle, "fasta")

	for reference_record in reference_records:
		reference_dict = {}
		reference_dict['Reference_name'] = reference_record.id
		reference_dict['Length'] = len(reference_record)
		DB.add_to_references_table(db_name, reference_dict)

	## Parameters_table
	parameters_dict = {}
	parameters_dict['DB_name'] = db_name
	parameters_dict['Job_name'] = job_name
	parameters_dict['Probe_file'] = probe_file_path
	parameters_dict['Reference_file'] = reference_file_path
	parameters_dict['Genome'] = genome_path
	parameters_dict['Probe_type'] = probe_type
	parameters_dict['Reference_type'] = reference_type
	parameters_dict['Min_hsp_distance'] = min_hsp_distance
	parameters_dict['E_value'] = e_value
	parameters_dict['Ext_buffer'] = ext_buffer
	parameters_dict['LTRharvest_sim'] = ltrharvest_sim
	parameters_dict['MinTSD'] = mintsd
	DB.add_to_parameters_table(db_name, parameters_dict) 

## Create a BLAST database of the target genome
def create_BLAST_db(blast_path, genome_path, reference_type, 
					reference_file_path, gag_file, env_file):

	## Make target genome and reference database
	makeblastdb_targets = blast_path + '/makeblastdb -in ' + genome_path + ' -dbtype nucl ' \
				+ '-title "Target genome BLAST database" ' + '-parse_seqids ' \
				+ '-out ' + job_path + '/blast_db/Target_genome_blastdb'
	os.system(makeblastdb_targets)


	if reference_type == "prot":
		makeblastdb_reference = blast_path + '/makeblastdb -in ' + reference_file_path + ' -dbtype prot ' \
					+ '-title "Reference screenset BLAST database" ' + '-parse_seqids ' \
					+ '-out ' + job_path + '/blast_db/Reference_screenset_blastdb'
		os.system(makeblastdb_reference)
	elif reference_type == "nucl":
		makeblastdb_reference = blast_path + '/makeblastdb -in ' + reference_file_path + ' -dbtype nucl ' \
					+ '-title "Reference screenset BLAST database" ' + '-parse_seqids ' \
					+ '-out ' + job_path + '/blast_db/Reference_screenset_blastdb'
		os.system(makeblastdb_reference)


	## Make annotation databases
	makeblastdb_gag = blast_path + '/makeblastdb -in ' + gag_file + ' -dbtype prot ' \
				+ '-title "Gag annotation BLAST database" ' + '-parse_seqids ' \
				+ '-out ' + job_path + '/blast_db/Gag_blastdb'
	os.system(makeblastdb_gag)

	makeblastdb_env = blast_path + '/makeblastdb -in ' + env_file + ' -dbtype prot ' \
				+ '-title "Env annotation BLAST database" ' + '-parse_seqids ' \
				+ '-out ' + job_path + '/blast_db/Env_blastdb'
	os.system(makeblastdb_env)

def run_first_blast(probe_type, blast_path, probe_file_path, threads,
					megablast):
	print "\n\n******************** First round of BLAST ********************"
	if probe_type == 'prot':
		print "- Running tblastn using protein probes against the genome database"
		tblastn = blast_path + '/tblastn -db ' + job_path + '/blast_db/Target_genome_blastdb' \
				+ ' -query ' + probe_file_path + ' -out ' + job_path + '/blast_results/tblastn_probes_to_target.xml' \
				+ ' -outfmt 5' + ' -num_threads ' + threads + ' -max_target_seqs 10000'
		os.system(tblastn)
		print "\n\n Done!"
	elif probe_type == 'nucl':
		if megablast == 'no':
			print "- Running blastn using nucleotide probes against the genome database"
			blastn = blast_path + '/blastn -db ' + job_path + '/blast_db/Target_genome_blastdb' \
					+ ' -query ' + probe_file_path + ' -out ' + job_path + '/blast_results/blastn_probes_to_target.xml' \
					+ ' -outfmt 5' + ' -num_threads ' + threads + ' -max_target_seqs 10000'
			os.system(blastn)
			print "\n\n Done!"
		elif megablast == 'yes':
			print "- Running blastn megablast using nucleotide probes against the genome database"
			blastn = blast_path + '/blastn -db ' + job_path + '/blast_db/Target_genome_blastdb' \
					+ ' -query ' + probe_file_path + ' -task megablast -out ' + job_path + '/blast_results/blastn_probes_to_target.xml' \
					+ ' -outfmt 5' + ' -num_threads ' + threads + ' -max_target_seqs 10000'
			os.system(blastn)
			print "\n\n Done!"

def run_second_blast(reference_type, blast_path, threads, megablast):
	print "\n\n******************** Second round of BLAST ********************"
	if reference_type == 'prot':
		print "- Running blastx using nucleotide extracted hits against the protein reference database"
		blastx = blast_path + '/blastx -db ' + job_path + '/blast_db/Reference_screenset_blastdb' \
				+ ' -query ' + job_path + '/blast_results/extracted_hits.fasta' + ' -out ' + job_path + '/blast_results/blastx_extracted_hits_to_reference.xml' \
				+ ' -outfmt 5' + ' -num_threads ' + threads + ' -max_target_seqs 10000'
		os.system(blastx)
	elif reference_type == 'nucl':
		if megablast == 'no':
			print "- Running blastn using nucleotide extracted hits against the nucleotide reference database"
			blastn = blast_path + '/blastn -db ' + job_path + '/blast_db/Reference_screenset_blastdb' \
					+ ' -query ' + job_path + '/blast_results/extracted_hits.fasta' + ' -out ' + job_path + '/blast_results/blastn_extracted_hits_to_reference.xml' \
					+ ' -outfmt 5' + ' -num_threads ' + threads + ' -max_target_seqs 10000'
			os.system(blastn)
		elif megablast == 'yes':
			print "- Running blastn megablast using nucleotide extracted hits against the nucleotide reference database"
			blastn = blast_path + '/blastn -db ' + job_path + '/blast_db/Reference_screenset_blastdb' \
					+ ' -query ' + job_path + '/blast_results/extracted_hits.fasta' + ' -task megablast -out ' + job_path + '/blast_results/blastn_extracted_hits_to_reference.xml' \
					+ ' -outfmt 5' + ' -num_threads ' + threads + ' -max_target_seqs 10000'
			os.system(blastn)









