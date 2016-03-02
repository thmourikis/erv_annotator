#!/usr/bin/python

#######################################################################
## Author: 	Thanos Mourikis
## Date: 	2014
## Purpose:	Python module to annotate ERVs based on NCBI's nr hits
##			and create a gff3 file for data visualisation			
#######################################################################


from Bio.Blast import NCBIXML
import DB
from blastdb_fetch import blastdb_fetch
from Bio import SeqIO
from StringIO import StringIO
import os
from itertools import *
import itertools


def prepare(job_name, master_path, db_name):

	global job_path
	job_path = master_path + '/scan_results/' + db_name + '/' + job_name


def bestScoreHsp(l):

	maxScoreSoFar = 0
	index = 0
	for no,hsp in enumerate(l):
		if hsp.score > maxScoreSoFar:
			maxScoreSoFar = hsp.score
			index = no
		else:
			continue

	return index


def best_alignment(handle): ## Return one hsp per query
	blast_records = NCBIXML.parse(handle)
	xml_dict = {}	
	for blast_record in blast_records:
		ids = blast_record.query.split()
		query_name = str(ids[0])
		record_dict = {}
		if len(blast_record.alignments) > 0:
			## Store the best HSP (bitscore) for each query
			bestHspIndex = bestScoreHsp(blast_record.alignments[0].hsps)			
			if blast_record.alignments[0].hsps[bestHspIndex].frame[0] > 0:
				frame = '+'
			elif blast_record.alignments[0].hsps[bestHspIndex].frame[0] < 0:
				frame = '-'
			identities = float(blast_record.alignments[0].hsps[bestHspIndex].identities)
			al_len = float(blast_record.alignments[0].hsps[bestHspIndex].align_length)
			perc_identity = (identities/al_len) * 100
			record_dict[query_name] =  {'hit'				:	blast_record.alignments[0].hit_id, 
										'score'				:	blast_record.alignments[0].hsps[bestHspIndex].score,
										'perc_identity'		:	perc_identity, 
										'alignment_length'	:	blast_record.alignments[0].hsps[bestHspIndex].align_length,
										'query_start'		:	blast_record.alignments[0].hsps[bestHspIndex].query_start,
										'query_end'			:	blast_record.alignments[0].hsps[bestHspIndex].query_end,
										'sbjct_start'		:	blast_record.alignments[0].hsps[bestHspIndex].sbjct_start,
										'sbjct_end'			:	blast_record.alignments[0].hsps[bestHspIndex].sbjct_end,
										'frame'				:	frame}

		for key in record_dict.keys():
			xml_dict[key] = record_dict[key]					
	return xml_dict

def return_Ns(fasta_file, key):
	N_list = []
	fasta_handle = open(fasta_file, "r")
	records = SeqIO.parse(fasta_handle, "fasta")
	for record in records:
		fasta_ids = record.description.split()
		fasta_name = str(fasta_ids[0])
		if fasta_name == key:			
			sequence = list(record.seq)
			for (bp_no,bp) in enumerate(sequence):
				if bp=='N':
					N_list.append(bp_no)
				else:
					continue
	return N_list

def ranges(L):
	range_list = []
	G=(list(x) for _,x in groupby(enumerate(L), lambda (i,x):i-x))
	for z in ("-".join(map(str,(g[0][1],g[-1][1])[:len(g)])) for g in G):
		range_list.append(z)
	return range_list

def annotate(organism, db_name, master_path, blast_path, threads): ## Run ERVs against three databases and creates a gff file 

	print "\n\n- Annotating ERVs"

	## Sometimes LTRharvest cannot detect LTRs and returns empty gff3
	## Check if this happens and control the error instead of break the code
	## This will help in batch jobs instead of error go to the next genome
	ltrharvest_gff = job_path + '/ltrharvest/extracted_hits_masked.fasta_ltrharvest.gff3'
	empty = os.stat(ltrharvest_gff)[6]==0
	if empty==True:
		print "\n\n LTRharvest could not detect LTRs. Genome LTR finding step skipped!"
		print " You can still find extracted blast hits in job's file."
		note = open("NO_LTRs_FOUND.txt", "w")
		note.write("LTRharvest could not detect LTRs. Genome LTR finding step skipped!" + "\n")
		note.write("You can still find extracted blast hits in job's file.")
	else:

		target_db = job_path + '/blast_db/Target_genome_blastdb'

		erv_file = job_path + '/annotation/' + organism + '_ERVs.fasta'
		erv_out_handle = open(erv_file, "w")

		## Log file just in case
		log_file = job_path + '/run.log'
		log_handle = open(log_file, "w")
		log_handle.write("-- \t Annotation step log \n")
	
		## Annotation dictionary 
		annotation_dict = {}
		N_regions_count = 1


		## Extract ERV sequences in a fasta file using blastdb_fetch
		res = DB.select_ervs(db_name, ['ERV_table.Sequence', 'Element_start', 'Element_end'], ['ERV_table', 'Extracted_hits'], organism) ## Here we take back all the elements

		for x in range(len(res)):
			gi = res[x][0].split("|")[1].split("_")[0]
			element_start = res[x][1]
			element_end = res[x][2]
			record = blastdb_fetch(blast_path, target_db, gi)
			record.id = "ref|" + str(gi) + "_" + str(element_start) + "_" + str(element_end) + "|"
			record.description = record.description + ' [' + str(element_start) + '-' + str(element_end) + ']'
			erv_out_handle.write(record[element_start:element_end+1].format('fasta'))

	
		## Run against GAG DB
		blastx_gag = blast_path + '/blastx -db ' + job_path + '/blast_db/Gag_blastdb' \
				+ ' -query ' + erv_file + ' -out ' + job_path + '/annotation/blastx_ervs_to_gag.xml' \
				+ ' -outfmt 5' + ' -num_threads ' + threads + ' -max_target_seqs 10000'
		os.system(blastx_gag)

		## Run against Pol DB
		blastx_pol = blast_path + '/blastx -db ' + job_path + '/blast_db/Reference_screenset_blastdb' \
				+ ' -query ' + erv_file + ' -out ' + job_path + '/annotation/blastx_ervs_to_pol.xml' \
				+ ' -outfmt 5' + ' -num_threads ' + threads + ' -max_target_seqs 10000'
		os.system(blastx_pol)

		## Run against Pol DB
		blastx_env = blast_path + '/blastx -db ' + job_path + '/blast_db/Env_blastdb' \
				+ ' -query ' + erv_file + ' -out ' + job_path + '/annotation/blastx_ervs_to_env.xml' \
				+ ' -outfmt 5' + ' -num_threads ' + threads + ' -max_target_seqs 10000'
		os.system(blastx_env)


		gff_file = job_path + '/annotation/' + organism + '_ERVs.gff3'
		gff_handle = open(gff_file, "w")
	
		blast_gag = job_path + '/annotation/blastx_ervs_to_gag.xml'
		gag_handle = open(blast_gag, "r")
		blast_pol = job_path + '/annotation/blastx_ervs_to_pol.xml'
		pol_handle = open(blast_pol, "r")
		blast_env = job_path + '/annotation/blastx_ervs_to_env.xml'
		env_handle = open(blast_env, "r")

		## Write header section in gff
		gff_handle.write("##gff-version\t3" + "\n")

		## Write target sequences in gff3
		res2 = DB.select_contigs(db_name, ['Record_ID', 'Length'], 'Target_genome', organism)
		for x in range(len(res2)):
			gff_handle.write("##sequence-region\t" + str(res2[x][0]) \
							+ "\t1\t" + str(res2[x][1]) + "\n")

		## Analyze each xml file
		gag_record = best_alignment(gag_handle)
		for key in gag_record:
			if key in annotation_dict.keys():
				annotation_dict[key]['gag'] = gag_record[key]
			else:
				annotation_dict[key] = {'gag' : gag_record[key]}

		pol_record = best_alignment(pol_handle)
		for key in pol_record:
			if key in annotation_dict.keys():
				annotation_dict[key]['pol'] = pol_record[key]
			else:
				annotation_dict[key] = {'pol' : pol_record[key]}

		env_record = best_alignment(env_handle)
		for key in env_record:
			if key in annotation_dict.keys():
				annotation_dict[key]['env'] = env_record[key]
			else:
				annotation_dict[key] = {'env' : env_record[key]}
	
		## Count the Ns
		for key in annotation_dict.keys():		
			N_list = return_Ns(erv_file,key)
			range_list = ranges(N_list)
			for n_range in range_list:
				n_key = "N" + str(N_regions_count)
				n_coords = n_range.split("-") 
				if len(n_coords) > 1:
					annotation_dict[key][n_key] = {'start' : n_coords[0], 'end' : n_coords[1]}
					N_regions_count += 1






		## Write gff3 and populate ERV_annotation table
		for key in annotation_dict.keys():
			seq_name = DB.select_pattern(db_name, 'Target_genome', 'Record_ID', key.split("|")[1].split("_")[0]) ## Get contig name from MySQL table
			annot_dict = {} ## For db insertion
			element_start = key.split("|")[1].split("_")[1]
			element_end = key.split("|")[1].split("_")[2]

			## Write parent element
			gff_handle.write(str(seq_name[0][0]) + "\tERVannotator\tERV\t" + str(element_start) \
							+ "\t" + str(element_end) +  "\t.\t.\t.\tID=" + str(key) \
							+ ";Name=" + str(key) + ";color=#FFFFFF" + "\n")

			## Write child annotation
			## gag -> red
			## pol -> green
			## env -> magenta
			## Ns -> grey
			## LTRs -> yellow
			## rest of element ->  blue
			annot_dict = {}
			for domain in annotation_dict[key].keys():
				if domain=='gag':
					gff_handle.write(str(seq_name[0][0]) + "\tERVannotator\t" + str(domain) + "\t" + str(int(element_start) + int(annotation_dict[key][domain]['query_start'])) \
						+ "\t" + str(int(element_start) + int(annotation_dict[key][domain]['query_end'])) + "\t.\t" + str(annotation_dict[key][domain]['frame']) \
						+ "\t.\tParent=" + str(key) + "\n")
				elif domain=='pol':
					gff_handle.write(str(seq_name[0][0]) + "\tERVannotator\t" + str(domain) + "\t" + str(int(element_start) + int(annotation_dict[key][domain]['query_start'])) \
						+ "\t" + str(int(element_start) + int(annotation_dict[key][domain]['query_end'])) + "\t.\t" + str(annotation_dict[key][domain]['frame']) \
						+ "\t.\tParent=" + str(key) + ";color=#00FF00" + "\n")
				elif domain=='env':
					gff_handle.write(str(seq_name[0][0]) + "\tERVannotator\t" + str(domain) + "\t" + str(int(element_start) + int(annotation_dict[key][domain]['query_start'])) \
						+ "\t" + str(int(element_start) + int(annotation_dict[key][domain]['query_end'])) + "\t.\t" + str(annotation_dict[key][domain]['frame']) \
						+ "\t.\tParent=" + str(key) + ";color=#FF00FF" + "\n")
				elif domain.startswith("N"):
					gff_handle.write(str(seq_name[0][0]) + "\tERVannotator\tNs\t" + str(int(element_start) + int(annotation_dict[key][domain]['start'])) \
						+ "\t" + str(int(element_start) + int(annotation_dict[key][domain]['end'])) + "\t.\t.\t.\tParent=" + str(key) \
						+ ";color=#7F7F7F" + "\n")
		
			## The key is in ref format but with tblastn coordinates
			## not with corrected ones according to the LTRdigest
			## just pass the gi and it should be ok		
			res3 = DB.select_exact_erv(db_name, key.split("|")[1].split("_")[0], element_start) ## Here we take back just one line

			## Check-point for returned values
			## If more than one found (distant hsps from tblastn)
			## choose the one with biggest length and write to log
			## for future reference
			if len(res3) > 1:
				for i in range(len(res3)):
					log_handle.write(str(res3[i]) + "\n")

				res3_cleaned = list(max(res3, key=lambda x: x[3]))

				log_handle.write("\nFound multiple ERVs all on the same location! \n")
				log_handle.write("Chosen: " + str(res3_cleaned) + "\n")
				log_handle.write("----- \n\n")

				res3 = [res3_cleaned]
			elif len(res3) == 0:
				print "ERROR: zero values returned from ERV_table"
				exit(1)

			gff_handle.write(str(seq_name[0][0]) + "\tERVannotator\tlLTR\t" + str(res3[0][4]) \
						+ "\t" + str(res3[0][5]) + "\t.\t.\t.\tParent=" + str(key) \
						+ ";color=#FFFF00"+ "\n")
			gff_handle.write(str(seq_name[0][0]) + "\tERVannotator\trLTR\t" + str(res3[0][7]) \
						+ "\t" + str(res3[0][8]) + "\t.\t.\t.\tParent=" + str(key) \
						+ ";color=#FFFF00"+ "\n")


			## MySQL DB insertion
			#seq_name = DB.select_pattern(db_name, 'Target_genome', 'Record_ID', key.split("|")[1].split("_")[0]) ## Get contig name from MySQL table
			annot_dict['Record_ID'] = seq_name[0][0]
			annot_dict['Element_start'] = res3[0][1] 
			annot_dict['Element_end'] = res3[0][2]
			annot_dict['lLTR_start'] = res3[0][4] 
			annot_dict['lLTR_end'] = res3[0][5] 
			annot_dict['rLTR_start'] = res3[0][7] 
			annot_dict['rLTR_end'] = res3[0][8]
			if 'gag' in annotation_dict[key].keys(): ## Probably some hits do not have gag, pol or env hits
				annot_dict['Gag_start'] = res3[0][1] + int(annotation_dict[key]['gag']['query_start'])
				annot_dict['Gag_end'] = res3[0][1] + int(annotation_dict[key]['gag']['query_end'])
				annot_dict['Gag_frame'] = annotation_dict[key]['gag']['frame']
				annot_dict['Gag_hit'] = annotation_dict[key]['gag']['hit']
				annot_dict['Gag_perc_identity'] = annotation_dict[key]['gag']['perc_identity']
			else:
				annot_dict['Gag_start'] = 'NA'
				annot_dict['Gag_end'] = 'NA'
				annot_dict['Gag_frame'] = 'NA'
				annot_dict['Gag_hit'] = 'NA'
				annot_dict['Gag_perc_identity'] = 'NA'

			if 'pol' in annotation_dict[key].keys():
				annot_dict['Pol_start'] = res3[0][1] + int(annotation_dict[key]['pol']['query_start'])
				annot_dict['Pol_end'] = res3[0][1] + int(annotation_dict[key]['pol']['query_end'])
				annot_dict['Pol_frame'] = annotation_dict[key]['pol']['frame']
				annot_dict['Pol_hit'] = annotation_dict[key]['pol']['hit']
				annot_dict['Pol_perc_identity'] = annotation_dict[key]['pol']['perc_identity']
			else:
				annot_dict['Pol_start'] = 'NA'
				annot_dict['Pol_end'] = 'NA'
				annot_dict['Pol_frame'] = 'NA'
				annot_dict['Pol_hit'] = 'NA'
				annot_dict['Pol_perc_identity'] = 'NA'

			if 'env' in annotation_dict[key].keys():
				annot_dict['Env_start'] = res3[0][1] + int(annotation_dict[key]['env']['query_start'])
				annot_dict['Env_end'] = res3[0][1] + int(annotation_dict[key]['env']['query_end'])
				annot_dict['Env_frame'] = annotation_dict[key]['env']['frame']
				annot_dict['Env_hit'] = annotation_dict[key]['env']['hit']
				annot_dict['Env_perc_identity'] = annotation_dict[key]['env']['perc_identity']
			else:
				annot_dict['Env_start'] = 'NA'
				annot_dict['Env_end'] = 'NA'
				annot_dict['Env_frame'] = 'NA'
				annot_dict['Env_hit'] = 'NA'
				annot_dict['Env_perc_identity'] = 'NA'


			DB.add_to_erv_annotation(db_name, annot_dict)




		print "\n\n Done!"

	## Reset path to the main executable so master path of the next step parsed correctly
	## Ugly. Probably in  a future version fix it along all submodules.
	os.chdir(master_path + '/ervannotator') 


	

	
		
			



