#!/usr/bin/python

#######################################################################
## Author: 	Thanos Mourikis
## Date: 	2014
## Purpose:	Python module for xml manipulation and BLAST interaction
#######################################################################


from Bio.Blast import NCBIXML
from Bio import SeqIO
from itertools import izip
from bisect import bisect_left
from segments import Segments
import sys
import os
import subprocess
from StringIO import StringIO
import DB
from blastdb_fetch import blastdb_fetch


def prepare(job_name, master_path, db_name):

	global job_path
	job_path = master_path + '/scan_results/' + db_name + '/' + job_name


def hsp_forward_coords(hsp): ## Return always forward coordinates
        forward_start = min(hsp.sbjct_start, hsp.sbjct_end)
        forward_end = max(hsp.sbjct_start, hsp.sbjct_end)
        return forward_start, forward_end

def blast_hit_extraction(probe_type, min_hsp_distance, e_value, ext_buffer, 
						blast_path, organism, db_name):
	
	"""Filter overlapping BLAST hits in a certain region."""

	min_hsp_distance = int(min_hsp_distance)
	e_value = float(e_value)
	ext_buffer = int(ext_buffer)

	if probe_type == 'prot':
		blast_xml = job_path + '/blast_results/tblastn_probes_to_target.xml'
	elif probe_type == 'nucl':
		blast_xml = job_path + '/blast_results/blastn_probes_to_target.xml'
	
	blast_handle = open(blast_xml, "r")
	blast_records = NCBIXML.parse(blast_handle)
	mysql_file = job_path + '/blast_results/extracted_hits.mysql'
	mysql_handle = open(mysql_file, "w")
	out_file = job_path + '/blast_results/extracted_hits.fasta'
	out_handle = open(out_file, "w")
	target_genome_blastdb = job_path + '/blast_db/Target_genome_blastdb'

	## If there are two pol hits in more than 5000kb distance
	## and the contig is small (e.g 7000) the region 0-7000 will
	## be exported twice and produce bugs in DB insertion and 
	## BLAST db construction (two gis the same)
	## So we need to keep a record of extracted seqs
	## Probably we will lose some pols but just a tiny fraction
	extracted_seqs = []
	
	sequence_dict = {} ## A dictionary per hit sequence
	

	## Iterate through xml file and extract unique hits based on parameters
	for blast_record in blast_records:
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				blast_hit_dict = {}
				blast_hit_dict['Record_ID'] =  alignment.hit_id
				blast_hit_dict['Probe_ID'] = blast_record.query				
				blast_hit_dict['Probe_type'] =  probe_type
				blast_hit_dict['Frame'] =  hsp.frame[1]
				blast_hit_dict['Bit_score'] = hsp.score
				blast_hit_dict['E_value'] = hsp.expect
				blast_hit_dict['Subject_start'] = hsp.sbjct_start
				blast_hit_dict['Subject_end'] = hsp.sbjct_end
				blast_hit_dict['Query_start'] = hsp.query_start
				blast_hit_dict['Query_end'] = hsp.query_end
				blast_hit_dict['Identities'] = hsp.identities
				blast_hit_dict['Positives'] = hsp.positives
				blast_hit_dict['Gaps'] = hsp.gaps
				blast_hit_dict['Align_len'] = hsp.align_length
				## Fill BLAST_hits table
				DB.add_to_blast_table(db_name, blast_hit_dict)
			if alignment.hit_id not in sequence_dict.keys():
				sequence_dict[alignment.hit_id] = []
				for hsp in alignment.hsps:
					if hsp.expect < e_value:
						forward_start, forward_end = hsp_forward_coords(hsp)
						sequence_dict[alignment.hit_id].append([forward_start, forward_end])
			else:
				for hsp in alignment.hsps:
					if hsp.expect < e_value:
						forward_start, forward_end = hsp_forward_coords(hsp)
						sequence_dict[alignment.hit_id].append([forward_start, forward_end])  

	correction_dict = {}
           
	for key in sequence_dict.keys():
			
		s = Segments(min_hsp_distance) ## Object to resolve overlapping regions

		for hsp in sequence_dict[key]:
			l = hsp[0]
			r = hsp[1]
			s.add(l,r)
		s.compact() ## Elimination of overlapping region
		
		
		gi = key.split("|")[1]
		record = blastdb_fetch(blast_path, target_genome_blastdb, gi)

		for (l,r) in s:
			l = l - ext_buffer
			r = r + ext_buffer

			## Check the length of the entry
			if l < 0:
				l = 0
			if r > len(record):
				r = len(record)


			subseq = record[l:r]
			subseq.id = "ref|" + gi + "_" + str(l) + "_" + str(r) + "|"
			subseq.description = " ".join(record.description.split()[1:])
			if subseq.id not in extracted_seqs:
				extracted_seqs.append(subseq.id)
				mysql_handle.write(gi + "\t" + subseq.id + "\t" + str(l) + "\t" + str(r) + "\n")
				out_handle.write(subseq.format('fasta'))
				DB.add_to_extracted_hits_table(db_name, key, subseq.id, organism, l, r)


			## Feed it back to Segements through a dictionary to correct overlapping region
			## Look below for explanation

			#if key not in correction_dict.keys(): ---- Activate if want to use below inactivated code!!!!
				#correction_dict[key] = []
			#else:
				#correction_dict[key].append([l,r])
			

		## Wipe s object
		s.wipe()



''' Not ideal solution to the problem
	## After the merging of hsps we need to re-search for overlapping regions
	## If the merging buffer is 3000 and an hsp exist at 3000 + 10 then two overlapping
	## regions will come up. We need to merge this regions as well!
	## Make a new dictionary and insert hits per contig after the min_hsp_buffer merging

	for key in correction_dict.keys():

		s = Segments(0) ## Object to resolve overlapping regions

		for extract in correction_dict[key]:

			l = extract[0]
			r = extract[1]
			s.add(l,r)
		s.compact()## Elimination of overlapping region

		for (l,r) in s:

			gi = key.split("|")[1]
			record = blastdb_fetch(blast_path, target_genome_blastdb, gi)
			subseq = record[l:r]
			subseq.id = "ref|" + gi + "_" + str(l) + "_" + str(r) + "|"
			subseq.description = " ".join(record.description.split()[1:])
			if subseq.id not in extracted_seqs:
				extracted_seqs.append(subseq.id)
				mysql_handle.write(gi + "\t" + subseq.id + "\t" + str(l) + "\t" + str(r) + "\n")
				out_handle.write(subseq.format('fasta'))
				DB.add_to_extracted_hits_table(db_name, key, subseq.id, organism, l, r)

		## Wipe s object
		s.wipe()
'''



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


def ref_similarity_extraction(reference_type, db_name):
	
	"""Filter best hits for each element against the reference database."""

	if reference_type == "prot":
		blast_xml = job_path + '/blast_results/blastx_extracted_hits_to_reference.xml'
	elif reference_type == "nucl":
		blast_xml = job_path + '/blast_results/blastn_extracted_hits_to_reference.xml'

	blast_handle = open(blast_xml, "r")
	blast_records = NCBIXML.parse(blast_handle)
	
	for blast_record in blast_records:
		blast_hit_dict = {}
		## Check if there are hits with the references
		if len(blast_record.alignments) > 0:
			blast_hit_dict['Sequence'] = blast_record.query.split()[0]
			blast_hit_dict['Reference'] = blast_record.alignments[0].hit_id
			
			## Find the HSP with the best score
			bestHspIndex = bestScoreHsp(blast_record.alignments[0].hsps)

			blast_hit_dict['Max_score'] = blast_record.alignments[0].hsps[bestHspIndex].score
			blast_hit_dict['Query_start'] = blast_record.alignments[0].hsps[bestHspIndex].query_start
			blast_hit_dict['Query_end'] = blast_record.alignments[0].hsps[bestHspIndex].query_end
			blast_hit_dict['Subject_start'] = blast_record.alignments[0].hsps[bestHspIndex].sbjct_start
			blast_hit_dict['Subject_end'] = blast_record.alignments[0].hsps[bestHspIndex].sbjct_end
			blast_hit_dict['E_value'] = blast_record.alignments[0].hsps[bestHspIndex].expect
			total_score = 0
			identities = float(blast_record.alignments[0].hsps[bestHspIndex].identities)
			align_len = float(blast_record.alignments[0].hsps[bestHspIndex].align_length)
			perc_identity = (identities/align_len) * 100
			blast_hit_dict['Perc_identity'] = perc_identity
			for hsp in blast_record.alignments[0].hsps:				
				total_score += hsp.score
			blast_hit_dict['Total_score'] = total_score
			
			## Push results to the MySQl database
			DB.add_to_extracted_hits_similarity_table(db_name, blast_hit_dict)
		else:
			blast_hit_dict['Sequence'] = blast_record.query.split()[0]
			blast_hit_dict['Reference'] = 'NA'
			blast_hit_dict['Max_score'] = 'NA'
			blast_hit_dict['Query_start'] = 'NA'
			blast_hit_dict['Query_end'] = 'NA'
			blast_hit_dict['Subject_start'] = 'NA'
			blast_hit_dict['Subject_end'] = 'NA'
			blast_hit_dict['E_value'] = 'NA'
			blast_hit_dict['Total_score'] = 'NA'
			blast_hit_dict['Perc_identity'] = 'NA'

			DB.add_to_extracted_hits_similarity_table(db_name, blast_hit_dict)










