#!/usr/bin/python
## (c) N.C. Kist (nico@kist.nl)
#from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from collections import defaultdict
import sys
def parse_exonerate(exonerate_file):
	hits_list=[]
	hit_dict={}
	for (counter,line) in enumerate(exonerate_file):
#		print line
		line=line.rstrip() ## remove whitespace
		if counter<=1:
			if line[:7] == 'Command' or line[:8]=='Hostname':
				continue
			else:
				print counter,line
				exit('ERROR: Not an exonerate file')
				
		if line[:12]=='-- completed': ## add the final alignment to the list.
			hit_dict=annotate_hitdict(hit_dict)
			hits_list.append(hit_dict)
			continue
		if(line[:7] == 'sugar: '):
			if( hit_dict != {} ): ## Don't add an empty dict on the first run through.
				hit_dict=annotate_hitdict(hit_dict)
				hits_list.append(hit_dict)
			hit_dict={}
			hit_dict['sugar_string']=line
			sugar_list=line.split(' ')
			hit_dict['query_id']=sugar_list[1]
			hit_dict['query_start']=int(sugar_list[2])
			hit_dict['query_end']=int(sugar_list[3])
			hit_dict['query_strand']=sugar_list[4]
			hit_dict['target_id']=sugar_list[5]
			hit_dict['target_start']=int(sugar_list[6])
			hit_dict['target_end']=int(sugar_list[7])
			hit_dict['target_strand']=sugar_list[8]
			hit_dict['score']=int(sugar_list[9])
			hit_dict['target_seq']=''
			hit_dict['query_length']=abs(hit_dict['query_end']-hit_dict['query_start'])
		else:
			hit_dict['target_seq']+=line
	return(hits_list)

def annotate_hitdict(hit_dict):
	hit_dict['target_length']=len(hit_dict['target_seq'])
	if(hit_dict['target_strand']=='-'):
		hit_dict['target_forward_start']=hit_dict['target_end']
		hit_dict['target_forward_end']=hit_dict['target_start']
	else:
		hit_dict['target_forward_start']=hit_dict['target_start']
		hit_dict['target_forward_end']=hit_dict['target_end']
	return(hit_dict)
def output_hit(hit):
	record=SeqRecord(Seq(hit['target_seq']))
	record.id="%s||%d||%d||%s" % (hit['target_id'],hit['target_start'], \
		hit['target_end'],hit['query_id'])
	record.description=hit['sugar_string']
	print(record.format('fasta'))


	#print(">%s::%d::%d::%s %s" % (hit['target_id'],hit['target_start'], \
		#hit['target_end'],hit['query_id'],hit['sugar_string']))
	#print(hit['target_seq'])
	

def best_alignment_per_target(exonerate_hitlist):
	''' For each target sequence used in the alignments, find the matching query sequence 
with the highest score. Returns a dictionary with the target ids as key and the
hit_dicts as value. Input is the output of parse_exonerate. '''
	hitsDict={}
	for hit in exonerate_hitlist:
		## Initialise a list for every target (marmoset) sequence id
		## This could be done much more efficiently
		hitsDict[hit['target_id']]=[]
		
	for hit in exonerate_hitlist:
		hitsDict[hit['target_id']].append(hit)
	
	besthitsDict={}
	for id in hitsDict.keys():
		scorelist=[hit['score'] for hit in hitsDict[id]]
		max_score=max(scorelist)
		max_score_index=scorelist.index(max_score)
		besthitsDict[id]=hitsDict[id][max_score_index]
	return(besthitsDict)	


def best_hit_per_target(hits,minscore=0, minlength=0):
	''' This function extracts the best exonerate hit per target sequence, so by-catch is not
included in downstream analyses (multiple alignment, treebuilding). It is the opposite
of the exonerate function -n 1, which makes exonerate output the best hit for each 
query. '''
	records=[]
	hits_by_contig=defaultdict(list)
	for hit in hits:
		hits_by_contig[hit['target_id']].append(hit)
	
	discarded=0	
	for key in hits_by_contig.keys(): ## For every target contig
		hits_in_single_contig=hits_by_contig[key]
		besthit=max(hits_in_single_contig,key=lambda hit: hit['score'])

		if besthit['score']>minscore and len(besthit['target_seq'])>minlength:
			thisRecord = SeqRecord(Seq(besthit['target_seq'],IUPAC.ambiguous_dna), \
				id=besthit['target_id'] + '_to_' + besthit['query_id'], description='')
			records.append(thisRecord)
		else:
			discarded+=1
	sys.stderr.write('Discarded ' + str(discarded) + ' sequences out of ' + str(len(hits_by_contig.keys())) + '\n')
	return(records)
