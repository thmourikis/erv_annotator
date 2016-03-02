#!/usr/bin/python
import sys
import subprocess
import os
import code
import argparse
import time
import mysql.connector
from mysql.connector import errorcode
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import exonerate_parser
from collections import defaultdict
from pprint import pprint


def exo_best_hit_per_target(exo_out):
	
	## Write to output

	out_best = open(exo_out + '_best_per_target', "w")

	exonerate_output_filename=exo_out
	exonerate_output = open(exonerate_output_filename,'r')
	hits = exonerate_parser.parse_exonerate(exonerate_output)			
	best_hits = exonerate_parser.best_hit_per_target(hits, minlength=0 , minscore=0)
	query_counter=defaultdict(lambda: 0)
	for record in best_hits:
		out_best.write(record.format('fasta'))
		query_counter[record]




def run(parser, args):

	## Clean terminal
	os.system('cls' if os.name == 'nt' else 'clear')

	print "\n\n##########################################################"
	print "##              ERVannotator export tool                ##"
	print "##########################################################"

	## Add running time indication and print basic info
	print "\n\nStart : \t\t%s" % time.ctime()
	print "Mysql DB: \t\t%s" % args.db_name
	print "Species name: \t\t%s" % args.species
	print "Threads: \t\t%s" % args.threads

	## Make the export path
	path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
	timestr = time.strftime("%Y%m%d-%H%M%S")
	
	export_path = path + '/export_results/' + args.db_name
	if not os.path.exists(export_path):
		os.makedirs(export_path)
	job_path = export_path + '/' + args.species
	if not os.path.exists(job_path):
		os.makedirs(job_path)

	os.chdir(job_path)
	query_fn = args.query.split('/')[-1]
	outfile = args.out + '_' + args.targets + '_' + query_fn + '_' + timestr
	outhandle = open(outfile, "w")  

	if args.targets == 'ervs':
		targets = path + '/scan_results/' + args.db_name + '/' + args.species + '/annotation/' + args.species + '_ERVs.fasta'
	elif args.targets == 'pols':
		targets = path + '/scan_results/' + args.db_name + '/' + args.species + '/blast_results/extracted_hits.fasta'



	## Define exonerate executable
	if args.os.startswith('linux'):
		args.os = 'linux'
	elif args.os.startswith('darwin'):
		args.os = 'mac'



	ryo_string = """sugar: %S\\n%tcs""" 
	pids=[]
	completed_output_handles=[]
	commands=[]
	for i in range(1,args.chunks+1):	## exonerate mandates 1-indexing
										## and python counts funny len(range(1,2))=1
		command_str=path + '/dist/' + args.os + '/exonerate/bin/' + 'exonerate -m protein2dna:bestfit --query ' +args.query+ ' --target '+targets+ ' --exhaustive --showalignment no --showvulgar no --targetchunktotal ' +str(args.chunks)+ " --targetchunkid " + str(i)
		command=command_str.split()+['--ryo',ryo_string] 
		output_fn='exonerate.chunk.'+str(i)+'.txt'
		commands.append((command,output_fn))

	for i in range(0,args.threads):
		if len(commands)!=0: ## still stuff left to execute:
			new_command,new_output_fn=commands.pop(0)
			new_output_handle=open(new_output_fn,'w')
			new_pid=subprocess.Popen(new_command,stdout=new_output_handle,stderr=subprocess.PIPE)
			pids.append((new_pid,new_output_handle))								
	
	while len(pids)!=0:
		time.sleep(1)
		for (pid,output_handle) in pids: ## check running processes
			pid_return=pid.poll()
			## poll returns None until the process finishes, then it returns the exit code
			if(pid_return!=0 and pid_return!=None):
				exit("Error: exonerate failed (code %d): %s" % (pid_return,pid.stderr))
			elif pid_return==0: ## Successful completion
				output_handle.close()
				completed_output_handles.append(output_handle)
				pids.remove((pid,output_handle))			
				if len(commands)!=0: ## still stuff left to execute:
					new_command,new_output_fn=commands.pop(0)
					new_output_handle=open(new_output_fn,'w')
					new_pid=subprocess.Popen(new_command,stdout=new_output_handle,stderr=subprocess.PIPE)
					pids.append((new_pid,new_output_handle))				

	assert(len(completed_output_handles)==args.chunks)
	string=''
	for file_num, file in enumerate(completed_output_handles):
	#	code.interact(local=locals())
		file.close()
		read_file=open(file.name,'r')
		for line in read_file:
			## include the header in the beginning
			if line.startswith('Command') and file_num!=0:
				continue

			if line.startswith('Hostname') and file_num!=0:
				continue
			## --completed only at the end of the file
			if line.startswith('-- completed') and file_num !=args.chunks:
				continue

			string+=line

		file.close()
 		os.unlink(file.name)

	outhandle.write(string)

	exo_best_hit_per_target(outfile)

	print "End : \t\t%s" % time.ctime()
