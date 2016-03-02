#!/usr/bin/python

#######################################################################
## Author: 	Thanos Mourikis
## Date: 	2014
## Purpose:	Python module for sequence fetching from BLAST database
#######################################################################



from StringIO import StringIO
from Bio import SeqIO
import subprocess
import sys
import os
import ConfigParser as cfg


def blastdb_fetch(blast_path, blastdb, accession):
	""" Fetch a sequence from a blastdatabase(eg: wgs, nr) using blastdbcmd, and return 
	it as a SeqRecord. Accession can be genbank or gi"""

	blastdbcmd = blast_path + '/blastdbcmd'

	args=[blastdbcmd, '-db', blastdb, '-entry', accession]
	a=subprocess.Popen(args,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	Popen_stdout,Popen_stderr=a.communicate()

	if len(Popen_stderr)!=0 or a.returncode!=0:
		sys.stderr.write(Popen_stderr)
		exit('Popen blastdbcmd error:')

	## Turn the fasta string from blastdbcmd into a temporary fake file handle.
	## and parse it with Biopython SeqIO, to then return a SeqRecord.
	fasta_handle=StringIO(Popen_stdout)
	record = SeqIO.read(fasta_handle,'fasta')
	return record