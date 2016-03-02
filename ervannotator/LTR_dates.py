#!/usr/bin/python

## A python script to estimate LTR dates

from __future__ import division
import code
import numpy
import sys
import os
import time
import csv
from collections import defaultdict,OrderedDict
from Bio import SeqIO
import math
from Bio import SeqIO,AlignIO
from StringIO import StringIO
from Bio.Emboss.Applications import NeedleCommandline
import argparse
import subprocess
from itertools import izip
import mysql.connector
from mysql.connector import errorcode


def get_os():
    
    if sys.platform.startswith('linux'):
        platform = 'linux'
    elif sys.platform.startswith('darwin'):
        platform = 'mac'

    return platform

        
## Return the alignment of the two sequences
def needle_wrapper(seq1,seq2,gapopen=10,gapextend=4):
	needle_input_1=open("needle-input-1.fa",'w')
	needle_input_2=open("needle-input-2.fa",'w')
#	code.interact(local=locals())
	
	SeqIO.write(seq1, needle_input_1, "fasta")
	SeqIO.write(seq2, needle_input_2, "fasta")	
	needle_input_1.close()
	needle_input_2.close()		
	needle_cline = NeedleCommandline(asequence="needle-input-1.fa", bsequence="needle-input-2.fa",
									  gapopen=gapopen, gapextend=gapextend, stdout=True, auto=True)# 
	stdout, stderr = needle_cline()
	if(stderr!=''):
		exit(stderr)
	align = AlignIO.read(StringIO(stdout), "emboss")
	return(align)

def estimate_nucleotide_frequencies(seq):
    seq = str(seq).replace('-','').upper() ## SeqRecord object has no replace method
    A = seq.count('A')
    C = seq.count('C')
    G = seq.count('G')
    T = seq.count('T')
    length = float(len(seq))
    return [ x/length for x in [A,C,G,T] ]

def pdistance(seq1, seq2):
    p = 0
    pairs = []
    for x in zip(seq1,seq2):
        if '-' not in x: pairs.append(x)
    #for (x,y) in zip(seq1,seq2):
    for (x,y) in pairs:
        if x != y:
            p += 1
    #length = (len(seq1) + len(seq2)) / 2
    length = len(pairs)
    return float(p) / length

def JCdistance(seq1, seq2):
    """ 
    distance = -b log(1 - p / b)
    where:
    b = 3/4
    and p = p-distance, i.e. uncorrected distance between seq1 and seq2
    """
    from math import log
    b = 0.75
    p = pdistance(seq1,seq2)
    try: d = -b * log( 1-p/b )
    except ValueError: 
        print "Tried to take log of a negative number"
        return None
    return d

def K2Pdistance(seq1,seq2):
    """
    Kimura 2-Parameter distance = -0.5 log( (1 - 2p -q) * sqrt( 1 - 2q ) )
    where:
    p = transition frequency
    q = transversion frequency
    """
    from math import log, sqrt
    pairs = []

    #collect ungapped pairs
    for x in zip(seq1,seq2):
        if '-' not in x: pairs.append(x)
        
    ts_count=0
    tv_count=0
    length = len(pairs)
    
    transitions = [ "AG", "GA", "CT", "TC"]
    transversions = [ "AC", "CA", "AT", "TA",
                      "GC", "CG", "GT", "TG" ]

    for (x,y) in pairs:
        if x+y in transitions: ts_count += 1 
        elif x+y in transversions: tv_count += 1
    
    p = float(ts_count) / length
    q = float(tv_count) / length
    try: d = -0.5 * log( (1 - 2*p - q) * sqrt( 1 - 2*q ) )
    except ValueError: 
        print "Tried to take log of a negative number"
        return None
    return d, p, q    


def distance_to_age(distance, rate):
    LTR_age=float((distance/2)/rate)
    LTR_age_mya=LTR_age/10**6
    return LTR_age_mya


## Read two LTR batch files, retrieve the sequences from the database and align them
## to calculate the distance
## Change batch files to fasta file and skip the blastcmd part
def run(parser, args):

        ## Clean terminal
    os.system('cls' if os.name == 'nt' else 'clear')

    print "\n\n##########################################################"
    print "##               ERVannotator date tool                 ##"
    print "##########################################################"

    ## Add running time indication and print basic info
    print "\n\nStart : \t\t%s" % time.ctime()
    print "Mysql DB: \t\t%s" % args.db_name
    print "Species name: \t\t%s" % args.species

    print("\n\nERV (gi_ERVstart_ERVend)" + "\t" + "Age (mya)" + "\n")

    ## Master path
    path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))

    platform = get_os()

    targetGenomeDbPath = path + '/scan_results/' + args.db_name + '/' + args.species + '/blast_db/Target_genome_blastdb'

    blast_path = path + '/dist/' + platform + '/blast+/bin/blastdbcmd'
    
    ## Connect to the database to export LTRs per ERV
    try:
        con = mysql.connector.connect(host=args.mysqldb_server, user=args.mysqldb_user, passwd=args.mysqldb_pass)
        cur = con.cursor()
        con.database = args.db_name
    except mysql.connector.Error as err:
        print "Database connection error: {}".format(err)


    select = "select ERV_annotation.Record_ID, Element_start, Element_end, lLTR_start, lLTR_end, rLTR_start, rLTR_end from ERV_annotation, Target_genome \
                where ERV_annotation.Record_ID=Target_genome.Record_ID and Organism=" + "'" + args.species + "';"


    ## Execute select statement
    try:
        cur.execute(select)
        res = cur.fetchall()
    except mysql.connector.Error as err:
        print ("ERROR - Fetching selecting items according to a pattern: {}".format(err))
        exit(1)

    ## Check point for res. If it is empty you probably have wrong species
    ## res is a list of tuples
    if len(res) == 0:
        print "No LTRs found in the database. Please check the species name."
        sys.exit(1)

    count = 0 
    for erv in res:
        gi = erv[0].split("|")[1]
        element_start = erv[1]
        element_end = erv[2]
        lLTR_start = erv[3]
        lLTR_end = erv[4]
        rLTR_start = erv[5]
        rLTR_end = erv[6]

        

        ## get record from the Target genome database
        cmd = [blast_path, '-db', targetGenomeDbPath, '-entry', gi]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        Popen_stdout,Popen_stderr = proc.communicate()

        if len(Popen_stderr)!=0 or proc.returncode!=0:
            sys.stderr.write(Popen_stderr)
            exit('Popen blastdbcmd error:')
        ## The Popen_stdout is a string. Let's have it as a fasta record!
        ## Turn the fasta string from blastdbcmd into a temporary fake file handle.
        ## and parse it with Biopython SeqIO, to then return a SeqRecord.
        fasta_handle=StringIO(Popen_stdout)
        record = SeqIO.read(fasta_handle, 'fasta')

        ##Seq1 5-LTR
        left_LTR = record[int(lLTR_start):int(lLTR_end)+1]
        left_LTR.description = record.description + ' [' + str(lLTR_start) + '-' + str(lLTR_end) + ']'
        ##Seq2 3-LTR
        right_LTR = record[int(rLTR_start):int(rLTR_end)+1]
        right_LTR.description = record.description + ' [' + str(rLTR_start) + '-' + str(rLTR_end) + ']'
        ## To be written to the dictionary as key
        ref_desc = record.description.split()[0].split("|")[1] + '_' + str(lLTR_start) + '_' + str(rLTR_end)






        ## Align them
        alignment = needle_wrapper(left_LTR,right_LTR)
        ## Access sequences
        ## As alignment[0,:],alignment[1,:]



        ## Transform similarity to age
        if args.correction == 'JC':
            distance = JCdistance(alignment[0,:],alignment[1,:])                   
            LTR_age_mya = distance_to_age(distance, args.rate)
            print(ref_desc + "\t" + str(LTR_age_mya))
        elif args.correction == 'K2P':
            distance, p, q = K2Pdistance(alignment[0,:],alignment[1,:])                   
            LTR_age_mya = distance_to_age(distance, args.rate)
            print(ref_desc + "\t" + str(LTR_age_mya))

    print "\n\nEnd : \t\t%s" % time.ctime()



''' Test code
    age_per_element = {}
    with open(ltr_batch, "r") as ltr_handle:
        
		for line in ltr_handle:
			gi, lltr_start, lltr_stop, rltr_start, rltr_stop = line.strip().split()
			cmd = [blast_path, '-db', args.blast_db, '-entry', gi]
			proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			Popen_stdout,Popen_stderr = proc.communicate()
        
			if len(Popen_stderr)!=0 or proc.returncode!=0:
				sys.stderr.write(Popen_stderr)
				exit('Popen blastdbcmd error:')
			## The Popen_stdout is a string. Let's have it as a fasta record!
			## Turn the fasta string from blastdbcmd into a temporary fake file handle.
			## and parse it with Biopython SeqIO, to then return a SeqRecord.
			fasta_handle=StringIO(Popen_stdout)
			record = SeqIO.read(fasta_handle, 'fasta')

			##Seq1 5-LTR
			left_LTR = record[int(lltr_start):int(lltr_stop)+1]
			left_LTR.description = record.description + ' [' + str(lltr_start) + '-' + str(lltr_stop) + ']'
			##Seq2 3-LTR
			right_LTR = record[int(rltr_start):int(rltr_stop)+1]
			right_LTR.description = record.description + ' [' + str(rltr_start) + '-' + str(rltr_stop) + ']'
			## To be written to the dictionary as key
			ref_desc = record.description.split()[0].split("|")[1] + '_' + str(lltr_start) + '_' + str(rltr_stop.strip())

			## Align them
			alignment = needle_wrapper(left_LTR,right_LTR)
			## Access sequences
			## As alignment[0,:],alignment[1,:]

			## Transform similarity to age
			if ref_desc not in age_per_element.keys():
				if args.distance_correction == 'JC':
					distance = JCdistance(alignment[0,:],alignment[1,:])                   
					LTR_age_mya = distance_to_age(distance)
					age_per_element[ref_desc] = LTR_age_mya
				elif args.distance_correction == 'K2P':
					distance, p, q = K2Pdistance(alignment[0,:],alignment[1,:])                   
					LTR_age_mya = distance_to_age(distance)
					age_per_element[ref_desc] = [LTR_age_mya, p, q]


    if args.distance_correction == 'K2P':
        print("ERV" + "\t" + "Age (mya)") #+ "\t" + "p" + "\t" + "q" + "\n") -- activate if you want to print p and q
        for key in age_per_element.keys():
            print(str(key) + "\t" + str(age_per_element[key][0])) #+ "\t" + str(age_per_element[key][1]) + "\t" + str(age_per_element[key][2])) 
    else:
        print("ERV" + "\t" + "Age (mya)" + "\n")
        for key in age_per_element.keys():
            print(str(key) + "\t" + str(age_per_element[key]))
'''




