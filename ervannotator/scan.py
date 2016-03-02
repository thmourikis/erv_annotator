#!/usr/bin/python

'''
This script created for a batch submission to ERVannotator.
The same as multiple ERVannotators on manually created configuration files.

'''


import os
import sys
import glob
import argparse
import config_maker
import scan_flow



def run(parser, args):

     ## Assert that at least the positional arguments are present
     if not args.directory or not args.db_name:
          print '''ERVannotator scan program needs at least 2 arguments. 
                    Please use "ervannotator scan -h" \
                    for the positional arguments.'''
          sys.exit(1)

     for genome_path in sorted(glob.glob(args.directory + "*.fasta")):

          ## Pass master path to scan_flow and then to Annotation
          ## master_path = os.getcwd()
          master_path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
          fasta_file = genome_path.split("/")[-1]
          species_name = fasta_file.split(".")[0]

          ## Make the configuration file for each of the genomes
          config_maker.configMaker(parser, args, master_path, genome_path, species_name)

          ## Pass the config file that was made in previous step to the scan_flow
          config_name = species_name + '.cfg'
          config_path = master_path + '/config/' + config_name

          ## Run the scan
          scan_flow.scanFlow(parser, args, config_path)





