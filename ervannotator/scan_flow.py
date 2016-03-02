#!/usr/bin/python

#######################################################################
## Genome scan workflow
#######################################################################



import os
import time
import configparser as cfg
import sys
import select
import DB
import Base as BS
import GenomeControl as GC
import BLASTxml as BXML
import RepeatMasker as RM
import LTRharvest as LH
import Annotation as AN


     



def scanFlow(parser, args, config_path):



     config = cfg.ConfigParser()
     config.read(config_path)

     ##General
     job_name = config.get('General', 'job_name')
     organism = config.get('General', 'organism')
     threads = config.get('General', 'threads')

     ##DatabaseInfo
     db_name = config.get('DatabaseInfo', 'db_name')
     mysql_server = config.get('DatabaseInfo', 'mysql_server')
     mysql_username = config.get('DatabaseInfo', 'mysql_username')
     mysql_password = config.get('DatabaseInfo', 'mysql_password')

     ##Paths
     master_path = config.get('Paths', 'master_path')
     config_path = config.get('Paths', 'config_path')
     blast_path = config.get('Paths', 'blast_path')
     rm_path = config.get('Paths', 'rm_path')
     gt_path = config.get('Paths', 'gt_path')
     probe_file_path = config.get('Paths', 'probe_file_path')
     reference_file_path = config.get('Paths', 'reference_file_path')
     genome_path = config.get('Paths', 'genome_path')
     trnas_path = config.get('Paths', 'trnas_path')
     hmms_path = config.get('Paths', 'hmms_path')
     gag_file = config.get('Paths', 'gag_file')
     env_file = config.get('Paths', 'env_file')

     ##Parameters
     megablast = config.get('Parameters', 'megablast')
     probe_type = config.get('Parameters', 'probe_type')
     reference_type = config.get('Parameters', 'reference_type')
     min_hsp_distance = config.get('Parameters', 'min_hsp_distance')
     e_value = config.get('Parameters', 'e_value')
     ext_buffer = config.get('Parameters', 'ext_buffer')
     ltrharvest_sim = config.get('Parameters', 'ltrharvest_sim')
     mintsd = config.get('Parameters', 'mintsd')

     ## Clean terminal
     os.system('cls' if os.name == 'nt' else 'clear')

     print "\n\n##########################################################"
     print "## Retroviral Detection using Database-Integrated Tools ##"
     print "##                     (ERVannotator)                   ##"
     print "##########################################################"

     ## Add running time indication and print basic info
     print "\n\nStart : %s" % time.ctime()
     print "\nJob name: \t\t%s" % job_name
     print "Mysql DB: \t\t%s" % db_name
     print "Genome: \t\t%s" % genome_path
     print "Probe file: \t\t%s" % probe_file_path
     print "Reference file: \t%s" % reference_file_path
     print "Threads: \t\t%s" % threads


     ## 1
     ## Create files, database and tables
     ## Report if database already exists
     BS.prepare(job_name, master_path, db_name)
     print "\n\n- Creating screening MySQL database\n"
     DB.prepare(job_name, mysql_server, mysql_username,
               mysql_password, master_path)
     DB.create_screening_db(db_name)
     DB.create_tables(db_name)

     ## 2
     ## Populate basic mysql tables
     GC.prepare(job_name, master_path, db_name)
     GC.insert_to_basic_tables(genome_path, organism, probe_file_path,
                              reference_file_path, db_name, job_name, probe_type,
                              reference_type, min_hsp_distance, e_value, ext_buffer, 
                              ltrharvest_sim, mintsd)

     ## 3
     ## Creating BLAST databases for the target genome and the references
     print "\n\n- Creating BLAST databases"
     GC.create_BLAST_db(blast_path, genome_path, reference_type, 
                    reference_file_path, gag_file, env_file)

     ## 4
     ## Run first round of BLAST with probes against the target genome
     GC.run_first_blast(probe_type, blast_path, probe_file_path, threads,
                    megablast)

     ## 5
     ## Extract BLAST hits based on parameters
     ## Populate BLAST_hits mysql table and Extracted_hits mysql table
     print "\n\n- Extracting buffered hits from BLAST xml"
     BXML.prepare(job_name, master_path, db_name)
     BXML.blast_hit_extraction(probe_type, min_hsp_distance, e_value, ext_buffer, 
                              blast_path, organism, db_name)
     print "\n\n Done!"

     ## 6 
     ## Run second BLAST round with extracted hits against the reference sequences
     ## Populate Extracted_hits_similarity mysql table
     GC.run_second_blast(reference_type, blast_path, threads, megablast)
     BXML.ref_similarity_extraction(reference_type, db_name)
     print "\n\n Done!"

     if args.ltr == 'yes':
          ## 7
          ## Mask LINEs, SINEs etc from extracted elements
          ## and extract masked elements to be imported by LTRharvest
          RM.prepare(job_name, master_path, db_name)
          RM.run_repeatmasker(rm_path, threads, blast_path)
          RM.mask_reps_not_LTRs(blast_path, db_name)

          ## 8
          ## Run LTRharvest on the masked extracted hits
          LH.prepare(job_name, master_path, db_name)
          LH.run_ltrharvest(gt_path, ltrharvest_sim, mintsd, threads, trnas_path, hmms_path)
          LH.ltrharvest_parsing(db_name)

          AN.prepare(job_name, master_path, db_name)
          AN.annotate(organism, db_name, master_path, blast_path, threads)


          DB.close_mysql_con()
          print "\n\nEnd : %s" % time.ctime()

     elif args.ltr == 'no':
          DB.close_mysql_con()
          print "\n\nEnd : %s" % time.ctime()




