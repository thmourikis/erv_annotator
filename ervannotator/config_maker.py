#!/usr/bin/python

"""
This script facilitates and automates the creation of ERVannotator configuration file.

"""
import ConfigParser

def configMaker(parser, args, master_path, genome_path, species_name):

     if args.os.startswith('linux'):
          args.os = 'linux'
     elif args.os.startswith('darwin'):
          args.os = 'mac'

     config_name = species_name + '.cfg'

     ## Path correspons to the absolute path of ervannotator folder
     config_path = master_path + '/config/' + config_name
     configfile = open(config_path, 'w')
     config = ConfigParser.RawConfigParser()

     ## General section
     config.add_section('General')
     config.set('General', 'job_name', species_name)
     config.set('General', 'organism', species_name)
     config.set('General', 'threads', args.threads)

     ## DatabaseInfo section
     config.add_section('DatabaseInfo')
     config.set('DatabaseInfo', 'db_name', args.db_name)
     config.set('DatabaseInfo', 'mysql_server', args.mysqldb_server)
     config.set('DatabaseInfo', 'mysql_username', args.mysqldb_user)
     config.set('DatabaseInfo', 'mysql_password', args.mysqldb_pass)

     ## Paths section
     config.add_section('Paths')
     config.set('Paths', 'master_path', master_path)
     config.set('Paths', 'config_path', config_path)
     config.set('Paths', 'blast_path', master_path + '/dist/' + args.os + '/blast+/bin')
     config.set('Paths', 'rm_path', master_path + '/dist/' + args.os + '/RepeatMasker/RepeatMasker')
     config.set('Paths', 'gt_path', master_path + '/dist/' + args.os + '/genometools/bin')
     config.set('Paths', 'exo_path', master_path + '/dist/' + args.os + '/exonerate/bin')
     config.set('Paths', 'probe_file_path', args.probe_fn)
     config.set('Paths', 'reference_file_path', args.reference_fn)
     config.set('Paths', 'genome_path', genome_path)
     config.set('Paths', 'trnas_path', args.trna_fn)
     config.set('Paths', 'hmms_path', args.hmm_fn)
     config.set('Paths', 'gag_file', args.gag_fn)
     config.set('Paths', 'env_file', args.env_fn)

     ## Parameters section
     config.add_section('Parameters')
     config.set('Parameters', 'megablast', args.mblast)
     config.set('Parameters', 'probe_type', args.probe_type)
     config.set('Parameters', 'reference_type', args.ref_type)
     config.set('Parameters', 'min_hsp_distance', args.min_hsp_dist)
     config.set('Parameters', 'e_value', args.evalue)
     config.set('Parameters', 'ext_buffer', args.exbuffer)
     config.set('Parameters', 'ltrharvest_sim', args.ltrsim)
     config.set('Parameters', 'mintsd', args.mintsd)



     # Writing our configuration file to 'example.cfg'  
     config.write(configfile)


