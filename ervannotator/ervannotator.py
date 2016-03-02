#!/usr/bin/python

#######################################################################
## Author: 	Thanos Mourikis
## Date: 	     04/2014
## Purpose:	ERVdb pipeline master wrapper script
## Notes:      Next step is to create a logger object for proper error
##             messaging
#######################################################################


import argparse
import sys
import os


def run_tool(parser, args):

     if args.tool == 'scan':
          import scan as tool
     elif args.tool == 'export':
          import fire_exonerate as tool
     elif args.tool == 'date':
          import LTR_dates as tool


     ##run the chosen subtool
     tool.run(parser, args)



def main():

     ## default paths for resources
     master_path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
     pf = master_path + '/resources/probes/retrovirus_pol_probes.fasta'
     rf = master_path + '/resources/probes/retrovirus_reference.fasta'
     tf = master_path + '/resources/trna/eukaryotic_tRNAs.fa'
     hf = master_path + '/resources/hmms/*.hmm'
     gf = master_path + '/resources/probes/retrovirus_gag_probes.fasta'
     ef = master_path + '/resources/probes/retrovirus_env_probes.fasta'


     #################################
     ## top-level parser
     #################################
     parser = argparse.ArgumentParser(prog='ervannotator', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

     ## subparser
     subparsers = parser.add_subparsers(title='Tools', dest='tool', metavar="")


     #################################
     ## subparsers
     #################################


     ## Genome scan
     parser_scan = subparsers.add_parser('scan', help='Scan genomes and annotate ERV sequences')
     parser_scan.add_argument('directory', help='Directory with genomes as FASTA files')
     parser_scan.add_argument('db_name', help='MySQL Database name')
     parser_scan.add_argument('-ds', '--db_server', dest='mysqldb_server', default='localhost', nargs='?', help='Specify MySQL Database server [localhost]')
     parser_scan.add_argument('-du', '--db_user', dest='mysqldb_user', default='root', nargs='?', help='Specify MySQL Database username [root]')
     parser_scan.add_argument('-dp', '--db_pass', dest='mysqldb_pass', default = '', nargs='?', help='Specify MySQL Database password [none]')
     parser_scan.add_argument('-t', '--threads', dest='threads', default='10', nargs='?', help='Specify the number of threads you want to run [10]')  
     parser_scan.add_argument('-mb', '--mblast', dest='mblast', default='no', nargs='?', help='Specify if you want to use megaBLAST utility [no]')
     parser_scan.add_argument('-pt', '--probe_type', dest='probe_type', default='prot', choices = ['prot', 'nucl'], nargs='?', help='Specify the type of probe sequences [prot]')
     parser_scan.add_argument('-rt', '--ref_type', dest='ref_type', default='prot', choices = ['prot', 'nucl'], nargs='?', help='Specify the type of reference sequences [prot]')
     parser_scan.add_argument('-pf', '--probe_fn', dest='probe_fn', help='Specify probe sequence file name [resources]', default=pf)
     parser_scan.add_argument('-rf', '--reference_fn', dest='reference_fn', help='Specify reference sequence file name [resources]', default=rf)
     parser_scan.add_argument('-tf', '--trna_fn', dest='trna_fn', help='Specify tRNA sequence file name [resources]', default=tf)
     parser_scan.add_argument('-hf', '--hmm_fn', dest='hmm_fn', help='Specify hmm directory [e.g path_to_hmms/*.hmm] [resources]', default=hf)
     parser_scan.add_argument('-gf', '--gag_fn', dest='gag_fn', help='Specify gag probe file name [resources]', default=gf)
     parser_scan.add_argument('-ef', '--env_fn', dest='env_fn', help='Specify env probe file name [resources]', default=ef)
     parser_scan.add_argument('-mhd', '--min_hsp_dist', dest='min_hsp_dist', default='5000', nargs='?', help='Specify minimum HSP distance to merge cancidate hsp hits [3000]')
     parser_scan.add_argument('-ev', '--evalue', dest='evalue', default='0.00001', nargs='?', help='Specify BLAST e-value for hit filtering [0.00001]')
     parser_scan.add_argument('-eb', '--exbuffer', dest='exbuffer', default='8000', nargs='?', help='Specify extension buffer for each side of a BLAST hit [8000]')
     parser_scan.add_argument('-ltr', '--ltr', dest='ltr', default='yes', choices=['yes', 'no'], nargs='?', help='Proceed to LTR identification using LTRharvest [yes]')
     parser_scan.add_argument('-ls', '--ltrsim', dest='ltrsim', default='70.0', nargs='?', help='Specify LTR similarity threshold for LTRharvest [70.0]')
     parser_scan.add_argument('-mt', '--mintsd', dest='mintsd', default='4', nargs='?', help='Specify minimum TSD length for LTRharvest[4]')


     parser_scan.set_defaults(func=run_tool)

     ## Exonerate tool parser
     parser_export = subparsers.add_parser('export', help='Use exonerate to export sequences from the database')
     parser_export.add_argument('db_name', help='MySQL Database name')
     parser_export.add_argument('query', help='File name of FASTA file of exonerate queries')
     parser_export.add_argument('species', help='Specify species name')
     parser_export.add_argument('-ds', '--db_server', dest='mysqldb_server', default='localhost', nargs='?', help='Specify MySQL Database server [localhost]')
     parser_export.add_argument('-du', '--db_user', dest='mysqldb_user', default='root', nargs='?', help='Specify MySQL Database username [root]')
     parser_export.add_argument('-dp', '--db_pass', dest='mysqldb_pass', default = '', nargs='?', help='Specify MySQL Database password [none]')
     parser_export.add_argument('-g','--targets', dest='targets', choices = ['pols', 'ervs'], help='Specify target level: Raw pol hits or ERVs [ervs]', default='ervs')
     parser_export.add_argument('-c','--chunks', dest='chunks', type=int, help='Divide the target FASTA file into N chunks [10]', default=10)
     parser_export.add_argument('-t','--threads', dest='threads', type=int, help='Use N threads [10]', default=10)
     parser_export.add_argument('-o','--out', dest='out', default='export_out', help='Output file name')
     
     
     parser_export.set_defaults(func=run_tool)

     ## Dating LTR tools
     parser_date = subparsers.add_parser('date', help='LTR dating')
     parser_date.add_argument('db_name', help='MySQL Database name')
     parser_date.add_argument('species', help='Specify species name')
     parser_date.add_argument('-ds', '--db_server', dest='mysqldb_server', default='localhost', nargs='?', help='Specify MySQL Database server [localhost]')
     parser_date.add_argument('-du', '--db_user', dest='mysqldb_user', default='root', nargs='?', help='Specify MySQL Database username [root]')
     parser_date.add_argument('-dp', '--db_pass', dest='mysqldb_pass', default = '', nargs='?', help='Specify MySQL Database password [none]')
     parser_date.add_argument('-r', '--rate', dest='rate', type=int, default=2.3*10**-9 , help='Specify the rate [2.3*10**-9]')
     parser_date.add_argument('-c', '--correction', dest='correction', choices=['JC', 'K2P'], default='K2P', help='Distance correction [K2P]')


     parser_date.set_defaults(func=run_tool)

     ## Basic Database stats
     parser_stats = subparsers.add_parser('stats', help='Basic database statistics')





     args = parser.parse_args()

     ## Pass the operating system to the argparse
     args.os= sys.platform


     if not (args.os.startswith('linux') or args.os.startswith('darwin')):
          print "Operating system was not recognised. Check manual for compatible platforms."
          sys.exit(1) 



     try:
          args.func(parser, args)
     except IOError, e:
          raise


if __name__ == "__main__":
     main()




