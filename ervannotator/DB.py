#!/usr/bin/python

#######################################################################
## Author: 	Thanos Mourikis
## Date: 	2014
## Purpose:	Python module for mysql database creation
#######################################################################

from __future__ import print_function
import MySQLdb as db
import os
import mysql.connector
from mysql.connector import errorcode

def prepare(job_name, mysql_server, mysql_username,
			mysql_password, master_path):

	global job_path
	job_path = master_path + '/' + job_name


	## Make connection global to avoid reopening a thousand times
	global con
	con = mysql.connector.connect(host=mysql_server, user=mysql_username, passwd=mysql_password)
	global cur
	cur = con.cursor()


def create_screening_db(db_name):

	try:
		cur.execute(
            "CREATE DATABASE IF NOT EXISTS {} DEFAULT CHARACTER SET 'utf8'".format(db_name))
	except mysql.connector.Error as err:
		print("Failed creating database: {}".format(err))
		exit(1)




TABLES = {}

## Table with all blast hits
TABLES['BLAST_hits'] = (
"CREATE TABLE `BLAST_hits` ("
"  `Record_ID`     	varchar(255) NOT NULL default ''," ## gi - alignment.hit_id
"  `Probe_ID`    	varchar(255) NOT NULL default ''," ## blast_record.query
"  `Probe_type`    	varchar(255) NOT NULL default ''," 
"  `Frame`   	   	int(11) NOT NULL default '0',"
"  `Bit_score`     	float NOT NULL default '0',"
"  `E_value`   	   	float  NOT NULL default '0',"
"  `Subject_start` 	int(11) NOT NULL default '0',"
"  `Subject_end`   	int(11) NOT NULL default '0',"
"  `Query_start`   	int(11) NOT NULL default '0',"
"  `Query_end`     	int(11) NOT NULL default '0',"
"  `Identities`     int(11) NOT NULL default '0',"
"  `Positives`      int(11) NOT NULL default '0',"
"  `Gaps`  		   	int(11) NOT NULL default '0',"
"  `Align_len`     	int(11) NOT NULL default '0',"
#"  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,"
"  PRIMARY KEY  (`Record_ID`, `Probe_ID`, `Subject_start`, `Subject_end`, `Query_start`, `Query_end`)" ## Composite primary key
") ENGINE=InnoDB")

## Table with extracted hits based on min distance and buffer
TABLES['Extracted_hits'] = (
"CREATE TABLE `Extracted_hits` ("
"  `Record_ID`          varchar(255) NOT NULL default '',"
"  `Sequence`          	varchar(255) NOT NULL default '',"
"  `Organism`   		varchar(255) NOT NULL default '',"
"  `Blast_start`        int(100) NOT NULL default '0',"
"  `Blast_end`         	int(100) NOT NULL default '0',"
#"  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,"
"  PRIMARY KEY  (`Sequence`), UNIQUE KEY `Sequence` (`Sequence`)"
") ENGINE=InnoDB")

## Table with extracted hits and their similarity with reference screenset
TABLES['Extracted_hits_similarity'] = (
"CREATE TABLE `Extracted_hits_similarity` ("
"  `Sequence`          		varchar(255) NOT NULL default '',"
"  `Reference`     			varchar(255) NOT NULL default '',"
"  `Max_score`     			float NOT NULL default '0',"
"  `Total_score`     		float  NOT NULL default '0',"
"  `Query_start`     		int(11) NOT NULL default '0',"
"  `Query_end`     			int(11) NOT NULL default '0',"
"  `Subject_start`     		int(11) NOT NULL default '0',"
"  `Subject_end`     		int(11) NOT NULL default '0',"
"  `E_value`   	   			float  NOT NULL default '0',"
"  `Perc_identity`   	   	int(11) NOT NULL default '0',"
#"  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,"
"  PRIMARY KEY  (`Sequence`), UNIQUE KEY `Sequence` (`Sequence`)"
") ENGINE=InnoDB")

## Basic info for target genome sequences
TABLES['Target_genome'] = (
"CREATE TABLE `Target_genome` ("
"  `Record_ID`          varchar(255) NOT NULL default '',"
"  `Organism`   		varchar(255) NOT NULL default '',"
"  `Record_desc`   		varchar(255) NOT NULL default '',"
"  `Length`             int(100) NOT NULL default '0',"
#"  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,"
"  PRIMARY KEY  (`Record_ID`)"
") ENGINE=InnoDB")

## Probe information
TABLES['Probes'] = (
"CREATE TABLE `Probes` ("
"  `Probe_name`        varchar(255) NOT NULL default '',"
"  `Length`            int(100) NOT NULL default '0',"
#"  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,"
"  PRIMARY KEY  (`Probe_name`), UNIQUE KEY `Probe_name` (`Probe_name`)"
") ENGINE=InnoDB")

## Reference information
TABLES['References'] = (
"CREATE TABLE `References` ("
"  `Reference_name`    varchar(255) NOT NULL default '',"
"  `Length`            int(100) NOT NULL default '0',"
#"  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,"
"  PRIMARY KEY  (`Reference_name`), UNIQUE KEY `Reference_name` (`Reference_name`)"
") ENGINE=InnoDB")

## Basic info for the run
TABLES['Parameters'] = (
"CREATE TABLE `Parameters` ("
"  `DB_name`        	varchar(255) NOT NULL default '',"
"  `Job_name`   		varchar(255) NOT NULL default '',"
"  `Probe_file`   		varchar(255) NOT NULL default '',"
"  `Reference_file`   	varchar(255) NOT NULL default '',"
"  `Genome`   			varchar(255) NOT NULL default '',"
"  `Probe_type`   		varchar(255) NOT NULL default '',"
"  `Reference_type`   	varchar(255) NOT NULL default '',"
"  `Min_hsp_distance`   int(100) NOT NULL default '0',"
"  `E_value`   			float NOT NULL default '0',"
"  `Ext_buffer`   		int(100) NOT NULL default '0',"
"  `LTRharvest_sim`   	int(100) NOT NULL default '0',"
"  `MinTSD`   			int(100) NOT NULL default '0',"
#"  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,"
"  PRIMARY KEY  (`Job_name`), UNIQUE KEY `Job_name` (`Job_name`)"
") ENGINE=InnoDB")

## RM info for each contig
TABLES['RepeatMasker'] = (
"CREATE TABLE `RepeatMasker` ("
"  `Sequence`          	varchar(255) NOT NULL default '',"
"  `Rep_start`          int(100) NOT NULL default '0',"
"  `Rep_end`            int(100) NOT NULL default '0',"
"  `Rep_type`   		varchar(255) NOT NULL default '',"
#"  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,"
"  PRIMARY KEY  (`Sequence`, `Rep_start`, `Rep_end`), UNIQUE KEY (`Sequence`, `Rep_start`, `Rep_end`)"
") ENGINE=InnoDB")

## ERV table
TABLES['ERV_table'] = (
"CREATE TABLE `ERV_table` ("
"  `Sequence`          	varchar(255) NOT NULL default '',"
"  `Element_start`      int(100) NOT NULL default '0',"
"  `Element_end`       	int(100) NOT NULL default '0',"
"  `Element_length`    	int(100) NOT NULL default '0',"
"  `lLTR_start`      	int(100) NOT NULL default '0',"
"  `lLTR_end`       	int(100) NOT NULL default '0',"
"  `lLTR_length`    	int(100) NOT NULL default '0',"
"  `rLTR_start`      	int(100) NOT NULL default '0',"
"  `rLTR_end`       	int(100) NOT NULL default '0',"
"  `rLTR_length`    	int(100) NOT NULL default '0',"
"  `lTSD_start`      	int(100) NOT NULL default '0',"
"  `lTSD_end`       	int(100) NOT NULL default '0',"
"  `lTSD_motif`    		varchar(100) NOT NULL default '',"
"  `rTSD_start`      	int(100) NOT NULL default '0',"
"  `rTSD_end`       	int(100) NOT NULL default '0',"
"  `rTSD_motif`    		varchar(100) NOT NULL default '',"
"  `PPT_start`      	int(100) NOT NULL default '0',"
"  `PPT_end`       		int(100) NOT NULL default '0',"
"  `PPT_motif`    		varchar(100) NOT NULL default '',"
"  `PPT_strand`    		varchar(100) NOT NULL default '',"
"  `PBS_start`      	int(100) NOT NULL default '0',"
"  `PBS_end`       		int(100) NOT NULL default '0',"
"  `PBS_strand`    		varchar(100) NOT NULL default '',"
"  `Domains`    		longtext NOT NULL default '',"
#"  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,"
"  PRIMARY KEY  (`Sequence`, `Element_start`, `Element_end`), UNIQUE KEY (`Sequence`, `Element_start`, `Element_end`)"
") ENGINE=InnoDB")

## ERV annotation table
TABLES['ERV_annotation'] = (
"CREATE TABLE `ERV_annotation` ("
"  `Record_ID`          varchar(255) NOT NULL default '',"
"  `Element_start`   	int(100) NOT NULL default '0',"
"  `Element_end`   		int(100) NOT NULL default '0',"
"  `lLTR_start`   		int(100) NOT NULL default '0',"
"  `lLTR_end`   		int(100) NOT NULL default '0',"
"  `rLTR_start`   		int(100) NOT NULL default '0',"
"  `rLTR_end`   		int(100) NOT NULL default '0',"
"  `Gag_start`   		int(100) NOT NULL default '0',"
"  `Gag_end`   			int(100) NOT NULL default '0',"
"  `Gag_frame`    		varchar(100) NOT NULL default '',"
"  `Gag_hit`    		varchar(100) NOT NULL default '',"
"  `Gag_perc_identity`  int(100) NOT NULL default '0',"
"  `Pol_start`   		int(100) NOT NULL default '0',"
"  `Pol_end`   			int(100) NOT NULL default '0',"
"  `Pol_frame`    		varchar(100) NOT NULL default '',"
"  `Pol_hit`    		varchar(100) NOT NULL default '',"
"  `Pol_perc_identity`  int(100) NOT NULL default '0',"
"  `Env_start`   		int(100) NOT NULL default '0',"
"  `Env_end`   			int(100) NOT NULL default '0',"
"  `Env_frame`    		varchar(100) NOT NULL default '',"
"  `Env_hit`    		varchar(100) NOT NULL default '',"
"  `Env_perc_identity`  int(100) NOT NULL default '0',"
#"  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,"
"  PRIMARY KEY  (`Record_ID`, `Element_start`, `Element_end`), UNIQUE KEY (`Record_ID`, `Element_start`, `Element_end`)"
") ENGINE=InnoDB")


## Create BLAST table
def create_tables(db_name):

	con.database = db_name

	for name, ddl in TABLES.iteritems():
		try:
			print ("Creating table {}: ".format(name), end='')
			cur.execute(ddl)
		except mysql.connector.Error as err:
			if err.errno == errorcode.ER_TABLE_EXISTS_ERROR:
				print("already exists.")
			else:
				print(err.msg)
		else:
			print("OK")


## Insert to mysql table

def add_to_target_table(db_name, target_dict):
	
	con.database = db_name

	add_targets = ("insert into Target_genome "
              "(Record_ID, Organism, Record_desc, Length) "
              "VALUES (%(Record_ID)s, %(Organism)s, %(Record_desc)s, %(Length)s)")

	try:
		cur.execute(add_targets, target_dict)
		# Make sure data is committed to the database
		con.commit()
		
	except mysql.connector.Error as err:
		print ("ERROR - Inserting into Target_genome table failed: {}".format(err))
		exit(1)

def add_to_probes_table(db_name, probe_dict):

	con.database = db_name

	add_probes = ("insert ignore into Probes "
              "(Probe_name, Length) "
              "VALUES (%(Probe_name)s, %(Length)s)")

	try:
		cur.execute(add_probes, probe_dict)
		con.commit()
	except mysql.connector.Error as err:
		print ("ERROR - Inserting into Probes table failes: {}".format(err))
		exit(1)

def add_to_references_table(db_name, reference_dict):

	con.database = db_name

	add_references = ("insert ignore into `References` "
              "(Reference_name, Length) "
              "VALUES (%(Reference_name)s, %(Length)s)")

	try:
		cur.execute(add_references, reference_dict)
		con.commit()
	except mysql.connector.Error as err:
		print ("ERROR - Inserting into References table failes: {}".format(err))
		exit(1)

def add_to_parameters_table(db_name, parameters_dict):

	con.database = db_name

	add_parameters = ("insert ignore into Parameters "
              "(DB_name, Job_name, Probe_file, Reference_file, Genome, Probe_type, \
              	Reference_type, Min_hsp_distance, E_value, Ext_buffer, LTRharvest_sim, MinTSD) "
              "VALUES (%(DB_name)s, %(Job_name)s, %(Probe_file)s, %(Reference_file)s, %(Genome)s, \
              	%(Probe_type)s, %(Reference_type)s, %(Min_hsp_distance)s, %(E_value)s, %(Ext_buffer)s, \
              	%(LTRharvest_sim)s, %(MinTSD)s )")

	try:
		cur.execute(add_parameters, parameters_dict)
		con.commit()
	except mysql.connector.Error as err:
		print ("ERROR - Inserting into Parameters table failes: {}".format(err))
		exit(1)

def add_to_extracted_hits_table(db_name, Record_ID, Sequence, Organism, Blast_start, Blast_end):
	
	con.database = db_name

	add_extracted_hits = ("insert ignore into Extracted_hits "
              	"(Record_ID, Sequence, Organism, Blast_start, Blast_end) "
              	"VALUES ('%s', '%s', '%s', '%s', '%s')") % (Record_ID, Sequence, Organism, Blast_start, Blast_end)


	try:
		cur.execute(add_extracted_hits)
		con.commit()
	except mysql.connector.Error as err:
		print ("ERROR - Inserting into Extracted_hits table failed: {}".format(err))
		exit(1)

def add_to_blast_table(db_name, blast_hit_dict):
	
	con.database = db_name

	add_blast_hits = ("insert into BLAST_hits "
              	"(Record_ID, Probe_ID, Probe_type, Frame, Bit_score, E_value, Subject_start, Subject_end, \
              		Query_start, Query_end, Identities, Positives, Gaps, Align_len) "
              	"VALUES (%(Record_ID)s, %(Probe_ID)s, %(Probe_type)s, %(Frame)s, \
              		%(Bit_score)s, %(E_value)s, %(Subject_start)s, %(Subject_end)s, %(Query_start)s, \
              		%(Query_end)s, %(Identities)s, %(Positives)s, %(Gaps)s, %(Align_len)s)")

	try:
		cur.execute(add_blast_hits, blast_hit_dict)
		con.commit()
	except mysql.connector.Error as err:
		print ("ERROR - Inserting into BLAST_hits table failed: {}".format(err))
		exit(1)

def add_to_extracted_hits_similarity_table(db_name, blast_hit_dict):
	
	con.database = db_name

	add_extracted_hits_similarity = ("insert ignore into Extracted_hits_similarity "
              	"(Sequence, Reference, Max_score, Total_score, Query_start, Query_end, \
              		Subject_start, Subject_end, E_value, Perc_identity) "
              	"VALUES (%(Sequence)s, %(Reference)s, %(Max_score)s, %(Total_score)s, %(Query_start)s, \
              		%(Query_end)s, %(Subject_start)s, %(Subject_end)s, %(E_value)s, %(Perc_identity)s)")


	try:
		cur.execute(add_extracted_hits_similarity, blast_hit_dict)
		con.commit()
	except mysql.connector.Error as err:
		print ("ERROR - Inserting into Extracted_hits_similarity table failed: {}".format(err))
		exit(1)

def add_to_erv_table(db_name, erv_dict):
	
	con.database = db_name

	add_erv_table = ("insert ignore into ERV_table "
              	"(Sequence, Element_start, Element_end, Element_length, lLTR_start, \
              		lLTR_end, lLTR_length, rLTR_start, rLTR_end, rLTR_length, lTSD_start, lTSD_end, lTSD_motif, \
              		rTSD_start, rTSD_end, rTSD_motif, PPT_start, PPT_end, PPT_motif, PPT_strand, PBS_start, PBS_end, \
              		PBS_strand, Domains) "
              	"VALUES (%(Sequence)s, %(Element_start)s, %(Element_end)s, %(Element_length)s, \
              		%(lLTR_start)s, %(lLTR_end)s, %(lLTR_length)s, %(rLTR_start)s, %(rLTR_end)s, %(rLTR_length)s, %(lTSD_start)s, \
              		%(lTSD_end)s, %(lTSD_motif)s, %(rTSD_start)s, %(rTSD_end)s, %(rTSD_motif)s, %(PPT_start)s, %(PPT_end)s, %(PPT_motif)s, \
              		%(PPT_strand)s, %(PBS_start)s, %(PBS_end)s, %(PBS_strand)s, %(Domains)s)")


	try:
		cur.execute(add_erv_table, erv_dict)
		con.commit()
	except mysql.connector.Error as err:
		print ("ERROR - Inserting into ERV table failed: {}.".format(err))
		exit(1)

def add_to_erv_annotation(db_name, annot_dict):
	
	con.database = db_name

	add_erv_annot = ("insert ignore into ERV_annotation "
              	"(Record_ID, Element_start, Element_end, lLTR_start, lLTR_end, \
              		rLTR_start, rLTR_end, Gag_start, Gag_end, Gag_frame, Gag_hit, Gag_perc_identity, Pol_start, \
              		Pol_end, Pol_frame, Pol_hit, Pol_perc_identity, Env_start, Env_end, Env_frame, Env_hit, Env_perc_identity) "
              	"VALUES (%(Record_ID)s, %(Element_start)s, %(Element_end)s, %(lLTR_start)s, \
              		%(lLTR_end)s, %(rLTR_start)s, %(rLTR_end)s, %(Gag_start)s, %(Gag_end)s, %(Gag_frame)s, \
              		%(Gag_hit)s, %(Gag_perc_identity)s, %(Pol_start)s, %(Pol_end)s, %(Pol_frame)s, %(Pol_hit)s, \
              		%(Pol_perc_identity)s, %(Env_start)s, %(Env_end)s, %(Env_frame)s, %(Env_hit)s, %(Env_perc_identity)s)")


	try:
		cur.execute(add_erv_annot, annot_dict)
		con.commit()
	except mysql.connector.Error as err:
		print ("ERROR - Inserting into ERV annotation failed: {}.".format(err))
		exit(1)

def add_to_rm_table(db_name, rm_dict):
	
	con.database = db_name

	add_rm = ("insert ignore into RepeatMasker "
              	"(Sequence, Rep_start, Rep_end, Rep_type) "
              	"VALUES (%(Sequence)s, %(Rep_start)s, %(Rep_end)s, %(Rep_type)s)")

	try:
		cur.execute(add_rm, rm_dict)
		con.commit()
	except mysql.connector.Error as err:
		print ("ERROR - Inserting into RM table failed: {}.".format(err))
		exit(1)

def select_contigs(db_name, field_list, table, organism):

	con.database = db_name
	columns = ",".join(field_list)

	select = "select " + columns + " from " + table + ' where Organism=' + "'" + organism + "'" + ';'

	try:
		cur.execute(select)
		res = cur.fetchall()
		return res
	except mysql.connector.Error as err:
		print ("ERROR - Select statement: {}".format(err))
		exit(1)

def select_ervs(db_name, field_list, table_list, organism):

	con.database = db_name
	columns = ",".join(field_list)
	tables = ",".join(table_list)

	select = "select " + columns + " from " + tables \
				+ " where ERV_table.Sequence=Extracted_hits.Sequence and Extracted_hits.Organism="\
				+ "'" + organism + "'" + ';'

	try:
		cur.execute(select)
		res = cur.fetchall()
		return res
	except mysql.connector.Error as err:
		print ("ERROR - ERV fetching: {}".format(err))
		exit(1)

def select_pattern(db_name, table, column, pattern):

	con.database = db_name

	select = 'select * from ' + table + ' where ' + column + " like '%" + pattern + "%';"

	try:
		cur.execute(select)
		res = cur.fetchall()
		return res
	except mysql.connector.Error as err:
		print ("ERROR - Fetching selecting items according to a pattern: {}".format(err))
		exit(1)

def select_exact_erv(db_name, ref_name, element_start):

	con.database = db_name

	select = 'select * from ERV_table where Sequence like ' + "'%" + ref_name + "%'" + ' and Element_start=' + "'" + element_start + "';"

	try:
		cur.execute(select)
		res = cur.fetchall()
		return res
	except mysql.connector.Error as err:
		print ("ERROR - Fetching selecting items according to a exact erv: {}".format(err))
		exit(1)

## More defs to be added for summarising tables and presenting results
'''
def erv_summary(): ## Not yet implemented in wrapper

	con.database = db_name

	summary = 'create view Summary as select ERV_table.Sequence, Reference, Perc_identity, Subject_start, Subject_end, Query_start, Query_end, Element_start, Element_end, Domains from ERV_table, Extracted_hits_similarity where ERV_table.Sequence=Extracted_hits_similarity.Sequence order by Reference;'

	try:
		cur.execute(summary)
		con.commit()
	except mysql.connector.Error as err:
		print ("ERROR - Creating Summary view failed: {}".format(err))
		exit(1)
'''
def close_mysql_con():

	if con:
		con.commit()
		cur.close()
		con.close()



