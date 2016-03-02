ERVannotator currently runs on Linux operating systems.
Main usage tips are summarized below.


**-- Dependencies**

	Our intention is to minimize dependencies in order to facilitate analaysis for empirical users.
	The only dependencies that ERVannotator asks for is MySQl and basic Python modules.

	MySQL:

	Mac:	http://dev.mysql.com/doc/refman/5.0/en/macosx-installation.html
	Linux:	http://dev.mysql.com/doc/refman/5.1/en/linux-installation.html


	Python modules required:

	biopython
	argparse
	configparser
	mysql.connector

	** You can check if you already have them installed by typing to the command line:
		$python
		> import module_name
		If there is no error you are good to go!


**-- Run**

	ERVannotator needs no installation. Just go to ERVannotator_master/ervannotator directory
	and type:
	$ ./ervannotator.py -h
	$ ./ervannotator.py tool -h 

	to find the help page for the parameters.


**-- Help**

	usage: ervannotator [-h]  ...

	optional arguments:
  		-h, --help  show this help message and exit

	Tools:
  
    	scan      Scan genomes and annotate ERV sequences
    	export    Use exonerate to export sequences from the database
    	date      LTR dating
    	stats     Basic database stats


**-- Example running instructions to be added here...**

	scan tool:

	./ervannotator.py scan -t 2 path_to_genomes name_of_db