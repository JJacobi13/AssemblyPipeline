The assembly snakemake pipeline
===============================

The assembly snakemake pipeline is a pipeline for the creation of an assembly 
from raw sequencing reads of any sequencing platform. 

Execution
---------
The pipeline can be run with the following command:

	snakemake --snakefile {/path/to/your/workflow} [{outputFile}]
	
Where the outputFile is required when the rule all is not defined in the workflow.

Rules
-----
The rules package only contains rules and functions required by these rules. The rules 
package contains multiple sup packages:

* assembly - all rules for assemblers are implemented in this package, each assembler is defined
in a different module
* fastqProcessing - all rules for the preprocessing of NGS data are stored in this package. For more information,
read the [README.md](rules/fastqProcessing/README.md) in this package.
* functions -  This package contains helper functions for creating the rules, this package only contains
functions which can be used in different packages!
* control - This package will contain the rules for 

Workflows
---------
The workflows implemented will be the following workflows

* Tomato150
	* Don't know preprocessing
	* Allpaths assembler and scaffolder
* IlluminaPacbioHybrid
	* Preprocessing with fastq-mcf
	* Assembly with WGS
	* Scaffolding with pbJelly
* assemblyBestPreprocess
	* adapter trimming with fastq-mcf
	* kmer correction with quake
	* PhiX contamination filtering
	* EColi contamination filtering
	* Assembly with WGS
	* Scaffolding with Scarpa
* assemblyFastPreprocess
	* Preprocessing with fastq-mcf
	* Assembly with WGS
	* Scaffolding with Scarpa
* Test
	* This workflow is changed all time for testing all new rules
	
Configuration
-------------
The configuration of the workflows is always stored in JSON format. The has to be called config.json,
this file has to be stored in the working directly. In each workflow the config has to be parsed by using
the load method of the json module, and has to be stored in a global variable called CONFIG.

The JSON configuration has to follow the schema described in configSchema.json.

	
