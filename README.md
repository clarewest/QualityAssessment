## Feature Gathering for Quality Assessment 

These are a collection of scripts I have been using to assess the quality of models output from SAINT2 structural prediction.

The directory structure is:
`$QA/data/DATASET/TARGET/MODELS`

Where DATASET is the name of a set of targets, and TARGET is a directory containing all the models for that target.

The following files are required in the data/DATASET directory: 
Many of these files are generated using the commands in $QA/bin/parallelize.sh 

- TARGET/MODELS            : Ensure they all have a chain id. 
  - You can add the chain id "A" using the script:
  - `$QA/bin/replace_chains.sh`

- TARGET.fasta             : sequence file
- TARGET.fasta.txt         : sequence file 
  - NB this is a duplicate as difference software may require either format
- TARGET.aln               : alignment file (from contact prediction using metaPSICOV)
- TARGET.metapsicov.stage1 : metapsicov output (predicted contacts)
- TARGET.fasta.ss2         : (for EigenTHREADER) psipred output file (predicted SS) 
  - not the one used by metapsicov!
  - see `$QA/bin/parallelize.sh`

- TARGET.lst               : A list of model names
- TARGET.pcons_lst         : A list of model names including the directory
- TARGET.proq_lst          : A list of model names including directory and .pdb extension
  - you can generate these using `$QA/bin/make_lsts.sh TARGET`

- TARGET.con               : predicted contact file for use in SAINT2 scoring
- config_TARGET file       : SAINT2 configuration file with scoring weights (ensure the paths are correct)


Run `get_features.sh` with the relevant arguments to calculate each stage.

## Running get_features.sh

Run with the flag `--help` to get help:

``` console 
$ bash get_features.sh --help

./get_features.sh [options] <NAME>

 Gather features for RFQAmodel

 Positional arguments:
 	-t | --target <NAME> 	prefix of the input files

 Optional arguments:
 	-h | --help		show this help text
 	-d | --dir		directory containing models (default: target)
 	-c | --chain		chain identifier of native structure (default: final character
 	-p | --ppv		calculate PPV of models
 	-m | --mapalign		run mapalign on models
 	-pe | --prepeigen	prepare files for EigenTHREADER
 	-s | --saint2		calculate SAINT2 score for all models (requires config file)
 	-e | --eigenthreader	run EigenTHREADER on models and gather scores
 	-p | --pcons		run Pcons for all models and gather scores
 	-q | --proq		run ProQ3D for all models and gather scores 
 	-m | --tmscore		calculate TM-score for all models
 	-l | --lddt		calculate lddt for all models
 	-b | --beff		calculate Beff for the target (requires TARGET.aln)
 	-ge | --gathereigen	gather EigenTHREADER scores from output files
 	-gp | --gatherpcons	gather Pcons scores from output files
 	-gq | --gatherproq	gather ProQ3D scores from output files
 	--scaffold <SEG> <TERM> flag as ScafFold and define segment length and terminal sampled
 	--gathersampled	gather ProQ3D, LDDT and PCons for sampled region only (requires --scaffold)


```

## Output 
During the process, the following files will be created:

- TARGET.out               : EigenTHREADER output
  - note that model names must be added with $QA/bin/paste_eig_lst.sh
  - as they are truncated in the output of EigenTHREADER 

- TARGET.pcons.txt         : Pcons output



