#########################################################################################################################
## MULTI-BLAST on a single node with multiple cpus. 																										#
# Description: Split one query.fasta, and blast it against a local db (built from a fasta file). Merge results.	    	#
#                                                                                       								#
## Usage                                                                                								#
# $ snakemake --snakefile multiblast.py --config fadb=base.fasta faqy=query.fasta dbtype=nucl -j 4 -np      			#
#           																											#
# -j: number of jobs run in paralel. 																					#
#																														#
## NO-SLURM VERSION (max, 24 cpus):																						#
# Run it from a working node: sbatch -p compute -c 24 -t 14-0 --mem 10GB --pty bash										#
# Possible configuration: -j 4 -blasthreads 6 OR -j 6 -blasthreads 4 OR -j 8 -blasthreads 3 OR -j 12 -blasthreads 2... 	#
#																														#
#																														#
# Output file:	query.blast.tsv                                                            								#
#########################################################################################################################

import subprocess,sys
from os.path import join,basename

# Globals ---------------------------------------------------------------------
#database fasta file
FASTA_DB = config["fadb"]
#query fasta file
FASTA_QUERY = config["faqy"]
#output prefix from query fasta file
prefix=basename(FASTA_QUERY).split(".")[0]
#nucl or prot
dbtype=config["dbtype"]
if dbtype=="nucl":
	rblast="blastn"
else:
	rblast="blastp"
#number of blast jobs to run (split the query.fasta file). Change line37 -> %01d with number of digits.
blastJobs=5
#threads for each blast run.
blasthreads=12

SAMPLESX=[]
for i in range(0,blastJobs):
	SAMPLESX.append(prefix+'.%01d'%i)

rule all:
	input:
		prefix+'.blast.tsv'	

rule splitFASTA:
	input:
		base=FASTA_DB,
		query=FASTA_QUERY,
	output:
		trimmed=expand('blastemp/{sample}.fasta', sample=SAMPLESX),
		datab='blastemp/blastDB.nsq'
	params: 
		prefix=prefix,
		dbtype=dbtype
	shell:"""
		pyfasta split -n {blastJobs} {input.query};
		makeblastdb -in {input.base} -dbtype {params.dbtype} -out blastemp/blastDB;
		mv {params.prefix}.*.fasta blastemp/ && mv {params.prefix}.fasta.* blastemp/;
	"""

FASTA_DIR = 'blastemp/'
SAMPLES = glob_wildcards(join(FASTA_DIR, '{sample,[^/]+}.fasta'))
PATTERN = '{sample}.fasta'

rule blast:
	input:
		fastas=join(FASTA_DIR, PATTERN),
		datab='blastemp/blastDB.nsq'
	output:
		'blastemp/{sample}.tsv'
	params: 
		blasthreads=blasthreads,
		blast_type=rblast
	shell:"""
		{params.blast_type} -db blastemp/blastDB -query {input.fastas} -outfmt 6 -out {output} -num_threads {params.blasthreads}
	"""
#add -evalue 0.01

rule catblast:
	input:
		expand("blastemp/{sample}.tsv", sample=SAMPLESX)
	output:
		prefix+'.blast.tsv'
	shell:"""
		cat {input} > {output}
		rm -r blastemp/
	"""
