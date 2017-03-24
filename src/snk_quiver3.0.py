##########################################################################################################
# snk.quiver3.0.py                                                                                       #
#                                                                                                        #
## Script to polish a PacBio assembly using Quiver with SMRTANALYSIS V3.0                                #
# 1. Convert bax.h5 files to bam                                                                         #
# 2. Run pbalign with each bam file.                                                                     #
# 3. Merge all the pbalign bam files output in a single bam file.                                        #
# 4. (sort/index bam and index fasta                                                                     #
# 5. run Quiver.                                                                                         #
#                                                                                                        #
## Requirements:                                                                                         #
# - pacbio assembly (from canu, falcon, etc)                                                             #
# - pacbio reads (.bax.h5 format)                                                                        #
#                                                                                                        #
## Example run:                                                                                          #
# one node, cpus=24  [1 run with 24threads]                                                              #
# (dry run) $ snakemake --snakefile snk.quiver3.0.py -j 1 --config rdir=raw assembly=assembly.fasta -np  #
#                                                                                                        #
# multi-node  [max 80 jobs at once, each one with threads=24]                                            #
# (dry run) $ snakemake -j 80 --snakefile snk.quiver3.0.py --cluster-config cluster.json                 #
#  --cluster "sbatch --partition=compute --cpus-per-task=1 --time=14-0 --job-name=snkmk --mem=10GB"      #
#  --config rdir=raw assembly=assembly.fasta -np                                                         #
#                                                                                                        #
##########################################################################################################

import subprocess,glob
from os.path import join
import os,re,sys
from Bio import SeqIO
from os.path import basename


# Globals ---------------------------------------------------------------------

#PATH of PacificBiosciences pitchfork SMRTANALYSIS3.0
SMRTloc="/apps/pitchfork"

# Full path to a folder that holds PacBio reads in bax.h5 format. Use symlinks in case separate folders.
BAX_DIR = config["rdir"]
# assemnbly folder path
ASBLY = config["assembly"]

# Regular expression matching the bax.h5 files.
SAMPLES, = glob_wildcards(join(BAX_DIR, '{sample,[^/]+}.1.bax.h5'))

# Patterns using the 'sample' wildcard.
PATTERN = '{sample}.1.bax.h5'

fnames=ffnames=[]
fnames=glob.glob(BAX_DIR+"/*.1.bax.h5")
for f in fnames:
  ffnames.append(os.path.basename(f).split(".1.bax.h5")[0])

# Rules -----------------------------------------------------------------------

rule all:
  input:
    fasta='quiverOut.fasta',
    fastq='quiverOut.fastq'

rule createBAM:
  input:
    files=join(BAX_DIR, PATTERN)
  output:
    compliantFiles='bax2bam/{sample}.subreads.bam'
  params: 
    outfiles='bax2bam/{sample}',
    SMRT=SMRTloc
  shell:"""
    source {params.SMRT}/deployment/setup-env.sh
    dname=$(dirname {input.files});
    fname=$(basename {input.files} .1.bax.h5);
    bax2bam $dname/$fname.1.bax.h5 $dname/$fname.2.bax.h5 $dname/$fname.3.bax.h5 -o {params.outfiles} --subread --pulsefeatures=DeletionQV,DeletionTag,InsertionQV,IPD,MergeQV,SubstitutionQV,PulseWidth,SubstitutionTag
  """

BAM_DIR = "bax2bam/"
SAMPLES, = glob_wildcards(join(BAM_DIR, '{sample,[^/]+}.subreads.bam'))
PATTERN = '{sample}.subreads.bam'

rule runPBALIGN:
  input:
    files=join(BAM_DIR, PATTERN)
  output:
    'bax2bam/{sample}_aligned.bam'
  params: 
    assembly=ASBLY
  shell:"""
    pbalign --nproc 12 {input.files} {params.assembly} {output}
  """

rule mergePBALIGN:
  input:
    expand("bax2bam/{sample}_aligned.bam", sample=ffnames)
  output:
    'bax2bam/all.bam'
  shell:"""
    din=$(echo {input} | sed 's/ / -in /g')
    bamtools merge -in $din -out {output};
  """

rule sortBAM:
  input:
    'bax2bam/all.bam'
  output:
    'bax2bam/all_sort.bam'
  shell:"""
    bamtools sort -in {input} -out {output}
  """

rule indexBAM:
  input:
    'bax2bam/all_sort.bam'
  output:
    'bax2bam/all_sort.bam.bai'
  shell:"""
    bamtools index -in {input} 
  """

rule indexFasta:
  input:
    ay=ASBLY,
    bai='bax2bam/all_sort.bam.bai'
  output:
    ASBLY+".fai"
  shell:"""
    samtools faidx {input.ay} 
  """
rule pbiBAM:
  input:
    bai='bax2bam/all_sort.bam.bai'
  output:
    ASBLY+".pbi"
  params: 
    assembly=ASBLY
  shell:"""
    python /apps/unit/MikheyevU/miquel/GenomicConsensus/bin/makePbi.py --referenceFasta {params.assembly} bax2bam/all_sort.bam
  """

#if pbi missing file error. http://pb-falcon.readthedocs.io/en/latest/quick_start.html
rule runQuiver:
  input:
    bam='bax2bam/all_sort.bam',
    bai='bax2bam/all_sort.bam.bai',
    fai=ASBLY+".fai"
  output:
    fasta='quiverOut.fasta',
    fastq='quiverOut.fastq'
  params: 
    assembly=ASBLY
  shell:"""
    module load gcc/4.9.2
    variantCaller -j 12 --algorithm=best {input.bam}  --referenceFilename {params.assembly} -o {output.fasta} -o {output.fastq}
  """
