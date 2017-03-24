##########################################################################################################
# snk.quiver2.3.py                                                                                       #
#                                                                                                        #
## Script to polish a PacBio assembly using Quiver with SMRTANALYSIS V2.3                                #
# 1. Create fofn files from the input bax.h5.                                                            #
# 2. Run pbalign on each fofn file.                                                                      #
# 3. Merge all the cmp.h5 files                                                                          #
# 4. Sort the cmp.h5 single file.                                                                        #
# 5. Run quiver with the assembly.                                                                       #
#                                                                                                        #
## Requirements:                                                                                         #
# - pacbio assembly (from canu, falcon, etc)                                                             #
# - pacbio reads (.bax.h5 format)                                                                        #
#                                                                                                        #
## Example run:                                                                                          #        
# one node, cpus=24  [1 run with 24threads]                                                              #
# (dry run) $ snakemake --snakefile snk.quiver2.3.py -j 1 --config rdir=raw assembly=assembly.fasta -np  #
#                                                                                                        #
# multi-node  [max 80 jobs at once, each one with threads=24]                                            #
# (dry run) $ snakemake -j 80 --snakefile snk.quiver2.3.py --cluster-config cluster.json                 #
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

#PATH of SMRTANALYSIS2.3
SMRTloc="/apps/SMRT/smrtanalysis"

# folder path that holds PacBio reads in bax.h5 format. Use symlinks in case separate folders.
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
#ffnames=list(set(ffnames))

# Rules -----------------------------------------------------------------------

rule all:
  input:
    gff='quiverOut.gff',
    fasta='quiverOut.fasta',
    fastq='quiverOut.fastq'
    

rule createFOFN:
  input:
    files=join(BAX_DIR, PATTERN)
  output:
    compliantFiles='fofn_files/{sample}.fofn'
  shell:"""
    dname=$(dirname {input.files});
    fname=$(basename {input.files} .1.bax.h5);
    echo -e $dname/$fname.1.bax.h5"\n"$dname/$fname.2.bax.h5"\n"$dname/$fname.3.bax.h5 > {output.compliantFiles};
  """

FOFN_DIR = "fofn_files/"
SAMPLES, = glob_wildcards(join(FOFN_DIR, '{sample,[^/]+}.fofn'))
PATTERN = '{sample}.fofn'

rule runPBALIGN:
  input:
    files=join(FOFN_DIR, PATTERN)
  output:
    'cmp_files/{sample}.cmp.h5'
  params: 
    assembly=ASBLY,
    SMRT=SMRTloc
  shell:"""
    source {params.SMRT}/current/etc/setup.sh;
    pbalign --forQuiver --nproc 8 {input.files} {params.assembly} {output}
  """

rule mergePBALIGN:
  input:
    expand("cmp_files/{sample}.cmp.h5", sample=ffnames),
  output:
    'cmp_files/all.cmp.h5'
  shell:"""
    cmph5tools.py merge --outFile {output} {input}
  """

rule sortCMP:
  input:
    'cmp_files/all.cmp.h5'
  output:
    'cmp_files/all.cmp_sort.h5'
  shell:"""
    cmph5tools.py sort --outFile {output} {input} 
  """

##create fasta index 'samtools faidx'

rule runQuiver:
 input:
   'cmp_files/all.cmp_sort.h5'
 output:
   gff='quiverOut.gff',
   fasta='quiverOut.fasta',
   fastq='quiverOut.fastq'
 params: 
   assembly=ASBLY
 shell:"""
   quiver -j 24 {input} -r {params.assembly} -o {output.gff} -o {output.fasta} -o {output.fastq} 
 """
