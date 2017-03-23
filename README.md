# Ecvotools. Ecology and Evolution Tools.

## Description 

Set of different tools for ecology and evolution analysis.

## Contents

* _src/snk_canupipe.py_ . Run correction of pacbio reads and performs assembly.
* _src/snk.quiver2.3.py_ . Polish canu assembly using Quiver 2.3.
* _src/snk.quiver3.0.py_ . Polish canu assembly using Quiver 3.0.

### snk_canupipe

First, it uses [colormap](https://github.com/cchauve/CoLoRMap) to run a correction of pacbio reads (long reads) using illumina reads (short reads). Since it takes long time (it runs bwa two times during the process), the script splits the pacbio reads (by default in 999 sub-files), it runs the correction in parallel and finally it merges the results.

Second, [canu](https://github.com/marbl/canu) performs the assembly.

_Considerations_

- Replace the _numJobs_-_threadsCorrection_ (colormap) and _genomeSize_ (canu) as needed.
- Script optimized to run on a slurm cluster.
- python/3.5.0 required to run the snakemake script.

_Usage_

(dry run) $ snakemake -j 60 -np --cluster "sbatch --partition=compute --cpus-per-task=12 --time=14-0 --job-name=snkmk --mem=20GB" -config pfasta=pacbio.fasta -config ifasta=illumina.fastq    

_Possible improvements_

- Add polishing step, after the canu assembly, using pacbio reads (quiver) or illumina reads (pilon).

### snk.quiver2.3

Once the assembly is obtained (from canu or falcon for example), the result can be polished using the PacBio reads and [Quiver](https://github.com/PacificBiosciences/GenomicConsensus). The _src/snk.quiver2.3.py_ script uses SMRTANALYSIS V2.3. The steps are:

 - Create fofn files from the input bax.h5.
 - Run pbalign on each fofn file.
 - Merge all the cmp.h5 files.
 - Sort the cmp.h5 single file.
 - Run quiver.
 
_Required files_

 - Raw reads in bax.h5 format
 - Canu/falcon assembly in fasta format.

_Usage_

In case we run the script in a single node, the available threads will be limited by the available cpus on the node:

(dry run) $ snakemake --snakefile snk.quiver2.3.py -j 1 --config rdir=raw_bax.h5_folder/ assembly=canu.fasta -np
 
It can be run also in multi-node mode (for example, 80 jobs at once, each one with 24 threads):

(dry run) $ snakemake -j 80 --snakefile snk.quiver2.3.py --cluster "sbatch --partition=compute --cpus-per-task=1 --time=14-0 --job-name=snkmk --mem=10GB" --config rdir=raw_bax.h5_folder/ assembly=canu.fasta -np
 
_Considerations_

Last step (the quiver run itself) has high memory demand. It took ~7 days for a ~450Mbps genome using, at least, 1T of memory.

Snakemake config file is attached.

### snk.quiver3.0


