# ecvotools

## Description 

Set of different tools for ecology and evolution analysis.

## Contents

* _src/snk_canupipe.py_ . Run correction of pacbio reads and performs assembly.

### snk_canupipe

First, it uses [colormap](https://github.com/cchauve/CoLoRMap) to run a correction of pacbio reads (long reads) using illumina reads (short reads). Since it takes long time (it runs bwa two times during the process), the script splits the pacbio reads (_by default_ in 999 sub-files), run the correction in parallel and finally merge the results.

Second, [canu](https://github.com/marbl/canu) performs the assembly.

_Considerations_

- Replace the numJobs-threadsCorrection (colormap) and genomeSize (canu) as needed.
- Script optimized to run on a slurm cluster.
- python/3.5.0 required to run the snakemake script.

_Usage_

(dry run) $ snakemake -j 60 -np --cluster "sbatch --partition=compute --cpus-per-task=12 --time=14-0 --job-name=snkmk --mem=20GB" -config pfasta=pacbio.fasta -config ifasta=illumina.fastq    


