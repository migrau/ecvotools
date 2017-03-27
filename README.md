# Ecvotools. Ecology and Evolution Tools.

## Description 

Set of different tools for ecology and evolution analysis.
 
## Contents

* _src/multiblast.py_ . Multi-core blast. [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow.
* _src/multiblast_slurm.py_ . Multi-core blast on slurm.  [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow.

### multiblast.py - multiblast_slurm.py

Speed up a blast run, trimming the query fasta file and running it on multiple cores/nodes with multi-threads.

 _Usage_
 
 ```{bash}
 (dry run) $ snakemake --snakefile multiblast.py --config fadb=base.fasta faqy=query.fasta dbtype=nucl -j 4 -np 
 ```
 
 or, in a cluster with slurm:
  ```{bash}
 (dry run) $ snakemake --snakefile multiblast_slurm.py --config fadb=base.fasta faqy=query.fasta dbtype=nucl -j 999 --cluster-config cluster.json --cluster "sbatch -c {cluster.cpus-per-task} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem}" -np
 ```
 A config file for snakemake is included (_cluster.json_)
 
 
