# Nextflow Implementation of the Blanco-Melo Lab Steady State (RNA-seq) Pipeline

### Usage

#### Rhino-Specific Usage Requirements
##### Primary Run Settings

This version of RNA-seq nextflow is implemented for Fred Hutch Scientific Computing server Rhino. The environment is configured on the `scratch` directory, under `/fh/scratch/delete30/blancomelo_d/RNAseq_pipeline`. This is where you would find all the dependencies and the nextflow pipeline configurations. Please do not edit files in this directory unless you know what you are doing.

To run the nextflow pipeline to fit your specific need, you might want to create a unique config file. You should adjust the config files based on the default `conf/rhino_grch38.config`, to ensure the proper paths are set for genome reference files and other executables. Variable names should hopefully be self-explanatory. It's generally a good idea to keep separate configuration files for samples using different reference genomes, and different organisms.

##### Setting up the environment

This pipeline requires a number of packages, which are already installed within a shared directory `/fh/scratch/delete30/blancomelo_d/RNAseq_pipeline`. You can configure the environment by running the following commands:

```
WD=/fh/scratch/delete30/blancomelo_d/RNAseq_pipeline
module load Python/3.8.2-GCCcore-9.3.0 nextflow
alias java=${WD}/opt/jdk-18.0.2.1/bin/java
export PATH=${WD}/opt/jdk-18.0.2.1/bin:$PATH
export JAVA_HOME=${WD}/opt/jdk-18.0.2.1/
export EBROOTJAVA=${WD}/opt/jdk-18.0.2.1/
export PATH=$PATH:${WD}/opt/STAR-2.7.10a/bin/Linux_x86_64_static
for d in $WD/opt; do export PATH="$PATH:$d"; done
for d in $WD/opt/*; do export PATH="$PATH:$d"; done
for d in $WD/opt/*/bin; do export PATH="$PATH:$d"; done
for d in $WD/opt/python_packages/*/bin; do export PATH="$PATH:$d"; done
for d in ${WD}/opt/python_packages/*; do export PYTHONPATH="$PYTHONPATH:$d"; done
```

##### Running Nextflow Using an sbatch script

Now with everything set up, you are ready to run Nextflow. The basic command for running the pipeline requires the pipeline file `main.nf`, the profile file (e.g. `rhino_hg38`) and your sequencing reads, in either the compressed fastq format (`fastq.gz`) or as a SRA accession number (SRR12345). To see a full list of options and pipeline version, enter:
```
nextflow run main.nf -profile slurm_grch38 --help
```

The best way to run Nextflow is using an sbatch command to submit the whole pipeline as a separate computing job with allocated node, CPUs and memory. Nextflow does a great job at keeping logs of every transaction, anyway, should you lose access to the console. The memory requirements do not exceed 8GB, so you do not need to request more RAM than this. Also 1 node with 16 CPUs will be more than sufficient to process the sequencing data. **Important: Nextflow put all of its intermediate files and output into your current directory, so create a working directory and call the nextflow pipeline from that directory.**

Example commands to be used on Fred Hutch Rhino node are listed below.
```
# To download and analyze single end RNAseq
SRR=SRR19572981
sbatch  -N 1 -n 1 -c 16 \
        --job-name="nextflow_se" \
        --error=./%x_%j.err --output=./%x_%j.out \
        --wrap="nextflow ${WD}/RNAseq_NextFlow/main.nf -profile rhino_hg38 --sras $SRR --singleEnd --salmon --count"

# To download and analyze paired end RNAseq
SRR=SRR19795679
sbatch  -N 1 -n 1 -c 16 \
        --job-name="nextflow_pe" \
        --error=./%x_%j.err --output=./%x_%j.out \
        --wrap="nextflow ${WD}/RNAseq_NextFlow/main.nf -profile rhino_hg38 --sras $SRR --salmon --count"
```

In addition, the following command can be used to process a list of SRA accession number in parallel. Simple provide the accession numbers, each number as a new line, in a text file (`srr.txt`). And run the following bash scripts:

```
SRR_FILE=srr.txt
for line in $(cat ${SRR_FILE}); do
  mkdir ${line}; cd ${line}
  if [ $(fastq-dump -X 1 -Z --split-spot $line | wc -l) -eq 4 ]
  then
    sbatch  -N 1 -n 1 -c 16 \
            --job-name="nextflow_rnaseq_${line}" \
            --error=./%x_%j.err --output=./%x_%j.out \
            --wrap="nextflow ${WD}/RNAseq_NextFlow/main.nf -profile rhino_hg38 --sras $line --singleEnd --salmon --count"
  else
    sbatch  -N 1 -n 1 -c 16 \
            --job-name="nextflow_rnaseq_${line}" \
            --error=./%x_%j.err --output=./%x_%j.out \
            --wrap="nextflow ${WD}/RNAseq_NextFlow/main.nf -profile rhino_hg38 --sras $line --salmon --count"
  fi
  cd ..
done
```

## Additional arguments for customization

**Required Arguments**

| Arugment  | Usage                            | Description                                                          |
|-----------|----------------------------------|----------------------------------------------------------------------|
| -profile  | \<base,slurm\>                    | Configuration profile to use.                                       |
| --fastqs  | \</project/\*\_{R1,R2}\*.fastq.gz\> | Directory pattern for fastq files (gzipped).                      |
| --sras    | \</project/\*.sra\>              | Directory pattern for sra files.                                     |
| --workdir | \</project/tmp/\>                | Nextflow working directory where all intermediate files are saved.   |
| --email   | \<EMAIL\>                        | Where to send workflow report email.                                 |

**Save Options**

| Arguments  | Usage         | Description                                               |
|------------|---------------|-----------------------------------------------------------|
| --outdir   | \</project/\> | Specifies where to save the output from the nextflow run. |
| --savefq   |               | Compresses and saves raw fastq reads.                     |
| --saveTrim |               | Compresses and saves trimmed fastq reads.                 |
| --saveAll  |               | Compresses and saves all fastq reads.                     |
| --skipBAM  |               | Skips saving BAM files (only save CRAM). Default=False    |
| --savebw   |               | Save normalized BigWig files for UCSC genome broswer.     |
| --savebg   |               | Saves concatenated pos/neg bedGraph file.                 |

**Input File Options**

| Arguments    | Usage       | Description                                                                  |
|--------------|-------------|------------------------------------------------------------------------------|
| --singleEnd  |             | Specifies that the input files are not paired reads (default is paired-end). |
| --flip       |             | Reverse complements each strand. Necessary for some library preps.           |
| --flipR2     |             | Reverse complements R2 only (will not work in singleEnd mode).               |

**QC Options**

| Arguments       | Usage       | Description                                             |
|-----------------|-------------|---------------------------------------------------------|
| --skipMultiQC   |             | Skip running MultiQC.                                   |
| --skipRSeQC     |             | Skip running RSeQC.                                     |

**Analysis Options**

| Arguments       | Usage       | Description                                                                         |
|-----------------|-------------|-------------------------------------------------------------------------------------|
| --count       |               | Count aligned reads over RefSeq gene file using FeatureCounts |
| --salmon       |               | Count trimmed reads over transcriptome using salmon |

### Author

* Qing Yang <qyang@fredhutch.org>
