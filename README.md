# ahcg_pipeline
Variant calling pipeline for genomic data analysis

## Requirements

1. [Python3 - version 3.4.1](https://www.python.org/download/releases/3.4.1/)
2. [Trimmomatic - version 0.36](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip)
3. [Bowtie2 - version 2.2.9](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/)
4. [Picard tools - version 2.6.0](https://github.com/broadinstitute/picard/releases/download/2.6.0/picard.jar)
5. [GATK - version 3.4](https://software.broadinstitute.org/gatk/download/)

## Reference genome

Reference genomes can be downloaded from [Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html)

## Test data

Use the following protocol to download and prepare test dataset from NIST sample NA12878

```{sh}
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
head -100000 NIST7035_TAAGGCGA_L001_R1_001.fastq > test_r1.fastq
head -100000 NIST7035_TAAGGCGA_L001_R2_001.fastq > test_r2.fastq
```

## Help

To access help use the following command:

```{sh}
python3 ahcg_pipeline.py -h
```
# **Mission**
The goal of this project is to standardize a variant calling pipeline for the detection of variants in the circulating tumor DNA (ctDNA).
ctDNA are recent biomarker which can help in deeper understanding of cancer genomics.


# **Data Acqusition Team** 
The data acqusition team used the data from the paper "Development and validation of a clinical cancer
genomic profiling test based on massively parallel DNA sequencing". The run numeber is SRR948994. This data was downloaded using the *fastq-dump* utility from the SRA Toolkit.
The command used is:
```{sh}
fastq-dump --split-files SRR948994
```
**Installing VirtualBox**

VirtualBox is useful to run more than one operating systems simultaneously. It makes software installations easy. It creates an isolated virtual environment from host OS. If anything goes wrong snapshot feature of VirtualBox can reset the VM to particular date/time. Thus, making it easier for testing and disaster recovery.

**Installation**
VirtualBox was installed on Ubuntu 16.04 operating system following instuctions on the link.
1. [VirtualBox](https://www.virtualbox.org/wiki/Linux_Downloads)

VirtualBox with name Ubuntu_VM_KP was installed on local machine.
