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
## **Mission**
The goal of this project is to standardize a variant calling pipeline for the detection of variants in the circulating tumor DNA (ctDNA).
ctDNA are recent biomarker which can help in deeper understanding of cancer genomics.


## **Data Acqusition Team** 
The data acqusition team used the data from the paper "Development and validation of a clinical cancer
genomic profiling test based on massively parallel DNA sequencing". The run numeber is SRR948994. 
SRAToolkit can be used to download data simply by using TOolkit command and accession number for the download. The following command downloads the paired end reads and splits them into two files storing them into *fastq* format.
This data was downloaded using the *fastq-dump* utility from the SRA Toolkit.
The command used is:
```{sh}
fastq-dump --split-files SRR948994
```
**Dataset2**
Second dataset was obtained from *Integrated digital error suppression for improved detection of circulating tumor DNA.* paper. The run number is *SRR3502999* This data was used for variant calling using pipeline version *ahcg_pipeline_v1.0.3.py*

## **Installing VirtualBox**

VirtualBox is useful to run more than one operating systems simultaneously. It makes software installations easy. It creates an isolated virtual environment from host OS. If anything goes wrong snapshot feature of VirtualBox can reset the VM to particular date/time. Thus, making it easier for testing and disaster recovery.

**Installation**
VirtualBox was installed on Ubuntu 16.04 operating system following instuctions on the link.
1. [VirtualBox](https://www.virtualbox.org/wiki/Linux_Downloads)

VirtualBox with name Ubuntu_VM_KP was installed on local machine.

```{sh}
VBoxManage import Ubuntu-64-DR-AHCG2017.ova -p 10022
VBoxmanage startvm Ubuntu-64-DR-AHCG2017 --type headless
```
## **Pipeline**

The updated version of pipeline *ahcg_pipeline_v1.0.1.py* was run using input of DNA tumor (NCI-H2126) and matched normal (NCI-H2126 BL) cell line DNA sample data obtained from run number SRR948994
```{sh}
python ahcg_pipeline_v1.0.1.py \
-t /data2/AHCG2017FALL/bin/Trimmomatic-0.36/trimmomatic-0.36.jar \
-b /data2/AHCG2017FALL/bin/bowtie2-2.2.9/bowtie2 \
-p /data2/AHCG2017FALL/bin/picard/picard.jar \
-g /data2/AHCG2017FALL/bin/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar \
-i /data2/AHCG2017FALL/data/SRR948994_1.fastq /data2/AHCG2017FALL/data/SRR948994_2.fastq \
-w /data2/AHCG2017FALL/reference_genome/Bowtie2Index/genome \
-r /data2/AHCG2017FALL/reference_genome/genome.fa \
-a /data2/AHCG2017FALL/bin/Trimmomatic-0.36/adapters/NexteraPE-PE.fa \
-o /data2/AHCG2017FALL/output \
-d /data2/AHCG2017FALL/reference_genome/GATKResourceBundle/dbsnp_146.hg38.vcf.gz
```


## **Calculate Coverage Per Gene**
Coverage per gene was calculated by GATK toolkit. A genelist was prepared by using UCSC table browser. The coverage was calculated for two genes BRAF and KRAS.

```{sh}
java -jar GenomeAnalysisTK.jar -T DepthOfCoverage --calculateCoverageOverGenes:REFSEQ /data2/users/kpatankar7/geneList.refSeq 
-R /data2/AHCG2017FALL/reference_genome/genome.fa -o doc_gene_summary -I /data2/AHCG2017FALL/output/SRR948994_1_trimmed_IR.bam 
```

## **Workflow**
Liquid Biopsy Workflow
Patient blood sample is first collected and subsequent isolation of exosome is done by differential ultracentrifugation protocols. High Molecular Weight (100 -200kb) DNA was extracted for downstream amplification and sequencing. Bioinformatics analysis was performed with a custom pipeline for variant calling. 
Exosome Isolation and HMW DNA Extraction
Blood sample from patients was collected using standard blood draw protocol. 30ml of blood was collected from each patient in 3 tubes of volume 10 ml each. Exosome from serum and plasma was isolated by performing a series of differential centrifugation cycle at centrifugal force of >120000g for more than 2 hrs. This step was repeated 3-4 times the supernatant containing the all cells and micro vesicles is discarded. The nanovesicles are pelleted down and the pellet contains nanovesicles with High Molecular Weight DNA(ng). This resulting pellet is used for HMW DNA extraction by MagAttract HMW DNA Kit by Quiagen.
Library Preparation and Sequencing
Library Preparation for downstream amplification and sequencing was done by using ILLUMINA Nextera Kit. Exome sequencing was performed by using Illumina HiSeq 2500 Rapid Run with genomic coverage of more than 200X for each sample using 100bp paired end reads. For sequencing the coverage for Illumina HiSeq 2500 rapid run is 200X with the 37 MB genome size. 6 samples can be sequenced at a time. 
Bioinformatics Analysis
Custom bioinformatics pipeline was used analyze patient exome. Raw Illumina HiSeq reads were used as input. A python script ahcg_pipeline_v1.0.4.py pipeline was used for calling variants. Trimmomatic was used to trim the adapter sequence from the Illumina HiSeq reads. The reads were aligned to human genome reference build hg38 using Bowtie. The SAM files were processed and manipulated to remove PCR duplicates using Picard. Genome Analysis Toolkit (GATK) was used to perform realignment to local sequences and indel realignment. Base recalibration is done to remove errors caused by sequencer. SNP and indel calling is done by using GATK haplotype caller. The variant calls are stored in vcf format. Coverage per gene is calculated for 70 actionable cancer gene using samtools and bcftools.




