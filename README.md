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

**Variant Calling Pipeline designed for detecting variants from exosomal DNA of cancer patients**

Blood samples were collected from cancer patients in 3 tubes (30 ml). The blood sample was centrifuged to obtain plasma/serum. The Extracellular Vesicles were isolated by differential ultracentrifugation at >120000g >2hrs. This centrifugation step removes all blood cells and microvesicles. The nanovesicles are pelleted out.  

**Purifying HMW DNA from exosomal vesicles
Isolation of Extracellular Vesicles from blood
liquid biopsy non invasive
blood/serum
3 blood tubes approx 3ml
derive serum/plasma
isolate nanovesicles/ exosomes by differential ultracentrifugation (>120000g >2hr) 3-4 steps remove all cells and microvesicles
nanovesicles are pelleted dowm. The pellet contains High Molecular Weight DNA (ng). To isolate the HMW DNA a QUIAGEN HMW isolation kit is used. The resulting exosome pellet was resuspended in buffer solution. 

*Library Preparation*
Library for sequencing was prepared using illumina NEXTERA kit for exosomal DNA enrichment.

Run illumina hiseq 2500 rapid run to obtain 100 bp paired end sequneces. >200X coverage. 


Purify exosomal DNA
library prep
Hiseq 
100bp 
exome kit to capture all genes
use 
