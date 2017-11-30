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

## Help

To access help use the following command:

```{sh}
python3 ahcg_pipeline.py -h
```
## **Mission**
The goal of this project is to standardize a variant calling pipeline for the detection of variants in the circulating tumor DNA (ctDNA).
ctDNA are recent biomarker which can help in deeper understanding of cancer genomics.


## **Data Acqusition Team** 
As part of the data acquisition team different dataset were obtained from various sources. The datasets were downloaded using SRAToolkit.
SRAToolkit was used to download data simply by using Toolkit command and accession number for the download. This data was downloaded using the *fastq-dump* utility from the SRA Toolkit. The aforementioned utility downloads the paired end reads and splits them into two files storing them into *fastq* format.

**Dataset1**
The first dataset belonged to the paper "Development and validation of a clinical cancer genomic profiling test based on massively parallel DNA sequencing". The run numeber is SRR948994. This data was processed using variant calling pipeline version *ahcg_pipeline_v1.0.3.py*

Command:
```{sh}
fastq-dump --split-files SRR948994
```
**Dataset2**
Second dataset was obtained from *Integrated digital error suppression for improved detection of circulating tumor DNA.* paper. The run number is *SRR3502999* This data was used for variant calling using pipeline version *ahcg_pipeline_v1.0.3.py*

Command:
```{sh}
fastq-dump --split-files SRR3502999
```

**Dataset3**
Third dataset was obtained from control and cancer patient. The sequence data is derived from exome sequencing of High Molecular Weight DNA isolated from nanovesicles of both control as well as cancer patient. MenCo002DNA and MenPa004DNA. The sequencing machine used was Illumina HiSeq 2500 with libraries for exome prepared using Nextera Rapid Capture Exome kit.
The dataset was provided in the form of two *fastq* files one for control and other for cancer

**Dataset4**
Dataset4 was obtained from the study *Analysis of ctDNA of the cerebrospinal fluid* The data was obtained from NCBI SRA project, The SRA run numbers SRR2530741 (BMLC2_Germline) and SRR2530742 (BMLC2_Tumor). The data was downlaoded using fastq-dump utility. The data is exome sequencing data with sequencing for samples done using Illumina HiSeq 2000 and exome capture was done using Nextera Rapid Capture Exome Kit. 

Command:
```{sh}
fastq-dump --split-files SRR2530741
fastq-dump --split-files SRR2530742
```
**Dataset5**
Dataset5 was obtained from the paper titled *Cerebrospinal fluid-derived circulating tumour DNA better represents the genomic alterations of brain tumours than plasma* published in Nature Communications journal. The dataset used wet biology protocol slightly similar to the dataset3 with focus on whole exome sequencing from the CSF ctDNA using Nextera Rapid Capture Exome Kit and Illumnia HiSeq 2000 as sequencing machine. The samples are exome sequencing data of GBM1 patient both tumor and germline DNA.

Command:
```{sh}
fastq-dump --split-files SRR1654210
fastq-dump --split-files SRR1654220
```



## **Pipeline**

The latest version of pipeline *ahcg_pipeline.py* is version 1.0.8 and can be found at the path /data2/AHCG2017FALL/bin/pipeline This version has features as follows:
1. Trimming Illumina sequence data with Trimmomatic
2. Aligning reads to the reference genome usign Bowtie2
3. Copy Number Alteration detection by Control-FREEC tool
4. Running pipeline using SRR accession number
5. Filtering Variants using Qual>30 and DP>=25
6. Calculating avrage, median and max coverage per gene


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

**filter variants based on DP and QUAL**
```{sh}
java -jar GenomeAnalysisTK.jar -T SelectVariants -R /data2/AHCG2017FALL/reference_genome/genome.fa -V /data2/AHCG2017FALL/output3/variants.vcf -select "DP >= 25 && QUAL >=30" -o /data2/AHCG2017FALL/gene_summary/filtered_vcf.vcf

wc -l filtered_vcf.vcf
50829 filtered_vcf.vcf
```
**ahcg_pipeline_v1.0.5.py**

This pipeline version which can be found on /data2/AHCG2017FALL/bin/pipeline is an updated version in which the variants are filtered on the basis of coverage per gene and quality scores. The pipeline was updated by including GATK SelectVariant command which filteres the final variants_vcf.vcf file to generate filtered variants with coverage per gene of more than 25 and quality score more than 30. 

The High Molecular Weight DNA was isolated from nanovesicles and sequenced using Illumina 2500 rapid run. This pipeline was used for variant calling for two samples MenCo002DNA and MenPa004DNA. The former sample is control and the later is from a cancer patient. This protocol of non-invasive liquid biopsy and isolation of nanovsicles to obtain High Molecular Weight DNA which is a biomarker of cancer incorporated with bioinformatics pipeline to detect variants can be used for early diagnosis of cancer.


## **Installing VirtualBox**

VirtualBox is useful to run more than one operating systems simultaneously. It makes software installations easy. It creates an isolated virtual environment from host OS. If anything goes wrong snapshot feature of VirtualBox can reset the VM to particular date/time. Thus, making it easier for testing and disaster recovery.

**Installation**
VirtualBox was installed on Ubuntu 16.04 operating system following instuctions on the link.
1. [VirtualBox](https://www.virtualbox.org/wiki/Linux_Downloads)

VirtualBox with name Ubuntu_VM_KP was installed on local machine.

```{sh}
VBoxManage import Ubuntu-64-DR-AHCG2017.ova -p 10027
VBoxmanage startvm Ubuntu-64-DR-AHCG2017 --type headless
```
Login into the virtualbox was done using the command:

```{sh}
ssh vannberglab@localhost -p 10027
```

The data from gpuvannberg server was transferred using command:

```{sh}
scp -r -P 10027 /data2/AHCG2017FALL vannberglab@localhost:~/
```



