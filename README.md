# CAST

## Introduction

**CAST** (**C**orrection **A**nd **S**caffolding **T**ool) is a tool to improve draft assembly by sequencing data of progeny or closely related strains. It exploits genetic information in progeny data to achieve better contiguity and correctness.

## Installation

CAST only needs Java(1.8 or above) and BLAST to run. It currently requires BLAST in system PATH so that CAST can call it directly.

Please go to [release](https://github.com/taoistly/CAST/releases) to download CAST.jar.

For easy configuration (e.g. on windows), CAST is also available on docker. Please check **Run on Docker** section below for step-by-step guide.

## Input

1. A parent draft assembly(fasta) from any platform/assembler
1. Over 10 progenies reads(fastq) with sequencing depth > 5

## Usage
First, please use [BWA](http://bio-bwa.sourceforge.net/) to map progenies reads to draft assembly and generate result in sorted bam format. (For details, please check **FAQ "How to prepare bam files and their indexes?"** below.)

Then:
```shell
java -jar CAST.jar <draft assembly> <list of bam files>
```
For example:
```shell
java -jar CAST.jar myassembly.fasta progeny1.bam progeny2.bam progeny3.bam
```
To specify many bam files, it is better to put all your bam files in a folder and then use wildcard like:
```shell
java -jar CAST.jar myassembly.fasta bamfiles/*.bam
```

## Output

* `report` contains candidate links inferred by genetic information. (intermediate output for check)
* `CAST.fasta` is the improved assembly. (final output)
* `heatmap` folder contains heatmaps for both ends of all contigs (by `--heatmap` option)

## FAQ

#### Q: What kind of data is appropriate?
Although CAST is designed to use progeny data, it can also improve a draft assembly by siblings or parents data. The sequencing platform for progenies is required to have a low sequencing error rate, e.g. Illumina.


#### Q: How to prepare bam files and their index?
If you have fastq file from progeny sequencing, you need to:
1. install BWA and samtools
2. align all progenies reads to your draft assembly
3. index generated bam files

For example:
```shell
bwa index draft.fasta
bwa mem draft.fasta p1_1.fastq p1_2.fastq | samtools sort - p1.bam
samtools index p1.bam
```

#### Q: How long does CAST take to run?

It basically depends on the amount of data as well as your disk access performance. You will see a progress bar, the elapsed time, and the estimated time remaining. A very rough estimate: 2 hours for 100G bam files. 

## Run on docker

1. Download and install docker. It's available for most platforms. (need root privilege)
2. Run docker to install CAST: `$docker pull taoistly/cast`
3. Place bam and fasta file in a folder
4. Run CAST: `docker run -tv /path/to/the/folder:/data taoistly/cast ` 

** Because of a technical limitation of docker, for Windows user, /path/to/files should be a folder under C:\users\yourUsername\ and type the path in the way like /c/users/yourUsername/somefolder **

## Citation
To be published

## Feedback
If you encounter any problem or have feedback, please feel free to contact me (x86@u.nus.edu).
I will try to reply ASAP.