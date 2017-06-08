## RNA Seq data analysis pipeline
This pipeline runs pre processing, splice junction alignment and gene expression quantification on paired end illumina RNASeq data.

### Requirements
The package has all tools for its functioning included within the git repo. Some prerequisites are necessary for running the pipeline.

1. Python3 : The pipeline is written with Python3 syntax, requires at least Python3.4 or above.
2. Bowtie2 : Tophat2 requires the path to bowtie2 to be added to the system path.
3. Illumina iGenomes reference : The Illumina iGenomes reference for the organism needs to be downloaded, specifically the pipeline needs the Bowtie2Index and genes.gtf file for the organism

### Running the pipeline
The program has a help included

```{sh}
python3 rnaSeq.py -h
usage: Nestv2 [-h] [-i INP_PATH] -r REF_PATH -g GTF_PATH -o OUT_PATH
              [-a ADP_PATH] [-t TPH_PATH] [-c CUF_PATH] [-b BBD_PATH]

optional arguments:
  -h, --help            show this help message and exit
  -i INP_PATH, --inp_path INP_PATH
                        Path to input directory
  -r REF_PATH, --ref REF_PATH
                        Path to Reference bowtie index
  -g GTF_PATH, --gtf GTF_PATH
                        Path to Reference GTF file
  -o OUT_PATH, --outpath OUT_PATH
                        Path where all outputs will be stored
  -a ADP_PATH, --adapter ADP_PATH
                        Path to Adpater fasta file
  -t TPH_PATH, --tophat TPH_PATH
                        Path to Tophat2 executable
  -c CUF_PATH, --cufflinks CUF_PATH
                        Path to Cufflinks executable
  -b BBD_PATH, --bbduk BBD_PATH
                        Path to BBDuk executable
```
#### Note :

 - To run the pipeline specify the path where the sample fastq files are stored. The pipeline has the ability to pick up the fastq files from the directory. Eg: If fastq files for the samples are stored in fq/sample1/*.fq then just provide fq/ to the -i argument

 - This enables the analysis of multiple samples if fq/ folder

 - The reference path should only contain the prefix of the bowtie index for eg: if the bowtie index is in the path ref/genome.1.bt2 then just provide ref/genome to the -r argument.
