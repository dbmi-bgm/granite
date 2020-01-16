# novoCaller 2

## About
novoCaller 2 is a Bayesian de novo variant calling algorithm that uses information from read-level data both in the pedigree and in unrelated samples.
The program represents an updated and improved implementation of the original algorithm described in [paper](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty749/5087716).

The program has also beeing extended with new functionalities:

  - Variants blacklist: read-level data in unrelated samples are used to filter out common shared variants and sequencing artifacts.
  - ..

## Requirements
novoCaller requires Python (version 3.X), *numpy* and *pysam* libraries.

To install [*numpy*](https://docs.scipy.org/doc/ "numpy documentation") under unix environment (linux, osx):

    sudo pip install numpy

To install [*pysam*](https://pysam.readthedocs.io/en/latest/ "pysam documentation") under unix environment (linux, osx):

    sudo pip install pysam

## Running novoCaller 2
To run novoCaller 2 under unix environment (linux, osx):

    # move to QPARSE folder
    cd PATH/TO/NOVOCALLER_2/FOLDER/

    # make QPARSE executable
    sudo chmod +x novoCaller_2.py

    # run QPARSE
    ./novoCaller_2.py [-h] -i INPUTFILE -o OUTPUTFILE [-u UNRELATEDBAMS]
                       [-t TRIOBAMS] [-p POSTPROBTHR] [-b BLACKLIST]
                       [--thr_bams THR_BAMS] [--thr_reads THR_READS]
                       [-a ALLELEFREQTHR]

    # !!! if the above is not working run
    python novoCaller_2.py [-h] -i INPUTFILE -o OUTPUTFILE [-u UNRELATEDBAMS]
                       [-t TRIOBAMS] [-p POSTPROBTHR] [-b BLACKLIST]
                       [--thr_bams THR_BAMS] [--thr_reads THR_READS]
                       [-a ALLELEFREQTHR]

## Input
The program accepts files in vcf format (VCFv4.x).

*De novo* analysis, input files must contain genotype information for the trio in addition to standard vcf columns. Columns IDs for trio must match the IDs provided together with the list of bam files (TRIOBAMS).

    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG002   HG003   HG004   ...

*Blacklist* analysis, input files must contain genotype information for the proband but do not require additional information for trio or family.

## Arguments
  - **-i**, **--inputfile**, *I/O*: path to input file in vcf format
  - **-o**, **--outputfile**, *I/O*: path to output file to write results, vcf format
  - [**-u**, **--unrelatedbams**], *DE NOVO*: path to tsv file listing `ID<TAB>Path/to/file` for unrelated bam files used to train the model
  - [**-t**, **--triobams**], *DE NOVO*: path to tsv file listing `ID<TAB>Path/to/file` for family bam files, **PROBAND must be listed as LAST**
  - [**-p**, **--postprobthr**], *DE NOVO*: threshold to filter by posterior probabilty for de novo calls [0]
  - [**-b**, **--blacklist**], *BLACKLIST*: path to tsv file listing `ID<TAB>Path/to/file` for unrelated bam files used to filter out shared variants/artifacts
  - [**--thr_reads**], *BLACKLIST*: minimum number of reads to count the bam file in "--thr_bams" [1]
  - [**--thr_bams**], *BLACKLIST*: minimum number of bam files with at least "--thr_reads" to blacklist the variant [2]
  - [**-a**, **--allelefreqthr**], threshold to filter by population allele frequency. If specified, the allele frequency for each variant must be provided as annotation in the INFO field using the tag `novoAF=<float>;` [1]

## Output
The program generates output in vcf format.
