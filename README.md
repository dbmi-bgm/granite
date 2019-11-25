# novoCaller 2

## About
novoCaller 2 is a Bayesian de novo variant calling algorithm that uses information from read-level data both in the pedigree and in unrelated samples.
The program represents an updated and improved implementation of the original algorithm described in [paper](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty749/5087716).

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
    ./novoCaller_2.py [-h] -i INPUTFILE -o OUTPUTFILE -u UNRELATEDBAMS -t TRIOBAMS [-a ALLELEFREQTHR]

    # !!! if the above is not working run
    python novoCaller_2.py [-h] -i INPUTFILE -o OUTPUTFILE -u UNRELATEDBAMS -t TRIOBAMS [-a ALLELEFREQTHR]

## Input
The program accepts files in vcf format (VCFv4.x).

## Arguments
  - **-i**, **--inputfile**, path to input file in vcf format
  - **-o**, **--outputfile**, path to output file to be used
  - **-u**, **--unrelatedbams**, path to tsv file containing `ID<TAB>Path/to/file` for unrelated bam files used to train the model
  - **-t**, **--triobams**, path to tsv file containing `ID<TAB>Path/to/file` for family bam files, the PROBAND must be listed as LAST
  - [**-a**, **--allelefreqthr**], threshold used to filter by population allele frequency. If specified, the allele frequency to be considered for each variant must be provided as annotation in the INFO field using the tag `novoAF=<float>;`.
