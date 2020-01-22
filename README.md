# Granite (in development)
Granite (Genomic vaRiANts calling and fIltering uTilitiEs) is a collection of software to call, filter and work with genomic variants.


## Requirements
A ready-to-use docker image is available to download.

    docker pull IMAGE

To run locally Python (3.x) is required together with the following libraries:

  - [*numpy*](https://docs.scipy.org/doc/ "numpy documentation")
  - [*pysam*](https://pysam.readthedocs.io/en/latest/ "pysam documentation")
  - [*bitarray*](https://pypi.org/project/bitarray/ "bitarray documentation")
  - [*pytabix*](https://pypi.org/project/pytabix/ "pytabix documentation")
  - [*h5py*](https://www.h5py.org/ "h5py documentation")

To install with pip:

    pip install numpy pysam bitarray h5py
    pip install --user pytabix

Additional software needs to be available in the environment:

  - [*samtools*](http://www.htslib.org/ "samtools documentation")
  - [*bgzip*](http://www.htslib.org/doc/bgzip.1.html "bgzip documentation")
  - [*tabix*](http://www.htslib.org/doc/tabix.1.html "tabix documentation")


# Callers

## novoCaller
novoCaller is a Bayesian variant calling algorithm for *de novo* mutations that uses information from read-level data both in pedigree and unrelated samples.
The software represents an updated and improved implementation of the original algorithm described in [paper](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty749/5087716).

### Arguments
  - **-i**, **--inputfile**, path to input file in VCF format
  - **-o**, **--outputfile**, path to output file to write results in VCF format
  - **-u**, **--unrelated**, path to TSV file listing `ID<TAB>Path/to/file` for unrelated files used to train the model
  - **-t**, **--trio**, path to TSV file listing `ID<TAB>Path/to/file` for trio files (pedigree), **PROBAND must be listed as LAST**
  - [**--bam**], by default the program expect family and unrelated files in XX format, if you want to use bam files instead use this flag
  - [**-p**, **--postprobthr**], threshold to filter calls by posterior probabilty [0]
  - [**-a**, **--allelefreqthr**], threshold to filter calls by population allele frequency. If specified, the allele frequency for each variant must be provided as annotation in the INFO field using the tag `novoAF=<float>;` [1]

### Input
The program accepts files in VCF format (VCFv4.x) as input. Files must contain genotype information for trio in addition to standard VCF columns. Columns IDs for trio must match the IDs provided together with the list of bam files (*--triobams*).

Required VCF format structure:

    #CHROM   POS   ID   REF   ALT   QUAL   FILTER   INFO   FORMAT   PROBAND_ID   MOTHER_ID   FATHER_ID   ...

### Output
The program generates output in VCF format.


# Utilities

## blackList
blackList allows to filter-out commonly shared variants and sequencing artifacts using read-level data or blacklisted positions in unrelated samples.

### Arguments
  - **-i**, **--inputfile**, path to input file in VCF format
  - **-o**, **--outputfile**, path to output file to write results in VCF format
  - **-u**, **--unrelated**,
  - **-b**, **--blacklist**,
  - [**-a**, **--allelefreqthr**], threshold to filter calls by population allele frequency

### Input
The program accepts files in VCF format (VCFv4.x) as input. Files must contain genotype information for the proband but do not require additional information for trio or family.

Required VCF format structure:

    #CHROM   POS   ID   REF   ALT   QUAL   FILTER   INFO   FORMAT   PROBAND_ID   ...


## mpileupCounts
mpileupCounts uses *samtools* to calculate statistics for reads pileup at each position in the specified region and returns counts in TSV format. Counts by strand (ForWard or ReVerse) are provided for reads that support REFerence allele, ALTernate alleles, INSertions or DELetions at POSition. The file can be further compressed with *bgzip* and indexed with *tabix* for storage, portability and speeding-up random access.

TSV format structure:

    #CHR   POS   COVERAGE   REF_FW   REF_RV   ALT_FW   ALT_RV   INS_FW   INS_RV   DEL_FW   DEL_REV

To compress and index the file:

    bgzip PATH/TO/FILE
    tabix -b 2 -s 1 -e 0 -c "#" PATH/TO/FILE.gz


## to_bitarray
to_bitarray calls positions by reads counts or allelic balance (TO DO) for single bam or multiple bams (joint calls) in specified region. Results are stored in a binary format where bits corresponding to called positions are set to 1.
