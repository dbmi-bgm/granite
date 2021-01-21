# granite

[![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://www.python.org/downloads/release/python-360/)

granite is a collection of software to work with genomic variants. The suite provides inheritance mode callers and utilities to filter and refine variants called by other methods in VCF format.

For more details, see granite [*documentation*](https://granite-suite.readthedocs.io/en/latest/ "granite documentation").

## Availability and requirements
A ready-to-use docker image is available to download.

    docker pull b3rse/granite:v0.1.7

To run locally, install the following libraries:

    pip install numpy pysam bitarray h5py matplotlib
    pip install --user pytabix

Additional software needs to be available in the environment:

  - [*samtools*](http://www.htslib.org/ "samtools documentation")
  - [*bgzip*](http://www.htslib.org/doc/bgzip.1.html "bgzip documentation")
  - [*tabix*](http://www.htslib.org/doc/tabix.1.html "tabix documentation")

To install the program from source, run the following commands:

    git clone https://github.com/dbmi-bgm/granite
    cd granite
    python setup.py install

To install the program with pip:

    pip install granite-suite

## File formats
The program is compatible with standard BED, BAM and VCF formats (`VCFv4.x`).

### ReadCountKeeper (.rck)
RCK is a tabular format that allows to efficiently store counts by strand (ForWard-ReVerse) for reads that support REFerence allele, ALTernate alleles, INSertions or DELetions at CHRomosome and POSition. RCK files can be further compressed with *bgzip* and indexed with *tabix* for storage, portability and faster random access. 1-based.

Tabular format structure:

    #CHR   POS   COVERAGE   REF_FW   REF_RV   ALT_FW   ALT_RV   INS_FW   INS_RV   DEL_FW   DEL_REV
    13     1     23         0        0        11       12       0        0        0        0
    13     2     35         18       15       1        1        0        0        0        0

Commands to compress and index files:

    bgzip PATH/TO/FILE
    tabix -b 2 -s 1 -e 0 -c "#" PATH/TO/FILE.gz

### BinaryIndexGenome (.big)
BIG is a hdf5-based binary format that stores boolean values for each genomic position as bit arrays. Each position is represented in three complementary arrays that account for SNVs (Single-Nucleotide Variants), insertions and deletions respectively. 1-based.

hdf5 format structure:

    e.g.
    chr1_snv: array(bool)
    chr1_ins: array(bool)
    chr1_del: array(bool)
    chr2_snv: array(bool)
    ...
    ...
    chrM_del: array(bool)

*note*: hdf5 keys are built as the chromosome name based on reference (e.g. chr1) plus the suffix specifying whether the array represents SNVs (_snv), insertions (_ins) or deletions (_del).

### Pedigree in JSON format
When the program requires pedigree information, the expected format is as follow:

    [
      {
        "individual": "NA12877",
        "sample_name": "NA12877_sample",
        "gender": "M",
        "parents": []
      },
      {
        "individual": "NA12878",
        "sample_name": "NA12878_sample",
        "gender": "F",
        "parents": []
      },
      {
        "individual": "NA12879",
        "sample_name": "NA12879_sample",
        "gender": "F",
        "parents": ["NA12878", "NA12877"]
      }
    ]

where `individual` is the unique identifier for member inside the pedigree, `sample_name` is the corresponding sample ID in VCF file, and `parents` is the list of unique identifiers for member parents if any.

## Usage
```text
    granite <command> ...

    positional arguments:
      <command>
        novoCaller   Bayesian de novo variant caller
        comHet       compound heterozygous variant caller
        mpileupCounts
                     samtools wrapper to calculate reads statistics for pileup at
                     each position
        blackList    utility to blacklist and filter out variants from input VCF
                     file based on positions set in BIG format file and/or
                     population allele frequency
        whiteList    utility to whitelist and select a subset of variants from
                     input VCF file based on specified annotations and positions
        cleanVCF     utility to clean INFO field of input VCF file
        geneList     utility to clean VEP annotations of input VCF file using a
                     list of genes
        toBig        utility that converts counts from bgzip and tabix indexed RCK
                     format into BIG format. Positions are "called" by reads
                     counts or allelic balance for single or multiple files (joint
                     calls) in specified regions
        rckTar       utility to create a tar archive from bgzip and tabix indexed
                     RCK files. Creates an index file for the archive
        qcVCF        utility to create a report of different metrics calculated
                     for input VCF file
        validateVCF  utility to calculate error models for input VCF file using
                     pedigree information
```

### novoCaller
novoCaller is a Bayesian calling algorithm for *de novo* mutations. The model uses read-level information both in pedigree (trio) and unrelated samples to rank and assign a probabilty to each call. The software represents an updated and improved implementation of the original algorithm described in [Mohanty et al. 2019](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty749/5087716).

### comHet
comHet is a calling algorithm for *compound heterozygous* mutations. The model uses genotype-level information in pedigree (trio) and VEP-based annotations to call possible compound heterozygous pairs. VEP annotations are used to assign variants to genes and transcripts, genotype information allows to refine calls based on inheritance mode. Calls are further flagged as "Phased" or "Unphased", where "Phased" means that genotype information supports in-trans inheritance for alternate alleles from parents.

### blackList
blackList allows to filter-out variants from input VCF file based on positions set in BIG format file and/or provided population allele frequency.

### whiteList
whiteList allows to select and filter-in a subset of variants from input VCF file based on specified annotations and positions. The software can use provided VEP, ClinVar or SpliceAI annotations. Positions can be also specfied as a BED format file.

### cleanVCF
cleanVCF allows to clean INFO field of input VCF file. The software can remove a list of TAG from INFO field, or can be used to clean VEP annotations.

### geneList
geneList allows to clean VEP annotations by applyng a list of genes. The software removes all the transcripts that do not map to a gene on the list.

### qcVCF
qcVCF produces a report in JSON format with different quality metrics calculated for input VCF file. Both single sample and family-based metrics are available.

### mpileupCounts
mpileupCounts uses *samtools* to access input BAM and calculates statistics for reads pileup at each position in the specified region, returns counts in RCK format.

### toBig
toBig converts counts from bgzip and tabix indexed RCK format into BIG format. Positions are "called" by read counts or allelic balance for single or multiple files (joint calls) in specified regions. Positions "called" are set to True (or 1) in BIG binary structure.

### rckTar
rckTar creates a tar archive from bgzip and tabix indexed RCK files. Creates an index file for the archive.

### validateVCF
validateVCF allows to calculate error models for different inheritance modes for input VCF file using pedigree information.
