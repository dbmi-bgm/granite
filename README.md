# granite

[![Python 3.9+](https://img.shields.io/badge/Python-3.9+-blue.svg)](https://www.python.org) ![License MIT](https://img.shields.io/badge/License-MIT-green)

granite is a software suite for working with genomic variants. The suite provides inheritance mode callers and utilities to filter and refine variants called by other methods in VCF format.

The granite library can also be used through a Python API to manipulate files in VCF format.

For more details, see granite [*documentation*](https://granite-suite.readthedocs.io/en/latest/ "granite documentation").

## Availability and requirements
A ready-to-use Docker image is available to download [here](https://hub.docker.com/r/b3rse/granite "granite Docker").

    docker pull b3rse/granite:<version>

Run granite using Docker as follows:

    docker run b3rse/granite:<version> granite <command> ...

If installed locally, additional software needs to be available in the environment to run certain commands:

  - [*Samtools*](http://www.htslib.org/ "Samtools documentation")
  - [*bgzip*](http://www.htslib.org/doc/bgzip.1.html "bgzip documentation")
  - [*tabix*](http://www.htslib.org/doc/tabix.1.html "tabix documentation")

To install the program from source, run the following commands:

    git clone https://github.com/dbmi-bgm/granite
    cd granite
    make configure
    make build

To install the program with pip:

    pip install granite-suite

### VEP annotations
The comHet command requires gene- and transcript-level annotations produced by the Ensembl Variant Effect Predictor (VEP). granite does not run VEP directly; instead, it consumes VEP annotations already present in the input VCF file and uses them to assign variants to genes and transcripts and to evaluate inheritance patterns.

VEP annotations are expected to follow the standard VCF INFO field format (e.g. the *CSQ* tag), as produced by default VEP configurations. For details on running VEP and generating compatible annotations, please refer to the official documentation:
https://github.com/Ensembl/ensembl-vep

## File formats
The program is compatible with standard BED, BAM and VCF formats (`VCFv4.x`).

### ReadCountKeeper (.rck)
RCK is a tabular format that allows to efficiently store counts by strand (ForWard-ReVerse) for reads that support REFerence allele, ALTernate alleles, INSertions or DELetions at CHRomosome and POSition. RCK files can be further compressed with *bgzip* and indexed with *tabix* for storage, portability and faster random access. 1-based.

Tabular format structure:

    #CHR   POS   COVERAGE   REF_FW   REF_RV   ALT_FW   ALT_RV   INS_FW   INS_RV   DEL_FW   DEL_RV
    13     1     23         0        0        11       12       0        0        0        0
    13     2     35         18       15       1        1        0        0        0        0

Commands to compress and index files:

    bgzip PATH/TO/FILE
    tabix -b 2 -s 1 -e 0 -c "#" PATH/TO/FILE.gz

### BinaryIndexGenome (.big)
BIG is a HDF5-based binary format that stores boolean values for each genomic position as bit arrays. Each position is represented in three complementary arrays that account for SNVs (Single-Nucleotide Variants), insertions and deletions respectively. 1-based.

HDF5 format structure:

    e.g.
    chr1_snv: array(bool)
    chr1_ins: array(bool)
    chr1_del: array(bool)
    chr2_snv: array(bool)
    ...
    ...
    chrM_del: array(bool)

*note*: HDF5 keys are built as the chromosome name based on reference (e.g. chr1) plus the suffix specifying whether the array represents SNVs (_snv), insertions (_ins) or deletions (_del).

### Pedigree in JSON format
When the program requires pedigree information, the expected format is as follows:

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

where `individual` is the unique identifier for a member within the pedigree, `sample_name` is the corresponding sample ID in VCF file, and `parents` is the list of unique identifiers for the parents, if any.

## Usage
```text
    granite <command> ...

    positional arguments:
      <command>
        novoCaller   Bayesian de novo variant caller
        comHet       compound heterozygous variant caller
        filterByTag  utility to filter variants from input VCF file
                     by INFO field tags and thresholds
        mpileupCounts
                     samtools wrapper to calculate read statistics for pileup at
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
        SVqcVCF      utility to create a report of different metrics calculated
                     for input SV VCF file
        validateVCF  utility to calculate error models for input VCF file using
                     pedigree information
```

### novoCaller
novoCaller is a Bayesian calling algorithm for *de novo* mutations. The model uses read-level information both in pedigree (trio) and unrelated samples to rank and assign a probability to each call. The software represents an updated and improved implementation of the original algorithm described in [Mohanty et al. 2019](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty749/5087716).

*WARNING: starting from version 0.1.12, novoCaller `--triofiles` expected order changed. Now PROBAND must be listed as first.*

### comHet
comHet is a calling algorithm for *compound heterozygous* mutations. The model uses genotype-level information in pedigree (trio) and VEP-based annotations to call possible compound heterozygous pairs. VEP annotations are used to assign variants to genes and transcripts, genotype information allows to refine calls based on inheritance mode. Calls are further flagged as "Phased" or "Unphased", where "Phased" means that genotype information supports in-trans inheritance for alternate alleles from parents.

### filterByTag
filterByTag allows to filter variants from input VCF file by INFO field tags and thresholds. Users can define one or more filters on numeric, string, or boolean annotations. Each filter specifies the tag or field name, a value, an operator, type, and logic for aggregation. Multiple filters can be combined with global across-tag logic (any or all).

### blackList
blackList allows to filter-out variants from input VCF file based on positions set in BIG format file and/or provided population allele frequency. Positions can be also specified as a BED format file.

### whiteList
whiteList allows to select and filter-in a subset of variants from input VCF file based on specified annotations and positions. The software can use provided VEP, ClinVar or SpliceAI annotations. Positions can be also specified as a BED format file.

### cleanVCF
cleanVCF allows to clean INFO field of input VCF file. The software can remove a list of TAG from INFO field, or can be used to clean VEP annotations.

### geneList
geneList allows to clean VEP annotations by applying a list of genes. The software removes all the transcripts that do not map to a gene on the list.

### qcVCF
qcVCF produces a report in JSON format with different quality metrics calculated for input VCF file. Both single sample and family-based metrics are available.

### SVqcVCF
SVqcVCF produces a report in JSON format with different quality metrics calculated for input SV VCF file. Currently, this function can count DEL and DUP SVs in single- and multi-sample SV VCF files. It reports the number of DEL, DUP, and total (the sum of only DEL and DUP) SVs for each sample provided in samples. Other SVTYPEs (INS, INV, CNV, BND) are currently ignored.

### mpileupCounts
mpileupCounts uses *Samtools* to access input BAM and calculates statistics for read pileups at each position in the specified region, returns read counts in RCK format.

### toBig
toBig converts counts from bgzip and tabix indexed RCK format into BIG format. Positions are "called" by read counts or allelic balance for single or multiple files (joint calls) in specified regions. Positions "called" are set to True (or 1) in BIG binary structure.

### rckTar
rckTar creates a tar archive from bgzip and tabix indexed RCK files. Creates an index file for the archive.

### validateVCF
validateVCF allows to calculate error models for different inheritance modes for input VCF file using pedigree information.
