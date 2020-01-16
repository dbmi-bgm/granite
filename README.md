# novoCaller 2


## About
novoCaller 2 is a Bayesian variant calling algorithm for *de novo* mutations that uses information from read-level data both in pedigree and unrelated samples.
The program represents an updated and improved implementation of the original algorithm described in [paper](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty749/5087716).

novoCaller 2 has been extended with new functionalities:

  - the program can now use a compressed tabular format () that encodes read-level information in a more accessible and portable way as compared with bam
  - the program can *blacklist* and filter-out commonly shared variants and sequencing artifacts using read-level data in unrelated samples

novoCaller 2 is accompanied by additional scripts to format data and create bitarrays used to fast blacklist commonly shared positions, *mpileup_parser* and *to_bitarray*.


## Requirements
novoCaller 2 and accompanying scripts are available as a docker image ready-to-use. Dockerfile is available inside docker folder.

To run locally novocaller 2 requires Python (version 3.X), *numpy*, *pysam* and *bitarray* libraries.

To install [*numpy*](https://docs.scipy.org/doc/ "numpy documentation") under unix environment (linux, osx):

    pip install numpy

To install [*pysam*](https://pysam.readthedocs.io/en/latest/ "pysam documentation") under unix environment (linux, osx):

    pip install pysam

To install [*bitarray*](https://pypi.org/project/bitarray/ "bitarray documentation") under unix environment (linux, osx):

    pip install bitarray

to_bitarray additionally requires *pytabix* and *h5py* libraries.

To install [*pytabix*](https://pypi.org/project/pytabix/ "pytabix documentation") under unix environment (linux, osx):

    pip install --user pytabix

To install [*h5py*](https://www.h5py.org/ "h5py documentation") under unix environment (linux, osx):

    pip install h5py


## novoCaller_2
To run novoCaller_2 locally under unix environment (linux, osx):

    # move to novoCaller_2 folder
    cd PATH/TO/NOVOCALLER_2/FOLDER/

    # make novoCaller_2 executable
    sudo chmod +x novoCaller_2.py

    # run novoCaller_2
    ./novoCaller_2.py [-h] -i INPUTFILE -o OUTPUTFILE [-u UNRELATEDBAMS]
                       [-t TRIOBAMS] [-p POSTPROBTHR] [-b BLACKLIST]
                       [--thr_bams THR_BAMS] [--thr_reads THR_READS]
                       [-a ALLELEFREQTHR]

    # !!! if the above is not working run
    python novoCaller_2.py [-h] -i INPUTFILE -o OUTPUTFILE [-u UNRELATEDBAMS]
                       [-t TRIOBAMS] [-p POSTPROBTHR] [-b BLACKLIST]
                       [--thr_bams THR_BAMS] [--thr_reads THR_READS]
                       [-a ALLELEFREQTHR]

### Input
The program accepts files in VCF format (VCFv4.x).

*De novo* analysis, input files must contain genotype information for the trio in addition to standard VCF columns. Columns IDs for trio must match the IDs provided together with the list of bam files (TRIOBAMS).

    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG002   HG003   HG004   ...

*Blacklist* analysis, input files must contain genotype information for the proband but do not require additional information for trio or family.

### Arguments
  - **-i**, **--inputfile**, *I/O*: path to input file in vcf format
  - **-o**, **--outputfile**, *I/O*: path to output file to write results, vcf format
  - [**-u**, **--unrelatedbams**], *DE NOVO*: path to tsv file listing `ID<TAB>Path/to/file` for unrelated bam files used to train the model
  - [**-t**, **--triobams**], *DE NOVO*: path to tsv file listing `ID<TAB>Path/to/file` for family bam files, **PROBAND must be listed as LAST**
  - [**-p**, **--postprobthr**], *DE NOVO*: threshold to filter by posterior probabilty for de novo calls [0]
  - [**-b**, **--blacklist**], *BLACKLIST*: path to tsv file listing `ID<TAB>Path/to/file` for unrelated bam files used to filter out shared variants/artifacts
  - [**--thr_reads**], *BLACKLIST*: minimum number of reads to count the bam file in "--thr_bams" [1]
  - [**--thr_bams**], *BLACKLIST*: minimum number of bam files with at least "--thr_reads" to blacklist the variant [2]
  - [**-a**, **--allelefreqthr**], threshold to filter by population allele frequency. If specified, the allele frequency for each variant must be provided as annotation in the INFO field using the tag `novoAF=<float>;` [1]

### Output
The program generates output in VCF format.


## mpileup_parser
mpileup_parser uses [*samtools*](http://www.htslib.org/ "samtools documentation") to calculate statistics for pileup at each position in specified region and returns counts in TSV format. Counts by strand (ForWard or ReVerse) are provided for reads that support REFerence allele, ALTernate alleles, INSertions or DELetions at POSition. The file can be compressed with bgzip and tabix indexed for storage, portability and speeding-up random access.

TSV format structure:

    #CHR    POS    COVERAGE    REF_FW    REF_RV    ALT_FW    ALT_RV    INS_FW    INS_RV    DEL_FW    DEL_REV

To compress and index the file:

    bgzip PATH/TO/FILE
    tabix -b 2 -s 1 -e 0 -c "#" PATH/TO/FILE.gz

mpileup_parser-parallel script based on [*parallel*](https://www.gnu.org/software/parallel/ 'parallel documentation') is available in docker folder to speed up computation by splitting the bam and parallelizing by regions.


## to_bitarray
to_bitarray calls positions by reads counts or allelic balance (TO DO) for single bam or multiple bams (joint calls) in the specified region. Results are stored in a binary format where bits corresponding to called positions are set to 1.
