## Utilities

### mpileupCounts
mpileupCounts uses *Samtools* to access input BAM and calculates statistics for read pileups at each position in the specified region, returns counts in RCK format.

#### Arguments
```text
    usage: granite mpileupCounts [-h] -i INPUTFILE -o OUTPUTFILE -r REFERENCE
                                 [--region REGION] [--MQthr MQTHR] [--BQthr BQTHR]

    optional arguments:
      -i INPUTFILE, --inputfile INPUTFILE
                            input file in BAM format
      -o OUTPUTFILE, --outputfile OUTPUTFILE
                            output file to write results as RCK format (TSV), use
                            .rck as extension
      -r REFERENCE, --reference REFERENCE
                            reference file in FASTA format
      --region REGION       region to be analyzed [e.g. chr1:1-10000000,
                            1:1-10000000, chr1, 1], chromosome name must match the
                            reference
      --MQthr MQTHR         minimum mapping quality for an alignment to be used
                            (>=) [0]
      --BQthr BQTHR         minimum base quality for a base to be considered (>=)
                            [13]
```

### toBig
toBig converts counts from bgzip and tabix indexed RCK format into BIG format. Positions are "called" by read counts or allelic balance for single or multiple files (joint calls) in specified regions. Positions "called" are set to True (or 1) in BIG binary structure.

#### Arguments
```text
    usage: granite toBig [-h] -f FILE [-f FILE ...] -o OUTPUTFILE -r
                         REGIONFILE -c CHROMFILE [--ncores NCORES] --fithr FITHR
                         [--rdthr RDTHR] [--abthr ABTHR]

    optional arguments:
      -f FILE, --file FILE  file to be used to call positions. To do joint calling
                            specify multiple files as: "-f file_1 -f file_2 -f ...".
                            Expected bgzip and tabix indexed RCK file
      -o OUTPUTFILE, --outputfile OUTPUTFILE
                            output file to write results as BIG format (binary
                            hdf5), use .big as extension
      -r REGIONFILE, --regionfile REGIONFILE
                            file containing regions to be used [e.g.
                            chr1:1-10000000, 1:1-10000000, chr1, 1] listed as a
                            column, chromosomes names must match the reference
      -c CHROMFILE, --chromfile CHROMFILE
                            chrom.sizes file containing chromosomes size
                            information
      --ncores NCORES       number of cores to be used if multiple regions are
                            specified [1]
      --fithr FITHR         minimum number of files with at least "--rdthr" for
                            the alternate allele or having the variant, "calls" by
                            allelic balance, to jointly "call" position (>=)
      --rdthr RDTHR         minimum number of alternate reads to count the file in
                            "--fithr", if not specified "calls" are made by
                            allelic balance (>=)
      --abthr ABTHR         minimum percentage of alternate reads compared to
                            reference reads to count the file in "--fithr" when
                            "calling" by allelic balance (>=) [15]
```

#### Examples
toBig can be used to calculate positions to blacklist for common variants by using unrelated samples. This command will set to `True` in BIG structure positions with allelic balance for alternate allele equal/higher than `--abthr` in more than `--fithr` samples (joint calling).

    granite toBig -f file -f file -f file -f file -f ... -o file.out.big -c file.chrom.sizes -r file.regions --fithr <int> --abthr <int>

Absolute read counts can be used instead of allelic balance to call positions. This command will set to `True` in BIG structure positions with read counts for alternate allele equal/higher than `--rdthr` in more than `--fithr` samples (joint calling).

    granite toBig -f file -f file -f file -f file -f ... -o file.out.big -c file.chrom.sizes -r file.regions --fithr <int> --rdthr <int>

### rckTar
rckTar creates a tar archive from bgzip and tabix indexed RCK files. Creates an index file for the archive.

#### Arguments
```text
    usage: granite rckTar [-h] -t TTAR -f FILE

    optional arguments:
      -t TTAR, --ttar TTAR  target tar to write results, use .rck.tar as extension
      -f FILE, --file FILE  file to be archived. Specify multiple files as: "-f
                            SampleID_1.rck.gz -f SampleID_2.rck.gz -f ...". Files
                            order is maintained while creating the index
```
