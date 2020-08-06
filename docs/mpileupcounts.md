### mpileupCounts
mpileupCounts uses *samtools* to access input BAM and calculates statistics for reads pileup at each position in the specified region, returns counts in RCK format.

#### Arguments
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
