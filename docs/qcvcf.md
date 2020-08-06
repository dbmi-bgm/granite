### qcVCF
qcVCF produces a report in JSON format with different quality metrics calculated for input VCF file. Both single sample and family-based metrics are available.

#### Arguments
```text
    usage: granite qcVCF [-h] -i INPUTFILE -o OUTPUTFILE -p PEDIGREE --samples
                         SAMPLES [SAMPLES ...] [--ti_tv] [--trio_errors]
                         [--het_hom]

    optional arguments:
      -i INPUTFILE, --inputfile INPUTFILE
                            input VCF file
      -o OUTPUTFILE, --outputfile OUTPUTFILE
                            output file to write results as JSON, use .json as
                            extension
      -p PEDIGREE, --pedigree PEDIGREE
                            pedigree information, either as JSON file or JSON
                            representation as string
      --samples SAMPLES [SAMPLES ...]
                            list of sample IDs to get stats for (e.g. --samples
                            SampleID_1 [SampleID_2] ...)
      --ti_tv               add transition-transversion ratio and statistics on
                            substitutions to report
      --trio_errors         add statistics on mendelian errors based on trio to
                            report
      --het_hom             add heterozygosity ratio and statistics on zygosity to
                            report
```
