## Quality Check

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

### validateVCF
validateVCF produces a report in JSON format with error models for different inheritance modes calculated for input VCF file. Additional supporting plots are generated in PNG format.

#### Arguments
```text
    usage: granite validateVCF [-h] -i INPUTFILE -o OUTPUTFILE -p PEDIGREE
                               [PEDIGREE ...] --anchor ANCHOR [ANCHOR ...]
                               [--het HET [HET ...]] [--novo NOVO] [--verbose]

    optional arguments:
      -i INPUTFILE, --inputfile INPUTFILE
                            input VCF file
      -o OUTPUTFILE, --outputfile OUTPUTFILE
                            output file to write results as JSON, use .json as
                            extension
      -p PEDIGREE [PEDIGREE ...], --pedigree PEDIGREE [PEDIGREE ...]
                            pedigree information, either as JSON file or JSON
                            representation as string. It is possible to specify
                            multiple pedigrees to load as list (e.g. --pedigree
                            pedigree_1 [pedigree_2] ...)
      --anchor ANCHOR [ANCHOR ...]
                            sample ID to be used as anchor in pedigree to build
                            family. It is possible to specify multiple sample IDs
                            as list. If multiple pedigrees are specified in "--
                            pedigree", anchors are positionally matched to
                            corresponding pedigree
      --het HET [HET ...]   sample ID to be used to calculate error model for
                            heterozygous mutations. It is possible to specify
                            multiple sample IDs as list. Each sample ID must
                            correspond to anchor specified in "--anchor"
      --novo NOVO           sample ID to be used to calculate error model for de
                            novo mutations. Must correspond to anchor specified in
                            "--anchor". Requires posterior probability from
                            novoCaller
      --type TYPE [TYPE ...]
                            by default error models are calculated only for SNV.
                            It is possible to specify different types of variant
                            to use (SNV, INS, DEL, MNV, MAV) as list (e.g. --type
                            SNV INS DEL)
```

### SVqcVCF
SVqcVCF produces a report in JSON format with different quality metrics calculated for input SV VCF file. Currently, this function can count DEL and DUP SVs in single- and multi-sample SV VCF files.  It reports the number of DEL, DUP, and total (the sum of only DEL and DUP) SVs for each sample provided in samples. Other SVTYPEs (INS, INV, CNV, BND) are currently ignored.

#### Arguments
```text
    usage: granite SVqcVCF [-h] -i INPUTFILE -o OUTPUTFILE
                         --samples SAMPLES [SAMPLES ...] [--verbose]

    optional arguments:
      -i INPUTFILE, --inputfile INPUTFILE
                            input SV VCF file
      -o OUTPUTFILE, --outputfile OUTPUTFILE
                            output file to write results as JSON, use .json as
                            extension
      --samples SAMPLES [SAMPLES ...]
                            list of sample IDs to get stats for (e.g. --samples
                            SampleID_1 [SampleID_2] ...)
```
