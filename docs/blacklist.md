### blackList
blackList allows to filter-out variants from input VCF file based on positions set in BIG format file and/or provided population allele frequency.

#### Arguments
    usage: granite blackList [-h] -i INPUTFILE -o OUTPUTFILE [-b BIGFILE]
                             [--aftag AFTAG] [--afthr AFTHR]

    optional arguments:
      -i INPUTFILE, --inputfile INPUTFILE
                            input VCF file
      -o OUTPUTFILE, --outputfile OUTPUTFILE
                            output file to write results as VCF, use .vcf as
                            extension
      -b BIGFILE, --bigfile BIGFILE
                            BIG format file with positions set for blacklist
      --aftag AFTAG         TAG (TAG=<float>) or TAG field to be used to filter by
                            population allele frequency
      --afthr AFTHR         threshold to filter by population allele frequency
                            (<=) [1]

#### Examples
Blacklist variants based on position set to `True` in BIG format file.

    granite blackList -i file.vcf -o file.out.vcf -b file.big

Blacklist variants based on population allele frequency. This filters out variants with allele frequency higher than `--afthr`. Allele frequency must be provided for each variant in INFO column.

    granite blackList -i file.vcf -o file.out.vcf --afthr <float> --aftag tag

Combine the two filters.

    granite blackList -i file.vcf -o file.out.vcf --afthr <float> --aftag tag -b file.big
