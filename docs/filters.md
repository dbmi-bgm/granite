## Variants Filtering

### whiteList
whiteList allows to select and filter-in a subset of variants from input VCF file based on specified annotations and positions. The software can use provided VEP, CLINVAR or SpliceAI annotations. Positions can be also specfied as a BED format file.

#### Arguments
```text
    usage: granite whiteList [-h] -i INPUTFILE -o OUTPUTFILE [--SpliceAI SPLICEAI]
                             [--SpliceAItag SPLICEAITAG] [--CLINVAR]
                             [--CLINVARonly CLINVARONLY [CLINVARONLY ...]]
                             [--CLINVARtag CLINVARTAG] [--VEP] [--VEPtag VEPTAG]
                             [--VEPrescue VEPRESCUE [VEPRESCUE ...]]
                             [--VEPremove VEPREMOVE [VEPREMOVE ...]]
                             [--VEPsep VEPSEP] [--BEDfile BEDFILE]

    optional arguments:
      -i INPUTFILE, --inputfile INPUTFILE
                            input VCF file
      -o OUTPUTFILE, --outputfile OUTPUTFILE
                            output file to write results as VCF, use .vcf as
                            extension
      --SpliceAI SPLICEAI   threshold to whitelist variants by SpliceAI value (>=)
      --SpliceAItag SPLICEAITAG
                            by default the program will search for "SpliceAI" TAG
                            (SpliceAI=<float>), use this parameter to specify a
                            different TAG | TAG field to be used (e.g. DS_DG)
      --CLINVAR             flag to whitelist all variants with a CLINVAR entry
                            [ALLELEID]
      --CLINVARonly CLINVARONLY [CLINVARONLY ...]
                            CLINVAR "CLNSIG" terms or keywords to be saved. Sets
                            for whitelist only CLINVAR variants with specified
                            terms or keywords
      --CLINVARtag CLINVARTAG
                            by default the program will search for CLINVAR
                            "ALLELEID" TAG, use this parameter to specify a
                            different TAG to be used
      --VEP                 use VEP "Consequence" annotations to whitelist exonic
                            and relevant variants (removed by default variants in
                            intronic, intergenic, or regulatory regions)
      --VEPtag VEPTAG       by default the program will search for "VEP" TAG
                            (VEP=<values>), use this parameter to specify a
                            different TAG to be used (e.g. CSQ)
      --VEPrescue VEPRESCUE [VEPRESCUE ...]
                            additional terms to overrule removed flags to rescue
                            and whitelist variants
      --VEPremove VEPREMOVE [VEPREMOVE ...]
                            additional terms to be removed
      --VEPsep VEPSEP       by default the program expects "&" as separator for
                            subfields in VEP (e.g.
                            intron_variant&splice_region_variant), use this
                            parameter to specify a different separator to be used
      --BEDfile BEDFILE     BED format file with positions to whitelist
```

#### Examples
Whitelists variants with CLINVAR entry. If available, CLINVAR annotation must be provided in INFO column.

    granite whiteList -i file.vcf -o file.out.vcf --CLINVAR

Whitelists only "Pathogenic" and "Likely_pathogenic" variants with CLINVAR entry. CLINVAR "CLNSIG" annotation must be provided in INFO column.

    granite whiteList -i file.vcf -o file.out.vcf --CLINVAR --CLINVARonly Pathogenic

Whitelists variants based on SpliceAI annotations. This filters in variants with SpliceAI score equal/higher than `--SpliceAI`. If available SpliceAI annotation must be provided in INFO column.

    granite whiteList -i file.vcf -o file.out.vcf --SpliceAI <float>

Whitelists variants based on VEP "Consequence" annotations. This withelists exonic and functional relevant variants by removing variants flagged as "intron_variant", "intergenic_variant", "downstream_gene_variant", "upstream_gene_variant", "regulatory_region_", "non_coding_transcript_". It is possible to specify additional [*terms*](https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html "VEP calculated consequences") to remove using `--VEPremove` and terms to rescue using `--VEPrescue`. To use VEP, annotation must be provided for each variant in INFO column.

    granite whiteList -i file.vcf -o file.out.vcf --VEP
    granite whiteList -i file.vcf -o file.out.vcf --VEP --VEPremove <str> <str>
    granite whiteList -i file.vcf -o file.out.vcf --VEP --VEPrescue <str> <str>
    granite whiteList -i file.vcf -o file.out.vcf --VEP --VEPrescue <str> <str> --VEPremove <str>

Whitelists variants based on positions specified as a BED format file.

    granite whiteList -i file.vcf -o file.out.vcf --BEDfile file.bed

Combine the above filters.

    granite whiteList -i file.vcf -o file.out.vcf --BEDfile file.bed --VEP --VEPrescue <str> <str> --CLINVAR --SpliceAI <float>

### blackList
blackList allows to filter-out variants from input VCF file based on positions set in BIG format file and/or provided population allele frequency.

#### Arguments
```text
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
```

#### Examples
Blacklist variants based on position set to `True` in BIG format file.

    granite blackList -i file.vcf -o file.out.vcf -b file.big

Blacklist variants based on population allele frequency. This filters out variants with allele frequency higher than `--afthr`. Allele frequency must be provided for each variant in INFO column.

    granite blackList -i file.vcf -o file.out.vcf --afthr <float> --aftag tag

Combine the two filters.

    granite blackList -i file.vcf -o file.out.vcf --afthr <float> --aftag tag -b file.big

### cleanVCF
cleanVCF allows to clean INFO field of input VCF file. The software can remove a list of TAG from INFO field, or can be used to clean VEP annotations.

#### Arguments
```text
    usage: granite cleanVCF [-h] -i INPUTFILE -o OUTPUTFILE [-t TAG] [--VEP]
                            [--VEPtag VEPTAG]
                            [--VEPrescue VEPRESCUE [VEPRESCUE ...]]
                            [--VEPremove VEPREMOVE [VEPREMOVE ...]]
                            [--VEPsep VEPSEP] [--SpliceAI SPLICEAI]
                            [--SpliceAItag SPLICEAITAG]

    optional arguments:
      -i INPUTFILE, --inputfile INPUTFILE
                            input VCF file
      -o OUTPUTFILE, --outputfile OUTPUTFILE
                            output file to write results as VCF, use .vcf as
                            extension
      -t TAG, --tag TAG     TAG to be removed from INFO field. Specify multiple
                            TAGs as: "-t TAG -t TAG -t ..."
      --VEP                 clean VEP "Consequence" annotations (removed by
                            default terms for intronic, intergenic, or regulatory
                            regions from annotations)
      --VEPtag VEPTAG       by default the program will search for "VEP" TAG
                            (VEP=<values>), use this parameter to specify a
                            different TAG to be used (e.g. CSQ)
      --VEPrescue VEPRESCUE [VEPRESCUE ...]
                            additional terms to overrule removed flags to rescue
                            annotations
      --VEPremove VEPREMOVE [VEPREMOVE ...]
                            additional terms to be removed from annotations
      --VEPsep VEPSEP       by default the program expects "&" as separator for
                            subfields in VEP (e.g.
                            intron_variant&splice_region_variant), use this
                            parameter to specify a different separator to be used
      --SpliceAI SPLICEAI   threshold to save intronic annotations, from VEP
                            "Consequence", for variants by SpliceAI value (>=)
      --SpliceAItag SPLICEAITAG
                            by default the program will search for "SpliceAI" TAG
                            (SpliceAI=<float>), use this parameter to specify a
                            different TAG | TAG field to be used (e.g. DS_DG)
```

#### Examples
Remove tag from INFO field.

    granite cleanVCF -i file.vcf -o file.out.vcf -t tag

Clean VEP based on VEP "Consequence" annotations. This remove annotations flagged as "intron_variant", "intergenic_variant", "downstream_gene_variant", "upstream_gene_variant", "regulatory_region_", "non_coding_transcript_". It is possible to specify additional [*terms*](https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html "VEP calculated consequences") to remove using `--VEPremove` and terms to rescue using `--VEPrescue`. VEP annotation must be provided for each variant in INFO column.

    granite cleanVCF -i file.vcf -o file.out.vcf --VEP
    granite cleanVCF -i file.vcf -o file.out.vcf --VEP --VEPremove <str> <str>
    granite cleanVCF -i file.vcf -o file.out.vcf --VEP --VEPrescue <str> <str>
    granite cleanVCF -i file.vcf -o file.out.vcf --VEP --VEPrescue <str> <str> --VEPremove <str>

The program also accepts a SpliceAI threshold that will rescue annotations for "intron_variant" by SpliceAI. SpliceAI annotation must be provided in INFO column.

    granite cleanVCF -i file.vcf -o file.out.vcf --VEP --SpliceAI <float>

Combine the above filters.

    granite cleanVCF -i file.vcf -o file.out.vcf -t tag --VEP --VEPrescue <str> <str> --SpliceAI <float>
