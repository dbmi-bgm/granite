## Variant Filtering

### filterByTag
filterByTag allows to filter variants from input VCF file by INFO field tags and thresholds. Users can define one or more filters on numeric, string, or boolean annotations. Each filter specifies the tag or field name, a value, an operator, type, and logic for aggregation. Multiple filters can be combined with global across-tag logic (`any` or `all`).

#### Arguments
```text
    usage: granite filterByTag [-h] -i INPUTFILE -o OUTPUTFILE
                               [-l {any, all}] -t TAG_FILTER [TAG_FILTER ...]
                               [--separator SEP] [-v]

    optional arguments:
      -i INPUTFILE, --inputfile INPUTFILE
                            input VCF file
      -o OUTPUTFILE, --outputfile OUTPUTFILE
                            output file to write results as VCF, use .vcf as extension
      -l {any, all}, --logic {any, all}
                            across-tag logic (combine multiple tag filters). 
                            Accept "any" or "all" [any]
      -t TAG_FILTER [TAG_FILTER ...], --tag TAG_FILTER [TAG_FILTER ...]
                            one or more tag filters. Quote each TAG_FILTER to protect
                            special characters

                            format:
                              'name/value/operator/type/logic[/entry=sep][/field=sep][/value=sep]'

                            components:
                              name       tag name (e.g. DP, CSQ) or
                                         field name (e.g. IMPACT, Consequence for VEP annotations)
                              value      threshold or string to compare against. For
                                         bool use placeholder "-"
                              operator   one of:
                                           ==      equal to
                                           !=      not equal to
                                           <       less than (int, float)
                                           >       greater than (int, float)
                                           <=      less than or equal to (int, float)
                                           >=      greater than or equal to (int, float)
                                           ~       substring contains (str)
                                           !~      substring does not contain (str)
                                           true    flag is set (bool)
                                           false   flag is unset (bool)
                              type       str | int | float | bool
                              logic      any | all (within-tag aggregation across entries)
                              entry=sep  entry separator within a tag, if tag has
                                         multiple entries (e.g. VEP transcripts)
                              field=sep  field separator within a tag, if tag/entry has
                                         embedded fields (e.g. VEP annotations)
                              value=sep  value separator within a field, if tag/entry
                                         has multiple values per field (e.g. VEP Consequence)

                            notes:
                              - if a numeric tag or embedded numeric value is missing in the
                                VCF INFO field, it is treated as 0
                              - if a string tag or embedded string value is missing in the
                                VCF INFO field, it is treated as empty string
                              - all string comparisons and tag matching are case-sensitive

      --separator SEP       tag separator within INFO field [;]
```

For complex annotations with multiple fields or multiple entries (e.g. transcript-level annotations from VEP), the program expects a VEP-like structure with proper field and entry definitions in the VCF header.  

**Format definition (example from VEP):**  

    ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature|...|gnomADg_AF|...">

#### Examples
Filter variants with depth >=10.

    granite filterByTag -i file.vcf -o file.out.vcf -t 'DP/10/>=/int/any'

Filter variants with gnomAD genome allele frequency <=0.01, evaluating all entries (transcripts) from VEP annotations.

    granite filterByTag -i file.vcf -o file.out.vcf -t 'gnomADg_AF/0.01/<=/float/all/field=|/entry=,/value=&'

Filter variants where a boolean PON (in panel of normal) flag is not set.

    granite filterByTag -i file.vcf -o file.out.vcf -t 'PON/-/false/bool/any'

Filter variants with an IMPACT value equal to HIGH in any entry (transcript) from VEP annotations.

    granite filterByTag -i file.vcf -o file.out.vcf -t 'IMPACT/HIGH/==/str/any/field=|/entry=,'

Combine filters with global across-tag logic to require `all` filters to be true.

    granite filterByTag -i file.vcf -o file.out.vcf -l all \
        -t 'DP/10/>=/int/any' \
           'gnomADg_AF/0.01/<=/float/all/field=|/entry=,/value=&' \
           'PON/-/false/bool/any' \
           'IMPACT/HIGH/==/str/any/field=|/entry=,'

### whiteList
whiteList allows to select and filter-in a subset of variants from input VCF file based on specified annotations and positions. The software can use provided VEP, ClinVar or SpliceAI annotations. Positions can be also specified as a BED format file.

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
      --SpliceAI SPLICEAI   threshold to whitelist variants by SpliceAI delta
                            scores value (>=)
      --SpliceAItag SPLICEAITAG
                            by default the program will search for SpliceAI delta
                            scores (DS_AG, DS_AL, DS_DG, DS_DL) to calculate the
                            max delta score for the variant. If a max value is
                            already defined, use this parameter to specify the TAG
                            | TAG field to be used
      --CLINVAR             flag to whitelist all variants with a ClinVar entry
                            [ALLELEID]
      --CLINVARonly CLINVARONLY [CLINVARONLY ...]
                            ClinVar "CLNSIG" terms or keywords to be saved. Sets
                            for whitelist only ClinVar variants with specified
                            terms or keywords
      --CLINVARtag CLINVARTAG
                            by default the program will search for ClinVar
                            "ALLELEID" TAG, use this parameter to specify a
                            different TAG to be used
      --VEP                 use VEP "Consequence" annotations to whitelist exonic
                            and relevant variants (removed by default variants in
                            intronic, intergenic, or regulatory regions)
      --VEPtag VEPTAG       by default the program will search for "CSQ" TAG
                            (CSQ=<values>), use this parameter to specify a
                            different TAG to be used (e.g. VEP)
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
Whitelists variants with ClinVar entry. If available, ClinVar annotation must be provided in INFO column.

    granite whiteList -i file.vcf -o file.out.vcf --CLINVAR

Whitelists only "Pathogenic" and "Likely_pathogenic" variants with ClinVar entry. ClinVar "CLNSIG" annotation must be provided in INFO column.

    granite whiteList -i file.vcf -o file.out.vcf --CLINVAR --CLINVARonly Pathogenic

Whitelists variants based on SpliceAI annotations. This filters in variants with SpliceAI score equal/higher than `--SpliceAI`. If available, SpliceAI annotation must be provided in INFO column.

    granite whiteList -i file.vcf -o file.out.vcf --SpliceAI <float>

Whitelists variants based on VEP "Consequence" annotations. This whitelists exonic and functional relevant variants by removing variants flagged as "intron_variant", "intergenic_variant", "downstream_gene_variant", "upstream_gene_variant", "regulatory_region_", "non_coding_transcript_". It is possible to specify additional [*terms*](https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html "VEP calculated consequences") to remove using `--VEPremove` and terms to rescue using `--VEPrescue`. To use VEP, annotation must be provided for each variant in INFO column.

    granite whiteList -i file.vcf -o file.out.vcf --VEP
    granite whiteList -i file.vcf -o file.out.vcf --VEP --VEPremove <str> <str>
    granite whiteList -i file.vcf -o file.out.vcf --VEP --VEPrescue <str> <str>
    granite whiteList -i file.vcf -o file.out.vcf --VEP --VEPrescue <str> <str> --VEPremove <str>

Whitelists variants based on positions specified as a BED format file.

    granite whiteList -i file.vcf -o file.out.vcf --BEDfile file.bed

Combine the above filters.

    granite whiteList -i file.vcf -o file.out.vcf --BEDfile file.bed --VEP --VEPrescue <str> <str> --CLINVAR --SpliceAI <float>

### blackList
blackList allows to filter-out variants from input VCF file based on positions set in BIG format file and/or provided population allele frequency. Positions can be also specified as a BED format file.

#### Arguments
```text
    usage: granite blackList [-h] -i INPUTFILE -o OUTPUTFILE [-b BIGFILE]
                             [--aftag AFTAG] [--afthr AFTHR] [--BEDfile BEDFILE]

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
      --BEDfile BEDFILE     BED format file with positions to blacklist
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
                            [--SpliceAItag SPLICEAITAG] [--filter_VEP]

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
      --VEPtag VEPTAG       by default the program will search for "CSQ" TAG
                            (CSQ=<values>), use this parameter to specify a
                            different TAG to be used (e.g. VEP)
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
                            "Consequence", for variants by SpliceAI delta scores
                            value (>=)
      --SpliceAItag SPLICEAITAG
                            by default the program will search for SpliceAI delta
                            scores (DS_AG, DS_AL, DS_DG, DS_DL) to calculate the
                            max delta score for the variant. If a max value is
                            already defined, use this parameter to specify the TAG
                            | TAG field to be used
      --filter_VEP          by default the program returns all variants in the input VCF file.
                            This flag will drop the variants with no VEP annotations after the
                            cleaning
```

#### Examples
Remove tag from INFO field.

    granite cleanVCF -i file.vcf -o file.out.vcf -t tag

Clean VEP based on VEP "Consequence" annotations. This removes annotations flagged as "intron_variant", "intergenic_variant", "downstream_gene_variant", "upstream_gene_variant", "regulatory_region_", "non_coding_transcript_". It is possible to specify additional [*terms*](https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html "VEP calculated consequences") to remove using `--VEPremove` and terms to rescue using `--VEPrescue`. VEP annotation must be provided for each variant in INFO column.

    granite cleanVCF -i file.vcf -o file.out.vcf --VEP
    granite cleanVCF -i file.vcf -o file.out.vcf --VEP --VEPremove <str> <str>
    granite cleanVCF -i file.vcf -o file.out.vcf --VEP --VEPrescue <str> <str>
    granite cleanVCF -i file.vcf -o file.out.vcf --VEP --VEPrescue <str> <str> --VEPremove <str>

The program also accepts a SpliceAI threshold that will rescue annotations for "intron_variant" by SpliceAI. SpliceAI annotation must be provided in INFO column.

    granite cleanVCF -i file.vcf -o file.out.vcf --VEP --SpliceAI <float>

Combine the above filters.

    granite cleanVCF -i file.vcf -o file.out.vcf -t tag --VEP --VEPrescue <str> <str> --SpliceAI <float>

### geneList
geneList allows to filter VEP annotations from input VCF file using a list of genes. If a transcript is not mapping to any of the genes in the list, the transcript is removed from VEP annotation in INFO field. If all transcripts are removed, the VEP tag is removed from INFO field for the variant.

#### Arguments
```text
    usage: granite geneList [-h] -i INPUTFILE -o OUTPUTFILE -g GENESLIST
                            [--VEPtag VEPTAG]

    optional arguments:
      -i INPUTFILE, --inputfile INPUTFILE
                            input VCF file
      -o OUTPUTFILE, --outputfile OUTPUTFILE
                            output file to write results as VCF, use .vcf as
                            extension
      -g GENESLIST, --geneslist GENESLIST
                            text file listing ensembl gene (ENSG) IDs for all
                            genes to save annotations for, IDs must be listed as a
                            column
      --VEPtag VEPTAG       by default the program will search for "CSQ" TAG
                            (CSQ=<values>), use this parameter to specify a
                            different TAG to be used (e.g. VEP)
```
