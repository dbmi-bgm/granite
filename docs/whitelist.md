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
