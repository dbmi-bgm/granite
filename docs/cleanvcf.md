### cleanVCF
cleanVCF allows to clean INFO field of input VCF file. The software can remove a list of TAG from INFO field, or can be used to clean VEP annotations.

#### Arguments
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
