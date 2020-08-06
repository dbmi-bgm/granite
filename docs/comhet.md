### comHet
comHet is a variant calling algorithm for compound heterozygous mutations. The model uses genotype-level information in pedigree (trio) and VEP-based annotations to call possible compound heterozygous pairs. VEP annotations are used to assign variants to genes and transcripts, genotype information allows to refine calls based on inheritance mode. Calls are further flagged as "Phased" or "Unphased", where "Phased" means that genotype information supports in-trans inheritance for alternate alleles from parents.

#### Arguments
```text
    usage: granite comHet [-h] -i INPUTFILE -o OUTPUTFILE --trio TRIO [TRIO ...]
                          [--VEPtag VEPTAG] [--sep SEP] [--filter_cmpHet]
                          [--allow_undef] [--SpliceAItag SPLICEAITAG] [--impact]

    optional arguments:
      -i INPUTFILE, --inputfile INPUTFILE
                            input VCF file
      -o OUTPUTFILE, --outputfile OUTPUTFILE
                            output file to write results as VCF, use .vcf as
                            extension
      --trio TRIO [TRIO ...]
                            list of sample IDs for trio, PROBAND is required and
                            must be listed FIRST (e.g. --trio PROBAND_ID
                            [PARENT_ID] [PARENT_ID])
      --VEPtag VEPTAG       by default the program will search for "VEP" TAG
                            (VEP=<values>), use this parameter to specify a
                            different TAG to be used (e.g. CSQ)
      --sep SEP             by default the program uses "&" as separator for
                            subfields in annotating VCF (e.g.
                            ENST00000643759&ENST00000643774), use this parameter
                            to specify a different separator to be used
      --filter_cmpHet       by default the program returns all variants in the
                            input VCF file. This flag will produce a shorter
                            output containing only variants that are potential
                            compound heterozygous
      --allow_undef         by default the program ignores variants with undefined
                            genotype in parents. This flag extends the output to
                            include these cases
      --SpliceAItag SPLICEAITAG
                            by default the program will search for "SpliceAI" TAG
                            (SpliceAI=<float>), use this parameter to specify a
                            different TAG | TAG field to be used (e.g. DS_DG)
      --impact              use VEP "IMPACT" or "Consequence" terms to assign an
                            impact to potential compound heterozygous. If
                            available, SpliceAI and CLINVAR "CLNSIG" information
                            is used together with VEP
```

#### Input
comHet accepts files in VCF format as input. Files must contain genotype information for trio members to be used in addition to standard VCF columns. Column IDs for trio must match the sample IDs provided as argument (`--trio`). Proband genotype information is mandatory. If available, parents information will be used to improve specificity by ruling-out false calls based on inheritance mode. VEP annotations for "Gene" and "Feature" are also required in INFO column for transcripts.

Required VCF format structure:

    #CHROM   POS   ID   REF   ALT   QUAL   FILTER   INFO   FORMAT   PROBAND_ID   [MOTHER_ID]   [FATHER_ID]

#### Output
comHet generates output in VCF format. The program adds a VEP-like tag to INFO field to report information for calls associated to each variant. *comHet* stores information for all the compound heterozygous pairs (cmpHet) that involve the variant.

*comHet* tag definition (INFO):

    ##INFO=<ID=comHet,Number=.,Type=String,Description="Putative compound heterozygous pairs. Subembedded:'cmpHet':Format:'phase|gene|transcript|mate_variant'">

*comHet* tag definition (INFO) with `--impact`:

    ##INFO=<ID=comHet,Number=.,Type=String,Description="Putative compound heterozygous pairs. Subembedded:'cmpHet':Format:'phase|gene|transcript|impact_gene|impact_transcript|mate_variant'">

#### Examples
Calls *compound heterozygous* variants.

    granite comHet -i file.vcf -o file.out.vcf --trio PROBAND_ID [PARENT_ID] [PARENT_ID]

It is possible to add impact information for gene (`impact_gene`) and for shared transcripts (`impact_transcript`). `impact_gene` is the worst impact calculated at gene level while considering all its associated transcripts. `impact_transcript` is the worst impact calculated considering only transcripts that are shared between the two mates, if any. VEP annotations for "IMPACT" or "Consequence" must be provided in INFO column in order to assign an impact. If available, SpliceAI and CLINVAR "CLNSIG" information is used together with VEP to refine the assignment.

    granite comHet -i file.vcf -o file.out.vcf --trio PROBAND_ID [PARENT_ID] [PARENT_ID] --impact

It is possible to reduce the output to only variants that are potential compound heterozygous.

    granite comHet -i file.vcf -o file.out.vcf --trio PROBAND_ID [PARENT_ID] [PARENT_ID] --filter_cmpHet
