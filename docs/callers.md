## Inheritance Mode Callers

### novoCaller
novoCaller is a Bayesian calling algorithm for *de novo* mutations. The model uses read-level information both in pedigree (trio) and unrelated samples to rank and assign a probabilty to each call. The software represents an updated and improved implementation of the original algorithm described in [Mohanty et al. 2019](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty749/5087716).

#### Arguments
```text
    usage: granite novoCaller [-h] -i INPUTFILE -o OUTPUTFILE -u UNRELATEDFILES -t
                              TRIOFILES [--ppthr PPTHR] [--afthr AFTHR]
                              [--aftag AFTAG] [--bam] [--MQthr MQTHR]
                              [--BQthr BQTHR] [--ADthr ADTHR]

    optional arguments:
      -i INPUTFILE, --inputfile INPUTFILE
                            input VCF file
      -o OUTPUTFILE, --outputfile OUTPUTFILE
                            output file to write results as VCF, use .vcf as
                            extension
      -u UNRELATEDFILES, --unrelatedfiles UNRELATEDFILES
                            TSV index file containing SampleID<TAB>Path/to/file
                            for unrelated files used to train the model (BAM or
                            bgzip and tabix indexed RCK)
      -t TRIOFILES, --triofiles TRIOFILES
                            TSV index file containing SampleID<TAB>Path/to/file
                            for family files, the PROBAND must be listed as LAST
                            (BAM or bgzip and tabix indexed RCK)
      --ppthr PPTHR         threshold to filter by posterior probabilty for de
                            novo calls (>=) [0]
      --afthr AFTHR         threshold to filter by population allele frequency
                            (<=) [1]
      --aftag AFTAG         TAG (TAG=<float>) or TAG field to be used to filter by
                            population allele frequency
      --bam                 by default the program expect bgzip and tabix indexed
                            RCK files for "--triofiles" and "--unrelatedfiles",
                            add this flag if files are in BAM format instead
                            (SLOWER)
      --MQthr MQTHR         (only with "--bam") minimum mapping quality for an
                            alignment to be used (>=) [0]
      --BQthr BQTHR         (only with "--bam") minimum base quality for a base to
                            be considered (>=) [0]
      --ADthr ADTHR         threshold to filter by alternate allele depth in
                            parents. This will ignore and set to "0" the posterior
                            probability for variants with a number of alternate
                            reads in parents higher than specified value
```

#### Input
novoCaller accepts files in VCF format as input. Files must contain genotype information for trio in addition to standard VCF columns. Column IDs for trio must match the sample IDs provided together with the list of RCK/BAM files (`--triofiles`).

Required VCF format structure:

    #CHROM   POS   ID   REF   ALT   QUAL   FILTER   INFO   FORMAT   PROBAND_ID   MOTHER_ID   FATHER_ID   ...

#### Trio and unrelated files
By default novoCaller expect bgzip and tabix indexed RCK files. To use BAM files directly specify `--bam` flag.

*note*: using BAM files directly will significantly slow down the software since pileup counts need to be calculated on the spot at each position and for each bam.

#### Output
novoCaller generates output in VCF format. Two new tags are used to report additional information for each call. *RSTR* stores reads counts by strand at position for reference and alternate alleles. *novoPP* stores posterior probabilty calculated for the call. Variants are sorted by posterior probability in desceding order.

*RSTR* tag definition (FORMAT):

    ##FORMAT=<ID=RSTR,Number=4,Type=Integer,Description="Read counts by strand for ref and alt alleles (Rf,Af,Rr,Ar)">

*novoPP* tag definition (INFO):

    ##INFO=<ID=novoPP,Number=1,Type=Float,Description="Posterior probability from novoCaller">

*note*: novoCaller model assumptions do not apply to unbalanced chromosomes (e.g. sex and mithocondrial chromosomes), therefore the model does not assign a posterior probabilty. When filtering by posterior probabilty (`--ppthr`), these variants are treated as if their posterior probabilty was 0.

#### Examples
Calls *de novo* variants. This will return the calls ranked and sorted by calculated posterior probabilty.

    granite novoCaller -i file.vcf -o file.out.vcf -u file.unrelatedfiles -t file.triofiles

It is possible to filter-out variants with posterior probabilty lower than `--ppthr`.

    granite novoCaller -i file.vcf -o file.out.vcf -u file.unrelatedfiles -t file.triofiles --ppthr <float>

It is possible to filter-out variants with population allele frequency higher than `--afthr`. Allele frequency must be provided for each variant in INFO column.

    granite novoCaller -i file.vcf -o file.out.vcf -u file.unrelatedfiles -t file.triofiles --afthr <float> --aftag tag

Filters can be combined.

    granite novoCaller -i file.vcf -o file.out.vcf -u file.unrelatedfiles -t file.triofiles --afthr <float> --aftag tag --ppthr <float>

### comHet
comHet is a calling algorithm for compound heterozygous mutations. The model uses genotype-level information in pedigree (trio) and VEP-based annotations to call possible compound heterozygous pairs. VEP annotations are used to assign variants to genes and transcripts, genotype information allows to refine calls based on inheritance mode. Calls are further flagged as "Phased" or "Unphased", where "Phased" means that genotype information supports in-trans inheritance for alternate alleles from parents.

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
      --VEPtag VEPTAG       by default the program will search for "CSQ" TAG
                            (CSQ=<values>), use this parameter to specify a
                            different TAG to be used (e.g. VEP)
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
comHet generates output in VCF format. The program adds a VEP-like tag to INFO field to report information for calls associated to each variant. *comHet* stores information for each compound heterozygous pair (cmpHet) that involves the variant.

*comHet* tag definition (INFO):

    ##INFO=<ID=comHet,Number=.,Type=String,Description="Putative compound heterozygous pairs. Subembedded:'cmpHet':Format:'phase|gene|transcript|mate_variant'">

*comHet* tag definition (INFO) with `--impact`:

    ##INFO=<ID=comHet,Number=.,Type=String,Description="Putative compound heterozygous pairs. Subembedded:'cmpHet':Format:'phase|gene|transcript|impact_gene|impact_transcript|mate_variant'">

A cmpHet is defined for each gene and for each possible mate variant. Multiple cmpHets are listed separated by comma.

Example:

    comHet=Phased|ENSG00000069424||STRONG_PAIR||chr1:6051661C>T,Phased|ENSG00000069424|ENST00000652845|STRONG_PAIR|STRONG_PAIR|chr1:6082358C>T,Phased|ENSG00000084636|ENST00000373672&ENST00000488897|STRONG_PAIR|STROING_PAIR|chr1:6051661G>A

All shared transcripts for a given pair are listed in `transcript` field. If the pair does not share any transcript, the field is empty.

#### Examples
Calls compound heterozygous variants.

    granite comHet -i file.vcf -o file.out.vcf --trio PROBAND_ID [PARENT_ID] [PARENT_ID]

It is possible to add impact information for gene (`impact_gene`) and for shared transcripts (`impact_transcript`). `impact_gene` is the worst impact calculated at gene level while considering all its associated transcripts. `impact_transcript` is the worst impact calculated considering only transcripts that are shared between the two mates, if any. VEP annotations for "IMPACT" or "Consequence" must be provided in INFO column in order to assign an impact. If available, SpliceAI and CLINVAR "CLNSIG" information is used together with VEP to refine the assignment.

    granite comHet -i file.vcf -o file.out.vcf --trio PROBAND_ID [PARENT_ID] [PARENT_ID] --impact

It is possible to reduce the output to only variants that are potential compound heterozygous.

    granite comHet -i file.vcf -o file.out.vcf --trio PROBAND_ID [PARENT_ID] [PARENT_ID] --filter_cmpHet

#### Impact
A variant is considered to have a potential STRONG impact if VEP impact is HIGH or MODERATE, spliceAI score is >= 0.8, or CLINVAR assignment is Pathogenic | Likely Pathogenic. If both variants are STRONG, the pair is assigned as a STRONG_PAIR. If only one of the two variants is STRONG, the pair is assigned as a MEDIUM_PAIR. If none of the variants is STRONG, the pair is assigned as a WEAK_PAIR.
