### novoCaller
novoCaller is a Bayesian variant calling algorithm for *de novo* mutations. The model uses read-level information both in pedigree (trio) and unrelated samples to rank and assign a probabilty to each call. The software represents an updated and improved implementation of the original algorithm described in [Mohanty et al. 2019](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty749/5087716).

#### Arguments
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
