# granite
granite (*genomic variants filtering utilities*) is a collection of software to call, filter and work with genomic variants.

- [granite](#granite)
  * [Availability and requirements](#availability-and-requirements)
  * [File formats](#file-formats)
  * [Tools](#tools)
    + [novoCaller](#novocaller)
      - [Arguments](#arguments)
      - [Input](#input)
      - [Trio and unrelated files](#trio-and-unrelated-files)
      - [Output](#output)
      - [Examples](#examples)
    + [comHet](#comhet)
      - [Arguments](#arguments-1)
      - [Input](#input-1)
      - [Output](#output-1)
      - [Examples](#examples-1)
    + [blackList](#blacklist)
      - [Arguments](#arguments-2)
      - [Examples](#examples-2)
    + [whiteList](#whitelist)
      - [Arguments](#arguments-3)
      - [Examples](#examples-3)
    + [cleanVCF](#cleanvcf)
      - [Arguments](#arguments-4)
      - [Examples](#examples-4)
    + [mpileupCounts](#mpileupcounts)
      - [Arguments](#arguments-5)
    + [toBig](#tobig)
      - [Arguments](#arguments-6)
      - [Examples](#examples-5)
    + [rckTar](#rcktar)
      - [Arguments](#arguments-7)

&nbsp;
## Availability and requirements
A ready-to-use docker image is available to download.

    docker pull IMAGE

To run locally, Python (3.x) is required together with the following libraries:

  - [*numpy*](https://docs.scipy.org/doc/ "numpy documentation")
  - [*pysam*](https://pysam.readthedocs.io/en/latest/ "pysam documentation")
  - [*bitarray*](https://pypi.org/project/bitarray/ "bitarray documentation")
  - [*pytabix*](https://pypi.org/project/pytabix/ "pytabix documentation")
  - [*h5py*](https://www.h5py.org/ "h5py documentation")

To install libraries with pip:

    pip install numpy pysam bitarray h5py
    pip install --user pytabix

Additional software needs to be available in the environment:

  - [*samtools*](http://www.htslib.org/ "samtools documentation")
  - [*bgzip*](http://www.htslib.org/doc/bgzip.1.html "bgzip documentation")
  - [*tabix*](http://www.htslib.org/doc/tabix.1.html "tabix documentation")

To install the program, run the following command inside granite folder:

    python setup.py install

&nbsp;
## File formats
The program is compatible with standard BED, BAM and VCF formats (VCFv4.x).

### ReadCountKeeper (.rck)
RCK is a tabular format that allows to efficiently store counts by strand (ForWard-ReVerse) for reads that support REFerence allele, ALTernate alleles, INSertions or DELetions at CHRomosome and POSition. RCK files can be further compressed with *bgzip* and indexed with *tabix* for storage, portability and faster random access. 1-based.

Tabular format structure:

    #CHR   POS   COVERAGE   REF_FW   REF_RV   ALT_FW   ALT_RV   INS_FW   INS_RV   DEL_FW   DEL_REV
    13     1     23         0        0        11       12       0        0        0        0
    13     2     35         18       15       1        1        0        0        0        0

Commands to compress and index files:

    bgzip PATH/TO/FILE
    tabix -b 2 -s 1 -e 0 -c "#" PATH/TO/FILE.gz

### BinaryIndexGenome (.big)
BIG is a hdf5-based binary format that stores boolean values for each genomic position as bit arrays. Each position is represented in three complementary arrays that account for SNVs (Single-Nucleotide Variants), insertions and deletions respectively. 1-based.

hdf5 format structure:

    e.g.
    chr1_snv: array(bool)
    chr1_ins: array(bool)
    chr1_del: array(bool)
    chr2_snv: array(bool)
    ...
    ...
    chrM_del: array(bool)

*note*: hdf5 keys are built as the chromosome name based on reference (e.g. chr1) plus the suffix specifying whether the array represents SNVs (_snv), insertions (_ins) or deletions (_del).

&nbsp;
## Tools
![tools chart](docs/chart.png)

&nbsp;
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

&nbsp;
### comHet
comHet is a variant calling algorithm for *compound heterozygous* mutations. The model uses genotype-level information in pedigree (trio) and VEP-based annotations to call possible compound heterozygous pairs. VEP annotations are used to assign variants to genes and transcripts, genotype information allows to refine calls based on inheritance mode. Calls are further flagged as "Phased" or "Unphased", where "Phased" means that genotype information supports in-trans inheritance for alternate alleles from parents.

#### Arguments
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

&nbsp;
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

&nbsp;
### whiteList
whiteList allows to select and filter-in a subset of variants from input VCF file based on specified annotations and positions. The software can use provided VEP, CLINVAR or SpliceAI annotations. Positions can be also specfied as a BED format file.

#### Arguments
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

&nbsp;
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

&nbsp;
### mpileupCounts
mpileupCounts uses *samtools* to access input BAM and calculates statistics for reads pileup at each position in the specified region, returns counts in RCK format.

#### Arguments
    usage: granite mpileupCounts [-h] -i INPUTFILE -o OUTPUTFILE -r REFERENCE
                                 [--region REGION] [--MQthr MQTHR] [--BQthr BQTHR]

    optional arguments:
      -i INPUTFILE, --inputfile INPUTFILE
                            input file in BAM format
      -o OUTPUTFILE, --outputfile OUTPUTFILE
                            output file to write results as RCK format (TSV), use
                            .rck as extension
      -r REFERENCE, --reference REFERENCE
                            reference file in FASTA format
      --region REGION       region to be analyzed [e.g. chr1:1-10000000,
                            1:1-10000000, chr1, 1], chromosome name must match the
                            reference
      --MQthr MQTHR         minimum mapping quality for an alignment to be used
                            (>=) [0]
      --BQthr BQTHR         minimum base quality for a base to be considered (>=)
                            [13]

&nbsp;
### toBig
toBig converts counts from bgzip and tabix indexed RCK format into BIG format. Positions are "called" by read counts or allelic balance for single or multiple files (joint calls) in specified regions. Positions "called" are set to True (or 1) in BIG binary structure.

#### Arguments
    usage: granite toBig [-h] [-i INPUTFILE [INPUTFILE ...]] -o OUTPUTFILE -r
                         REGIONFILE -f CHROMFILE [--ncores NCORES] --fithr FITHR
                         [--rdthr RDTHR] [--abthr ABTHR]

    optional arguments:
      -f FILE, --file FILE  file to be used to call positions. To do joint calling
                            specify multiple files as: "-f file_1 -f file_2 -f ...".
                            Expected bgzip and tabix indexed RCK file
      -o OUTPUTFILE, --outputfile OUTPUTFILE
                            output file to write results as BIG format (binary
                            hdf5), use .big as extension
      -r REGIONFILE, --regionfile REGIONFILE
                            file containing regions to be used [e.g.
                            chr1:1-10000000, 1:1-10000000, chr1, 1] listed as a
                            column, chromosomes names must match the reference
      -c CHROMFILE, --chromfile CHROMFILE
                            chrom.sizes file containing chromosomes size
                            information
      --ncores NCORES       number of cores to be used if multiple regions are
                            specified [1]
      --fithr FITHR         minimum number of files with at least "--rdthr" for
                            the alternate allele or having the variant, "calls" by
                            allelic balance, to jointly "call" position (>=)
      --rdthr RDTHR         minimum number of alternate reads to count the file in
                            "--fithr", if not specified "calls" are made by
                            allelic balance (>=)
      --abthr ABTHR         minimum percentage of alternate reads compared to
                            reference reads to count the file in "--fithr" when
                            "calling" by allelic balance (>=) [15]

#### Examples
toBig can be used to calculate positions to blacklist for common variants by using unrelated samples. This command will set to `True` in BIG structure positions with allelic balance for alternate allele equal/higher than `--abthr` in more that `--fithr` samples (joint calling).

    granite toBig -f file -f file -f file -f file -f ... -o file.out.big -c file.chrom.sizes -r file.regions --fithr <int> --abthr <int>

Absolute reads count can be used instead of allelic balance to call positions. This command will set to `True` in BIG structure positions with reads count for alternate allele equal/higher than `--rdthr` in more that `--fithr` samples (joint calling).

    granite toBig -f file -f file -f file -f file -f ... -o file.out.big -c file.chrom.sizes -r file.regions --fithr <int> --rdthr <int>

&nbsp;
### rckTar
rckTar creates a tar archive from bgzip and tabix indexed RCK files. Creates an index file for the archive.

#### Arguments
    usage: granite rckTar [-h] -t TTAR -f FILE

    optional arguments:
      -t TTAR, --ttar TTAR  target tar to write results, use .rck.tar as extension
      -f FILE, --file FILE  file to be archived. Specify multiple files as: "-f
                            SampleID_1.rck.gz -f SampleID_2.rck.gz -f ...". Files
                            order is maintained while creating the index
