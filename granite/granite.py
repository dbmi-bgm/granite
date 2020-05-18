#!/usr/bin/env python

#################################################################
#
#    granite
#        Entry point for granite from the command line
#
#        Michele Berselli
#        Harvard Medical School
#        berselli.michele@gmail.com
#
#################################################################


#################################################################
#
#    LIBRARIES
#
#################################################################
import sys, os
import argparse
# Tools
from granite import novoCaller
from granite import comHet
from granite import blackList
from granite import whiteList
from granite import mpileupCounts
from granite import toBig
from granite import rckTar
from granite import cleanVCF


#################################################################
#
#    FUNCTIONS
#
#################################################################
#################################################################
#    runner
#################################################################
def main():
    ''' command line wrapper around available tools '''
    # Adding parser and subparsers
    parser = argparse.ArgumentParser(prog='granite', description='granite is a collection of software to call, filter and work with genomic variants')
    subparsers = parser.add_subparsers(dest='func')

    # Add novoCaller to subparsers
    novoCaller_parser = subparsers.add_parser('novoCaller', description='Bayesian de novo variant caller',
                                                help='Bayesian de novo variant caller')

    novoCaller_parser.add_argument('-i', '--inputfile', help='input VCF file', type=str, required=True)
    novoCaller_parser.add_argument('-o', '--outputfile', help='output file to write results as VCF, use .vcf as extension', type=str, required=True)
    novoCaller_parser.add_argument('-u', '--unrelatedfiles', help='TSV index file containing SampleID<TAB>Path/to/file for unrelated files used to train the model (BAM or bgzip and tabix indexed RCK)', type=str, required=True)
    novoCaller_parser.add_argument('-t', '--triofiles', help='TSV index file containing SampleID<TAB>Path/to/file for family files, the PROBAND must be listed as LAST (BAM or bgzip and tabix indexed RCK)', type=str, required=True)
    novoCaller_parser.add_argument('--ppthr', help='threshold to filter by posterior probabilty for de novo calls (>=) [0]', type=float, required=False)
    novoCaller_parser.add_argument('--afthr', help='threshold to filter by population allele frequency (<=) [1]', type=float, required=False)
    novoCaller_parser.add_argument('--afthr_unrelated', help='threshold to filter by allele frequency calculated among unrelated (<=) [1]', type=float, required=False)
    novoCaller_parser.add_argument('--aftag', help='TAG (TAG=<float>) to be used to filter by population allele frequency', type=str, required=False)
    novoCaller_parser.add_argument('--bam', help='by default the program expect bgzip and tabix indexed RCK files for "--triofiles" and "--unrelatedfiles", add this flag if files are in BAM format instead (SLOWER)', action='store_true', required=False)
    novoCaller_parser.add_argument('--MQthr', help='(only with "--bam") minimum mapping quality for an alignment to be used (>=) [0]', type=int, required=False)
    novoCaller_parser.add_argument('--BQthr', help='(only with "--bam") minimum base quality for a base to be considered (>=) [0]', type=int, required=False)
    novoCaller_parser.add_argument('--ADthr', help='threshold to filter by alternate allele depth in parents. This will ignore and set to "0" the posterior probability for variants with a number of alternate reads in parents higher than specified value', type=int, required=False)
    novoCaller_parser.add_argument('--verbose', help='show progress status in terminal', action='store_true', required=False)

    # Add comHet to subparsers
    comHet_parser = subparsers.add_parser('comHet', description='',
                                                help='')

    comHet_parser.add_argument('-i', '--inputfile', help='input VCF file', type=str, required=True)
    comHet_parser.add_argument('-o', '--outputfile', help='output file to write results as VCF, use .vcf as extension', type=str, required=True)
    comHet_parser.add_argument('--trio', help='list of sample IDs for trio, PROBAND is required and must be listed FIRST (e.g. --trio PROBAND_ID [PARENT_ID] [PARENT_ID])', nargs='+', required=True)
    comHet_parser.add_argument('--VEPtag', help='by default the program will search for "VEP" TAG (VEP=<values>), use this parameter to specify a different TAG to be used', type=str, required=False)
    comHet_parser.add_argument('--sep', help='by default the program uses "&" as separator for subfields in annotating VCF (e.g. ENST00000643759&ENST00000643774), use this parameter to specify a different separator to be used', type=str, required=False)
    comHet_parser.add_argument('--filter_cmpHet', help='by default the program returns all variants in the input VCF file. This flag will produce a shorter output containing only variants that are potential compound heterozygous', action='store_true', required=False)
    comHet_parser.add_argument('--allow_undef', help='by default the program ignores variants with undefined genotype in parents. This flag extends the output to include these cases', action='store_true', required=False)
    comHet_parser.add_argument('--verbose', help='show progress status in terminal', action='store_true', required=False)

    # Add mpileupCounts to subparsers
    mpileupCounts_parser = subparsers.add_parser('mpileupCounts', description='samtools wrapper to calculate reads statistics for pileup at each position',
                                                    help='samtools wrapper to calculate reads statistics for pileup at each position')

    mpileupCounts_parser.add_argument('-i', '--inputfile', help='input file in BAM format', type=str, required=True)
    mpileupCounts_parser.add_argument('-o', '--outputfile', help='output file to write results as RCK format (TSV), use .rck as extension', type=str, required=True)
    mpileupCounts_parser.add_argument('-r', '--reference', help='reference file in FASTA format', type=str, required=True)
    mpileupCounts_parser.add_argument('--region', help='region to be analyzed [e.g. chr1:1-10000000, 1:1-10000000, chr1, 1], chromosome name must match the reference', type=str, required=False)
    mpileupCounts_parser.add_argument('--MQthr', help='minimum mapping quality for an alignment to be used (>=) [0]', type=int, required=False)
    mpileupCounts_parser.add_argument('--BQthr', help='minimum base quality for a base to be considered (>=) [13]', type=int, required=False)

    # Add blackList to subparsers
    blackList_parser = subparsers.add_parser('blackList', description='utility to blacklist and filter out variants from input VCF file based on positions set in BIG format file and/or population allele frequency',
                                                help='utility to blacklist and filter out variants from input VCF file based on positions set in BIG format file and/or population allele frequency')

    blackList_parser.add_argument('-i', '--inputfile', help='input VCF file', type=str, required=True)
    blackList_parser.add_argument('-o', '--outputfile', help='output file to write results as VCF, use .vcf as extension', type=str, required=True)
    blackList_parser.add_argument('-b', '--bigfile', help='BIG format file with positions set for blacklist', type=str, required=False)
    blackList_parser.add_argument('--aftag', help='TAG (TAG=<float>) to be used to filter by population allele frequency', type=str, required=False)
    blackList_parser.add_argument('--afthr', help='threshold to filter by population allele frequency (<=) [1]', type=float, required=False)
    blackList_parser.add_argument('--verbose', help='show progress status in terminal', action='store_true', required=False)

    # Add whiteList to subparsers
    whiteList_parser = subparsers.add_parser('whiteList', description='utility to whitelist and select a subset of variants from input VCF file based on specified annotations and positions',
                                                help='utility to whitelist and select a subset of variants from input VCF file based on specified annotations and positions')

    whiteList_parser.add_argument('-i', '--inputfile', help='input VCF file', type=str, required=True)
    whiteList_parser.add_argument('-o', '--outputfile', help='output file to write results as VCF, use .vcf as extension', type=str, required=True)
    whiteList_parser.add_argument('--SpliceAI', help='threshold to whitelist variants by SpliceAI value (>=)', type=float, required=False)
    whiteList_parser.add_argument('--CLINVAR', help='flag to whitelist all variants with a CLINVAR entry', action='store_true', required=False)
    whiteList_parser.add_argument('--CLINVARonly', help='CLINVAR "CLINSIG" terms or keywords to be saved. Sets for whitelist only CLINVAR variants with specified terms or keywords', nargs='+', required=False)
    whiteList_parser.add_argument('--CLINVARtag', help='by default the program will search for "CLINVAR" TAG (CLINVAR=<values>), use this parameter to specify a different TAG to be used', type=str, required=False)
    whiteList_parser.add_argument('--VEP', help='use VEP "Consequence" annotations to whitelist exonic and relevant variants (removed by default variants in intronic, intergenic, or regulatory regions)', action='store_true', required=False)
    whiteList_parser.add_argument('--VEPtag', help='by default the program will search for "VEP" TAG (VEP=<values>), use this parameter to specify a different TAG to be used', type=str, required=False)
    whiteList_parser.add_argument('--VEPrescue', help='additional terms to overrule removed flags to rescue and whitelist variants', nargs='+', required=False)
    whiteList_parser.add_argument('--VEPremove', help='additional terms to be removed', nargs='+', required=False)
    whiteList_parser.add_argument('--VEPsep', help='by default the program expects "&" as separator for subfields in VEP (e.g. intron_variant&splice_region_variant), use this parameter to specify a different separator to be used', type=str, required=False)
    whiteList_parser.add_argument('--BEDfile', help='BED format file with positions to whitelist', type=str, required=False)
    whiteList_parser.add_argument('--verbose', help='show progress status in terminal', action='store_true', required=False)

    # Add cleanVCF to subparsers
    cleanVCF_parser = subparsers.add_parser('cleanVCF', description='utility to clean INFO field of input VCF file',
                                                help='utility to clean INFO field of input VCF file')

    cleanVCF_parser.add_argument('-i', '--inputfile', help='input VCF file', type=str, required=True)
    cleanVCF_parser.add_argument('-o', '--outputfile', help='output file to write results as VCF, use .vcf as extension', type=str, required=True)
    cleanVCF_parser.add_argument('-t', '--tag', help='TAG to be removed from INFO field. Specify multiple TAGs as: "-t TAG -t TAG -t ..."', action='append', required=False)
    cleanVCF_parser.add_argument('--VEP', help='clean VEP "Consequence" annotations (removed by default terms for intronic, intergenic, or regulatory regions from annotations)', action='store_true', required=False)
    cleanVCF_parser.add_argument('--VEPtag', help='by default the program will search for "VEP" TAG (VEP=<values>), use this parameter to specify a different TAG to be used', type=str, required=False)
    cleanVCF_parser.add_argument('--VEPrescue', help='additional terms to overrule removed flags to rescue annotations', nargs='+', required=False)
    cleanVCF_parser.add_argument('--VEPremove', help='additional terms to be removed from annotations', nargs='+', required=False)
    cleanVCF_parser.add_argument('--VEPsep', help='by default the program expects "&" as separator for subfields in VEP (e.g. intron_variant&splice_region_variant), use this parameter to specify a different separator to be used', type=str, required=False)
    cleanVCF_parser.add_argument('--SpliceAI', help='threshold to save intronic annotations, from VEP "Consequence", for variants by SpliceAI value (>=)', type=float, required=False)
    cleanVCF_parser.add_argument('--verbose', help='show progress status in terminal', action='store_true', required=False)

    # Add toBig to subparsers
    toBig_parser = subparsers.add_parser('toBig', description='utility that converts counts from bgzip and tabix indexed RCK format into BIG format. Positions are "called" by reads counts or allelic balance for single or multiple files (joint calls) in specified regions',
                                                help='utility that converts counts from bgzip and tabix indexed RCK format into BIG format. Positions are "called" by reads counts or allelic balance for single or multiple files (joint calls) in specified regions')

    toBig_parser.add_argument('-f', '--file', help='file to be used to call positions. To do joint calling specify multiple files as: "-f file_1 -f file_2 -f ...". Expected bgzip and tabix indexed RCK file', action='append', required=True)
    toBig_parser.add_argument('-o', '--outputfile', help='output file to write results as BIG format (binary hdf5), use .big as extension', type=str, required=True)
    toBig_parser.add_argument('-r', '--regionfile', help='file containing regions to be used [e.g. chr1:1-10000000, 1:1-10000000, chr1, 1] listed as a column, chromosomes names must match the reference', type=str, required=True)
    toBig_parser.add_argument('-c', '--chromfile', help='chrom.sizes file containing chromosomes size information', type=str, required=True)
    toBig_parser.add_argument('--ncores', help='number of cores to be used if multiple regions are specified [1]', type=int, required=False)
    toBig_parser.add_argument('--fithr', help='minimum number of files with at least "--rdthr" for the alternate allele or having the variant, "calls" by allelic balance, to jointly "call" position (>=)', type=int, required=True)
    toBig_parser.add_argument('--rdthr', help='minimum number of alternate reads to count the file in "--fithr", if not specified "calls" are made by allelic balance (>=)', type=int, required=False)
    toBig_parser.add_argument('--abthr', help='minimum percentage of alternate reads compared to reference reads to count the file in "--fithr" when "calling" by allelic balance (>=) [15]', type=int, required=False)

    # Add rckTar to subparsers
    rckTar_parser = subparsers.add_parser('rckTar', description='utility to create a tar archive from bgzip and tabix indexed RCK files. Create an index file for the archive',
                                                help='utility to create a tar archive from bgzip and tabix indexed RCK files. Create an index file for the archive')

    rckTar_parser.add_argument('-t', '--ttar', help='target tar to write results, use .rck.tar as extension', type=str, required=True)
    rckTar_parser.add_argument('-f', '--file', help='file to be archived. Specify multiple files as: "-f SampleID_1.rck.gz -f SampleID_2.rck.gz -f ...". Files order is maintained while creating the index', action='append', required=True)

    # Subparsers map
    subparser_map = {
                    'novoCaller': novoCaller_parser,
                    'comHet': comHet_parser,
                    'blackList': blackList_parser,
                    'whiteList': whiteList_parser,
                    'cleanVCF': cleanVCF_parser,
                    'mpileupCounts': mpileupCounts_parser,
                    'toBig': toBig_parser,
                    'rckTar': rckTar_parser
                    }

    # Checking arguments
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    elif len(sys.argv) == 2:
        if sys.argv[1] in subparser_map:
            subparser_map[sys.argv[1]].print_help(sys.stderr)
            sys.exit(1)
        else:
            parser.print_help(sys.stderr)
            sys.exit(1)
        #end if
    #end if
    args = vars(parser.parse_args())

    # Call the right tool
    if args['func'] == 'novoCaller':
        novoCaller.main(args)
    elif args['func'] == 'comHet':
        comHet.main(args)
    elif args['func'] == 'blackList':
        blackList.main(args)
    elif args['func'] == 'whiteList':
        whiteList.main(args)
    elif args['func'] == 'mpileupCounts':
        mpileupCounts.main(args)
    elif args['func'] == 'toBig':
        toBig.main(args)
    elif args['func'] == 'rckTar':
        rckTar.main(args)
    elif args['func'] == 'cleanVCF':
        cleanVCF.main(args)
    #end if
#end def


#################################################################
#
#    MAIN
#
#################################################################
if __name__ == "__main__":

    main()

#end if
