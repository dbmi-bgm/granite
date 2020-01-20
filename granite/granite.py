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


#################################################################
#
#    FUNCTIONS
#
#################################################################
def main():
    ''' wrapper around programs from command line '''
    # Adding parser and subparsers
    parser = argparse.ArgumentParser(prog='granite')
    subparsers = parser.add_subparsers(dest='func')

    # Add novoCaller to subparsers
    novoCaller_parser = subparsers.add_parser('novoCaller', description='Bayesian de novo variant caller')

    novoCaller_parser.add_argument('-i', '--inputfile', help='input VCF file, must contain novoAF=<float> in INFO field to filter by allele frequency', required=True)
    novoCaller_parser.add_argument('-o', '--outputfile', help='output file to write results as VCF', required=True)
    novoCaller_parser.add_argument('-u', '--unrelatedbams', help='TSV file containing ID<TAB>Path/to/file for unrelated bam files used to train the model', required=True)
    novoCaller_parser.add_argument('-t', '--triobams', help='TSV file containing ID<TAB>Path/to/file for family bam files, the PROBAND must be listed as LAST', required=True)
    novoCaller_parser.add_argument('-p', '--postprobthr', help='threshold to filter by posterior probabilty for de novo calls [0]', required=False)
    novoCaller_parser.add_argument('-a', '--allelefreqthr', help='threshold to filter by population allele frequency [1]', required=False)

    # Add mpileupCounts to subparsers
    mpileupCounts_parser = subparsers.add_parser('mpileupCounts', description='The program uses samtools to calculate statistics for pileup at each position in the specified region. Returns counts in a tsv format')

    mpileupCounts_parser.add_argument('-i', '--inputfile', help='input file in BAM format', required=True)
    mpileupCounts_parser.add_argument('-o', '--outputfile', help='output file to write results as TSV', required=True)
    mpileupCounts_parser.add_argument('-r', '--reference', help='reference file in FASTA format', required=True)
    mpileupCounts_parser.add_argument('--region', help='region to be analyzed [e.g chr1:1-10000000, 1:1-10000000, chr1, 1], chromosome name must match the reference', required=False)
    mpileupCounts_parser.add_argument('--MQthr', help='minimum mapping quality for an alignment to be used [0]', required=False)
    mpileupCounts_parser.add_argument('--BQthr', help='minimum base quality for a base to be considered [13]', required=False)

    # Add blackList to subparsers
    # blackList_parser = subparsers.add_parser('blackList', description='')
    # parser.add_argument('-b', '--blacklist', help='BLACKLIST: tsv file containing ID<TAB>Path/to/file for bam files \
    #                                             used to filter out shared variants/artifacts', required=False)
    # parser.add_argument('--thr_bams', help='BLACKLIST: minimum number of bam files with at least "--thr_reads" to blacklist the variant [2]', required=False)
    # parser.add_argument('--thr_reads', help='BLACKLIST: minimum number of reads to count the bam file in "--thr_bams" [1]', required=False)

    # Add toBGI to so subparsers

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()

#end def


#################################################################
#
# MAIN
#
#################################################################
if __name__ == "__main__":

    main()

#end if
