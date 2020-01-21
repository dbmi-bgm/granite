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
    ''' command line wrapper around available tools '''
    # Adding parser and subparsers
    parser = argparse.ArgumentParser(prog='granite')
    subparsers = parser.add_subparsers(dest='func')

    # Add novoCaller to subparsers
    novoCaller_parser = subparsers.add_parser('novoCaller', description='Bayesian de novo variant caller')

    novoCaller_parser.add_argument('-i', '--inputfile', help='input VCF file, must contain novoAF=<float> in INFO field to filter by allele frequency', required=True)
    novoCaller_parser.add_argument('-o', '--outputfile', help='output file to write results as VCF, use .vcf as extension', required=True)
    novoCaller_parser.add_argument('-u', '--unrelatedfiles', help='TSV file containing ID<TAB>Path/to/file for unrelated files used to train the model (BAM or bgzip and tabix indexed RCK)', required=True)
    novoCaller_parser.add_argument('-t', '--triofiles', help='TSV file containing ID<TAB>Path/to/file for family files, the PROBAND must be listed as LAST (BAM or bgzip and tabix indexed RCK)', required=True)
    novoCaller_parser.add_argument('--ppthr', help='threshold to filter by posterior probabilty for de novo calls [0]', required=False)
    novoCaller_parser.add_argument('--afthr', help='threshold to filter by population allele frequency [1]', required=False)

    # Add mpileupCounts to subparsers
    mpileupCounts_parser = subparsers.add_parser('mpileupCounts', description='Samtools wrapper to calculate reads statistics for pileup at each position in the specified region')

    mpileupCounts_parser.add_argument('-i', '--inputfile', help='input file in BAM format', required=True)
    mpileupCounts_parser.add_argument('-o', '--outputfile', help='output file to write results as RCK format (TSV), use .rck as extension', required=True)
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

    # Add rckTobgi to subparsers
    rckTobgi_parser = subparsers.add_parser('rckTobgi', description='Utility that converts counts from bgzip and tabix indexed RCK format into BGI format. Positions are called by reads counts or allelic balance for single or multiple files (joint calls) in the specified region')

    rckTobgi_parser.add_argument('-i', '--inputfiles', help='list of files to be used for the single/joint calling [e.g -i file_1 file_2 ...], expected bgzip and tabix indexed RCK files', nargs='+')
    rckTobgi_parser.add_argument('-o', '--outputfile', help='output file to write results as BGI format (binary hdf5), use .bgi as extension', required=True)
    rckTobgi_parser.add_argument('-r', '--regionfile', help='file containing regions to be used [e.g chr1:1-10000000, 1:1-10000000, chr1, 1] listed as a column, chromosomes names must match the reference', required=True)
    rckTobgi_parser.add_argument('-f', '--chromfile', help='chrom.sizes file containing chromosomes size information', required=True)
    rckTobgi_parser.add_argument('--ncores', help='number of cores to be used if multiple regions are specified [1]', required=False)
    rckTobgi_parser.add_argument('--bmthr', help='minimum number of bam files with at least "--rdthr" for the alternate allele or having the variant, if calls by allelic balance, to jointly call position', required=True)
    rckTobgi_parser.add_argument('--rdthr', help='minimum number of reads to count the bam file in "--bmthr", if not specified calls are made by allelic balance', required=False)

    # Checking arguments
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
