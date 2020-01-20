#!/usr/bin/env python

#################################################################
#
#    mpileupCounts
#        note: require samtools available in environment
#        Michele Berselli, berselli.michele@gmail.com
#        Harvard Medical School
#
#################################################################


#################################################################
#
#    LIBRARIES
#
#################################################################
import sys, os
import argparse
import subprocess
# mpileup_parser
try: from .lib import mpileup_parser
except Exception: from lib import mpileup_parser
#end try
# fasta_parser
try: from .lib import fasta_parser
except Exception: from lib import fasta_parser
#end try


#################################################################
#
#    FUNCTIONS
#
#################################################################
def run_mpileupParser(fi, fo, ref_dict):
    ''' '''
    mP = mpileup_parser.mpileupParser()
    first = True
    for mC in mP.generator(fi):
        try: mC.get_AD_noreference(ref_dict[mC.chr][mC.pos-1])
        except Exception:
            sys.exit('ERROR in reading position information: chr format ({0}) is not matching reference\n'
                      .format(mC.chr))
        #end try
        if first:
            first = False
            mC.write_AD(fo, header=True)
        else:
            mC.write_AD(fo)
        #end if
    #end for
#end def

def main(args):
    ''' '''
    # Initialize objects and variables
    handler = fasta_parser.FastaHandler()
    ref_dict = {} # {chr: seq, ...}
    is_chr, chr = False, ''

    # Parsing region if available
    if args['region']:
        is_chr = True
        if ':' in args['region']:
            try:
                chr, region = args['region'].split(':')
                strt, end = map(int, region.split('-'))
                if strt >= end:
                    sys.exit('ERROR in parsing region argument: start index is larger than end index\n')
                #end if
            except Exception:
                sys.exit('ERROR in parsing region argument: the format is not recognized\n')
            #end try
        else:
            try:
                chr = args['region']
            except Exception:
                sys.exit('ERROR in parsing region argument: the format is not recognized\n')
            #end try
        #end if
    #end if

    # Building command line
    command_line = ['samtools', 'mpileup', '-a']
    if is_chr: command_line += ['-r', args['region']]
    #end if
    if args['MQthr']: command_line += ['--min-MQ', args['MQthr']]
    #end if
    if args['BQthr']: command_line += ['--min-BQ', args['BQthr']]
    #end if
    command_line += [args['inputfile']]

    # Reading reference into iterator
    IT = handler.parse_generator(args['reference'])

    # Output
    fo = open(args['outputfile'], 'w')

    # Loading reference
    if is_chr:
        for header, seq in IT:
            if header.split()[0] == chr:
                ref_dict.setdefault(chr, seq)
                break
            #end if
        #end for
    else:
        for header, seq in IT:
            chr = header.split()[0]
            ref_dict.setdefault(chr, seq)
        #end for
    #end if

    # Running samtools
    pipe_in = subprocess.Popen(command_line, stdout=subprocess.PIPE)

    # Parsing mpileup
    run_mpileupParser(pipe_in.stdout, fo, ref_dict)

    # Closing files
    fo.close()
#end def


#################################################################
#
# MAIN
#
#################################################################
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='The program uses samtools to calculate statistics for pileup at each position in the specified region. Returns counts in a tsv format')

    parser.add_argument('-i', '--inputfile', help='I/O: input file in bam format', required=True)
    parser.add_argument('-o', '--outputfile', help='I/O: output file', required=True)
    parser.add_argument('-r', '--reference', help='OTHER: reference file', required=False)
    parser.add_argument('--region', help='OTHER: region to be analyzed [e.g chr1:1-10000000, 1:1-10000000, chr1, 1], chromosome name must match the reference', required=False)
    parser.add_argument('--MQthr', help='OTHER: minimum mapping quality for an alignment to be used [0]', required=False)
    parser.add_argument('--BQthr', help='OTHER: minimum base quality for a base to be considered [13]', required=False)

    args = vars(parser.parse_args())

    main(args)

#end if
