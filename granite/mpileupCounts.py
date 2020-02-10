#!/usr/bin/env python

#################################################################
#
#    mpileupCounts
#        note: require samtools available in environment
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
import subprocess
# mpileup_parser
from granite.lib import mpileup_parser
# fasta_parser
from granite.lib import fasta_parser



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
            sys.exit('\nERROR in reading position information: chr format ({0}) is not matching reference\n'
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

#################################################################
#    runner
#################################################################
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
                    sys.exit('\nERROR in parsing region argument: start index is larger than end index\n')
                #end if
            except Exception:
                sys.exit('\nERROR in parsing region argument: the format is not recognized\n')
            #end try
        else:
            try:
                chr = args['region']
            except Exception:
                sys.exit('\nERROR in parsing region argument: the format is not recognized\n')
            #end try
        #end if
    #end if

    # Building command line
    command_line = ['samtools', 'mpileup', '-aa']
    if is_chr: command_line += ['-r', args['region']]
    #end if
    if args['MQthr']: command_line += ['--min-MQ', str(args['MQthr'])]
    #end if
    if args['BQthr']: command_line += ['--min-BQ', str(args['BQthr'])]
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
#    MAIN
#
#################################################################
if __name__ == "__main__":

    main()

#end if
