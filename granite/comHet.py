#!/usr/bin/env python

#################################################################
#
#    comHet
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
# shared_functions as *
from granite.lib.shared_functions import *
# vcf_parser
from granite.lib import vcf_parser


#################################################################
#
#    FUNCTIONS
#
#################################################################
def add_vnt_to_trscrpt():
    pass
#end def

#################################################################
#    runner
#################################################################
def main(args):
    ''' '''
    # Variables

    # Buffers
    fo = open(args['outputfile'], 'w')

    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf(args['inputfile'])

    # Data structures
    trscrpts_dict = {} # {trscrptID: [vnt_obj1, vnt_obj2], ...}


    # Reading variants
    analyzed, CHROM = 0, ''
    for i, vnt_obj in enumerate(vcf_obj.parse_variants(args['inputfile'])):
        sys.stderr.write('\rAnalyzing variant... ' + str(i + 1))
        sys.stderr.flush()

        if vnt_obj.CHROM != CHROM and trscrpts_dict:
            # DO MAGIC AND PRINT OUT
            pass
        else:

        #end if
    #end for

    # Writing output
    sys.stderr.write('\n\n...Writing results for ' + str(analyzed) + ' analyzed variants out of ' + str(i + 1) + ' total variants\n')
    sys.stderr.flush()

#end def


#################################################################
#
#    MAIN
#
#################################################################
if __name__ == "__main__":

    main()

#end if
