#!/usr/bin/env python

#################################################################
#
#    geneList
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
#################################################################
#    runner
#################################################################
def main(args):
    ''' '''
    # Variables
    ENSG_set, ENSG_idx = set(), 0
    VEPtag = args['VEPtag'] if args['VEPtag'] else 'CSQ'
    is_verbose = True if args['verbose'] else False

    # Read genes list
    with open(args['geneslist']) as fi:
        for line in fi:
            line = line.rstrip()
            if line: ENSG_set.add(line)
            #end if
        #end for
    #end with

    # Buffers
    fo = open(args['outputfile'], 'w')

    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf(args['inputfile'])

    # Writing header
    vcf_obj.write_header(fo)

    # Get ENSG (Gene) index in VEP
    ENSG_idx = vcf_obj.header.get_tag_field_idx(VEPtag, 'Gene')

    # Reading variants and writing passed
    analyzed = 0
    for i, vnt_obj in enumerate(vcf_obj.parse_variants()):
        if is_verbose:
            sys.stderr.write('\rAnalyzing variant... ' + str(i + 1))
            sys.stderr.flush()
        #end if

        # # Check if chromosome is canonical and in valid format
        # if not check_chrom(vnt_obj.CHROM):
        #     continue
        # #end if
        analyzed += 1

        # Apply genes list and clean VEP
        VEP_clean = clean_VEP_byfield(vnt_obj, ENSG_idx, ENSG_set, VEPtag)

        # Remove old VEP
        vnt_obj.remove_tag_info(VEPtag)

        # Add cleaned VEP if any
        if VEP_clean:
            vnt_obj.add_tag_info('{0}={1}'.format(VEPtag, VEP_clean))
        #end if

        # Write variant
        vcf_obj.write_variant(fo, vnt_obj)
    #end for
    sys.stderr.write('\n\n...Wrote results for ' + str(analyzed) + ' analyzed variants out of ' + str(i + 1) + ' total variants\n')
    sys.stderr.flush()

    # Closing buffers
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
