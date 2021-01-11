#!/usr/bin/env python

#################################################################
#
#    blackList
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
    afthr, aftag = 0., ''
    big_dict = {}
    is_afthr = True if args['afthr'] else False
    is_bigfile = True if args['bigfile'] else False
    is_verbose = True if args['verbose'] else False

    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf(args['inputfile'])

    # Check arguments
    if is_afthr:
        afthr = float(args['afthr'])
        if args['aftag']:
            aftag, aftag_idx = vcf_obj.header.check_tag_definition(args['aftag'])
        else:
            sys.exit('\nERROR in parsing arguments: to filter by population allele frequency please specify the TAG to use\n')
        #end if
    else:
        if not is_bigfile:
            sys.exit('\nERROR in parsing arguments: to blacklist specify a BIG file and/or a threshold for population allele frequency and the TAG to use\n')
        #end if
    #end if

    # Buffers
    fo = open(args['outputfile'], 'w')

    # Loading big if specified
    if is_bigfile: big_dict = load_big(args['bigfile'])
    #end if

    # Writing header
    vcf_obj.write_header(fo)

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

        # Get allele frequency from aftag tag if requested
        if is_afthr:
            af = allele_frequency(vnt_obj, aftag, aftag_idx)
            # Check allele frequency
            if af > afthr:
                continue
            #end if
        #end if

        if is_bigfile:
            vtype = variant_type(vnt_obj.REF, vnt_obj.ALT)
            try:
                key = vnt_obj.CHROM + '_' + vtype
                is_blacklist = big_dict[key][vnt_obj.POS]
            except Exception:
                sys.exit('\nERROR in blacklist check: {0}:{1} missing in BIG file'.format(key, vnt_obj.POS))
            #end try
            if is_blacklist:
                continue
            #end if
        #end if

        # All good, pass and write variant
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
