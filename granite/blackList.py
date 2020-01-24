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
try: from .lib.shared_functions import *
except Exception: from lib.shared_functions import *
#end try
# vcf_parser
try: from .lib import vcf_parser
except Exception: from lib import vcf_parser
#end try


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
    bgi_dict = {}
    is_afthr = True if args['afthr'] else False
    is_bgifile = True if args['bgifile'] else False

    # Check arguments
    if is_afthr:
        afthr = float(args['afthr'])
        if args['aftag']: aftag = args['aftag']
        else:
            sys.exit('ERROR in parsing arguments: to filter by population allele frequency please specify the TAG to use\n')
        #end if
    else:
        if not is_bgifile:
            sys.exit('ERROR in parsing arguments: to blacklist specify a BGI file and/or a threshold for population allele frequency and the TAG to use\n')
        #end if
    #end if

    # Buffers
    fo = open(args['outputfile'], 'w')

    # Loading bgi if specified
    if is_bgifile: bgi_dict = load_bgi(args['bgifile'])
    #end if

    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf(args['inputfile'])

    # Writing header
    fo.write(vcf_obj.header.definitions)
    fo.write(vcf_obj.header.columns)

    # Reading variants and writing passed
    analyzed = 0
    for i, vnt_obj in enumerate(vcf_obj.parse_variants(args['inputfile'])):
        sys.stderr.write('Analyzing variant... ' + str(i + 1) + '\n')
        sys.stderr.flush()

        # # Check if chromosome is canonical and in valid format
        # if not check_chrom(vnt_obj.CHROM):
        #     continue
        # #end if
        analyzed += 1

        # Get allele frequency from aftag tag if requested
        if is_afthr:
            try:
                af = float(vnt_obj.get_tag_value(aftag))
            except Exception:
                sys.exit('ERROR in parsing VCF: TAG is missing or in the wrong format for variant:\n\t{0}\n'
                            .format(vnt_obj.to_string()))
            #end try
            # Check allele frequency
            if af > afthr:
                continue
            #end if
        #end if

        if is_bgifile:
            vtype = variant_type(vnt_obj.REF, vnt_obj.ALT)
            try:
                key = vnt_obj.CHROM + '_' + vtype
                is_blacklist = bgi_dict[key][vnt_obj.POS]
            except:
                sys.exit('ERROR in blacklist check: {0} missing in BGI file'.format(key))
            #end try
            if is_blacklist:
                continue
            #end if
        #end if

        # All good, pass and write variant
        fo.write(vnt_obj.to_string())
    #end for
    sys.stderr.write('\n...Writing results for ' + str(analyzed) + ' analyzed variants out of ' + str(i + 1) + ' total variants\n')
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
