#!/usr/bin/env python

#################################################################
#
#    cleanVCF
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
# shared_vars
from granite.lib.shared_vars import VEPremove
from granite.lib.shared_vars import VEPSpliceAI
from granite.lib.shared_vars import DStags


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
    VEPrescue, consequence_idx = set(), 0
    # VEPremove = {...} -> import from shared_vars
    # VEPSpliceAI = {...} -> import from shared_vars
    # DStags = {...} -> import from shared_vars
    SpAItag_list, SpAI_idx_list = [], []
    is_VEP = True if args['VEP'] else False
    is_filter_VEP = True if args['filter_VEP'] else False
    VEPtag = args['VEPtag'] if args['VEPtag'] else 'CSQ'
    VEPsep = args['VEPsep'] if args['VEPsep'] else '&'
    SpliceAI_thr = float(args['SpliceAI']) if args['SpliceAI'] else 0.
    SpliceAItag = args['SpliceAItag'] # default None
    is_SpAI = False
    is_verbose = True if args['verbose'] else False

    # Buffers
    fo = open(args['outputfile'], 'w')

    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf(args['inputfile'])

    # Clean header definitions in INFO block for specified tags
    if args['tag']:
        for tag in args['tag']:
            vcf_obj.header.remove_tag_definition(tag, 'INFO')
        #end for
    #end if

    # Writing header
    vcf_obj.write_header(fo)

    # VEP
    if is_VEP:
        consequence_idx = vcf_obj.header.get_tag_field_idx(VEPtag, 'Consequence')
        if args['VEPrescue']: VEPrescue = {term for term in args['VEPrescue']}
        #end if
        if args['VEPremove']: VEPremove.update({term for term in args['VEPremove']})
        #end if
    elif args['VEPrescue'] or args['VEPremove']:
        sys.exit('\nERROR in parsing arguments: specify the flag "--VEP" to filter by VEP annotations to apply rescue terms or remove additional terms\n')
    #end if

    # SpliceAI
    if SpliceAI_thr:
        if SpliceAItag: # single tag has been specified
            tag, idx = vcf_obj.header.check_tag_definition(SpliceAItag)
            SpAItag_list.append(tag)
            SpAI_idx_list.append(idx)
        else: # search for delta scores as default
            for DStag in DStags:
                tag, idx = vcf_obj.header.check_tag_definition(DStag)
                SpAItag_list.append(tag)
                SpAI_idx_list.append(idx)
            #end for
        #end if
    #end if

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

        # Clean tags if specified
        if args['tag']:
            for tag in args['tag']:
                vnt_obj.remove_tag_info(tag)
            #end for
        #end if

        # is_SpAI reset
        is_SpAI = False

        # Check SpliceAI
        if SpliceAI_thr:
            if check_spliceAI(vnt_obj, SpAI_idx_list, SpAItag_list, SpliceAI_thr):
                is_SpAI = True
            #end if
        #end if

        # Clean VEP
        if is_VEP:
            # Get cleaned VEP
            if is_SpAI:
                VEP_clean = clean_VEP(vnt_obj, consequence_idx, VEPremove, VEPrescue.union(VEPSpliceAI), VEPtag, VEPsep)
            else:
                VEP_clean = clean_VEP(vnt_obj, consequence_idx, VEPremove, VEPrescue, VEPtag, VEPsep)
            #end if
            # Remove old VEP
            vnt_obj.remove_tag_info(VEPtag)
            # Add cleaned VEP if any
            if VEP_clean:
                vnt_obj.add_tag_info('{0}={1}'.format(VEPtag, VEP_clean))
            else: # check if is filter_VEP
                if is_filter_VEP:
                    continue
                #end if
            #end if
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
