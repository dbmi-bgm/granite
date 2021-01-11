#!/usr/bin/env python

#################################################################
#
#    whiteList
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
    VEPrescue, consequence_idx = {}, 0
    # VEPremove = {...} -> import from shared_vars
    # DStags = {...} -> import from shared_vars
    CLINVARonly = {}
    CLNtag = ''
    CLNSIGtag, CLNSIG_idx = '', 0
    SpAItag_list, SpAI_idx_list = [], []
    BED_bitarrays = {}
    is_VEP = True if args['VEP'] else False
    is_CLINVAR = True if args['CLINVAR'] else False
    SpliceAI_thr = float(args['SpliceAI']) if args['SpliceAI'] else 0.
    is_BEDfile = True if args['BEDfile'] else False
    VEPtag = args['VEPtag'] if args['VEPtag'] else 'CSQ'
    CLINVARtag = args['CLINVARtag'] if args['CLINVARtag'] else 'ALLELEID'
    SpliceAItag = args['SpliceAItag'] # default None
    VEPsep = args['VEPsep'] if args['VEPsep'] else '&'
    is_verbose = True if args['verbose'] else False

    # Buffers
    fo = open(args['outputfile'], 'w')

    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf(args['inputfile'])

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

    #CLINVAR
    if is_CLINVAR:
        CLNtag, CLN_idx = vcf_obj.header.check_tag_definition(CLINVARtag)
        if args['CLINVARonly']:
            CLINVARonly = {term for term in args['CLINVARonly']}
            CLNSIGtag, CLNSIG_idx = vcf_obj.header.check_tag_definition('CLNSIG')
        #end if
    elif args['CLINVARonly']:
        sys.exit('\nERROR in parsing arguments: specify the flag "--CLINVAR" to filter by CLINVAR annotations to specify tags or keywords to whitelist\n')
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

    # BED
    if is_BEDfile:
        BED_bitarrays = bed_to_bitarray(args['BEDfile'])
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

        # Check BED
        if is_BEDfile:
            try: # CHROM and POS can miss in the BED file, if that just pass to next checks
                if BED_bitarrays[vnt_obj.CHROM][vnt_obj.POS]:
                    vcf_obj.write_variant(fo, vnt_obj)
                    continue
                #end if
            except Exception: pass
            #end try
        #end if

        # Check VEP
        if is_VEP:
            if check_VEP(vnt_obj, consequence_idx, VEPremove, VEPrescue, VEPtag, VEPsep):
                vcf_obj.write_variant(fo, vnt_obj)
                continue
            #end if
        #end if

        # Check SpliceAI
        if SpliceAI_thr:
            if check_spliceAI(vnt_obj, SpAI_idx_list, SpAItag_list, SpliceAI_thr):
                vcf_obj.write_variant(fo, vnt_obj)
                continue
            #end if
        #end if

        # Check CLINVAR
        if is_CLINVAR:
            if check_CLINVAR(vnt_obj, CLN_idx, CLNtag, CLNSIG_idx, CLNSIGtag, CLINVARonly):
                vcf_obj.write_variant(fo, vnt_obj)
                continue
            #end if
        #end if
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
