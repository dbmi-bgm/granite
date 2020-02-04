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


#################################################################
#
#    FUNCTIONS
#
#################################################################
def check_tags(vnt_obj, tags_dict):
    ''' check if variant has allowed values for the tags
    specified in tags_dict = {tag: [val, ...], ...} '''
    for tag, val_list in tags_dict.items():
        try: val_get = vnt_obj.get_tag_value(tag)
        except Exception:
            sys.exit('ERROR in whitelist check: {0} tag missing for variant:\n\t{1}\n'
                        .format(tag, vnt_obj.to_string()))
        #end try
        for val in val_list:
            if val in val_get:
                return True
            #end if
        #end for
    #end for
    return False
#end def

def check_spliceAI(vnt_obj, thr=0.8):
    ''' '''
    try: val_get = float(vnt_obj.get_tag_value('SpliceAI'))
    except Exception:
        sys.exit('ERROR in whitelist check: SpliceAI tag missing or in wrong format for variant:\n\t{0}\n'
                    .format(vnt_obj.to_string()))
    #end try
    if val_get >= thr:
        return True
    #end if
    return False
#end def

def check_CLINVAR(vnt_obj):
    ''' '''
    try: val_get = vnt_obj.get_tag_value('CLINVAR')
    except Exception: return False
    #end try
    return True
#end def

#################################################################
#    runner
#################################################################
def main(args):
    ''' '''
    # Variables
    tags_dict = {}
    is_CLINVAR = True if args['CLINVAR'] else False
    SpliceAI_thr = float(args['SpliceAI']) if args['SpliceAI'] else 0.

    # Populate tags_dict
    tags_dict = {
        'VEP': [
            'transcript_ablation',
            'splice_acceptor_variant',
            'splice_donor_variant',
            'stop_gained',
            'frameshift_variant',
            'stop_lost',
            'start_lost',
            'transcript_amplification',
            'inframe_insertion',
            'inframe_deletion',
            'missense_variant',
            'protein_altering_variant',
            'splice_region_variant',
            'incomplete_terminal_codon_variant',
            'start_retained_variant',
            'stop_retained_variant',
            'synonymous_variant',
            'coding_sequence_variant',
            'mature_miRNA_variant',
            '5_prime_UTR_variant',
            '3_prime_UTR_variant',
            'non_coding_transcript_exon_variant',
            'NMD_transcript_variant',
            'non_coding_transcript_variant',
            'TFBS_ablation',
            'TFBS_amplification',
            'TF_binding_site_variant',
            'regulatory_region_ablation',
            'regulatory_region_amplification',
            'feature_elongation',
            'regulatory_region_variant',
            'feature_truncation'
        ]
    }

    # Buffers
    fo = open(args['outputfile'], 'w')

    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf(args['inputfile'])

    # Writing header
    fo.write(vcf_obj.header.definitions)
    fo.write(vcf_obj.header.columns)

    # Reading variants and writing passed
    analyzed = 0
    for i, vnt_obj in enumerate(vcf_obj.parse_variants(args['inputfile'])):
        sys.stderr.write('\rAnalyzing variant... ' + str(i + 1))
        sys.stderr.flush()

        # # Check if chromosome is canonical and in valid format
        # if not check_chrom(vnt_obj.CHROM):
        #     continue
        # #end if
        analyzed += 1

        # Check SpliceAI
        if SpliceAI_thr:
            if check_spliceAI(vnt_obj, SpliceAI_thr)):
                fo.write(vnt_obj.to_string())
                continue
            #end if
        #end if

        # Check CLINVAR
        if is_CLINVAR:
            if check_CLINVAR(vnt_obj):
                fo.write(vnt_obj.to_string())
                continue
            #end if
        #end if

        # Check tags
        if tags_dict:
            if check_tags(vnt_obj, tags_dict):
                fo.write(vnt_obj.to_string())
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
