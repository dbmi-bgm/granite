#!/usr/bin/env python

#################################################################
#
#    SVqcVCF
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
import json
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
#    General stats for variants
#################################################################
def get_stats(vnt_obj, stat_dict, ID_list):
    ''' extract information from variant for single samples,
    update counts for each sample in ID_list '''
    var_type = variant_type_sv(vnt_obj)
    if var_type in ['del', 'dup']:
        for ID in ID_list:
            _genotype(vnt_obj, ID, var_type, stat_dict)
        #end for
    #end if
#end def

def _genotype(vnt_obj, ID, var_type, stat_dict):
    ''' genotype information, update counts for ID '''
    GT = vnt_obj.get_genotype_value(ID, 'GT').replace('|', '/')
    if GT not in ['0/0', './.']: # sample has variant
        stat_dict[ID][var_type]['total'] += 1
    #end if
#end def

def to_json(stat_dict):
    ''' '''
    stat_json = {
        'total variants': []
    }

    for ID in stat_dict:
        tmp_total = {
            'name': ID,
            'total': 0
        }
        for k, v in stat_dict[ID].items():
            tmp_dict = {}
            # total variants
            if k in ['del', 'dup']:
                tmp_total.setdefault(k.upper(), v['total'])
                tmp_total['total'] += v['total']
            #end if
        #end for
        stat_json['total variants'].append(tmp_total)
    #end for
    return stat_json
#end def

#################################################################
#    runner
#################################################################
def main(args):
    ''' '''
    # Variables
    is_verbose = True if args['verbose'] else False
    stat_dict = {}

    # Buffers
    fo = open(args['outputfile'], 'w')

    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf(args['inputfile'])

    # Get list of sample IDs to use
    ID_list = args['samples'] # list of sample IDs

    # Initializing stat_dict
    for ID in ID_list:
        stat_dict.setdefault(ID, {
                            'del': {'total': 0},
                            'dup': {'total': 0}
                            })
    #end for

    # Reading variants
    analyzed = 0
    for i, vnt_obj in enumerate(vcf_obj.parse_variants()):
        if is_verbose:
            sys.stderr.write('\rAnalyzing variant... ' + str(i + 1))
            sys.stderr.flush()
        #end if
        analyzed += 1
        # Getting and updating stats
        get_stats(vnt_obj, stat_dict, ID_list)
    #end for

    # Writing output
    sys.stderr.write('\n\n...Writing results for ' + str(analyzed) + ' analyzed variants out of ' + str(i + 1) + ' total variants\n')
    sys.stderr.flush()

    # Create json
    stat_json = to_json(stat_dict)

    # Write json to file
    json.dump(stat_json, fo, indent=2, sort_keys=False)

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
