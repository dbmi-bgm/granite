#!/usr/bin/env python

#################################################################
#
#    mergeVCF
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
def definitions(vcf_obj_list):
    ''' '''
    #TODO, for now definitions of the first file will be used
    return vcf_obj_list[0].header.definitions
#end def

def columns(vcf_obj_list):
    ''' '''
    #TODO, for now columns of the first file will be used
    return vcf_obj_list[0].header.columns
#end def

#################################################################
#    runner
#################################################################
def main(args):
    ''' '''
    # Variables
    buffers = []

    # Check arguments
    files = args['file']
    if len(files) < 2:
        sys.exit('\nERROR in parsing arguments: specify at least two files to merge\n')
    #end if

    # Buffer output
    fo = open(args['outputfile'], 'w')

    # Creating Vcf objects
    vcf_obj_list = [vcf_parser.Vcf(filename) for filename in files]

    # Creating variants generators
    vnt_gen_list = []
    for i, vcf_obj in enumerate(vcf_obj_list):
        vnt_gen_list.append(vcf_obj.parse_variants(files[i]))
    #end for

    # Writing header
    fo.write(definitions(vcf_obj_list))
    fo.write(columns(vcf_obj_list))

    # Merge and write variants
    vnt_obj_list = [next(vnt_gen) for vnt_gen in vnt_gen_list]

    while True:
        # Check duplicates
        for i, vnt_obj in vnt_obj_list:
        # Check lowest

        # Write variant

        # Increase generator

    #end while

#end def


#################################################################
#
#    MAIN
#
#################################################################
if __name__ == "__main__":

    main()

#end if
