#!/usr/bin/env python

#################################################################
#
#   TITLE
#       contact: berselli.michele@gmail.com
#
#################################################################


#################################################################
#   Libraries
#################################################################
import sys, os
import argparse
import bitarray
import tabix


#################################################################
#   Functions
#################################################################
def tabix_IT(filename, region):
    ''' '''
    tb = tabix.open(filename)
    return tb.querys(region)
#ene def

def main(args):
    ''' '''
    # Variables
    region = args['region']

    # Check if region format is valid
    if ':' in region:
        try:
            _, strt_end = region.split(':')
            strt, end = map(int, strt_end.split('-'))
            if strt >= end:
                sys.exit('ERROR in parsing region argument: start index is larger than end index\n')
        except Exception:
            sys.exit('ERROR in parsing region argument: the format is not recognized\n')
        #end try
    #end if

    # Opening buffers to read
    buffers = [tabix_IT(filename, region) for filename in args['inputfiles']]

    p, tmp_chr, tmp_pos = 0, 0, '' # p is a pointer to current position
    while True:
        try:
            chr, pos, cov, ref_fw, ref_rv, alt_fw, alt_rv, \
                ins_fw, ins_rv, del_fw, del_rv = next(buffers[0])
        except: break
        #end try
        tmp_chr, tmp_pos = chr, pos

        # print(0, '--', chr, pos, cov, ref_fw, ref_rv, alt_fw, alt_rv, ins_fw, ins_rv, del_fw, del_rv)

        for i, buffer in enumerate(buffers[1:]):
            chr, pos, cov, ref_fw, ref_rv, alt_fw, alt_rv, \
                ins_fw, ins_rv, del_fw, del_rv = next(buffer)
            #Check consistency among the files
            if tmp_chr != chr or tmp_pos != pos:
                sys.exit('ERROR in file: position {0}:{1} in file {2} is not consistent with other input files\n'
                        .format(chr, pos, args['inputfiles'][i+1]))
            #end if
        # print(i + 1, '--', chr, pos, cov, ref_fw, ref_rv, alt_fw, alt_rv, ins_fw, ins_rv, del_fw, del_rv)
        #end for
        p += 1
    #end while
#end def


#################################################################
# MAIN
#################################################################
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-i', '--inputfiles', help='INPUT: bam files used to calculate the blacklist positions', nargs='+')
    parser.add_argument('--bamsthr', '--outputfile', help='OTHER: ', required=False)
    parser.add_argument('--readsthr', '--reference', help='OTHER: ', required=False)
    parser.add_argument('--region', help='OTHER: region to be used [e.g chr1:1-10000000, 1:1-10000000, chr1, 1], chromsome name have to match the reference', required=True)

    args = vars(parser.parse_args())

    main(args)
