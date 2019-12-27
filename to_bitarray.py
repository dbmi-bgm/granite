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

def check_pos(thr_reads, cov, ref_fw, ref_rv, alt_fw, alt_rv, ins_fw, ins_rv, del_fw, del_rv):
    ''' '''
    is_snv, is_ins, is_del = False, False, False
    if not thr_reads:
        is_snv = calc_allbal(ref_fw, ref_rv, alt_fw, alt_rv)
        is_ins = calc_allbal(ref_fw, ref_rv, ins_fw, ins_rv)
        is_del = calc_allbal(ref_fw, ref_rv, del_fw, del_rv)
    else:
        is_snv = calc_reads(thr_reads, alt_fw, alt_rv)
        is_ins = calc_reads(thr_reads, ins_fw, ins_rv)
        is_del = calc_reads(thr_reads, del_fw, del_rv)
    #end if

    return is_snv, is_ins, is_del
#end def

def calc_allbal(ref_fw, ref_rv, alt_fw, alt_rv):
    ''' '''
    pass
#end def

def calc_reads(thr_reads, alt_fw, alt_rv):
    ''' '''
    if alt_fw + alt_rv >= thr_reads:
        return False
    #end if

    return True
#end def

def bitarray_tofile(bit_array, filename):
    ''' convert bitarray to bytes and write to filename,
    positions are indexed by 1 in the bitarray '''
    with open(filename, 'wb') as fo:
        bit_array.tofile(fo)
    #end with
#end def

def check_region(region):
    ''' check if region format is valid,
    if valid return length of the region '''
    if ':' in region:
        try:
            chr, strt_end = region.split(':')
            strt, end = map(int, strt_end.split('-'))
            if strt >= end:
                sys.exit('ERROR in parsing region argument: start index is larger than end index\n')
            #end if
            region_len = end - start + 1
        except Exception:
            sys.exit('ERROR in parsing region argument: the format is not recognized\n')
        #end try
    else:
        pass #need to figure out length of the chromosome here
    #end if

    return region_len
#end def

def main(args):
    ''' '''
    # Variables
    region = args['region']
    thr_reads = int(args['thr_reads']) if args['thr_reads'] else 0
    thr_bams = int(args['thr_bams'])

    # Check region
    region_len =  check_region(region)

    # Initialize bitarray
    bit_array_snv = bitarray.bitarray(region_len + 1) # +1 to index positions in bitarray by 1
    bit_array_ins = bitarray.bitarray(region_len + 1)
    bit_array_del = bitarray.bitarray(region_len + 1)

    # Set all bits to 0
    bit_array_snv.setall(False)
    bit_array_ins.setall(False)
    bit_array_del.setall(False)

    # Opening buffers to read
    buffers = [tabix_IT(filename, region) for filename in args['inputfiles']]

    bams_snv, bams_ins, bams_del = 0, 0, 0
    p, tmp_chr, tmp_pos = 1, 0, '' # p is a pointer to current position
                                   # indexing positions by 1
    while True:
        bams_snv, bams_ins, bams_del = 0, 0, 0 # new position
                                               # reset bams counts
        # Check first bam
        try: chr, pos, cov, ref_fw, ref_rv, alt_fw, alt_rv, \
                 ins_fw, ins_rv, del_fw, del_rv = next(buffers[0])
        except: break
        #end try
        tmp_chr, tmp_pos = chr, pos
        # Check position and update bams counts
        is_snv, is_ins, is_del = \
            check_pos(thr_reads, cov, ref_fw, ref_rv, alt_fw, alt_rv, ins_fw, ins_rv, del_fw, del_rv)
        bams_snv += int(is_snv); bams_ins += int(is_ins); bams_del += int(is_del)
        # Check ramaining bams
        for i, buffer in enumerate(buffers[1:]):
            chr, pos, cov, ref_fw, ref_rv, alt_fw, alt_rv, \
                ins_fw, ins_rv, del_fw, del_rv = next(buffer)
            #Check consistency among the files
            if tmp_chr != chr or tmp_pos != pos:
                sys.exit('ERROR in file: position {0}:{1} in file {2} is not consistent with other input files\n'
                        .format(chr, pos, args['inputfiles'][i+1]))
            #end if
            # Check position and update bams counts
            is_snv, is_ins, is_del = \
                check_pos(thr_reads, cov, ref_fw, ref_rv, alt_fw, alt_rv, ins_fw, ins_rv, del_fw, del_rv)
            bams_snv += int(is_snv); bams_ins += int(is_ins); bams_del += int(is_del)
        #end for
        # Check thresholds to blacklist
        if bams_snv >= thr_bams:
            bit_array_snv[p] = 1
        #end if
        if bams_ins >= thr_bams:
            bit_array_ins[p] = 1
        #end if
        if bam_del >= thr_bams:
            bit_array_del[p] = 1
        #end if
        p += 1
    #end while
#end def


#################################################################
# MAIN
#################################################################
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-i', '--inputfiles', help='I/O: consecutive list of bam files used to calculate the blacklist positions [e.g -i file_1 file_2 ...]', nargs='+')
    parser.add_argument('-r', '--region', help='OTHER: region to be used [e.g chr1:1-10000000, 1:1-10000000, chr1, 1], chromsome name have to match the reference', required=True)
    parser.add_argument('--thr_bams', help='THRESHOLD: minimum number of bam files with at least "--thr_reads" for the alternate allele or having the variant (allelic balance call, default) to blacklist the variant [2] ', required=True)
    parser.add_argument('--thr_reads', help='THRESHOLD: minimum number of reads to count the bam file in "--thr_bams" [1]', required=False)


    args = vars(parser.parse_args())

    main(args)
