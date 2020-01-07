#!/usr/bin/env python

#################################################################
#
#   to_bitarray
#       contact: berselli.michele@gmail.com
#
#################################################################


#################################################################
#   Libraries
#################################################################
import sys, os
import argparse
import tabix
import ctypes
import bitarray
import numpy as np
import multiprocessing
from multiprocessing import Process, Pool
from multiprocessing.sharedctypes import Array


#################################################################
#   Functions
#################################################################
def tabix_IT(filename, region):
    ''' open buffer to bgzip indexed file using tabix,
    return an iterator to file content (tsv rows as lists) '''
    tb = tabix.open(filename)
    return tb.querys(region)
#ene def

def check_pos(thr_reads, cov, ref_fw, ref_rv, alt_fw, alt_rv, ins_fw, ins_rv, del_fw, del_rv):
    ''' check if position can be called as snv, insertion or delition.
    absolute reads counts or allelic balance can be used alternatively '''
    is_snv, is_ins, is_del = False, False, False
    if not thr_reads:
        is_snv = __routine_allbal(ref_fw, ref_rv, alt_fw, alt_rv)
        is_ins = __routine_allbal(ref_fw, ref_rv, ins_fw, ins_rv)
        is_del = __routine_allbal(ref_fw, ref_rv, del_fw, del_rv)
    else:
        is_snv = __routine_reads(thr_reads, alt_fw, alt_rv)
        is_ins = __routine_reads(thr_reads, ins_fw, ins_rv)
        is_del = __routine_reads(thr_reads, del_fw, del_rv)
    #end if
    return is_snv, is_ins, is_del
#end def

def __routine_allbal(ref_fw, ref_rv, alt_fw, alt_rv):
    ''' check if position can be called as alternate based on allelic balance '''
    pass
#end def

def __routine_reads(thr_reads, alt_fw, alt_rv):
    ''' check if position can be called as alternate based on absolute
    read counts '''
    if alt_fw + alt_rv >= thr_reads:
        return True
    #end if
    return False
#end def

def bitarray_tofile(bit_array, filename):
    ''' convert bitarray to bytes and write to filename,
    positions are indexed by 1 in the bitarray '''
    with open(filename, 'wb') as fo:
        bit_array.tofile(fo)
    #end with
#end def

def check_region(region, chr_length):
    ''' check if chromosme and region format are valid '''
    # Parse and check if region is valid
    if ':' in region:
        try:
            chr, strt_end = region.split(':')
            strt, end = map(int, strt_end.split('-'))
            if strt >= end:
                sys.exit('ERROR in parsing region: in region {0} starting index is larger than ending index\n'
                        .format(region))
            #end if
        except Exception:
            sys.exit('ERROR in parsing region: region {0} format is not recognized\n'
                    .format(region))
        #end try
    else:
        chr = region
    #end if
    # Check if chr is valid
    if not chr in chr_length:
        sys.exit('ERROR in parsing region: {0} is not a valid chromosome format\n'
                .format(chr))
    #end if
#end def

def run_region(shared_arrays, files, region, thr_bams, thr_reads):
    ''' '''
    # Variables
    snv, ins, dele = [], [], []
    # Opening buffers to region
    buffers = [tabix_IT(filename, region) for filename in files]
    # Reading from buffers and processing
    bams_snv, bams_ins, bams_del = 0, 0, 0
    tmp_chr, tmp_pos = '', 0
    while True:
        bams_snv, bams_ins, bams_del = 0, 0, 0 # new position
                                               # reset bams counts
        # Check first bam
        try:
            line_split = next(buffers[0])
            chr = line_split[0]
            pos, cov, ref_fw, ref_rv, alt_fw, alt_rv, \
                ins_fw, ins_rv, del_fw, del_rv = map(int, line_split[1:])
        except Exception: break
        #end try
        tmp_chr, tmp_pos = chr, pos
        # Check position and update bams counts
        is_snv, is_ins, is_del = \
            check_pos(thr_reads, cov, ref_fw, ref_rv, alt_fw, alt_rv, ins_fw, ins_rv, del_fw, del_rv)
        bams_snv += int(is_snv); bams_ins += int(is_ins); bams_del += int(is_del)
        # Check ramaining bams
        for i, buffer in enumerate(buffers[1:]):
            line_split = next(buffer)
            chr = line_split[0]
            pos, cov, ref_fw, ref_rv, alt_fw, alt_rv, \
                ins_fw, ins_rv, del_fw, del_rv = map(int, line_split[1:])
            # Check consistency among the files
            if tmp_chr != chr or tmp_pos != pos:
                sys.exit('ERROR in file: position {0}:{1} in file {2} is not consistent with other input files\n'
                        .format(chr, pos, args['inputfiles'][i+1]))
            #end if
            # Check position and update bams counts
            is_snv, is_ins, is_del = \
                check_pos(thr_reads, cov, ref_fw, ref_rv, alt_fw, alt_rv, ins_fw, ins_rv, del_fw, del_rv)
            bams_snv += int(is_snv); bams_ins += int(is_ins); bams_del += int(is_del)
        #end for
        # Check thresholds
        if bams_snv >= thr_bams:
            snv.append(tmp_pos)
        #end if
        if bams_ins >= thr_bams:
            ins.append(tmp_pos)
        #end if
        if bams_del >= thr_bams:
            dele.append(tmp_pos)
        #end if
    #end while
    # Setting bits in shared bitarray
    for idx in snv:
        shared_arrays[tmp_chr]['snv'][idx] = 1
    #end for
    for idx in ins:
        shared_arrays[tmp_chr]['ins'][idx] = 1
    #end for
    for idx in dele:
        shared_arrays[tmp_chr]['del'][idx] = 1
    #end for
#end def

def main(args):
    ''' '''
    # Variables
    thr_bams = int(args['thr_bams'])
    thr_reads = int(args['thr_reads']) if args['thr_reads'] else 0
    nthreads = int(args['nthreads']) if args['nthreads'] else 1
    files = args['inputfiles']

    # Getting regions
    regions = []
    with open(args['regionfile']) as fi:
        for line in fi:
            line = line.rstrip()
            if line: regions.append(line)
            #end if
        #end for
    #end with

    # Initializing bitarray data structure
    # with open() as fi:
    #     pass

    # Initialize bitarrays
    bit_array_snv = Array(ctypes.c_bool, 115169878 + 1) # +1 to index positions in bitarray by 1
    bit_array_ins = Array(ctypes.c_bool, 115169878 + 1)
    bit_array_del = Array(ctypes.c_bool, 115169878 + 1)

    shared_arrays = {}
    shared_arrays.setdefault('13', {'snv': bit_array_snv, 'ins': bit_array_ins, 'del': bit_array_del})
    #end with

    # Multithreading


    # proc = []
    # for region in regions:
    #     # Check region is valid
    #     check_region(region, shared_arrays)
    #     # Create and start process
    #     p = multiprocessing.Process(target=run_region, args=(shared_arrays, files, region, thr_bams, thr_reads, ))
    #     p.start()
    #     proc.append(p)
    # for p in proc:
    #     p.join()

    # Writing bitarrays to files
    filename = 'blacklist_' + 'TEST' + '_bamsthr-{0}'.format(thr_bams)
    if thr_reads: filename += '_readsthr-{0}'.format(thr_reads)
    else: filename += '_allelebalance'
    #end if
    # a = np.array(shared_arrays['13']['snv'][:], dtype=np.bool)

    bitarray_tofile(bitarray.bitarray(shared_arrays['13']['snv'][:]), filename + '_snv.bin')
    bitarray_tofile(bitarray.bitarray(shared_arrays['13']['ins'][:]), filename + '_ins.bin')
    bitarray_tofile(bitarray.bitarray(shared_arrays['13']['del'][:]), filename + '_del.bin')

#end def


#################################################################
# MAIN
#################################################################
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='The program calls variants by reads counts or allelic balance for single bams or multiple bams (joint calls) in the specified region. Results are stored as binary files where bits corresponding to called positions are set to 1. Three different files are created for SNV, INSERTIONS and DELITIONS respectively')

    parser.add_argument('-i', '--inputfiles', help='I/O: list of files to be used for the joint calling [e.g -i file_1 file_2 ...]. Files need to follow mpileup_parser output format', nargs='+')
    parser.add_argument('-r', '--regionfile', help='OTHER: file containing regions to be used [e.g chr1:1-10000000, 1:1-10000000, chr1, 1], chromosome names must match the reference. Regions must be listed as a column', required=True)
    parser.add_argument('-f', '--chromfile', help='OTHER: file containing chromosomes information', required=False)
    parser.add_argument('--nthreads', help='OTHER: number of threads to be used [1]', required=False)
    parser.add_argument('--thr_bams', help='THRESHOLD: minimum number of bam files with at least "--thr_reads" for the alternate allele or having the variant (default call by allelic balance) to jointly call position', required=True)
    parser.add_argument('--thr_reads', help='THRESHOLD: minimum number of reads to count the bam file in "--thr_bams", if not specified calls are made by allelic balance', required=False)

    args = vars(parser.parse_args())

    main(args)
