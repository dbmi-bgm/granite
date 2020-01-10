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
import h5py
# import bitarray
import numpy
import multiprocessing
from multiprocessing import Process, Pool
from multiprocessing.sharedctypes import Array
from functools import partial


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
    sys.exit('TODO, allele balance calling not yet implemented\n')
#end def

def __routine_reads(thr_reads, alt_fw, alt_rv):
    ''' check if position can be called as alternate based on absolute
    read counts '''
    if alt_fw + alt_rv >= thr_reads:
        return True
    #end if
    return False
#end def

# def bitarray_tofile(bit_array, filename):
#     ''' convert bitarray to bytes and write to filename,
#     positions are indexed by 1 in the bitarray '''
#     with open(filename, 'wb') as fo:
#         bit_array.tofile(fo)
#     #end with
# #end def

def bitarrays_toHDF5(filename):
    ''' write bitarrays to file in HDF5 format '''
    fo = h5py.File(filename, 'w')
    # shared_arrays is in global scope
    for chr in shared_arrays:
        # Packing and writing snv
        tmp_arr = numpy.array(shared_arrays[chr]['snv'][:], dtype=bool)
        fo[chr + '_snv'] = numpy.packbits(tmp_arr.view(numpy.uint8))
        # Packing and writing ins
        tmp_arr = numpy.array(shared_arrays[chr]['ins'][:], dtype=bool)
        fo[chr + '_ins'] = numpy.packbits(tmp_arr.view(numpy.uint8))
        # Packing and writing del
        tmp_arr = numpy.array(shared_arrays[chr]['del'][:], dtype=bool)
        fo[chr + '_del'] = numpy.packbits(tmp_arr.view(numpy.uint8))
    #end for
    fo.close()
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

def run_region(files, thr_bams, thr_reads, region):
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
                        .format(chr, pos, files[i+1]))
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
    # Setting bits in shared bitarrays
    # shared_arrays is in global scope
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


#################################################################
#   main (runner)
#################################################################
def main(args):
    ''' '''
    # Global variables
    global shared_arrays

    # Variables
    thr_bams = int(args['thr_bams'])
    thr_reads = int(args['thr_reads']) if args['thr_reads'] else 0
    nthreads = int(args['nthreads']) if args['nthreads'] else 1
    files = args['inputfiles']

    # Data structures
    chr_length, shared_arrays = {}, {}
    regions = []

    # Reading chrom.sizes file
    with open(args['chromfile']) as fi:
        for line in fi:
            line = line.rstrip()
            if line:
                chr, length = line.split('\t')
                chr_length.setdefault(chr, int(length))
            #end if
        #end for
    #end with

    # Getting regions
    with open(args['regionfile']) as fi:
        for line in fi:
            line = line.rstrip() # line is a region
            if line:
                check_region(line, chr_length)
                regions.append(line)
            #end if
        #end for
    #end with

    # Initializing bitarrays data structure
    for chr, length in chr_length.items(): # +1 to index positions in bitarrays by 1
        shared_arrays.setdefault(chr, {'snv': Array(ctypes.c_bool, length + 1),
                                       'ins': Array(ctypes.c_bool, length + 1),
                                       'del': Array(ctypes.c_bool, length + 1)})
    #end for

    # Multiprocessing
    with Pool(nthreads) as pool:
        results = pool.map(partial(run_region, files, thr_bams, thr_reads), regions)
    #end with

    # Writing bitarrays to files
    filename = 'bitarrays_bamsthr-{0}'.format(thr_bams)
    if thr_reads: filename += '_readsthr-{0}.hdf5'.format(thr_reads)
    else: filename += '_allelebalance.hdf5'
    #end if
    bitarrays_toHDF5(filename)
#end def


#################################################################
# MAIN
#################################################################
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='The program calls variants by reads counts or allelic balance for single bams or multiple bams (joint calls) in the specified region. Results are stored in a binary format where bits corresponding to called positions are set to 1')

    parser.add_argument('-i', '--inputfiles', help='list of files to be used for the joint calling [e.g -i file_1 file_2 ...]. Files must follow mpileup_parser output format. Files must also have bgzip compression and a tabix generated index file', nargs='+')
    parser.add_argument('-r', '--regionfile', help='file containing regions to be used [e.g chr1:1-10000000, 1:1-10000000, chr1, 1], chromosome names must match the reference. Regions must be listed as a column', required=True)
    parser.add_argument('-f', '--chromfile', help='chrom.sizes file containing chromosomes information', required=True)
    parser.add_argument('--nthreads', help='number of threads to be used if multiple regions are specified [1]', required=False)
    parser.add_argument('--thr_bams', help='minimum number of bam files with at least "--thr_reads" for the alternate allele or having the variant (default call by allelic balance) to jointly call position', required=True)
    parser.add_argument('--thr_reads', help='minimum number of reads to count the bam file in "--thr_bams", if not specified calls are made by allelic balance', required=False)

    args = vars(parser.parse_args())

    main(args)

#end if
