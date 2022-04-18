#!/usr/bin/env python

#################################################################
#
#    toBig
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
import h5py
import multiprocessing
import bitarray
from multiprocessing import Pool
# shared_functions as *
from granite.lib.shared_functions import *


#################################################################
#
#    FUNCTIONS
#
#################################################################
def check_pos(rdthr, abthr, cov, ref_fw, ref_rv, alt_fw, alt_rv, ins_fw, ins_rv, del_fw, del_rv):
    ''' check if position can be called as snv, insertion or delition.
    absolute reads counts or allelic balance can be used alternatively '''
    is_snv, is_ins, is_del = False, False, False
    if not rdthr:
        is_snv = __routine_allbal(abthr, cov, alt_fw, alt_rv)
        is_ins = __routine_allbal(abthr, cov, ins_fw, ins_rv)
        is_del = __routine_allbal(abthr, cov, del_fw, del_rv)
    else:
        is_snv = __routine_reads(rdthr, alt_fw, alt_rv)
        is_ins = __routine_reads(rdthr, ins_fw, ins_rv)
        is_del = __routine_reads(rdthr, del_fw, del_rv)
    #end if
    return is_snv, is_ins, is_del
#end def

def __routine_allbal(abthr, cov, alt_fw, alt_rv):
    ''' check if position can be called as alternate based on allelic balance '''
    alt_tot = alt_fw + alt_rv
    prc_alt = alt_tot / cov * 100
    if prc_alt >= abthr:
        return True
    #end if
    return False
#end def

def __routine_reads(rdthr, alt_fw, alt_rv):
    ''' check if position can be called as alternate based on absolute
    read counts '''
    if alt_fw + alt_rv >= rdthr:
        return True
    #end if
    return False
#end def

def bitarrays_toHDF5(filename):
    ''' write bitarrays to file in HDF5 format '''
    fo = h5py.File(filename, 'w')
    # bitarrays is in global scope
    for chr in bitarrays:
        # Packing and writing snv
        fo[chr + '_snv'] = bitarrays[chr]['snv']
        # Packing and writing ins
        fo[chr + '_ins'] = bitarrays[chr]['ins']
        # Packing and writing del
        fo[chr + '_del'] = bitarrays[chr]['del']
    #end for
    fo.close()
#end def

def run_region(files, fithr, rdthr, abthr, region):
    ''' '''
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
            check_pos(rdthr, abthr, cov, ref_fw, ref_rv, alt_fw, alt_rv, ins_fw, ins_rv, del_fw, del_rv)
        bams_snv += int(is_snv); bams_ins += int(is_ins); bams_del += int(is_del)
        # Check ramaining bams
        for i, buffer in enumerate(buffers[1:]):
            line_split = next(buffer)
            chr = line_split[0]
            pos, cov, ref_fw, ref_rv, alt_fw, alt_rv, \
                ins_fw, ins_rv, del_fw, del_rv = map(int, line_split[1:])
            # Check consistency among the files
            if tmp_chr != chr or tmp_pos != pos:
                raise IndexError('\nERROR in file: position {0}:{1} in file {2} is not consistent with other files\n'
                        .format(chr, pos, files[i+1]))
            #end if
            # Check position and update bams counts
            is_snv, is_ins, is_del = \
                check_pos(rdthr, abthr, cov, ref_fw, ref_rv, alt_fw, alt_rv, ins_fw, ins_rv, del_fw, del_rv)
            bams_snv += int(is_snv); bams_ins += int(is_ins); bams_del += int(is_del)
        #end for
        # Check thresholds
        if bams_snv >= fithr:
            snv.append((tmp_chr, tmp_pos))
        #end if
        if bams_ins >= fithr:
            ins.append((tmp_chr, tmp_pos))
        #end if
        if bams_del >= fithr:
            dele.append((tmp_chr, tmp_pos))
        #end if
    #end while
    return [snv, ins, dele]
#end def

def update_bitarrays(snv_ins_dele):
    ''' '''
    # bitarrays is in global scope
    for tmp_chr, tmp_pos in snv_ins_dele[0]:
        bitarrays[tmp_chr]['snv'][tmp_pos] = True
    #end for
    for tmp_chr, tmp_pos in snv_ins_dele[1]:
        bitarrays[tmp_chr]['ins'][tmp_pos] = True
    #end for
    for tmp_chr, tmp_pos in snv_ins_dele[2]:
        bitarrays[tmp_chr]['del'][tmp_pos] = True
    #end for

#################################################################
#    runner
#################################################################
def main(args):
    ''' '''

    global bitarrays

    # Variables
    fithr = int(args['fithr'])
    rdthr = int(args['rdthr']) if args['rdthr'] else 0
    abthr = int(args['abthr']) if args['abthr'] else 15
    ncores = int(args['ncores']) if args['ncores'] else 1
    files = args['file']

    # Data structures
    chr_length = {}
    regions = []
    bitarrays = {}

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
        bitarrays.setdefault(chr, {'snv': bitarray.bitarray(length + 1),
                                   'ins': bitarray.bitarray(length + 1),
                                   'del': bitarray.bitarray(length + 1)})
        bitarrays[chr]['snv'].setall(False)
        bitarrays[chr]['ins'].setall(False)
        bitarrays[chr]['del'].setall(False)
    #end for

    # Multiprocessing
    pool = Pool(ncores)
    for region in regions:
        p = pool.apply_async(run_region, (files, fithr, rdthr, abthr, region, ), callback=update_bitarrays)
        p.get()
    pool.close()
    pool.join()

    # Writing bitarrays to files
    bitarrays_toHDF5(args['outputfile'])
#end def


#################################################################
#
#    MAIN
#
#################################################################
if __name__ == "__main__":

    main()

#end if
