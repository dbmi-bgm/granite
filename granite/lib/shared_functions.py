#!/usr/bin/env python

#################################################################
#
#    shared_functions
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
import bitarray
import tabix
import h5py
import numpy


#################################################################
#
#    FUNCTIONS
#
#################################################################
#################################################################
#    Functions to buffer
#################################################################
def tabix_IT(filename, region):
    ''' open buffer to bgzip indexed filename using tabix,
    return an iterator to file content (tsv rows as lists)
    for region '''
    tb = tabix.open(filename)
    return tb.querys(region)
#end def

#################################################################
#    Functions to load
#################################################################
def load_big(filename):
    ''' read big filename into bitarrays_dict with the following
    structure {key: bitarray, ...} '''
    big = h5py.File(filename, 'r')
    bitarrays_dict = {k: bitarray.bitarray() for k in big.keys()}
    for k in big.keys():
        bitarrays_dict[k].frombytes(big[k][:].tostring())
    #end for
    big.close()
    return bitarrays_dict
#end def

def bed_to_bitarray(filename):
    ''' read bed filename into bitarrays_dict with the following
    structure {chrID: bitarray, ...} '''
    chr_tmp, pos_tmp = '', set()
    bitarrays_dict = {}
    with open(filename) as fi:
        for line in fi:
            line = line.rstrip().split('\t')
            if len(line) >= 3: # valid line for position in bed format
                chr, start, end = line[0], int(line[1]), int(line[2])
                if not chr_tmp: chr_tmp = chr
                #end if
                # Set positions in bitarray for current chromosome and reset data structures for next chromosome
                if chr_tmp != chr: # next chromosome
                    bitarrays_dict.setdefault(chr_tmp, bitarray.bitarray(max(pos_tmp) + 1))
                    bitarrays_dict[chr_tmp].setall(False)
                    for i in pos_tmp:
                        bitarrays_dict[chr_tmp][i] = True
                    #end for
                    chr_tmp, pos_tmp = chr, set()
                #end if
                # Adding new positions
                i = 0
                while (start + i) < end:
                    pos_tmp.add(start + i + 1) # +1 to index by one
                                               # 'The first base in a chromosome is numbered 0'
                                               # 'The end position in each BED feature is one-based'
                    i += 1
                #end while
            #end if
        #end for
        # Set positions into bitarray for last chromosome
        bitarrays_dict.setdefault(chr_tmp, bitarray.bitarray(max(pos_tmp) + 1))
        bitarrays_dict[chr_tmp].setall(False)
        for i in pos_tmp:
            bitarrays_dict[chr_tmp][i] = True
        #end for
    #end with
    return bitarrays_dict
#end def

#################################################################
#    Functions to write
#################################################################
def bitarray_tofile(bit_array, filename):
    ''' convert bit_array (bitarray) to bytes and write to filename '''
    with open(filename, 'wb') as fo:
        bit_array.tofile(fo)
    #end with
#end def

#################################################################
#    Functions to check
#################################################################
def check_region(region, chr_dict):
    ''' check if chromosme and region format are valid,
    chr_dict follow the structure {chrID: ..., ...} '''
    # Parse and check if region is valid
    if ':' in region:
        try:
            chr, strt_end = region.split(':')
            strt, end = map(int, strt_end.split('-'))
            if strt >= end:
                raise ValueError('\nERROR in parsing region, in region {0} starting index is larger than ending index\n'
                        .format(region))
            #end if
        except Exception:
            raise ValueError('\nERROR in parsing region, region {0} format is not recognized\n'
                    .format(region))
        #end try
    else:
        chr = region
    #end if
    # Check if chr is valid
    if not chr in chr_dict:
        raise ValueError('\nERROR in parsing region, {0} is not a valid chromosome format\n'
                .format(chr))
    #end if
#end def

def check_chrom(chrom):
    ''' check if chromosome is canonical and in a valid format '''
    chrom_repl = chrom.replace('chr', '')
    if chrom_repl in {'M', 'MT', 'X', 'Y'}:
        return True
    else:
        try:
            int_chrom_repl = int(chrom_repl)
        except Exception:
            return False
        #end try
        if int_chrom_repl > 0 and int_chrom_repl < 23:
            return True
        #end if
    #end if
    return False
#end de

#################################################################
#    Functions to get info
#################################################################
def variant_type(REF, ALT):
    ''' return variant type as snv, ins, del '''
    if len(ALT.split(',')) > 1:
        return 'snv' # TO DECIDE WHAT TO DO, as snv for now
    elif len(REF) > 1:
        return 'del'
    elif len(ALT) > 1:
        return 'ins'
    #end if
    return 'snv'
#end def

def allele_frequency(vnt_obj, aftag, idx=0):
    ''' return allele frequency for variant from aftag in INFO,
    set allele frequency to 0. if aftag missing or value not a float '''
    try:
        return float(vnt_obj.get_tag_value(aftag).split('|')[idx])
    except:
        return 0.
    #end try
#end def
