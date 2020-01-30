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
def load_bgi(filename):
    ''' read bgi filename into bitarrays_dict with the following
    structure {key: bitarray, ...} '''
    bgi = h5py.File(filename, 'r')
    bitarrays_dict = {k: bitarray.bitarray() for k in bgi.keys()}
    for k in bgi.keys():
        bitarrays_dict[k].frombytes(bgi[k][:].tostring())
    #end for
    bgi.close()
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
                raise ValueError('ERROR in parsing region, in region {0} starting index is larger than ending index\n'
                        .format(region))
            #end if
        except Exception:
            raise ValueError('ERROR in parsing region, region {0} format is not recognized\n'
                    .format(region))
        #end try
    else:
        chr = region
    #end if
    # Check if chr is valid
    if not chr in chr_dict:
        raise ValueError('ERROR in parsing region, {0} is not a valid chromosome format\n'
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
