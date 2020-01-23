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
    ''' read bgi filename into dict_bitarrays with the following
    structure {key: bitarray, ...} '''
    bgi = h5py.File(filename)
    dict_bitarrays = {k: bitarray.bitarray() for k in bgi.keys()}
    for k in bgi.keys():
        dict_bitarrays[k].frombytes(bgi[k][:].tostring())
    #end for
    bgi.close()
    return dict_bitarrays
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
    if not chr in chr_dict:
        sys.exit('ERROR in parsing region: {0} is not a valid chromosome format\n'
                .format(chr))
    #end if
#end def
