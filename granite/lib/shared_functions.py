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
#end def

def check_VEP(vnt_obj, idx, VEPremove, VEPrescue, VEPtag, sep='&'):
    ''' check VEP annotations from VEPtag '''
    try: val_get = vnt_obj.get_tag_value(VEPtag)
    except Exception: return False
    #end try
    trscrpt_list = val_get.split(',')
    # Get all terms
    for trscrpt in trscrpt_list:
        trscrpt_terms = set(trscrpt.split('|')[idx].split(sep))
        if trscrpt_terms.intersection(VEPrescue):
            return True
        elif trscrpt_terms.difference(VEPremove):
            return True
        #end if
    #end for
    return False
#end def

def check_spliceAI(vnt_obj, idx_list, SpAItag_list, thr=0.8):
    ''' check if SpliceAI max delta score value is over threshold thr,
    either access the max delta score directly from a tag,
    or calculate as the max of the different delta scores '''
    # if SpliceAI is within VEP annotations,
    # this fetch only the first transcript
    # and assumes all the transcripts share the same scores
    # Get max delta score
    if len(SpAItag_list) == 1: # max score from tag
        try:
            val_get = vnt_obj.get_tag_value(SpAItag_list[0])
            max_DS = float(val_get.split('|')[idx_list[0]])
        except Exception: return False
        #end try
    else: # default, multiple delta scores
        DS_list = []
        for i, SpAItag in enumerate(SpAItag_list):
            try:
                # if SpliceAI is within VEP and is at the end of Format
                # need to remove , that separate next transcript
                val_get = vnt_obj.get_tag_value(SpAItag)
                DS = val_get.split('|')[idx_list[i]]
                DS_list.append(float(DS.split(',')[0]))
            except Exception: return False
            #end try
        #end for
        max_DS = max(DS_list)
    #end if
    # Check max delta score
    if max_DS >= thr:
        return True
    #end if
    return False
#end def

def check_CLINVAR(vnt_obj, CLN_idx, CLNtag, CLNSIG_idx=0, CLNSIGtag='', CLINVARonly=[]):
    ''' check if CLNtag is present, if CLINVARonly
    check if variant has specified tags or keywords in CLNSIGtag '''
    # if clinvar is within VEP annotations,
    # this fetch only the first transcript
    # and assumes all the transcripts share the same clinvar
    try: val_get = vnt_obj.get_tag_value(CLNtag)
    except Exception: return False
    #end try
    # Check CLN
    CLN = val_get.split('|')[CLN_idx]
    if not CLN: return False
    #end if
    # Check CLNSIG
    if CLINVARonly:
        val_get = vnt_obj.get_tag_value(CLNSIGtag)
        CLNSIG = val_get.split('|')[CLNSIG_idx]
        for term in CLINVARonly:
            if term.lower() in CLNSIG.lower():
                return True
            #end if
        #end for
        return False
    #end if
    return True
#end def

#################################################################
#    Functions to get info
#################################################################
def variant_type(REF, ALT):
    ''' return variant type as snv, ins, del '''
    if len(ALT.split(',')) > 1:
        return 'snv' # TO DECIDE WHAT TO DO, as snv for now
    elif len(REF) > len(ALT) or ALT == '*':
        return 'del'
    elif len(REF) < len(ALT) or REF == '*':
        return 'ins'
    #end if
    return 'snv'
#end def

def variant_type_ext(REF, ALT):
    ''' return variant type as snv, ins, del,
    mnv (multi-nucleotide), mav (multi-allelic) '''
    # multi-allelic variant
    if len(ALT.split(',')) > 1:
        return 'mav'
    #end if
    # normal variant
    if len(REF) < len(ALT) or REF == '*':
        return 'ins'
    elif len(REF) > len(ALT) or ALT == '*':
        return 'del'
    elif len(REF) == len(ALT) and len(REF) != 1:
        return 'mnv'
    #end if
    return 'snv'
#end def

def allele_frequency(vnt_obj, aftag, idx=0):
    ''' return allele frequency for variant from aftag in INFO,
    return 0. if tag is missing or value is not a float '''
    try:
        return float(vnt_obj.get_tag_value(aftag).split('|')[idx])
    except Exception:
        return 0.
    #end try
#end def

def VEP_field(vnt_obj, idx, VEPtag):
    ''' return list of annotations at idx across all transcripts from VEPtag '''
    try: val_get = vnt_obj.get_tag_value(VEPtag)
    except Exception: return []
    #end try
    trscrpt_list = val_get.split(',')
    return [trscrpt.split('|')[idx] for trscrpt in trscrpt_list]
#end def

def get_tag_idx(vnt_obj, tag, idx=0):
    ''' return value from tag at idx in INFO for variant,
    return None if tag is missing '''
    try:
        return vnt_obj.get_tag_value(tag).split('|')[idx]
    except Exception:
        return None
    #end try
#end def

#################################################################
#    Functions to modify
#################################################################
def clean_VEP(vnt_obj, idx, VEPremove, VEPrescue, VEPtag, sep='&'):
    ''' clean VEP annotations from VEPtag '''
    try: val_get = vnt_obj.get_tag_value(VEPtag)
    except Exception: return False
    #end try
    trscrpt_clean = []
    trscrpt_list = val_get.split(',')
    # Check Consequence terms and clean transcripts
    for trscrpt in trscrpt_list:
        cnsquence_terms = []
        trscrpt_values = trscrpt.split('|')
        for term in trscrpt_values[idx].split(sep):
            if term in VEPrescue:
                cnsquence_terms.append(term)
            elif term not in VEPremove:
                cnsquence_terms.append(term)
            #end if
        #end for
        if cnsquence_terms:
            trscrpt_values[idx] = sep.join(cnsquence_terms)
            trscrpt_clean.append('|'.join(trscrpt_values))
        #end if
    #end for
    return ','.join(trscrpt_clean)
#end def

def clean_VEP_byfield(vnt_obj, idx, values, VEPtag):
    ''' clean VEP annotations from VEPtag using a set of values,
    if VEPtag field at idx not in values remove transcript '''
    try: val_get = vnt_obj.get_tag_value(VEPtag)
    except Exception: return False
    #end try
    trscrpt_clean = []
    trscrpt_list = val_get.split(',')
    # Check field and clean transcripts
    for trscrpt in trscrpt_list:
        if trscrpt.split('|')[idx] in values:
            trscrpt_clean.append(trscrpt)
        #end if
    #end for
    return ','.join(trscrpt_clean)
#end def
