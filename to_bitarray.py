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

def check_region(region):
    ''' check if region format is valid,
    if valid return length of the region '''
    hg_19 = {'1': 249250621, '2': 243199373, '3': 198022430, '4': 191154276,
             '5': 180915260, '6': 171115067, '7': 159138663, 'X': 155270560,
             '8': 146364022, '9': 141213431, '10': 135534747, '11': 135006516,
             '12': 133851895, '13': 115169878, '14': 107349540, '15': 102531392,
             '16': 90354753, '17': 81195210, '18': 78077248, '20': 63025520,
             'Y': 59373566, '19': 59128983, '22': 51304566, '21': 48129895,
             'MT': 16569}
    build_38 = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
                'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
                'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
                'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
                'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
                'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415,
                'chrM': 16569}
    # Parse and check region to return the length
    if ':' in region:
        try:
            chr, strt_end = region.split(':')
            strt, end = map(int, strt_end.split('-'))
            if strt >= end:
                sys.exit('ERROR in parsing region argument: start index is larger than end index\n')
            #end if
        except Exception:
            sys.exit('ERROR in parsing region argument: the format is not recognized\n')
        #end try
    else:
        chr = region
        #end if
    #end if

    if chr in build_38:
        return build_38[chr]
    elif chr in hg_19:
        return hg_19[chr]
    else:
        sys.exit('ERROR in parsing region argument: not a valid chromosome format\n')
    #end if
#end def

def main(args):
    ''' '''
    # Variables
    region = args['region']
    thr_reads = int(args['thr_reads']) if args['thr_reads'] else 0
    thr_bams = int(args['thr_bams'])

    # Check region
    region_len = check_region(region)
    print(region_len)

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
    tmp_chr, tmp_pos = 0, 0
    while True:
        bams_snv, bams_ins, bams_del = 0, 0, 0 # new position
                                               # reset bams counts
        # Check first bam
        try: chr, pos, cov, ref_fw, ref_rv, alt_fw, alt_rv, \
                 ins_fw, ins_rv, del_fw, del_rv = map(int, next(buffers[0]))
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
                ins_fw, ins_rv, del_fw, del_rv = map(int, next(buffer))
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
            bit_array_snv[pos] = 1
        #end if
        if bams_ins >= thr_bams:
            bit_array_ins[pos] = 1
        #end if
        if bams_del >= thr_bams:
            bit_array_del[pos] = 1
        #end if
    #end while

    # TEST
    bit_9996 = bitarray.bitarray(9996)
    bit_9996.setall(False)
    if bit_array_snv[:9996] == bit_9996:
        print('YASS')
    if bit_array_ins[:9996] == bit_9996:
        print('YASS')
    if bit_array_del[:9996] == bit_9996:
        print('YASS')
    #end if block
    print('bit_array_snv', bit_array_snv[9996:10015])
    print('bit_array_ins', bit_array_ins[9996:10015])
    print('bit_array_del', bit_array_del[9996:10015])
    bit_249240606 = bitarray.bitarray(249240607)
    bit_249240606.setall(False)
    if bit_array_snv[10015:] == bit_249240607:
        print('YASS')
    if bit_array_ins[10015:] == bit_249240607:
        print('YASS')
    if bit_array_del[10015:] == bit_249240607:
        print('YASS')
    #end if block

    # Writing bitarrays to files
    filename = 'blacklist_' + region + '_bams{0}'.format(thr_bams)
    if thr_reads: filename += '_reads{0}'.format(thr_reads)
    else: filename += '_allelebalance'
    #end if
    with open(filename + '_snv.bin', 'wb') as fo:
        bit_array_snv.tofile(fo)
    #end with
    with open(filename + '_ins.bin', 'wb') as fo:
        bit_array_ins.tofile(fo)
    #end with
    with open(filename + '_del.bin', 'wb') as fo:
        bit_array_del.tofile(fo)
    #end with
#end def


#################################################################
# MAIN
#################################################################
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-i', '--inputfiles', help='I/O: consecutive list of bam files used to calculate the blacklist positions [e.g -i file_1 file_2 ...]', nargs='+')
    parser.add_argument('-r', '--region', help='OTHER: region to be used [e.g chr1:1-10000000, 1:1-10000000, chr1, 1], chromsome name have to match the reference', required=True)
    parser.add_argument('--thr_bams', help='THRESHOLD: minimum number of bam files with at least "--thr_reads" for the alternate allele or having the variant (allelic balance call, default) to blacklist the variant', required=True)
    parser.add_argument('--thr_reads', help='THRESHOLD: minimum number of reads to count the bam file in "--thr_bams"', required=False)


    args = vars(parser.parse_args())

    main(args)
