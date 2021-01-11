#!/usr/bin/env python

#################################################################
#
#    novoCaller
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
import pysam
import tabix
import numpy as np
# shared_functions as *
from granite.lib.shared_functions import *
# vcf_parser
from granite.lib import vcf_parser
# shared_vars
from granite.lib.shared_vars import real_NA_chroms
from granite.lib.shared_vars import test_NA_chroms

# suppress warnings
import warnings
warnings.filterwarnings("ignore")


#################################################################
#
#    FUNCTIONS
#
#################################################################
#################################################################
#    GT_ordering_alternate (original)
#################################################################
def GT_ordering_alternate(ALT_count):
    ''' '''
    combos = (ALT_count + 1) * (ALT_count + 2) // 2
    ordering = np.empty([combos, 2])
    count = 0
    for a1 in range(0, ALT_count + 1):
        for a2 in range(a1, ALT_count + 1):
            ordering[count, 0] = a1
            ordering[count, 1] = a2
            count += 1
        #end for
    #end for

    return ordering
#end def

#################################################################
#    row_gen (original)
#################################################################
def row_gen(GT1, GT2, alt_count, mut_rate):
    ''' '''
    N = alt_count
    combos = (N + 1) * (N + 2) // 2
    row = np.zeros(combos)
    count = 0
    for a1 in range(N + 1):
        for a2 in range(N + 1):
            for a3 in range(N + 1):
                for a4 in range(N + 1):
                    P = 1.0
                    if a1 == GT1[0]:
                        P = P * (1 - mut_rate)
                    else:
                        P = P * mut_rate / N
                    #end if
                    if a2 == GT1[1]:
                        P = P * (1 - mut_rate)
                    else:
                        P = P * mut_rate / N
                    #end if
                    if a3 == GT2[0]:
                        P = P * (1 - mut_rate)
                    else:
                        P = P * mut_rate / N
                    #end if
                    if a4 == GT2[1]:
                        P = P * (1 - mut_rate)
                    else:
                        P = P * mut_rate / N
                    #end if
                    count += 1

                    for b1 in [a1, a2]:
                        for b2 in [a3, a4]:
                            gt_work = np.sort([b1, b2])
                            index = (2 * N + 3 - gt_work[0]) * gt_work[0] // 2 + gt_work[1] - gt_work[0]
                            row[index] = row[index] + 0.25 * P
                        #end for
                    #end for
                #end for
            #end for
        #end for
    #end for

    return row
#end def

#################################################################
#    table_gen (original)
#################################################################
def table_gen(alt_count, mut_rate):
    ''' '''
    N = alt_count
    II_prev = -1

    combos = (N + 1) * (N + 2) // 2
    table = np.zeros([combos ** 2, combos])
    for a1 in range(N + 1):
        for a2 in range(a1, N + 1):
            for a3 in range(N + 1):
                for a4 in range(a3, N + 1):
                    GT1 = [a1, a2]
                    GT2 = [a3, a4]
                    I1 = (2 * N + 3 - GT1[0]) * GT1[0] // 2 + GT1[1] - GT1[0]
                    I2 = (2 * N + 3 - GT2[0]) * GT2[0] // 2 + GT2[1] - GT2[0]
                    II = I1 * combos + I2

                    if II <= II_prev:
                        sys.exit("\nERROR in II calc!!!\n")
                    #end if

                    row = row_gen(GT1, GT2, alt_count, mut_rate)
                    table[II, :] = row
                #end for
            #end for
        #end for
    #end for

    return table
#end def

#################################################################
#    GT_likelihood_wrt_allele_calc (original)
#################################################################
def GT_likelihood_wrt_allele_calc(ALT_count):
    ''' '''
    ordering = GT_ordering_alternate(ALT_count)

    combos = (ALT_count + 1) * (ALT_count + 2) // 2
    table = np.zeros([combos, ALT_count + 1]) * 1.

    for i in range(combos):
        a1 = int(ordering[i, 0]) # casting to int to be used as indexes
        a2 = int(ordering[i, 1])
        table[i, a1] = table[i, a1] + 0.5
        table[i, a2] = table[i, a2] + 0.5
    #end for

    return table
#end def

#################################################################
#    buffering_bams
#################################################################
def buffering_bams(bams_infofile):
    ''' return a list containing reading buffers to bam files (pysam),
    return also a list with the corrisponding IDs associated
    to the bam files '''
    bamfiles, IDs  = [], []
    with open(bams_infofile) as fi:
        for line in fi:
            line_strip = line.rstrip()
            if line_strip:
                try:
                    ID, filepath = line_strip.split('\t') #ID    path/to/file
                    bamfile = pysam.AlignmentFile(filepath, "rb" )
                    IDs.append(ID)
                    bamfiles.append(bamfile)
                except Exception:
                    sys.exit('\nERROR in parsing BAMs info file\n')
                #end try
            #end if
        #end for
    #end with

    return bamfiles, IDs
#end def

#################################################################
#    buffering_rcks
#################################################################
def buffering_rcks(rcks_infofile):
    ''' return a list containing reading buffers to rck files (tabix),
    return also a list with the corrisponding IDs associated to files'''
    rckfiles, IDs  = [], []
    with open(rcks_infofile) as fi:
        for line in fi:
            line_strip = line.rstrip()
            if line_strip:
                try:
                    ID, filepath = line_strip.split('\t') #ID    path/to/file
                    rckfile = tabix.open(filepath)
                    IDs.append(ID)
                    rckfiles.append(rckfile)
                except Exception:
                    sys.exit('\nERROR in parsing RCKs info file\n')
                #end try
            #end if
        #end for
    #end with

    return rckfiles, IDs
#end def

#################################################################
#    get_ADs_bam (original refactored)
#################################################################
def get_ADs_bam(bamfile, chrom, pos, REF, MQthr, BQthr, deletion=False, insertion=False):
    ''' access bam file and return AD counts by strand for reference and alternate allele for variant,
    the way pileup is used depend on variant type (insertion, delition or snv) '''
    # coordinates in pysam are always 0-based
    #   -> current position - 1 to access the right position from pysam
    position = pos - 1
    REF_upper = REF.upper()
    ADf, ADr = np.array([0., 0]), np.array([0., 0])

    # Getting pileup info
    try: SP = bamfile.pileup(chrom, position, position + 1)
    except Exception:
        sys.exit('\nERROR in accessing BAM file: variant and file chromosome formats are not matching\n')
    #end try

    # Getting AD info
    for pileupcolumn in SP:
        if pileupcolumn.pos == position:
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    qp = pileupread.query_position
                    MQ = pileupread.alignment.mapping_quality
                    BQ = pileupread.alignment.query_qualities[qp]
                    if MQ >= MQthr and BQ >= BQthr:
                        seq = pileupread.alignment.query_sequence[qp].upper()
                        # DELETION
                        if deletion:
                            indel_val = pileupread.indel
                            if seq == REF_upper and indel_val >= 0:
                                if pileupread.alignment.is_reverse: ADr[0] += 1
                                else: ADf[0] += 1
                                #end if
                            elif seq == REF_upper and indel_val < 0:
                                if pileupread.alignment.is_reverse: ADr[1] += 1
                                else: ADf[1] += 1
                                #end if
                            #end if
                        # INSERTION
                        elif insertion:
                            indel_val = pileupread.indel
                            if seq == REF_upper and indel_val <= 0:
                                if pileupread.alignment.is_reverse: ADr[0] += 1
                                else: ADf[0] += 1
                                #end if
                            elif seq == REF_upper and indel_val > 0:
                                if pileupread.alignment.is_reverse: ADr[1] += 1
                                else: ADf[1] += 1
                                #end if
                            #end if
                        # SNV
                        else:
                            if seq == REF_upper:
                                if pileupread.alignment.is_reverse: ADr[0] += 1
                                else: ADf[0] += 1
                                #end if
                            else:
                                if pileupread.alignment.is_reverse: ADr[1] += 1
                                else: ADf[1] += 1
                                #end if
                            #end if
                        #end if
                    #end if
                #end if
            #end for
        #end if
    #end for

    return ADf, ADr
#end def

#################################################################
#    get_ADs_rck
#################################################################
def get_ADs_rck(rckfile, chrom, pos, deletion=False, insertion=False):
    ''' '''
    # rck format is 1-based following mpileup standard (1-based position on the chromosome)
    region = chrom + ':' + str(pos) + '-' + str(pos)
    chr, rec_pos, cov, ref_fw, ref_rv, alt_fw, alt_rv, \
        ins_fw, ins_rv, del_fw, del_rv = next(rckfile.querys(region))
    if pos != int(rec_pos):
        raise IndexError('\nERROR in RCK file indexing: position received is not consistent with position called\n')
    #end if
    if deletion:
        ADf = np.array([int(ref_fw), int(del_fw)])
        ADr = np.array([int(ref_rv), int(del_rv)])
    elif insertion:
        ADf = np.array([int(ref_fw), int(ins_fw)])
        ADr = np.array([int(ref_rv), int(ins_rv)])
    else:
        ADf = np.array([int(ref_fw), int(alt_fw)])
        ADr = np.array([int(ref_rv), int(alt_rv)])
    #end if

    return ADf, ADr
#end def

#################################################################
#    get_ADs_caller
#################################################################
def get_ADs_caller(file, chrom, pos, REF, ALT, MQthr, BQthr, is_bam):
    ''' define variant type and run the appropriate function to access file
    to retrieve AD counts by strand for reference and alternate allele for variant'''
    split_ALT = ALT.split(',')
    if len(split_ALT) > 1:
        if is_bam:
            return get_ADs_bam(file, chrom, pos, REF[0], MQthr, BQthr)
        else:
            return get_ADs_rck(file, chrom, pos)
        #end if
    elif len(REF) > 1:
        if is_bam:
            return get_ADs_bam(file, chrom, pos, REF[0], MQthr, BQthr, deletion=True)
        else:
            return get_ADs_rck(file, chrom, pos, deletion=True)
        #end if
    elif len(ALT) > 1:
        if is_bam:
            return get_ADs_bam(file, chrom, pos, REF[0], MQthr, BQthr, insertion=True)
        else:
            return get_ADs_rck(file, chrom, pos, insertion=True)
        #end if
    #end if

    if is_bam:
        return get_ADs_bam(file, chrom, pos, REF[0], MQthr, BQthr)
    else:
        return get_ADs_rck(file, chrom, pos)
    #end if
#end def

#################################################################
#    get_all_ADs
#################################################################
def get_all_ADs(files, chrom, pos, REF, ALT, MQthr, BQthr, is_bam):
    ''' return the AD counts by strand for reference and alternate allele for variant
    in all files '''
    ADfs, ADrs = [], []
    for file in files:
        ADf, ADr = get_ADs_caller(file, chrom, pos, REF, ALT, MQthr, BQthr, is_bam)
        ADfs.append(ADf)
        ADrs.append(ADr)
    #end for

    return np.array(ADfs), np.array(ADrs)
#end def


#################################################################
#    M1_L_calc_aux (original)
#################################################################
def M1_L_calc_aux(rho, k):
    ''' '''
    ALT_count = 1
    M1_L_k = np.zeros(ALT_count + 1)
    default = np.log((1. - rho) / ALT_count)
    M1_L_k = M1_L_k + default
    M1_L_k[k] = np.log(rho)

    return M1_L_k
#end def

#################################################################
#    M1_L_calc (original)
#################################################################
def M1_L_calc(AD, rho):
    ''' '''
    ALT_count = 1
    if AD.size - 1 != ALT_count:
        sys.exit("\nERROR in M1_L_calc\n")
    #end if
    M1_L = []
    for k in range(ALT_count + 1):
        M1_L_k = M1_L_calc_aux(rho,k)
        M1_L.append(M1_L_k)
    #end for

    return M1_L
#end def

#################################################################
#    M2_L_calc_aux (original)
#################################################################
def M2_L_calc_aux(M1_L_k, GT_likelihood_wrt_allele_L):
    ''' '''
    ALT_count = 1
    if M1_L_k.size - 1 != ALT_count:
        sys.exit("\nERROR in M2_L_calc_aux\n")
    #end if
    combos = (ALT_count + 1) * (ALT_count + 2) // 2
    temp_table = GT_likelihood_wrt_allele_L + np.tile(M1_L_k.reshape([1, ALT_count + 1]),[combos, 1])
    M2_L_k = np.zeros(combos)
    for i in range(combos):
        row = temp_table[i, :]
        row_max = np.max(row)
        row = row - row_max
        M2_L_k[i] = np.log(np.sum(np.exp(row))) + row_max
    #end for

    return M2_L_k
#end def

#################################################################
#    M2_L_calc (original)
#################################################################
def M2_L_calc(M1_L, GT_likelihood_wrt_allele_L):
    ''' '''
    ALT_count = 1
    if (M1_L[0]).size - 1 != ALT_count:
        sys.exit("\nERROR in M2_L_calc\n")
    #end if
    M2_L = []
    for k in range(ALT_count + 1):
        M1_L_k = M1_L[k]
        M2_L_k = M2_L_calc_aux(M1_L_k,GT_likelihood_wrt_allele_L)
        M2_L.append(M2_L_k)
    #end for

    return M2_L
#end def

#################################################################
#    GT_marg_L_calc (original)
#################################################################
def GT_marg_L_calc(M2_L_f, M2_L_r, ADf, ADr, prior_L):
    ''' '''
    GT_marg_L = prior_L
    ALT_count = 1
    if ADf.size - 1 != ALT_count or ADr.size - 1 != ALT_count:
        sys.exit("\nERROR in GT_marg_L_calc\n")
    #end if
    for k in range(ALT_count + 1):
        M2_L_k = M2_L_f[k]
        GT_marg_L = GT_marg_L + ADf[k] * M2_L_k
    #end for
    for k in range(ALT_count + 1):
        M2_L_k = M2_L_r[k]
        GT_marg_L = GT_marg_L + ADr[k] * M2_L_k
    #end for

    return GT_marg_L
#end def

#################################################################
#    M3_L_calc_aux (original)
#################################################################
def M3_L_calc_aux(GT_marg_L, M2_L_k):
    ''' '''
    M3_L_k = GT_marg_L - M2_L_k

    return M3_L_k
#end def

#################################################################
#    M3_L_calc (original)
#################################################################
def M3_L_calc(GT_marg_L, M2_L):
    ''' '''
    ALT_count = 1
    if len(M2_L) - 1 != ALT_count:
        sys.exit("\nERROR in M3_L_calc\n")
    #end if
    M3_L = []
    for k in range(ALT_count + 1):
        M2_L_k = M2_L[k]
        M3_L_k = M3_L_calc_aux(GT_marg_L, M2_L_k)
        M3_L.append(M3_L_k)
    #end for

    return M3_L
#end def

#################################################################
#    M4_L_calc_aux (original)
#################################################################
def M4_L_calc_aux(M3_L_k, GT_likelihood_wrt_allele_L):
    ''' '''
    ALT_count = 1
    if (GT_likelihood_wrt_allele_L.shape)[1] - 1 != 1:
        sys.exit("\nERROR in M4_L_calc_aux\n")
    #end if
    combos = (ALT_count + 1) * (ALT_count + 2) // 2
    temp_table = GT_likelihood_wrt_allele_L + np.tile(M3_L_k.reshape([combos, 1]), [1, ALT_count + 1])
    M4_L_k = np.zeros(ALT_count + 1)
    for i in range(ALT_count + 1):
        column = temp_table[:, i]
        column_max = np.max(column)
        column = column - column_max
        M4_L_k[i] = np.log(np.sum(np.exp(column))) + column_max
    #end for

    return M4_L_k
#end def

#################################################################
#    M4_L_calc (original)
#################################################################
def M4_L_calc(M3_L, GT_likelihood_wrt_allele_L):
    ''' '''
    ALT_count = 1
    if (GT_likelihood_wrt_allele_L.shape)[1] - 1 != ALT_count:
        sys.exit("\nERROR in M4_L_calc\n")
    #end if
    M4_L = []
    for k in range(ALT_count + 1):
        M3_L_k = M3_L[k]
        M4_L_k = M4_L_calc_aux(M3_L_k, GT_likelihood_wrt_allele_L)
        M4_L.append(M4_L_k)
    #end for

    return M4_L
#end def

#################################################################
#    A_marg_L_calc (original)
#################################################################
def A_marg_L_calc(M1_L, M4_L):
    ''' '''
    ALT_count = 1
    if len(M1_L) - 1 != ALT_count:
        sys.exit("\nERROR in A_marg_L_calc\n")
    #end if
    A_marg_L = []
    for k in range(ALT_count + 1):
        M1_L_k = M1_L[k]
        M4_L_k = M4_L[k]
        A_marg_L_k = M1_L_k + M4_L_k
        A_marg_L.append(A_marg_L_k)
    #end for

    return A_marg_L
#end def

#################################################################
#    T_term_calc_for_rho (original)
#################################################################
def T_term_calc_for_rho(A_marg_L, AD):
    ''' '''
    if len(A_marg_L) != AD.size:
        sys.exit("\nERROR in T_term_calc\n")
    #end if
    ALT_count = AD.size - 1
    T1_term = 0.
    T2_term = 0.
    for k in range(ALT_count + 1):
        A_marg_L_k = A_marg_L[k]
        A_marg_temp = np.exp(A_marg_L_k - np.max(A_marg_L_k))
        A_marg = A_marg_temp / np.sum(A_marg_temp)
        T1_term = T1_term + A_marg[k] * AD[k]
        T2_term = T2_term + (1. - A_marg[k]) * AD[k]
    #end for

    return T1_term, T2_term
#end def

#################################################################
#    GT_marg_L_to_GT_marg (original)
#################################################################
def GT_marg_L_to_GT_marg(GT_marg_L):
    ''' '''
    M = np.max(GT_marg_L)
    GT_marg_L = GT_marg_L - M
    GT_marg = np.exp(GT_marg_L)
    S = np.sum(GT_marg)
    GT_marg = GT_marg / S
    joint_probty_term = np.log(S) + M

    return GT_marg, joint_probty_term
#end def

#################################################################
#    EM_step (original)
#################################################################
def EM_step(ADf_list, ADr_list, rho_f_old, rho_r_old, prior_L_old, GT_likelihood_wrt_allele_L, a, b, D_original, allele_freq):
    ''' '''
    D = np.zeros(3)
    D[0] = D_original[0]
    D[1] = D_original[1]
    D[2] = D_original[2]

    if allele_freq <= 0.:
        AF = 0.
    else:
        AF = allele_freq
    #end if

    f0 = (1. - AF) ** 2.
    f2 = AF ** 2.
    f1 = 1. - f0 - f2
    D = np.array([f0, f1, f2]) * 1000. + 2.

    T1_f = a - 1.
    T2_f = b - 1.
    T1_r = a - 1.
    T2_r = b - 1.
    T_for_prior = D - 1.
    joint_probty =                (a - 1.) * np.log(rho_f_old) + (b - 1.) * np.log(1. - rho_f_old)
    joint_probty = joint_probty + (a - 1.) * np.log(rho_r_old) + (b - 1.) * np.log(1. - rho_r_old)
    for i in range(3):
        joint_probty = joint_probty + (D[i] - 1) * prior_L_old[i]
    #end for

    if len(ADf_list) != len(ADr_list):
        sys.exit("\nERROR1 in EM_step\n")
    #end if

    for i in range(len(ADf_list)):
        ADf = ADf_list[i]
        ADr = ADr_list[i]
        M1_L_f = M1_L_calc(ADf, rho_f_old)
        M1_L_r = M1_L_calc(ADr, rho_r_old)
        M2_L_f = M2_L_calc(M1_L_f, GT_likelihood_wrt_allele_L)
        M2_L_r = M2_L_calc(M1_L_r, GT_likelihood_wrt_allele_L)
        GT_marg_L = GT_marg_L_calc(M2_L_f, M2_L_r, ADf, ADr, prior_L_old)
        M3_L_f = M3_L_calc(GT_marg_L, M2_L_f)
        M3_L_r = M3_L_calc(GT_marg_L, M2_L_r)
        M4_L_f = M4_L_calc(M3_L_f, GT_likelihood_wrt_allele_L)
        M4_L_r = M4_L_calc(M3_L_r, GT_likelihood_wrt_allele_L)
        A_marg_L_f = A_marg_L_calc(M1_L_f, M4_L_f)
        A_marg_L_r = A_marg_L_calc(M1_L_r, M4_L_r)

        T1_term_f, T2_term_f = T_term_calc_for_rho(A_marg_L_f, ADf)
        T1_term_r, T2_term_r = T_term_calc_for_rho(A_marg_L_r, ADr)

        T1_f = T1_f + T1_term_f
        T2_f = T2_f + T2_term_f
        T1_r = T1_r + T1_term_r
        T2_r = T2_r + T2_term_r

        GT_marg,joint_probty_term = GT_marg_L_to_GT_marg(GT_marg_L)
        joint_probty = joint_probty + joint_probty_term

        T_for_prior = T_for_prior + GT_marg
    #end for

    rho_f_new = 1. / (1. + T2_f / T1_f)
    rho_r_new = 1. / (1. + T2_r / T1_r)
    prior_new = T_for_prior / np.sum(T_for_prior)
    prior_L_new = np.log(prior_new)

    return rho_f_new, rho_r_new, prior_L_new, joint_probty
#end def

#################################################################
#    EM_full (original)
#################################################################
def EM_full(ADfs, ADrs, rho_f_old, rho_r_old, prior_L_old, GT_likelihood_wrt_allele_L, a, b, D, allele_freq):
    ''' '''
    joint_probty_s = []
    joint_probty_new = np.nan
    for i in range(3):
        joint_probty_old = joint_probty_new
        rho_f_new, rho_r_new, prior_L_new, joint_probty_new = \
            EM_step(ADfs, ADrs, rho_f_old, rho_r_old, prior_L_old, GT_likelihood_wrt_allele_L, a, b, D, allele_freq)
        rho_f_old = rho_f_new
        rho_r_old = rho_r_new
        prior_L_old = prior_L_new
        joint_probty_s.append(joint_probty_new)
    #end for

    while np.abs(joint_probty_old - joint_probty_new) > 10 ** -7:
        joint_probty_old = joint_probty_new
        rho_f_new, rho_r_new, prior_L_new, joint_probty_new = \
            EM_step(ADfs, ADrs, rho_f_old, rho_r_old, prior_L_old, GT_likelihood_wrt_allele_L, a, b, D, allele_freq)
        rho_f_old = rho_f_new
        rho_r_old = rho_r_new
        prior_L_old = prior_L_new
        joint_probty_s.append(joint_probty_new)
    #end while

    return rho_f_new, rho_r_new, prior_L_new, joint_probty_s
#end def

#################################################################
#    GTL_L_calc (original)
#################################################################
def GTL_L_calc(ADf, ADr, rho_f, rho_r, GT_likelihood_wrt_allele_L):
    ''' '''
    M1_L_f = M1_L_calc(ADf, rho_f)
    M1_L_r = M1_L_calc(ADr, rho_r)
    M2_L_f = M2_L_calc(M1_L_f, GT_likelihood_wrt_allele_L)
    M2_L_r = M2_L_calc(M1_L_r, GT_likelihood_wrt_allele_L)
    prior_L = np.zeros(3)
    GTL_L = GT_marg_L_calc(M2_L_f, M2_L_r, ADf, ADr, prior_L)
    GTL_L = GTL_L - np.max(GTL_L)

    return GTL_L
#end def

#################################################################
#    posterior_probty_calc_exact (original)
#################################################################
def posterior_probty_calc_exact(prior_L, table_L, C_GL_L, M_GL_L, D_GL_L):
    ''' '''
    combos = 3
    work_column = np.empty(combos ** 2)
    for I1 in range(combos):
        for I2 in range(combos):
            II = I1 * combos + I2
            work_column[II] = prior_L[I1] + prior_L[I2] + M_GL_L[I1] + D_GL_L[I2]
        #end for
    #end for

    work_table = table_L + np.tile(C_GL_L, [combos ** 2, 1]) + np.tile(np.reshape(work_column, [combos ** 2, 1]), [1, combos])
    work_table = work_table - np.max(work_table)
    work_table = np.exp(work_table)
    work_table = work_table / np.sum(work_table)
    PP = np.max(np.array([work_table[0][1], work_table[0][2]]))

    return PP, work_table
#end def

#################################################################
#    denovo_P_calc (original)
#################################################################
def denovo_P_calc(ADfs, ADrs, rho_f, rho_r, GT_likelihood_wrt_allele_L, table_L, prior_L):
    ''' '''
    M_GL_L = GTL_L_calc(ADfs[0], ADrs[0], rho_f, rho_r, GT_likelihood_wrt_allele_L)
    D_GL_L = GTL_L_calc(ADfs[1], ADrs[1], rho_f, rho_r, GT_likelihood_wrt_allele_L)
    C_GL_L = GTL_L_calc(ADfs[2], ADrs[2], rho_f, rho_r, GT_likelihood_wrt_allele_L) # child is the last one
    PP, work_table = posterior_probty_calc_exact(prior_L, table_L, C_GL_L, M_GL_L, D_GL_L)

    return PP, work_table
#end def

#################################################################
#    PP_calc (original)
#################################################################
def PP_calc(trio_files, unrelated_files, chrom, pos, REF, ALT, allele_freq, MQthresh, BQthresh, is_bam):
    ''' '''
    ADfs_U, ADrs_U = get_all_ADs(unrelated_files, chrom, pos, REF, ALT, MQthresh, BQthresh, is_bam)
    rho_f_old, rho_r_old = 0.8, 0.8
    prior_old = np.array([1. / 3, 1. / 3, 1. / 3])
    prior_old = prior_old / np.sum(prior_old)
    prior_L_old = np.log(prior_old)
    GT_likelihood_wrt_allele = GT_likelihood_wrt_allele_calc(1)
    GT_likelihood_wrt_allele_L = np.log(GT_likelihood_wrt_allele)
    a, b, D = 2., 2., np.array([2., 2, 2])

    rho_f_new, rho_r_new, prior_L_new, joint_probty_s = \
        EM_full(ADfs_U, ADrs_U, rho_f_old, rho_r_old, prior_L_old, GT_likelihood_wrt_allele_L, a, b, D, allele_freq)

    AF_unrel = 0.
    for i in range(ADfs_U.shape[0]):
        temp1 = GTL_L_calc(ADfs_U[i], ADrs_U[i], rho_f_new, rho_r_new, GT_likelihood_wrt_allele_L)
        temp = temp1 + prior_L_new
        temp = temp - np.max(temp)
        temp = np.exp(temp)
        temp = temp / np.sum(temp)

        AF_unrel = AF_unrel + temp[1] + temp[2] * 2.
    #end for

    AF_unrel = AF_unrel / 2. / ADfs_U.shape[0]

    ADfs, ADrs = get_all_ADs(trio_files, chrom, pos, REF, ALT, MQthresh, BQthresh, is_bam)

    table = table_gen(1, 1e-8)
    table_L = np.log(table)
    PP, work_table = denovo_P_calc(ADfs, ADrs, rho_f_new, rho_r_new, GT_likelihood_wrt_allele_L, table_L, prior_L_new)

    return PP, ADfs, ADrs, ADfs_U, ADrs_U, rho_f_new, rho_r_new, prior_L_new, AF_unrel
#end def

#################################################################
#    ALT_count_check_parents
#################################################################
def ALT_count_check_parents(ADfs, ADrs, thr):
    ''' check if total alternate reads count in parents is over threshold '''
    if len(ADfs) != 3 or len(ADrs) != 3:
        sys.exit("\nERROR in retrieving stranded AD counts: missing information for trio\n")
    #end if
    alt_count = ADfs[0][1] + ADfs[1][1] + ADrs[0][1] + ADrs[1][1]

    if alt_count > thr: return True
    else: return False
    #end if
#end def

#################################################################
#    runner
#################################################################
def main(args, test=False):
    ''' '''
    # Variables
    is_bam = True if args['bam'] else False
    is_afthr = True if args['afthr'] else False
    afthr, aftag, aftag_idx = 1., 'novoAF', 0 # novoAF as aftag placeholder if not is_afthr
    ppthr = float(args['ppthr']) if args['ppthr'] else 0.
    afthr_unrelated = float(args['afthr_unrelated']) if args['afthr_unrelated'] else 1.
    MQthr = int(args['MQthr']) if args['MQthr'] else 0
    BQthr = int(args['BQthr']) if args['BQthr'] else 0
    ADthr = int(args['ADthr']) if args['ADthr'] else 0
    RSTR_def = '##FORMAT=<ID=RSTR,Number=4,Type=Integer,Description="Read counts by strand for ref and alt alleles (Rf,Af,Rr,Ar)">'
    novoCaller_def = '##INFO=<ID=novoPP,Number=1,Type=Float,Description="Posterior probability from novoCaller">'
    # NA chromosomes set -> import from shared_vars
    if test: NA_chroms = test_NA_chroms
    else: NA_chroms = real_NA_chroms
    #end if
    is_NA = False
    is_verbose = True if args['verbose'] else False

    # Buffers
    fo = open(args['outputfile'], 'w')

    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf(args['inputfile'])

    # Check arguments
    if is_afthr:
        afthr = float(args['afthr'])
        if args['aftag']:
            aftag, aftag_idx = vcf_obj.header.check_tag_definition(args['aftag'])
        else:
            sys.exit('\nERROR in parsing arguments: to filter by population allele frequency please specify the TAG to use\n')
        #end if
    #end if

    # Data structures
    variants_passed = []

    # Getting files and associated IDs
    sys.stderr.write('Getting unrelated and trio files...\n')
    sys.stderr.flush()

    if is_bam: # if bam files
        unrelated_files, IDs_unrelated = buffering_bams(args['unrelatedfiles'])
        trio_files, IDs_trio = buffering_bams(args['triofiles']) # [parent, parent, child]
    else:
        unrelated_files, IDs_unrelated = buffering_rcks(args['unrelatedfiles'])
        trio_files, IDs_trio = buffering_rcks(args['triofiles']) # [parent, parent, child]
    #end if

    # Checking info files for trio is complete
    if len(trio_files) != 3:
        sys.exit('\nERROR in BAMs info file for trio: missing information for some family member\n')
    #end if

    # Checking information for trio is complete in the vcf
    for ID in IDs_trio:
        if ID not in vcf_obj.header.IDs_genotypes:
            sys.exit('\nERROR in VCF file: missing information for some family member\n')
        #end if
    #end for

    # Reading variants
    analyzed = 0
    for i, vnt_obj in enumerate(vcf_obj.parse_variants()):
        if is_verbose:
            sys.stderr.write('\rAnalyzing variant... ' + str(i + 1))
            sys.stderr.flush()
        #end if

        # # Check if chromosome is canonical and in valid format
        # if not check_chrom(vnt_obj.CHROM): # skip variant if not
        #     continue
        # #end if

        # Getting allele frequency from novoAF tag
        af = allele_frequency(vnt_obj, aftag, aftag_idx)

        # is_NA reset
        is_NA = False

        # Calculate statistics
        if af <= afthr: # hard filter on allele frequency
            analyzed += 1
            PP, ADfs, ADrs, ADfs_U, ADrs_U, _, _, _, AF_unrel = \
                PP_calc(trio_files, unrelated_files, vnt_obj.CHROM, vnt_obj.POS, vnt_obj.REF, vnt_obj.ALT, af, MQthr, BQthr, is_bam)
            if vnt_obj.CHROM.replace('chr', '') in NA_chroms:
            # model assumptions does not apply to sex and mithocondrial chromosomes, PP -> NA
                PP = 0.
                is_NA = True
            #end if
            if ADthr and ALT_count_check_parents(ADfs, ADrs, ADthr):
            # AD in parents over ADthr, PP -> 0
                PP = 0.
                is_NA = False
            #end if
            if AF_unrel <= afthr_unrelated and PP >= ppthr: # hard filter on AF_unrel, PP
                variants_passed.append([PP, ADfs, ADrs, ADfs_U, ADrs_U, AF_unrel, is_NA, vnt_obj])
            #end if
        #end if
    #end for

    # Writing output
    sys.stderr.write('\n\n...Writing results for ' + str(analyzed) + ' analyzed variants out of ' + str(i + 1) + ' total variants\n')
    sys.stderr.flush()

    # Header definitions
    is_RSTR = 'RSTR' in vcf_obj.header.definitions
    is_novoCaller = 'novoPP' in vcf_obj.header.definitions

    if not is_RSTR:
        vcf_obj.header.add_tag_definition(RSTR_def, 'FORMAT')
    #end if
    if not is_novoCaller:
        vcf_obj.header.add_tag_definition(novoCaller_def, 'INFO')
    #end if
    vcf_obj.write_definitions(fo)

    # Adding to header columns unrelated samples missing IDs
    fo.write(vcf_obj.header.columns.rstrip())
    for ID in IDs_unrelated:
        if ID not in vcf_obj.header.IDs_genotypes:
            fo.write('\t' + ID)
        #end if
    #end for
    fo.write('\n')

    # Variants passed
    for variant in sorted(variants_passed, key=lambda x: x[0], reverse=True):
        PP, ADfs, ADrs, ADfs_U, ADrs_U, AF_unrel, is_NA, vnt_obj = variant

        # Removing older tags fields if present
        if is_RSTR:
            vnt_obj.remove_tag_genotype('RSTR')
        #end if
        if is_novoCaller:
            vnt_obj.remove_tag_info('novoPP')
        #end if

        # Adding new tag
        if not is_NA:
            vnt_obj.add_tag_info('novoPP={0}'.format(PP))
        #end if

        # Fill the trailing fields dropped in genotypes
        vnt_obj.complete_genotype()

        # Updating genotypes trio
        for i, ID in enumerate(IDs_trio):
            values = '{0},{1},{2},{3}'.format(int(ADfs[i][0]), int(ADfs[i][1]), int(ADrs[i][0]), int(ADrs[i][1]))
            vnt_obj.add_values_genotype(ID, values)
        #end for

        # Updating genotypes unrelated
        unrelated_genotypes = []
        for i, ID in enumerate(IDs_unrelated):
            values = '{0},{1},{2},{3}'.format(int(ADfs_U[i][0]), int(ADfs_U[i][1]), int(ADrs_U[i][0]), int(ADrs_U[i][1]))
            if ID in vnt_obj.GENOTYPES:
                vnt_obj.add_values_genotype(ID, values)
            else:
                unrelated_genotypes.append(vnt_obj.empty_genotype() + ':' + values)
            #end if
        #end for

        # Updating FORMAT
        vnt_obj.add_tag_format('RSTR')

        # Writing output
        if unrelated_genotypes:
            fo.write(vnt_obj.to_string().rstrip() + '\t' + '\t'.join(unrelated_genotypes) + '\n')
        else:
            vcf_obj.write_variant(fo, vnt_obj)
        #end if
    #end for

    # Closing files buffers
    fo.close()
    if is_bam:
        for buffer in unrelated_files:
            buffer.close()
        #end for
        for buffer in trio_files:
            buffer.close()
        #end for
    #end if
#end def


#################################################################
#
#    MAIN
#
#################################################################
if __name__ == "__main__":

    main()

#end if
