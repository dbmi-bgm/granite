#################################################################
#   Libraries
#################################################################
import sys, os
import pytest
from granite.validateVCF import (
                            main as main_validateVCF
                            )

#################################################################
#   Functions
#################################################################
def _do_cmp(f1, f2):
    with open(f1, 'rb') as fp1, open(f2, 'rb') as fp2:
        while True:
            b1 = fp1.read()
            b2 = fp2.read()
            if b1 != b2:
                return False
            #end if
            if not b1:
                return True
            #end if
        #end while
    #end with
#end def

#################################################################
#   Tests
#################################################################
# def test_run_validateVCF_het():
#     ''' '''
#     # Variables
#     args = {'inputfile': 'tests/files/GAPFIR4IUBTH.merged.vcf.gz', 'outputfile': 'tests/files/main_test.out',
#             'pedigree': ['tests/files/pedigree.UGRP1_ind_5.json', 'tests/files/pedigree.UGRP1_ind_6.json'], 'type': None,
#             'anchor': ['UGRP1_ind_5_GAPFIR4IUBTH', 'UGRP1_ind_6_GAPFIR4IUBTH'], 'novo': None,
#             'het': ['UGRP1_ind_5_GAPFIR4IUBTH', 'UGRP1_ind_6_GAPFIR4IUBTH'], 'verbose': None}
#     # Run
#     main_validateVCF(args)
#     # Tests
#     assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/GAPFIR4IUBTH.het.json')]
#     f1 = 'autosomal_heterozygous_accuracy_family_UGRP1_ind_5_GAPFIR4IUBTH-UGRP1_ind_6_GAPFIR4IUBTH.png'
#     f2 = 'tests/files/autosomal_heterozygous_accuracy_family_UGRP1_ind_5_GAPFIR4IUBTH-UGRP1_ind_6_GAPFIR4IUBTH.out'
#     assert _do_cmp(f1, f2)
#     f1 = 'autosomal_heterozygous_distribution_family_UGRP1_ind_5_GAPFIR4IUBTH-UGRP1_ind_6_GAPFIR4IUBTH.png'
#     f2 = 'tests/files/autosomal_heterozygous_distribution_family_UGRP1_ind_5_GAPFIR4IUBTH-UGRP1_ind_6_GAPFIR4IUBTH.out'
#     assert _do_cmp(f1, f2)
#     f1 = 'autosomal_heterozygous_distribution_trio_UGRP1_ind_5_GAPFIR4IUBTH-UGRP1_ind_6_GAPFIR4IUBTH.png'
#     f2 = 'tests/files/autosomal_heterozygous_distribution_trio_UGRP1_ind_5_GAPFIR4IUBTH-UGRP1_ind_6_GAPFIR4IUBTH.out'
#     assert _do_cmp(f1, f2)
#     # Clean
#     os.remove('tests/files/main_test.out')
#     os.remove('autosomal_heterozygous_accuracy_family_UGRP1_ind_5_GAPFIR4IUBTH-UGRP1_ind_6_GAPFIR4IUBTH.png')
#     os.remove('autosomal_heterozygous_distribution_family_UGRP1_ind_5_GAPFIR4IUBTH-UGRP1_ind_6_GAPFIR4IUBTH.png')
#     os.remove('autosomal_heterozygous_distribution_trio_UGRP1_ind_5_GAPFIR4IUBTH-UGRP1_ind_6_GAPFIR4IUBTH.png')
# #end def
