#################################################################
#   Libraries
#################################################################
import sys, os
import pytest
from granite.novoCaller import (
                            main as main_novoCaller
                            )


#################################################################
#   Tests
#################################################################
# this tests may fail because of differences in float handling
# e.g. novoCaller=1.4030838153054451e-05 Vs novoCaller=1.403083815305445e-05
# depending on the machine
def test_run_novoCaller_rck():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_RCK.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_rck.tsv', 'triofiles':'tests/files/trio_rck.tsv',
            'ppthr': None, 'afthr': None, 'aftag': None, 'bam': None,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': None}
    # Run
    main_novoCaller(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_novoCaller_RCK.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_noUNR():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_noUNR.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': '0.01', 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': None}
    # Run
    main_novoCaller(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_novoCaller_BAM_noUNR.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_noUNR_noafthr():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_noUNR.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': None, 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': None}
    # Run
    main_novoCaller(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_novoCaller_BAM_noUNR_noafthr.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_withUNR():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_withUNR.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': '0.01', 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': None}
    # Run
    main_novoCaller(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_novoCaller_BAM_withUNR.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_withUNR_asSibling_ALLBAMS():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_withUNR_asSibling.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': '0.01', 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': None}
    # Run
    main_novoCaller(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_novoCaller_BAM_withUNR_asSibling_ALLBAMS.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_withUNR_asSibling_missingSibling():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_withUNR_asSibling.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam_PSC-01-003_asSibling.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': '0.01', 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': None}
    # Run
    main_novoCaller(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_novoCaller_BAM_withUNR_asSibling_missingSibling.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_withUNR_asSibling_SWAP_SON_missingSibling():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_withUNR_asSibling_SWAP_SON.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam_PSC-01-003_asSibling.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': '0.01', 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': None}
    # Run
    main_novoCaller(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_novoCaller_BAM_withUNR_asSibling_SWAP_SON_missingSibling.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_rerun_noUNR():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_noUNR.out', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': '0.01', 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': None}
    # Run
    main_novoCaller(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_novoCaller_BAM_noUNR.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_rerun_withUNR_asSibling_missingSibling():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_withUNR_asSibling_missingSibling.out', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam_PSC-01-003_asSibling.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': '0.01', 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': None}
    # Run
    main_novoCaller(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_novoCaller_BAM_withUNR_asSibling_missingSibling.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_different_tag():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_noUNR_AlleleFreq.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': '0.01', 'aftag': 'AlleleFreq', 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': None}
    # Run
    main_novoCaller(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_novoCaller_BAM_noUNR_AlleleFreq.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_jc50_wgenome_plus_indels():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': '0.01', 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': None}
    # Run
    main_novoCaller(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_jc50_wgenome_plus_indels_11():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': None, 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': None}
    # Run
    main_novoCaller(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels_11.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_jc50_wgenome_plus_indels_9():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': '7.579382496057039e-05', 'afthr': None, 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': None}
    # Run
    main_novoCaller(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels_9.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_jc50_wgenome_plus_indels_8():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': '0.9999998', 'afthr': None, 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': None}
    # Run
    main_novoCaller(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels_8.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_jc50_wgenome_plus_indels_7():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': '0.9999999', 'afthr': None, 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': None}
    # Run
    main_novoCaller(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels_7.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_jc50_wgenome_plus_indels_4():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': '0.999999999999', 'afthr': None, 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': None}
    # Run
    main_novoCaller(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels_4.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_jc50_wgenome_plus_indels_2():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': '1', 'afthr': None, 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': None}
    # Run
    main_novoCaller(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels_2.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def


#################################################################
#   Errors
#################################################################
# def test_run_novoCaller_bam_annot_missing_novoAF():
#     ''' '''
#     # Variables
#     args = {'inputfile': 'tests/files/input_novoCaller_BAM_annot_missing_novoAF.vcf', 'outputfile': 'tests/files/main_test.out',
#             'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
#             'ppthr': None, 'afthr': '0.01', 'aftag': None, 'bam': True,
#             'MQthr': None, 'BQthr': None}
#     # Run and Tests
#     with pytest.raises(SystemExit) as e:
#         assert main_novoCaller(args)
#     assert '\nERROR in variant parsing: novoAF allele frequency tag in INFO field is missing for variant:\n' in str(e.value)
#     assert '16805'
#     # Clean
#     os.remove('tests/files/main_test.out')
# #end def
#
# def test_run_novoCaller_bam_annot_wrong_novoAF():
#     ''' '''
#     # Variables
#     args = {'inputfile': 'tests/files/input_novoCaller_BAM_annot_wrong_novoAF.vcf', 'outputfile': 'tests/files/main_test.out',
#             'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
#             'ppthr': None, 'afthr': '0.01', 'aftag': None, 'bam': True,
#             'MQthr': None, 'BQthr': None}
#     # Run and Tests
#     with pytest.raises(SystemExit) as e:
#         assert main_novoCaller(args)
#     assert '\nERROR in variant parsing: novoAF allele frequency tag in INFO field is in the wrong format for variant:\n' in str(e.value)
#     assert '16805'
#     # Clean
#     os.remove('tests/files/main_test.out')
# #end def
#
# def test_run_novoCaller_bam_missing_aftag():
#     ''' '''
#     # Variables
#     args = {'inputfile': 'tests/files/input_novoCaller_BAM_annot_missing_novoAF.vcf', 'outputfile': 'tests/files/main_test.out',
#             'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
#             'ppthr': None, 'afthr': '0.01', 'aftag': None, 'bam': True,
#             'MQthr': None, 'BQthr': None}
#     # Run and Tests
#     with pytest.raises(SystemExit) as e:
#         assert main_novoCaller(args)
#     assert '\nERROR in variant parsing: novoAF allele frequency tag in INFO field is missing for variant:\n' in str(e.value)
#     assert '14699'
#     # Clean
#     os.remove('tests/files/main_test.out')
# #end def

def test_run_novoCaller_bam_missing_trio():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_annot_missing_novoAF.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam_missingOne.tsv',
            'ppthr': None, 'afthr': '0.01', 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': None}
    # Run and Tests
    with pytest.raises(SystemExit) as e:
        assert main_novoCaller(args)
    assert '\nERROR in BAMs info file for trio: missing information for some family member\n' == str(e.value)
    # Clean
    os.remove('tests/files/main_test.out')
#end def
