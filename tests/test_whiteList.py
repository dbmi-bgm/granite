#################################################################
#   Libraries
#################################################################
import sys, os
import pytest
from granite.whiteList import (
                            main as main_whiteList
                            )


#################################################################
#   Tests
#################################################################
def test_run_whiteList():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'CLINVAR': None, 'VEP': None, 'VEPrescue': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_whiteList_CLINVAR():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'CLINVAR': True, 'VEP': None, 'VEPrescue': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_CLINVAR.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_SpliceAI_001():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': '0.01', 'CLINVAR': None, 'VEP': None, 'VEPrescue': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_SpliceAI_0-01.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_SpliceAI_002():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': '0.02', 'CLINVAR': None, 'VEP': None, 'VEPrescue': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_SpliceAI_0-02.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_SpliceAI_008():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': '0.08', 'CLINVAR': None, 'VEP': None, 'VEPrescue': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_SpliceAI_0-08.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_SpliceAI_0081():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': '0.081', 'CLINVAR': None, 'VEP': None, 'VEPrescue': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_SpliceAI_008_CLINVAR():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': '0.08', 'CLINVAR': True, 'VEP': None, 'VEPrescue': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_SpliceAI_0-08_CLINVAR.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_SpliceAI_0081_CLINVAR():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': '0.081', 'CLINVAR': True, 'VEP': None, 'VEPrescue': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_CLINVAR.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_whiteList_VEP():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList_VEP.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'CLINVAR': None, 'VEP': True, 'VEPrescue': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_VEP.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_whiteList_VEP_splice_region_variant():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList_VEP.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'CLINVAR': None, 'VEP': True, 'VEPrescue': ['splice_region_variant']}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_VEP_splice_region_variant.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_whiteList_VEP_non_coding_transcript_variant():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList_VEP.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'CLINVAR': None, 'VEP': True, 'VEPrescue': ['non_coding_transcript_variant']}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_VEP_non_coding_transcript_variant.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_whiteList_VEP_rescue():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList_VEP.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'CLINVAR': None, 'VEP': True, 'VEPrescue': ['splice_region_variant', 'non_coding_transcript_variant']}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_VEP_rescue.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_whiteList_VEP_CLINVAR():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList_VEP.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'CLINVAR': True, 'VEP': True, 'VEPrescue': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_VEP_CLINVAR.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def


#################################################################
#   Errors
#################################################################
def test_args_VEP_conflict():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'CLINVAR': None, 'VEP': None, 'VEPrescue': ['splice_region_variant', 'non_coding_transcript_variant']}
    # Run and Tests
    with pytest.raises(SystemExit) as e:
        assert main_whiteList(args)
    assert str(e.value) == 'ERROR in parsing arguments: specify the flag "--VEP" to filter by VEP annotations to apply rescue terms\n'
#end def

def test_missing_VEP():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList_VEP_missing.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'CLINVAR': None, 'VEP': True, 'VEPrescue': None}
    # Run and Tests
    with pytest.raises(SystemExit) as e:
        assert main_whiteList(args)
    assert '942451' in str(e.value)
    # Clean
    os.remove('tests/files/main_test.out')
#end def
