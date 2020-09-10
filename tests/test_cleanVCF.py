#################################################################
#   Libraries
#################################################################
import sys, os
import pytest
from granite.cleanVCF import (
                            main as main_cleanVCF
                            )


#################################################################
#   Tests
#################################################################
def test_run_cleanVCF_void():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_cleanVCF.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'VEP': None, 'VEPrescue': None, 'VEPremove': None,
            'VEPtag': 'VEP', 'tag': [], 'verbose': None, 'VEPsep': None, 'SpliceAItag': None}
    # Run
    main_cleanVCF(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_cleanVCF.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_cleanVCF_VEP():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_cleanVCF.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'VEP': True, 'VEPrescue': None, 'VEPremove': None,
            'VEPtag': 'VEP', 'tag': [], 'verbose': None, 'VEPsep': None, 'SpliceAItag': None}
    # Run
    main_cleanVCF(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_cleanVCF_VEP.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_cleanVCF_VEP_as_ANNOTTAG():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_cleanVCF_ANNOTTAG.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'VEP': True, 'VEPrescue': None, 'VEPremove': None,
            'VEPtag': 'ANNOTTAG', 'tag': [], 'verbose': None, 'VEPsep': None, 'SpliceAItag': None}
    # Run
    main_cleanVCF(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_cleanVCF_ANNOTTAG_VEP.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_cleanVCF_VEP_SpAI():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_cleanVCF_SpAI.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': 0.9, 'VEP': True, 'VEPrescue': None, 'VEPremove': None,
            'VEPtag': 'VEP', 'tag': [], 'verbose': None, 'VEPsep': None, 'SpliceAItag': None}
    # Run
    main_cleanVCF(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_cleanVCF_SpAI.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_cleanVCF_VEP_rescue():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_cleanVCF.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'VEP': True, 'VEPrescue': ['intron_variant'], 'VEPremove': None,
            'VEPtag': 'VEP', 'tag': [], 'verbose': None, 'VEPsep': None, 'SpliceAItag': None}
    # Run
    main_cleanVCF(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_cleanVCF_VEP_rescue.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_cleanVCF_VEP_SpAI_remove_rescue_tag_bis():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_cleanVCF_SpAI.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': 0.9, 'VEP': True, 'VEPrescue': ['downstream_gene_variant', 'NMD_transcript_variant'],
            'VEPremove': None, 'VEPtag': 'VEP', 'tag': ['BaseQRankSum'], 'verbose': None, 'VEPsep': None, 'SpliceAItag': None}
    # Run
    main_cleanVCF(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_cleanVCF_VEP_SpAI_remove_rescue_tag_bis.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_cleanVCF_VEP_SpAI_remove_rescue_tag_bis_sep():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_cleanVCF_SpAI_sep.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': 0.9, 'VEP': True, 'VEPrescue': ['downstream_gene_variant', 'NMD_transcript_variant'],
            'VEPremove': None, 'VEPtag': 'VEP', 'tag': ['BaseQRankSum'], 'verbose': None, 'VEPsep': '~', 'SpliceAItag': None}
    # Run
    main_cleanVCF(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_cleanVCF_VEP_SpAI_remove_rescue_tag_bis_sep.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_cleanVCF_VEP_remove():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_cleanVCF.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'VEP': True, 'VEPrescue': None, 'VEPremove': ['3_prime_UTR_variant'],
            'VEPtag': 'VEP', 'tag': [], 'verbose': None, 'VEPsep': None, 'SpliceAItag': None}
    # Run
    main_cleanVCF(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_cleanVCF_VEP_remove.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_cleanVCF_VEP_remove_rescue():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_cleanVCF.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'VEP': True, 'VEPrescue': ['downstream_gene_variant', 'NMD_transcript_variant'],
            'VEPremove': ['5_prime_UTR_variant', '3_prime_UTR_variant'], 'VEPtag': 'VEP', 'tag': [], 'verbose': None, 'VEPsep': None, 'SpliceAItag': None}
    # Run
    main_cleanVCF(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_cleanVCF_VEP_remove_rescue.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_cleanVCF_VEP_SpAI_remove_rescue():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_cleanVCF_SpAI.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': 0.9, 'VEP': True, 'VEPrescue': ['downstream_gene_variant', 'NMD_transcript_variant'],
            'VEPremove': ['5_prime_UTR_variant', '3_prime_UTR_variant'], 'VEPtag': 'VEP', 'tag': [], 'verbose': None, 'VEPsep': None, 'SpliceAItag': None}
    # Run
    main_cleanVCF(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_cleanVCF_VEP_SpAI_remove_rescue.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_cleanVCF_VEP_SpAI_remove_rescue_tag():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_cleanVCF_SpAI.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': 0.9, 'VEP': True, 'VEPrescue': ['downstream_gene_variant', 'NMD_transcript_variant'],
            'VEPremove': ['5_prime_UTR_variant', '3_prime_UTR_variant'], 'VEPtag': 'VEP', 'tag': ['BaseQRankSum'], 'verbose': None, 'VEPsep': None, 'SpliceAItag': None}
    # Run
    main_cleanVCF(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_cleanVCF_VEP_SpAI_remove_rescue_tag.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_cleanVCF_VEP_SpAI_remove_rescue_tag_2():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_cleanVCF_SpAI.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': 0.9, 'VEP': True, 'VEPrescue': ['downstream_gene_variant', 'NMD_transcript_variant'],
            'VEPremove': ['5_prime_UTR_variant', '3_prime_UTR_variant'], 'VEPtag': 'VEP', 'tag': ['BaseQRankSum', 'novoCaller'], 'verbose': None, 'VEPsep': None, 'SpliceAItag': None}
    # Run
    main_cleanVCF(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_cleanVCF_VEP_SpAI_remove_rescue_tag_2.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_cleanVCF_VEP_SpAI_SpAItag_remove_rescue_tag_2():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_cleanVCF_SpAI_SpAItag.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': 0.9, 'VEP': True, 'VEPrescue': ['downstream_gene_variant', 'NMD_transcript_variant'],
            'VEPremove': ['5_prime_UTR_variant', '3_prime_UTR_variant'], 'VEPtag': 'VEP', 'tag': ['BaseQRankSum', 'novoCaller'], 'verbose': None, 'VEPsep': None, 'SpliceAItag': 'VALUE_MAX'}
    # Run
    main_cleanVCF(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_cleanVCF_VEP_SpAI_SpAItag_remove_rescue_tag_2.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

#################################################################
#   Errors
#################################################################
def test_args_conflict():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_cleanVCF.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'VEP': None, 'VEPrescue': ['splice_region_variant', 'non_coding_transcript_variant'], 'VEPremove': None,
            'VEPtag': 'VEP', 'tag': [], 'verbose': None, 'VEPsep': None, 'SpliceAItag': None}
    # Run and Tests
    with pytest.raises(SystemExit) as e:
        assert main_cleanVCF(args)
    assert str(e.value) == '\nERROR in parsing arguments: specify the flag "--VEP" to filter by VEP annotations to apply rescue terms or remove additional terms\n'
    # Clean
    os.remove('tests/files/main_test.out')
#end def
