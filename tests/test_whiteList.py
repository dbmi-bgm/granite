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
            'SpliceAI': None, 'CLINVAR': None, 'CLINVARonly': None, 'VEP': None, 'VEPrescue': None, 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
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
            'SpliceAI': None, 'CLINVAR': True, 'CLINVARonly': None, 'VEP': None, 'VEPrescue': None, 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_CLINVAR.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_whiteList_CLINVAR_as_ANNOT():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList_as_ANNOT.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'CLINVAR': True, 'CLINVARonly': None, 'VEP': None, 'VEPrescue': None, 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': 'ANNOT'}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_CLINVAR_as_ANNOT.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_whiteList_CLINVARonly():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList_CLINVARonly.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'CLINVAR': True, 'CLINVARonly': ['Whatever', 'Pathogenic'], 'VEP': None, 'VEPrescue': None, 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_CLINVARonly.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_whiteList_CLINVARonly_bis():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList_CLINVARonly.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'CLINVAR': True, 'CLINVARonly': ['Pathogenic'], 'VEP': None, 'VEPrescue': None, 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_CLINVARonly_bis.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_SpliceAI_001():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': '0.01', 'CLINVAR': None, 'CLINVARonly': None, 'VEP': None, 'VEPrescue': None, 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
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
            'SpliceAI': '0.02', 'CLINVAR': None, 'CLINVARonly': None, 'VEP': None, 'VEPrescue': None, 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
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
            'SpliceAI': '0.08', 'CLINVAR': None, 'CLINVARonly': None, 'VEP': None, 'VEPrescue': None, 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
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
            'SpliceAI': '0.081', 'CLINVAR': None, 'CLINVARonly': None, 'VEP': None, 'VEPrescue': None, 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
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
            'SpliceAI': '0.08', 'CLINVAR': True, 'CLINVARonly': None, 'VEP': None, 'VEPrescue': None, 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
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
            'SpliceAI': '0.081', 'CLINVAR': True, 'CLINVARonly': None, 'VEP': None, 'VEPrescue': None, 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
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
            'SpliceAI': None, 'CLINVAR': None, 'CLINVARonly': None, 'VEP': True, 'VEPrescue': ['TF_binding_site_variant'], 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_VEP.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_whiteList_VEP_as_CSQ():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList_VEP_as_CSQ.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'CLINVAR': None, 'CLINVARonly': None, 'VEP': True, 'VEPrescue': ['TF_binding_site_variant'], 'VEPremove': None, 'BEDfile': None,
            'VEPtag': 'CSQ', 'CLINVARtag': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_VEP_as_CSQ.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_whiteList_VEP_splice_region_variant():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList_VEP.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'CLINVAR': None, 'CLINVARonly': None, 'VEP': True, 'VEPrescue': ['TF_binding_site_variant', 'splice_region_variant'], 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
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
            'SpliceAI': None, 'CLINVAR': None, 'CLINVARonly': None, 'VEP': True, 'VEPrescue': ['TF_binding_site_variant', 'non_coding_transcript_variant'], 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
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
            'SpliceAI': None, 'CLINVAR': None, 'CLINVARonly': None, 'VEP': True, 'VEPrescue': ['TF_binding_site_variant', 'splice_region_variant', 'non_coding_transcript_variant'], 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
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
            'SpliceAI': None, 'CLINVAR': True, 'CLINVARonly': None, 'VEP': True, 'VEPrescue': None, 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_VEP_CLINVAR.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_whiteList_microannot_VEP_CLINVAR_SPLICEAI():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList_microannot.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': '0.8', 'CLINVAR': True, 'CLINVARonly': None, 'VEP': True, 'VEPrescue': ['TF_binding_site_variant', 'non_coding_transcript_exon_variant'], 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_microannot_VEP_CLINVAR_SPLICEAI.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_whiteList_microannot_VEP_CLINVAR():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList_microannot.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'CLINVAR': True, 'CLINVARonly': None, 'VEP': True, 'VEPrescue': ['TF_binding_site_variant', 'non_coding_transcript_exon_variant'], 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_microannot_VEP_CLINVAR.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_whiteList_microannot_VEP_SPLICEAI():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList_microannot.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': '0.8', 'CLINVAR': None, 'CLINVARonly': None, 'VEP': True, 'VEPrescue': ['TF_binding_site_variant', 'non_coding_transcript_exon_variant'], 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_microannot_VEP_SPLICEAI.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_whiteList_microannot_VEP():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList_microannot.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'CLINVAR': None, 'CLINVARonly': None, 'VEP': True, 'VEPrescue': ['TF_binding_site_variant', 'non_coding_transcript_exon_variant'], 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_microannot_VEP.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_whiteList_microannot_VEP_remove():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList_microannot.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'CLINVAR': None, 'CLINVARonly': None, 'VEP': True, 'VEPrescue': ['TF_binding_site_variant'], 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_microannot_VEP_remove.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_whiteList_microannot_VEP_remove_bis():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList_microannot.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'CLINVAR': None, 'CLINVARonly': None, 'VEP': True, 'VEPrescue': ['TF_binding_site_variant', 'non_coding_transcript_exon_variant'], 'VEPremove': ['mature_miRNA_variant'], 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_microannot_VEP_remove_bis.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_whiteList_microannot_VEP_remove_save():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList_microannot.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'CLINVAR': None, 'CLINVARonly': None, 'VEP': True, 'VEPrescue': ['splice_region_variant', 'non_coding_transcript_exon_variant', 'TF_binding_site_variant'], 'VEPremove': ['mature_miRNA_variant'], 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_microannot_VEP_remove_save.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_whiteList_microannot_VEP_remove_double():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList_microannot.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'CLINVAR': None, 'CLINVARonly': None, 'VEP': True, 'VEPrescue': ['non_coding_transcript_exon_variant'], 'VEPremove': ['mature_miRNA_variant'], 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_microannot_VEP_remove_double.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_whiteList_microannot_VEP_remove_double_BED():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList_microannot_BED.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'CLINVAR': None, 'CLINVARonly': None, 'VEP': True, 'VEPrescue': ['non_coding_transcript_exon_variant'], 'VEPremove': ['mature_miRNA_variant'],
            'BEDfile': 'tests/files/input_BED_whiteList.bed', 'VEPtag': None, 'CLINVARtag': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_microannot_VEP_remove_double_BED.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_whiteList_VEP_missing():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList_VEP_missing.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'CLINVAR': True, 'CLINVARonly': None, 'VEP': True, 'VEPrescue': None, 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
    # Run
    main_whiteList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_whiteList_VEP_missing.out')]
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
            'SpliceAI': None, 'CLINVAR': None, 'CLINVARonly': None, 'VEP': None, 'VEPrescue': ['splice_region_variant', 'non_coding_transcript_variant'], 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
    # Run and Tests
    with pytest.raises(SystemExit) as e:
        assert main_whiteList(args)
    assert str(e.value) == '\nERROR in parsing arguments: specify the flag "--VEP" to filter by VEP annotations to apply rescue terms or remove additional terms\n'
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_args_CLINVAR_conflict():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_whiteList.vcf', 'outputfile': 'tests/files/main_test.out',
            'SpliceAI': None, 'CLINVAR': None, 'CLINVARonly': ['Pathogenic'], 'VEP': None, 'VEPrescue': None, 'VEPremove': None, 'BEDfile': None,
            'VEPtag': None, 'CLINVARtag': None}
    # Run and Tests
    with pytest.raises(SystemExit) as e:
        assert main_whiteList(args)
    assert str(e.value) == '\nERROR in parsing arguments: specify the flag "--CLINVAR" to filter by CLINVAR annotations to specify tags or keywords to whitelist\n'
    # Clean
    os.remove('tests/files/main_test.out')
#end def
