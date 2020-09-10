#################################################################
#   Libraries
#################################################################
import sys, os
import pytest
from granite.comHet import (
                            main as main_comHet
                            )


#################################################################
#   Tests
#################################################################
def test_run_comHet_proband(): # NA12879
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_comHet.vcf', 'outputfile': 'tests/files/main_test.out',
            'trio': ['NA12879_sample'], 'VEPtag': 'VEP', 'allow_undef': None, 'filter_cmpHet': True, 'verbose': None, 'sep': None,
            'impact': None, 'SpliceAItag': None}
    # Run
    main_comHet(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_comHet.out')]
    # Clean
    os.remove('tests/files/main_test.out')
    os.remove('tests/files/main_test.out.summary')
    os.remove('tests/files/main_test.out.json')
#end def

def test_run_comHet_proband_plus_NA12877():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_comHet.vcf', 'outputfile': 'tests/files/main_test.out',
            'trio': ['NA12879_sample', 'NA12877_sample'], 'VEPtag': 'VEP', 'allow_undef': None, 'filter_cmpHet': True, 'verbose': None, 'sep': None,
            'impact': None, 'SpliceAItag': None}
    # Run
    main_comHet(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_comHet_plus_NA12877.out')]
    # Clean
    os.remove('tests/files/main_test.out')
    os.remove('tests/files/main_test.out.summary')
    os.remove('tests/files/main_test.out.json')
#end def

def test_run_comHet_proband_plus_NA12878():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_comHet.vcf', 'outputfile': 'tests/files/main_test.out',
            'trio': ['NA12879_sample', 'NA12878_sample'], 'VEPtag': 'VEP', 'allow_undef': None, 'filter_cmpHet': True, 'verbose': None, 'sep': None,
            'impact': None, 'SpliceAItag': None}
    # Run
    main_comHet(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_comHet_plus_NA12878.out')]
    # Clean
    os.remove('tests/files/main_test.out')
    os.remove('tests/files/main_test.out.summary')
    os.remove('tests/files/main_test.out.json')
#end def

def test_run_comHet_proband_plus_NA12877_NA12878():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_comHet.vcf', 'outputfile': 'tests/files/main_test.out',
            'trio': ['NA12879_sample', 'NA12877_sample', 'NA12878_sample'], 'VEPtag': 'VEP', 'allow_undef': None, 'filter_cmpHet': True, 'verbose': None, 'sep': None,
            'impact': None, 'SpliceAItag': None}
    # Run
    main_comHet(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_comHet_plus_NA12877_NA12878.out')]
    # Clean
    os.remove('tests/files/main_test.out')
    os.remove('tests/files/main_test.out.summary')
    os.remove('tests/files/main_test.out.json')
#end def

def test_run_comHet_proband_plus_NA12877_NA12878_GT():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_comHet_GT.vcf', 'outputfile': 'tests/files/main_test.out',
            'trio': ['NA12879_sample', 'NA12877_sample', 'NA12878_sample'], 'VEPtag': 'VEP', 'allow_undef': None, 'filter_cmpHet': True, 'verbose': None, 'sep': None,
            'impact': None, 'SpliceAItag': None}
    # Run
    main_comHet(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_comHet_GT_plus_NA12877_NA12878.out')]
    # Clean
    os.remove('tests/files/main_test.out')
    os.remove('tests/files/main_test.out.summary')
    os.remove('tests/files/main_test.out.json')
#end def

def test_run_comHet_proband_VEP_as_CSQ_plus_NA12877_NA12878():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_comHet_VEP_as_CSQ.vcf', 'outputfile': 'tests/files/main_test.out',
            'trio': ['NA12879_sample', 'NA12877_sample', 'NA12878_sample'], 'VEPtag': None, 'allow_undef': None, 'filter_cmpHet': True, 'verbose': None, 'sep': None,
            'impact': None, 'SpliceAItag': None}
    # Run
    main_comHet(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_comHet_VEP_as_CSQ_plus_NA12877_NA12878.out')]
    # Clean
    os.remove('tests/files/main_test.out')
    os.remove('tests/files/main_test.out.summary')
    os.remove('tests/files/main_test.out.json')
#end def

def test_run_comHet_proband_plus_NA12877_allow_undef():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_comHet.vcf', 'outputfile': 'tests/files/main_test.out',
            'trio': ['NA12879_sample', 'NA12877_sample'], 'VEPtag': 'VEP', 'allow_undef': True, 'filter_cmpHet': True, 'verbose': None, 'sep': None,
            'impact': None, 'SpliceAItag': None}
    # Run
    main_comHet(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_comHet_plus_NA12877_allow_undef.out')]
    # Clean
    os.remove('tests/files/main_test.out')
    os.remove('tests/files/main_test.out.summary')
    os.remove('tests/files/main_test.out.json')
#end def

def test_run_comHet_proband_plus_NA12878_allow_undef():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_comHet.vcf', 'outputfile': 'tests/files/main_test.out',
            'trio': ['NA12879_sample', 'NA12878_sample'], 'VEPtag': 'VEP', 'allow_undef': True, 'filter_cmpHet': True, 'verbose': None, 'sep': None,
            'impact': None, 'SpliceAItag': None}
    # Run
    main_comHet(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_comHet_plus_NA12878_allow_undef.out')]
    # Clean
    os.remove('tests/files/main_test.out')
    os.remove('tests/files/main_test.out.summary')
    os.remove('tests/files/main_test.out.json')
#end def

def test_run_comHet_proband_plus_NA12878_NA12877_allow_undef():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_comHet.vcf', 'outputfile': 'tests/files/main_test.out',
            'trio': ['NA12879_sample', 'NA12878_sample', 'NA12877_sample'], 'VEPtag': 'VEP', 'allow_undef': True, 'filter_cmpHet': True, 'verbose': None, 'sep': None,
            'impact': None, 'SpliceAItag': None}
    # Run
    main_comHet(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_comHet_plus_NA12878_NA12877_allow_undef.out')]
    # Clean
    os.remove('tests/files/main_test.out')
    os.remove('tests/files/main_test.out.summary')
    os.remove('tests/files/main_test.out.json')
#end def

def test_run_comHet_NA12878():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_comHet.vcf', 'outputfile': 'tests/files/main_test.out',
            'trio': ['NA12878_sample'], 'VEPtag': 'VEP', 'allow_undef': None, 'filter_cmpHet': True, 'verbose': None, 'sep': None,
            'impact': None, 'SpliceAItag': None}
    # Run
    main_comHet(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_comHet_NA12878.out')]
    # Clean
    os.remove('tests/files/main_test.out')
    os.remove('tests/files/main_test.out.summary')
    os.remove('tests/files/main_test.out.json')
#end def

def test_run_comHet_proband_plus_NA12878_NA12877_denovo():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_comHet_denovo.vcf', 'outputfile': 'tests/files/main_test.out',
            'trio': ['NA12879_sample', 'NA12878_sample', 'NA12877_sample'], 'VEPtag': 'VEP', 'allow_undef': None, 'filter_cmpHet': True, 'verbose': None, 'sep': '~',
            'impact': None, 'SpliceAItag': None}
    # Run
    main_comHet(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_comHet_plus_NA12878_NA12877_denovo.out')]
    assert [row for row in open('tests/files/main_test.out.summary')] == [row for row in open('tests/files/input_comHet_plus_NA12878_NA12877_denovo.out.summary')]
    assert [row for row in open('tests/files/main_test.out.json')] == [row for row in open('tests/files/input_comHet_plus_NA12878_NA12877_denovo.out.json')]
    # Clean
    os.remove('tests/files/main_test.out')
    os.remove('tests/files/main_test.out.summary')
    os.remove('tests/files/main_test.out.json')
#end def

def test_run_comHet_proband_plus_NA12878_NA12877_denovo_impact_test():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_comHet_denovo_impact.vcf', 'outputfile': 'tests/files/main_test.out',
            'trio': ['NA12879_sample', 'NA12878_sample', 'NA12877_sample'], 'VEPtag': 'VEP', 'allow_undef': None, 'filter_cmpHet': True, 'verbose': None,
            'sep': '~', 'impact': True, 'SpliceAItag': None}
    # Run
    main_comHet(args, test=True)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_comHet_plus_NA12878_NA12877_denovo_impact.out')]
    assert [row for row in open('tests/files/main_test.out.summary')] == [row for row in open('tests/files/input_comHet_plus_NA12878_NA12877_denovo_impact.out.summary')]
    assert [row for row in open('tests/files/main_test.out.json')] == [row for row in open('tests/files/input_comHet_plus_NA12878_NA12877_denovo_impact.out.json')]
    # Clean
    os.remove('tests/files/main_test.out')
    os.remove('tests/files/main_test.out.summary')
    os.remove('tests/files/main_test.out.json')
#end def

def test_run_comHet_proband_plus_NA12878_NA12877_denovo_impact_dict_test():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_comHet_denovo_impact_dict.vcf', 'outputfile': 'tests/files/main_test.out',
            'trio': ['NA12879_sample', 'NA12878_sample', 'NA12877_sample'], 'VEPtag': 'VEP', 'allow_undef': None, 'filter_cmpHet': True, 'verbose': None, 'sep': '~',
            'impact': True, 'SpliceAItag': None}
    # Run
    main_comHet(args, test=True)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_comHet_plus_NA12878_NA12877_denovo_impact_dict.out')]
    assert [row for row in open('tests/files/main_test.out.summary')] == [row for row in open('tests/files/input_comHet_plus_NA12878_NA12877_denovo_impact_dict.out.summary')]
    assert [row for row in open('tests/files/main_test.out.json')] == [row for row in open('tests/files/input_comHet_plus_NA12878_NA12877_denovo_impact_dict.out.json')]
    # Clean
    os.remove('tests/files/main_test.out')
    os.remove('tests/files/main_test.out.summary')
    os.remove('tests/files/main_test.out.json')
#end def

#################################################################
#   Errors
#################################################################
def test_run_trio_error(): # NA12879
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_comHet.vcf', 'outputfile': 'tests/files/main_test.out',
            'trio': ['NA12879_sample', 'NA12878_sample', 'NA12877_sample', 'extra_sample'], 'VEPtag': 'VEP', 'allow_undef': None, 'filter_cmpHet': True, 'verbose': None, 'sep': None,
            'impact': None, 'SpliceAItag': None}
    # Run and Tests
    with pytest.raises(SystemExit) as e:
        assert main_comHet(args)
    assert str(e.value) == '\nERROR in parsing arguments: too many sample IDs provided for trio\n'
    # Clean
    os.remove('tests/files/main_test.out')
#end def
