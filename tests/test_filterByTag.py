#################################################################
#   Libraries
#################################################################
import sys, os
import pytest
from granite.filterByTag import (
                            main as main_filterByTag
                            )

#################################################################
#   Tests
#################################################################
def test_run_gnomADc_AF_joint_0_01():
    ''' '''
    # Variables
    inputfile = 'input_filterByTag_gnomADc_AF_joint.vcf'
    compare_outputfile = 'input_filterByTag_gnomADc_AF_joint_0.01.out'
    args = {'inputfile': f'tests/files/{inputfile}', 'outputfile': 'tests/files/main_test.out',
            'logic': 'any', 'verbose': None, 'separator': ';',
            'tag': [
                'gnomADc_AF_joint/0.01/</float/any/field=|/entry=,'
            ]}
    # Run
    main_filterByTag(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open(f'tests/files/{compare_outputfile}')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_gnomADc_AF_joint_0_1():
    ''' '''
    # Variables
    inputfile = 'input_filterByTag_gnomADc_AF_joint.vcf'
    compare_outputfile = 'input_filterByTag_gnomADc_AF_joint_0.1.out'
    args = {'inputfile': f'tests/files/{inputfile}', 'outputfile': 'tests/files/main_test.out',
            'logic': 'any', 'verbose': None, 'separator': ';',
            'tag': [
                'gnomADc_AF_joint/0.1/>=/float/all/field=|/entry=,'
            ]}
    # Run
    main_filterByTag(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open(f'tests/files/{compare_outputfile}')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_gnomADc_AF_joint_0_01_0_1():
    ''' '''
    # Variables
    inputfile = 'input_filterByTag_gnomADc_AF_joint.vcf'
    compare_outputfile = 'input_filterByTag_gnomADc_AF_joint_0.01_0.1.out'
    args = {'inputfile': f'tests/files/{inputfile}', 'outputfile': 'tests/files/main_test.out',
            'logic': 'any', 'verbose': None, 'separator': ';',
            'tag': [
                'gnomADc_AF_joint/0.01/</float/any/field=|/entry=,',
                'gnomADc_AF_joint/0.1/>=/float/all/field=|/entry=,'
            ]}
    # Run
    main_filterByTag(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open(f'tests/files/{compare_outputfile}')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_gnomADc_AF_joint_AN_grpmax_genomes_any():
    ''' '''
    # Variables
    inputfile = 'input_filterByTag_gnomADc_AF_joint.vcf'
    compare_outputfile = 'input_filterByTag_gnomADc_AF_joint_AN_grpmax_genomes_any.out'
    args = {'inputfile': f'tests/files/{inputfile}', 'outputfile': 'tests/files/main_test.out',
            'logic': 'any', 'verbose': None, 'separator': ';',
            'tag': [
                'gnomADc_AF_joint/0.1/>=/float/all/field=|/entry=,',
                'gnomADc_AN_grpmax_genomes/10/>/int/any/field=|/entry=,',
                'gnomADc_AN_grpmax_genomes/9/==/int/any/field=|/entry=,'
            ]}
    # Run
    main_filterByTag(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open(f'tests/files/{compare_outputfile}')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_gnomADc_AF_joint_AN_grpmax_genomes_all():
    ''' '''
    # Variables
    inputfile = 'input_filterByTag_gnomADc_AF_joint.vcf'
    compare_outputfile = 'input_filterByTag_gnomADc_AF_joint_AN_grpmax_genomes_all.out'
    args = {'inputfile': f'tests/files/{inputfile}', 'outputfile': 'tests/files/main_test.out',
            'logic': 'all', 'verbose': None, 'separator': ';',
            'tag': [
                'gnomADc_AF_joint/0.1/>/float/all/field=|/entry=,',
                'gnomADc_AN_grpmax_genomes/9/>=/int/all/field=|/entry=,'
            ]}
    # Run
    main_filterByTag(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open(f'tests/files/{compare_outputfile}')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_IMPACT_FLAG_HIGH_true_any():
    ''' '''
    # Variables
    inputfile = 'input_filterByTag_IMPACT_FLAG.vcf'
    compare_outputfile = 'input_filterByTag_IMPACT_FLAG_HIGH_true_any.out'
    args = {'inputfile': f'tests/files/{inputfile}', 'outputfile': 'tests/files/main_test.out',
            'logic': 'any', 'verbose': None, 'separator': ';',
            'tag': [
                'IMPACT/HIGH/==/str/all/field=|/entry=,',
                'PIPPO/-/true/bool/any'
            ]}
    # Run
    main_filterByTag(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open(f'tests/files/{compare_outputfile}')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_IMPACT_FLAG_not_MODIFIER_true_any():
    ''' '''
    # Variables
    inputfile = 'input_filterByTag_IMPACT_FLAG.vcf'
    compare_outputfile = 'input_filterByTag_IMPACT_FLAG_HIGH_true_any.out'
    args = {'inputfile': f'tests/files/{inputfile}', 'outputfile': 'tests/files/main_test.out',
            'logic': 'any', 'verbose': None, 'separator': ';',
            'tag': [
                'IMPACT/MODIFIER/!=/str/all/field=|/entry=,',
                'PIPPO/-/true/bool/any'
            ]}
    # Run
    main_filterByTag(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open(f'tests/files/{compare_outputfile}')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_IMPACT_FLAG_HIGH_true_all():
    ''' '''
    # Variables
    inputfile = 'input_filterByTag_IMPACT_FLAG.vcf'
    compare_outputfile = 'input_filterByTag_IMPACT_FLAG_HIGH_true_all.out'
    args = {'inputfile': f'tests/files/{inputfile}', 'outputfile': 'tests/files/main_test.out',
            'logic': 'all', 'verbose': None, 'separator': ';',
            'tag': [
                'IMPACT/HIG/~/str/any/field=|/entry=,',
                'PIPPO/-/true/bool/all'
            ]}
    # Run
    main_filterByTag(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open(f'tests/files/{compare_outputfile}')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_IMPACT_FLAG_MODIF_notinstring_true_all():
    ''' '''
    # Variables
    inputfile = 'input_filterByTag_IMPACT_FLAG.vcf'
    compare_outputfile = 'input_filterByTag_IMPACT_FLAG_HIGH_true_all.out'
    args = {'inputfile': f'tests/files/{inputfile}', 'outputfile': 'tests/files/main_test.out',
            'logic': 'all', 'verbose': None, 'separator': ';',
            'tag': [
                'IMPACT/MODIF/!~/str/any/field=|/entry=,',
                'PIPPO/-/true/bool/all'
            ]}
    # Run
    main_filterByTag(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open(f'tests/files/{compare_outputfile}')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_IMPACT_FLAG_false():
    ''' '''
    # Variables
    inputfile = 'input_filterByTag_IMPACT_FLAG.vcf'
    compare_outputfile = 'input_filterByTag_IMPACT_FLAG_false.out'
    args = {'inputfile': f'tests/files/{inputfile}', 'outputfile': 'tests/files/main_test.out',
            'logic': 'any', 'verbose': None, 'separator': ';',
            'tag': [
                'PIPPO/-/false/bool/all'
            ]}
    # Run
    main_filterByTag(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open(f'tests/files/{compare_outputfile}')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_missing_tag_value():
    ''' '''
    # Variables
    inputfile = 'input_filterByTag_missing_tag_value.vcf'
    compare_outputfile = 'input_filterByTag_missing_tag_value.out'
    args = {'inputfile': f'tests/files/{inputfile}', 'outputfile': 'tests/files/main_test.out',
            'logic': 'any', 'verbose': None, 'separator': ';',
            'tag': [
                'gnomADc_AF_joint/0.01/</float/all/field=|/entry=,',
                'DP/100/<=/int/any',
                'IMPACT/HIGH/==/str/all/field=|/entry=,'
            ]}
    # Run
    main_filterByTag(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open(f'tests/files/{compare_outputfile}')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_missing_tag_value_missing_field():
    ''' '''
    # Variables
    inputfile = 'input_filterByTag_missing_tag_value.vcf'
    compare_outputfile = 'input_filterByTag_missing_tag_value.out'
    args = {'inputfile': f'tests/files/{inputfile}', 'outputfile': 'tests/files/main_test.out',
            'logic': 'any', 'verbose': None, 'separator': ';',
            'tag': [
                'gnomADc_AF_joint/0.01/</float/all/entry=,',
                'DP/100/<=/int/any',
                'IMPACT/HIGH/==/str/all/entry=,'
            ]}
    # Run
    main_filterByTag(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open(f'tests/files/{compare_outputfile}')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_errors():
    ''' '''
    # Variables
    inputfile = 'input_filterByTag_missing_tag_value.vcf'
    compare_outputfile = 'input_filterByTag_missing_tag_value.out'
    args = {'inputfile': f'tests/files/{inputfile}', 'outputfile': 'tests/files/main_test.out',
            'logic': 'any', 'verbose': None, 'separator': ';',
            'tag': []}
    # Run
    args['tag'] = ['IMPACT/HIGH/==/str/all/entry-,']
    with pytest.raises(SystemExit) as e:
        assert main_filterByTag(args)
    assert '\nERROR in tag filter: malformed option "entry-," in IMPACT/HIGH/==/str/all/entry-,. Expected key=value\n' == str(e.value)

    args['tag'] = ['DP/100/<=/str/any']
    with pytest.raises(SystemExit) as e:
        assert main_filterByTag(args)
    assert '\nERROR in tag filter type for operator <=: str\n' == str(e.value)

    args['tag'] = ['IMPACT/HIGH/==/float/all/field=|/entry=,']
    with pytest.raises(SystemExit) as e:
        assert main_filterByTag(args)
    assert '\nERROR in tag filter value: cannot convert "HIGH" to float\n' == str(e.value)

    args['tag'] = ['IMPACT/HIGH/==/all/field=|/entry=,']
    with pytest.raises(SystemExit) as e:
        assert main_filterByTag(args)
    assert '\nERROR in tag filter type: all. Accepted: str, int, float, bool\n' == str(e.value)

    args['tag'] = ['IMPACT/HIGH/==/str']
    with pytest.raises(SystemExit) as e:
        assert main_filterByTag(args)
    assert '\nERROR in tag filter format: IMPACT/HIGH/==/str. Expected format: name/value/operator/type/logic[/field=sep/entry=sep]\n' == str(e.value)
#end def