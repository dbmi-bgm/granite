#################################################################
#   Libraries
#################################################################
import sys, os
import pytest
from granite.blackList import (
                            main as main_blackList
                            )


#################################################################
#   Tests
#################################################################
def test_run_blackList_big():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_blackList.vcf', 'outputfile': 'tests/files/main_test.out',
            'aftag': None, 'afthr': None, 'bigfile': 'tests/files/input_blackList.big'}
    # Run
    main_blackList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_blackList_big.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_blackList_032_big():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_blackList.vcf', 'outputfile': 'tests/files/main_test.out',
            'aftag': 'novoAF', 'afthr': '0.32', 'bigfile': 'tests/files/input_blackList.big'}
    # Run
    main_blackList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_blackList_0-32_big.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_blackList_032():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_blackList.vcf', 'outputfile': 'tests/files/main_test.out',
            'aftag': 'novoAF', 'afthr': '0.32', 'bigfile': None}
    # Run
    main_blackList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_blackList_0-32.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_blackList_01():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_blackList.vcf', 'outputfile': 'tests/files/main_test.out',
            'aftag': 'novoAF', 'afthr': '0.1', 'bigfile': None}
    # Run
    main_blackList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_blackList_0-1.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_blackList_001():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_blackList.vcf', 'outputfile': 'tests/files/main_test.out',
            'aftag': 'novoAF', 'afthr': '0.01', 'bigfile': None}
    # Run
    main_blackList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_blackList_0-01.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def


#################################################################
#   Errors
#################################################################
def test_args_afthr_bigfile_missing():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_blackList.vcf', 'outputfile': 'tests/files/main_test.out',
            'aftag': None, 'afthr': None, 'bigfile': None}
    # Run and Tests
    with pytest.raises(SystemExit) as e:
        assert main_blackList(args)
    assert str(e.value) == '\nERROR in parsing arguments: to blacklist specify a BIG file and/or a threshold for population allele frequency and the TAG to use\n'
#end def

def test_args_afthr_conflict():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_blackList.vcf', 'outputfile': 'tests/files/main_test.out',
            'aftag': None, 'afthr': '0.2', 'bigfile': None}
    # Run and Tests
    with pytest.raises(SystemExit) as e:
        assert main_blackList(args)
    assert str(e.value) == '\nERROR in parsing arguments: to filter by population allele frequency please specify the TAG to use\n'
#end def

def test_aftag_missing():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_blackList_novoAF_missing.vcf', 'outputfile': 'tests/files/main_test.out',
            'aftag': 'novoAF', 'afthr': '0.2', 'bigfile': None}
    # Run and Tests
    with pytest.raises(SystemExit) as e:
        assert main_blackList(args)
    assert '27100779' in str(e.value)
    # Clean
    os.remove('tests/files/main_test.out')
#end def
