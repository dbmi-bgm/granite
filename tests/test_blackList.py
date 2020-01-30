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
def test_run_blackList_bgi():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_blackList.vcf', 'outputfile': 'tests/files/main_test.out',
            'aftag': None, 'afthr': None, 'bgifile': 'tests/files/input_blackList.bgi'}
    # Run
    main_blackList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_blackList_bgi.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_blackList_032_bgi():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_blackList.vcf', 'outputfile': 'tests/files/main_test.out',
            'aftag': 'novoAF', 'afthr': '0.32', 'bgifile': 'tests/files/input_blackList.bgi'}
    # Run
    main_blackList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_blackList_0-32_bgi.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_blackList_032():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_blackList.vcf', 'outputfile': 'tests/files/main_test.out',
            'aftag': 'novoAF', 'afthr': '0.32', 'bgifile': None}
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
            'aftag': 'novoAF', 'afthr': '0.1', 'bgifile': None}
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
            'aftag': 'novoAF', 'afthr': '0.01', 'bgifile': None}
    # Run
    main_blackList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_blackList_0-01.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def
