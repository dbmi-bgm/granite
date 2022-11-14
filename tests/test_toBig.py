#################################################################
#   Libraries
#################################################################
import sys, os
import pytest
import bitarray
from granite.toBig import (
                            main as main_toBig
                            )
from granite.lib.shared_functions import *


#################################################################
#   Tests
#################################################################
def test_run_toBig_rdthr_2_all():
    ''' '''
    # Variables
    args = {'file': ['tests/files/input_toBig_1.rck.gz', 'tests/files/input_toBig_2.rck.gz',
            'tests/files/input_toBig_3.rck.gz'], 'outputfile': 'tests/files/main_test.out',
            'fithr': '3', 'rdthr': '2', 'ncores': '2', 'abthr': None,
            'regionfile': 'tests/files/input_toBig.regions',
            'chromfile': 'tests/files/input_toBig.chrom.size'}
    # Run
    main_toBig(args)
    # Expected
    snv_expect = [11001, 11007, 11010]
    ins_expect = [11005, 11022]
    del_expect = [11017, 11030]
    # Tests
    bit = bitarray.bitarray(11030 + 1)
    big_dict = load_big('tests/files/main_test.out')
    # Check snv
    bit.setall(False)
    for i in snv_expect:
        bit[i] = True
    #end for
    assert big_dict['13_snv'][:11031] == bit
    # Check ins
    bit.setall(False)
    for i in ins_expect:
        bit[i] = True
    #end for
    assert big_dict['13_ins'][:11031] == bit
    # Check del
    bit.setall(False)
    for i in del_expect:
        bit[i] = True
    #end for
    assert big_dict['13_del'][:11031] == bit
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_toBig_rdthr_2_2():
    ''' '''
    # Variables
    args = {'file': ['tests/files/input_toBig_1.rck.gz', 'tests/files/input_toBig_2.rck.gz',
            'tests/files/input_toBig_3.rck.gz'], 'outputfile': 'tests/files/main_test.out',
            'fithr': '2', 'rdthr': '2', 'ncores': '2', 'abthr': None,
            'regionfile': 'tests/files/input_toBig.regions',
            'chromfile': 'tests/files/input_toBig.chrom.size'}
    # Run
    main_toBig(args)
    # Expected
    snv_expect = [11001, 11002, 11007, 11010, 11013, 11023]
    ins_expect = [11005, 11022]
    del_expect = [11017, 11030]
    # Tests
    bit = bitarray.bitarray(11030 + 1)
    big_dict = load_big('tests/files/main_test.out')
    # Check snv
    bit.setall(False)
    for i in snv_expect:
        bit[i] = True
    #end for
    assert big_dict['13_snv'][:11031] == bit
    # Check ins
    bit.setall(False)
    for i in ins_expect:
        bit[i] = True
    #end for
    assert big_dict['13_ins'][:11031] == bit
    # Check del
    bit.setall(False)
    for i in del_expect:
        bit[i] = True
    #end for
    assert big_dict['13_del'][:11031] == bit
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_toBig_rdthr_17_all():
    ''' '''
    # Variables
    args = {'file': ['tests/files/input_toBig_1.rck.gz', 'tests/files/input_toBig_2.rck.gz',
            'tests/files/input_toBig_3.rck.gz'], 'outputfile': 'tests/files/main_test.out',
            'fithr': '3', 'rdthr': '17', 'ncores': None, 'abthr': None,
            'regionfile': 'tests/files/input_toBig.regions',
            'chromfile': 'tests/files/input_toBig.chrom.size'}
    # Run
    main_toBig(args)
    # Expected
    snv_expect = [11007]
    ins_expect = []
    del_expect = []
    # Tests
    bit = bitarray.bitarray(11030 + 1)
    big_dict = load_big('tests/files/main_test.out')
    # Check snv
    bit.setall(False)
    for i in snv_expect:
        bit[i] = True
    #end for
    assert big_dict['13_snv'][:11031] == bit
    # Check ins
    bit.setall(False)
    for i in ins_expect:
        bit[i] = True
    #end for
    assert big_dict['13_ins'][:11031] == bit
    # Check del
    bit.setall(False)
    for i in del_expect:
        bit[i] = True
    #end for
    assert big_dict['13_del'][:11031] == bit
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_toBig_rdthr_17_2():
    ''' '''
    # Variables
    args = {'file': ['tests/files/input_toBig_1.rck.gz', 'tests/files/input_toBig_2.rck.gz',
            'tests/files/input_toBig_3.rck.gz'], 'outputfile': 'tests/files/main_test.out',
            'fithr': '2', 'rdthr': '17', 'ncores': '1', 'abthr': None,
            'regionfile': 'tests/files/input_toBig.regions',
            'chromfile': 'tests/files/input_toBig.chrom.size'}
    # Run
    main_toBig(args)
    # Expected
    snv_expect = [11007, 11010]
    ins_expect = []
    del_expect = []
    # Tests
    bit = bitarray.bitarray(11030 + 1)
    big_dict = load_big('tests/files/main_test.out')
    # Check snv
    bit.setall(False)
    for i in snv_expect:
        bit[i] = True
    #end for
    assert big_dict['13_snv'][:11031] == bit
    # Check ins
    bit.setall(False)
    for i in ins_expect:
        bit[i] = True
    #end for
    assert big_dict['13_ins'][:11031] == bit
    # Check del
    bit.setall(False)
    for i in del_expect:
        bit[i] = True
    #end for
    assert big_dict['13_del'][:11031] == bit
    # Clean
    os.remove('tests/files/main_test.out')
#end def


def test_run_toBig_abthr_15_all():
    ''' '''
    # Variables
    args = {'file': ['tests/files/input_toBig_1.rck.gz', 'tests/files/input_toBig_2.rck.gz',
            'tests/files/input_toBig_3.rck.gz'], 'outputfile': 'tests/files/main_test.out',
            'fithr': '3', 'rdthr': None, 'ncores': '2', 'abthr': None,
            'regionfile': 'tests/files/input_toBig.regions',
            'chromfile': 'tests/files/input_toBig.chrom.size'}
    # Run
    main_toBig(args)
    # Expected
    snv_expect = [11007, 11010, 11013]
    ins_expect = [11022]
    del_expect = [11030]
    # Tests
    bit = bitarray.bitarray(11030 + 1)
    big_dict = load_big('tests/files/main_test.out')
    # Check snv
    bit.setall(False)
    for i in snv_expect:
        bit[i] = True
    #end for
    assert big_dict['13_snv'][:11031] == bit
    # Check ins
    bit.setall(False)
    for i in ins_expect:
        bit[i] = True
    #end for
    assert big_dict['13_ins'][:11031] == bit
    # Check del
    bit.setall(False)
    for i in del_expect:
        bit[i] = True
    #end for
    assert big_dict['13_del'][:11031] == bit
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_toBig_abthr_25_all():
    ''' '''
    # Variables
    args = {'file': ['tests/files/input_toBig_1.rck.gz', 'tests/files/input_toBig_2.rck.gz',
            'tests/files/input_toBig_3.rck.gz'], 'outputfile': 'tests/files/main_test.out',
            'fithr': '3', 'rdthr': None, 'ncores': '2', 'abthr': '25',
            'regionfile': 'tests/files/input_toBig.regions',
            'chromfile': 'tests/files/input_toBig.chrom.size'}
    # Run
    main_toBig(args)
    # Expected
    snv_expect = [11010]
    ins_expect = []
    del_expect = []
    # Tests
    bit = bitarray.bitarray(11030 + 1)
    big_dict = load_big('tests/files/main_test.out')
    # Check snv
    bit.setall(False)
    for i in snv_expect:
        bit[i] = True
    #end for
    assert big_dict['13_snv'][:11031] == bit
    # Check ins
    bit.setall(False)
    for i in ins_expect:
        bit[i] = True
    #end for
    assert big_dict['13_ins'][:11031] == bit
    # Check del
    bit.setall(False)
    for i in del_expect:
        bit[i] = True
    #end for
    assert big_dict['13_del'][:11031] == bit
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_toBig_abthr_25_2():
    ''' '''
    # Variables
    args = {'file': ['tests/files/input_toBig_1.rck.gz', 'tests/files/input_toBig_2.rck.gz',
            'tests/files/input_toBig_3.rck.gz'], 'outputfile': 'tests/files/main_test.out',
            'fithr': '2', 'rdthr': None, 'ncores': '1', 'abthr': '25',
            'regionfile': 'tests/files/input_toBig.regions',
            'chromfile': 'tests/files/input_toBig.chrom.size'}
    # Run
    main_toBig(args)
    # Expected
    snv_expect = [11007, 11010]
    ins_expect = []
    del_expect = []
    # Tests
    bit = bitarray.bitarray(11030 + 1)
    big_dict = load_big('tests/files/main_test.out')
    # Check snv
    bit.setall(False)
    for i in snv_expect:
        bit[i] = True
    #end for
    assert big_dict['13_snv'][:11031] == bit
    # Check ins
    bit.setall(False)
    for i in ins_expect:
        bit[i] = True
    #end for
    assert big_dict['13_ins'][:11031] == bit
    # Check del
    bit.setall(False)
    for i in del_expect:
        bit[i] = True
    #end for
    assert big_dict['13_del'][:11031] == bit
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_toBig_rdthr_2_1_single():
    ''' '''
    # Variables
    args = {'file': ['tests/files/input_toBig_1.rck.gz'],
            'outputfile': 'tests/files/main_test.out',
            'fithr': '1', 'rdthr': '2', 'ncores': '2', 'abthr': None,
            'regionfile': 'tests/files/input_toBig.regions',
            'chromfile': 'tests/files/input_toBig.chrom.size'}
    # Run
    main_toBig(args)
    # Expected
    snv_expect = [11001, 11002, 11007, 11010, 11013, 11023]
    ins_expect = [11005, 11022]
    del_expect = [11017, 11030]
    # Tests
    bit = bitarray.bitarray(11030 + 1)
    big_dict = load_big('tests/files/main_test.out')
    # Check snv
    bit.setall(False)
    for i in snv_expect:
        bit[i] = True
    #end for
    assert big_dict['13_snv'][:11031] == bit
    # Check ins
    bit.setall(False)
    for i in ins_expect:
        bit[i] = True
    #end for
    assert big_dict['13_ins'][:11031] == bit
    # Check del
    bit.setall(False)
    for i in del_expect:
        bit[i] = True
    #end for
    assert big_dict['13_del'][:11031] == bit
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_toBig_rdthr_2_2_single():
    ''' '''
    # Variables
    args = {'file': ['tests/files/input_toBig_1.rck.gz'],
            'outputfile': 'tests/files/main_test.out',
            'fithr': '2', 'rdthr': '2', 'ncores': '2', 'abthr': None,
            'regionfile': 'tests/files/input_toBig.regions',
            'chromfile': 'tests/files/input_toBig.chrom.size'}
    # Run
    main_toBig(args)
    # Expected
    snv_expect = []
    ins_expect = []
    del_expect = []
    # Tests
    bit = bitarray.bitarray(11030 + 1)
    big_dict = load_big('tests/files/main_test.out')
    # Check snv
    bit.setall(False)
    for i in snv_expect:
        bit[i] = True
    #end for
    assert big_dict['13_snv'][:11031] == bit
    # Check ins
    bit.setall(False)
    for i in ins_expect:
        bit[i] = True
    #end for
    assert big_dict['13_ins'][:11031] == bit
    # Check del
    bit.setall(False)
    for i in del_expect:
        bit[i] = True
    #end for
    assert big_dict['13_del'][:11031] == bit
    # Clean
    os.remove('tests/files/main_test.out')
#end def


#################################################################
#   Errors
#################################################################
def test_run_toBig_rdthr_2_all_miss_pos():
    ''' '''
    # Variables
    args = {'file': ['tests/files/input_toBig_1.rck.gz', 'tests/files/input_toBig_2.rck.gz',
            'tests/files/input_toBig_3.rck.gz', 'tests/files/input_toBig_miss_pos.rck.gz'], 'outputfile': 'tests/files/main_test.out',
            'fithr': '4', 'rdthr': '2', 'ncores': '2', 'abthr': None,
            'regionfile': 'tests/files/input_toBig.regions',
            'chromfile': 'tests/files/input_toBig.chrom.size'}
    # Run and Tests
    with pytest.raises(IndexError) as e:
        assert main_toBig(args)
    assert str(e.value) == '\nERROR in file: position 13:11006 in file tests/files/input_toBig_miss_pos.rck.gz is not consistent with other files\n'
#end def
