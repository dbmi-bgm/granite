#################################################################
#   Libraries
#################################################################
import sys, os
import pytest
import bitarray
from granite.lib.shared_functions import *


#################################################################
#   Tests
#################################################################
def test_bed_to_bitarray():
    ''' '''
    # Run
    test_bit = bed_to_bitarray('tests/files/input_BED.bed')
    # Checks
    assert test_bit['chr1'][0] == False
    assert test_bit['chr1'][1:100].count(False) == 0
    assert test_bit['chr1'][100] == True
    assert test_bit['chr1'][101:106].count(True) == 0
    assert test_bit['chr1'][106] == True
    assert test_bit['chr2'][:54].count(True) == 0
    assert test_bit['chr2'][54] == True
    assert test_bit['X'][:2000].count(True) == 0
    assert test_bit['X'][2000:2010].count(False) == 0
    assert test_bit['X'][2010] == True
#end def

def test_bed_to_bitarray_short():
    ''' '''
    # Run
    test_bit = bed_to_bitarray('tests/files/input_BED_short.bed')
    # Checks
    assert test_bit['chr1'][0] == False
    assert test_bit['chr1'][1:100].count(False) == 0
    assert test_bit['chr1'][100] == True
#end def
