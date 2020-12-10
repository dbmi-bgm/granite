#################################################################
#   Libraries
#################################################################
import sys, os
import pytest
from granite.geneList import (
                            main as main_geneList
                            )


#################################################################
#   Tests
#################################################################
def test_run_geneList():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_geneList.vcf', 'outputfile': 'tests/files/main_test.out',
            'geneslist': 'tests/files/ENSG_list.txt', 'VEPtag': 'VEP', 'verbose': None}
    # Run
    main_geneList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_geneList.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_geneList_CSQ():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_geneList_CSQ.vcf', 'outputfile': 'tests/files/main_test.out',
            'geneslist': 'tests/files/ENSG_list_CSQ.txt', 'VEPtag': None, 'verbose': None}
    # Run
    main_geneList(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_geneList_CSQ.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def
