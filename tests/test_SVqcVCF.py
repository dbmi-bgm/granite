#################################################################
#   Libraries
#################################################################
import sys, os
import pytest
import json
from granite.SVqcVCF import (
                            main as main_SVqcVCF
                            )


#################################################################
#   Tests
#################################################################

def test_fail_SVqcVCF_no_SVTYPE():
    #Variables
    args = {'inputfile': 'tests/files/SVqcVCF_no_SVTYPE.vcf.gz','outputfile':'tests/files/main_test.out','samples':["TESTSAMPLE", "TESTSAMPLE2"], 'verbose': None}
    # Run and Tests
    with pytest.raises(Exception, match = "ERROR in parsing vcf, variant at position 1:46000 does not contain SVTYPE in INFO"):
        main_SVqcVCF(args)

def test_fail_SVqcVCF_wrong_SVTYPE():
    #Variables
    args = {'inputfile': 'tests/files/SVqcVCF_wrong_SVTYPE.vcf.gz','outputfile':'tests/files/main_test.out','samples':["TESTSAMPLE", "TESTSAMPLE2"], 'verbose': None}
    # Run and Tests
    with pytest.raises(Exception, match = "ERROR in parsing vcf, variant at position 1:46000 contains unexpected SVTYPE \"SNV\" in INFO"):
        main_SVqcVCF(args)

def test_success_SVqcVCF_twoSamples():
    #Variables
    args = {'inputfile': 'tests/files/SVqcVCF_success.vcf.gz','outputfile':'tests/files/main_test.out','samples':["TESTSAMPLE", "TESTSAMPLE2"], 'verbose': None}
    # Run
    main_SVqcVCF(args)
    # Tests
    with open('tests/files/main_test.out') as fi:
        d1 = json.load(fi)
    with open('tests/files/SVqcVCF_twoSamples.json') as fi:
        d2 = json.load(fi)
    assert d1 == d2
    # Clean
    os.remove('tests/files/main_test.out')

def test_success_SVqcVCF_oneSample():
    #Variables
    args = {'inputfile': 'tests/files/SVqcVCF_success.vcf.gz','outputfile':'tests/files/main_test.out','samples':["TESTSAMPLE"], 'verbose': None}
    # Run
    main_SVqcVCF(args)
    # Tests
    with open('tests/files/main_test.out') as fi:
        d1 = json.load(fi)
    with open('tests/files/SVqcVCF_oneSample.json') as fi:
        d2 = json.load(fi)
    assert d1 == d2
    # Clean
    os.remove('tests/files/main_test.out')
