#################################################################
#   Libraries
#################################################################
import sys, os
import pytest
from granite.qcVCF import (
                            main as main_qcVCF
                            )


#################################################################
#   Tests
#################################################################
def test_run_qcVCF_string():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/NA_12879.vcf.gz', 'outputfile': 'tests/files/main_test.out',
            'samples': ["NA12879_sample", "NA12878_sample", "NA12877_sample"],
            'trio_errors': True, 'het_hom': False, 'ti_tv': False, 'verbose': None,
            'pedigree': '[{"gender": "M", "individual": "GAPIDERSF9QH", "parents": [], "sample_name": "NA12877_sample"}, {"gender": "F", "individual": "GAPIDGE5KLIH", "parents": [], "sample_name": "NA12878_sample"}, {"gender": "F", "individual": "GAPID832VXDG", "parents": ["GAPIDGE5KLIH", "GAPIDERSF9QH"], "sample_name": "NA12879_sample"}]'}
    # Run
    main_qcVCF(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/NA_12879.qcvcf.json')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_qcVCF_file():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/NA_12879.vcf.gz', 'outputfile': 'tests/files/main_test.out',
            'samples': ["NA12879_sample", "NA12878_sample", "NA12877_sample"],
            'trio_errors': True, 'het_hom': False, 'ti_tv': False, 'verbose': None,
            'pedigree': 'tests/files/ped_trio.json'}
    # Run
    main_qcVCF(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/NA_12879.qcvcf.json')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def
