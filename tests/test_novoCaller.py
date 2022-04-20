#################################################################
#   Libraries
#################################################################
import sys, os
import pytest
from granite.lib import vcf_parser
from granite.novoCaller import (
                            main as main_novoCaller
                            )

# define flexible assert function to compare vcf files
def assert_(tfile, ofile, test=False):
    ''' '''
    # Creating Vcf object
    vcf_obj_t = vcf_parser.Vcf(tfile)
    vcf_obj_o = vcf_parser.Vcf(ofile)
    # Check headers
    assert vcf_obj_t.header.definitions == vcf_obj_o.header.definitions
    # Get variants
    vnt_obj_t_list, vnt_obj_o_list = [], []
    for vnt_obj in vcf_obj_t.parse_variants(): vnt_obj_t_list.append(vnt_obj)
    for vnt_obj in vcf_obj_o.parse_variants(): vnt_obj_o_list.append(vnt_obj)
    # Check variants number is consistent
    assert len(vnt_obj_t_list) == len(vnt_obj_o_list)
    # Check variants
    for i, vnt_obj_t in enumerate(vnt_obj_t_list):
        vnt_obj_o = vnt_obj_o_list[i]
        # Check basic fields are matching
        assert vnt_obj_t.CHROM == vnt_obj_o.CHROM
        assert vnt_obj_t.POS == vnt_obj_o.POS
        assert vnt_obj_t.ID == vnt_obj_o.ID
        assert vnt_obj_t.REF == vnt_obj_o.REF
        assert vnt_obj_t.ALT == vnt_obj_o.ALT
        assert vnt_obj_t.QUAL == vnt_obj_o.QUAL
        assert vnt_obj_t.FILTER == vnt_obj_o.FILTER
        assert vnt_obj_t.INFO.split(';')[:-2] == vnt_obj_o.INFO.split(';')[:-2]
        assert vnt_obj_t.FORMAT == vnt_obj_o.FORMAT
        # novoAF
        #   since there may be minor differences in rounding
        #   we check that the results are close enough to the expected
        try:
            af_t = float(vnt_obj_t.get_tag_value('novoAF'))
            af_o = float(vnt_obj_o.get_tag_value('novoAF'))
        except ValueError:
            af_t = float(vnt_obj_t.get_tag_value('AlleleFreq'))
            af_o = float(vnt_obj_o.get_tag_value('AlleleFreq'))
        assert abs(af_t - af_o) < 0.0000000001
        # novoPP
        #   since there may be minor differences in rounding
        #   we check that the results are close enough to the expected
        try:
            pp_t = float(vnt_obj_t.get_tag_value('novoPP'))
            pp_o = float(vnt_obj_o.get_tag_value('novoPP'))
            if pp_t > 0.1: assert abs(pp_t - pp_o) < 0.0000000001
            else: assert abs(pp_t - pp_o) < 0.0001
        except ValueError as e:
            assert str(e) == '\nERROR in variant INFO field, novoPP tag is missing\n'
        # genotypes
        for id in vnt_obj_t.IDs_genotypes:
            try:
                GT_t, RSTR_t = vnt_obj_t.get_genotype_value(id, 'GT'), list(map(int, vnt_obj_t.get_genotype_value(id, 'RSTR').split(',')))
                GT_o, RSTR_o = vnt_obj_o.get_genotype_value(id, 'GT'), list(map(int, vnt_obj_o.get_genotype_value(id, 'RSTR').split(',')))
                ref_t, alt_t = RSTR_t[0] + RSTR_t[2], RSTR_t[1] + RSTR_t[3]
                ref_o, alt_o = RSTR_o[0] + RSTR_o[2], RSTR_o[1] + RSTR_o[3]
                assert GT_t == GT_o
                assert ref_t == ref_o
                assert alt_t == alt_o
            except ValueError as e:
                assert str(e) == '\nERROR in GENOTYPES identifiers, PSC-01-003 identifier is missing in VCF\n'

#################################################################
#   Tests
#################################################################
# this tests may fail because of differences in float handling
# e.g. novoCaller=1.4030838153054451e-05 Vs novoCaller=1.403083815305445e-05
# depending on the machine
def test_run_novoCaller_rck_all():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_RCK.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_rck.tsv', 'triofiles':'tests/files/trio_rck.tsv',
            'ppthr': None, 'afthr': None, 'aftag': None, 'bam': None,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': None, 'ADthr': 3, 'verbose': None}
    # Run
    main_novoCaller(args)
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_novoCaller_RCK_all.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_noUNR():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_noUNR.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': '0.01', 'aftag': 'novoAF', 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': 0.01, 'ADthr': 3, 'verbose': None}
    # Run
    main_novoCaller(args)
    # Tests
    testfile = 'tests/files/main_test.out'
    outfile = 'tests/files/input_novoCaller_BAM_noUNR.out'
    assert_(testfile, outfile)
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_noUNR_NA():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_noUNR.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': '0.01', 'aftag': 'novoAF', 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': 0.01, 'ADthr': 3, 'verbose': None}
    # Run
    main_novoCaller(args, test=True)
    # Tests
    testfile = 'tests/files/main_test.out'
    outfile = 'tests/files/input_novoCaller_BAM_noUNR_NA.out'
    assert_(testfile, outfile, test=True)
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_noUNR_noafthr():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_noUNR.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': None, 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': 0.01, 'ADthr': 3, 'verbose': None}
    # Run
    main_novoCaller(args)
    # Tests
    testfile = 'tests/files/main_test.out'
    outfile = 'tests/files/input_novoCaller_BAM_noUNR_noafthr.out'
    assert_(testfile, outfile)
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_withUNR():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_withUNR.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': '0.01', 'aftag': 'novoAF', 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': 0.01, 'ADthr': 3, 'verbose': None}
    # Run
    main_novoCaller(args)
    # Tests
    testfile = 'tests/files/main_test.out'
    outfile = 'tests/files/input_novoCaller_BAM_withUNR.out'
    assert_(testfile, outfile)
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_withUNR_asSibling_ALLBAMS():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_withUNR_asSibling.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': '0.01', 'aftag': 'novoAF', 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': 0.01, 'ADthr': 3, 'verbose': None}
    # Run
    main_novoCaller(args)
    # Tests
    testfile = 'tests/files/main_test.out'
    outfile = 'tests/files/input_novoCaller_BAM_withUNR_asSibling_ALLBAMS.out'
    assert_(testfile, outfile)
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_withUNR_asSibling_missingSibling():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_withUNR_asSibling.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam_PSC-01-003_asSibling.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': '0.01', 'aftag': 'novoAF', 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': 0.01, 'ADthr': 3, 'verbose': None}
    # Run
    main_novoCaller(args)
    # Tests
    testfile = 'tests/files/main_test.out'
    outfile = 'tests/files/input_novoCaller_BAM_withUNR_asSibling_missingSibling.out'
    assert_(testfile, outfile)
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_withUNR_asSibling_SWAP_SON_missingSibling():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_withUNR_asSibling_SWAP_SON.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam_PSC-01-003_asSibling.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': '0.01', 'aftag': 'novoAF', 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': 0.01, 'ADthr': 3, 'verbose': None}
    # Run
    main_novoCaller(args)
    # Tests
    testfile = 'tests/files/main_test.out'
    outfile = 'tests/files/input_novoCaller_BAM_withUNR_asSibling_SWAP_SON_missingSibling.out'
    assert_(testfile, outfile)
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_rerun_noUNR():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_noUNR.out', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': '0.01', 'aftag': 'novoAF', 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': 0.01, 'ADthr': 3, 'verbose': None}
    # Run
    main_novoCaller(args)
    # Tests
    testfile = 'tests/files/main_test.out'
    outfile = 'tests/files/input_novoCaller_BAM_noUNR.out'
    assert_(testfile, outfile)
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_rerun_withUNR_asSibling_missingSibling():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_withUNR_asSibling_missingSibling.out', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam_PSC-01-003_asSibling.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': '0.01', 'aftag': 'novoAF', 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': 0.01, 'ADthr': 3, 'verbose': None}
    # Run
    main_novoCaller(args)
    # Tests
    testfile = 'tests/files/main_test.out'
    outfile = 'tests/files/input_novoCaller_BAM_withUNR_asSibling_missingSibling.out'
    assert_(testfile, outfile)
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_different_tag():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_noUNR_AlleleFreq.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': '0.01', 'aftag': 'AlleleFreq', 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': 0.01, 'ADthr': 3, 'verbose': None}
    # Run
    main_novoCaller(args)
    # Tests
    testfile = 'tests/files/main_test.out'
    outfile = 'tests/files/input_novoCaller_BAM_noUNR_AlleleFreq.out'
    assert_(testfile, outfile)
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_jc50_wgenome_plus_indels():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': '0.01', 'aftag': 'novoAF', 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': 0.01, 'ADthr': 3, 'verbose': None}
    # Run
    main_novoCaller(args)
    # Tests
    testfile = 'tests/files/main_test.out'
    outfile = 'tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels.out'
    assert_(testfile, outfile)
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_jc50_wgenome_plus_indels_11():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': None, 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': 0.01, 'ADthr': 3, 'verbose': None}
    # Run
    main_novoCaller(args)
    # Tests
    testfile = 'tests/files/main_test.out'
    outfile = 'tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels_11.out'
    assert_(testfile, outfile)
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_jc50_wgenome_plus_indels_9():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': '7.579382496057039e-05', 'afthr': None, 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': 0.01, 'ADthr': 3, 'verbose': None}
    # Run
    main_novoCaller(args)
    # Tests
    testfile = 'tests/files/main_test.out'
    outfile = 'tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels_9.out'
    assert_(testfile, outfile)
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_jc50_wgenome_plus_indels_9_NA():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': '7.579382496057039e-05', 'afthr': None, 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': 0.01, 'ADthr': 3, 'verbose': None}
    # Run
    main_novoCaller(args, test=True)
    # Tests
    testfile = 'tests/files/main_test.out'
    outfile = 'tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels_9_NA.out'
    assert_(testfile, outfile)
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_jc50_wgenome_plus_indels_8():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': '0.9999998', 'afthr': None, 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': 0.01, 'ADthr': 3, 'verbose': None}
    # Run
    main_novoCaller(args)
    # Tests
    testfile = 'tests/files/main_test.out'
    outfile = 'tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels_8.out'
    assert_(testfile, outfile)
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_jc50_wgenome_plus_indels_7():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': '0.9999999', 'afthr': None, 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': 0.01, 'ADthr': 3, 'verbose': None}
    # Run
    main_novoCaller(args)
    # Tests
    testfile = 'tests/files/main_test.out'
    outfile = 'tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels_7.out'
    assert_(testfile, outfile)
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_jc50_wgenome_plus_indels_4():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': '0.999999999999', 'afthr': None, 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': 0.01, 'ADthr': 3, 'verbose': None}
    # Run
    main_novoCaller(args)
    # Tests
    testfile = 'tests/files/main_test.out'
    outfile = 'tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels_4.out'
    assert_(testfile, outfile)
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_jc50_wgenome_plus_indels_2():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': '1', 'afthr': None, 'aftag': None, 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': 0.01, 'ADthr': 3, 'verbose': None}
    # Run
    main_novoCaller(args)
    # Tests
    testfile = 'tests/files/main_test.out'
    outfile = 'tests/files/input_novoCaller_BAM_jc50_wgenome_plus_indels_2.out'
    assert_(testfile, outfile)
    # Clean
    os.remove('tests/files/main_test.out')
#end def


#################################################################
#   Errors
#################################################################
# def test_run_novoCaller_bam_annot_missing_novoAF():
#     ''' '''
#     # Variables
#     args = {'inputfile': 'tests/files/input_novoCaller_BAM_annot_missing_novoAF.vcf', 'outputfile': 'tests/files/main_test.out',
#             'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
#             'ppthr': None, 'afthr': '0.01', 'aftag': None, 'bam': True,
#             'MQthr': None, 'BQthr': None}
#     # Run and Tests
#     with pytest.raises(SystemExit) as e:
#         assert main_novoCaller(args)
#     assert '\nERROR in variant parsing: novoAF allele frequency tag in INFO field is missing for variant:\n' in str(e.value)
#     assert '16805'
#     # Clean
#     os.remove('tests/files/main_test.out')
# #end def
#
# def test_run_novoCaller_bam_annot_wrong_novoAF():
#     ''' '''
#     # Variables
#     args = {'inputfile': 'tests/files/input_novoCaller_BAM_annot_wrong_novoAF.vcf', 'outputfile': 'tests/files/main_test.out',
#             'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
#             'ppthr': None, 'afthr': '0.01', 'aftag': None, 'bam': True,
#             'MQthr': None, 'BQthr': None}
#     # Run and Tests
#     with pytest.raises(SystemExit) as e:
#         assert main_novoCaller(args)
#     assert '\nERROR in variant parsing: novoAF allele frequency tag in INFO field is in the wrong format for variant:\n' in str(e.value)
#     assert '16805'
#     # Clean
#     os.remove('tests/files/main_test.out')
# #end def
#
# def test_run_novoCaller_bam_missing_aftag():
#     ''' '''
#     # Variables
#     args = {'inputfile': 'tests/files/input_novoCaller_BAM_annot_missing_novoAF.vcf', 'outputfile': 'tests/files/main_test.out',
#             'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
#             'ppthr': None, 'afthr': '0.01', 'aftag': None, 'bam': True,
#             'MQthr': None, 'BQthr': None}
#     # Run and Tests
#     with pytest.raises(SystemExit) as e:
#         assert main_novoCaller(args)
#     assert '\nERROR in variant parsing: novoAF allele frequency tag in INFO field is missing for variant:\n' in str(e.value)
#     assert '14699'
#     # Clean
#     os.remove('tests/files/main_test.out')
# #end def

def test_run_novoCaller_bam_missing_trio():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_annot_missing_novoAF.vcf', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam_missingOne.tsv',
            'ppthr': None, 'afthr': '0.01', 'aftag': 'novoAF', 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': 0.01, 'ADthr': 3, 'verbose': None}
    # Run and Tests
    with pytest.raises(SystemExit) as e:
        assert main_novoCaller(args)
    assert '\nERROR in BAMs info file for trio: missing information for some family member\n' == str(e.value)
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_run_novoCaller_bam_rerun_noUNR_NORSTR():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_novoCaller_BAM_noUNR_NORSTR.out', 'outputfile': 'tests/files/main_test.out',
            'unrelatedfiles':'tests/files/unrelated_bam.tsv', 'triofiles':'tests/files/trio_bam.tsv',
            'ppthr': None, 'afthr': '0.01', 'aftag': 'novoAF', 'bam': True,
            'MQthr': None, 'BQthr': None, 'afthr_unrelated': 0.01, 'ADthr': 3, 'verbose': None}
    # Run and Tests
    with pytest.raises(ValueError) as e:
        assert main_novoCaller(args)
    assert '\nERROR in variant FORMAT field, RSTR tag is missing\n' == str(e.value)
    # Clean
    os.remove('tests/files/main_test.out')
#end def
