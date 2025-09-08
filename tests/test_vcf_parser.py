#################################################################
#   Libraries
#################################################################
import sys, os
import pytest
from granite.lib import vcf_parser


#################################################################
#   Tests
#################################################################
def test_Header__get_tag_field_idx():
    ''' '''
    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf('tests/files/input_vcf_parser.vcf')
    # Tests
    assert vcf_obj.header.get_tag_field_idx('VEP', 'Feature_type') == 2
    assert vcf_obj.header.get_tag_field_idx('SpliceAI', 'MAXDS') == 0
    assert vcf_obj.header.get_tag_field_idx('CLINVAR', 'CLNSIG') == 1
#end def

def test_Header__check_tag_definition():
    ''' '''
    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf('tests/files/input_vcf_parser.vcf')
    # Tests
    assert vcf_obj.header.check_tag_definition('PID', 'FORMAT') == ('PID', 0)
    assert vcf_obj.header.check_tag_definition('AN', 'INFO') == ('AN', 0)
    assert vcf_obj.header.check_tag_definition('Consequence') == ('VEP', 3)
    assert vcf_obj.header.check_tag_definition('SYMBOL') == ('VEP', 4)
    assert vcf_obj.header.check_tag_definition('VEP') == ('VEP', 0)
    assert vcf_obj.header.check_tag_definition('ALLELEID') == ('CLINVAR', 0)
    assert vcf_obj.header.check_tag_definition('CLNSIG') == ('CLINVAR', 1)
    assert vcf_obj.header.check_tag_definition('MAXDS') == ('SpliceAI', 0)
    assert vcf_obj.header.check_tag_definition('NEGATIVE_TRAIN_SITE') == ('NEGATIVE_TRAIN_SITE', 0)
#end def

def test_Header__attributes():
    ''' '''
    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf('tests/files/input_vcf_parser.vcf')
    # Tests
    assert vcf_obj.header.columns == '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\tHG003\tHG004\n'
    assert vcf_obj.header.IDs_genotypes == ['HG002', 'HG003', 'HG004']
#end def

def test_Header__add_tag_definition():
    ''' '''
    # Variables
    FORMAT_tag = '##FORMAT=<ID=TEST_TAG,Number=1,Type=Test,Description="I\'m a tag to test the add_tag_definition applied to FORMAT">'
    INFO_tag = '##INFO=<ID=TestTag,Number=2,Type=Float,Description="I\'m a tag to test the add_tag_definition applied to INFO. Format:\'XX|YY\' ">'
    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf('tests/files/input_vcf_parser.vcf')
    # Open buffer
    fo = open('tests/files/main_test.out', 'w')
    # Add tags
    vcf_obj.header.add_tag_definition(FORMAT_tag, 'FORMAT')
    vcf_obj.header.add_tag_definition(INFO_tag, 'INFO')
    # Write to buffer
    fo.write(vcf_obj.header.definitions)
    fo.write(vcf_obj.header.columns)
    for vnt_obj in vcf_obj.parse_variants():
        fo.write(vnt_obj.to_string())
    #end for
    # Close buffer
    fo.close()
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_vcf_parser_TAG.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def

def test_read_vcf():
    ''' '''
    vcflines = vcf_parser.Vcf.read_vcf('tests/files/input_vcf_parser.vcf')
    assert next(vcflines) == '##fileformat=VCFv4.1\n'
    assert next(vcflines) == '##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">\n'
#end def

def test_read_vcf_gz():
    ''' '''
    vcflines = vcf_parser.Vcf.read_vcf('tests/files/input_vcf_parser.vcf')
    vcflines_gz = vcf_parser.Vcf.read_vcf('tests/files/input_vcf_parser.vcf.gz')
    assert next(vcflines) == next(vcflines_gz)
    assert next(vcflines) == next(vcflines_gz)
    assert next(vcflines) == next(vcflines_gz)
    assert next(vcflines) == next(vcflines_gz)
#end def

def test_write_header():
    ''' '''
    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf('tests/files/input_vcf_parser.vcf')
    # Write
    with open('tests/files/test_write_header.out', 'w') as fo:
        vcf_obj.write_definitions(fo)
    #end with
    # Test
    with open('tests/files/test_write_header.out', 'r') as f:
        output_headers = f.read()
    #end with
    input_headers = vcf_obj.header.definitions
    assert output_headers == input_headers
    # Clean up
    os.remove('tests/files/test_write_header.out')
#end def

def test_write_variant():
    ''' '''
    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf('tests/files/input_vcf_parser.vcf')
    # Write
    with open('tests/files/test_write_variants.out', 'w') as fo:
        for rec in vcf_obj.parse_variants():
            vcf_obj.write_variant(fo, rec)
        #end for
    #end with
    # Test
    recs = vcf_obj.parse_variants()
    with open('tests/files/test_write_variants.out', 'r') as f:
        for line in f:
            assert line == next(recs).to_string()
        #end for
    #end with
    # Clean up
    os.remove('tests/files/test_write_variants.out')
#end def

def test_minimal_vcf():
    ''' '''
    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf('tests/files/input_vcf_parser_minimal.vcf')
    # Writing test file
    with open('tests/files/test_write_variants.out', 'w') as fo:
        vcf_obj.write_header(fo)
        for rec in vcf_obj.parse_variants():
            assert rec.FORMAT == ''
            assert rec.GENOTYPES == {}
            vcf_obj.write_variant(fo, rec)
        #end for
    #end with
    # Tests
    assert vcf_obj.header.columns.rstrip().split('\t') == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    assert vcf_obj.header.IDs_genotypes == []
    assert [row for row in open('tests/files/test_write_variants.out')] == [row for row in open('tests/files/input_vcf_parser_minimal.vcf')]
    # Clean up
    os.remove('tests/files/test_write_variants.out')
#end def

def test_minimal_flag_emptyINFO():
    ''' '''
    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf('tests/files/input_vcf_parser_flag.vcf')
    # Parsing variants
    for i, rec in enumerate(vcf_obj.parse_variants()):
        if i == 0:
            rec.add_tag_info('TAG1=VALUE1')
            assert rec.INFO == 'TAG1=VALUE1'
        elif i == 1:
            assert rec.get_tag_value('FLAG1', is_flag=True) == True
            assert rec.get_tag_value('FLAG2', is_flag=True) == False
            assert rec.get_tag_value('AC', is_flag=True) == True
            assert rec.get_tag_value('AC', is_flag=False) == '1'
        elif i == 2:
            rec.remove_tag_info('FLAG1')
            assert rec.INFO == 'AC=1;FLAG2'
            rec.remove_tag_info('FLAG2')
            assert rec.INFO == 'AC=1'
            rec.remove_tag_info('AC')
            assert rec.INFO == '.'
            with pytest.raises(vcf_parser.MissingTag) as e:
                assert rec.get_tag_value('AC')
            assert '\nERROR in variant INFO field, AC tag is missing\n' == str(e.value)
            rec.add_tag_info('AC')
            assert rec.get_tag_value('AC', is_flag=True) == True
        #end if
    #end for
#end def
