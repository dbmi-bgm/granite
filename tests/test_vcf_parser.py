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
    for vnt_obj in vcf_obj.parse_variants('tests/files/input_vcf_parser.vcf'):
        fo.write(vnt_obj.to_string())
    #end for
    # Close buffer
    fo.close()
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_vcf_parser_TAG.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def
