#################################################################
#   Libraries
#################################################################
import sys, os
import pytest
import json
from granite.lib import pedigree_parser


#################################################################
#   Tests
#################################################################
def test_pedigree_parser():
    ''' '''
    samples = ['NA12881_smpl', 'NA12879_sample', 'NA12878_sample', 'NA12877_sample',
               'NA12876_sample', 'NA12875_sample', 'NA12873_sample']
    results = {
        'NA12873' : {'parents': [], 'children': ['NA12878'],
                     'sample': 'NA12873_sample', 'gender': 'M',
                     'siblings': []},
        'NA12874' : {'parents': [], 'children': ['NA12878'],
                     'sample': None, 'gender': 'F',
                     'siblings': []},
        'NA12875' : {'parents': [], 'children': ['NA12877'],
                     'sample': 'NA12875_sample', 'gender': 'M',
                     'siblings': []},
        'NA12876' : {'parents': [], 'children': ['NA12877'],
                     'sample': 'NA12876_sample', 'gender': 'F',
                     'siblings': []},
        'NA12877' : {'parents': ['NA12875', 'NA12876'], 'children': ['NA12879', 'NA12880', 'NA12881'],
                     'sample': 'NA12877_sample', 'gender': 'M',
                     'siblings': []},
        'NA12878' : {'parents': ['NA12873', 'NA12874'], 'children': ['NA12879', 'NA12880', 'NA12881'],
                     'sample': 'NA12878_sample', 'gender': 'F',
                     'siblings': []},
        'NA12879' : {'parents': ['NA12878', 'NA12877'], 'children': [],
                     'sample': 'NA12879_sample', 'gender': 'F',
                     'siblings': ['NA12880', 'NA12881']},
        'NA12880' : {'parents': ['NA12877', 'NA12878'], 'children': [],
                     'sample': None, 'gender': 'U',
                     'siblings': ['NA12879', 'NA12881']},
        'NA12881' : {'parents': ['NA12877', 'NA12878'], 'children': [],
                     'sample': 'NA12881_smpl', 'gender': 'U',
                     'siblings': ['NA12879', 'NA12880']}
    }
    # Loading pedigree
    with open('tests/files/pedigree.json') as fi:
        pedigree = json.load(fi)
    #end with
    # Creating Pedigree object
    pedigree_obj = pedigree_parser.Pedigree(pedigree)
    # Test
    for name, member_obj in pedigree_obj.members.items():
        assert name == member_obj.name
        assert [obj.name for obj in member_obj.get_parents()] == results[name]['parents']
        assert [obj.name for obj in member_obj.get_children()] == results[name]['children']
        assert member_obj.sample == results[name]['sample']
        assert member_obj.gender == results[name]['gender']
        assert [obj.name for obj in member_obj.get_siblings()] == results[name]['siblings']
    #end for
    for sample in samples:
        sample_obj = pedigree_obj.get_member_by_sample(sample)
        assert sample_obj.sample == sample
    #end for
#end def

def test_pedigree_parser_spouses():
    ''' '''
    sample = 'NA12877_sample'
    results = {
        'spouses': ['NA12878', 'NA12896', 'NA12895'],
        'NA12878': {'NA12880', 'NA12881', 'NA12882'},
        'NA12896': {'NA12883', 'NA12884'},
        'NA12895': {'NA12879'}
    }
    # Loading pedigree
    with open('tests/files/pedigree_spouses.json') as fi:
        pedigree = json.load(fi)
    #end with
    # Creating Pedigree object
    pedigree_obj = pedigree_parser.Pedigree(pedigree)
    # Test
    sample_obj = pedigree_obj.get_member_by_sample(sample)
    NA12878_obj = pedigree_obj.get_member_by_sample('NA12878_sample')
    NA12896_obj = pedigree_obj.get_member_by_sample('NA12896_sample')
    NA12895_obj = pedigree_obj.get_member_by_sample('NA12895_sample')
    assert [obj.name for obj in sample_obj.get_spouses()] == results['spouses']
    assert set([obj.name for obj in sample_obj.common_children(NA12878_obj)]) == results['NA12878']
    assert set([obj.name for obj in sample_obj.common_children(NA12896_obj)]) == results['NA12896']
    assert set([obj.name for obj in sample_obj.common_children(NA12895_obj)]) == results['NA12895']
#end def

#################################################################
#   Errors
#################################################################
def test_pedigree_parser_missing_sample():
    ''' '''
    # Loading pedigree
    with open('tests/files/pedigree.json') as fi:
        pedigree = json.load(fi)
    #end with
    # Creating Pedigree object
    pedigree_obj = pedigree_parser.Pedigree(pedigree)
    # Test
    with pytest.raises(ValueError) as e:
        sample_obj = pedigree_obj.get_member_by_sample('NA12880_sample')
    assert str(e.value) == '\nERROR in pedigree structure, missing sample NA12880_sample in pedigree\n'
#end def

def test_pedigree_parser_missing_member_info():
    ''' '''
    member = {
        "sample_name": "NA12878_sample",
        "gender": "F",
        "parents": ["NA12873", "NA12874"]
    }
    # Loading pedigree
    with open('tests/files/pedigree.json') as fi:
        pedigree = json.load(fi)
    #end with
    # Creating Pedigree object
    pedigree_obj = pedigree_parser.Pedigree(pedigree)
    # Test
    with pytest.raises(ValueError) as e:
        pedigree_obj.add_member(member)
    assert str(e.value) == '\nERROR in pedigree structure, missing name for pedigree member\n'
#end def
