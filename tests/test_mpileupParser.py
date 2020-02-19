#################################################################
#   Libraries
#################################################################
import sys, os
import pytest
from granite.lib import mpileup_parser
from granite.mpileupCounts import (
                            run_mpileupParser
                            )
from granite.lib import fasta_parser


#################################################################
#   Tests
#################################################################
def test_mpileupColumn__parser_reads():
    ''' '''
    reads = ['..+4ACAC.+4ACAC.+2AC', '.,,-7tttttgtt', ',$.,,^>.',
             ',.,.*.**.^*.', ',..A$..,Tn.t',
             'A-1N', ',,....*-1A.-1A.^g.-1A']
    parsed_reads = [['.', '.+ACAC', '.+ACAC', '.+AC'], ['.', ',', ',-tttttgt', 't'],
                    [',$', '.', ',', ',', '^.'], [',', '.', ',', '.', '*', '.', '*', '*', '.', '^.'],
                    [',', '.', '.', 'A$', '.', '.', ',', 'T', 'n', '.', 't'], ['A-N'],
                    [',', ',', '.', '.', '.', '.', '*-A', '.-A', '.', '^.-A']]
    parsed_reads_basic = [['.', '.+', '.+', '.+'], ['.', ',', ',-', 't'],
                          [',', '.', ',', ',', '.'], [',', '.', ',', '.', '*', '.', '*', '*', '.', '.'],
                          [',', '.', '.', 'A', '.', '.', ',', 'T', 'n', '.', 't'], ['A-'],
                          [',', ',', '.', '.', '.', '.', '*-', '.-', '.', '.-']]
    mP = mpileup_parser.mpileupParser()
    # Tests
    for i, read in enumerate(reads):
        mc = mP.mpileupColumn('chr', 0, 'ref', 'cov', read, 'BQs')
        assert mc._mpileupColumn__parser_reads(basic=False) == parsed_reads[i]
        assert mc._mpileupColumn__parser_reads() == parsed_reads_basic[i]
    #end for
#end def

def test_run_mpileupParser():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_noREF.mpileup', 'region': '1',
            'reference': 'tests/files/ref_37_chr1_50Mb.fa', 'outputfile': 'tests/files/main_test.out',
            'MQthr': None, 'BQthr': None}
    ref_dict = {} # {chr: seq, ...}
    # Fasta reader init
    handler = fasta_parser.FastaHandler()
    IT = handler.parse_generator(args['reference'])
    # Output file
    fi = open(args['inputfile'])
    fo = open(args['outputfile'], 'w')
    # Load reference
    for header, seq in IT:
        chr = header.split()[0]
        ref_dict.setdefault(chr, seq)
    #end for
    # Parsing mpileup
    run_mpileupParser(fi, fo, ref_dict)
    # Closing files
    fo.close()
    # Tests
    assert [row for row in open('tests/files/main_test.out')] == [row for row in open('tests/files/input_noREF.out')]
    # Clean
    os.remove('tests/files/main_test.out')
#end def


#################################################################
#   Errors
#################################################################
def test_run_mpileupParser_wrong_chr():
    ''' '''
    # Variables
    args = {'inputfile': 'tests/files/input_noREF_wrong_chr.mpileup', 'region': '1',
            'reference': 'tests/files/ref_37_chr1_50Mb.fa', 'outputfile': 'tests/files/main_test.out',
            'MQthr': None, 'BQthr': None}
    ref_dict = {} # {chr: seq, ...}
    # Fasta reader init
    handler = fasta_parser.FastaHandler()
    IT = handler.parse_generator(args['reference'])
    # Output file
    fi = open(args['inputfile'])
    fo = open(args['outputfile'], 'w')
    # Load reference
    for header, seq in IT:
        chr = header.split()[0]
        ref_dict.setdefault(chr, seq)
    #end for
    # Run and Tests
    with pytest.raises(SystemExit) as e:
        assert run_mpileupParser(fi, fo, ref_dict)
    assert str(e.value) == '\nERROR in reading position information: chr format (chr1) is not matching reference\n'
    # Closing files
    fo.close()
    # Clean
    os.remove('tests/files/main_test.out')
#end def
