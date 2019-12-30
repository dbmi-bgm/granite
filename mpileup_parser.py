#!/usr/bin/env python

#################################################################
#
#   mpileupParser
#       requires samtools
#       contact: berselli.michele@gmail.com
#
#################################################################


#################################################################
#   Libraries
#################################################################
import sys, os
import argparse
import re, subprocess
import fasta_parser


#################################################################
#   Objects
#################################################################
class mpileupParser(object):
    ''' class to manipulate mpileup format from samtools '''
    def __init__(self):
        self.mpileupColumns = []
    #end def

    class mpileupColumn(object):
        ''' class to manipulate mpileup column at position '''
        def __init__(self, chr, pos, ref, cov, reads, BQs):
            ''' '''
            self.chr = chr
            self.pos = int(pos)
            self.ref = ref
            self.cov = cov
            self.reads = reads
            self.BQs = BQs
            self.counts = {
                    'ref_fw': 0, 'alt_fw': 0, 'ref_rv': 0, 'alt_rv': 0,
                    'ins_fw': 0, 'ins_rv': 0, 'del_fw': 0, 'del_rv': 0
                    }
        #end def

        def __parser_reads(self, basic=True):
            ''' parse the reads field and return a list of the reads.
            can provide only basic info (e.g. match/mismatch/indel),
            or be more comprehensive and return all info for each read '''
            reads_list, reads = [], self.reads
        # def test_parser_reads(reads, basic=True):
        #     reads_list = []
            i, l_i = 0, 0 # l_i is current index in reads_list
            while i < len(reads):
                if reads[i] in "ACGTNacgtn.,*":
                    reads_list.append(reads[i])
                    i += 1; l_i += 1
                elif reads[i] == '$':
                    if not basic: reads_list[l_i-1] += '$'
                    #end if
                    i += 1
                elif reads[i] == '^':
                    if not basic: reads_list.append(reads[i] + reads[i+2])
                    else: reads_list.append(reads[i+2])
                    #end if
                    i += 3; l_i += 1
                elif reads[i] in '+-':
                    re_match = re.search('[\+-](\d+)([ACGTNacgtn*]+)', reads[i:])
                    indl_len = int(re_match.group(1))
                    if not basic: reads_list[l_i-1] += reads[i] + re_match.group(2)[:indl_len]
                    else : reads_list[l_i-1] += reads[i]
                    #end if
                    i += 1 + len(re_match.group(1)) + indl_len
                else:
                    sys.exit('ERROR in mpileup parser: Unknown char {0} in {1}\n'
                              .format(reads[i], reads))
                #end if
            #end while
            return reads_list
        #end def

        def get_AD_noreference(self, ref):
            ''' '''
            encode = {
                    ref.upper(): 'ref_fw', ref.lower(): 'ref_rv',
                    'A+': 'ins_fw', 'A-': 'del_fw', 'a+': 'ins_rv', 'a-': 'del_rv',
                    'C+': 'ins_fw', 'C-': 'del_fw', 'c+': 'ins_rv', 'c-': 'del_rv',
                    'T+': 'ins_fw', 'T-': 'del_fw', 't+': 'ins_rv', 't-': 'del_rv',
                    'G+': 'ins_fw', 'G-': 'del_fw', 'g+': 'ins_rv', 'g-': 'del_rv',
                    'N+': 'ins_fw', 'N-': 'del_fw', 'n+': 'ins_rv', 'n-': 'del_rv'
                    }
            for read in self.__parser_reads():
                if '*' not in read:
                    try:
                        if read in encode:
                            self.counts[encode[read]] += 1
                        else:
                            if read.isupper(): self.counts['alt_fw'] += 1
                            else: self.counts['alt_rv'] += 1
                            #end if
                        #end if
                    except Exception:
                        sys.exit('ERROR in mpileup parser: Unknown char {0} in reads\n'
                                  .format(read))
                    #end try
                #end if
            #end for
        #end def

        def __write_AD_header(self, fo):
            ''' '''
            fo.write('#chr\tpos\tcov\t')
            fo.write('ref_fw\tref_rv\talt_fw\talt_rv\t')
            fo.write('ins_fw\tins_rv\tdel_fw\tdel_rv\n')
        #end def

        def write_AD(self, fo, header=False):
            ''' '''
            if header:
                self.__write_AD_header(fo)
            #end if
            fo.write('{0}\t{1}\t{2}\t'.format(self.chr, self.pos, self.cov))
            fo.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'
            .format(
              self.counts['ref_fw'], self.counts['ref_rv'],
              self.counts['alt_fw'], self.counts['alt_rv'],
              self.counts['ins_fw'], self.counts['ins_rv'],
              self.counts['del_fw'], self.counts['del_rv']
              ))
        #end def

    #end class

    def generator(self, fi):
        ''' '''
        for line in fi:
            try: line = str(line, 'utf-8').rstrip()
            except Exception: line = line.rstrip()
            #end try
            chr, pos, ref, cov, reads, BQs = line.split()
            yield self.mpileupColumn(chr, pos, ref, cov, reads, BQs)
        #end for
    #end def

#end class


#################################################################
#   Functions
#################################################################
def run_mpileupParser(fi, fo, ref_dict):
    ''' '''
    mP = mpileupParser()
    first = True
    for mC in mP.generator(fi):
        try: mC.get_AD_noreference(ref_dict[mC.chr][mC.pos-1])
        except Exception:
            sys.exit('ERROR in reading position information: chr format ({0}) is not matching reference\n'
                      .format(mC.chr))
        #end try
        if first:
            first = False
            mC.write_AD(fo, header=True)
        else:
            mC.write_AD(fo)
        #end if
    #end for
#end def

def main(args):
    ''' '''
    # Initialize objects and variables
    handler = fasta_parser.FastaHandler()
    ref_dict = {} # {chr#: {pos#: REF, ...}, ...}
    is_chr, chr = False, ''

    # Parsing region if available
    if args['region']:
        is_chr = True
        if ':' in args['region']:
            try:
                chr, region = args['region'].split(':')
                strt, end = map(int, region.split('-'))
                if strt >= end:
                    sys.exit('ERROR in parsing region argument: start index is larger than end index\n')
                #end if
            except Exception:
                sys.exit('ERROR in parsing region argument: the format is not recognized\n')
            #end try
        else:
            try:
                chr = args['region']
            except Exception:
                sys.exit('ERROR in parsing region argument: the format is not recognized\n')
            #end try
        #end if
    #end if

    # Building command line
    command_line = ['samtools', 'mpileup']
    if is_chr: command_line += ['-r', args['region']]
    #end if
    if args['MQthr']: command_line += ['--min-MQ', args['MQthr']]
    #end if
    if args['BQthr']: command_line += ['--min-BQ', args['MQthr']]
    #end if
    command_line += [args['inputfile']]

    # Reading reference into iterator
    IT = handler.parse_generator(args['reference'])

    # Output
    fo = open(args['outputfile'], 'w')

    # Reference headers example GRCh38 and 37
    # >chr1  AC:CM000663.2  gi:568336023  LN:248956422  rl:Chromosome  M5:6aef897c3d6ff0c78aff06ac189178dd  AS:GRCh38
    # >1 dna:chromosome chromosome:GRCh37:1:1:249250621:1

    # Loading reference
    if is_chr:
        for header, seq in IT:
            if header.split()[0] == chr:
                ref_dict = {chr: seq}
                break
            #end if
        #end for
    else:
        for header, seq in IT:
            chr = header.split()[0]
            ref_dict.setdefault(chr, seq)
        #end for
    #end if

    # Running samtools
    pipe_in = subprocess.Popen(command_line, stdout=subprocess.PIPE)

    # Parsing mpileup
    run_mpileupParser(pipe_in.stdout, fo, ref_dict)

    # Closing files
    fo.close()
#end def


#################################################################
# MAIN
#################################################################
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Use samtools to calculate statistics for pileup at each position in the specified region')

    parser.add_argument('-i', '--inputfile', help='I/O: input file', required=True)
    parser.add_argument('-o', '--outputfile', help='I/O: output file', required=True)
    parser.add_argument('-r', '--reference', help='OTHER: reference file', required=False)
    parser.add_argument('--region', help='OTHER: region to be analyzed [e.g chr1:1-10000000, 1:1-10000000, chr1, 1], chromsome name have to match the reference', required=False)
    parser.add_argument('--MQthr', help='OTHER: minimum mapping quality for an alignment to be used [0]', required=False)
    parser.add_argument('--BQthr', help='OTHER: minimum base quality for a base to be considered [13]', required=False)

    args = vars(parser.parse_args())

    main(args)


#################################################################
#   Test
#################################################################
# a = '..+4ACAC.+4ACAC.+2AC' #=> ['.', '.+ACAC', '.+ACAC', '.+AC']
# b = '.,,-7tttttgtt' #=> ['.', ',', ',-tttttgt', 't']
# c = ',$.,,^>.' #=> [',$', '.', ',', ',', '^.']
# d = ',.,.*.**.^*.' #=> [',', '.', ',', '.', '*', '.', '*', '*', '.', '^.']
# e = ',..A$..,Tn.t' #=> [',', '.', '.', 'A$', '.', '.', ',', 'T', 'n', '.', 't']
# f = 'A-1N' #=> ['A-N']
# g = ',,....*-1A.-1A.^g.-1A' # => [',', ',', '.', '.', '.', '.', '*-A', '.-A', '.', '^.-A']
#
# print(test_parser_reads(a, basic=False) == ['.', '.+ACAC', '.+ACAC', '.+AC'])
# print(test_parser_reads(a) == ['.', '.+', '.+', '.+'])
# print(test_parser_reads(b, basic=False) == ['.', ',', ',-tttttgt', 't'])
# print(test_parser_reads(b) == ['.', ',', ',-', 't'])
# print(test_parser_reads(c, basic=False) == [',$', '.', ',', ',', '^.'])
# print(test_parser_reads(c) == [',', '.', ',', ',', '.'])
# print(test_parser_reads(d, basic=False) == [',', '.', ',', '.', '*', '.', '*', '*', '.', '^.'])
# print(test_parser_reads(d) == [',', '.', ',', '.', '*', '.', '*', '*', '.', '.'])
# print(test_parser_reads(e, basic=False) == [',', '.', '.', 'A$', '.', '.', ',', 'T', 'n', '.', 't'])
# print(test_parser_reads(e) == [',', '.', '.', 'A', '.', '.', ',', 'T', 'n', '.', 't'])
# print(test_parser_reads(f, basic=False) == ['A-N'])
# print(test_parser_reads(f) == ['A-'])
# print(test_parser_reads(g, basic=False) == [',', ',', '.', '.', '.', '.', '*-A', '.-A', '.', '^.-A'])
# print(test_parser_reads(g) == [',', ',', '.', '.', '.', '.', '*-', '.-', '.', '.-'])
