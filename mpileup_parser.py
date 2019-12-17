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
        def __init__(self, chr, pos, ref, cov, reads, BQs, *, MQs=None):
            ''' '''
            self.chr = chr
            self.pos = int(pos)
            self.ref = ref
            self.cov = cov
            self.reads = reads
            self.BQs = BQs
            self.MQs = MQs
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
                if not '*' in read:
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

    def generator(self, fi, is_MQs=False):
        for line in fi:
            if not is_MQs:
                chr, pos, ref, cov, reads, BQs = str(line, 'utf-8').rstrip().split()
                yield self.mpileupColumn(chr, pos, ref, cov, reads, BQs)
            else:
                chr, pos, ref, cov, reads, BQs, MQs = str(line, 'utf-8').rstrip().split()
                yield self.mpileupColumn(chr, pos, ref, cov, reads, BQs, MQs=MQs)
            #end if
        #end for
    #end def

#end class


#################################################################
#   Functions
#################################################################
def main(args):

    mP = mpileupParser()
    handler = fasta_parser.FastaHandler()

    # Reading reference into generator
    IT = handler.parse_generator(args['referencefile'])

    ref_dict = {x: 'N' for x in range(10000010)}

    fo = open(args['outputfile'], 'w')

    h, s = next(IT)

    for i, s in enumerate(s[:10000009]):
        ref_dict[i] = s
    #end for

    pipe_in = subprocess.Popen(['samtools', 'mpileup', '-r', '1:1-10000000', args['inputfile']], stdout=subprocess.PIPE)
    #pipe_in = subprocess.Popen(['samtools', 'mpileup', args['inputfile']], stdout=subprocess.PIPE)

    first = True
    for mC in mP.generator(pipe_in.stdout):
        mC.get_AD_noreference(ref_dict[mC.pos])
        if first:
            first = False
            mC.write_AD(fo, header=True)
        else:
            mC.write_AD(fo)
        #end if
    #end for

    fo.close()
#end def


#################################################################
# MAIN
#################################################################
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Use samtools to calculate statistics for pileup at each position')

    parser.add_argument('-i', '--inputfile', help='I/O: input file', required=True)
    parser.add_argument('-o', '--outputfile', help='I/O: output file', required=True)
    parser.add_argument('-r', '--referencefile', help='OTHER: reference file', required=False)

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


#################################################################
#   Old
#################################################################
# def get_AD_noreference_fast(self, ref):
#     ''' '''
#     encode = {
#             ref.upper(): 'ref_fw', ref.lower(): 'ref_rv',
#             'A+': 'ins_fw', 'A-': 'del_fw', 'a+': 'ins_rv', 'a-': 'del_rv',
#             'C+': 'ins_fw', 'C-': 'del_fw', 'c+': 'ins_rv', 'c-': 'del_rv',
#             'T+': 'ins_fw', 'T-': 'del_fw', 't+': 'ins_rv', 't-': 'del_rv',
#             'G+': 'ins_fw', 'G-': 'del_fw', 'g+': 'ins_rv', 'g-': 'del_rv',
#             'N+': 'ins_fw', 'N-': 'del_fw', 'n+': 'ins_rv', 'n-': 'del_rv'
#             }
#     i, j = 0, 0
#     reads, indls, r_tmp = self.reads, [], ''
#     if '^' in reads or '$' in reads:
#         reads = re.sub(r'\^.?|\$', '', reads)
#     #end if
#     if '+' in reads or '-' in reads:
#         indls = re.findall('[\+-](\d+)[ACGTNacgtn*]+', reads)
#     #end if
#     while i < len(reads):
#         r_tmp = reads[i]
#         try:
#             r_nex = reads[i+1]
#             if r_nex in '+-':
#                 r_tmp += r_nex
#                 i += 1 + len(indls[j]) + int(indls[j]) + 1
#                 j += 1
#             else:
#                 i += 1
#             #end if
#         except Exception:
#             i += 1
#         #end try
#         if '*' not in r_tmp:
#             try:
#                 if r_tmp in encode:
#                     self.counts[encode[r_tmp]] += 1
#                 else:
#                     if r_tmp.isupper(): self.counts['alt_fw'] += 1
#                     else: self.counts['alt_rv'] += 1
#                     #end if
#                 #end if
#             except Exception:
#                 sys.exit('ERROR in mpileup parser: Unknown char {0} in reads\n'
#                           .format(read))
#             #end try
#         #end if
#     #end while
# #end def
