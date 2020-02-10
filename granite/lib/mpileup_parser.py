#!/usr/bin/env python

#################################################################
#
#    mpileup_parser
#        Michele Berselli
#        Harvard Medical School
#        berselli.michele@gmail.com
#
#################################################################


#################################################################
#
#    LIBRARIES
#
#################################################################
import sys, os
import re


#################################################################
#
#    mpileupParser
#      -> mpileupColumn
#
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
                    raise ValueError('\nERROR in mpileup parser, unknown char {0} in {1}\n'
                              .format(reads[i], reads))
                #end if
            #end while
            return reads_list
        #end def

        def get_AD_noreference(self, ref):
            ''' return counts from reads '''
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
                        raise ValueError('\nERROR in mpileup parser, unknown char {0} in reads\n'
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

    #end class mpileupColumn

    def generator(self, fi):
        ''' parse mpileup format and return mpileupColumn objects as generator '''
        for line in fi:
            try: line = str(line, 'utf-8').rstrip()
            except Exception: line = line.rstrip()
            #end try
            chr, pos, ref, cov, reads, BQs = line.split()
            yield self.mpileupColumn(chr, pos, ref, cov, reads, BQs)
        #end for
    #end def

#end class mpileupParser
