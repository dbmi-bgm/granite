#!/usr/bin/env python

#################################################################
#
#    comHet
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
# shared_functions as *
from granite.lib.shared_functions import *
# vcf_parser
from granite.lib import vcf_parser


#################################################################
#
#    OBJECTS
#
#################################################################
class VariantHet(object):
    ''' '''

    def __init__(self, vnt_obj, ENSG, ENST_list):
        ''' '''
        self.comHet = []
        self.vnt_obj = vnt_obj
        self.ENSG = ENSG
        self.ENST_set = set(ENST_list)
    #end def

    def add_pair(self, vntHet_obj):
        comHet_pair = [self.ENSG]
        common_ENST = self.ENST_set.intersection(vntHet_obj.ENST_set)
        if common_ENST:
            comHet_pair.append('&'.join(common_ENST))
        else: comHet_pair.append('')
        #end if
        comHet_pair.append('{0}:{1}{2}>{3}'.format(vntHet_obj.vnt_obj.CHROM,
                                  vntHet_obj.vnt_obj.POS,
                                  vntHet_obj.vnt_obj.REF,
                                  vntHet_obj.vnt_obj.ALT))
        self.comHet.append('|'.join(comHet_pair))
    #end def

    def to_string(self):
        ''' '''
        if self.comHet:
            self.vnt_obj.add_tag_info('comHet=' + ','.join(self.comHet))
        #end if
        return self.vnt_obj.to_string()
    #end def

#end class VariantHet


#################################################################
#
#    FUNCTIONS
#
#################################################################
#################################################################
#    runner
#################################################################
def main(args):
    ''' '''
    # Variables
    VEPtag = 'VEP'

    # Buffers
    fo = open(args['outputfile'], 'w')

    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf(args['inputfile'])

    # Writing header
    fo.write(vcf_obj.header.definitions)
    fo.write(vcf_obj.header.columns)

    # Data structures
    ENSG_dict = {} # {ENSG: [vntHet_obj1, vntHet_obj2], ...}

    # Get idx for ENST and ENSG
    ENSG_idx = vcf_obj.header.get_tag_field_idx(VEPtag, 'Gene')
    ENST_idx = vcf_obj.header.get_tag_field_idx(VEPtag, 'Feature')

    # Reading variants
    analyzed = 0
    for i, vnt_obj in enumerate(vcf_obj.parse_variants(args['inputfile'])):
        sys.stderr.write('\rAnalyzing variant... ' + str(i + 1))
        sys.stderr.flush()

        # # Check if chromosome is canonical and in valid format
        # if not check_chrom(vnt_obj.CHROM):
        #     continue
        # #end if
        analyzed += 1

        # Get transcript and gene infos from VEP
        ENSG_list = VEP_field(vnt_obj, ENSG_idx, VEPtag)
        ENST_list = VEP_field(vnt_obj, ENST_idx, VEPtag)

        for ENSG in set(ENSG_list):
            # set to remove duplicated ENSG
            ENSG_dict.setdefault(ENSG, [])
            ENSG_dict[ENSG].append(VariantHet(vnt_obj, ENSG, ENST_list))
        #end for
    #end for

    # Pairing variants
    for ENSG, vntHet_list in ENSG_dict.items():
        p, l = 0, len(vntHet_list)
        while p < l:
            vntHet_obj = vntHet_list[p]
            for i, vntHet_obj_i in enumerate(vntHet_list):
                if i != p: vntHet_obj.add_pair(vntHet_obj_i)
            #end for
            fo.write(vntHet_obj.to_string())
            p += 1
        #end while
    #end for

    # Writing output
    sys.stderr.write('\n\n...Writing results for ' + str(analyzed) + ' analyzed variants out of ' + str(i + 1) + ' total variants\n')
    sys.stderr.flush()

#end def


#################################################################
#
#    MAIN
#
#################################################################
if __name__ == "__main__":

    main()

#end if

## comHet=ENSG|ENST0&ENST1&ENST2|HGNC, ENSG|ENST3|HGNC,

#'chr%s:%s %s/%s' % (CHROM, POS, REF, ALT)
