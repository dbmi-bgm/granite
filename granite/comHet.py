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
    ''' object to extend a Variant object to store additional
    information on genes, transcripts and compound heterozygous pairs '''

    def __init__(self, vnt_obj, i):
        ''' initialize VariantHet object '''
        self.comHet = []
        self.vnt_obj = vnt_obj
        self.ENST_dict = {} # {ENSG: ENST_set, ...}
        self.i = i # variant index
    #end def

    def add_ENST(self, ENSG, ENST_set):
        ''' add gene and its associated transcripts information '''
        self.ENST_dict.setdefault(ENSG, ENST_set)
    #end def

    def add_pair(self, vntHet_obj, ENSG):
        ''' add information for compound heterozygous pair with vntHet_obj '''
        comHet_pair = [ENSG]
        common_ENST = self.ENST_dict[ENSG].intersection(vntHet_obj.ENST_dict[ENSG])
        if common_ENST:
            comHet_pair.append('&'.join(sorted(common_ENST)))
        else: comHet_pair.append('')
        #end if
        comHet_pair.append('{0}:{1}{2}>{3}'.format(vntHet_obj.vnt_obj.CHROM,
                                  vntHet_obj.vnt_obj.POS,
                                  vntHet_obj.vnt_obj.REF,
                                  vntHet_obj.vnt_obj.ALT))
        self.comHet.append('|'.join(comHet_pair))
    #end def

    def to_string(self):
        ''' return variant as a string after adding comHet information to INFO field '''
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
def is_comHet(vntHet_obj_1, vntHet_obj_2, ID_list, allow_undef=False):
    ''' check genotypes combination for parents if available and
    refine the assignment for the pair as a compound heterozygous or not '''
    for ID in ID_list[1:]:
        GT_1 = vntHet_obj_1.vnt_obj.get_genotype_value(ID, 'GT').replace('|', '/')
        GT_2 = vntHet_obj_2.vnt_obj.get_genotype_value(ID, 'GT').replace('|', '/')
        if GT_1 == '1/1' or GT_2 == '1/1':
            return False
        elif GT_1 in ['0/1', '1/0'] and GT_2 in ['0/1', '1/0']:
            return False
        #end if
        if not allow_undef:
            if GT_1 == './.' or GT_2 == './.':
                return False
            elif GT_1 == '0/0' and GT_2 == '0/0': # this could be a potential de novo
                                                  # in a compound het
                return False
            #end if
        #end if
    #end for
    return True
#end def

#################################################################
#    runner
#################################################################
def main(args):
    ''' '''
    # Variables
    VEPtag = args['VEPtag'] if args['VEPtag'] else 'VEP'
    allow_undef = True if args['allow_undef'] else False
    comHet_def = '##INFO=<ID=comHet,Number=3,Type=String,Description="Putative compound heterozygous pairs. Format:\'ENSG_ID|ENST_ID|VARIANT\'">'

    # Buffers
    fo = open(args['outputfile'], 'w')

    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf(args['inputfile'])

    # Add definition to header
    vcf_obj.header.add_tag_definition(comHet_def, 'INFO')

    # Writing header
    fo.write(vcf_obj.header.definitions)
    fo.write(vcf_obj.header.columns)

    # Data structures
    ENSG_dict = {} # {ENSG: [vntHet_obj1, vntHet_obj2], ...}
    ENST_dict_tmp = {} # {ENSG: ENST_set, ...}
    vntHet_set = set() # variants to write in output -> {(i, vntHet_obj), ...}
                       # i is to track and keep variants order as they are read from input

    # Get idx for ENST and ENSG
    ENSG_idx = vcf_obj.header.get_tag_field_idx(VEPtag, 'Gene')
    ENST_idx = vcf_obj.header.get_tag_field_idx(VEPtag, 'Feature')

    # Get trio IDs
    if len(args['trio']) > 3:
        sys.exit('\nERROR in parsing arguments: too many sample IDs provided for trio\n')
    #end if
    ID_list = args['trio'] # [proband_ID, parent_ID, parent_ID]

    # Reading variants
    analyzed = 0
    for c, vnt_obj in enumerate(vcf_obj.parse_variants(args['inputfile'])):
        sys.stderr.write('\rAnalyzing variant... ' + str(c + 1))
        sys.stderr.flush()

        # # Check if chromosome is canonical and in valid format
        # if not check_chrom(vnt_obj.CHROM):
        #     continue
        # #end if
        analyzed += 1

        # Reset data structures
        ENST_dict_tmp = {}

        # Check proband_ID genotype
        if vnt_obj.get_genotype_value(ID_list[0], 'GT').replace('|', '/') not in ['0/1', '1/0']:
            continue # go next if is not 0/1
        #end if

        # Get transcripts and genes information from VEP
        ENSG_list = VEP_field(vnt_obj, ENSG_idx, VEPtag)
        ENST_list = VEP_field(vnt_obj, ENST_idx, VEPtag)

        # Assign transcripts to genes
        for ENSG, ENST in zip(ENSG_list, ENST_list):
            if ENSG and ENST:
                ENST_dict_tmp.setdefault(ENSG, set())
                ENST_dict_tmp[ENSG].add(ENST)
            #end if
        #end for

        # Assign variants to genes if VEP
        if ENST_dict_tmp:
            vntHet_obj = VariantHet(vnt_obj, c)
            # Assign variant to genes and update transcripts for variant
            for ENSG, ENST_set in ENST_dict_tmp.items():
                ENSG_dict.setdefault(ENSG, [])
                ENSG_dict[ENSG].append(vntHet_obj)
                vntHet_obj.add_ENST(ENSG, ENST_set)
            #end for
        #end if
    #end for

    # Pairing variants
    sys.stderr.write('\n')
    n = len(ENSG_dict)
    for n_i, (ENSG, vntHet_list) in enumerate(ENSG_dict.items()):
        sys.stderr.write('\rPairing variants... {:.0f}%'.format(float(n_i)/n*100))
        sys.stderr.flush()
        p, l = 0, len(vntHet_list)
        while p < l:
            vntHet_obj = vntHet_list[p]
            for i, vntHet_obj_i in enumerate(vntHet_list):
                if i != p:
                    # if parents information,
                    # check genotypes to confirm is compound het or not
                    if is_comHet(vntHet_obj, vntHet_obj_i, ID_list, allow_undef):
                        vntHet_obj.add_pair(vntHet_obj_i, ENSG)
                        # Add vntHet to set to write since there is at least one pair
                        vntHet_set.add((vntHet_obj.i, vntHet_obj))
                    #end if
                #end if
            #end for
            p += 1
        #end while
    #end for
    sys.stderr.write('\rPairing variant... {0}%'.format(100))
    sys.stderr.flush()

    # Writing output
    sys.stderr.write('\n\n...Writing results for ' + str(analyzed) + ' analyzed variants out of ' + str(c + 1) + ' total variants\n')
    sys.stderr.flush()

    # Order and write variants to output file
    for _, vntHet_obj in sorted(vntHet_set, key=lambda x: x[0]):
        fo.write(vntHet_obj.to_string())
    #end for
#end def


#################################################################
#
#    MAIN
#
#################################################################
if __name__ == "__main__":

    main()

#end if
