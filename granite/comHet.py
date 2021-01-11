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
import json
# shared_functions as *
from granite.lib.shared_functions import *
# vcf_parser
from granite.lib import vcf_parser
# shared_vars
from granite.lib.shared_vars import VEP_encode
from granite.lib.shared_vars import DStags


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
        self.ENSG_IMPCT_dict = {} # {ENSG: IMPACT, ...}
        self.ENST_IMPCT_dict = {} # {ENST: IMPACT, ...}
        self.i = i # variant index
        self.SpAI = ''
        self.CLINVAR = ''
    #end def

    def add_SpAI(self, SpAI_val, SpAI_thr=0.8):
        ''' add SpliceAI, encoded as S or s based on SpAI_thr '''
        if SpAI_val >= SpAI_thr: self.SpAI = 'S'
        else: self.SpAI = 's'
        #end if
    #end def

    def add_CLINVAR(self, CLNSIG_val, CLNSIG_encode):
        ''' add CLINVAR, encoded as C or c based on CLNSIG value '''
        for clnsig, encd in CLNSIG_encode:
            if clnsig in CLNSIG_val:
                self.CLINVAR = encd
                break
            #end if
        #end for
    #end def

    def add_ENST(self, ENSG, ENST_set):
        ''' add gene and its associated transcripts information '''
        self.ENST_dict.setdefault(ENSG, ENST_set)
    #end def

    def add_ENST_IMPCT(self, ENST, IMPCT):
        ''' add transcript and its associated IMPACT,
        IMPACT encoded as number '''
        self.ENST_IMPCT_dict.setdefault(ENST, IMPCT)
    #end def

    def add_ENSG_IMPCT(self, ENSG, IMPCT_set):
        ''' add gene and worst IMPACT in its associated transcripts,
        IMPACT encoded as number '''
        self.ENSG_IMPCT_dict.setdefault(ENSG, sorted(list(IMPCT_set))[0])
    #end def

    def get_gene_impct(self, ENSG):
        ''' return combined impact information for gene '''
        return self.ENSG_IMPCT_dict[ENSG], self.SpAI + self.CLINVAR
    #end def

    def get_trscrpt_impct(self, ENST):
        ''' return combined impact information for transcript '''
        return self.ENST_IMPCT_dict[ENST], self.SpAI + self.CLINVAR
    #end def

    def map_impct(self, impct, impct_, IMPCT_decode, test):
        ''' '''
        if test:
            return IMPCT_decode[impct[0]] + impct[1] + '/' + IMPCT_decode[impct_[0]] + impct_[1]
        else:
            STRONG_set = {'H', 'M', 'S', 'C'}
            impct_set = set(list(IMPCT_decode[impct[0]] + impct[1]))
            impct_set_ = set(list(IMPCT_decode[impct_[0]] + impct_[1]))
            if impct_set.intersection(STRONG_set) and \
               impct_set_.intersection(STRONG_set):
                return 'STRONG_PAIR'
            elif impct_set.intersection(STRONG_set) or \
                 impct_set_.intersection(STRONG_set):
                return 'MEDIUM_PAIR'
            else:
                return 'WEAK_PAIR'
            #end if
        #end if
    #end def

    def worse_ENST(self, vntHet_obj, common_ENST):
        ''' return the shared transcript with worse IMPACT '''
        score, tmp_score, tmp_ENST = 0, 0, ''
        for ENST in common_ENST:
            score = self.ENST_IMPCT_dict[ENST] + vntHet_obj.ENST_IMPCT_dict[ENST]
            if not tmp_score:
                tmp_ENST = ENST
                tmp_score = score
            else:
                if score < tmp_score:
                    tmp_ENST = ENST
                #end if
            #end if
        #end for
        return tmp_ENST
    #end def

    def add_pair(self, vntHet_obj, ENSG, phase, sep, is_impct, IMPCT_decode, test):
        ''' add information for compound heterozygous pair with vntHet_obj '''
        comHet_pair = [phase, ENSG]
        common_ENST = self.ENST_dict[ENSG].intersection(vntHet_obj.ENST_dict[ENSG])
        if common_ENST:
            comHet_pair.append(sep.join(sorted(common_ENST)))
        else: comHet_pair.append('')
        #end if
        if is_impct:
            # add gene impact
            impct, impct_ = sorted([self.get_gene_impct(ENSG), vntHet_obj.get_gene_impct(ENSG)])
            comHet_pair.append(self.map_impct(impct, impct_, IMPCT_decode, test))
            # add transcript impact
            if common_ENST:
                ENST = self.worse_ENST(vntHet_obj, sorted(common_ENST))
                impct, impct_ = sorted([self.get_trscrpt_impct(ENST), vntHet_obj.get_trscrpt_impct(ENST)])
                comHet_pair.append(self.map_impct(impct, impct_, IMPCT_decode, test))
            else: comHet_pair.append('')
            #end if
        #end if
        comHet_pair.append(vntHet_obj.vnt_obj.repr())
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
            #end if
        #end if
    #end for
    return True
#end def

def phase(vntHet_obj_1, vntHet_obj_2, ID_list):
    ''' check genotypes combination for parents if available and
    refine the assignment for the pair as Phased or Unphased '''
    if len(ID_list[1:]) < 2:
        return 'Unphased'
    #end if
    for ID in ID_list[1:]:
        GT_1 = vntHet_obj_1.vnt_obj.get_genotype_value(ID, 'GT').replace('|', '/')
        GT_2 = vntHet_obj_2.vnt_obj.get_genotype_value(ID, 'GT').replace('|', '/')
        if GT_1 == '0/0' and GT_2 == '0/0': # this could be a potential de novo
                                            # in a compound het
                return 'Unphased'
        elif GT_1 == './.' or GT_2 == './.':
                return 'Unphased'
        #end if
    #end for
    return 'Phased'
#end def

def encode_IMPACT(IMPCT_list, VEP_encode, IMPCT_encode, is_IMPACT, sep):
    ''' return IMPCT_list encoded as number '''
    IMPCT_encoded = []
    if is_IMPACT: # encode IMPACT
        IMPCT_encoded = [IMPCT_encode[IMPCT] for IMPCT in IMPCT_list]
    else: # encode terms from Consequence
        for terms in IMPCT_list:
            impct_set = set()
            for term in terms.split(sep):
                try: impct_set.add(VEP_encode[term])
                except Exception: impct_set.add(VEP_encode['MODIFIER'])
                #end try
            #end for
            IMPCT_encoded.append(sorted(list(impct_set))[0])
        #end for
    #end if
    return IMPCT_encoded
#end def

def update_stats(vntHet_obj, stat_dict, sep, is_impct):
    ''' '''
    try: val_get = vntHet_obj.vnt_obj.get_tag_value('comHet')
    except Exception: return
    #end try
    var_repr = vntHet_obj.vnt_obj.repr()
    var_phase = 'Unphased'
    for cmpHet in val_get.split(','):
        if is_impct:
            PHASE, ENSG_ID, ENST_ID, IMPCT_G, IMPCT_T, VARIANT = cmpHet.split('|')
        else: PHASE, ENSG_ID, ENST_ID, VARIANT = cmpHet.split('|')
        #end if
        if PHASE == 'Phased': var_phase = PHASE
        #end if
        # global
        if (var_repr, VARIANT) not in stat_dict['pairs']['pairs_set']:
            # increase count
            stat_dict['pairs'][PHASE] += 1
            # add pair to set
            stat_dict['pairs']['pairs_set'].add((var_repr, VARIANT))
            stat_dict['pairs']['pairs_set'].add((VARIANT, var_repr))
        #end if
        # by gene
        stat_dict['genes'].setdefault(ENSG_ID, {
                                            'vntHet_set': {},
                                            'pairs_set': set(),
                                            'transcript_set': {},
                                            'Phased': 0,
                                            'Unphased': 0
                                            })
        stat_dict['genes'][ENSG_ID]['vntHet_set'].setdefault(vntHet_obj, 0)
        if PHASE == 'Phased':
            stat_dict['genes'][ENSG_ID]['vntHet_set'][vntHet_obj] = 1
        #end if
        if (var_repr, VARIANT) not in stat_dict['genes'][ENSG_ID]['pairs_set']:
            # increase count
            stat_dict['genes'][ENSG_ID][PHASE] += 1
            # add pair to set
            stat_dict['genes'][ENSG_ID]['pairs_set'].add((var_repr, VARIANT))
            stat_dict['genes'][ENSG_ID]['pairs_set'].add((VARIANT, var_repr))
        #end if
        # by transcripts
        if ENST_ID:
            for trscrpt in ENST_ID.split(sep):
                stat_dict['trscrpts'].setdefault(trscrpt, {
                                                    'gene': ENSG_ID,
                                                    'vntHet_set': {},
                                                    'pairs_set': set(),
                                                    'Phased': 0,
                                                    'Unphased': 0
                                                    })
                stat_dict['genes'][ENSG_ID]['transcript_set'].setdefault(trscrpt, 0)
                stat_dict['trscrpts'][trscrpt]['vntHet_set'].setdefault(vntHet_obj, 0)
                if PHASE == 'Phased':
                    stat_dict['genes'][ENSG_ID]['transcript_set'][trscrpt] = 1
                    stat_dict['trscrpts'][trscrpt]['vntHet_set'][vntHet_obj] = 1
                #end if
                if (var_repr, VARIANT) not in stat_dict['trscrpts'][trscrpt]['pairs_set']:
                    # increase count
                    stat_dict['trscrpts'][trscrpt][PHASE] += 1
                    # add pair to set
                    stat_dict['trscrpts'][trscrpt]['pairs_set'].add((var_repr, VARIANT))
                    stat_dict['trscrpts'][trscrpt]['pairs_set'].add((VARIANT, var_repr))
                #end if
            #end for
        #end if
        # by impact
        if is_impct:
            stat_dict['impact'].setdefault(IMPCT_G, {
                                                'vntHet_set': {},
                                                'pairs_set': set(),
                                                'gene_set': {},
                                                'transcript_set': {},
                                                'Phased': 0,
                                                'Unphased': 0
                                                    })
            stat_dict['impact'][IMPCT_G]['vntHet_set'].setdefault(vntHet_obj, 0)
            stat_dict['impact'][IMPCT_G]['gene_set'].setdefault(ENSG_ID, 0)
            if ENST_ID:
                for trscrpt in ENST_ID.split(sep):
                    stat_dict['impact'][IMPCT_G]['transcript_set'].setdefault(trscrpt, 0)
                    if PHASE == 'Phased':
                        stat_dict['impact'][IMPCT_G]['transcript_set'][trscrpt] = 1
                    #end if
                #end for
            #end if
            if PHASE == 'Phased':
                stat_dict['impact'][IMPCT_G]['vntHet_set'][vntHet_obj] = 1
                stat_dict['impact'][IMPCT_G]['gene_set'][ENSG_ID] = 1
            #end if
            if (var_repr, VARIANT) not in stat_dict['impact'][IMPCT_G]['pairs_set']:
                # increase count
                stat_dict['impact'][IMPCT_G][PHASE] += 1
                # add pair to set
                stat_dict['impact'][IMPCT_G]['pairs_set'].add((var_repr, VARIANT))
                stat_dict['impact'][IMPCT_G]['pairs_set'].add((VARIANT, var_repr))
            #end if
        #end if
    #end for
    stat_dict['vnts'][var_phase] += 1
#end def

def to_json(stat_dict, is_impct):
    ''' '''
    stat_json = {}
    phased_genes, total_genes = 0, 0
    phased_trscrpts, total_trscrpts = 0, 0
    # phased genes
    for ENSG_ID in stat_dict['genes']:
        if stat_dict['genes'][ENSG_ID]['Phased']:
            phased_genes += 1
        #end if
        total_genes += 1
    #end for
    # phased transcripts
    for trscrpt in stat_dict['trscrpts']:
        if stat_dict['trscrpts'][trscrpt]['Phased']:
            phased_trscrpts += 1
        #end if
        total_trscrpts += 1
    #end for
    # global stats
    stat_json.setdefault('general', {})
    stat_json['general'].setdefault('genes', {'phased': phased_genes, 'total': total_genes})
    stat_json['general'].setdefault('transcripts', {'phased': phased_trscrpts, 'total': total_trscrpts})
    stat_json['general'].setdefault('variants', {'phased': stat_dict['vnts']['Phased'],
                                                 'total': stat_dict['vnts']['Phased'] + stat_dict['vnts']['Unphased']})
    stat_json['general'].setdefault('pairs', {'phased': stat_dict['pairs']['Phased'],
                                              'total': stat_dict['pairs']['Phased'] + stat_dict['pairs']['Unphased']})
    # by genes
    stat_json.setdefault('by_genes', [])
    for ENSG_ID in sorted(stat_dict['genes']):
        variants_phased, variants_total = 0, 0
        for v, c in stat_dict['genes'][ENSG_ID]['vntHet_set'].items():
            if c == 1: variants_phased += 1
            #end if
            variants_total += 1
        #end for
        transcripts_phased, transcripts_total = 0, 0
        for v, c in stat_dict['genes'][ENSG_ID]['transcript_set'].items():
            if c == 1: transcripts_phased += 1
            #end if
            transcripts_total += 1
        #end for
        tmp_dict = {}
        tmp_dict.setdefault('name', ENSG_ID)
        tmp_dict.setdefault('transcripts', {'phased': transcripts_phased,
                                            'total': transcripts_total
                                            })
        tmp_dict.setdefault('variants', {'phased': variants_phased,
                                         'total': variants_total
                                        })
        tmp_dict.setdefault('pairs', {'phased': stat_dict['genes'][ENSG_ID]['Phased'],
                                      'total': stat_dict['genes'][ENSG_ID]['Phased'] + stat_dict['genes'][ENSG_ID]['Unphased']
                                     })
        # append to list
        stat_json['by_genes'].append(tmp_dict)
    #end for
    # by transcripts
    stat_json.setdefault('by_transcripts', [])
    for trscrpt in sorted(stat_dict['trscrpts']):
        variants_phased, variants_total = 0, 0
        for v, c in stat_dict['trscrpts'][trscrpt]['vntHet_set'].items():
            if c == 1: variants_phased += 1
            #end if
            variants_total += 1
        #end for
        tmp_dict = {}
        tmp_dict.setdefault('name', trscrpt)
        tmp_dict.setdefault('gene', stat_dict['trscrpts'][trscrpt]['gene'])
        tmp_dict.setdefault('variants', {'phased': variants_phased,
                                         'total': variants_total
                                        })
        tmp_dict.setdefault('pairs', {'phased': stat_dict['trscrpts'][trscrpt]['Phased'],
                                      'total': stat_dict['trscrpts'][trscrpt]['Phased'] + stat_dict['trscrpts'][trscrpt]['Unphased']
                                     })
        # append to list
        stat_json['by_transcripts'].append(tmp_dict)
    #end for
    # by impact
    if is_impct:
        stat_json.setdefault('by_impact', [])
        for impact in sorted(stat_dict['impact']):
            genes_phased, genes_total = 0, 0
            for v, c in stat_dict['impact'][impact]['gene_set'].items():
                if c == 1: genes_phased += 1
                #end if
                genes_total += 1
            #end for
            variants_phased, variants_total = 0, 0
            for v, c in stat_dict['impact'][impact]['vntHet_set'].items():
                if c == 1: variants_phased += 1
                #end if
                variants_total += 1
            #end for
            transcripts_phased, transcripts_total = 0, 0
            for v, c in stat_dict['impact'][impact]['transcript_set'].items():
                if c == 1: transcripts_phased += 1
                #end if
                transcripts_total += 1
            #end for
            tmp_dict = {}
            tmp_dict.setdefault('name', impact)
            tmp_dict.setdefault('genes', {'phased': genes_phased,
                                          'total': genes_total
                                         })
            tmp_dict.setdefault('transcripts', {'phased': transcripts_phased,
                                                'total': transcripts_total
                                                })
            tmp_dict.setdefault('variants', {'phased': variants_phased,
                                             'total': variants_total
                                            })
            tmp_dict.setdefault('pairs', {'phased': stat_dict['impact'][impact]['Phased'],
                                          'total': stat_dict['impact'][impact]['Phased'] + stat_dict['impact'][impact]['Unphased']
                                         })
            # append to list
            stat_json['by_impact'].append(tmp_dict)
        #end for
    #end if
    return stat_json
#end def

def print_stats(stat_json, fo, is_impct):
    ''' '''
    # global stats
    fo.write('##general stats\n')
    fo.write('#category\tphased\ttotal\n')
    fo.write('genes\t{0}\t{1}\n'.format(stat_json['general']['genes']['phased'],
                                        stat_json['general']['genes']['total']))
    fo.write('transcripts\t{0}\t{1}\n'.format(stat_json['general']['transcripts']['phased'],
                                              stat_json['general']['transcripts']['total']))
    fo.write('variants\t{0}\t{1}\n'.format(stat_json['general']['variants']['phased'],
                                           stat_json['general']['variants']['total']))
    fo.write('pairs\t{0}\t{1}\n'.format(stat_json['general']['pairs']['phased'],
                                        stat_json['general']['pairs']['total']))
    # by genes
    fo.write('\n##stats by genes\n')
    fo.write('#ENSG_ID\ttranscripts_phased\ttranscripts_total\tvariants_phased\tvariants_total\tpairs_phased\tpairs_total\n')
    for gene_dict in stat_json['by_genes']:
        fo.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(
                                            gene_dict['name'],
                                            gene_dict['transcripts']['phased'],
                                            gene_dict['transcripts']['total'],
                                            gene_dict['variants']['phased'],
                                            gene_dict['variants']['total'],
                                            gene_dict['pairs']['phased'],
                                            gene_dict['pairs']['total']
                                            ))
    #end for
    # by transcripts
    fo.write('\n##stats by transcripts\n')
    fo.write('#ENST_ID\tENSG_ID\tvariants_phased\tvariants_total\tpairs_phased\tpairs_total\n')
    for trscrpt_dict in stat_json['by_transcripts']:
        fo.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(
                                            trscrpt_dict['name'],
                                            trscrpt_dict['gene'],
                                            trscrpt_dict['variants']['phased'],
                                            trscrpt_dict['variants']['total'],
                                            trscrpt_dict['pairs']['phased'],
                                            trscrpt_dict['pairs']['total']
                                            ))
    #end for
    # by impact
    if is_impct:
        fo.write('\n##stats by impact\n')
        fo.write('#impact\tgenes_phased\tgenes_total\ttranscripts_phased\ttranscripts_total\tvariants_phased\tvariants_total\tpairs_phased\tpairs_total\n')
        for impact_dict in stat_json['by_impact']:
            fo.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(
                                                impact_dict['name'],
                                                impact_dict['genes']['phased'],
                                                impact_dict['genes']['total'],
                                                impact_dict['transcripts']['phased'],
                                                impact_dict['transcripts']['total'],
                                                impact_dict['variants']['phased'],
                                                impact_dict['variants']['total'],
                                                impact_dict['pairs']['phased'],
                                                impact_dict['pairs']['total']
                                                ))
        #end for
    #end if
#end def

#################################################################
#    runner
#################################################################
def main(args, test=False):
    ''' '''
    # Definitions
    CLNSIG_encode = [
                      # high impact
                      ('Pathogenic', 'C'), ('Likely_pathogenic', 'C'),
                      # moderate/low impact
                      ('Conflicting_interpretations', 'c'), ('Uncertain_significance', 'c'), ('risk_factor', 'c')
                    ]
    # VEP_encode = {...} -> import from shared_vars
    # DStags = {...} -> import from shared_vars
    IMPCT_encode = {'HIGH': 1, 'MODERATE': 2, 'LOW': 3, 'MODIFIER': 4}
    IMPCT_decode = {1: 'H', 2: 'M', 3: 'L', 4: 'm'}

    # Variables
    is_impct = True if args['impact'] else False
    is_IMPACT = False # True if IMPACT field is in VEP
    CLNSIGtag, CLNSIG_idx, is_CLNSIG = '', 0, False
    SpAItag_list, SpAI_idx_list, is_SpAI = [], [], False
    ENSG_idx, ENST_idx, IMPCT_idx = 0, 0, 0
    SpliceAItag = args['SpliceAItag'] # default None
    VEPtag = args['VEPtag'] if args['VEPtag'] else 'CSQ'
    sep = args['sep'] if args['sep'] else '&'
    allow_undef = True if args['allow_undef'] else False
    filter_cmpHet = True if args['filter_cmpHet'] else False
    granite_def = '##GRANITE=<ID=comHet>'
    comHet_def = '##INFO=<ID=comHet,Number=.,Type=String,Description="Putative compound heterozygous pairs. Subembedded:\'cmpHet\':Format:\'phase|gene|transcript|mate_variant\'">'
    comHet_impct_def = '##INFO=<ID=comHet,Number=.,Type=String,Description="Putative compound heterozygous pairs. Subembedded:\'cmpHet\':Format:\'phase|gene|transcript|impact_gene|impact_transcript|mate_variant\'">'
    is_verbose = True if args['verbose'] else False

    # Buffers
    fo = open(args['outputfile'], 'w')

    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf(args['inputfile'])

    # Add definition to header
    if is_impct:
        vcf_obj.header.add_tag_definition(granite_def + '\n' + comHet_impct_def, 'INFO')
    else:
        vcf_obj.header.add_tag_definition(granite_def + '\n' + comHet_def, 'INFO')
    #end if

    # Writing header
    vcf_obj.write_header(fo)

    # Data structures
    stat_dict = {'genes': {},
                 'trscrpts': {},
                 'pairs': {
                    'pairs_set': set(),
                    'Phased': 0,
                    'Unphased': 0
                    },
                 'vnts': {
                    'Phased': 0,
                    'Unphased': 0
                    },
                 'impact' : {}
                }
    ENSG_dict = {} # {ENSG: [vntHet_obj1, vntHet_obj2], ...}
    ENST_dict_tmp = {} # {ENSG: ENST_set, ...}
    vntHet_set = set() # variants to write in output -> {(i, vntHet_obj), ...}
                       # i is to track and keep variants order as they are read from input

    # Get idx for ENST and ENSG
    ENSG_idx = vcf_obj.header.get_tag_field_idx(VEPtag, 'Gene')
    ENST_idx = vcf_obj.header.get_tag_field_idx(VEPtag, 'Feature')

    # Get idx for SpliceAI, CLINVAR and VEP IMPACT
    if is_impct:
        try:
            IMPCT_idx = vcf_obj.header.get_tag_field_idx(VEPtag, 'IMPACT')
            is_IMPACT = True
        except Exception: # missing IMPACT, IMPCT_idx will point to Consequence
            try:
                IMPCT_idx = vcf_obj.header.get_tag_field_idx(VEPtag, 'Consequence')
            except Exception:
                sys.exit('\nERROR in VCF structure: either IMPACT or Consequence field in VEP is necessary to assign "--impact"\n')
            #end try
        #end try
        try:
            if SpliceAItag: # single tag has been specified
                tag, idx = vcf_obj.header.check_tag_definition(SpliceAItag)
                SpAItag_list.append(tag)
                SpAI_idx_list.append(idx)
            else: # search for delta scores as default
                for DStag in DStags:
                    tag, idx = vcf_obj.header.check_tag_definition(DStag)
                    SpAItag_list.append(tag)
                    SpAI_idx_list.append(idx)
                #end for
            #end if
            is_SpAI = True
        except Exception: is_SpAI = False
        #end try
        try:
            CLNSIGtag, CLNSIG_idx = vcf_obj.header.check_tag_definition('CLNSIG')
            is_CLNSIG = True
        except Exception: is_CLNSIG = False
        #end try
    #end if

    # Get trio IDs
    if len(args['trio']) > 3:
        sys.exit('\nERROR in parsing arguments: too many sample IDs provided for trio\n')
    #end if
    ID_list = args['trio'] # [proband_ID, parent_ID, parent_ID]

    # Reading variants
    analyzed = 0
    for c, vnt_obj in enumerate(vcf_obj.parse_variants()):
        if is_verbose:
            sys.stderr.write('\rAnalyzing variant... ' + str(c + 1))
            sys.stderr.flush()
        #end if

        # # Check if chromosome is canonical and in valid format
        # if not check_chrom(vnt_obj.CHROM):
        #     continue
        # #end if
        analyzed += 1

        # Reset data structures
        ENST_dict_tmp = {}
        IMPCT_dict_tmp = {}

        # Creating VariantHet object
        vntHet_obj = VariantHet(vnt_obj, c)
        if not filter_cmpHet: # if not filter, all variants are added to vntHet_set here
                              # if filter, no variant is added here to vntHet_set,
                              # compound heterozygous variants will be added after pairing
            vntHet_set.add((vntHet_obj.i, vntHet_obj))
        #end if

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

        # Assign variant to genes if VEP
        if ENST_dict_tmp:
            # Assign variant to genes and update transcripts for variant
            for ENSG, ENST_set in ENST_dict_tmp.items():
                ENSG_dict.setdefault(ENSG, [])
                ENSG_dict[ENSG].append(vntHet_obj)
                vntHet_obj.add_ENST(ENSG, ENST_set)
            #end for
        #end if

        # Add impact information if is_impct
        if is_impct:
            # VEP
            IMPCT_list = VEP_field(vnt_obj, IMPCT_idx, VEPtag)
            IMPCT_encoded = encode_IMPACT(IMPCT_list, VEP_encode, IMPCT_encode, is_IMPACT, sep)
            for i, (ENSG, IMPCT) in enumerate(zip(ENSG_list, IMPCT_encoded)):
                if ENSG and IMPCT:
                    IMPCT_dict_tmp.setdefault(ENSG, set())
                    IMPCT_dict_tmp[ENSG].add(IMPCT)
                    vntHet_obj.add_ENST_IMPCT(ENST_list[i], IMPCT)
                #end if
            #end for
            if IMPCT_dict_tmp:
                for ENSG, IMPCT_set in IMPCT_dict_tmp.items():
                    vntHet_obj.add_ENSG_IMPCT(ENSG, IMPCT_set)
                #end for
            #end if
            # SpliceAI
            if is_SpAI:
                # if SpliceAI is within VEP
                # fetching only the first transcript
                # expected the same scores for all transcripts
                SpAI_vals = []
                for i, SpAItag in enumerate(SpAItag_list):
                    SpAI_val = get_tag_idx(vnt_obj, SpAItag, SpAI_idx_list[i])
                    # if SpliceAI is with VEP and is at the end of Format
                    # need to remove , that separate next transcript
                    try: SpAI_vals.append(float(SpAI_val.split(',')[0]))
                    except Exception:
                        break
                    #end try
                #end for
                if SpAI_vals:
                    max_SpAI_vals = max(SpAI_vals)
                    if max_SpAI_vals >= 0.2:
                        vntHet_obj.add_SpAI(max_SpAI_vals)
                    #end if
                #end if
            #end if
            # CLINVAR
            if is_CLNSIG:
                # if CLNSIG is within VEP
                # fetching only the first transcript
                # expected the same CLNSIG for all transcripts
                CLNSIG_val = get_tag_idx(vnt_obj, CLNSIGtag, CLNSIG_idx)
                if CLNSIG_val:
                    vntHet_obj.add_CLINVAR(CLNSIG_val, CLNSIG_encode)
                #end if
            #end if
        #end if
    #end for

    # Pairing variants
    sys.stderr.write('\n')
    n = len(ENSG_dict)
    for n_i, (ENSG, vntHet_list) in enumerate(ENSG_dict.items()):
        if is_verbose:
            sys.stderr.write('\rPairing variants... {:.0f}%'.format(float(n_i)/n*100))
            sys.stderr.flush()
        #end if
        p, l = 0, len(vntHet_list)
        while p < l:
            vntHet_obj = vntHet_list[p]
            for i, vntHet_obj_i in enumerate(vntHet_list):
                if i != p:
                    # if parents information,
                    # check genotypes to confirm is compound het or not
                    if is_comHet(vntHet_obj, vntHet_obj_i, ID_list, allow_undef):
                        vntHet_obj.add_pair(vntHet_obj_i, ENSG, phase(vntHet_obj, vntHet_obj_i, ID_list), sep, is_impct, IMPCT_decode, test)
                        # Add vntHet to set to write since there is at least one pair
                        vntHet_set.add((vntHet_obj.i, vntHet_obj))
                    #end if
                #end if
            #end for
            p += 1
        #end while
    #end for
    if is_verbose:
        sys.stderr.write('\rPairing variants... {0}%'.format(100))
        sys.stderr.flush()
    #end if

    # Writing output
    sys.stderr.write('\n\n...Writing results for ' + str(analyzed) + ' analyzed variants out of ' + str(c + 1) + ' total variants\n')
    sys.stderr.flush()

    # Order and write variants to output file
    for _, vntHet_obj in sorted(vntHet_set, key=lambda x: x[0]):
        fo.write(vntHet_obj.to_string())
        update_stats(vntHet_obj, stat_dict, sep, is_impct)
    #end for

    # Print summary
    fs = open(args['outputfile'] + '.summary', 'w')
    fj = open(args['outputfile'] + '.json', 'w')

    # Get stats as json
    stat_json = to_json(stat_dict, is_impct)

    # Write to file
    print_stats(stat_json, fs, is_impct)
    json.dump(stat_json, fj, indent=2, sort_keys=True)

    # Close buffers
    fo.close()
    fs.close()
    fj.close()
#end def


#################################################################
#
#    MAIN
#
#################################################################
if __name__ == "__main__":

    main()

#end if
