#!/usr/bin/env python

#################################################################
#
#    validateVCF
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
import numpy
import copy
# matplotlib
from matplotlib import pyplot
# shared_functions as *
from granite.lib.shared_functions import *
# vcf_parser
from granite.lib import vcf_parser
# pedigree_parser
from granite.lib import pedigree_parser
# shared_vars
from granite.lib.shared_vars import real_NA_chroms


#################################################################
#
#    FUNCTIONS
#
#################################################################
#################################################################
#    Heterozygosity stats
#################################################################
def get_stats_het(vnt_obj, stat_dict, family, NA_chroms):
    ''' extract heterozygosity information from variant using family '''
    var_type = variant_type_ext(vnt_obj.REF, vnt_obj.ALT)
    if vnt_obj.CHROM.replace('chr', '') in NA_chroms:
        return # skip unbalanced chromosomes
    #end if
    is_ID = True if vnt_obj.ID != '.' else False
    if var_type == 'snv':
        _error_het(vnt_obj, stat_dict, family, is_ID)
        _error_het_family(vnt_obj, stat_dict, family)
    #end if
#end def

def _error_het_family(vnt_obj, stat_dict, family):
    ''' calculate error model for heterozygous calls using family, MODEL '''
    # Counting children that are heterozygous
    cnt_children = 0
    for child in family['children']:
        try:
            GT_ = vnt_obj.get_genotype_value(child.sample, 'GT').replace('|', '/')
        except Exception: continue
        #end try
        if GT_ == './.' or GT_ == '1/1': return
        elif GT_ in ['0/1', '1/0']: cnt_children += 1
        #end if
    #end for
    # Parents samples
    sample_0 = family['parents'][0].sample
    sample_1 = family['parents'][1].sample
    # Update stat_dict for cnt_children
    stat_dict['error_het_family'].setdefault(cnt_children, {})
    for sample in [sample_0, sample_1]:
        stat_dict['error_het_family'][cnt_children].setdefault(sample, {
                                                            'ref': 0, # undercalled het
                                                            'hom': 0, # overcalled het
                                                            'het': 0, # properly called het
                                                            'missing': 0, # missing genotype
                                                            'total': 0
                                                        })
    #end for
    # Get trio genotypes for parents
    trio_0 = _GT_trio(vnt_obj, family['parents'][0])
    trio_1 = _GT_trio(vnt_obj, family['parents'][1])
    # Check trio genotypes combination
    if trio_0 and trio_1: # trio genotypes are complete
        if trio_0 == ['ref', {'ref'}] and trio_1[1] == {'ref', 'het'}:
            stat_dict['error_het_family'][cnt_children][sample_1][trio_1[0]] += 1
            stat_dict['error_het_family'][cnt_children][sample_1]['total'] += 1
        elif trio_1 == ['ref', {'ref'}] and trio_0[1] == {'ref', 'het'}:
            stat_dict['error_het_family'][cnt_children][sample_0][trio_0[0]] += 1
            stat_dict['error_het_family'][cnt_children][sample_0]['total'] += 1
        #end if
    #end if
#end def

def _error_het(vnt_obj, stat_dict, family, is_ID):
    ''' calculate error model for heterozygous calls using family, GENERAL '''
    # Counting children that are heterozygous or homozygous alternate
    cnt_children = 0
    for child in family['children']:
        try:
            GT_ = vnt_obj.get_genotype_value(child.sample, 'GT').replace('|', '/')
        except Exception: continue
        #end try
        if GT_ == './.': return
        # if GT_ == './.' or GT_ == '1/1': return
        elif GT_ in ['0/1', '1/0', '1/1']: cnt_children += 1
        # elif GT_ in ['0/1', '1/0']: cnt_children += 1
        #end if
    #end for
    # Parents samples
    sample_0 = family['parents'][0].sample
    sample_1 = family['parents'][1].sample
    # Update stat_dict for cnt_children
    stat_dict['error_het'].setdefault(cnt_children, {})
    for sample in [sample_0, sample_1]:
        stat_dict['error_het'][cnt_children].setdefault(sample, {
                                                    'no_spouse': {
                                                        'is_ID': {
                                                            'in_gparents': 0,
                                                            'no_gparents': 0,
                                                            'total': 0
                                                        },
                                                        'no_ID': {
                                                            'in_gparents': 0,
                                                            'no_gparents': 0,
                                                            'total': 0
                                                        }
                                                    },
                                                    'in_spouse': {
                                                        'is_ID': {
                                                            'in_gparents': 0,
                                                            'no_gparents': 0,
                                                            'total': 0
                                                        },
                                                        'no_ID': {
                                                            'in_gparents': 0,
                                                            'no_gparents': 0,
                                                            'total': 0
                                                        }
                                                    }
                                                })
    #end for
    # Get trio genotypes for parents
    trio_0 = _GT_trio(vnt_obj, family['parents'][0])
    trio_1 = _GT_trio(vnt_obj, family['parents'][1])
    key_ID = 'is_ID' if is_ID else 'no_ID'
    if trio_0 and trio_1:
        if trio_0[0] == 'het':
            key_sps = 'no_spouse' if trio_1 == ['ref', {'ref'}] else 'in_spouse'
            stat_dict['error_het'][cnt_children][sample_0][key_sps][key_ID]['total'] += 1
            if trio_0[1] != {'ref'}:
                stat_dict['error_het'][cnt_children][sample_0][key_sps][key_ID]['in_gparents'] += 1
            else:
                stat_dict['error_het'][cnt_children][sample_0][key_sps][key_ID]['no_gparents'] += 1
            #end if
        #end if
        if trio_1[0] == 'het':
            key_sps = 'no_spouse' if trio_0 == ['ref', {'ref'}] else 'in_spouse'
            stat_dict['error_het'][cnt_children][sample_1][key_sps][key_ID]['total'] += 1
            if trio_1[1] != {'ref'}:
                stat_dict['error_het'][cnt_children][sample_1][key_sps][key_ID]['in_gparents'] += 1
            else:
                stat_dict['error_het'][cnt_children][sample_1][key_sps][key_ID]['no_gparents'] += 1
            #end if
        #end if
    #end if
#end def

#################################################################
#    De novo stats
#################################################################
def get_stats_novo(vnt_obj, stat_dict, family, NA_chroms, sample_novo, novotag):
    ''' extract de novo information from variant using family '''
    var_type = variant_type_ext(vnt_obj.REF, vnt_obj.ALT)
    if vnt_obj.CHROM.replace('chr', '') in NA_chroms:
        return # skip unbalanced chromosomes
    #end if
    if var_type == 'snv':
        try:
            novoPP = float(vnt_obj.get_tag_value(novotag))
        except Exception: return
        #end try
        if novoPP > 0:
            _error_novo_family(vnt_obj, stat_dict, family, sample_novo, novoPP)
        #end if
    #end if
#end def

def _error_novo_family(vnt_obj, stat_dict, family, sample_novo, novoPP):
    ''' calculate error model for de novo calls using family '''
    # Counting children that are heterozygous or homozygous alternate
    cnt_children, is_child = 0, False
    for child in family['children']:
        if child.sample != sample_novo:
            try:
                GT_ = vnt_obj.get_genotype_value(child.sample, 'GT').replace('|', '/')
            except Exception: continue
            #end try
            if GT_ == './.': return
            elif GT_ in ['0/1', '1/0', '1/1']: cnt_children += 1
            #end if
        else: is_child = True
        #end if
    #end for
    if novoPP >= 0.9: bin = 0.9
    elif novoPP >= 0.1: bin = 0.1
    else: bin = 0
    #end if
    # Update stat_dict for cnt_children
    stat_dict['error_novo_family'].setdefault(bin, {})
    stat_dict['error_novo_family'][bin].setdefault(cnt_children, {
                                                            'in_gparents': 0,
                                                            'no_gparents': 0,
                                                            'total': 0,
                                                            'no_gparents_vnt': [],
                                                            'total_vnt': [],
                                                            'no_gparents_AD-DP': []
                                                        })
    # Check if sample_novo is child or parent
    if is_child:
        _error_novo_child(vnt_obj, stat_dict, family, bin, cnt_children, sample_novo)
    else:
        _error_novo_parent(vnt_obj, stat_dict, family, bin, cnt_children, sample_novo)
    #end if
#end def

def _error_novo_child(vnt_obj, stat_dict, family, bin, cnt_children, sample_novo):
    ''' sample to calculate error for is child '''
    # Get trio genotypes for parents
    trio_0 = _GT_trio(vnt_obj, family['parents'][0])
    trio_1 = _GT_trio(vnt_obj, family['parents'][1])
    # Check trio genotypes combination
    if trio_0 and trio_1: # trio genotypes are complete
        stat_dict['error_novo_family'][bin][cnt_children]['total'] += 1
        stat_dict['error_novo_family'][bin][cnt_children]['total_vnt'].append(vnt_obj)
        if trio_0[1] != {'ref'} or trio_1[1] != {'ref'}:
            stat_dict['error_novo_family'][bin][cnt_children]['in_gparents'] += 1
        else:
            stat_dict['error_novo_family'][bin][cnt_children]['no_gparents'] += 1
            stat_dict['error_novo_family'][bin][cnt_children]['no_gparents_vnt'].append(vnt_obj)
            # Get AD-DP ratio
            AD = sum(map(int, vnt_obj.get_genotype_value(sample_novo, 'AD').split(',')[1:]))
            DP = int(vnt_obj.get_genotype_value(sample_novo, 'DP'))
            try:
                stat_dict['error_novo_family'][bin][cnt_children]['no_gparents_AD-DP'].append(AD / DP)
            except Exception:
                stat_dict['error_novo_family'][bin][cnt_children]['no_gparents_AD-DP'].append(0)
            #end try
        #end if
    #end if
#end def

def _error_novo_parent(vnt_obj, stat_dict, family, bin, cnt_children, sample_novo):
    ''' sample to calculate error for is parent '''
    # Get trio for parent only
    for parent in family['parents']:
        if parent.sample == sample_novo:
            trio = _GT_trio(vnt_obj, parent)
            break
        #end if
    #end for
    if trio: # ignoring spouse family, checking parent only
        stat_dict['error_novo_family'][bin][cnt_children]['total'] += 1
        stat_dict['error_novo_family'][bin][cnt_children]['total_vnt'].append(vnt_obj)
        if trio[1] != {'ref'}:
            stat_dict['error_novo_family'][bin][cnt_children]['in_gparents'] += 1
        else:
            stat_dict['error_novo_family'][bin][cnt_children]['no_gparents'] += 1
            stat_dict['error_novo_family'][bin][cnt_children]['no_gparents_vnt'].append(vnt_obj)
            # Get AD-DP ratio
            AD = sum(map(int, vnt_obj.get_genotype_value(sample_novo, 'AD').split(',')[1:]))
            DP = int(vnt_obj.get_genotype_value(sample_novo, 'DP'))
            stat_dict['error_novo_family'][bin][cnt_children]['no_gparents_AD-DP'].append(AD / DP)
        #end if
    #end if
#end def

#################################################################
#    Routines for stats
#################################################################
def _GT_trio(vnt_obj, member_obj):
    ''' return genotype information for member_obj and parents,
    [GT_member_obj, set([GT_parent_0, GT_parent_1])] '''
    encode = {'0/0': 'ref', '1/1': 'hom',
              '1/0': 'het', '0/1': 'het',
              './.': 'missing'}
    GT_trio, GT_parents = [], set()
    GT = vnt_obj.get_genotype_value(member_obj.sample, 'GT').replace('|', '/')
    for parent in member_obj.get_parents():
        GT_ = vnt_obj.get_genotype_value(parent.sample, 'GT').replace('|', '/')
        if GT_ == './.': return
        else: GT_parents.add(encode[GT_])
        #end if
    #end for
    GT_trio.append(encode[GT])
    GT_trio.append(GT_parents)
    return GT_trio
#end def

#################################################################
#
#################################################################
def to_json(stat_dict_list, sample_het_list, sample_novo):
    ''' '''
    stat_json = {}
    # Error heterozygous calls, model
    if sample_het_list:
        stat_json.setdefault('autosomal heterozygous calls error, MODEL', [])
        # Calculate error
        tmp_dict = {}
        for stat_dict in stat_dict_list:
            if stat_dict['error_het_family']:
                for cnt_children in sorted(stat_dict['error_het_family']):
                    tmp_dict.setdefault(cnt_children, {
                                                        'children': cnt_children
                                                        })
                    for sample in stat_dict['error_het_family'][cnt_children]:
                        if sample in sample_het_list:
                            tmp_dict[cnt_children].setdefault(sample, {})
                            _counts = {
                                'total': stat_dict['error_het_family'][cnt_children][sample]['total']
                            }
                            for k, v in stat_dict['error_het_family'][cnt_children][sample].items():
                                # Counts
                                if k == 'ref': _counts.setdefault('undercall', v)
                                elif k == 'hom': _counts.setdefault('overcall', v)
                                elif k == 'het': _counts.setdefault('correct', v)
                                elif k == 'missing': _counts.setdefault('missing', v)
                                #end if
                            #end for
                            # Ratios
                            if _counts['total']:
                                tmp_dict[cnt_children][sample].setdefault('undercall ratio', round(_counts['undercall'] / _counts['total'], 3))
                                tmp_dict[cnt_children][sample].setdefault('overcall ratio', round(_counts['overcall'] / _counts['total'], 3))
                                tmp_dict[cnt_children][sample].setdefault('correct ratio', round(_counts['correct'] / _counts['total'], 3))
                                tmp_dict[cnt_children][sample].setdefault('missing ratio', round(_counts['missing'] / _counts['total'], 3))
                            #end if
                            tmp_dict[cnt_children][sample].setdefault('counts', _counts)
                        #end if
                    #end for
                #end for
            #end if
        #end for
        for cnt_children in sorted(tmp_dict):
            stat_json['autosomal heterozygous calls error, MODEL'].append(tmp_dict[cnt_children])
        #end for
        _extend_json_het(stat_json)
    #end if
    # Error heterozygous calls, general
    if sample_het_list:
        stat_json.setdefault('autosomal heterozygous calls error, GENERAL', [])
        # Calculate error
        tmp_dict = {}
        for sample in sample_het_list:
            tmp_dict.setdefault(sample, {
                                'sample': sample,
                                'by_siblings': {}
                                })
        #end for
        for stat_dict in stat_dict_list:
            if stat_dict['error_het']:
                for cnt_children in sorted(stat_dict['error_het']):
                    for sample in stat_dict['error_het'][cnt_children]:
                        if sample in tmp_dict:
                            tmp_dict[sample]['by_siblings'].setdefault(cnt_children, {
                                                                            'in_gparents': 0,
                                                                            'no_gparents': 0,
                                                                            'total': 0,
                                                                            'is_ID': {
                                                                                'in_gparents': 0,
                                                                                'no_gparents': 0,
                                                                                'total': 0
                                                                            },
                                                                            'no_spouse': {
                                                                                'is_ID': {
                                                                                    'in_gparents': 0,
                                                                                    'no_gparents': 0,
                                                                                    'total': 0
                                                                                },
                                                                                'in_gparents': 0,
                                                                                'no_gparents': 0,
                                                                                'total': 0
                                                                            }
                                                                        })
                            # In spouse
                            for k, v in stat_dict['error_het'][cnt_children][sample]['in_spouse']['is_ID'].items():
                                tmp_dict[sample]['by_siblings'][cnt_children][k] += v
                                tmp_dict[sample]['by_siblings'][cnt_children]['is_ID'][k] += v
                            #end for
                            for k, v in stat_dict['error_het'][cnt_children][sample]['in_spouse']['no_ID'].items():
                                tmp_dict[sample]['by_siblings'][cnt_children][k] += v
                            #end for
                            # No spouse
                            for k, v in stat_dict['error_het'][cnt_children][sample]['no_spouse']['is_ID'].items():
                                tmp_dict[sample]['by_siblings'][cnt_children][k] += v
                                tmp_dict[sample]['by_siblings'][cnt_children]['is_ID'][k] += v
                                tmp_dict[sample]['by_siblings'][cnt_children]['no_spouse'][k] += v
                                tmp_dict[sample]['by_siblings'][cnt_children]['no_spouse']['is_ID'][k] += v
                            #end for
                            for k, v in stat_dict['error_het'][cnt_children][sample]['no_spouse']['no_ID'].items():
                                tmp_dict[sample]['by_siblings'][cnt_children][k] += v
                                tmp_dict[sample]['by_siblings'][cnt_children]['no_spouse'][k] += v
                            #end for
                        #end if
                    #end for
                #end for
            #end if
        #end for
        for sample in sorted(tmp_dict):
            stat_json['autosomal heterozygous calls error, GENERAL'].append(tmp_dict[sample])
        #end for
    #end if
    # Error de novo calls
    if sample_novo:
        stat_json.setdefault('de novo calls error', [])
        for stat_dict in stat_dict_list: # only one stat_dict at most
            if stat_dict['error_novo_family']:
                tmp_dict = {
                        'novoPP': '',
                        'total variants': 0,
                        'not in grandparents': 0,
                        'not in grandparents, in siblings': {}
                        }
                for bin in sorted(stat_dict['error_novo_family'], reverse=True):
                    if bin != 0: tmp_dict['novoPP'] = '>= ' + str(bin)
                    else: tmp_dict['novoPP'] = '> ' + str(bin)
                    #end if
                    for cnt_children in sorted(stat_dict['error_novo_family'][bin]):
                        if stat_dict['error_novo_family'][bin][cnt_children]['no_gparents']:
                            tmp_dict['not in grandparents, in siblings'].setdefault(cnt_children, 0)
                        #end if
                        for k, v in stat_dict['error_novo_family'][bin][cnt_children].items():
                            if k == 'total':
                                tmp_dict['total variants'] += v
                            elif k == 'no_gparents' and v:
                                tmp_dict['not in grandparents'] += v
                                tmp_dict['not in grandparents, in siblings'][cnt_children] += v
                            #end if
                        #end for
                    #end for
                    stat_json['de novo calls error'].append(copy.deepcopy(tmp_dict))
                #end for
            #end if
        #end for
    #end if
    return stat_json
#end def

def _extend_json_het(stat_json):
    ''' '''
    # Error heterozygous calls in cnt_children range
    cnt_children = stat_json['autosomal heterozygous calls error, MODEL'][-1]['children']
    median = cnt_children // 2 + 1
    range = cnt_children * 25 // 100
    min = median - range
    max = median + range
    dict_range = {
        'children': '[' + str(min) + ',' + str(max) + ']',
        'counts': {
          "undercall": 0, "overcall": 0,
          "correct": 0, "missing": 0, "total": 0
        }
    }
    dict_median = {
        'children': '[' + str(median) + ',' + str(cnt_children) + ']',
        'counts': {
          "undercall": 0, "overcall": 0,
          "correct": 0, "missing": 0, "total": 0
        }
    }
    for tmp_dict in stat_json['autosomal heterozygous calls error, MODEL']:
        cnt_children_ = tmp_dict['children']
        if cnt_children_ >= min:
            for k, v in tmp_dict.items():
                if k != 'children':
                    for k_v, v_v in v['counts'].items():
                        if cnt_children_ <= max:
                            dict_range['counts'][k_v] += v_v
                        #end if
                        if cnt_children_ >= median:
                            dict_median['counts'][k_v] += v_v
                        #end if
                    #end for
                #end if
            #end for
        #end if
    #end for
    for dict_ in [dict_range, dict_median]:
        if dict_['counts']['total']:
            dict_.setdefault('undercall ratio', round(dict_['counts']['undercall'] / dict_['counts']['total'], 3))
            dict_.setdefault('overcall ratio', round(dict_['counts']['overcall'] / dict_['counts']['total'], 3))
            dict_.setdefault('correct ratio', round(dict_['counts']['correct'] / dict_['counts']['total'], 3))
            dict_.setdefault('missing ratio', round(dict_['counts']['missing'] / dict_['counts']['total'], 3))
        #end if
        stat_json['autosomal heterozygous calls error, MODEL'].append(dict_)
    #end for
#end def

def novo_variants(stat_dict_list, filename):
    ''' '''
    for stat_dict in stat_dict_list: # only one stat_dict at most
        if stat_dict['error_novo_family']:
            fv = open(filename.split('.')[0] + '_total_variants.txt', 'w')
            for bin in sorted(stat_dict['error_novo_family'], reverse=True):
                for cnt_children in sorted(stat_dict['error_novo_family'][bin]):
                    if stat_dict['error_novo_family'][bin][cnt_children]['total_vnt']:
                        with open(filename.split('.')[0] + '_total_' + str(bin) + '-' + str(cnt_children) + '.txt', 'w') as fo:
                            for vnt_obj in stat_dict['error_novo_family'][bin][cnt_children]['total_vnt']:
                                fo.write('{0}\t{1}\t{2}\t{3}\n'.format(vnt_obj.CHROM, vnt_obj.POS, vnt_obj.REF, vnt_obj.ALT))
                                fv.write(vnt_obj.to_string())
                            #end for
                        #end with
                    #end if
                    if stat_dict['error_novo_family'][bin][cnt_children]['no_gparents_vnt']:
                        with open(filename.split('.')[0] + '_no-gprnts_' + str(bin) + '-' + str(cnt_children) + '.txt', 'w') as fo:
                            for vnt_obj in stat_dict['error_novo_family'][bin][cnt_children]['no_gparents_vnt']:
                                fo.write('{0}\t{1}\t{2}\t{3}\n'.format(vnt_obj.CHROM, vnt_obj.POS, vnt_obj.REF, vnt_obj.ALT))
                            #end for
                        #end with
                    #end if
                #end for
            #end for
            fv.close()
        #end if
    #end for
#end def

#################################################################
#    Plots
#################################################################
def plot_error_het_family(stat_dict_list, sample_het_list):
    ''' plot accuracy for autosomal heterozygous variants in parents,
    family '''
    # Labels
    sample_0, sample_1 = sorted(sample_het_list)
    # Data
    data = {
        'bins': [],
        sample_0: [],
        sample_1: []
    }
    # Get data points
    for stat_dict in stat_dict_list:
        if stat_dict['error_het_family']:
            for cnt_children in sorted(stat_dict['error_het_family']):
                if cnt_children not in data['bins']:
                    data['bins'].append(cnt_children)
                #end if
                for sample in stat_dict['error_het_family'][cnt_children]:
                    if sample in sample_het_list:
                        try:
                            correct = stat_dict['error_het_family'][cnt_children][sample]['het'] / \
                                      stat_dict['error_het_family'][cnt_children][sample]['total']
                            data[sample].append(correct * 100)
                        except Exception: data[sample].append(0)
                        #end try
                    #end if
                #end for
            #end for
        #end if
    #end for
    filename = 'autosomal_heterozygous_accuracy_family_{0}-{1}.png'.format(sample_0, sample_1)
    xlabel = 'Number of heterozygous children'
    ylabel = '% of correct calls in parent'
    title = 'Accuracy of autosomal heterozygous variant calls (family)'
    _plot_hist_2(data['bins'], data[sample_0], data[sample_1],
                        sample_0, sample_1, filename, xlabel, ylabel, title, percent=True)
#end def

def plot_distr_het_family(stat_dict_list, sample_het_list):
    ''' plot distribution for autosomal heterozygous variants in parents,
    family '''
    # Labels
    sample_0, sample_1 = sorted(sample_het_list)
    # Data
    data = {
        'bins': [],
        sample_0: [],
        sample_1: []
    }
    # Get data points
    for stat_dict in stat_dict_list:
        if stat_dict['error_het_family']:
            for cnt_children in sorted(stat_dict['error_het_family']):
                if cnt_children not in data['bins']:
                    data['bins'].append(cnt_children)
                #end if
                for sample in stat_dict['error_het_family'][cnt_children]:
                    if sample in sample_het_list:
                        try:
                            correct = stat_dict['error_het_family'][cnt_children][sample]['het']
                            data[sample].append(correct)
                        except Exception: data[sample].append(0)
                        #end try
                    #end if
                #end for
            #end for
        #end if
    #end for
    filename = 'autosomal_heterozygous_distribution_family_{0}-{1}.png'.format(sample_0, sample_1)
    xlabel = 'Number of heterozygous children'
    ylabel = 'Number of variants heterozygous in parent'
    title = 'Distribution of autosomal heterozygous variant calls (family)'
    _plot_hist_2(data['bins'], data[sample_0], data[sample_1],
                        sample_0, sample_1, filename, xlabel, ylabel, title)
#end def

def plot_distr_het(stat_dict_list, sample_het_list):
    ''' plot distribution for autosomal heterozygous variants in parents,
    trio '''
    # Labels
    sample_0, sample_1 = sorted(sample_het_list)
    # Data
    data = {
        'bins': [],
        sample_0: [],
        sample_1: []
    }
    # Get data points
    for stat_dict in stat_dict_list:
        if stat_dict['error_het']:
            for cnt_children in sorted(stat_dict['error_het']):
                if cnt_children not in data['bins']:
                    data['bins'].append(cnt_children)
                #end if
                for sample in stat_dict['error_het'][cnt_children]:
                    if sample in sample_het_list:
                        try:
                            data[sample].append(stat_dict['error_het'][cnt_children][sample]['no_spouse']['is_ID']['in_gparents'])
                        except Exception: data[sample].append(0)
                        #end try
                    #end if
                #end for
            #end for
        #end if
    #end for
    filename = 'autosomal_heterozygous_distribution_trio_{0}-{1}.png'.format(sample_0, sample_1)
    xlabel = 'Number of heterozygous children'
    ylabel = 'Number of variants heterozygous in parent'
    title = 'Distribution of autosomal heterozygous variant calls (trio)'
    _plot_hist_2(data['bins'], data[sample_0], data[sample_1],
                        sample_0, sample_1, filename, xlabel, ylabel, title)
#end def

def plot_AD_DP_ratio(stat_dict, sample_novo):
    ''' '''
    data_0, data_no_0 = [], []
    for cnt_children in stat_dict['error_novo_family'][0.9]:
        if cnt_children == 0:
            data_0 += stat_dict['error_novo_family'][0.9][cnt_children]['no_gparents_AD-DP']
        else:
            data_no_0 += stat_dict['error_novo_family'][0.9][cnt_children]['no_gparents_AD-DP']
        #end if
    #end for
    # Plotting
    pyplot.hist(numpy.array(data_0), alpha=0.5, color='#67a9cf', label='in 0 children')
    pyplot.hist(numpy.array(data_no_0), alpha=0.5, color='#ef8a62', label='in 1+ children')
    # Format and save
    pyplot.gca().set(title='AD-DP ratio distribution (novoPP >= 0.9)', ylabel='Frequency')
    pyplot.legend()
    filename = 'AD-DP_distribution_denovo_{0}.png'.format(sample_novo)
    pyplot.savefig(filename)
#end def

def _plot_hist_2(bins, data_0, data_1, label_0, label_1, filename, xlabel, ylabel, title, percent=False):
    ''' '''
    # Create plot
    fig, ax = pyplot.subplots()
    index = numpy.arange(len(bins))
    bar_width = 0.35
    # Plot data
    rects1 = pyplot.bar(index, data_0, bar_width, color='#ef8a62', label=label_0, zorder=3)
    rects2 = pyplot.bar(index + bar_width + 0.02, data_1, bar_width, color='#67a9cf', label=label_1, zorder=3)
    # Grid and axis
    pyplot.grid(axis='y', alpha=0.8, zorder=0, linestyle='--')
    pyplot.xticks(index + bar_width / 2 + 0.01, bins)
    if percent: pyplot.ylim(0, 120)
    #end if
    yticks = ax.yaxis.get_major_ticks()
    yticks[-1].set_visible(False)
    # Labels and layout
    pyplot.xlabel(xlabel)
    pyplot.ylabel(ylabel)
    pyplot.title(title)
    pyplot.legend()
    pyplot.tight_layout()
    # Save
    pyplot.savefig(filename)
#end def

#################################################################
#    runner
#################################################################
def main(args):
    ''' '''
    # Variables
    NA_chroms = real_NA_chroms
    is_verbose = True if args['verbose'] else False
    sample_novo = args['novo'] if args['novo'] else ''
    sample_het_list = args['het'] if args['het'] else []
    anchor_list = args['anchor']

    # Check pedigree and anchor args
    if len(anchor_list) != len(args['pedigree']) and \
       len(args['pedigree']) != 1:
        sys.exit('\nERROR in parsing arguments: number of pedigrees and anchors is different\n')
    #end if

    # Buffers
    fo = open(args['outputfile'], 'w')

    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf(args['inputfile'])

    # Get pedigree / pedigrees information
    pedigree_list = []
    for pedigree in args['pedigree']:
        # Loading pedigree
        if os.path.isfile(pedigree):
            with open(pedigree) as fi:
                pedigree_list.append(json.load(fi))
            #end with
        else:
            try: pedigree_list.append(json.loads(pedigree))
            except Exception:
                sys.exit('\nERROR in parsing arguments: {0} must be either a json file or a string representing a json\n'
                        .format(pedigree))
            #end try
        #end if
    #end for

    # Check novoPP
    if sample_novo:
        try: novotag, _ = vcf_obj.header.check_tag_definition('novoPP')
        except Exception:
                sys.exit('\nERROR in parsing arguments: novoCaller information missing in VCF file\n')
        #end try
    #end if

    # Creating Pedigree object / objects
    pedigree_obj_list = []
    for pedigree in pedigree_list:
        pedigree_obj_list.append(pedigree_parser.Pedigree(pedigree))
    #end for

    # Initializing stat_dict for each anchor
    stat_dict_list = []
    for anchor in anchor_list:
        stat_dict_list.append({
                            'error_het_family': {},
                            'error_het': {},
                            'error_novo_family': {}
                             })
    #end for

    # Building family / families
    family_list = []
    if len(pedigree_obj_list) == 1:
        for anchor in anchor_list:
            family_list.append(pedigree_obj_list[0].get_family(anchor))
        #end for
    else:
        for i, pedigree_obj in enumerate(pedigree_obj_list):
            family_list.append(pedigree_obj.get_family(anchor_list[i]))
        #end for
    #end if

    # Reading variants
    analyzed = 0
    vnt_obj_ = None
    for i, vnt_obj in enumerate(vcf_obj.parse_variants()):
        if is_verbose:
            sys.stderr.write('\rAnalyzing variant... ' + str(i + 1))
            sys.stderr.flush()
        #end if

        # # Check if chromosome is canonical and in valid format
        # if not check_chrom(vnt_obj.CHROM):
        #     continue
        # #end if
        analyzed += 1

        # Skip MAV that are redundant
        if vnt_obj_:
            if vnt_obj_.CHROM == vnt_obj.CHROM and vnt_obj_.POS == vnt_obj.POS:
                continue
            #end if
        #end if
        vnt_obj_ = vnt_obj

        # Getting and updating stats for each stat_dict / family
        for l, family in enumerate(family_list):
            if sample_novo:
                if anchor_list[l] == sample_novo:
                    get_stats_novo(vnt_obj, stat_dict_list[l], family, NA_chroms, sample_novo, novotag)
                #end if
            #end if
            if sample_het_list:
                if anchor_list[l] in sample_het_list:
                    get_stats_het(vnt_obj, stat_dict_list[l], family, NA_chroms)
                #end if
            #end if
        #end for
    #end for

    # Writing output
    sys.stderr.write('\n\n...Writing results for ' + str(analyzed) + ' analyzed variants out of ' + str(i + 1) + ' total variants\n')
    sys.stderr.flush()

    # Create json
    stat_json = to_json(stat_dict_list, sample_het_list, sample_novo)

    # Write variants
    if sample_novo:
        novo_variants(stat_dict_list, args['outputfile'])
        for l, family in enumerate(family_list):
            if anchor_list[l] == sample_novo:
                try: plot_AD_DP_ratio(stat_dict_list[l], sample_novo)
                except Exception: # no variants with novoPP >= 0.9, skip
                    pass
                #end try
            #end if
        #end for
    #end if

    # Plots
    if len(sample_het_list) == 2:
        plot_error_het_family(stat_dict_list, sample_het_list)
        plot_distr_het_family(stat_dict_list, sample_het_list)
        plot_distr_het(stat_dict_list, sample_het_list)
    #end if

    # Write json to file
    json.dump(stat_json, fo, indent=2, sort_keys=False)

    # Closing buffers
    fo.close()
#end def


#################################################################
#
#    MAIN
#
#################################################################
if __name__ == "__main__":

    main()

#end if
