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
# shared_functions as *
from granite.lib.shared_functions import *
# vcf_parser
from granite.lib import vcf_parser
# pedigree_parser
from granite.lib import pedigree_parser
# shared_vars
from granite.lib.shared_vars import real_NA_chroms
# matplotlib
from matplotlib import pyplot


#################################################################
#
#    FUNCTIONS
#
#################################################################
#################################################################
#    Stats
#################################################################
def get_stats(vnt_obj, stat_dict, family, NA_chroms):
    ''' extract information from variant using family '''
    var_type = variant_type_ext(vnt_obj.REF, vnt_obj.ALT)
    if vnt_obj.CHROM.replace('chr', '') in NA_chroms:
        return # skip unbalanced chromosomes
    #end if
    if var_type == 'snv':
        _error_het(vnt_obj, stat_dict, family)
    #end if
#end def

def _error_het(vnt_obj, stat_dict, family):
    ''' calculate error model for heterozygous calls using family '''
    # Counting children that are heterozygous
    cnt_children = 0
    for child in family['children']:
        try:
            GT_ = vnt_obj.get_genotype_value(child.sample, 'GT').replace('|', '/')
        except Exception: continue
        #end try
        if GT_ == './.': return
        elif GT_ in ['0/1', '1/0']: cnt_children += 1
        #end if
    #end for
    if not cnt_children: return
    #end if
    # Parents samples
    sample_0 = family['parents'][0].sample
    sample_1 = family['parents'][1].sample
    # Update stat_dict for cnt_children
    stat_dict['error_het'].setdefault(cnt_children, {})
    for sample in [sample_0, sample_1]:
        stat_dict['error_het'][cnt_children].setdefault(sample, {
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
            stat_dict['error_het'][cnt_children][sample_1][trio_1[0]] += 1
            stat_dict['error_het'][cnt_children][sample_1]['total'] += 1
        elif trio_1 == ['ref', {'ref'}] and trio_0[1] == {'ref', 'het'}:
            stat_dict['error_het'][cnt_children][sample_0][trio_0[0]] += 1
            stat_dict['error_het'][cnt_children][sample_0]['total'] += 1
        #end if
    #end if
#end def

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
#    Build family structure from pedigree
#################################################################
def get_family(pedigree_obj, anchor):
    ''' given anchor in pedigree,
    build family structure around it '''
    family = {
        'children': [],
        'parents': []
    }
    anchor_obj = pedigree_obj.get_member_by_sample(anchor)
    if anchor_obj.get_children(): # anchor is the center of the family,
                                  # check for parents and children
        if _check_parents(anchor_obj): # parents information complete
            spouse = None
            for spouse_ in anchor_obj.get_spouses():
                if _check_parents(spouse_):
                    spouse = spouse_
                    break
                #end if
            #end for
            if not spouse:
                sys.exit('\nERROR in building family from pedigree: missing parents information for sample {0} spouse\n'
                        .format(anchor))
            #end if
            # Create family
            family['children'] = anchor_obj.common_children(spouse)
            family['parents'] = [obj for obj in sorted([anchor_obj, spouse], key=lambda x: x.name)]
        else: # missing parents
            sys.exit('\nERROR in building family from pedigree: missing parents information for sample {0}\n'
                    .format(anchor))
        #end if
    else: # anchor is the newest generation,
          # check for siblings, parents and grandparents
        if _check_parents(anchor_obj): # parents information complete, check grandparents
            parents = anchor_obj.get_parents()
            for parent in parents:
                if not _check_parents(parent):
                    sys.exit('\nERROR in building family from pedigree: missing family information for sample {0}\n'
                            .format(anchor))
                #end if
            #end for
            # Create family
            family['children'] = parents[0].common_children(parents[1])
            family['parents'] = [obj for obj in sorted(parents, key=lambda x: x.name)]
        else:
            sys.exit('\nERROR in building family from pedigree: missing parents information for sample {0}\n'
                    .format(anchor))
        #end if
    #end if
    return family
#end def

def _check_parents(member_obj):
    ''' check if member_object has two parents with samples '''
    parents = member_obj.get_parents()
    if len(parents) < 2: return False # missing parent
    #end if
    for parent in parents:
        if not parent.is_sample(): # missing parent sample information
            return False
        #end if
    #end for
    return True
#end def

#################################################################
#
#################################################################
def to_json(stat_dict):
    ''' '''
    stat_json = {
        'autosomal heterozygous calls error': []
    }
    # Error heterozygous call
    for cnt_children in sorted(stat_dict['error_het']):
        tmp_dict = {
            'children': cnt_children
        }
        for sample in stat_dict['error_het'][cnt_children]:
            tmp_dict.setdefault(sample, {})
            _counts = {
                'total': stat_dict['error_het'][cnt_children][sample]['total']
            }
            for k, v in stat_dict['error_het'][cnt_children][sample].items():
                # Counts
                if k == 'ref': _counts.setdefault('undercall', v)
                elif k == 'hom': _counts.setdefault('overcall', v)
                elif k == 'het': _counts.setdefault('correct', v)
                elif k == 'missing': _counts.setdefault('missing', v)
                #end if
            #end for
            # Ratios
            if _counts['total']:
                tmp_dict[sample].setdefault('undercall ratio', round(_counts['undercall'] / _counts['total'], 3))
                tmp_dict[sample].setdefault('overcall ratio', round(_counts['overcall'] / _counts['total'], 3))
                tmp_dict[sample].setdefault('correct ratio', round(_counts['correct'] / _counts['total'], 3))
                tmp_dict[sample].setdefault('missing ratio', round(_counts['missing'] / _counts['total'], 3))
            #end if
            tmp_dict[sample].setdefault('counts', _counts)
        #end for
        stat_json['autosomal heterozygous calls error'].append(tmp_dict)
    #end for
    _extend_json(stat_json)
    return stat_json
#end def

def _extend_json(stat_json):
    ''' '''
    # Error heterozygous call in cnt_children range
    cnt_children = stat_json['autosomal heterozygous calls error'][-1]['children']
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
    for tmp_dict in stat_json['autosomal heterozygous calls error']:
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
        stat_json['autosomal heterozygous calls error'].append(dict_)
    #end for
#end def

#################################################################
#   Plots
#################################################################
def plot_error_het(stat_dict, family):
    ''' '''
    # Labels
    sample_0 = family['parents'][0].sample
    sample_1 = family['parents'][1].sample
    # Data
    data = {
        'bins': [],
        sample_0: [],
        sample_1: []
    }
    # Get data points
    for cnt_children in sorted(stat_dict['error_het']):
        data['bins'].append(cnt_children)
        for sample in stat_dict['error_het'][cnt_children]:
            try:
                correct = stat_dict['error_het'][cnt_children][sample]['het'] / stat_dict['error_het'][cnt_children][sample]['total']
                data[sample].append(correct * 100)
            except Exception: data[sample].append(0)
            #end try
        #end for
    #end for
    filename = 'autosomal_heterozygous_accuracy_{0}-{1}.png'.format(sample_0, sample_1)
    xlabel = 'Number of heterozygous children'
    ylabel = '% of correct calls in parent'
    title = 'Accuracy of autosomal heterozygous variant calls'
    _plot_hist_2_percent(data['bins'], data[sample_0], data[sample_1],
                        sample_0, sample_1, filename, xlabel, ylabel, title)
#end def

def _plot_hist_2_percent(bins, data_0, data_1, label_0, label_1, filename, xlabel, ylabel, title):
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
    pyplot.ylim(0, 120)
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
    anchor = args['sample']

    # Buffers
    fo = open(args['outputfile'], 'w')

    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf(args['inputfile'])

    # Loading pedigree
    if os.path.isfile(args['pedigree']):
        with open(args['pedigree']) as fi:
            pedigree = json.load(fi)
        #end with
    else:
        try: pedigree = json.loads(args['pedigree'])
        except Exception:
            sys.exit('\nERROR in parsing arguments: pedigree must be either a json file or a string representing a json\n')
        #end try
    #end if

    # Creating Pedigree object
    pedigree_obj = pedigree_parser.Pedigree(pedigree)

    # Initializing stat_dict
    stat_dict = {
                'error_het': {}
    }

    # Building family
    family = get_family(pedigree_obj, anchor)

    # Reading variants
    analyzed = 0
    for i, vnt_obj in enumerate(vcf_obj.parse_variants(args['inputfile'])):
        if is_verbose:
            sys.stderr.write('\rAnalyzing variant... ' + str(i + 1))
            sys.stderr.flush()
        #end if

        # # Check if chromosome is canonical and in valid format
        # if not check_chrom(vnt_obj.CHROM):
        #     continue
        # #end if
        analyzed += 1

        # Getting and updating stats
        get_stats(vnt_obj, stat_dict, family, NA_chroms)
    #end for

    # Writing output
    sys.stderr.write('\n\n...Writing results for ' + str(analyzed) + ' analyzed variants out of ' + str(i + 1) + ' total variants\n')
    sys.stderr.flush()

    # Plots
    plot_error_het(stat_dict, family)

    # Create json
    stat_json = to_json(stat_dict)

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
