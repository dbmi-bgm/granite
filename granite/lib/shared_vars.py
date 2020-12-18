#!/usr/bin/env python

#################################################################
#
#    shared_vars
#        Michele Berselli
#        Harvard Medical School
#        berselli.michele@gmail.com
#
#################################################################

#################################################################
#
#   SpliceAI delta scores fields
#
#################################################################
DStags = {'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL'}


#################################################################
#
#   VEP terms corresponding to intronic, intergenic or
#       regulatory regions to be removed or ignored
#
#################################################################
VEPremove = {
            # intronic and intergenic features tags
            'intron_variant', 'intergenic_variant',
            'downstream_gene_variant', 'upstream_gene_variant',
            'NMD_transcript_variant', 'non_coding_transcript_variant',
            'non_coding_transcript_exon_variant',
            # regulatory features tags
            'feature_elongation', 'feature_truncation',
            'regulatory_region_variant', 'regulatory_region_amplification',
            'regulatory_region_ablation', 'splice_region_variant',
            'TFBS_amplification', 'TFBS_ablation', 'TF_binding_site_variant'
            }


#################################################################
#
#   VEP terms to rescue only for SpliceAI variants
#
#################################################################
VEPSpliceAI = {'intron_variant'}


#################################################################
#
#   Unbalanced chromosomes to ignore in novoCaller, PP -> NA
#
#################################################################
real_NA_chroms = {'M', 'MT', 'X', 'Y'}

# for automated tests
test_NA_chroms = {'2'}


#################################################################
#
#   VEP Consequence to impact
#
#################################################################
VEP_encode = {
                # HIGH
                'transcript_ablation': 1,
                'splice_acceptor_variant': 1,
                'splice_donor_variant': 1,
                'stop_gained': 1,
                'frameshift_variant': 1,
                'stop_lost': 1,
                'start_lost': 1,
                'transcript_amplification': 1,
                # MODERATE
                'inframe_insertion': 2,
                'inframe_deletion': 2,
                'missense_variant': 2,
                'protein_altering_variant': 2,
                # LOW
                'splice_region_variant': 3,
                'incomplete_terminal_codon_variant': 3,
                'start_retained_variant': 3,
                'stop_retained_variant': 3,
                'synonymous_variant': 3,
                # MODIFIER
                'MODIFIER': 4
            }
