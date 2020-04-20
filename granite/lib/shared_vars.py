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
