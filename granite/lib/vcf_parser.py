#!/usr/bin/env python

#################################################################
#
#    vcf_parser
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
#    Vcf
#      -> Header
#      -> Variant
#
#################################################################
class Vcf(object):
    ''' object to read and manipulate vcf file format '''

    def __init__(self, inputfile):
        ''' open input vcf, read header lines and save
        information as Header object to initialize Vcf object'''
        self.header = self.parse_header(inputfile)
    #end def

    class Header(object):
        ''' object to store vcf header information '''

        def __init__(self, definitions, columns, IDs_genotypes):
            ''' initialize Header object '''
            self.definitions = definitions
            self.columns = columns
            self.IDs_genotypes = IDs_genotypes
        #end def

        def add_tag_definition(self, tag_definition, tag_type='FORMAT'):
            ''' add tag_definition to the header on top
            of the block specified by tag_type (e.g. FORMAT, INFO) '''
            added_tag, new_definitions = False, ''
            for line in self.definitions.split('\n')[:-1]:
                if line.startswith('##' + tag_type) and not added_tag:
                    added_tag = True
                    new_definitions += tag_definition + '\n'
                #end if
                new_definitions += line + '\n'
            #end for
            self.definitions = new_definitions
        #end def

    #end class Header

    class Variant(object):
        ''' object to store information for variant in vcf format '''

        def __init__(self, line_strip, IDs_genotypes):
            ''' initialize Variant object '''
            line_split = line_strip.split('\t')
            self.CHROM = line_split[0]
            self.POS = int(line_split[1])
            self.ID = line_split[2]
            self.REF = line_split[3]
            self.ALT = line_split[4]
            self.QUAL = line_split[5]
            self.FILTER = line_split[6]
            self.INFO = line_split[7]
            self.FORMAT = line_split[8]
            self.IDs_genotypes = IDs_genotypes
            self.GENOTYPES = {k: v for k, v in zip(IDs_genotypes, line_split[9:])}
        #end def

        def to_string(self):
            ''' variant as string rapresentation '''
            genotypes_as_list = []
            variant_as_string = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t'.format(self.CHROM,
                                                                                        self.POS,
                                                                                        self.ID,
                                                                                        self.REF,
                                                                                        self.ALT,
                                                                                        self.QUAL,
                                                                                        self.FILTER,
                                                                                        self.INFO,
                                                                                        self.FORMAT)
            for IDs_genotype in self.IDs_genotypes:
                genotypes_as_list.append(self.GENOTYPES[IDs_genotype])
            #end for

            return variant_as_string + '\t'.join(genotypes_as_list) + '\n'
        #end def

        def remove_tag_genotype(self, tag_to_remove, sep=':'):
            ''' remove tag field from FORMAT and GENOTYPES '''
            idx_tag_to_remove, new_format = 0, []
            # Removing tag field from FORMAT
            for i, tag in enumerate(self.FORMAT.split(sep)):
                if tag_to_remove == tag:
                    idx_tag_to_remove = i
                else:
                    new_format.append(tag)
                #end if
            #end for
            self.FORMAT = sep.join(new_format)
            # Removing tag field from GENOTYPES
            for ID_genotype, genotype in self.GENOTYPES.items():
                genotype_as_list = genotype.split(sep)
                del genotype_as_list[idx_tag_to_remove]
                self.GENOTYPES[ID_genotype] = sep.join(genotype_as_list)
            #end for
        #end def

        def remove_tag_info(self, tag_to_remove, sep=';'):
            ''' remove tag field from INFO '''
            self.INFO = re.sub(r'{0}=.*?{1}'.format(tag_to_remove, sep), '', self.INFO)
        #end def

        def add_tag_format(self, tag_to_add, sep=':'):
            ''' add tag field to FORMAT '''
            self.FORMAT += sep + tag_to_add
        #end def

        def add_values_genotype(self, ID_genotype, values, sep=':'):
            ''' add values field to genotype specified by corresponding ID '''
            self.GENOTYPES[ID_genotype] += sep + values
        #end def

        def add_tag_info(self, tag_to_add, sep=';'):
            ''' add tag field and value (tag_to_add) to INFO '''
            # tag_to_add -> tag=<value>
            self.INFO += tag_to_add + sep
        #end def

        def get_tag_value(self, tag_to_get, sep=';'):
            ''' get value from tag (tag_to_get) in INFO '''
            for tag in self.INFO.split(sep):
                if tag.startswith(tag_to_get):
                    try:
                        return tag.split('=')[1]
                    except Exception: # tag field is in a wrong format
                        raise ValueError('ERROR in variant INFO field, {0} tag is in the wrong format\n'
                                    .format(tag_to_get))
                    #end try
                #end if
            #end for

            # tag_to_get not found
            raise ValueError('ERROR in variant INFO field, {0} tag is missing\n'.format(tag_to_get))
        #end def

    #end class Variant

    def parse_header(self, inputfile):
        ''' read header and save information as Header object '''
        definitions, columns, IDs_genotypes = '', '', ''
        with open(inputfile) as fi:
            for line in fi:
                if line.startswith('#'): # reading a header line
                    line_strip = line.rstrip()
                    if line_strip.startswith('##'): # header definition line
                        definitions += line_strip + '\n'
                    elif line_strip.startswith('#CHROM'): # header columns line
                        columns += line_strip + '\n'
                        IDs_genotypes = line_strip.split('\t')[9:]
                    #end if
                else: # finished to read the header
                    break # exit and close buffer
                #end if
            #end for
        #end with

        # Checking header is correct
        if definitions and columns and IDs_genotypes:
            return self.Header(definitions, columns, IDs_genotypes)
        else:
            raise ValueError('ERROR in VCF header structure, missing essential lines\n')
        #end if
    #end def

    def parse_variants(self, inputfile): # generator
        ''' return a generator to variants stored as Variant objects '''
        with open(inputfile) as fi:
            for line in fi:
                if not line.startswith('#'):
                    line_strip = line.rstrip()
                    if line_strip:
                        try:
                            yield self.Variant(line_strip, self.header.IDs_genotypes)
                        except Exception:
                            raise ValueError('ERROR in variant VCF structure, missing essential columns\n')
                        #end try
                    #end if
                #end if
            #end for
        #end with
    #end def

#end class Vcf
