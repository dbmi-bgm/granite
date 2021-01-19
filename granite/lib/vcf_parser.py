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
import gzip


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
        information as Header object to initialize Vcf object '''
        self.inputfile = inputfile
        self.header = self.parse_header()
    #end def

    class Header(object):
        ''' object to store vcf header information '''

        def __init__(self, definitions, columns, IDs_genotypes):
            ''' initialize Header object '''
            self.definitions = definitions
            self.columns = columns
            self.IDs_genotypes = IDs_genotypes
        #end def

        def add_tag_definition(self, tag_definition, tag_type='INFO'):
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

        def remove_tag_definition(self, tag, tag_type='INFO'):
            ''' remove tag definition from header,
            block specified by tag_type (e.g. FORMAT, INFO) '''
            new_definitions = ''
            for line in self.definitions.split('\n')[:-1]:
                if line.startswith('##' + tag_type + '=<ID=' + tag + ','): ##<tag_type>=<ID=<tag>,...
                    continue
                #end if
                new_definitions += line + '\n'
            #end for
            self.definitions = new_definitions
        #end def

        def get_tag_field_idx(self, tag, field, tag_type='INFO', sep='|'):
            ''' get idx for value field in tag from definition,
            block specified by tag_type (e.g. FORMAT, INFO) '''
            for line in self.definitions.split('\n')[:-1]:
                if line.startswith('##' + tag_type + '=<ID=' + tag + ','):
                    try:
                        format = line.split('Format:')[1]
                        # Cleaning format
                        format = format.replace('\'', '')
                        format = format.replace('\"', '')
                        format = format.replace('>', '')
                    except Exception:
                        raise ValueError('\nERROR in VCF header structure, {0} tag definition has no format specification\n'
                                            .format(tag))
                    #end try
                    # Search exact match
                    # if not exact match, search for partial match (included in field name)
                    for i, field_i in enumerate(format.split(sep)):
                        if field == field_i.strip(): return i # exact match
                        #end if
                    #end for
                    for i, field_i in enumerate(format.split(sep)):
                        if field in field_i.strip(): return i # partial match
                        #end if
                    #end for
                #end if
            #end for
            raise ValueError('\nERROR in VCF header structure, {0} tag definition is missing\n'
                                .format(tag))
        #end def

        def check_tag_definition(self, tag, tag_type='INFO', sep='|'):
            ''' check if tag is standalone or field of another leading tag,
            return leading tag and field index, if any, to acces requested tag '''
            for line in self.definitions.split('\n')[:-1]:
                if line.startswith('##' + tag_type):
                    if ('=<ID=' + tag + ',') in line: ##<tag_type>=<ID=<tag>,..
                        # tag is already a standalone tag
                        return tag, 0
                    elif tag in line and 'Format:' in line: ##<tag_type>=<ID=<lead_tag>,...,Description="... Format:<tag>">
                        # tag is a field, get leading tag and field index
                        lead_tag = line.split('=<ID=')[1].split(',')[0]
                        idx = self.get_tag_field_idx(lead_tag, tag, tag_type, sep)
                        return lead_tag, idx
                    #end if
                #end if
            #end for
            raise ValueError('\nERROR in VCF header structure, {0} tag definition is missing\n'
                                .format(tag))
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

        def repr(self):
            ''' variant representation as CHROM:POSREF>ALT'''
            return '{0}:{1}{2}>{3}'.format(self.CHROM,
                                      self.POS,
                                      self.REF,
                                      self.ALT)
        #end def

        def remove_tag_genotype(self, tag_to_remove, sep=':'):
            ''' remove tag field from FORMAT and GENOTYPES '''
            idx_tag_to_remove, new_format = -1, []
            # Removing tag field from FORMAT
            for i, tag in enumerate(self.FORMAT.split(sep)):
                if tag_to_remove == tag:
                    idx_tag_to_remove = i
                else:
                    new_format.append(tag)
                #end if
            #end for
            # Error if tag_to_remove not found in FORMAT
            if idx_tag_to_remove == -1:
                raise ValueError('\nERROR in variant FORMAT field, {0} tag is missing\n'
                            .format(tag_to_remove))
            #end if
            # Updating FORMAT
            self.FORMAT = sep.join(new_format)
            # Removing tag field from GENOTYPES
            for ID_genotype, genotype in self.GENOTYPES.items():
                genotype_as_list = genotype.split(sep)
                try:
                    del genotype_as_list[idx_tag_to_remove]
                except Exception: # del will fail for trailing fields that are dropped
                                  # field to remove is missing already
                    pass
                #end try
                self.GENOTYPES[ID_genotype] = sep.join(genotype_as_list)
            #end for
        #end def

        def complete_genotype(self, sep=':'):
            ''' fill the trailing fields dropped in GENOTYPES,
            based on FORMAT structure '''
            len_FORMAT = len(self.FORMAT.split(sep))
            for ID_genotype, genotype in self.GENOTYPES.items():
                genotype_as_list = genotype.split(sep)
                for i in range(len_FORMAT - len(genotype_as_list)):
                    genotype_as_list.append('.')
                #end for
                self.GENOTYPES[ID_genotype] = sep.join(genotype_as_list)
            #end for
        #end def

        def empty_genotype(self, sep=':'):
            ''' return a empty genotype based on FORMAT structure '''
            len_FORMAT = len(self.FORMAT.split(sep))
            return './.' + (sep + '.') * (len_FORMAT - 1)
        #end def

        def remove_tag_info(self, tag_to_remove, sep=';'):
            ''' remove tag field from INFO '''
            new_INFO = []
            for tag in self.INFO.split(sep):
                if tag.startswith(tag_to_remove + '='):
                    continue
                #end if
                new_INFO.append(tag)
            #end for
            self.INFO = sep.join(new_INFO)
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
            if self.INFO.endswith(sep): # if INFO ending is wrongly formatted
                self.INFO += tag_to_add
            else:
                self.INFO += sep + tag_to_add
            #end if
        #end def

        def get_tag_value(self, tag_to_get, sep=';'):
            ''' get value from tag (tag_to_get) in INFO '''
            for tag in self.INFO.split(sep):
                if tag.startswith(tag_to_get + '='):
                    try:
                        return tag.split(tag_to_get + '=')[1]
                    except Exception: # tag field is in a wrong format
                        raise ValueError('\nERROR in variant INFO field, {0} tag is in the wrong format\n'
                                    .format(tag_to_get))
                    #end try
                #end if
            #end for

            # tag_to_get not found
            raise ValueError('\nERROR in variant INFO field, {0} tag is missing\n'.format(tag_to_get))
        #end def

        def get_genotype_value(self, ID_genotype, tag_to_get, sep=':'):
            ''' get value from tag (tag_to_get) in genotype specified by corresponding ID '''
            # Get index from FORMAT
            idx_tag_to_get = -1
            for i, tag in enumerate(self.FORMAT.split(sep)):
                if tag_to_get == tag:
                    idx_tag_to_get = i
                    break
                #end if
            #end for
            # Error if tag_to_get not found in FORMAT
            if idx_tag_to_get == -1:
                raise ValueError('\nERROR in variant FORMAT field, {0} tag is missing\n'
                            .format(tag_to_get))
            #end if
            # Get value from index in genotype by ID
            try:
                return self.GENOTYPES[ID_genotype].split(sep)[idx_tag_to_get]
            except Exception:
                raise ValueError('\nERROR in GENOTYPES identifiers, {0} identifier is missing in VCF\n'
                            .format(ID_genotype))
            #end try
        #end def

    #end class Variant

    @staticmethod
    def read_vcf(inputfile):
        ''' read vcf file, gzipped or ungzipped,
        return a generator '''
        if inputfile.endswith('.gz'):
            with gzip.open(inputfile, 'rb') as fz:
                for byteline in fz:
                    yield byteline.decode()
                #end for
            #end with
        else:
            with open(inputfile) as fi:
                for line in fi:
                    yield line
                #end for
            #end with
        #end if
    #end def

    def parse_header(self):
        ''' read header and save information as Header object '''
        definitions, columns, IDs_genotypes = '', '', ''
        for line in self.read_vcf(self.inputfile):
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

        # Checking header is correct
        if definitions and columns and IDs_genotypes:
            return self.Header(definitions, columns, IDs_genotypes)
        else:
            raise ValueError('\nERROR in VCF header structure, missing essential lines\n')
        #end if
    #end def

    def parse_variants(self): # generator
        ''' return a generator to variants stored as Variant objects '''
        for line in self.read_vcf(self.inputfile):
            if not line.startswith('#'):
                line_strip = line.rstrip()
                if line_strip:
                    try:
                        yield self.Variant(line_strip, self.header.IDs_genotypes)
                    except Exception:
                        raise ValueError('\nERROR in variant VCF structure, missing essential columns\n')
                    #end try
                #end if
            #end if
        #end for
    #end def

    def write_definitions(self, outputfile_obj):
        ''' write header definitions to outputfile_obj buffer '''
        outputfile_obj.write(self.header.definitions)
    #end def

    def write_columns(self, outputfile_obj):
        ''' write header columns to outputfile_obj buffer,
        #CHROM ... '''
        outputfile_obj.write(self.header.columns)
    #end def

    def write_header(self, outputfile_obj):
        ''' write header definitions and columns to outputfile_obj buffer '''
        self.write_definitions(outputfile_obj)
        self.write_columns(outputfile_obj)
    #end def

    def write_variant(self, outputfile_obj, parsed_variant):
        ''' write parsed_variant (Variant object generated by parse_variants) to outputfile_obj buffer '''
        outputfile_obj.write(parsed_variant.to_string())
    #end def

#end class Vcf
