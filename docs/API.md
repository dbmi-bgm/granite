## VCF Parser
granite library can be used directly to access and manipulate information in VCF format.

### Import the library

    from granite.lib import vcf_parser

### Usage
The library implements the objects [*Vcf*](#vcf), [*Header*](#header) and [*Variant*](#variant).

#### Vcf
This is the main object and has methods to read and write VCF format.

##### Initialize the object

    vcf_obj = vcf_parser.Vcf('inputfile.vcf')

This will automatically read the file header into a *Header* object.

##### Read and access variants
The method *parse_variants()* will read the file and return a generator to *Variant* objects that store variants information.

    for vnt_obj in vcf_obj.parse_variants():
        # do something with vnt_obj

##### Write to file
The method *write_header(fo)* allows to write header definitions and columns to specified buffer (fo).

    with open('outputfile.vcf', 'w') as fo:
        vcf_obj.write_header(fo)

It is possible to write only definitions or columns respectively with the methods *write_definitions(fo)* and *write_columns(fo)*.

    with open('outputfile.vcf', 'w') as fo:
        vcf_obj.write_definitions(fo)
        vcf_obj.write_columns(fo)

The method *write_variant(fo, Variant_obj)* allows to write information from *Variant* object to specified buffer (fo).

    with open('outputfile.vcf', 'w') as fo:
        vcf_obj.write_variant(fo, vnt_obj)

#### Header
This is the object used to store information for the header in VCF format. Methods are available to extract and modify information in the header.

##### Attributes
###### definitions *\<str\>*
Stores the full header information minus the last line where columns are defined.

    vcf_obj.header.definitions

###### columns *\<str\>*
Stores the last header line where columns are defined.

    # columns example
    # #CHROM    POS    ID    REF    ALT    ...
    vcf_obj.header.columns

###### IDs_genotypes *\<list\>*
Stores sample ID(s) available in the VCF as list. If multiple samples, the order from the VCF is maintained.

    vcf_obj.header.IDs_genotypes

##### Add or remove definitions
The method *add_tag_definition(tag_definition, tag_type='INFO')* allows to add tag_definition to the header on top of the block specified by tag_type (e.g. FORMAT, INFO).

    tag_definition = '##INFO=<ID=tag,Number=.,Type=.,Description="INFO tag definition example">'
    vcf_obj.header.add_tag_definition(tag_definition)

The method *remove_tag_definition(tag, tag_type='INFO')* allows to remove tag definition from the header block specified by tag_type (e.g. FORMAT, INFO).

    # remove CSQ definition (VEP) from the header INFO block
    tag = 'CSQ'
    vcf_obj.header.remove_tag_definition(tag)

##### Extract information
The method *get_tag_field_idx(tag, field, tag_type='INFO', sep='|')* allows to get the index corresponding to value field in tag from definition, block specified by tag_type (e.g. FORMAT, INFO). sep is the fields separator used in the tag definition.

    # return the index corresponding to 'Consequence' field
    # from CSQ definition (VEP) in the header INFO block
    # ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Gene|Feature|Consequence|IMPACT|...">
    tag, field = 'CSQ', 'Consequence'
    idx <int> = vcf_obj.header.get_tag_field_idx(tag, field) # idx -> 3

The method *check_tag_definition(tag, tag_type='INFO', sep='|')* allows to check if a tag is in the header and if is standalone or field of another leading tag. Returns the leading tag and the field corresponding index, if any, to acces the tag. sep is the fields separator used in the tag definition.

    # return the leading tag and index corresponding to 'Consequence' field
    # from CSQ definition (VEP) in the header INFO block
    # ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Gene|Feature|Consequence|IMPACT|...">
    tag = 'Consequence'
    lead_tag <str>, idx <int> = vcf_obj.header.check_tag_definition(tag) # lead_tag -> CSQ
                                                                         # idx -> 3

*note: tag and field are case sensitive.*

#### Variant
This is the object used to store information for variants in VCF format.

##### Attributes
###### CHROM *\<str\>*
Stores chromosome name (e.g. 1, chr1), as in the VCF file.

    vnt_obj.CHROM

###### POS *\<int\>*
Stores variant position.

    vnt_obj.POS

###### ID *\<str\>*
Stores variant ID(s), as in the VCF file.

    vnt_obj.ID

###### REF *\<str\>*
Stores reference allele at position.

    vnt_obj.REF

###### ALT *\<str\>*
Stores alternate allele(s) at position.

    vnt_obj.ALT

###### QUAL *\<str\>*
Stores phred-scaled quality score for the assertion made in ALT.

    vnt_obj.QUAL

###### FILTER *\<str\>*
Stores filter status.

    vnt_obj.FILTER

###### INFO *\<str\>*
Additional information for the variant.

    vnt_obj.INFO

###### FORMAT *\<str\>*
Stores specification for the genotype column(s) structure.

    vnt_obj.FORMAT

###### IDs_genotypes *\<list\>*
Stores sample ID(s) available in the VCF as list. If multiple samples, the order from the VCF is maintained.

    vnt_obj.IDs_genotypes

###### GENOTYPES *\<dict\>*
Stores a dictionary linking genotype(s) for the variant to corresponding sample ID(s).

    # {ID_genotype: genotype, ...}
    vnt_obj.GENOTYPES

##### Format variants
The method *to_string()* returns the variant representation in VCF format.

    vnt_vcf <str> = vnt_obj.to_string()

The method *repr()* returns the variant representation as *CHROM:POSREF>ALT*.

    vnt_repr <str> = vnt_obj.repr()

##### Manipulate genotype(s)
The method *remove_tag_genotype(tag, sep=':')* allows to remove a tag from FORMAT and GENOTYPES. sep is the tags separator used in format definition and genotype(s).

    # remove AD tag from format definition and genotype(s)
    tag = 'AD'
    vnt_obj.remove_tag_genotype(tag)

The method *complete_genotype(sep=':')* fills in the trailing fields that are missing and by default dropped in GENOTYPES. sep is the tags separator used in format definition and genotype(s).

    vnt_obj.complete_genotype()

The method *empty_genotype(sep=':')* returns a empty genotype based on FORMAT structure. sep is the tags separator used in format definition and genotype(s).

    empty <str> = vnt_obj.empty_genotype()

The method *add_tag_format(tag, sep=':')* allows to add a tag at the end of FORMAT structure. sep is the tags separator used in format definition and genotype(s).

    # add RSTR tag to format
    tag = 'RSTR'
    vnt_obj.add_tag_format(tag)

The method *add_values_genotype(ID_genotype, values, sep=':')* allows to add values at the end of the genotype specified by corresponding ID. sep is the tags separator used in format definition and genotype(s).

    vnt_obj.add_values_genotype(ID_genotype, values)

The method *get_genotype_value(ID_genotype, tag, complete_genotype=False, sep=':')* returns value for tag from the genotype specified by corresponding ID. sep is the tags separator used in format definition and genotype(s).
If complete_genotype=True, return '.' if tag is missing. If complete_genotype=False (default) raise exception for the missing tag.

    tag_val <str> = vnt_obj.get_genotype_value(ID_genotype, tag)

##### Manipulate INFO
The method *remove_tag_info(tag, sep=';')* allows to remove a tag from INFO. sep is the tags separator used in INFO.

    vnt_obj.remove_tag_info(tag)

The method *add_tag_info(tag_value, sep=';')* allows to add a tag and its value at the end of INFO. sep is the tags separator used in INFO.

    # add tag and value to INFO
    tag_value = 'tag=value'
    vnt_obj.add_tag_info(tag_value)

The method *get_tag_value(tag, sep=';')* returns the value from tag in INFO. sep is the tags separator used in INFO.

    tag_val <str> = vnt_obj.get_tag_value(tag)

*note: tag and ID are case sensitive.*

### Custom error classes

*MissingTag* describes a missing tag or tag value.

*MissingTagDefinition* describes a missing tag definition.

*TagDefinitionError* describes a format error for a tag definition.

*TagFormatError* describes a format error for a tag.

*MissingIdentifier* describes a missing genotype identifier in the VCF file.

*VcfFormatError* describes an error in the VCF format.
