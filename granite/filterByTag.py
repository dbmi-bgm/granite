#!/usr/bin/env python

#################################################################################
#
#    filterByTag
#        Michele Berselli
#        Harvard Medical School
#        berselli.michele@gmail.com
#
#################################################################################


#################################################################
#
#    LIBRARIES
#
#################################################################
import sys
# vcf_parser
from granite.lib import vcf_parser


#################################################################################
#
#    OBJECTS
#
#################################################################################
class Tag:
    def __init__(self, name, idx, value, operator, type_, logic, field_sep, entry_sep, value_sep):
        self.name = name
        self.idx = idx
        self.value = value
        self.operator = operator
        self.type_ = type_
        self.logic = logic
        self.field_sep = field_sep
        self.entry_sep = entry_sep
        self.value_sep = value_sep


#################################################################
#
#    FUNCTIONS
#
#################################################################
def get_tag_value(vnt_obj, tag_obj, info_sep):
    ''' get the value(s) of a specific tag from the variant object '''
    if tag_obj.type_ == 'bool':
        # Return as a 1-element list so filter_check can iterate uniformly
        return [vnt_obj.get_tag_value(tag_obj.name, is_flag=True, sep=info_sep)]

    try:
        raw = vnt_obj.get_tag_value(tag_obj.name, sep=info_sep)
    except vcf_parser.MissingTag:
        raw = ''

    # Split into entries if requested (e.g., per-transcript lists like CSQ)
    entries = raw.split(tag_obj.entry_sep) if tag_obj.entry_sep else [raw]

    values = []
    for entry in entries:
        if tag_obj.field_sep:
            parts = entry.split(tag_obj.field_sep)
            value = parts[tag_obj.idx] if tag_obj.idx < len(parts) else ''
        else:
            value = entry

        # Split into multiple values if needed
        split_vals = value.split(tag_obj.value_sep) if tag_obj.value_sep else [value]

        # Normalize numeric missing *per element*
        for value_ in split_vals:
            if (value_ == '' or value_ == '.') and tag_obj.type_ in ['int', 'float']:
                value_ = '0'
            values.append(value_)

    # Apply type conversion (per-element so we can pinpoint failures)
    if tag_obj.type_ in ('int', 'float'):
        cast = int if tag_obj.type_ == 'int' else float
        casted = []
        for idx, v in enumerate(values):
            try:
                casted.append(cast(v))
            except ValueError:
                # Try to include some helpful context
                vnt_repr = vnt_obj.repr()
                # Mention separators and field index used to extract this value
                sys.exit(
                    f'\nERROR in tag value: cannot convert token "{v}" to {tag_obj.type_} for tag "{tag_obj.name}"\n'
                    f'  Value extracted using field index: {tag_obj.idx}\n'
                    f'  Token index within flattened values: {idx}\n'
                    f'  Separators: entry_sep="{tag_obj.entry_sep}", field_sep="{tag_obj.field_sep}", value_sep="{tag_obj.value_sep}"\n'
                    f'  Variant: {vnt_repr}\n'
                    'Please confirm the separators and field index are correct\n'
                )
        values = casted
    # else: str -> leave as-is; bool handled above

    return values

def _eval_one(value, op, target):
    ''' evaluate a single value against the target using the specified operator '''
    if op == '==':   return value == target
    if op == '!=':   return value != target
    if op == '<':    return value < target
    if op == '>':    return value > target
    if op == '<=':   return value <= target
    if op == '>=':   return value >= target
    if op == '~':    return str(target) in str(value)
    if op == '!~':   return str(target) not in str(value)
    if op == 'true': return value is True
    if op == 'false': return value is False
    sys.exit(f'\nERROR in tag evaluation: unknown operator {op}\n')

def filter_check(vnt_obj, tag_obj, info_sep):
    ''' check if entry satisfies the filter logic (any or all) '''
    values = get_tag_value(vnt_obj, tag_obj, info_sep)
    op = tag_obj.operator
    target = tag_obj.value

    results = [_eval_one(v, op, target) for v in values]
    return any(results) if tag_obj.logic == 'any' else all(results)

def main(args):
    ''' '''
    # Variables
    info_sep = args['separator']
    tags_logic = args['logic']

    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf(args['inputfile'])

    # Parsing tag filters
    tag_filters = []
    for tag in args['tag']:
        # Splitting the tag into its components
        tag_ = tag.split('/')
        if len(tag_) < 5:
            sys.exit(f'\nERROR in tag filter format: {tag}. Expected format: name/value/operator/type/logic[/entry=sep][/field=sep][/value=sep]\n')
        name_, value, operator, type_, logic = tag_[:5]
        # Sanitize components
        operator = operator.lower()
        type_    = type_.lower()
        logic    = logic.lower()
        # Get separators
        field_sep = None
        entry_sep = None
        value_sep = None
        for item in tag_[5:]:
            try:
                key, sep = item.split('=', 1)
            except ValueError:
                sys.exit(f'\nERROR in tag filter: malformed option "{item}" in {tag}. Expected key=value\n')
            if key == 'field':
                if sep == '':
                    sys.exit(f'\nERROR in tag filter: empty field separator in {tag}\n')
                field_sep = sep
            elif key == 'entry':
                if sep == '':
                    sys.exit(f'\nERROR in tag filter: empty entry separator in {tag}\n')
                entry_sep = sep
            elif key == 'value':
                if sep == '':
                    sys.exit(f'\nERROR in tag filter: empty value separator in {tag}\n')
                value_sep = sep
            else:
                sys.exit(f'\nERROR in tag filter: unknown extra option "{key}" in {tag}\n')
        # Validate components
        if type_ not in ['str', 'int', 'float', 'bool']:
            sys.exit(f'\nERROR in tag filter type: {type_}. Accepted: str, int, float, bool\n')
        if type_ == 'bool':
            if operator not in ['true', 'false']:
                sys.exit(f'\nERROR in tag filter operator for bool: {operator}. Use "true" or "false"\n')
        if logic not in ['any', 'all']:
            sys.exit(f'\nERROR in tag filter logic: {logic}. Accepted: any, all\n')
        if operator not in ['==', '!=', '<', '>', '<=', '>=', '~', '!~', 'true', 'false']:
            sys.exit(f'\nERROR in tag filter operator: {operator}. Accepted: ==, !=, <, >, <=, >=, ~, !~, true, false\n')
        if (
            (operator in ['true', 'false'] and type_ != 'bool') or
            (operator in ['~', '!~'] and type_ != 'str') or
            (operator in ['<', '<=', '>', '>='] and type_ not in ['int', 'float'])
        ):
            sys.exit(f'\nERROR in tag filter type for operator {operator}: {type_}\n')
        if value == '':
            sys.exit('\nERROR in tag filter value: empty string not allowed for value\n')
        # Need a field_sep for check_tag_definition, the function will take care of everything else
        # If field_sep is not provided, it will default to '|' (for VEP-like annotations)
        field_sep_ = field_sep if field_sep is not None else '|'
        name, idx = vcf_obj.header.check_tag_definition(name_, sep=field_sep_)
        # If check_tag_definition changed the name (e.g., CSQ -> ANN), the tag is actually a subfield
        # Use default field_sep_
        if name != name_ and field_sep is None:
            field_sep = field_sep_
        # Converting value
        if type_ == 'int':
            try:
                value = int(value)
            except ValueError: sys.exit(f'\nERROR in tag filter value: cannot convert "{value}" to int\n')
        elif type_ == 'float':
            try:
                value = float(value)
            except ValueError: sys.exit(f'\nERROR in tag filter value: cannot convert "{value}" to float\n')
        # Create Tag object
        tag_filters.append(Tag(name, idx, value, operator, type_, logic, field_sep, entry_sep, value_sep))

    # Parsing variants and write output
    with open(args['outputfile'], 'w', encoding='utf-8') as fo:
        vcf_obj.write_header(fo)

        kept = 0
        analyzed = 0
        for i, vnt_obj in enumerate(vcf_obj.parse_variants()):
            analyzed += 1
            if args['verbose']:
                sys.stderr.write('\rAnalyzing variant... ' + str(i + 1))
                sys.stderr.flush()

            # Check filters
            if tags_logic == 'any':
                for t in tag_filters:
                    if filter_check(vnt_obj, t, info_sep):
                        vcf_obj.write_variant(fo, vnt_obj)
                        kept += 1
                        break
            else:
                check_results = [filter_check(vnt_obj, t, info_sep) for t in tag_filters]
                if all(check_results):
                    vcf_obj.write_variant(fo, vnt_obj)
                    kept += 1

        sys.stderr.write(f'\n\nWrote {kept} variants out of {analyzed} analyzed\n')
        sys.stderr.flush()


#################################################################
#
#    MAIN
#
#################################################################
if __name__ == '__main__':

    main()
