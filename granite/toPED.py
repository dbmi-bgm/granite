#!/usr/bin/env python

#################################################################
#
#    toPED
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
# pedigree_parser
from granite.lib import pedigree_parser


#################################################################
#
#    FUNCTIONS
#
#################################################################
#################################################################
#    runner
#################################################################
def main(args):
    ''' '''
    # Variables
    encode_gender = {'M': 1, 'F': 2}
    familyID = args['family'] if args['family'] else 'FAM'

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

    # Buffers
    fo = open(args['outputfile'], 'w')

    # Writing PED
    for _, member_obj in pedigree_obj.members.items():
        if member_obj.is_sample():
            memberID = member_obj.sample
            try: sex = encode_gender[member_obj.gender]
            except Exception: sex = 'U'
            #end try
            paternalID, maternalID = 0, 0
            if member_obj.has_parents():
                for parent in member_obj.get_parents():
                    if parent.gender == 'M':
                        paternalID = parent.sample
                    elif parent.gender == 'F':
                        maternalID = parent.sample
                    #end if
                #end for
            #end if
            fo.write('{0}\t{1}\t{2}\t{3}\t{4}\t-9\n'.format(
                                                    familyID,
                                                    memberID,
                                                    paternalID,
                                                    maternalID,
                                                    sex
                                                    ))
        #end if
    #end for

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
