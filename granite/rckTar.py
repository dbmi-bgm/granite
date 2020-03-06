#!/usr/bin/env python

#################################################################
#
#    rckTar
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
import tarfile


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
    ttar = args['ttar']
    files = args['file']

    # Create index file
    with open(ttar + '.index', 'w') as fo:
        for file in files:
            fo.write(file.split('.')[0] + '\t' + file + '\n')
        # end for
    #end with

    # Create tar file
    tar_file = tarfile.open(ttar, 'w')

    for file in files:
        tar_file.add(file)
        tar_file.add(file + '.tbi')
    #end for

    tar_file.close()
#end def


#################################################################
#
#    MAIN
#
#################################################################
if __name__ == "__main__":

    main()

#end if
