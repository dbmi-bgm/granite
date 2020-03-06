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
import ntpath


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
            filename = ntpath.basename(file)
            fo.write(filename.split('.')[0] + '\t' + filename + '\n')
        # end for
    #end with

    # Create tar file
    tar_file = tarfile.open(ttar, 'w')

    for file in files:
        filename = ntpath.basename(file)
        tar_file.add(file, filename)
        tar_file.add(file + '.tbi', filename + '.tbi')
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
