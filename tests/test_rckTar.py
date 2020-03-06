#################################################################
#   Libraries
#################################################################
import sys, os
import pytest
import shutil
import filecmp
import tarfile
from granite.rckTar import (
                            main as main_rckTar
                            )


#################################################################
#   Tests
#################################################################
def test_run_rckTar():
    ''' '''
    # Variables
    args = {'ttar': 'tests/files/main_test.rck.tar',
            'file': ['tests/files/rck_tar/GAPFI9BBS1MK.rck.gz',
                     'tests/files/rck_tar/GAPFIZWAW2KS.rck.gz',
                     'tests/files/rck_tar/GAPFI487H4XM.rck.gz']}
    # Run
    main_rckTar(args)
    # Tests
    tf = tarfile.open('tests/files/main_test.rck.tar')
    tf.extractall(path='tests/files/rck_tar_out')
    for file in os.listdir('tests/files/rck_tar_out'):
        assert filecmp.cmp('tests/files/rck_tar_out/' + file, 'tests/files/rck_tar/' + file, shallow=False)
    #end for
    assert [row for row in open('tests/files/main_test.rck.tar.index')] == [row for row in open('tests/files/rck_tar.index')]
    # Clean
    shutil.rmtree('tests/files/rck_tar_out')
    os.remove('tests/files/main_test.rck.tar')
    os.remove('tests/files/main_test.rck.tar.index')
#end def
