## Install

A ready-to-use docker image is available to download.

    docker pull b3rse/granite:v0.1.1

To run locally, Python 3.6+ is required together with the following libraries:

  - [*numpy*](https://docs.scipy.org/doc/ "numpy documentation")
  - [*pysam*](https://pysam.readthedocs.io/en/latest/ "pysam documentation")
  - [*bitarray*](https://pypi.org/project/bitarray/ "bitarray documentation")
  - [*pytabix*](https://pypi.org/project/pytabix/ "pytabix documentation")
  - [*h5py*](https://www.h5py.org/ "h5py documentation")
  - [*matplotlib*](https://matplotlib.org/ "matplotlib documentation")

To install libraries with pip:

    pip install numpy pysam bitarray h5py matplotlib
    pip install --user pytabix

Additional software needs to be available in the environment:

  - [*samtools*](http://www.htslib.org/ "samtools documentation")
  - [*bgzip*](http://www.htslib.org/doc/bgzip.1.html "bgzip documentation")
  - [*tabix*](http://www.htslib.org/doc/tabix.1.html "tabix documentation")

To install the program from source:

    git clone https://github.com/dbmi-bgm/granite
    cd granite
    python setup.py install

To install the program with pip:

    pip install granite-suite
