## Install

A ready-to-use docker image is available to download.

    docker pull b3rse/granite:<version>

Run granite using docker as follow:

    docker run b3rse/granite:<version> granite <command> ...

If installed locally, additional software needs to be available in the environment to run certain commands:

  - [*samtools*](http://www.htslib.org/ "samtools documentation")
  - [*bgzip*](http://www.htslib.org/doc/bgzip.1.html "bgzip documentation")
  - [*tabix*](http://www.htslib.org/doc/tabix.1.html "tabix documentation")

To install the program from source:

    git clone https://github.com/dbmi-bgm/granite
    cd granite
    make configure && make build

To install the program with pip:

    pip install granite-suite
