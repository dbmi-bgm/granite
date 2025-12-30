## Install

A ready-to-use Docker image is available to download [here](https://hub.docker.com/r/b3rse/granite "granite Docker").

    docker pull b3rse/granite:<version>

Run granite using Docker as follows:

    docker run b3rse/granite:<version> granite <command> ...

If installed locally, additional software needs to be available in the environment to run certain commands:

  - [*Samtools*](http://www.htslib.org/ "Samtools documentation")
  - [*bgzip*](http://www.htslib.org/doc/bgzip.1.html "bgzip documentation")
  - [*tabix*](http://www.htslib.org/doc/tabix.1.html "tabix documentation")

To install the program from source:

    git clone https://github.com/dbmi-bgm/granite
    cd granite
    make configure && make build

To install the program with pip:

    pip install granite-suite

### VEP

The comHet command requires gene- and transcript-level annotations produced by the Ensembl Variant Effect Predictor (VEP). granite does not run VEP directly; instead, it consumes VEP annotations already present in the input VCF file and uses them to assign variants to genes and transcripts and to evaluate inheritance patterns.

VEP annotations are expected to follow the standard VCF INFO field format (e.g. the *CSQ* tag), as produced by default VEP configurations. For details on running VEP and generating compatible annotations, please refer to the official documentation:
https://github.com/Ensembl/ensembl-vep
