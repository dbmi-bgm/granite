
0.4.0
=====
*First Official Release*
* Updated Dockerfile
* Added BED support for blackList module
* Improved error handling with custom error objects
  - *MissingTag* describes a missing tag or tag value
  - *MissingTagDefinition* describes a missing tag definition
  - *TagDefinitionError* describes a format error for a tag definition
  - *TagFormatError* describes a format error for a tag
  - *MissingIdentifier* describes a missing genotype identifier in the VCF file
  - *VcfFormatError* describes an error in the VCF format

0.3.0
=====
* Added this CHANGELOG.rst
* Updated to support Python 3.12
  - updated bitarray: ">=1.2.0" -> "^2.9.2"
  - updated h5py: ">=2.10.0" -> "^3.11.0"
  - updated matplotlib: "==3.3.4" -> "^3.9.0"
  - udpated numpy: ">=1.18.0" -> "^1.26.4"
  - updated pysam ">=0.15.0" -> "^0.22.1"
  - updated pytabix: ">=0.0.2" -> "^0.1"
* Supported Python versions now are: 3.9, 3.10, 3.11, 3.12
* Note specifically that Python 3.8 is no longer supported (end of life: October 2024)

0.2.0
=====
* Original
