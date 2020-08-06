## File formats
The program is compatible with standard BED, BAM and VCF formats (VCFv4.x).

### ReadCountKeeper (.rck)
RCK is a tabular format that allows to efficiently store counts by strand (ForWard-ReVerse) for reads that support REFerence allele, ALTernate alleles, INSertions or DELetions at CHRomosome and POSition. RCK files can be further compressed with *bgzip* and indexed with *tabix* for storage, portability and faster random access. 1-based.

Tabular format structure:

    #CHR   POS   COVERAGE   REF_FW   REF_RV   ALT_FW   ALT_RV   INS_FW   INS_RV   DEL_FW   DEL_REV
    13     1     23         0        0        11       12       0        0        0        0
    13     2     35         18       15       1        1        0        0        0        0

Commands to compress and index files:
```text
    bgzip PATH/TO/FILE
    tabix -b 2 -s 1 -e 0 -c "#" PATH/TO/FILE.gz
```

### BinaryIndexGenome (.big)
BIG is a hdf5-based binary format that stores boolean values for each genomic position as bit arrays. Each position is represented in three complementary arrays that account for SNVs (Single-Nucleotide Variants), insertions and deletions respectively. 1-based.

hdf5 format structure:

    e.g.
    chr1_snv: array(bool)
    chr1_ins: array(bool)
    chr1_del: array(bool)
    chr2_snv: array(bool)
    ...
    ...
    chrM_del: array(bool)

*note*: hdf5 keys are built as the chromosome name based on reference (e.g. chr1) plus the suffix specifying whether the array represents SNVs (_snv), insertions (_ins) or deletions (_del).

### Pedigree in JSON format
When the program requires pedigree information, the expected format is as follow:

    [
      {
        "individual": "NA12877",
        "sample_name": "NA12877_sample",
        "gender": "M",
        "parents": []
      },
      {
        "individual": "NA12878",
        "sample_name": "NA12878_sample",
        "gender": "F",
        "parents": []
      },
      {
        "individual": "NA12879",
        "sample_name": "NA12879_sample",
        "gender": "F",
        "parents": ["NA12878", "NA12877"]
      }
    ]

where `individual` is the unique identifier for member inside the pedigree, `sample_name` is the corresponding sample ID in VCF file, and `parents` is the list of unique identifiers for member parents if any.
