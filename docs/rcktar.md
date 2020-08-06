### rckTar
rckTar creates a tar archive from bgzip and tabix indexed RCK files. Creates an index file for the archive.

#### Arguments
```text
    usage: granite rckTar [-h] -t TTAR -f FILE

    optional arguments:
      -t TTAR, --ttar TTAR  target tar to write results, use .rck.tar as extension
      -f FILE, --file FILE  file to be archived. Specify multiple files as: "-f
                            SampleID_1.rck.gz -f SampleID_2.rck.gz -f ...". Files
                            order is maintained while creating the index
```
