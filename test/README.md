# Test Data

The example reads are the first 40000 lines of SRR6357070 (yeast). The genome and gtf are from nf-core rnaseq.
CI tests need to be run locally as we did not manage to build a STAR index small enough to fit into GitHub.
Use the below command for the indexing.

**Indexing**

```bash

STAR-2.7.5a/bin/Linux_x86_64/STAR --runMode genomeGenerate --runThreadN 1 --genomeDir index --genomeFastaFiles genome.fa --sjdbGTFfile genes.gtf --sjdbOverhang 149

```

**Testing**

Run this test after `cd`ing into the `$baseDir` of this repo:

```bash

nextflow run -profile docker,test main.nf

```
