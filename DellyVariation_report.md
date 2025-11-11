# DellyVariation Trio SV Report

## Summary

- Total variants: **17311**

  - DEL: 11823
  - INS: 5315
  - DUP: 173

## Variants by chromosome (top 10)

- chr2: 1393
- chr1: 1310
- chr4: 1236
- chr3: 1210
- chr5: 1146
- chr6: 1119
- chr7: 982
- chr10: 849
- chr8: 847
- chr12: 834

## Trio genotype pattern counts (top 10)

- HG00512:0/0|HG00513:0/0|HG00514:0/0: 9136
- HG00512:0/0|HG00513:0/1|HG00514:0/0: 1047
- HG00512:0/0|HG00513:0/1|HG00514:0/1: 1042
- HG00512:0/1|HG00513:0/0|HG00514:0/0: 1036
- HG00512:0/1|HG00513:0/0|HG00514:0/1: 994
- HG00512:0/1|HG00513:0/1|HG00514:0/1: 685
- HG00512:1/1|HG00513:1/1|HG00514:1/1: 595
- HG00512:1/1|HG00513:0/0|HG00514:0/1: 391
- HG00512:0/0|HG00513:1/1|HG00514:0/1: 381
- HG00512:0/1|HG00513:0/1|HG00514:0/0: 370

## De novo candidates assuming **HG00512** is the child

- Count: **159**

  - chr1:40501100-40504900 DEL (3800 bp), HG00512 0/1 (DV=92.0, GQ=99.0); parents: HG00513 0/0, HG00514 0/0
  - chr1:63239685-63242476 DEL (2791 bp), HG00512 0/1 (DV=69.0, GQ=99.0); parents: HG00513 0/0, HG00514 0/0
  - chr1:107860349-107862863 DEL (2514 bp), HG00512 0/1 (DV=91.0, GQ=99.0); parents: HG00513 0/0, HG00514 0/0
  - chr1:153071218-153094204 DUP (22986 bp), HG00512 0/1 (DV=29.0, GQ=99.0); parents: HG00513 0/0, HG00514 0/0
  - chr1:242964469-242965914 DEL (1445 bp), HG00512 0/1 (DV=13.0, GQ=99.0); parents: HG00513 0/0, HG00514 0/0

## De novo candidates assuming **HG00513** is the child

- Count: **178**

  - chr1:9786596-9796609 DEL (10013 bp), HG00513 0/1 (DV=49.0, GQ=99.0); parents: HG00512 0/0, HG00514 0/0
  - chr1:13384061-13385121 DEL (1060 bp), HG00513 0/1 (DV=37.0, GQ=99.0); parents: HG00512 0/0, HG00514 0/0
  - chr1:20605315-20606204 DEL (889 bp), HG00513 0/1 (DV=14.0, GQ=99.0); parents: HG00512 0/0, HG00514 0/0
  - chr1:54626600-54630294 DEL (3694 bp), HG00513 0/1 (DV=53.0, GQ=99.0); parents: HG00512 0/0, HG00514 0/0
  - chr1:56365434-56369290 DEL (3856 bp), HG00513 0/1 (DV=81.0, GQ=99.0); parents: HG00512 0/0, HG00514 0/0

## De novo candidates assuming **HG00514** is the child

- Count: **2**

  - chr13:105744786-105747125 DEL (2339 bp), HG00514 0/1 (DV=24.0, GQ=99.0); parents: HG00512 0/0, HG00513 0/0
  - chr14:105745717-105859942 DEL (114225 bp), HG00514 0/1 (DV=3.0, GQ=70.2); parents: HG00512 0/0, HG00513 0/0

## Notes on functional importance

- These are **structural variants (SVs)** called by DELLY (DEL/INS). Without a gene annotation file, we cannot assign exonic/intronic status here.

- **Deletions** can disrupt coding regions or regulatory elements; larger events have higher likelihood of functional impact.

- **Insertions** (especially long) may affect splicing or gene regulation.

- To prioritize, intersect these SVs with a gene annotation (e.g., GTF/BED) and known disease genes for your phenotype.

- Consider validating de novo candidates with an orthogonal method (e.g., PCR/long-read, or a second caller) and visually inspect in IGV.
