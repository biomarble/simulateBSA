# simulateBSA

A simulate script to generate an artificial SNP-index result  for Bulk Segregate Sequencing analysis.

Usage:
```bash
w=2000000   #sliding window size
s=100000    #sliding step

Rscript simulateBSA.r   Chr05  25000000   index.raw.txt
Rscript sliding.r index.raw.txt $w  $s 5  index.windowed.txt
```


