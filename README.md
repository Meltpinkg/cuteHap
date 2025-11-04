# cuteHap

A haplotype-resolved SV detector in phased long read sequencing data.

---
### Installation
We recommand the installation via conda or pip:
```
    The installation via conda and pypi is now in construction.
```
Alternatively, installing from source with:
```
    git clone https://github.com/Meltpinkg/cuteHap.git
    cd cuteHap
    CFLAGS="-std=c99" python setup.py install
```

---	
### Introduction
cuteHap detects germline or mosaic structural variations (SVs) through phased long reads sequencing alignments. 

---
### Dependence
	
	1. python3
 	2. scipy
	2. pysam
	3. Biopython
	4. cigar
	5. numpy
	6. Cython

---
### Quick Start
```
cuteHap <phased.bam> <reference.fa> <output.vcf> <workdir>
```

---
### Changelog
    cuteHap (v1.0.0)
    1. the initial version of cuteHap

---
### Contact
For advising, bug reporting and requiring help, please post on [Github Issue](https://github.com/Meltpinkg/cuteHap/issues) or contact sqcao@stu.hit.edu.cn.
