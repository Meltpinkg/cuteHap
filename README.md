# cuteHap

A haplotype-resolved SV detector in phased long read sequencing data.

---
### Installation
We recommand the installation via conda or pip:
```
    $ conda install -c bioconda cutehap
	or
	$ pip install cutehap
```
Alternatively, installing from source with:
```
    git clone https://github.com/Meltpinkg/cuteHap.git
    cd cuteHap
    CFLAGS="-std=c99" python setup.py install
```

---	
### Introduction
Structural variations (SVs), as a major category of genomic rearrangements, are capable of altering millions of nucleotides within the human genome. The detection of germline SVs and somatic mosaicism has emerged as a critical frontier in genomic research. Long-read sequencing technologies have demonstrated transformative potential in characterizing these variants. cuteHap is designed to produce high-quality, phased call sets for germline SV detection while simultaneously identifying low-frequency somatic mosaic events. The method delivers high-performance, haplotype-resolved SV detection and comprehensive detection of low-frequency mosaicism. A detailed usage guide documentation is available at [here].

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
