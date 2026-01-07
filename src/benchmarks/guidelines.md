# Usage Guide

cuteHap requires phased alignment bam file and a corresponding reference to apply SV calling. Here, we provide a demo for implementing cuteHap, including sequencing acquiring, reads alignment, read phasing, and SV calling.

# Get tools

The installation of Clair3 (implement SNV calling) is referred [here](https://github.com/HKU-BAL/Clair3?tab=readme-ov-file#installation).
The installation of cuteHap can be applied through conda (see example below), pip, or source code.
```sh
conda create -n test_env python=3.10
conda activate test_env
conda install -c bioconda cutehap
```

# Get data
1) Download the alignment files:
```sh
wget -c https://downloads.pacbcloud.com/public/revio/2022Q4/HG002-rep1/analysis/HG002.m84011_220902_175841_s1.GRCh38.bam
samtools index HG002.m84011_220902_175841_s1.GRCh38.bam
```

2) Download GRCh38 reference:
```sh
wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
```

# Run read phasing
1) Run Clair3
```sh
run_clair3.sh --bam_fn=HG002.m84011_220902_175841_s1.GRCh38.bam --ref_fn=GRCh38_full_analysis_set_plus_decoy_hla.fa --threads=16 --platform=hifi --model_path=~/clair3/bin/models/hifi --output=clair3_output
```

2) Run Longphase
```sh
longphase phase -s clair3_output/merge_output.vcf.gz -b HG002.m84011_220902_175841_s1.GRCh38.bam -r GRCh38_full_analysis_set_plus_decoy_hla.fa -t 16 -o longphase -pb
longphase haplotag -r GRCh38_full_analysis_set_plus_decoy_hla.fa -s longphase.vcf -b HG002.m84011_220902_175841_s1.GRCh38.bam -t 16 -o longphase.tag
samtools index longphase.tag.bam
```

# Run cuteHap
```sh
cuteHap longphase.tag.bam GRCh38_full_analysis_set_plus_decoy_hla.fa cuteHap.vcf ./
```