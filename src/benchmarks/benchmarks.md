# Benchmark example

We provide a demo to help users of performing SV detection and evaluating performance. Here, we apply various SV detection methods on alignments from PacBio HiFi sequencing, and regard the "High confidence variant" of HG002 from T2T assembly as ground truth.

# Get tools

Information about how to install `conda` and add the `bioconda` channel is available on https://bioconda.github.io/.
Instructions of installing other state-of-the-art methods are on their GitHub page.

```sh
conda create -n test_env python=3.11
conda activate test_env
conda install -c bioconda sniffles==2.2 cuteSV==2.1.2 cutehap==1.0.2
```

# Get data
1) Activate conda environment:
```sh
conda activate test_env
```

2) Download T2T assembly ground truth:
```sh
wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/GRCh38_HG2-T2TQ100-V1.1_stvar.vcf.gz
wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed
```

3) Download GRCh38 reference:
```sh
wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
```

4) Download the alignment files:
```sh
wget -c https://downloads.pacbcloud.com/public/revio/2022Q4/HG002-rep1/analysis/HG002.m84011_220902_175841_s1.GRCh38.bam
samtools index HG002.m84011_220902_175841_s1.GRCh38.bam
```

5) Perform read phasing
```sh
run_clair3.sh --bam_fn=HG002.m84011_220902_175841_s1.GRCh38.bam --ref_fn=GRCh38_full_analysis_set_plus_decoy_hla.fa --threads=16 --platform=hifi --model_path=~/miniconda3/envs/clair3/bin/models/hifi --output=clair3_output
longphase phase -s clair3_output/merge_output.vcf.gz -b HG002.m84011_220902_175841_s1.GRCh38.bam -r GRCh38_full_analysis_set_plus_decoy_hla.fa -t 16 -o $clair3_output/longphase --pb
longphase haplotag -r GRCh38_full_analysis_set_plus_decoy_hla.fa -s clair3_output/longphase.vcf -b HG002.m84011_220902_175841_s1.GRCh38.bam -t 16 -o HG002.phased
samtools index HG002.phased.bam
```

# Run Tools

6) Run Sniffles2 (v2.2):
```sh
sniffles --input HG002.phased.bam --vcf sniffles2.vcf --reference GRCh38_full_analysis_set_plus_decoy_hla.fa
```

# Run cuteSV

7) Run cuteSV (v2.1.2):
```sh
cuteSV HG002.phased.bam GRCh38_full_analysis_set_plus_decoy_hla.fa cutesv.vcf ./ --genotype --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 -S HG002 -s 5 -l 50
```

# Run cuteHap

8) Run cuteHap (v1.0.2):
```sh
cuteHap HG002.phased.bam GRCh38_full_analysis_set_plus_decoy_hla.fa cutehap.vcf ./
```

# Final comparison

9) Compare to ground truth via truvari (v5.1.0):
```sh
for tools in {sniffles2,cutesv,cutehap}
do
        bcftools sort ${tools}.vcf -Oz -o ${tools}.vcf.gz
        tabix ${tools}.vcf.gz
        truvari bench -b GRCh38_HG2-T2TQ100-V1.1_stvar.vcf.gz -c ${tools}.vcf.gz --includebed GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed -o cmp --pick ac --passonly -r 2000 -C 5000 -f GRCh38_full_analysis_set_plus_decoy_hla.fa
        truvari refine --use-original-vcfs --align mafft -f GRCh38_full_analysis_set_plus_decoy_hla.fa -r cmp/candidate.refine.bed cmp -t 16 --mafft-params "--auto --thread 16"
done
```