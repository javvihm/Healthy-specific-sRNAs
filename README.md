# Healthy-specific-sRNAs
Scripts used in my project Discovery of Healthy-Specific Small  RNAs Using Reference-Free  Computational Approaches. They include scripts used for the downstream analysis of a novel miRNA.

- FindHealthySpecificRNAs.r. Last step of the SURFR pipeline. Adaptation for healthy-specific sRNA discovery.
- novelmiRNA_distribution.R. It builds a count distribution after grepping the forward and reverse sequence on raw sequencing data.
- cor_test.R. It downloads bulk transcriptomic data and perform Spearman correlation tests between all genes and the novel miRNA.
- get_coverage_scheme.py. It build the coverage scheme from bowtie output.
