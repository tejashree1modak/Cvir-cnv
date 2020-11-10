## Exceptional genome wide copy number variation in the eastern oyster (*C. virginica*). 
### Modak et al., 2020

#### Reference genome and resequencing data:
- [Reference genome](https://www.ncbi.nlm.nih.gov/genome/398)
- Resequence data: 60 individuals of *C.virginica* across 9 locations along the eastern coast of the United States were sampled. Details about processing Illumina reads from Whole Genome Sequencing of these samples to generate BAM files used in CNV discovery are available in dir [WGS-BAMs](https://github.com/tejashree1modak/Cvir-cnv/tree/master/WGS-BAMs).
#### CNV discovery
- **Software:** CNVs were called using 'germline SV calling' with default parameters in [DELLY2](https://github.com/dellytools/delly) (v0.7.8).
- **Requirements:** BAM files per sample and reference genome fasta. 
- All steps performed as described in the [Germline SV calling section](https://github.com/dellytools/delly#germline-sv-calling).
- Convert output file format from .bcf to .vcf using [bcftools](http://samtools.github.io/bcftools/bcftools.html#view). 

#### R scripts for data processing and analysis 
**Requirements:** Combined vcf output from DELLY for all samples. 

#### 1. Filtration

This step filters duplications identified by delly using the following criteria:

1. Duplications that pass the quality filter as applied in DELLY are kept.
2. Duplications that are present in >90% of samples hence likely fixed in the population are filtered out.   
3. Duplications that overlap >10% with a repeat region identified in the reference genome are filtered out. 

Script: filter_dups.R

We identified 19,315 duplications in total. Filters 1, 2 and 3 identified 6815, 288 and 1641 duplications to be removed respectively. Thus, we removed 8,716 unique duplications resulting in 10,599 duplications that were then analyzed for all downstream analysis as shown in the next sections.

#### 2. Characterization

This step is used to characterize duplications remaining after filtration per step 1 above. 

1. across locations
2. across genome

Script 1: post_filtration_analyses.R 

Script 2: cnv_analyses.R

#### 3. Annotation

This step annotates filtered duplications using the reference genome annotation.

1. GO and KEGG annotation
2. Mapping duplications to different features in the genome

Script 1: annotation_dups.R

Script 2: Oyster_Annotation_Processing.R
