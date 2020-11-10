## Exceptional genome wide copy number variation in the eastern oyster (*C. virginica*). 
### Modak et al., 2020

#### Duplication filtration:
#### Filter 1: Duplications that pass the quality filter as applied in DELLY are kept.
#### Filter 2: Duplications that are present in >90% of samples hence likely fixed in the population are filtered out.
#### Filter 3: Duplications overlapping with repeat regions identified in the reference genome
##### Duplications that overlap >10% with a repeat region were filtered from the analyses
- [Repeats identified in C.virginica genome](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/022/765/GCA_002022765.4_C_virginica-3.0)
Repeats file from NCBI was preprocessed to make a bedfile of repeats as follows:
Modify chromosome names to match gff3 to get Cvir_genome_repeats_mod.bed

```shell
 awk -v OFS='\t' 'BEGIN {
     change["CM008241.1"] = "NC_035780.1"
     change["CM008242.1"] = "NC_035781.1"
     change["CM008243.1"] = "NC_035782.1"
     change["CM008244.1"] = "NC_035783.1"
     change["CM008245.1"] = "NC_035784.1"
     change["CM008246.1"] = "NC_035785.1"
     change["CM008247.1"] = "NC_035786.1"
     change["CM008248.1"] = "NC_035787.1"
     change["CM008249.1"] = "NC_035788.1"
     change["CM008250.1"] = "NC_035789.1"
 }
 
 NR > 3 {
    if ($5 in change)  {
        $5 = change[$5]
    }
    print $5,$6,$7,$5"_"$6
 }' GCA_002022765.4_C_virginica-3.0_rm.out  > Cvir_genome_repeats_mod.bed
```

##### BEDTools used to obtain overlap between duplications and repeats
Use bedtools merge to merge repeats as follows:

```shell
bedtools merge -i Cvir_genome_repeats_mod.bed -c 1,5 -o count,collapse > Cvir_genome_repeats_merged.bed
```

Use bedtools using the -wo flag to obtain overlaps between merged repeats and duplications

```shell
bedtools intersect -a oysterduplicate_sort.bed -b Cvir_repeats_merged.bed -wo > dup_repeat_merged_overlap_mod.bed
```

This file was used to filter duplications that overlap >10% with a repeat region.
Use filter_repeats.R for this step.
This script produces the following files of filtered duplications that are used in ALL further analysis
- **FILE 1: BEDFILE of filtered duplications**
- **FILE 2: VCF of filtered duplications**
- **FILE 3: population counts of filtered duplications**
