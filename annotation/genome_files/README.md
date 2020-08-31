### To generate ref_C_virginica-3.0_top_level.fa
To conserve space the genome fasta file was split in two compressed files. For use do the following:

```shell
zcat cvir_chr1-chr5.fa.gz cvir_chr6-mito.fa.gz > ref_C_virginica-3.0_top_level.fa
```

Keep the unzipped file in this dir and it will be used by Oyster_Annotation_Processing.R as one of the input files. 
