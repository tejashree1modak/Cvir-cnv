Oyster Genome Bioinformatics
============================
Jon Puritz
12/7/2019


This will walk through all code and outputs for generating BAM files for the Eastern Oyster Genome Project.

``` r
setwd("~/Oyster_Genome_Project")
```

Initial read trimming and mapping
---------------------------------

For this, I will be using a modified version of `dDocent` [LINK](./Scripts/dDocent_ngs.sh).  This version of `dDocent`utilized [Trimmomatic](https://github.com/timflutre/trimmomatic) for read trimming, [BWA](https://github.com/lh3/bwa) for read mapping, [PICARD](https://broadinstitute.github.io/picard/) for duplicate dmarking, and [SAMTOOLS](https://github.com/samtools/samtools) for BAM file manipulation.  

Here are the first parameters used:

``` bash
cat config
```

    ## Number of Processors
    ## 40
    ## Maximum Memory
    ## 0
    ## Trimming
    ## no
    ## Assembly?
    ## no
    ## Type_of_Assembly
    ## 
    ## Clustering_Similarity%
    ## 
    ## Mapping_Reads?
    ## yes
    ## Mapping_Match_Value
    ## 1
    ## Mapping_MisMatch_Value
    ## 3
    ## Mapping_GapOpen_Penalty
    ## 5
    ## Calling_SNPs?
    ## no
    ## Email
    ## jpuritz@gmail.com

``` bash
./dDocent_ngs.sh config
```

This code is not run directly in the markdown notebook, but here is what the output would look like:

    dDocent run started Tue Apr 17 11:02:56 EDT 2018

    At this point, all configuration information has been entered and dDocent may take several hours to run.
    It is recommended that you move this script to a background operation and disable terminal input and output.
    All data and logfiles will still be recorded.
    To do this:
    Press control and Z simultaneously
    Type 'bg' without the quotes and press enter
    Type 'disown -h' again without the quotes and press enter

    Now sit back, relax, and wait for your analysis to finish.
    Trimming reads

    dDocent run started Tue Apr 17 11:17:32 EDT 2018

    At this point, all configuration information has been entered and dDocent may take several hours to run.
    It is recommended that you move this script to a background operation and disable terminal input and output.
    All data and logfiles will still be recorded.
    To do this:
    Press control and Z simultaneously
    Type 'bg' without the quotes and press enter
    Type 'disown -h' again without the quotes and press enter

    Now sit back, relax, and wait for your analysis to finish.
    Trimming reads
    Using BWA to map reads.
    [bam_sort_core] merging from 0 files and 40 in-memory blocks...
    [bam_sort_core] merging from 0 files and 40 in-memory blocks...
    [bam_sort_core] merging from 0 files and 40 in-memory blocks...
    [bam_sort_core] merging from 0 files and 40 in-memory blocks...
    [bam_sort_core] merging from 0 files and 40 in-memory blocks...
    [bam_sort_core] merging from 0 files and 40 in-memory blocks...
    [bam_sort_core] merging from 0 files and 40 in-memory blocks...
    [bam_sort_core] merging from 0 files and 40 in-memory blocks...
    [bam_sort_core] merging from 0 files and 40 in-memory blocks...

Mark Duplicates
---------------

These are PCR-free libraries, but many process generate duplicate sequences. This step is likely overly conservative, so I am open to suggestions for altering it, such as only removing optical duplicates. For Picard, I am using the recommended distance of 2500 for a patterned flowcell, like the HiSeq X 10.

``` bash
cat picard.sh
```

    ## #!/usr/bin/env bash
    ## 
    ## cat namelist | parallel -j 8 "java -Xms4g -jar /usr/local/bin/picard.jar MarkDuplicates I={}-RG.bam O={}-RGmd.bam M={}_dup_metrics.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 TAGGING_POLICY=OpticalOnly &> md.{}.log"
    ## 
    ## echo -e "Picard has finished  in" `pwd` | mailx -s "Analysis has finished" jpuritz@uri.edu

``` bash
bash picard.sh
```

Now that duplicates are marked. Duplicate sequences and secondary alignments are removed.

``` bash
filter_bam(){
samtools view -@32 -h -F 0x100 -F 0x400 $1-RGmd.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/'| samtools view -@ 32 -b
}

export -f filter_bam
cat namelist | parallel -j 8 "filter_bam {} > {}.F.bam"
```

Now, merge all the individual bam files into a single bam file for easier parallelized variant calling

``` bash
ls *.F.bam > bam.list
samtools merge -@64 -f filter.merged.bam -b bam.list 
samtools index -@64 filter.merged.bam
```

