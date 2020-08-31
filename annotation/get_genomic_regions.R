#This script pulls out different genomic regions from the reference genome of C.virginica.
# Author: Robert Literman

rm(list=ls())
library(biomaRt)
library(bedr)
library(parallel)
library(testit)
library(reticulate)
library(tidyverse)

cdsIntron <- function(gene_id){
  library(tidyverse)
  library(bedr)
  library(testit)
  temp_cds_range <- cds_range_bed %>% filter(Description == gene_id)
  temp_cds <- cds_bed %>% filter(Description == gene_id)
  temp_gene <- gene_bed %>% filter(Description == gene_id)
  temp_genic <- genic_stuff_bed %>% filter(Description == gene_id)
  
  intron_bed <- temp_genic[0,]
  genic_other_bed <- temp_genic[0,]
  
  intron_bed$Annotation <- character()
  genic_other_bed$Annotation <- character()
  
  if(!(has_error(bedr(input = list(a = temp_cds_range, b = temp_cds), method = "subtract", params="",check.chr = FALSE,verbose = FALSE,check.zero.based = FALSE,check.valid = FALSE,check.merge = FALSE,check.sort = FALSE), silent = !interactive()))){
    intron_bed <- bedr(input = list(a = temp_cds_range, b = temp_cds), method = "subtract", params="",check.chr = FALSE,verbose = FALSE,check.zero.based = FALSE,check.valid = FALSE,check.merge = FALSE,check.sort = FALSE) %>%
      `row.names<-`(NULL) %>%
      mutate(Annotation = 'intron')
  }
  
  if(!(has_error(bedr(input = list(a = temp_gene, b = temp_genic), method = "subtract", params="",check.chr = FALSE,verbose = FALSE,check.zero.based = FALSE,check.valid = FALSE,check.merge = FALSE,check.sort = FALSE), silent = !interactive()))){
    genic_other_bed <- bedr(input = list(a = temp_gene, b = temp_genic), method = "subtract", params="",check.chr = FALSE,verbose = FALSE,check.zero.based = FALSE,check.valid = FALSE,check.merge = FALSE,check.sort = FALSE) %>%
      `row.names<-`(NULL) %>%
      mutate(Annotation = 'genic_other')
  }
  return(rbind(intron_bed,genic_other_bed))
}

firstSemiColonNoTag <- function(string_to_split){
  first_element <- str_split(string_to_split,";")[[1]][1]
  return(gsub(".*:","",first_element))
}

mergeBED <- function(anno){
  library(tidyverse)
  library(bedr)
  library(testit)
  anno_bed <- gene_id_bed %>% filter(Annotation == anno)
  anno_table <- table(anno_bed$Description)
  
  temp_merged <- anno_bed[0,]
  
  # Find Gene IDs with 1 entry (no need to merge/sort)
  solo_genes <- names(anno_table[anno_table == 1])
  
  # Find Gene IDs with > 1 entry (need to merge/sort)
  multi_genes <- names(anno_table[anno_table > 1])
  
  if(length(solo_genes)>0){
    temp_merged <- rbind(temp_merged,anno_bed %>% filter(Description %in% solo_genes)) # Add solo Gene IDs to merged BED dataset
  }
  
  if(length(multi_genes)>0){
    newly_merged <- anno_bed %>% filter(Description %in% multi_genes) %>%
      mutate(Description = paste(Description,Strand,Annotation,sep = ":")) %>% # BEDTools merge removes these columns, so we can store them here and retrieve them later
      select(chr,start,end,Description) %>%
      as.data.frame() %>%
      bedr.sort.region(check.chr = FALSE,verbose=FALSE,check.merge = FALSE,check.zero.based = FALSE,check.valid = FALSE) %>%
      bedr.merge.region(check.chr = FALSE,verbose = FALSE,check.zero.based = FALSE,check.valid = FALSE,check.sort = FALSE,stratify.by = 'Description') %>%
      `row.names<-`(NULL) %>%
      rename(Description=names) %>%
      mutate(Score=1000) %>%
      separate(Description,into = c('Description','Strand','Annotation'),sep = ":",remove = TRUE) %>%
      select(chr,start,end,Description,Score,Strand,Annotation)
    temp_merged <- rbind(temp_merged,newly_merged)
  }
  return(temp_merged)
}

intergenicBED <- function(chromosome){
  library(tidyverse)
  library(bedr)
  library(testit)
  temp_chr <- chr_bed %>% filter(chr == chromosome)
  temp_features <- all_features %>% filter(chr == chromosome) %>% select(-names)
  temp_intergenic <- bedr(input = list(a = temp_chr, b = temp_features), method = "subtract", params="",check.chr = FALSE,verbose = FALSE,check.zero.based = FALSE,check.valid = FALSE,check.merge = FALSE,check.sort = FALSE) %>%
    bedr.sort.region(check.chr = FALSE,verbose = FALSE) %>%
    `row.names<-`(NULL) %>%
    rownames_to_column(var='count') %>%
    mutate(Description = paste('intergenic_chr',chr,count,sep = "_")) %>%
    mutate(Score = 1000,Strand = "+",Annotation = 'intergenic') %>%
    select(chr,start,end,Description,Score,Strand,Annotation)
  return(temp_intergenic)
}

getGene <- function(string){
  string_list <- str_split(string,pattern = ";") %>% unlist()
  return(string_list[grep("gene=", string_list)])
}

getGene_smRNA <- function(string){
  string_list <- str_split(string,pattern = ";") %>% unlist()
  return(string_list[grep("product=", string_list)])
}

setwd('/home/ralubuntu/Work/Tejashree_Intron/')

cores <- parallel::detectCores()

if(cores == 1){
  cores <- 1
} else{
  cores <- cores-1
}

fasta <- 'data/ref_C_virginica-3.0_top_level.fa'
gff3 <- 'data/ref_C_virginica-3.0_top_level.gff3'

# Remove header rows and read in GFF3
raw_gff3 <- read_tsv(system(paste("grep -E -v '##|#!'",gff3,sep = " "),intern = TRUE),
                     col_types = 'ccciicccc',
                     col_names = c('Chr','Build','Annotation','Start','End','Dot1','Strand','Dot3','Description'))
all_anno <- unique(raw_gff3$Annotation)

# HAVE TO GENERATE CHROMOSOME LENGTHS FROM FASTA BECAUSE GFF MAY NOT CONTAIN ENTIRE CHROMOSOME
source_python('scripts/seq_length.py')
chr_bed <- get_length(fasta) %>%
  mutate(chr = as.character(chr),start=as.integer(start),end=as.integer(end)) %>%
  as.data.frame() %>%
  bedr.sort.region(check.chr = FALSE,verbose = FALSE,check.zero.based = FALSE,check.merge = FALSE) %>%
  `row.names<-`(NULL)

# Set desired annotation list, as well as those to be merged downstream as smRNA
desired_list <- c('gene','lnc_RNA','CDS','rRNA','tRNA')
smRNA_list <- c('rRNA','tRNA')

# Set file ID
file_id <- 'Oyster'

# Extract annotations from desired list
smRNA_gff3 <- raw_gff3 %>% filter(Annotation=='rRNA' | Annotation =='tRNA') %>%
  mutate(Description = getGene_smRNA(Description)) %>%
  mutate(Description = str_replace(Description,pattern = "product=",replacement = "gene=")) %>%
  mutate(Annotation = "smRNA")

filtered_gff3 <- raw_gff3 %>% filter(Annotation %in% desired_list & !Annotation %in% smRNA_list)  %>%
  mutate(Description = getGene(Description)) %>%
  rbind(smRNA_gff3)

gene_id_bed <- filtered_gff3 %>%
  mutate(start=as.integer(Start - 1)) %>% # BED files are 0-based, GFF files are 1-based
  rename(end = End,chr=Chr) %>% # For BEDTools
  mutate(Score=1000) %>% # Score column required to have ID column; Dummy value
  filter(!(is.na(Description))) %>% # Remove transcipts or proteins that have deprecated genes
  select(chr,start,end,Description,Score,Strand,Annotation) %>% as.data.frame()

# Create non-coding gene fragment

cds_list <- gene_id_bed %>% filter(Annotation == 'CDS') %>% pull(Description) %>% unique() # Get genes with CDS

nc_genes <- gene_id_bed %>% filter(Annotation == 'gene') %>% filter(!(Description %in% cds_list)) %>% pull(Description) %>% unique() # Get genes without CDS

gene_id_bed <- gene_id_bed %>%
  filter(Annotation == 'gene' & Description %in% nc_genes) %>%
  mutate(Annotation = 'nc_gene') %>%
  rbind(gene_id_bed %>% filter(!(Annotation == 'gene' & Description %in% nc_genes))) %>%
  as.data.frame() %>%
  bedr.sort.region(check.chr = FALSE, verbose=FALSE)  %>%
  `row.names<-`(NULL)

# Write clean BED file for each annotation type prior to merging

output_dir <- paste('output/unmerged/',file_id,sep = "")

annos <- unique(gene_id_bed$Annotation)
for(anno in annos){
  anno_bed <- gene_id_bed %>%
    filter(Annotation == anno) %>%
    mutate(Description = paste(Description,Annotation,sep="_")) %>% # BED Files lack the columns necessary to separate
    select(-Annotation) %>%
    as.data.frame() %>%
    bedr.sort.region(check.chr = FALSE,verbose = FALSE,check.zero.based = FALSE,check.valid = FALSE,check.merge = FALSE) %>%
    write_tsv(path = paste(output_dir,anno,'Unmerged.bed',sep = "_"),col_names = FALSE)
}

# Create merged BED file
# In this pipeline, all features of a type (ie. alternative transcripts of CDS or UTR) are merged within a gene_ID, by type.

clust_merge <- makeCluster(cores)
clusterExport(clust_merge,c("gene_id_bed","mergeBED"))
merged_bed <- parLapply(clust_merge,annos,function(x) mergeBED(x))
merged_bed <- bind_rows(merged_bed)
stopCluster(clust_merge)

# Generate intronic + 'Other' genic regions

cds_bed <- merged_bed %>% filter(Annotation=='CDS') %>% select(-Annotation) %>% as.data.frame() %>%
  bedr.sort.region(check.chr = FALSE,verbose = FALSE) %>%
  `row.names<-`(NULL)

cds_range_bed <- merged_bed %>% filter(Annotation=='CDS') %>% select(-Annotation) %>%
  mutate(grouper = Description) %>%
  group_by(grouper) %>%
  summarise(chr=unique(chr),start=min(start),end=max(end),Description=unique(Description),Score=1000,Strand=unique(Strand)) %>%
  select(chr,start,end,Description,Score,Strand) %>%
  as.data.frame() %>%
  bedr.sort.region(check.chr = FALSE,verbose = FALSE) %>%
  `row.names<-`(NULL)

gene_bed <- merged_bed %>% filter(Annotation=='gene') %>% select(-Annotation) %>% as.data.frame() %>%
  bedr.sort.region(check.chr = FALSE,verbose = FALSE) %>%
  `row.names<-`(NULL)

coding_genes <- gene_bed %>% pull(Description) %>% unique()

genic_stuff_bed <- rbind(cds_range_bed) %>%
  bedr.sort.region(check.chr = FALSE,verbose = FALSE,check.zero.based = FALSE,check.valid = FALSE,check.merge = FALSE) %>%
  `row.names<-`(NULL)

clust_intronic <- makeCluster(cores)
clusterExport(clust_intronic,c("cds_range_bed","cds_bed",'gene_bed','genic_stuff_bed','cdsIntron'))
intron_bed <- parLapply(clust_intronic,coding_genes,function(x) cdsIntron(x))

non_empty <- intron_bed[sapply(intron_bed, function(x) dim(x)[1]) > 0]

intron_bed <- bind_rows(non_empty)
stopCluster(clust_intronic)

merged_bed <- rbind(merged_bed,intron_bed) # Add introns to DF

# # Generate intronic (OLD)
#
# genic_stuff_bed <- rbind(cds_bed,threeUTR_bed,fiveUTR_bed) %>%
#   bedr.sort.region(check.chr = FALSE,verbose = FALSE,check.zero.based = FALSE,check.valid = FALSE,check.merge = FALSE) %>%
#   `row.names<-`(NULL)
#
# clust_intronic2 <- makeCluster(cores)
# clusterExport(clust_intronic2,c("gene_bed","genic_stuff_bed",'geneIntron'))
# old_intron_bed <- parLapply(clust_intronic2,coding_genes,function(x) geneIntron(x))
# old_intron_bed <- bind_rows(old_intron_bed)
# stopCluster(clust_intronic2)
# gc()
#
# merged_bed <- rbind(merged_bed,old_intron_bed) # Add introns to DF

# Generate intergenic regions

all_features <- merged_bed %>%
  select(-Annotation) %>%
  as.data.frame() %>%
  bedr.sort.region(check.zero.based = FALSE,check.chr = FALSE,verbose = FALSE,check.valid = FALSE,check.merge = FALSE) %>%
  bedr.merge.region(check.zero.based = FALSE,check.chr = FALSE,check.valid = FALSE,check.sort = FALSE,verbose = FALSE)

chrs <- unique(chr_bed$chr)

clust_intergenic <- makeCluster(cores)
clusterExport(clust_intergenic,c("all_features","chr_bed",'intergenicBED'))
intergenic_bed <- parLapply(clust_intergenic,chrs,function(x) intergenicBED(x))
intergenic_bed <- bind_rows(intergenic_bed)
stopCluster(clust_intergenic)

# Write clean BED file for each merged annotation type

merged_bed <- rbind(merged_bed,intergenic_bed) %>% as.data.frame() %>% # Add intergenic to DF
  `row.names<-`(NULL) %>%
  bedr.sort.region(check.chr = FALSE,verbose = FALSE,check.merge = FALSE,check.zero.based = FALSE) %>%
  `row.names<-`(NULL)

output_dir <- paste('output/merged/',file_id,sep = "")

# Convert chromosome names in Intergenic description
merged_bed <- merged_bed %>%
  mutate(Description = str_replace(Description,'NC_035780.1','1')) %>%
  mutate(Description = str_replace(Description,'NC_035781.1','2')) %>%
  mutate(Description = str_replace(Description,'NC_035782.1','3')) %>%
  mutate(Description = str_replace(Description,'NC_035783.1','4')) %>%
  mutate(Description = str_replace(Description,'NC_035784.1','5')) %>%
  mutate(Description = str_replace(Description,'NC_035785.1','6')) %>%
  mutate(Description = str_replace(Description,'NC_035786.1','7')) %>%
  mutate(Description = str_replace(Description,'NC_035787.1','8')) %>%
  mutate(Description = str_replace(Description,'NC_035788.1','9')) %>%
  mutate(Description = str_replace(Description,'NC_035789.1','10')) %>%
  mutate(Description = str_replace(Description,'NC_007175.2','MT'))


for(anno in unique(merged_bed$Annotation)){
  anno_bed <- merged_bed %>%
    filter(Annotation == anno) %>%
    mutate(Description = paste(Description,Annotation,sep="_")) %>%
    select(chr,start,end,Description,Score,Strand) %>%
    as.data.frame() %>%
    bedr.sort.region(check.chr = FALSE,verbose = FALSE,check.zero.based = FALSE,check.valid = FALSE,check.merge = FALSE) %>%
    write_tsv(path = paste(output_dir,anno,'Merged.bed',sep = "_"),col_names = FALSE)
}