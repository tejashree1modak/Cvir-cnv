# This script uses filtered duplication data obtained from filtration for annotation of duplications

# Dependencies installation and loading
for (i in c("tidyverse","here")) {
  if(!require(i, character.only=TRUE)) {
    install.packages(i, dependencies=TRUE)
    require(i, character.only=TRUE)
  }
}

# Read in file with filtered duplications. 
oysterdup_fil <- read.table(here("filtration/oysterdup_fil"), 
                            sep="\t" , stringsAsFactors = FALSE, header = TRUE)


#### Getting annotation from reference genome and mapping DUPs ####
# Make a bedfile of reference annotation
# GFF3 to BED
#Grep for F3 = gene and then cut F9 to keep just the gene= LOC. and created a bed file. 
#Oyster_gene_sort.bed with chr, st, stop and gene_name. 
# Annotation of duplications using bedtools 
# ` bedtools intersect -wa -wb ` > Oyster_Dup_gene

# Extract LOC and annotation from C. virginica gff3 file 
# awk -v FS=';' -v OFS='\t' ' /gene=/ && /product=/ { gene = product = ""; for (i=1; i <= NF ;i++) { if ($i ~ /^gene=/) { gene = substr($i, 6); } else if ($i ~ /product=/) { product = substr($i, 9); } } if (length(gene) && length(product)) { print gene, product } } ' ref_C_virginica-3.0_top_level.gff3 > ref_annot
# Annotating dups 
# Read in annotations from ref genome
ref_annot <- read.table("annotation/ref_annot", 
                        sep="\t" , quote="", fill=FALSE, stringsAsFactors = FALSE)
colnames(ref_annot) <- c("LOC", "annot")
# Read in bed file of dups mapped to LOCs from ref genome
map_dup <- read.table("annotation/Oyster_Dup_gene", sep="\t" , stringsAsFactors = FALSE)
colnames(map_dup) <- c("ID", "LOC")
# Keep those duplications that passed the filter
map_dup_fil <- map_dup %>% filter(map_dup$ID %in% oysterdup_fil$ID)
# Get annotation for filtered duplications
dup_annot <- left_join(map_dup_fil, ref_annot, by = "LOC") 
# write annotation file
dup_annot %>%  
  write.table(here("annotation/dup_annot"), append = FALSE, sep = "\t",quote = TRUE,
              row.names = F, col.names = TRUE)

#### KEGG annotation for duplications ####
#Extract LOC and XP_ID (protein_id) from gff3
# Annotation for C.virginica genome in the GFF3 format can be found at https://www.ncbi.nlm.nih.gov/genome/398
#awk -v FS=';' -v OFS='\t' ' /gene=/ && /protein_id=/ { gene = protein_id = ""; for (i=1; i <= NF ;i++) { if ($i ~ /^gene=/) { gene = substr($i, 6); } else if ($i ~ /protein_id=/) { protein_id = substr($i, 12); } } if (length(gene) && length(protein_id)) { print gene, protein_id } } ' ref_C_virginica-3.0_top_level.gff3 | uniq > ref_annot_prot
ref_annot_prot <- read.table("annotation/ref_annot_prot", 
                            sep="\t" , quote="", fill=FALSE, stringsAsFactors = FALSE)
colnames(ref_annot_prot) <- c("LOC", "Sequence_name")
#Join this with DUP_IDs (map_dup has DUP_ID and LOC)
dup_loc_xp <- left_join(map_dup_fil, ref_annot_prot, by="LOC") 
#Join that with XP_sequences_Cvirginica_GCF_002022765.2_GO.tab by XP 
ref_annot_go_kegg <- read.table(here("annotation/XP_sequences_Cvirginica_GCF_002022765.2_GO.tab"), 
                                sep="\t" , quote="", fill=FALSE, stringsAsFactors = FALSE, header = TRUE)
colnames(ref_annot_go_kegg) <- c("Sequence_name","Sequence_length","Sequence_description","GO_ID","Enzyme_code","Enzyme_name")
#What % genes (LOCs) are mapped to kegg or GO ID.
ref_annot_go_kegg %>% filter(!is.na(Enzyme_name)) %>% filter(Enzyme_name != "") %>% nrow() #10.4 % (100*6273)/60213
ref_annot_go_kegg %>% filter(!is.na(GO_ID)) %>% filter(GO_ID != "") %>% nrow() #58.2% (100*35081)/30213

dup_kegg <- left_join(dup_loc_xp, ref_annot_go_kegg, by="Sequence_name") %>% select(ID, Enzyme_code, Enzyme_name) %>% unique() 
# write KEGG annotation file
dup_kegg %>%  
  write.table(here("annotation/dup_kegg"), append = FALSE, sep = "\t",quote = TRUE,
              row.names = F, col.names = TRUE)
#What % dups mapped to an EC number via kegg
dup_kegg %>% filter(!is.na(Enzyme_name)) %>% filter(Enzyme_name != "") %>% nrow() # 10.8% (1144*100)/10599
#separate the enzyme names and get count for each
kegg_vector <- as.data.frame(table(unlist(strsplit(as.character(dup_kegg$Enzyme_name), ";"))))
colnames(kegg_vector) <- c("Enzyme Name", "Number of genes mapped")
kegg_vector_sorted <-  kegg_vector[order(kegg_vector$`Number of genes mapped`, decreasing=TRUE),] 
kegg_vector2 <- dup_kegg %>% filter(Enzyme_name == "Acting on peptide bonds (peptidases)") 
kegg_vector2 <- as.data.frame(table(unlist(strsplit(as.character(kegg_vector2$Enzyme_code), ";"))))
colnames(kegg_vector2) <- c("Enzyme Code", "Number of genes mapped")
kegg_vector2 <-  kegg_vector2[order(kegg_vector2$`Number of genes mapped`, decreasing=TRUE),] 
#top3 enzymes EC:3.4.24 = metalloendopeptidases,EC:3.4.22=cysteine endopeptidases ,EC:3.4.21 = serine endopeptidases

#make a csv file for paper
write.table(kegg_vector_sorted, here("annotation/dup_kegg_freq"), append = FALSE, sep = ",", quote = FALSE,
            row.names = F, col.names = TRUE)


#### GO annotation for duplications ####
#Extract DUP_IDs, GO_IDs
dup_go <- left_join(dup_loc_xp, ref_annot_go_kegg, by="Sequence_name") %>% select(ID, GO_ID) %>% unique() 
#What % dups mapped to GO terms
dup_go %>% filter(!is.na(GO_ID)) %>% filter(GO_ID != "") %>% nrow() # 63.45% (6726*100)/10599
# write GO annotation file
dup_go %>%  
  write.table(here("annotation/dup_go"), append = FALSE, sep = "\t",quote = TRUE,
              row.names = F, col.names = TRUE)
#separate the GO_IDs and get count for each
go_vector <- as.data.frame(table(unlist(strsplit(as.character(dup_go$GO_ID), ";"))))
go_vector_sorted <-  go_vector[order(go_vector$Freq, decreasing=TRUE),] 
#write
#write.table(go_vector_sorted, (here("characterization/go_vector_sorted.txt"), append = FALSE, sep = " ",quote = FALSE,
#            row.names = F, col.names = TRUE)

# Gene Ontology enrichment analysis for all duplications was done using topGO in topGO.R script
#write files for topGO: 
#interestinggenes.txt
left_join(dup_loc_xp, ref_annot_go_kegg, by="Sequence_name") %>% 
  select(Sequence_name) %>% unique() %>% 
  write.table(here("annotation/interestinggenes.txt"), append = FALSE, 
              sep = "\t",quote = F,row.names = F, col.names = F)
#Gene Universe
ref_annot_go_kegg %>% select(Sequence_name,GO_ID) %>% 
  write.table(here("annotation/annotations.txt"), append = FALSE, sep = "\t",quote = F,row.names = F, col.names = F)
#Use sed 's/;/,/g' annotations.txt > annotations_topgo.txt to convert it to format required by topGO
# USe topgo_pipe.R to get enrichment results
