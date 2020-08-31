#This script filters duplicatons identified by delly using the following criteria:
  # Duplications that pass quality filter as applied in DELLY.
  # Duplications that are likely fixed in the population are filtered out.   
  # Duplications that overlap >10% with a repeat region are filtered out. 

# Dependencies installation and loading
for (i in c("tidyverse" , "here")) {
  if(!require(i, character.only=TRUE)) {
    install.packages(i, dependencies=TRUE)
    require(i, character.only=TRUE)
  }
}

##################################################################################################
#Repeatmasker output was obtained from the ftp server where C.virginica genome files are at
#ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/022/765/GCA_002022765.4_C_virginica-3.0
#The file is GCA_002022765.4_C_virginica-3.0_rm.out.gz
#This was preprocessed to make a bedfile of repeats as mentioned in README.md
  
#Get overlap 
#Use bedtools merge to merge repeats as mentioned in README.md
##################################################################################################

#### Read in merged vcf output file from delly ####
#all vcf data for each individual for each duplication obtained from DELLY in a vcf format
oysterdup <- read.table(here("filtration/germline_nohead_dup.vcf"),stringsAsFactors = FALSE)
#Add header with samples names
header <- strsplit("CHROM POS ID      REF     ALT     QUAL    FILTER  
                   INFO    FORMAT  CL_1    CL_2    CL_3    CL_4    CL_5    CL_6    CLP_1   CLP_2   
                   CLP_3   CLP_4   CLP_5   CLP_6   CS_1    CS_2    CS_3    CS_5    CS_6    CS_7    
                   HC_1    HC_3    HC_4    HC_5    HC_6    HC_7    HCVA_1 HCVA_2 HCVA_3 HCVA_4 HCVA_5 HCVA_6 HG_HG0F2       
                   HG_HG2F1        HG_HG2M5
                   HI_1    HI_2    HI_3    HI_4    HI_5    HI_6    LM_1_pool       LM_3    LM_4    
                   LM_7    LM_8    NG_NH0H4        NG_NH2F6        NG_NH2F8        
                   NG_NH2M1        SL_1    SL_2    SL_3    SL_4    SL_5    SL_6
                   SM_10   SM_11   SM_12   SM_7    SM_8    SM_9", "\\s+")[[1]]
colnames(oysterdup)<-header

#### Filter 1: Keep duplications that pass quality criteria in DELLY ####
oysterdup <-dplyr::filter(oysterdup,FILTER=="PASS")
#Get length of each duplication
oysterdup$end <- str_split(oysterdup$INFO, ';') %>%
  map_chr(5) %>% str_split('=') %>% map_chr(2) %>% as.integer()
oysterdup$length <- oysterdup$end - oysterdup$POS

#### Extracting genotypes for duplications per sample ####
# Genotype 0/0 = homozygous for absence of duplication, genotype 0/1 or 1/0 = heterozygous for presence of duplication
# genotype 1/1 = homozygous for presence of duplication

#Function to pull out genotype from a col in the vcf for a sample
getg <- function(bedout_col){
  str_split( bedout_col, ':') %>% map_chr(1)
}

gtypes_only <- map_dfr(select(oysterdup,CL_1:SM_9),getg)
gtypes_only$ID <- oysterdup$ID
gtypes_long <- gather(gtypes_only,key=sample,value=gtype,-ID)
gtypes_long$pop <- str_split(gtypes_long$sample,'_') %>% map(1) %>% as.character()
gtypes_long$pop <- as.vector(gtypes_long$pop)
gtypes_long$num_alts <- str_split(gtypes_long$gtype,'/') %>% 
  map(as.integer) %>% 
  map_int(sum)

#adding dups in all individuals of same pop to give pop count
pop_num_alts <- gtypes_long %>% filter(!is.na(num_alts)) %>%
  group_by(pop,ID) %>% summarize(num_alts = sum(num_alts)) 
#join to get duplication information per population
pop_num_alts <- left_join(pop_num_alts,select(oysterdup,ID,length) )
#Pulling out duplications present in all populations
#when num_alts = genotype 0/0 meaning duplication is absent in the sample at that location
#when num_alts > 0 genptype 0/1 or 1/0 or 1/1 meaning duplication is present in the sample at that location
pop_num_alts_present <- filter(pop_num_alts,num_alts >0)

#### Filter 2: Fixed duplications ####
#Get dups common to all populations since they are likely artifacts 
common_dups <- pop_num_alts_present %>% group_by(ID) %>% tally(sort = TRUE) %>% head(1046) %>% select(ID)
#Get dups common to all populations but present in all samples of each population
sample_num_alts <- gtypes_long %>% filter(!is.na(num_alts)) %>% filter(num_alts >0)
#Criteria for filteration of common dups: 
#Out of the 1046 dups that are present in all populations how many are present in >90% samples (i.e. in 54 samples) 
common_filter_dups <- 
  semi_join(sample_num_alts,common_dups, by="ID") %>% group_by(ID) %>% summarize(count=n()) %>% filter(count > 54) %>% select("ID")

#### Filter 3: Duplications in repeat regions ####
# Read in bedtools Ouput of intersect between repeat regions in reference genome and duplications  
dup_repeat_overlap <- read.table(here("filtration/dup_repeat_merged_overlap_mod.bed"), 
                                 sep="\t" , stringsAsFactors = FALSE)
colnames(dup_repeat_overlap) <- c("CHROM", "POS","end","ID","R_POS","R_end","R_ID","l")
#Number of repeats mapped to each duplicate
dup_repeat_overlap %>% select("ID","l") %>% group_by(ID) %>% tally() %>% View()
#Total len of repeats included in dups
overlap_total_len <- dup_repeat_overlap %>% select("ID","l") %>% group_by(ID) %>% summarise(sum(l))
colnames(overlap_total_len) <- c("ID","total_len")
#Get % of dup len covered by repeats for each dup that overlaps with repeats
percent_overlap <- oysterdup %>% select(ID,length) %>% left_join(overlap_total_len,by = 'ID') %>% na.omit() 
percent_overlap$percent <- (percent_overlap$total_len/percent_overlap$length)*100
#dups with >10% repeat coverage
percent_overlap %>% filter(percent > 10) %>% nrow() #filter out 1641 dups
repeat_filter_dups <- percent_overlap %>% filter(percent > 10) %>% select("ID")

#### Get final set of duplications post filtration ####
#Combining list of dups to be filtered because they are shared among >90% samples or have >10% repeat coverage.
filter_dups <- rbind(common_filter_dups, repeat_filter_dups) %>% distinct() 
# Make a bedfile for filtered duplications for further analysis
cvir_dup_bed <- oysterdup %>% select(CHROM,POS,end,ID)
colnames(cvir_dup_bed) <- c("CHROM","start","stop","ID")
# Number of duplications post filtration
anti_join(cvir_dup_bed,filter_dups) %>% group_by(ID) %>% summarize(count=n()) %>% nrow() #10599 
cvir_dups_fil_bed <- anti_join(cvir_dup_bed,filter_dups)

# Output files:
# These files are input for further characterization of duplications
# by scripts in characterization dir. 
#Write the files if needed. Files already available in the dir for use.
# FILE 1: BEDFILE of filtered duplications
cvir_dups_fil_bed %>%
 write.table(here("filtration/cvir_filtered_dups.bed"), append = FALSE, sep = "\t",quote = FALSE,
             row.names = F, col.names = FALSE)
# FILE 2: VCF of filtered duplications
oysterdup_fil <- anti_join(oysterdup, filter_dups)
oysterdup_fil %>%
 write.table(here("filtration/oysterdup_fil"), append = FALSE, sep = "\t",quote = TRUE,
                                            row.names = F, col.names = TRUE)

# FILE 3: population counts of filtered duplications
pop_num_alts_present_fil <- anti_join(pop_num_alts_present,filter_dups)
pop_num_alts_present_fil %>%
  write.table(here("filtration/pop_num_alts_present_fil"), append = FALSE, sep = "\t",quote = FALSE,
                                         row.names = F, col.names = TRUE)
