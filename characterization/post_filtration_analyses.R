# This script uses filtered duplication data obtained from filtration for further characterization


# Dependencies installation and loading
for (i in c("UpSetR","tidyverse","here")) {
  if(!require(i, character.only=TRUE)) {
    install.packages(i, dependencies=TRUE)
    require(i, character.only=TRUE)
  }
}

#### Genome coverage ####
# Merged duplications to avoid counting overlapping duplication multiple times
# The output of filter_dups.R is used as input for bedtools
# bedtools merge -i cvir_filtered_dups.bed -c 1,2,3 -o count,collapse,collapse  > characterization/cvir_filtered_dups_merged.bed
dups_fil_merged <- read.table(here("characterization/cvir_filtered_dups_merged.bed"), 
                              sep="\t" , stringsAsFactors = FALSE)
colnames(dups_fil_merged) <- c("CHROM", "POS","end","count","POS_collapse","end_collapse")
# Length of merged dups
dups_fil_merged$len <- (dups_fil_merged$end - dups_fil_merged$POS) + 1
#Total number of bases of all duplications
sum(dups_fil_merged$len) #112809562
# % of genome covered by duplications
# Genome size of C.virginica genome can be found at https://www.ncbi.nlm.nih.gov/genome/398
(sum(dups_fil_merged$len)/684000000)*100 #16.5% 
# basic stats
median(dups_fil_merged$len) #863
mean(dups_fil_merged$len) #16205.94
summary(dups_fil_merged$len)
# % below 1000bp
dups_fil_merged %>% filter(len < 1000) %>% nrow() #3715 ie 3715/10599 = 35% below 1000bp
#basic stats for LSDs (dups>1kb)
dups_fil_merged %>% filter(len > 1000) %>% nrow()
dups_fil_merged %>% filter(len > 1000) %>% select(len) %>% summary()


#### Genotypes for Duplication  ####
#Extracting genotype will inform us of presence/absence of a duplication per sample
# Genotype 0/0 = homozygous for absence of duplication, genotype 0/1 or 1/0 = heterozygouse for presence of duplication
# genotype 1/1 = homozygous for presence of duplication
gett <- function(bedout_col){
  sep_out <- str_split( bedout_col, ':') %>% as.data.frame()
  g <- sep_out[1,] %>% unname() %>% t()
  g2<-g[,1]
  table(g2)
}
#run above function on all samples (columns)
gtypesa <- select(oysterdup_fil,CL_1:HCVA_6) %>% map_dfr(gett)
gtypesb <- select(oysterdup_fil,HG_HG2F1:SM_9) %>% map_dfr(gett)
gtypesc <- gett(oysterdup_fil$HG_HG0F2)
gtypesc <- as.data.frame(gtypesc)
gtypesc <- as.data.frame(gtypesc$Freq)
absent_count <- c(as.integer(0)) # freq of genptype "./." for HG_HG0F2 is 0
gtypesc <- rbind(absent_count,gtypesc)
gtypes <- cbind(gtypesa, HG_HG0F2=gtypesc$`gtypesc$Freq` ,gtypesb)
gtypes2 <-as.data.frame(t(gtypes))
colnames(gtypes2) <- c("./.", "0/0",  "0/1",  "1/1")
gtypes2$pop <- map_chr(str_split(rownames(gtypes2),"_"),1) #add pop info (no sample num)
gtypes3 <- gather(gtypes2,genotype,number,-pop) #make long
#proportions
gtypes_p <- gtypes2 %>% mutate(sum=rowSums(select(gtypes2,("0/0":"1/1")))) %>%
  mutate(p0 = gtypes2$"0/0"/sum) %>% mutate(p01 = gtypes2$"0/1"/sum) %>%
  mutate(p1 = gtypes2$"1/1"/sum) %>% select(pop,p0,p01,p1)
gtypesp2 <- gather(gtypes_p,genotype,number,-pop)
gtypesp2$pop <- as.factor(gtypesp2$pop)
levels(gtypesp2$pop) <- c("HI","SM","HC","HCVA",  "CS", "CLP",
                          "SL","CL","LM",
                          "HG","NG")  #reorder pops
# Proportion of genotypes for duplications per population
#Plot shows heterozygous genotype (p01) for presence of duplications is more in proportion than homogzygous (p1) 
ggplot(gtypesp2,aes(genotype,number,color=pop))+geom_boxplot(outlier.size = 0.5) + labs(x="Genotype", y="Proportion") + theme_classic() +
  theme(axis.text.x  = element_text(size=12), axis.text.y  = element_text(size=12), 
        axis.title.x  = element_text(face = "bold", size=12), axis.title.y  = element_text(face = "bold", size=12),
        legend.title = element_blank())
ggsave(filename = "FigS2.png", 
       path = here("characterization/figures"),
       width = 6, height = 5, units = "in" )


#### Duplication lengths ####
oysterdup_fil <- read.table(here("filtration/oysterdup_fil"), 
              sep="\t" , stringsAsFactors = FALSE, header = TRUE)
# Fig S1 from paper: Frequency distribution of duplication lengths
# Figure shows shorter duplications more common than longer ones
ggplot(oysterdup_fil, aes(length))+geom_histogram(binwidth = 60,fill="steelblue")+ylim(c(0,1000))+
  xlim(c(0,5000)) + 
  labs(x="Length of duplications", y="Frequency") + theme_classic() +
  theme(axis.text.x  = element_text(size=12), axis.text.y  = element_text(size=12), axis.title.x  = element_text(face = "bold", size=12), axis.title.y  = element_text(face = "bold", size=12)) 
ggsave(filename = "FigS1.png", 
       path = here("characterization/figures"),
       width = 4, height = 3, units = "in" )
# Distribution of duplication lengths per population
pop_num_alts_present_fil <- read.table(here("filtration/pop_num_alts_present_fil"), 
              sep="\t" , stringsAsFactors = FALSE, header = TRUE)
#plot shows distribution of lengths is similar in all sampling locations
ggplot(pop_num_alts_present_fil, aes(pop,length)) +geom_violin() + ylim(c(0,2500))
# Mean lengths of duplications compared across populations
meanl <- group_by(pop_num_alts_present_fil,pop) %>% summarize(mean_len = mean(length),sd = sd(length))
ggplot(meanl,aes(pop,mean_len))+geom_point()+
  geom_errorbar(aes(ymin=mean_len+sd,ymax=mean_len-sd))


#### Duplication comparison across locations ####
## Upset plot of duplications across locations POST FILTERATION ##
# get a de-duplicated list of locus id's
ids <- unique(pop_num_alts_present_fil$ID)
# for each id, get a binary indicator of whether it is in a pop and bind to one dataframe
pops <- unique(pop_num_alts_present_fil$pop)
binaries <- pops %>% 
  map_dfc(~ ifelse(ids %in% filter(pop_num_alts_present_fil, pop == .x)$ID, 1, 0) %>% 
            as.data.frame) # UpSetR doesn't like tibbles
# set column names
names(binaries) <- pops
# have a look at the data
head(binaries)  
# how many duplications are present in more than 3 locations
filter(binaries,rowSums(binaries)>3) %>% nrow() #5954 ie 56%
# Fig 2b from paper: UpSet plot of the intersected duplications across locations
upset(binaries, nsets = length(pops), main.bar.color = "SteelBlue", sets.bar.color = "DarkCyan", 
      sets.x.label = "Number duplicate loci", text.scale = c(rep(1.4, 5), 2), order.by = "freq")

# Dups per location POST FILTRATION #
pop_sum_fil <- as.data.frame(colSums(binaries))
pop_sum_fil <- data.frame(pop = names(binaries),total_dups=colSums(binaries))
#get proportion of duplications by:
# dividing total duplications per location by total number of duplications across locations
pop_sum_fil$prop <- pop_sum_fil$total_dups/length(oysterdup_fil$ID)  #number of filtered dups are 10599 
# Fig: Proportion of duplications per location
ggplot(pop_sum_fil, aes(x=pop,y=prop, color=pop)) + geom_bar(stat = "identity", fill="white") + 
  labs(x="Populations", y="Proportion of total duplications per location") 


## Upset plot of duplications to compare hightemp and lowtemp among themselves POST FILTERATION ##
# For h_temp samples
h_temp <- c("CL", "LM", "SL")
pop_num_alts_present_fil_h_temp <- pop_num_alts_present_fil %>% 
  filter(pop %in% h_temp)
ids_h_temp <- unique(pop_num_alts_present_fil_h_temp$ID)
#Count number of duplications present in h_temp samples
length(ids_h_temp)#7731
binaries_h_temp <- h_temp %>% 
  map_dfc(~ ifelse(ids_h_temp %in% filter(pop_num_alts_present_fil_h_temp, pop == .x)$ID, 1, 0) %>% 
            as.data.frame) # UpSetR doesn't like tibbles
# set column names
names(binaries_h_temp) <- h_temp
# have a look at the data
head(binaries_h_temp)  
# UpSet plot of the intersected duplications across locations
upset(binaries_h_temp, nsets = length(h_temp), main.bar.color = "SteelBlue", sets.bar.color = "DarkCyan", 
      sets.x.label = "Number duplicate loci", text.scale = c(rep(1.4, 5), 1.4), order.by = "freq")

l_temp <- c("CLP","CS","HC","HCVA")
pop_num_alts_present_fil_l_temp <- pop_num_alts_present_fil %>% 
  filter(pop %in% l_temp)
ids_l_temp <- unique(pop_num_alts_present_fil_l_temp$ID)
#Count number of duplications present in l_temp samples
length(ids_l_temp)#7684 #8243 with CS
binaries_l_temp <- l_temp %>% 
  map_dfc(~ ifelse(ids_l_temp %in% filter(pop_num_alts_present_fil_l_temp, pop == .x)$ID, 1, 0) %>% 
            as.data.frame) # UpSetR doesn't like tibbles
# set column names
names(binaries_l_temp) <- l_temp
# have a look at the data
head(binaries_l_temp)  
# UpSet plot of the intersected duplications across locations
upset(binaries_l_temp, nsets = length(l_temp), main.bar.color = "SteelBlue", sets.bar.color = "DarkCyan", 
      sets.x.label = "Number duplicate loci", text.scale = c(rep(1.4, 5), 1.4), order.by = "freq")

## How mnay dups are present in both low and high temp samples? 
length(Reduce(intersect,list(ids_h_temp,ids_l_temp))) #5778 #5986 with CS
common_hl_temp <- Reduce(intersect,list(ids_h_temp,ids_l_temp))
grid.newpage()
draw.pairwise.venn(7731, 8243, 5986, category = c("High temp samples", "Low temp samples"), 
                   lty = rep("blank", 2), fill = c("sky blue", "orange"), 
                   alpha = rep(0.5, 2), cat.pos = c(340, 45), cat.dist = rep(0.025, 2), scaled = TRUE)
#ids in SL,SM (wild) and HG,NG (inbred)
rest_temp <- c("HG","NG","HI","SM")
pop_num_alts_present_fil_rest_temp <- pop_num_alts_present_fil %>% 
  filter(pop %in% rest_temp)
ids_rest_temp <- unique(pop_num_alts_present_fil_rest_temp$ID) 
length(ids_rest_temp) 
#rest and h_temp
length(Reduce(intersect,list(ids_rest_temp,ids_h_temp)))
#rest and l_temp
length(Reduce(intersect,list(ids_rest_temp,ids_l_temp)))
# rest_h_l_temp
length(Reduce(intersect,list(ids_rest_temp,ids_h_temp,ids_l_temp)))
#Triple Venn between high temp, low temp and the rest of the locations
grid.newpage()
venn.plot <- draw.triple.venn(
  #area1 = 7731,
  area1 = length(ids_h_temp),
  # area2 = 8243,
  area2 = length(ids_l_temp),
  # area3 = 7576,
  area3 = length(ids_rest_temp),
  # n12 = 5986,
  n12 = length(Reduce(intersect,list(ids_h_temp,ids_l_temp))),
  # n23 = 6274,
  n23 = length(Reduce(intersect,list(ids_rest_temp,ids_l_temp))),
  # n13 = 6501,
  n13 = length(Reduce(intersect,list(ids_rest_temp,ids_h_temp))),
  # n123 = 5575,
  n123 = length(Reduce(intersect,list(ids_rest_temp,ids_h_temp,ids_l_temp))),
  category = c("High", "Low", "Rest"),
  fill = c("dodgerblue", "goldenrod1", "seagreen3"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("dodgerblue", "goldenrod1", "seagreen3"))
#Get IDs that are unique to high and low temp
high_df <- as.data.frame(ids_h_temp)
colnames(high_df) <- c("ID")
low_df <- as.data.frame(ids_l_temp)
colnames(low_df) <- c("ID")
rest_df <- as.data.frame(ids_rest_temp)
colnames(rest_df) <- c("ID")
#get IDs unique to high temp
high_uniq1 <- anti_join(high_df,low_df)
high_uniq2 <- anti_join(high_uniq1,rest_df)
nrow(high_uniq2)  #1459
#get IDs unique to low temp
low_uniq1 <- anti_join(low_df,high_df)
low_uniq2 <- anti_join(low_uniq1,rest_df)
nrow(low_uniq2) #1322

## Comparison of FREQUENCY of dups between high and low temp
# High temp dups present in all high temp loc AND > 90% individuals
common_dups_htemp <- pop_num_alts_present_fil_h_temp %>% group_by(ID) %>% 
                          tally(sort = TRUE) %>% filter(n == 3) %>% select(ID) 
#Get dups common to all populations but present in all samples of each population
gtypes_long_htemp <- gtypes_long_fil %>% filter(!is.na(num_alts)) %>% filter(num_alts > 0) %>% 
                          filter(pop %in% h_temp)
#Out of the 3730 dups that are present in all high temp locations  how many are present in 
# Most frequent: >90% samples (i.e. in 15 samples) #325dups
most_freq_htemp <- semi_join(gtypes_long_htemp,common_dups_htemp, by="ID") %>% group_by(ID) %>% 
                      summarize(count=n()) %>% filter(count > 15) %>% select("ID")
# Least frequent: <10% samples (ie 1.7=2 individuals)
# no dups present in <10%, the least number of individuals is 4 hence using <5 as criteria, ie <30% #4dups
least_freq_htemp <- gtypes_long_fil %>% filter(!is.na(num_alts)) %>% filter(pop %in% h_temp) %>% group_by(ID) %>% 
                        summarize(count=n()) %>% filter(count < 5) %>% select("ID")
# Low temp 
common_dups_ltemp <- pop_num_alts_present_fil_l_temp %>% group_by(ID) %>% 
      tally(sort = TRUE) %>% filter(n == 4) %>% select(ID) 
#Get dups common to all populations but present in all samples of each population
gtypes_long_ltemp <- gtypes_long_fil %>% filter(!is.na(num_alts)) %>% filter(num_alts > 0) %>% 
      filter(pop %in% l_temp)
#Out of the 3806 dups that are present in all populations how many are present in 
# Most frequent: >90% samples (i.e. in 22 samples) #513dups
most_freq_ltemp <- semi_join(gtypes_long_ltemp,common_dups_ltemp, by="ID") %>% group_by(ID) %>% 
      summarize(count=n()) %>% filter(count > 22) %>% select("ID")
# Least frequent: <10% samples (ie 3 individuals)
# no dups present in <10%, the least number of individuals is 12 with 1 dup and 13 with 3 dups hence using <14 as criteria, ie <54% #4dups
least_freq_ltemp <- gtypes_long_fil %>% filter(!is.na(num_alts)) %>% filter(pop %in% l_temp) %>% group_by(ID) %>% 
      summarize(count=n()) %>% filter(count < 14) %>% select("ID")
# Most freq dups in high temp that are least freq in low temp
semi_join(most_freq_htemp,least_freq_ltemp) #1dup (not annotated)
# Most freq dups in low temp that are least freq in high temp
semi_join(most_freq_ltemp,least_freq_htemp) #1dup (not annotated)

#### Frequency of Duplications per chromosome comparison across locations ####
# get chromosome locations
chrom_pos <- oysterdup_fil %>% select(CHROM, POS)
# get chromosome information including lengths (can be obtained from https://www.ncbi.nlm.nih.gov/genome/398)
chrom_len <- data.frame(CHROM=c("NC_035780.1","NC_035781.1","NC_035782.1","NC_035783.1","NC_035784.1","NC_035785.1",
                                "NC_035786.1", "NC_035787.1","NC_035788.1","NC_035789.1"), 
                        start=c(1,1,1,1,1,1,1,1,1,1), 
                        end=c(65668440,61752955,77061148,59691872,98698416,51258098,57830854,75944018,104168038,32650045))
chrom_len$len <- chrom_len$end - chrom_len$start
#get frequency of duplications per location
gtypes_pos_fil <- map_dfr(select(oysterdup_fil,CL_1:SM_9),getg)
gtypes_pos_fil$POS <- oysterdup_fil$POS
gtypes_pos_long_fil <- gather(gtypes_pos_fil,key=sample,value=gtype,-POS)
gtypes_pos_long_fil$pop <- str_split(gtypes_pos_long_fil$sample,'_') %>% map(1) %>% as.character()
gtypes_pos_long_fil$pop <- as.vector(gtypes_pos_long_fil$pop)
gtypes_pos_long_fil$num_alts <- str_split(gtypes_pos_long_fil$gtype,'/') %>% 
  map(as.integer) %>% 
  map_int(sum)
#adding dups in all individuals of same pop to give pop count
pop_num_pos_alts_fil <- gtypes_pos_long_fil %>% filter(!is.na(num_alts)) %>%
  group_by(pop,POS) %>%
  summarize(num_alts = sum(num_alts))
pop_num_pos_alts_present_fil <- filter(pop_num_pos_alts_fil,num_alts >0)
#join to get chromosome information 
pop_num_pos_alts_present_chrom_fil <- left_join(pop_num_pos_alts_present_fil, chrom_pos, by = "POS")
pop_alts_per_chrom_fil <- pop_num_pos_alts_present_chrom_fil %>% group_by(pop,CHROM) %>% 
  summarize(num_alts = sum(num_alts)) 
pop_alts_per_chrom_len_fil <- left_join(pop_alts_per_chrom_fil, chrom_len, by = "CHROM")
pop_alts_per_chrom_len_fil$pop <- factor (as.character(pop_alts_per_chrom_len_fil$pop), 
                                          levels=c("HI","SM","CS","HC","HCVA","CLP","CL","SL","LM","HG","NG"))
# ANOVA for frequency of CNVs per chromosome
res_aov <- aov(num_alts ~ CHROM, data = pop_alts_per_chrom_fil)
summary(res_aov)
res_tuk <- TukeyHSD(res_aov)
#shows chr5 has significantly higher freq of cnv than any other chr
res_tuk_df <- as.data.frame(res_tuk$CHROM) 

# Fig2 in manuscript: Frequency of duplications per chromosome across locations normalized by chromosome length
freq_cnv <- pop_alts_per_chrom_len_fil %>% mutate(cnv_freq_norm = (num_alts/len)) %>% select(pop, CHROM, cnv_freq_norm)
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}
freq_cnv2 <- freq_cnv %>% group_by(CHROM) %>% 
  mutate(outlier = ifelse(is_outlier(cnv_freq_norm), cnv_freq_norm, as.numeric(NA))) 
freq_cnv2$pop[which(is.na(freq_cnv2$outlier))] <- as.numeric(NA)
freq_cnv2 <- freq_cnv2 %>%
  mutate(inbred_st = case_when(pop == 'HG' ~ 'inbred',
                               pop == 'NG' ~ 'inbred',
                               pop == 'CL' ~ 'nt_inbred',
                               pop == 'LM' ~ 'nt_inbred',
                               pop == 'SL' ~ 'nt_inbred'))
freq_cnv2$inbred_st <- as.factor(freq_cnv2$inbred_st)
freq_cnv3 <- freq_cnv2 %>% 
  mutate(outlier_inbred = case_when(inbred_st == 'inbred' ~ outlier, TRUE ~ NA_real_), 
         outlier_nt_inbred = case_when(inbred_st == 'nt_inbred' ~ outlier, TRUE ~ NA_real_))
ggplot(freq_cnv3, aes(x=CHROM,y=cnv_freq_norm)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(x=CHROM,y=outlier_inbred), shape=17, size=2)+
  geom_jitter(aes(x=CHROM,y=outlier_nt_inbred), shape=15, size=2)+
  labs(x="Chromosome Number", y="Frequency of CNVs") + theme_classic() +
  theme(axis.text.x  = element_text(size=12), axis.text.y  = element_text(size=12), 
        axis.title.x  = element_text(face = "bold", size=12), 
        axis.title.y  = element_text(face = "bold", size=12)) + 
  scale_x_discrete(labels=c("NC_035780.1"= "1","NC_035781.1"="2","NC_035782.1"="3","NC_035783.1"="4",
                            "NC_035784.1"="5","NC_035785.1"="6", "NC_035786.1"="7", "NC_035787.1"="8",
                            "NC_035788.1"="9","NC_035789.1"="10"))
ggsave(filename = "Fig2.png", 
       path = here("characterization/figures"),
       width = 4, height = 4, units = "in" )
