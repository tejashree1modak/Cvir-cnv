# This script uses filtered duplication data obtained from filtration for cnv analyses

# Dependencies installation and loading
for (i in c("tidyverse","here")) {
  if(!require(i, character.only=TRUE)) {
    install.packages(i, dependencies=TRUE)
    require(i, character.only=TRUE)
  }
}

#### Copy Number comparison of Duplications across locations ####
# Read in file with filtered duplications. 
oysterdup_fil <- read.table(here("filtration/oysterdup_fil"), 
                            sep="\t" , stringsAsFactors = FALSE, header = TRUE)
#Getting the genotype and copy number for all pop
#function to pull out copy num from a col in the vcf for a sample
getcn <- function(bedout_col){
  str_split( bedout_col, ':') %>% map_chr(8)
}
cn_only <- map_dfr(select(oysterdup_fil,CL_1:SM_9),getcn)
cn_only$ID <- oysterdup_fil$ID
cn_only$POS <- oysterdup_fil$POS
cn_only$CHROM <- oysterdup_fil$CHROM
#Function to pull out genotype from a col in the vcf for a sample
getg <- function(bedout_col){
  str_split( bedout_col, ':') %>% map_chr(1)
}
gtypes_only <- map_dfr(select(oysterdup_fil,CL_1:SM_9),getg)
gtypes_only$ID <- oysterdup_fil$ID
#pulling out gtype and cn for each pop side by side to visually compare
gtypes_cn <- left_join(cn_only, gtypes_only, by = 'ID') 
gtypes_cn <- gtypes_cn[,order(colnames(gtypes_cn))]
#long table with both cn and gtype info
cn_long <- gather(cn_only,key=sample,value=cn,-ID, -POS, -CHROM)
cn_long$pop <- str_split(cn_long$sample,'_') %>% map(1) %>% as.character()
cn_long$pop <- as.vector(cn_long$pop)
#removing 0/0 and ./. genotypes since the cn is not real for those
#cn_gtypes_long <- left_join(cn_long, gtypes_long) %>% filter(gtype != "0/0" & gtype != "./.") %>% select(ID, CHROM, POS, pop, sample, cn, gtype, num_alts) 
cn_gtypes_long <- left_join(cn_long, gtypes_long)
#Converting cn value to 0 for genotypes 0/0 and ./. because they are assumed homologous to reference. 
cn_gtypes_long <- within(cn_gtypes_long, cn[gtype == '0/0'] <- 0)
cn_gtypes_long <- within(cn_gtypes_long, cn[gtype == './.'] <- 0)
cn_gtypes_long$cn <- as.numeric(as.character(cn_gtypes_long$cn))

# Plot cnv per location per individual per chromosome in one plot
#Chromosome 1
cn_gtypes_long_chr1_fil <- filter(cn_gtypes_long, CHROM == "NC_035780.1") %>% select(POS, sample, cn)
# Fig 3 from paper: Genome wide copy number variation profiles across 11 locations and 60 samples
cn_chr1_hmap_fil <- ggplot(data = cn_gtypes_long_chr1_fil, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 1") +
  scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10)) +
  theme(axis.text.x  = element_text(size=12), axis.text.y  = element_text(size=11), 
        axis.title.x  = element_text(face = "bold", size=16), 
        axis.title.y  = element_text(face = "bold", size=16))
cn_chr1_hmap_fil + geom_vline(xintercept = (65668439/2), color = "red", size=0.3) 
ggsave(filename = "Fig3.1.png", 
       path = here("characterization/figures"),
       width = 9, height = 8, units = "in" )

#Chromosome 2
cn_gtypes_long_chr2_fil <- filter(cn_gtypes_long, CHROM == "NC_035781.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr2_fil$cn <- as.numeric(as.character(cn_gtypes_long_chr2_fil$cn))
cn_chr2_hmap_fil <- ggplot(data = cn_gtypes_long_chr2_fil, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 2") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10)) +
  theme(axis.text.x  = element_text(size=12), axis.text.y  = element_text(size=11), 
        axis.title.x  = element_text(face = "bold", size=16), 
        axis.title.y  = element_text(face = "bold", size=16))
cn_chr2_hmap_fil + geom_vline(xintercept = (61752954/2), color = "red", size=0.3)
ggsave(filename = "Fig3.2.png", 
       path = here("characterization/figures"),
       width = 9, height = 8, units = "in" )

#Chromosome 3
cn_gtypes_long_chr3_fil <- filter(cn_gtypes_long, CHROM == "NC_035782.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr3_fil$cn <- as.numeric(as.character(cn_gtypes_long_chr3_fil$cn))
cn_chr3_hmap_fil <- ggplot(data = cn_gtypes_long_chr3_fil, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 3") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10)) +
  theme(axis.text.x  = element_text(size=12), axis.text.y  = element_text(size=11), 
        axis.title.x  = element_text(face = "bold", size=16), 
        axis.title.y  = element_text(face = "bold", size=16))
cn_chr3_hmap_fil + geom_vline(xintercept = (77061147/2), color = "red", size=0.3)
ggsave(filename = "Fig3.3.png", 
       path = here("characterization/figures"),
       width = 9, height = 8, units = "in" )

#Chromosome 4
cn_gtypes_long_chr4_fil <- filter(cn_gtypes_long, CHROM == "NC_035783.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr4_fil$cn <- as.numeric(as.character(cn_gtypes_long_chr4_fil$cn))
cn_chr4_hmap_fil <- ggplot(data = cn_gtypes_long_chr4_fil, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 4") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10)) +
  theme(axis.text.x  = element_text(size=12), axis.text.y  = element_text(size=11), 
        axis.title.x  = element_text(face = "bold", size=16), 
        axis.title.y  = element_text(face = "bold", size=16))
cn_chr4_hmap_fil + geom_vline(xintercept = (59691871/2), color = "red", size=0.3)
ggsave(filename = "Fig3.4.png", 
       path = here("characterization/figures"),
       width = 9, height = 8, units = "in" )

#Chromosome 5
cn_gtypes_long_chr5_fil <- filter(cn_gtypes_long, CHROM == "NC_035784.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr5_fil$cn <- as.numeric(as.character(cn_gtypes_long_chr5_fil$cn))
cn_chr5_hmap_fil <- ggplot(data = cn_gtypes_long_chr5_fil, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 5") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10)) +
  theme(axis.text.x  = element_text(size=12), axis.text.y  = element_text(size=11), 
        axis.title.x  = element_text(face = "bold", size=16), 
        axis.title.y  = element_text(face = "bold", size=16))
cn_chr5_hmap_fil + geom_vline(xintercept = (98698415/2), color = "red", size=0.3)
ggsave(filename = "Fig3.5.png", 
       path = here("characterization/figures"),
       width = 9, height = 8, units = "in" )

#Chromosome 6
cn_gtypes_long_chr6_fil <- filter(cn_gtypes_long, CHROM == "NC_035785.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr6_fil$cn <- as.numeric(as.character(cn_gtypes_long_chr6_fil$cn))
cn_chr6_hmap_fil <- ggplot(data = cn_gtypes_long_chr6_fil, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 6") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10)) +
  theme(axis.text.x  = element_text(size=12), axis.text.y  = element_text(size=11), 
        axis.title.x  = element_text(face = "bold", size=16), 
        axis.title.y  = element_text(face = "bold", size=16))
cn_chr6_hmap_fil + geom_vline(xintercept = (51258097/2), color = "red", size=0.3)
ggsave(filename = "Fig3.6.png", 
       path = here("characterization/figures"),
       width = 9, height = 8, units = "in" )

#Chromosome 7
cn_gtypes_long_chr7_fil <- filter(cn_gtypes_long, CHROM == "NC_035786.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr7_fil$cn <- as.numeric(as.character(cn_gtypes_long_chr7_fil$cn))
cn_chr7_hmap_fil <- ggplot(data = cn_gtypes_long_chr7_fil, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 7") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10)) +
  theme(axis.text.x  = element_text(size=12), axis.text.y  = element_text(size=11), 
        axis.title.x  = element_text(face = "bold", size=16), 
        axis.title.y  = element_text(face = "bold", size=16))
cn_chr7_hmap_fil + geom_vline(xintercept = (57830853/2), color = "red", size=0.3)
ggsave(filename = "Fig3.7.png", 
       path = here("characterization/figures"),
       width = 9, height = 8, units = "in" )

#Chromosome 8
cn_gtypes_long_chr8_fil <- filter(cn_gtypes_long, CHROM == "NC_035787.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr8_fil$cn <- as.numeric(as.character(cn_gtypes_long_chr8_fil$cn))
cn_chr8_hmap_fil <- ggplot(data = cn_gtypes_long_chr8_fil, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 8") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10)) +
  theme(axis.text.x  = element_text(size=12), axis.text.y  = element_text(size=11), 
        axis.title.x  = element_text(face = "bold", size=16), 
        axis.title.y  = element_text(face = "bold", size=16))
cn_chr8_hmap_fil + geom_vline(xintercept = (75944017/2), color = "red", size=0.3)
ggsave(filename = "Fig3.8.png", 
       path = here("characterization/figures"),
       width = 9, height = 8, units = "in" )

#Chromosome 9
cn_gtypes_long_chr9_fil <- filter(cn_gtypes_long, CHROM == "NC_035788.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr9_fil$cn <- as.numeric(as.character(cn_gtypes_long_chr9_fil$cn))
cn_chr9_hmap_fil <- ggplot(data = cn_gtypes_long_chr9_fil, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 9") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10)) +
  theme(axis.text.x  = element_text(size=12), axis.text.y  = element_text(size=11), 
        axis.title.x  = element_text(face = "bold", size=16), 
        axis.title.y  = element_text(face = "bold", size=16))
cn_chr9_hmap_fil + geom_vline(xintercept = (104168037/2), color = "red", size=0.3)
ggsave(filename = "Fig3.9.png", 
       path = here("characterization/figures"),
       width = 9, height = 8, units = "in" )

#Chromosome 10
cn_gtypes_long_chr10_fil <- filter(cn_gtypes_long, CHROM == "NC_035789.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr10_fil$cn <- as.numeric(as.character(cn_gtypes_long_chr10_fil$cn))
cn_chr10_hmap_fil <- ggplot(data = cn_gtypes_long_chr10_fil, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 10") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10)) +
  theme(axis.text.x  = element_text(size=12), axis.text.y  = element_text(size=11), 
        axis.title.x  = element_text(face = "bold", size=16), 
        axis.title.y  = element_text(face = "bold", size=16))
cn_chr10_hmap_fil + geom_vline(xintercept = (32650044/2), color = "red", size=0.3)
ggsave(filename = "Fig3.10.png", 
       path = here("characterization/figures"),
       width = 9, height = 8, units = "in" )


