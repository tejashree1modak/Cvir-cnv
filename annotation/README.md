## Exceptional genome wide copy number variation in the eastern oyster (*C. virginica*). 
### Modak et al.2020

### Duplication annotation:
#### Annotation files for the reference genome can be found here:
- [*C.virginica* genome](https://www.ncbi.nlm.nih.gov/genome/398)

oyster_genes.bed: BEDfile of genes obtained from the genome annotation 
oyster_dups.bed: BEDfile of duplications identified by DELLY
#### BEDTools used to annotate duplications 
Bedtools intersect was used to annotate duplications  
`bedtools intersect -wa -wb oyster_dups.bed oyster_genes.bed > Oyster_Dup_gene.bed`

This file is available in annotation dir. 
annotation_dups.R produces file (dup_annot) for duplication with their mapped annotation.

### KEGG and GO annotation using annotation_dups.R 
GO and KEGG mapping from the reference genome was used to map duplications to these terms. 

### Input files for  annot_dups.R
These files were made from reference annotation in GFF3 format as follows: 

Extract LOC and annotation

`awk -v FS=';' -v OFS='\t' ' /gene=/ && /product=/ { gene = product = ""; for (i=1; i <= NF ;i++) { if ($i ~ /^gene=/) { gene = substr($i, 6); } else if ($i ~ /product=/) { product = substr($i, 9); } } if (length(gene) && length(product)) { print gene, product } } ' ref_C_virginica-3.0_top_level.gff3 > ref_annot`

Extract LOC and protein ID

`awk -v FS=';' -v OFS='\t' ' /gene=/ && /protein_id=/ { gene = protein_id = ""; for (i=1; i <= NF ;i++) { if ($i ~ /^gene=/) { gene = substr($i, 6); } else if ($i ~ /protein_id=/) { protein_id = substr($i, 12); } } if (length(gene) && length(protein_id)) { print gene, protein_id } } ' ref_C_virginica-3.0_top_level.gff3 | uniq > ref_annot_prot`

Output: annotation_dups.R produces files (dup_kegg, dup_go) for duplications with their mapped KEGG or GO annotation.

### Mapping duplications to different features in the genome

Bedtools with -wo flag was used to obtain overlap of duplications with each type of genomic feature.
The bedtools output was processed as follows to obtain the number of duplicaitons completely within the genome feature:

`cut -f 4 -d ' ' <bedtools_ouput> | sort -n | uniq | wc -l`

#### Overlap with genes

Genes were pulled out from the reference.gff3 file:

`awk -F'\t' -v OFS='\t' '$3 == "gene" {print $1, $4, $5, $3, $9}' ref_C_virginica-3.0_top_level.gff3 | awk -F'\t' -v OFS='\t' '/gene=/ { split($NF, a, /;/); gene = ""; for (i in a)  { if (a[i] ~ /^gene=/) { gene = substr(a[i], 6) } }  print $1, $2, $3, $4, gene; next }  { print $1, $2, $3, $4, "NA" } ' > gene_annot.bed`

Result: 4270 (40.3%) duplications lie completely within genes

#### Overlap with exons

Exons were pulled out from the reference.gff3 file:

`awk -F'\t' -v OFS='\t' '$3 == "exon" {print $1, $4, $5, $3, $9}' ref_C_virginica-3.0_top_level.gff3 | awk -F'\t' -v OFS='\t' '/gene=/ && /product=/ { split($NF, a, /;/); gene = product = ""; for (i in a)  { if (a[i] ~ /^gene=/) { gene = substr(a[i], 6) } else if (a[i] ~ /^product=/) { product = substr(a[i], 9); } } print $1, $2, $3, $4, gene, product; next }  { print $1, $2, $3, $4, "NA", "NA" } '  > exon_annot.bed`

Result: 535 (5.05%) duplications lie completely within exons

#### Overlap with introns and intergenic regions separately

Introns and intergenic regions were pulled out separately for C.virginica reference genome. This was done using the script Oyster_Annotation_Processing.R present in this directory.
Bedtools intersect with -f1 flag was used to get duplications completely within introns or intergenic region.
Result: 2306 (21.8%) duplications overlap introns completely and 3105 (29.3%) overlap intergenic regions completely.

### GO enrichment 
GO enrichment was performed using topGO using script topGO.R present in this dir. This script performs enrichment using GO terms mapped to duplications as compared to GO mapping of the reference genome.  

