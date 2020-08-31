## Exceptional genome wide copy number variation in the eastern oyster (*C. virginica*). 
### Modak et al., 2020

#### Steps to get GO and KEGG annotations for reference genome
- Author: [Kevin Johnson](https://github.com/drkmj)

#### 1. Download NCBI amino acid sequences

```shell
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/022/765/GCF_002022765.2_C_virginica-3.0/GCF_002022765.2_C_virginica-3.0_protein.faa.gz
gunzip GCF_002022765.2_C_virginica-3.0_protein.faa.gz
```

#### 2. split sequences into 10,000 sequence bins

File - [split_multifasta.pl](split_multifasta.pl)

```shell
perl split_multifasta.pl --input_file GCF_002022765.2_C_virginica-3.0_protein.faa \
                         --output_dir ./ -f CV_prot_id_aa  --seqs_per_file=10000
```

#### 3. Run each smaller file through interproscan

interproscan version 5.27-66.0 should be used, installed with the recommended Panther member database analysis package (v12.0).

```shell
interproscan-5.27-66.0/interproscan.sh -i CV_prot_id_aa1.fsa -b ./ -goterms -f TSV,XML,GFF3
interproscan-5.27-66.0/interproscan.sh -i CV_prot_id_aa2.fsa -b ./ -goterms -f TSV,XML,GFF3
interproscan-5.27-66.0/interproscan.sh -i CV_prot_id_aa3.fsa -b ./ -goterms -f TSV,XML,GFF3
interproscan-5.27-66.0/interproscan.sh -i CV_prot_id_aa4.fsa -b ./ -goterms -f TSV,XML,GFF3
interproscan-5.27-66.0/interproscan.sh -i CV_prot_id_aa5.fsa -b ./ -goterms -f TSV,XML,GFF3
interproscan-5.27-66.0/interproscan.sh -i CV_prot_id_aa6.fsa -b ./ -goterms -f TSV,XML,GFF3
interproscan-5.27-66.0/interproscan.sh -i CV_prot_id_aa7.fsa -b ./ -goterms -f TSV,XML,GFF3
```

- The XML files were subsequently imported into BLAST2GO (V5.0), annotations were merged with those found in GCF_002022765.2_C_virginica-3.0_genomic.gff.
- GO SLIM terms were added in Blast2GO using the run GO-Slim feature
- KEGG codes were added in Blast2GO using the GO-EnzymeCode Mapping feature.
- These final annotations were exported into a tab deliminated file.

