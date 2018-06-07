# Analysis of archaeal amoA gene amplicons

Written by alexandre.bagnoud@gmail.com and henri.siljanen@uef.fi in May 2018

This script was used for analyzing archaeal amoA amplicon sequences generated with the primers pair CamoaA-19F and TamoA-629R-2 and sequenced with Illumina MiSeq with a 2x300 bp configuration.
These data were published here:

*add here the reference of the article once published*

#### Software used:
* DADA2 v1.6.0 (https://benjjneb.github.io/dada2/index.html)
* QIIME1 v1.9.1 (http://qiime.org/)
* R v.3.4.4 and packages

#### Database used:
* Alves et al., 2018, Nature Communications (Supplementary Information)

#### 1) Preparing the data

* Place all the fastq.bz2 files in one folder '0/raw_data'. Used the the reads trimmed from any non-biological sequences.
* Rename them if needed.
* Convert the fastq files into fastq (on the terminal):

```
cd /Users/siljanen/Documents/MiSeq_data_LCG/LGC_G20002861_part1/PrimerClipped/Henri_AOA_amoA/Peat_soil_DNA/Raw_R1R2
bzip2 -d *.bz2
```
#### 2) DADA2 pipeline

This part is was written based on the [DADA2 tutorial](https://benjjneb.github.io/dada2/tutorial.html)

##### 2.1) Set-up the R environment
* Load DADA2
```
library("dada2")
```
* Set the working directory
```
setwd("/Users/siljanen/Documents/AA_MiSeq_data_LCG/LGC_G20002861_part1/PrimerClipped/Henri_AOA_amoA/Peat_soil_DNA/Raw_R1R2") 
```
* Define path variable for the fastq files
```
path <- "/Users/siljanen/Documents/AA_MiSeq_data_LCG/LGC_G20002861_part1/PrimerClipped/Henri_AOA_amoA/Peat_soil_DNA/Raw_R1R2"
list.files(path)
```

##### 2.2) Filtering and trimming
* Extracting sample names
```
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names
```
* Read quality vizualisation of some samples
```
plotQualityProfile(fnFs[1:12]) # Visualize the quality of other samples by modifying the numbers in bracket.
```
* Read filtering and trimming
```
filt_path <- file.path("1-filtered_reads")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))


out <- filterAndTrim(fnFs, filtFs,  truncLen=c(200),
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # Modify accordingly

out
```

##### 2.3) Learn the error rates
```
errF <- learnErrors(filtFs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
```

##### 2.4) Dereplication
```
derepFs <- derepFastq(filtFs, verbose=TRUE)
```
* Name the derep-class objects by the sample names
```
names(derepFs) <- sample.names
```
##### 2.5) Sample inferrence
```
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
```

## CONSTRUCT SEQUENCE TABLE

seqtab <- makeSequenceTable(dadaFs, derepFs)
dim(seqtab)
# Surface without Taz  = 12 207
# New: with Taz samples  = 18 265

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# 200 
#    207 

## With Taz:
#  200 
#      265 

## REMOVE CHIMERAS

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# Identified 47 bimeras out of 207 input sequences.
# With Taz  = Identified 55 bimeras out of 265 input sequences.
dim(seqtab.nochim)
# [1]  12 160
#  With Taz:  18 210

sum(seqtab.nochim)/sum(seqtab)
# 1] 0.9866975
# With Taz = 0.9888035



## TRACK THE READS THROUGH THE PIPELINE

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "tabled", "nonchim")
rownames(track) <- sample.names
track

dir.create("2-dada2")
write.table(track, "2-dada2/1-track.txt", quote = FALSE, sep = "\t", col.names = NA)

## EXPORT DATA

# fasta of uniques non-chimeric reads
uniquesToFasta(getUniques(seqtab.nochim), "2-dada2/2-uniques_nochim.fasta")

# ASV table
write.table(t(seqtab.nochim), "2-dada2/3-asv_table.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

## ANNOTATION OF UNIQUES READS WITH DATABASE (on the terminal)

#################################################################################################

# First: Install qiime1  (This protocol is for MacBook, PC might have different commands.)
## Install before miniconda 2.7.

# Download miniconda3 from website: 
##    https://conda.io/miniconda.html

cd /Users/siljanen/Downloads
bash Miniconda2-latest-MacOSX-x86_64.sh

# accept Licence agreement and answer "yes" to following question (this keeps miniconda3 active all the time at background):
# "Do you wish the installer to prepend the Miniconda3 install location
#  to PATH in your /Users/siljanen/.bash_profile ? [yes|no]"
#  If you answer "no", then you need activate miniconda each time you want use it. Like this: 
#  export PATH=/Users/siljanen/miniconda3/bin:$PATH


##  Now open new Terminal window, because this is need for activation of miniconda
# Then just install qiime1



conda create -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda

# To activate this environment, use:
# > source activate qiime1
#
# To deactivate an active environment, use:
# > source deactivate


# Optional way to install qiime1:

conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda

conda install qiime

#################################################################################################


"""
cd /Users/siljanen/Documents/AA_MiSeq_data_LCG/LGC_G20002861_part1/PrimerClipped/Henri_AOA_amoA/Peat_soil_DNA/Raw_R1R2/2-dada2/


# with usearch -pipeline: (from usearch_v5_aoa_18-04-2018.sh)
# Change to path to the database file accordingly
db_seq="/Users/siljanen/Documents/AA_MiSeq_data_LCG/ricardos_database/Supplementary_Data_1_Databases_Trees/d_AamoA.db_nr_aln.fasta"
qiime_tax="/Users/siljanen/Documents/AA_MiSeq_data_LCG/ricardos_database/Supplementary_Data_1_Databases_Trees/e_AamoA.db_nr_aln_taxonomy_qiime.txt"
mothur_tax="/Users/siljanen/Documents/AA_MiSeq_data_LCG/ricardos_database/Supplementary_Data_1_Databases_Trees/f_AamoA.db_nr_aln_taxonomy_mothur.txt"
chimera_db="/Users/siljanen/Documents/AA_MiSeq_data_LCG/ricardos_database/Supplementary_Data_1_Databases_Trees/j_AamoA_chimera.ref.db_aln.trim.fasta"


### 7) Discard all sequences that share less than 55% identity with any reference sequences


usearch8 -usearch_global 2-uniques_nochim.fasta -db $db_seq -id 0.55 -strand plus \
-uc 4a-uclust_report.txt -matched 4b-uniques_nochim_match.fasta -notmatched 4c-uniques_nochim_nomatch.fasta
# 78.8% matched (dada2, F=200bp R=200bp)
# 49.3%  matched (dada2, F=200bp R=150bp)
# 2.5% matched (dada2, F=150bp R=150bp)
#  58.1% matched (dada2, F only 200bp)
# With Taz =  55.7% matched (dada2, F only 200bp)

### 8) UCHIME chimera filtration (using parameters defined by Alves et al., 2018)

usearch8 -uchime_ref  4b-uniques_nochim_match.fasta -db $chimera_db \
-nonchimeras 5a-uniques_nochim_match_uchimed.fasta -strand plus -mindiv 1.7 -minh 0.1 -uchimeout 5b-uchime_report.txt
# Found 19/128 chimeras (14.8%), 42 not classified (32.8%) (min1) 
# Found 12/327 chimeras (3.7%), 304 not classified (93.0%) (dada2,  F=200bp R=200bp)
# Found 64/145 chimeras (44.1%), 70 not classified (48.3%) (dada2,  F=200bp R=150bp)
#  Found 0/7 chimeras (0.0%), 7 not classified (100.0%) (dada2,  F=150bp R=150bp)
# Found 2/93 chimeras (2.2%), 36 not classified (38.7%) (dada2, F only 200bp)
## With Taz:  Found 2/117 chimeras (1.7%), 46 not classified (39.3%)

### 9) Annotation of the OTUs

# 9.1) With UCLUST implemented in QIIME1

source activate qiime1
assign_taxonomy.py -i 5a-uniques_nochim_match_uchimed.fasta -t $qiime_tax -r $db_seq --similarity 0.8 \
-o 6-uniques_nochim_match_uchimed_uclust_annotation/
source deactivate qiime1

# Number of uniques annotations
less 6-uniques_nochim_match_uchimed_uclust_annotation/5a-uniques_nochim_match_uchimed_tax_assignments.txt | cut -f2 | sort -u | wc -l
#  6 (dada2,  F=200bp R=200bp)
#  6 (dada2,  F=200bp R=150bp)
#  8 (dada2, F only 200bp)
# With Taz:  9  (dada2, F only 200bp)

# Number of Unassigned OTUs
grep -c Unassigned 6-uniques_nochim_match_uchimed_uclust_annotation/5a-uniques_nochim_match_uchimed_tax_assignments.txt
#  11 (dada2,  F=200bp R=200bp)
#  10 (dada2,  F=200bp R=150bp)
# 0 (dada2, F only 200bp)
# With Taz:  0  (dada2, F only 200bp)

### 10) Check the chimeric nature of the unassigned OTUs.

#make amoA databases free of gaps in alignment:
less $chimera_db | sed 's/-//g' > temp2-amoa_chimera_db_nogaps.fasta
less $db_seq | sed 's/-//g' > temp-amoa_ref_db_nogaps.fasta


# Make blast database 
makeblastdb -in temp2-amoa_chimera_db_nogaps.fasta -out temp5-amoa_chimera_db -dbtype nucl
makeblastdb -in temp-amoa_ref_db_nogaps.fasta -out temp5-amoa_db -dbtype nucl


# Make the fast flat (for easier selection of OTUs) 
less 5a-uniques_nochim_match_uchimed.fasta | \
awk -v RS='>'     -v FS="\n"     -v OFS=""     -v ORS="" '{ if (NR > 1) { printf ">%s\n",$1; $1=""; printf "%s\n",$0 } }' \
> 5c-uniques_nochim_match_uchimed_flat.fasta

# List of unassigned OTUs from UCLUST annotation
grep Unassigned 6-uniques_nochim_match_uchimed_uclust_annotation/5a-uniques_nochim_match_uchimed_tax_assignments.txt | cut -f1 > temp6a-unassigned_uclust_dada2_otu_list.txt


# Extract the unassigned sequences from the flat fast file for UCLUST annotations
grep -A1 -f temp6a-unassigned_uclust_dada2_otu_list.txt 5c-uniques_nochim_match_uchimed_flat.fasta \
| grep -v "^--" > 7-unassigned_uclust_dada2_nochim_amoa_otus.fasta

# Run blast (UCLUST)

blastn -query 7-unassigned_uclust_dada2_nochim_amoa_otus.fasta -db temp5-amoa_chimera_db -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq" -out 8a-blast_outfmt6_report_unassigned.txt
echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqseq" >> 8a-blast_outfmt6_report_unassigned.txt

blastn -query 7-unassigned_uclust_dada2_nochim_amoa_otus.fasta -db temp5-amoa_chimera_db -outfmt 0 -out 8b-blast_report_unassigned_uclust_dada2_otus.txt

## Outcome all unassigned have got blast result only for another end of the sequence. So, these unassigned ones we can exclude. 


blastn -query 5c-uniques_nochim_match_uchimed_flat.fasta  -db temp5-amoa_chimera_db -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -out 8c-blast_outfmt6_report_assigned.txt
echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" >> 8c-blast_outfmt6_report_assigned.txt

blastn -query 5c-uniques_nochim_match_uchimed_flat.fasta -db temp5-amoa_chimera_db -outfmt 0 -out 8d-blast_report_assigned_uclust_dada2_otus.txt

## Outcome all assigned have got blast result only for R1 sequence. Maybe E-value for R2 is lower because the length of the seq is smaller for R2. 


# Carefully review the BLAST alignements and determine if the OTUs sequences are artefacts or not.
# If some are, manually remove them from the file '5c-uniques_nochim_match_uchimed_flat.fasta'
# and save it as '9-clean_asvs_dada2.fasta'

"""

# Import annotation file
annot <- read.table("/Users/siljanen/Documents/AA_MiSeq_data_LCG/LGC_G20002861_part1/PrimerClipped/Henri_AOA_amoA/Peat_soil_DNA/Raw_R1R2/2-dada2/
                    6-uniques_nochim_match_uchimed_uclust_annotation/5a-uniques_nochim_match_uchimed_tax_assignments_NVcleaned.txt",
                    header = FALSE)
annot$V4 <- NULL
names(annot) <- c("asv", "annotation", "confidence")

# Retain the clade annotation in 1 column
library(stringr)
annot_l2 <- str_split_fixed(annot$annotation, ";", 11)[,2]
annot_l1 <- str_split_fixed(annot$annotation, ";", 11)[,1]
annot_l2[annot_l2 == ""] <- annot_l1[annot_l2 == ""]

annot$annotation_l2 <- annot_l2

head(annot)
# Import fasta file (uniques and chemira-free reads) (here we have dada2 produced non-chemric uniques which have been uchimed with usearch8)
library("Biostrings")
citation("Biostrings")

fasta_file <- readDNAStringSet("/Users/siljanen/Documents/AA_MiSeq_data_LCG/LGC_G20002861_part1/PrimerClipped/Henri_AOA_amoA/Peat_soil_DNA/Raw_R1R2/2-dada2/9-clean_asvs_dada2.fasta")
head(fasta_file)
asv = names(fasta_file)
seq = paste(fasta_file)
fasta_df <- data.frame(asv, seq)

head(fasta_file)
head(asv)
head(seq)
head(fasta_df)

# Import ASV table (dada2 generated asv-table)
asv_tab <- read.table("/Users/siljanen/Documents/AA_MiSeq_data_LCG/LGC_G20002861_part1/PrimerClipped/Henri_AOA_amoA/Peat_soil_DNA/Raw_R1R2/2-dada2/3-asv_table.txt", header = T, row.names = NULL)
names(asv_tab)[1] <- "seq"

colnames(asv_tab)
colnames(afa_aggr_rel)

# Merge all three data frames
af <- merge(fasta_df, annot)
afa_abs <- merge(af, asv_tab)

# Compute relative abundances of ASVs
asv_counts <- colSums(afa_abs[,6:ncol(afa_abs)])
afa_rel <- afa_abs
afa_rel[,6:ncol(afa_abs)] <- sweep(afa_abs[,6:ncol(afa_abs)], 2, asv_counts, `/`)
afa_rel[is.na(afa_rel)] <- 0


# Aggregate all identical phylogenetic annotations (full annotations)
afa_aggr <- aggregate(afa_abs[,6:ncol(afa_abs)], by=list(annotation=afa_abs$annotation), FUN=sum)
afa_aggr_rel <- aggregate(afa_rel[,6:ncol(afa_rel)], by=list(annotation=afa_rel$annotation), FUN=sum)

dim(afa_aggr_rel)
# 7  13
# With Taz:  8 19


# Aggregate all identical phylogenetic annotations (level-2 annotations)
afa_aggr_l2 <- aggregate(afa_abs[,6:ncol(afa_abs)], by=list(annotation=afa_abs$annotation_l2), FUN=sum)
afa_aggr_l2_rel <- aggregate(afa_rel[,6:ncol(afa_rel)], by=list(annotation=afa_rel$annotation_l2), FUN=sum)

head(afa_aggr_l2)
head(afa_aggr_l2_rel)

# export the annotation asv tables
"""
mkdir 3-annotation
"""

write.table(afa_abs, "2-dada2/3-annotation/2a-asv_tab_tax.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(afa_rel, "2-dada2/3-annotation/2b-asv_tab_rel_tax.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(afa_aggr, "2-dada2/3-annotation/3a-asv_tab_tax_merged.txt", sep="\t", quote = FALSE, row.names = FALSE)
write.table(afa_aggr_rel, "2-dada2/3-annotation/3b-asv_tab_rel_tax_merged.txt", sep="\t", quote = FALSE, row.names = FALSE)
write.table(afa_aggr_l2, "2-dada2/3-annotation/4a-asv_tab_tax_merged_l2.txt", sep="\t", quote = FALSE, row.names = FALSE)
write.table(afa_aggr_l2_rel, "2-dada2/3-annotation/4b-asv_tab_rel_tax_merged_l2.txt", sep="\t", quote = FALSE, row.names = FALSE)




## 12-top OTUS BARPLOTS for unaggregated data:



## NMDS for simplified asv list
library(vegan)
asv_table_simple <- read.csv("2-dada2/3-annotation/2b-asv_tab_rel_tax_manual_heatmap_simple3pros_for_NMDS.csv",header=T,sep=";",row.names=1)
asv_table_simple_bare <- read.csv("2-dada2/3-annotation/2b-asv_tab_rel_tax_manual_heatmap_simple3pros_onlyBARE_for_NMDS.csv",header=T,sep=";",row.names=1)
asv_table_NMDS <- read.table("2b-asv_tab_rel_tax_forNMDS_noseq copy.txt",header=T,sep="\t",row.names=1)

t_asv_table_simple_bare <- t(asv_table_simple_bare)
t_asv_table_NMDS <- t(asv_table_NMDS)

asv_table_NMDS_withTaz <- read.table("2b-asv_tab_rel_tax_NMDS_withTaz.txt",header=T,sep="\t",row.names=1)

t_asv_table_NMDS_withTaz <- t(asv_table_NMDS_withTaz)


head(t_asv_table_NMDS_withTaz)

test_NMDS <- metaMDS(t_asv_table_NMDS_withTaz)
stressplot(test_NMDS)
ordiplot(test_NMDS,type="n")
orditorp(test_NMDS,display="species",cex=0.3,col="red",air=0.01)
orditorp(test_NMDS,display="sites",cex=0.3,air=0.01)



replicates_list_NEW <- c("peat_002", "peat_005", "peat_008", "peat_017",  
                          "peat_030", "peat_031", "peat_032", "peat_033",
                         "peat_038", "peat_039", "peat_064", 
                         "peat_065", "peat_066", "peat_063", "peat_067", "peat_069", "peat_071", "peat_077")


#surface only:
replicates_groups_NEW <- c("Kev_BS", "Kev_BS", "Kev_BS", "Kev_VS",
                           "Taz_PP_VS", "Taz_PP_VS", "Taz_PP_BS", "Taz_PP_BS", "Taz_PB_BS", "Taz_PB_BS",
                           "Sei_BS", "Sei_BS", "Sei_BS", "Sei_VS",
                           "Tay_BS", "Tay_BS", "Tay_BS", 
                           "Tay_VS")


ps<-c(rep(0,3), rep(0,1), rep(1,2), rep(1,2), rep(2,2), rep(5,3), rep(5,1), rep(6,3), rep(6,1))
cl<-c(rep("black",3), rep("white",1), rep("white",2), rep("black",2), rep("black",2), rep("black",3), rep("white",1), rep("black",3), rep("white",1))

ps_b<-c(rep(21,3),  rep(24,3), rep(25,3))
cl_b<-c(rep(3,3), rep(2,3), rep(6,3))

ordiplot(test_NMDS,type="n")

?pch
# Make points:
points(test_NMDS, pch=ps, col="black", bg=cl, cex=1)
points(test_NMDS_b, pch=ps_b, col="black", bg=cl_b, cex=1.2)
orditorp(test_NMDS_b,display="species",cex=0.5,col="red",air=0.01)

legend(x="topleft", legend=c("Kevo", "Tazorvsky PP","Tazorvsky PB" "Seida", "Taymyr"),
       pch=c(0,1,2, 5, 6), col="black", bty="n", title="Sites", cex=0.6)

legend(x="topright", legend=c("Kevo bare", "Kevo veg","Seida veg", "Seida bare", "Taymyr bare",  "Taymyr veg"),
       pch=c(21,22,23,24,25,3), col="black", pt.bg=c(3,4,5,2,6,7), bty="n", title="Sites", cex=0.6)

#########################################################################################
##  Bubbleplot script modification:

# written by alexandre.bagnoud@gmail.com, April 2018

setwd("/Users/siljanen/Documents/AA_MiSeq_data_LCG/LGC_G20002861_part1/PrimerClipped/Henri_AOA_amoA/Peat_soil_DNA/Raw_R1R2/2-dada2")

library(stringr)
citation("stringr")


# files to import

annot_file <- "/Users/siljanen/Documents/AA_MiSeq_data_LCG/LGC_G20002861_part1/PrimerClipped/Henri_AOA_amoA/Peat_soil_DNA/Raw_R1R2/2-dada2/6-uniques_nochim_match_uchimed_uclust_annotation/5a-uniques_nochim_match_uchimed_tax_assignments_NVcleaned_DADA+USEARCH.txt"

annot_type <- "uclust" # either "uclust" or "mothur"
otu_tab_file <- "/Users/siljanen/Documents/AA_MiSeq_data_LCG/LGC_G20002861_part1/PrimerClipped/Henri_AOA_amoA/Peat_soil_DNA/Raw_R1R2/2-dada2/3-annotation/2a-asv_tab_tax_mod_for_bubble.txt"

## Other variables to define
output_folder <- "10_plots_and_tables"
plot_title <- "dada2+usearch_R1_pipeline"

# Number of taxons to display in the full-length taxonomy bar plots (the low abundance one will be pool together in a bin nameed "Other")
tax_number <- 6

#Surface only
replicates_list <- c("peat_2", "peat_5", "peat_8", "peat_17",  "peat_64", 
                    "peat_65", "peat_66", "peat_63", "peat_67", "peat_69", "peat_71", "peat_77")

#surface only:
replicates_groups <- c("Kev_BS", "Kev_BS", "Kev_BS", "Kev_VS",
                        "Sei_BS", "Sei_BS", "Sei_BS", "Sei_VS",
                       "Tay_BS", "Tay_BS", "Tay_BS", 
                       "Tay_VS")

dir.create(output_folder)

# Import annotation file
if (annot_type == "uclust") {
  annot <- read.table(annot_file)
  annot$V3 <- NULL
  annot$V4 <- NULL
  names(annot) <- c("otu", "annotation")
  annot$annot_l2 <- str_split_fixed(annot$annotation, ";", 11)[,2]
  annot$annotation_short <- gsub("^.*;", "", annot$annotation)
  annot$otu <- sub(";.*;", "", annot$otu)
} else if (annot_type == "mothur") {
  annot_temp <- read.table(annot_file)
  otu.vec <- gsub(";.*;", "", annot_temp$V1)
  annot.vec <- gsub("\\([0-9]*\\)", "", sub("_unclassified.*$","" , annot_temp$V2))
  annot_l2.vec <- str_split_fixed(annot.vec, ";", 11)[,2]
  annot_short.vec <- sub("^.*;", "", annot.vec)
  annot <- data.frame(otu = otu.vec, annotation = annot.vec, annot_l2=annot_l2.vec, annot_short=annot_short.vec)
}


# Import ASV table
otu_tab <- read.table(otu_tab_file, header = T, row.names = NULL, sep = "\t", comment.char = "@")
names(otu_tab)[1] <- "otu"

# Compute relative abundances of ASVs
otu_counts <- colSums(otu_tab[,-1])
otu_tab_rel <- otu_tab
otu_tab_rel[,-1] <- sweep(otu_tab[,-1], 2, otu_counts, `/`)
otu_tab_rel[is.na(otu_tab_rel)] <- 0

# Merge annotation and otu table
otutab_tax_abs <- merge(annot, otu_tab)
otutab_tax_rel <- merge(annot, otu_tab_rel)

# Aggregate all identical phylogenetic annotations (full annotations)
otutab_tax_aggr <- aggregate(otutab_tax_abs[,5:ncol(otutab_tax_abs)], by=list(annotation=otutab_tax_abs$annotation), FUN=sum)
otutab_tax_aggr_rel <- aggregate(otutab_tax_rel[,5:ncol(otutab_tax_rel)], by=list(annotation=otutab_tax_rel$annotation), FUN=sum)

# Aggregate all identical phylogenetic annotations (level-2 annotations)
otutab_tax_aggr_l2 <- aggregate(otutab_tax_abs[,5:ncol(otutab_tax_abs)], by=list(annotation=otutab_tax_abs$annot_l2), FUN=sum)
otutab_tax_aggr_l2_rel <- aggregate(otutab_tax_rel[,5:ncol(otutab_tax_rel)], by=list(annotation=otutab_tax_rel$annot_l2), FUN=sum)

paste0(output_folder, "/", "1a-otu_tab_tax.txt")

# export the annotation asv tables
write.table(otutab_tax_abs, paste0(output_folder, "/", "1a-otu_tab_tax.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(otutab_tax_rel, paste0(output_folder, "/", "1b-otu_tab_rel_tax.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(otutab_tax_aggr, paste0(output_folder, "/", "2a-otu_tab_tax_merged.txt"), sep="\t", quote = FALSE, row.names = FALSE)
write.table(otutab_tax_aggr_rel, paste0(output_folder, "/", "2b-otu_tab_rel_tax_merged.txt"), sep="\t", quote = FALSE, row.names = FALSE)
write.table(otutab_tax_aggr_l2, paste0(output_folder, "/", "3a-otu_tab_tax_merged_l2.txt"), sep="\t", quote = FALSE, row.names = FALSE)
write.table(otutab_tax_aggr_l2_rel, paste0(output_folder, "/", "3b-otu_tab_rel_tax_merged_l2.txt"), sep="\t", quote = FALSE, row.names = FALSE)



############################################################
  
library("RColorBrewer")
library("reshape2")
library("ggplot2")

citation("RColorBrewer")
citation("reshape2")
citation("ggplot2")
citation("vegan")

## This was made for previsou version of the data set without Taz:
dat_3pros_lim <- "/Users/siljanen/Documents/AA_MiSeq_data_LCG/LGC_G20002861_part1/PrimerClipped/Henri_AOA_amoA/Peat_soil_DNA/Raw_R1R2/2-dada2/3-annotation/data_3pros_lim_mod_molten.txt"

setwd("/Users/siljanen/Documents/AA_MiSeq_data_LCG/LGC_G20002861_part1/PrimerClipped/Henri_AOA_amoA/Peat_soil_DNA/Raw_R1R2/2-dada2/3-annotation/")

otu_3pros_lim_molten <- read.table("data_3pros_lim_mod_molten.txt" , header = T, row.names = NULL, sep = "\t", comment.char = "@")
otu_3pros_lim_molten$X <- NULL

otutab_Top40asvs_molten <- read.table("Final_otutab_Top40asvs.txt" , header = T, row.names = NULL, sep = "\t", comment.char = "@")

## This was made for previous version of the data set without Taz:

setwd("/Users/siljanen/Documents/AA_MiSeq_data_LCG/LGC_G20002861_part1/PrimerClipped/Henri_AOA_amoA/Peat_soil_DNA/Raw_R1R2/2-dada2/3-annotation/")

#Taxonomy merged to OTU-level 
# Files re-organized in Excel: mean of replicates calculated and dataframe molten manually,
##  In Text-editor Varible names are shorten: e.g. NS;NS-Gamma;NS-Gamma-2;NS-Gamma-2,3;NS-Gamma-2,3,2;NS-Gamma-2,3,2_OTU5 
####  to NS-Gamma-2,3,2_OTU5

otutab_with_Taz_molten <- read.table("3b-asv_tab_rel_tax_merged_molten.txt" , header = T, row.names = NULL, sep = "\t", comment.char = "@")
otutab_Top40asvs_with_Taz_molten <- read.table("With_Taz_otutab_top40_molten_final.txt" , header = T, row.names = NULL, sep = "\t", comment.char = "@")

otutab_Top40asvs_with_Taz_molten_B_V_order <- read.table("With_Taz_otutab_top40_molten_final_Bare_Veg_order.txt" , header = T, row.names = NULL, sep = "\t", comment.char = "@")

otutab_g_merged_with_Taz_molten <- read.table("4b-asv_rel_g-merged_molten.txt" , header = T, row.names = NULL, sep = "\t", comment.char = "@")

otutab_g_merged_with_Taz_molten_B_V_order <- read.table("4b-asv_rel_g-merged_molten_Bare_Veg_order.txt" , header = T, row.names = NULL, sep = "\t", comment.char = "@")

plot_title <- "DADA2+usearch pipeline"
bubble_plot_file_name <- "5-barplots/3-bubbleplot_dada2nochim.fasta"


bubble_plot_Top40asvs_with_Taz_molten_B_V_order <- ggplot(otutab_Top40asvs_with_Taz_molten_B_V_order,aes(sample,variable)) +
  #geom_point(aes(size=value+sd),shape=16, color = "red") +
  geom_point(aes(size=value, fill=clade),shape=21,color="black") +
  #geom_point(aes(size=value, color=otu_molten2$clade)) +
  theme(panel.grid.major=element_line(linetype=1,color="grey"),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0),
        panel.background = element_blank()) +
  ylab("Thaumarchaeota OTUs") +
  xlab("Samples") +
  scale_fill_brewer(palette="Set1", name="Taxonomic\nclade") +
  scale_size(name = "Relative\nabundance") +
  ggtitle(plot_title)

bubble_plot_Top40asvs_with_Taz_molten_B_V_order

svg(paste0(output_folder, "/7-bubbleplot.svg"), width = 6, height = 6)
bubble_plot
dev.off()

############################################################

# export the annotation asv tables
write.table(otutab_tax_aggr, paste0(output_folder, "/", "2a-otu_tab_tax_merged.txt"), sep="\t", quote = FALSE, row.names = FALSE)
write.table(otutab_tax_aggr_rel, paste0(output_folder, "/", "2b-otu_tab_rel_tax_merged.txt"), sep="\t", quote = FALSE, row.names = FALSE)

## In 2a-asv_tab_tax_30top.xlsx file relative abundance was recalculated and TOP40 asvs were selected.
##  From these new otutab_tax_aggr_rel object was created to make Bubbleplot. 


dir()
otutab_tax_aggr_rel <- read.table("2a-asv_tab_tax_30top.txt" , header = T, row.names = NULL, sep = "\t", comment.char = "@")

#Surface only
replicates_list <- c("Peat_2", "Peat_5", "Peat_8", "Peat_17", "Peat_64", 
                     "Peat_65", "Peat_66", "Peat_63", "Peat_67", "Peat_69", "Peat_71", "Peat_77")

#surface only:
replicates_groups <- c("Kev_BS", "Kev_BS", "Kev_BS", "Kev_VS",
                       "Sei_BS", "Sei_BS", "Sei_BS", "Sei_VS",
                       "Tay_BS", "Tay_BS", "Tay_BS", 
                       "Tay_VS")


# Transpose the ASV table
n <- otutab_tax_aggr_rel$annotation2
t.otu_aggr_rel <- as.data.frame(t(otutab_tax_aggr_rel[,c(-1,-nrow(otutab_tax_aggr_rel))]))
colnames(t.otu_aggr_rel) <- n
t.otu_aggr_rel$sample <- rownames(t.otu_aggr_rel)
rownames(t.otu_aggr_rel) <- NULL
t.otu_aggr_rel <- t.otu_aggr_rel[1:ncol(otu_tab)-1,]

# Calculate average of each replicate
t.otu_aggr_rel$replicate <- rep(NA, nrow(t.otu_aggr_rel))
for (line in 1:(nrow(t.otu_aggr_rel)-1)){
  t.otu_aggr_rel$replicate[line] <- replicates_groups[t.otu_aggr_rel$sample[line] == replicates_list]
}

t.otu_replicates <- aggregate(t.otu_aggr_rel[,1:(ncol(t.otu_aggr_rel)-2)],
                              by = list(t.otu_aggr_rel$replicate),
                              FUN = "mean")

names(t.otu_replicates)[1] <- "sample"


# Melt the dataframe
otu_replicates_molten <- melt(t.otu_replicates, id.vars = "sample")
otu_replicates_molten$clade <- str_split_fixed(otu_replicates_molten$variable, ";", 11)[,2]
otu_replicates_molten$clade[otu_replicates_molten$clade == ""] <- "Unassigned"
otu_replicates_molten$variable <- gsub("^.*;", "", otu_replicates_molten$variable)
otu_replicates_molten$sample <- factor(otu_replicates_molten$sample, levels = unique(replicates_groups))
otu_replicates_molten$row.id <- paste0(otu_replicates_molten$sample, otu_replicates_molten$variable)

# Calculate the standard deviation for each replicate
t.otu_sd <- aggregate(t.otu_aggr_rel[,1:(ncol(t.otu_aggr_rel)-2)],
                      by = list(t.otu_aggr_rel$replicate),
                      FUN = "sd")

names(t.otu_sd)[1] <- "sample"

# Melt the sd dataframe
otu_sd_molten <- melt(t.otu_sd, id.vars = "sample")
otu_sd_molten$variable <- gsub("^.*;", "", otu_sd_molten$variable)
names(otu_sd_molten)[3] <- "sd"
otu_sd_molten$row.id <- paste0(otu_sd_molten$sample, otu_sd_molten$variable)
otu_sd_molten$sample <- NULL
otu_sd_molten$variable <- NULL

# Merge both dataframe (mean and sd)
otu_molten <- merge(otu_replicates_molten, otu_sd_molten)

# Reorder the taxonomic annotation for the plot
tax_levels <- sort(unique(list(sub("^.*;", "", otutab_tax_aggr_rel$annotation))[[1]]), decreasing = TRUE)
otu_molten$variable <- factor(otu_molten$variable, levels = tax_levels)
otu_molten$variable <- factor(otu_molten$variable)

# Remove null values
otu_molten2 <- otu_molten[otu_molten$value > 0,]

View(otu_molten2)

bubble_plot <- ggplot(otu_molten2,aes(sample,variable)) +
  #geom_point(aes(size=value+sd),shape=16, color = "red") +
  geom_point(aes(size=value, fill=otu_molten2$clade),shape=21,color="black") +
  #geom_point(aes(size=value, color=otu_molten2$clade)) +
  theme(panel.grid.major=element_line(linetype=1,color="grey"),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0),
        panel.background = element_blank()) +
  ylab("Thaumarchaeota OTUs") +
  xlab("Samples") +
  scale_fill_brewer(palette="Set1", name="Taxonomic\nclade") +
  scale_size(name = "Relative\nabundance") +
  ggtitle(plot_title)

bubble_plot

svg(paste0(output_folder, "/7-bubbleplot.svg"), width = 6, height = 6)
bubble_plot
dev.off()

