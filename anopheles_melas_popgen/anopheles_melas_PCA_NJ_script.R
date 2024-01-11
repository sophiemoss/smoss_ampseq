
## Shift + Options + R to open a new terminal
## Ctrl + Option + Enter to send commands to terminal 
## ssh sophie@s11
## conda activate renv

## attach all required R packages
library(dplyr)
library(scales)
library(data.table)
library(amap)
library(RColorBrewer)
library(colorRamps)
library(scales)
library(unikn)
library(ggplot2)
library(showtext)
showtext_auto()

##setwd
setwd('/mnt/storage11/sophie/bijagos_mosq_wgs/2019_melas_fq2vcf/genomics_database_melas2019plusglobal/genomics_database_melas2019plusglobal_vcf/melas2019plusglobal_vcf_filtering')
getwd()

## load binary snp matrix (biallelic nucleotide matrix for making the PCA plots, ending mat.bin)(this binary matrix was made from the filtered VCF file)
snp<-fread("mt_melas_plus_global_subset_filtered.mat.bin",sep="\t",header=T)
snp<-as.data.frame(snp)
count(snp)
nrow(snp)
ncol(snp)

## read in metadata

met <- read.table("metadata_melasplusglobal.csv", sep=",", stringsAsFactors = FALSE, header = TRUE)
count(met)

## check that metadata has loaded correctly
# ensure the order of the samples in the metadata matches that of the SNP matrix
# then selects relevant columns from the metadata for downstream analysis

head(met,15)
met<-met[match(colnames(snp)[-(1:3)], met$sample),]
count(met)

## allow metadata fields for plotting (aka plot by country, region, site, sample, etc.)

met <- met %>% select(c(sample,year,country,species,island))

# inspect metadata
colnames(met)
table(met$"country")
unique(met[c("island")])
unique(met[c("species")])
count(met)

# Select which countries to include in analysis

met_subset<-met[which(met$country %in% c("Guinea-Bissau")),]
met_sample<-met_subset$sample
snp_subset <- snp[,which(colnames(snp) %in% c("chr", "pos", "ref", met_sample))]

count(snp_subset)

#Inspect final metadata

colnames(met_subset)
table(met_subset$"country")
unique(met_subset[c("country")])
unique(met_subset[c("species")])
count(met_subset)

##Remove outliers: bu1003_Combined, bu1034_Combined, so1023_Combined

met_ro<-met_subset[which(met_subset$sample=="bu1003_Combined"),]
met_ro<-met_subset$sample[which(met_subset$sample=="bu1003_Combined")]
met_ro<-as.vector(met_ro)
snp_subset<-snp_subset[,-which(colnames(snp_subset)%in%met_ro)]
met_subset<-met_subset[-which(met_subset$sample=="bu1003_Combined"),]
count(met_subset)
# 
# met_ro<-met_subset[which(met_subset$sample=="bu1034_Combined"),]
# met_ro<-met_subset$sample[which(met_subset$sample=="bu1034_Combined")]
# met_ro<-as.vector(met_ro)
# snp_subset<-snp_subset[,-which(colnames(snp_subset)%in%met_ro)]
# met_subset<-met_subset[-which(met_subset$sample=="bu1034_Combined"),]
# count(met_subset)
# 
# met_ro<-met_subset[which(met_subset$sample=="so1023_Combined"),]
# met_ro<-met_subset$sample[which(met_subset$sample=="so1023_Combined")]
# met_ro<-as.vector(met_ro)
# snp_subset<-snp_subset[,-which(colnames(snp_subset)%in%met_ro)]
# met_subset<-met_subset[-which(met_subset$sample=="so1023_Combined"),]
# count(met_subset)
# 
## Choose plot type (do this before generating the distance matrix)

#run by island
met_subset$island<-factor(met_subset$island)

#run by country
met_subset$country<-factor(met_subset$country)

# calculate distance matrix within this script if small dataset
# filter snp to contain only the samples in met_subset
# snp_c <- snp[, c(1:3, which(colnames(snp[-c(1:3)]) %in% met_subset$sample))]

# calculate the distance matrix used to make the plots
# dist<-Dist(t(snp_c), method = "manhattan")
# summary(dist)

# Or you can create distance matrix using plink on the server instead (not in R).
# Then load it in and combine with slimmed metadata.

dist <- read.table("mito_only_melas_plusglobal_renamed.dist", header = F)
id <- read.table("mito_only_melas_plusglobal_renamed.dist.id", header = F)

ncol(dist)
ncol(id)

names(id)
names(met)

# need to subset id variable so that it only has the sample_ids that are in met_subset

samplestoinclude <- which(id[,1] %in% met_subset$sample) #takes sample IDs from met_subset and matches with id variable
length(samplestoinclude) #prints the variable so you can check it's the right number of samples to be in the distance matrix
dist_filt <- dist[samplestoinclude,samplestoinclude] #filters the distance matrix to only include the coordinates which contain distances for these samples
dim(dist_filt) #check that it contains the same numbers of rows and columns, and that the number of distances = the number of samples

id_filt <-id[samplestoinclude,] #filters the id variable to only include the sample id's we care about
dim(id_filt) #check that it's the right number of samples
data <- id_filt %>% left_join(met_subset, by = c("V1" = "sample")) #data is produced by joining up the id_filt variable with the rest of the metadata information for the samples
head(data) #check and see that the metadata has joined

dist_m <- as.matrix(dist_filt)
colnames(dist_m) <- data$V1
rownames(dist_m) <- data$V1
dim(dist_m)

# can save dist_m to start from next time, as this takes a while to produce
# saveRDS(dist_m, file = "/mnt/storage9/sophie/basebijagospopgen/wgs_global_database/filtered/pca/distancematrix_20222311.rds")

# to read this RDS file in future to start from this point, make sure you're in the right directory and then use this code:
# readRDS(distancematrix.rds, refhook = NULL)

# Perform principal components analysis with cmdscale on the distance matrix. k = 10, so the first 10 principal components are calculated. 

cmd_m<-cmdscale(dist_m, k = 10, eig = TRUE, x.ret = TRUE)

# Calculate the proportion of variance explained by the first two components
# include the first and second values in the axes titles, first value = PC1, second value = PC2, etc...

#vars <- round(cmd_pv$eig / sum(cmd_pv$eig)*100, 1)
variance <- (cmd_m$eig / sum(cmd_m$eig)*100)
var1<-variance[[1]]
var2<-variance[[2]]

island<-met_subset$island
island<-as.character(island)

country<-met_subset$country
country<-as.character(country)

year<-met_subset$year
year<-as.character(year)

sample<-met_subset$sample
sample<-as.character(sample)

# PCA plotting

# Selecting first 6 pairs of colors from pal_unikn_pair
colours <- as.vector(unikn::pal_unikn_pref[1:6])

# Create a ggplot object
p <- ggplot(data, aes(x = cmd_m$points[, 1], y = cmd_m$points[, 2])) +
  geom_point(aes(colour = as.factor(island)), pch = 20, size = 3) +
  scale_color_manual(values = colours) +
  xlab(paste("PC1 ", round(var1, digits = 2), "%")) +
  ylab(paste("PC2 ", round(var2, digits = 2), "%")) +
  theme_minimal()

# Add sample names as labels
#text(
# cmd_m$points[, 1], cmd_m$points[, 2], labels = sample,
#pos = 3, col = "black", cex = 0.7
#)

#legend("bottomright", inset = c(-0.43, 0), col = unique(pal[as.factor(island)]), pch = 20, leg = unique(as.factor(island)))
#dev.off()

# Save the plot as a PDF file
pdf("mito_only_melas_pca_plot.pdf", width = 8, height = 6)  # Specify the filename and dimensions
print(p)  # Print the plot to the PDF file
dev.off()  # Close the PDF device

# generate a tree file that corresponds to the PCA plots, can be visualised in iTOL or in R
# install.packages("ape")
# library(ape)
# tree_m<-nj(dist_m)

# Writing to a newick tree file
# write.tree(phy=tree_m, file="20230809_tree.newick")

# Create PCA using old script (not ggplot)

pal<-c("red", "cyan", "purple", "violet", "green", "blue")

# generate plots 
#PCA by island
pdf("mito_only_melas_smalltext.pdf",width=8)
par(mar=c(5.1, 4.1, 4.1, 11.1), xpd=TRUE)
plot(cmd_m$points[,1], cmd_m$points[,2],col=pal[as.numeric(as.factor(island))],pch=25,cex=1.25,xlab=paste("PC1 (",round(var1,digits=2),")"),ylab=paste("PC2 (",round(var2,digits=2),")"))

text(
 cmd_m$points[, 1], cmd_m$points[, 2], labels = sample,
pos = 3, col = "black", cex = 0.2
)

legend("bottomright",inset=c(-0.3,0), col=unique(pal[as.factor(country)]),pch=25,leg=unique(as.factor(country)))
dev.off()
