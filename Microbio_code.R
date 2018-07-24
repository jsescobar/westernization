# Code used in the statistical analyses of the paper:
# de la Cuesta-Zuluaga J, Corrales-Agudelo V, Velásquez-Mejía EP, Carmona JA, Abad JM, Escobar JS. 2018.
# Gut microbiota is associated with obesity and cardiometabolic disease in a population in the midst of Westernization. 
# Scientific Reports. doi: 10.1038/s41598-018-29687-x.
# www.nature.com/articles/s41598-018-29687-x

# Cleans the workspace
#rm(list = ls())

# Libraries
library(GUniFrac) # UniFracs
library(ade4) # s.class
library(phytools) # Load tree - read.newick
library(cluster) # Clustering
library(clusterSim) # Clustering
library(ggplot2) # Graphs
library(car) # Anova
library(NMF) # Heatmap
library(qvalue) # False discovery rate
library(cowplot) # plot_grid
library(ggdendro) # Dendrogram
library(BiodiversityR) # Alpha diversity
library(reshape2) # melt
#library(remotes) # for intalling github packages (useful for Tax4Fun installation)
library(Tax4Fun) # Metagenome prediction
library(ecodist) # Mantel test
library(plot3D) # 3D graphs

# Loading data ----
# Change directory name to where the data were extracted
setwd(dir = "XXXX")

# Metadata table
microbio.meta = read.table(file = "./files/microbio_selected.meta", header = T, 
                           sep = "\t", dec = ",", row.names = 1)

# OTU table
microbio.otus = read.table(file = "./files/microbio_selected.otus", header = T, 
                           sep = "\t", row.names = 1)

# Phylogenetic tree
microbio.tree = read.newick(file = "./files/microbio_selected.tre")

# Taxonomy
microbio.taxonomy = read.table("./files/microbio_selected.taxonomy", sep = "\t",
                               row.names = 1, header = T)

# OTU rarefaction
# By default, it uses the minimum number. Verify with rowSums(microbio.rare)
microbio.rare = Rarefy(microbio.otus)$otu.tab.rff

# Calculate OTU relative frecuencies
microbio.relative = t(microbio.otus/rowSums(microbio.otus)) 

# In the files, there are five samples that were sequenced twice
replicate_samples = c("MI_008_H2", "MI_093_H12", "MI_130_H2", "MI_198_H2", "MI_458_H2")
replicate_positions = c(9, 95, 132, 201, 445)

# Calculate the difference between replicate samples
# Rarefied counts
dif_MI008 = summary(unlist(microbio.rare[8,] - microbio.rare[9,]))
dif_MI093 = summary(unlist(microbio.rare[94,] - microbio.rare[95,]))
dif_MI130 = summary(unlist(microbio.rare[131,] - microbio.rare[132,]))
dif_MI198 = summary(unlist(microbio.rare[200,] - microbio.rare[201,]))
dif_MI458 = summary(unlist(microbio.rare[444,] - microbio.rare[445,]))

# Relative frequencies
dif_MI008 = summary(unlist(microbio.relative[8,] - microbio.relative[9,]))
dif_MI093 = summary(unlist(microbio.relative[94,] - microbio.relative[95,]))
dif_MI130 = summary(unlist(microbio.relative[131,] - microbio.relative[132,]))
dif_MI198 = summary(unlist(microbio.relative[200,] - microbio.relative[201,]))
dif_MI458 = summary(unlist(microbio.relative[444,] - microbio.relative[445,]))

# Remove replicates for the final dataset
microbio.meta = microbio.meta[-replicate_positions,]
microbio.otus = microbio.otus[-replicate_positions,]
microbio.rare = microbio.rare[-replicate_positions,]
microbio.relative = microbio.relative[,-replicate_positions]

# To create phylotypes for each taxonomic level ----
# Sums all the OTUs with the same taxonomy

# The following code was taken from a tutorial found here:
# http://www.r-bloggers.com/from-otu-table-to-heatmap/
# https://learningomics.wordpress.com/2013/02/23/from-otu-table-to-heatma/
#
# Functions
# Function to separate taxonomies
extract.name.level = function(x, level){
  a=c(unlist(strsplit(x,';')),'Other')
  paste(a[1:min(level,length(a))],collapse=';')
}

# Function to sumarize the data at different taxonomic levels
otu2taxonomy = function(x, level, taxa=NULL){
  if(is.null(taxa)){
    taxa = colnames(x)
  }
  if(length(taxa)!=dim(x)[2]){
    print("ERROR: taxonomy should have the same length as the number of columns in the OTU table")
    return;
  }
  level.names = sapply(as.character(taxa), 
                       function(x)
                         extract.name.level(x,level=level))
  t(apply(x, 1, 
          function(y) 
            tapply(y,level.names,sum)))
}

# Separate the abundance data from the taxonomy
taxa.names = microbio.taxonomy$Taxonomy
dat2 = t(microbio.otus)

# Remove samples with low sequence count
s_abundances = apply(dat2,2,sum)

# Separate the data that are above and below the threshold (1000 reads in this case)
bads = dat2[,s_abundances<1000]
goods = dat2[,s_abundances>1000]

# Number of samples that are above and below the threshold
ncol(goods)
ncol(bads)

# Keep only 'good' samples
dat2 = goods

dat2 = scale(dat2, center=F, scale=colSums(dat2))
dat2 <-t(dat2)

# Separate objects for each taxonomic level
# Greengenes taxonomy has 7 levels, starting from 1 (Kingdom), 2 (Phylum), ... 
# k__Kingdom;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Speciess;
d.phylum = otu2taxonomy(dat2,level=2,taxa=taxa.names)
d.class = otu2taxonomy(dat2,level=3,taxa=taxa.names)
d.order = otu2taxonomy(dat2,level=4,taxa=taxa.names)
d.family = otu2taxonomy(dat2,level=5,taxa=taxa.names)
d.genus = otu2taxonomy(dat2,level=6,taxa=taxa.names)
d.species = otu2taxonomy(dat2,level=7,taxa=taxa.names)

# # Transpose the tables and export the files
# phylum2 <-t(d.phylum)
# class2 <-t(d.class)
# order2 <-t(d.order)
# family2 <-t(d.family)
# genus2 <-t(d.genus)
# species2 <-t(d.species)
# 
# # Create a new folder and save taxonomy (phylotype) tables
# dir.create(path = "./phylotypes/")
# 
# write.table(phylum2, file="phylotypes/phyla.txt", col.names=NA,row.names=TRUE,
# sep="\t", quote=FALSE)
# write.table(class2, file="phylotypes/classes.txt", col.names=NA,row.names=TRUE,
# sep="\t", quote=FALSE)
# write.table(order2, file="phylotypes/orders.txt", col.names=NA,row.names=TRUE,
# sep="\t", quote=FALSE)
# write.table(family2, file="phylotypes/families.txt", col.names=NA,row.names=TRUE,
# sep="\t", quote=FALSE)
# write.table(genus2, file="phylotypes/genera.txt", col.names=NA,row.names=TRUE,
# sep="\t", quote=FALSE)
# write.table(species2, file="phylotypes/species.txt", col.names=NA,row.names=TRUE,
# sep="\t", quote=FALSE)
# 
# # Using absolute OTU frecuencies
# abs_phylum = otu2taxonomy(t(goods),level=2,taxa=taxa.names)
# abs_class = otu2taxonomy(t(goods),level=3,taxa=taxa.names)
# abs_order = otu2taxonomy(t(goods),level=4,taxa=taxa.names)
# abs_family = otu2taxonomy(t(goods),level=5,taxa=taxa.names)
# abs_genus = otu2taxonomy(t(goods),level=6,taxa=taxa.names)
# abs_species = otu2taxonomy(t(goods),level=7,taxa=taxa.names)
 
# Descriptive statistics by phylum and OTU ----

# By phylum
# Melt the phylum table
phylum = t(d.phylum)
phylum_melt = melt(phylum)

# Compute the mean and standard deviation
mean_abund_phylum = aggregate(value ~ Var1, data = phylum_melt, FUN = mean)
sd_abund_phylum = aggregate(value ~ Var1, data = phylum_melt, FUN = sd)

# Ordered mean and SD table
mean_abund_phylum = cbind(mean_abund_phylum, sd = sd_abund_phylum$value)
mean_abund_phylum = mean_abund_phylum[order(mean_abund_phylum[,2], decreasing = T),]


# Complete OTU table
# Empty OTUs (that only appeared in the replicates) should be removed
otu_summary = data.frame(mean.abundance = round(rowMeans(microbio.relative[rowSums(microbio.relative) > 0,])*100, 2),
                         Taxonomy = microbio.taxonomy[rowSums(microbio.relative) > 0, 2])

# By OTUs
# Melt the OTU table
otus_melt = melt(microbio.relative)

# Compute mean and standard deviation
mean_abund_otus = aggregate(value ~ Var1, data = otus_melt, FUN = mean)
sd_abund_otus = aggregate(value ~ Var1, data = otus_melt, FUN = sd)

# Ordered mean and SD table
mean_abund_otus = cbind(mean_abund_otus, sd = sd_abund_otus$value)
mean_abund_otus = mean_abund_otus[order(mean_abund_otus[,2], decreasing = T),]
top_ten_otus = head(mean_abund_otus, 10)

# Boxplot of phyla and top OTUs
# Boxplot of phyla

# Combine phyla with very low abundance
phyla_median = aggregate(value ~ Var1, data = phylum_melt, FUN = median)
top_phyla = phyla_median$Var1[phyla_median$value > 0]
bottom_phyla = phyla_median$Var1[phyla_median$value == 0]

top_bottom_phyla = rbind(phylum[top_phyla, ], "Other" = colSums(phylum[bottom_phyla, ]))

phylum_melt = melt(top_bottom_phyla)


#phylum_melt$value[phylum_melt$value < 0.00005] = 0.00005
phyla_labels = c("Firmicutes", "Bacteroidetes", "Actinobacteria",
                 "Proteobacteria", "Verrucomicrobia", "Tenericutes", 
                 "Cyanobacteria", "Euryarchaeota", "Fusobacteria", 
                 "Other phyla")

box_phylum = ggplot(phylum_melt, aes(x = reorder(Var1, -value, FUN=median), y = value)) +
  geom_boxplot() + 
  labs(x = "", y = "Relative abundance") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = phyla_labels)

# Boxplot of OTUs
top_otus_melt = melt(microbio.relative[top_ten_otus$Var1,])
OTU_labels = c("Akkermansia muciniphila", "Prevotella copri", 
               "Escherichia coli", "Faecalibacterium prausnitzii", 
               "Bifidobacterium adolescentis", "Enterobacter hormaechei", 
               "Gemmiger formicilis", "Ruminococcus bromii",
               "Methanobrevibacter", "Oscillospira")

box_otus = ggplot(top_otus_melt, aes(x = reorder(Var1, -value, FUN=median), y = value)) +
  geom_boxplot() + 
  labs(x = "", y = "Relative abundance") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        axis.text.y = element_text(size = 8), 
        axis.title=element_text(size=9)) +
  scale_x_discrete(labels = OTU_labels)

# Figure 1
ggdraw() +
  draw_plot(box_phylum, 0, 0, 1, 1, scale = 1) +
  draw_plot(box_otus,  0.45, 0.35, 0.60, 0.73, scale = 0.8) +
  draw_plot_label(c("A", "B"), c(0, 0.50), c(1, 1), size = 15)


# Obtention of UniFrac distance matrices ----
 
unifracs <- GUniFrac(microbio.rare, microbio.tree, alpha=c(0, 0.5, 1))$unifracs

dw <- unifracs[, , "d_1"]   # Weighted UniFrac
du <- unifracs[, , "d_UW"]  # Unweighted UniFrac    
#d5 <- unifracs[, , "d_0.5"] # Generalized UniFrac with alpha 0.5


# Enterotype analysis ----
# Clustering analysis and evaluation of the enterotyping
# Modified from Arumugam et al., 2011 (doi:10.1038/nature09944)
# Code available at enterotyping.embl.de

# K-medoids function
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

# Weighted UniFrac for use with enterotypes
data.dist = as.dist(dw)

# PCoA analysis and computation of the proportion of the variance of each axis
e.pcoa = cmdscale(data.dist, k=5, eig = T)
e.PC1 = round(e.pcoa$eig[1]/sum(e.pcoa$eig), 4)* 100
e.PC2 = round(e.pcoa$eig[2]/sum(e.pcoa$eig), 4)* 100
e.PC3 = round(e.pcoa$eig[3]/sum(e.pcoa$eig), 4)* 100

# Optimum number of clusters using the Silhouette index instead of Calinski-Harabasz index
nclusters=NULL

for (k in 1:20) {
  if (k==1) {
    nclusters[k]=NA
  } else {
    data.cluster_temp=pam.clustering(data.dist, k)
    nclusters[k] = mean(silhouette(data.cluster_temp, data.dist)[,3])
    #nclusters[k]=index.G1(t(data),data.cluster_temp,  d = data.dist, centrotypes = "medoids")
  }
}

# Grouping with k clusters
data.cluster=pam.clustering(data.dist, k=2)

# SI index for the above clustering
# -1 <= S(i) <= 1
# A sample which is closer to its own cluster than to any other cluster has a high S(i),
# S(i) close to 0 implies that the given sample lies somewhere between two clusters.
# Large negative S(i) values indicate that the sample was assigned to the wrong cluster.
obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])
cat(obs.silhouette)

# plots SI index ~ number of clusters
ggplot(data.frame(n = factor(2:20), nclusters[-1]), aes(n, nclusters[-1])) + 
  geom_bar(stat = "identity") + labs(x = "Number of clusters", y = "Average silhouette index") + 
  theme_gray()

# PCoA graph of the enterotyping with k clusters
clusters = as.factor(data.cluster)
w = as.data.frame(cmdscale(dw, k = 3))
ggplot(w, aes(V1, V2)) + geom_point(aes(color = clusters)) + 
  stat_ellipse(type = "t", mapping = aes(color = clusters), level = 0.75) + 
  labs(x = paste("PC1",e.PC1, "%"),
       y = paste("PC2",e.PC2, "%"), 
       title = "PCoA Weighted UniFrac")

# Enterogradient analysis ----
# Selection of OTUs with median relative abundance >= 0.01%
microbio.relative_tax = cbind(microbio.taxonomy[,2], microbio.relative)

# Compute the median of each OTU and filter
median_abund = apply(microbio.relative, MARGIN = 1, FUN = median)
abundant_otus = microbio.relative[median_abund >= 0.0001, ]
nrow(abundant_otus)

# Descriptive statistics of the subset of most abundant OTUs
mean(colSums(abundant_otus))
sd(colSums(abundant_otus))
range(range(colSums(abundant_otus)))
summary(colSums(abundant_otus))
rownames(abundant_otus)

# Correlation analysis between the first 3 axis of weighted UniFrac PCoA and the most abundant OTUs
# For each axis
# 1. Spearman's correlation tests are done between each axis and OTU
# 2. A table is generated with rho and p-value of each result of point 1
# 3. p-value correction is made using qvalue
# 4. Significant correlations are selected if q<0.05
# 5. A table is generated with the taxonomy, rho, p-value and q-value of the significant OTUs

# PC1
cor_entero_PC1 = apply(X = abundant_otus, MARGIN = 1, 
                       FUN = function(x) cor.test(x, e.pcoa$points[, 1], m = "s"))

cor_PC1 = sapply(X = cor_entero_PC1,
                 FUN = function(x) cbind(rho = as.vector(x$estimate),
                                         p =as.vector(x$p.value)))

cor_PC1 = as.data.frame(t(cor_PC1))
colnames(cor_PC1) = c("rho", "p")
q_value = as.vector(qvalue(cor_PC1$p)$qvalue)
cor_PC1$q = q_value
sig_cor_PC1 = cor_PC1[cor_PC1$q <= 0.05,]
result_cor_PC1 = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% row.names(sig_cor_PC1),2], sig_cor_PC1)
result_cor_PC1 = result_cor_PC1[abs(result_cor_PC1$rho) >= 0.3,]
sig_otus_PC1 = abundant_otus[cor_PC1[,3] <= 0.05,]
nrow(result_cor_PC1)

# PC2
cor_entero_PC2 = apply(X = abundant_otus, MARGIN = 1, 
                       FUN = function(x) cor.test(x, e.pcoa$points[, 2], m = "s"))
cor_PC2 = sapply(X = cor_entero_PC2, 
                 FUN = function(x) cbind(rho = as.vector(x$estimate), 
                                         p =as.vector(x$p.value)))

cor_PC2 = as.data.frame(t(cor_PC2))
colnames(cor_PC2) = c("rho", "p")
q_value = as.vector(qvalue(cor_PC2$p)$qvalue)
cor_PC2$q = q_value
sig_cor_PC2 = cor_PC2[cor_PC2$q <= 0.05,]
result_cor_PC2 = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% row.names(sig_cor_PC2),2], sig_cor_PC2)
result_cor_PC2 = result_cor_PC2[abs(result_cor_PC2$rho) >= 0.3,]
sig_otus_PC2 = abundant_otus[cor_PC2[,3] <= 0.05,]
nrow(result_cor_PC2)

# PC3
cor_entero_PC3 = apply(X = abundant_otus, MARGIN = 1, 
                       FUN = function(x) cor.test(x, e.pcoa$points[, 3],
                                                  m = "s"))
cor_PC3 = sapply(X = cor_entero_PC3, 
                 FUN = function(x) cbind(rho = as.vector(x$estimate),
                                         p =as.vector(x$p.value)))

cor_PC3 = as.data.frame(t(cor_PC3))
colnames(cor_PC3) = c("rho", "p")
q_value = as.vector(qvalue(cor_PC3$p)$qvalue)
cor_PC3$q = q_value
sig_cor_PC3 = cor_PC3[cor_PC3$q <= 0.05,]
result_cor_PC3 = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% row.names(sig_cor_PC3),2], sig_cor_PC3)
result_cor_PC3 = result_cor_PC3[abs(result_cor_PC3$rho) >= 0.3,]
sig_otus_PC3 = abundant_otus[cor_PC3[,3] <= 0.05,]
nrow(result_cor_PC3)

# Table with the OTUs correlated with at least one PCoA axis (Table S4)
all_axes_otus_names = sort(unique(c(rownames(result_cor_PC1), 
                                    rownames(result_cor_PC2), 
                                    rownames(result_cor_PC3))))

all_axes_otus = abundant_otus[all_axes_otus_names,]
cor_all_PC1 = apply(X = all_axes_otus, MARGIN = 1, 
                    FUN = function(x) cor(x, e.pcoa$points[, 1], m = "s"))

cor_all_PC2 = apply(X = all_axes_otus, MARGIN = 1, 
                    FUN = function(x) cor(x, e.pcoa$points[, 2], m = "s"))

cor_all_PC3 = apply(X = all_axes_otus, MARGIN = 1, 
                    FUN = function(x) cor(x, e.pcoa$points[, 3], m = "s"))

result_cor_all = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% row.names(all_axes_otus),2], 
                            round(cor_all_PC1, 2), 
                            round(cor_all_PC2, 2), 
                            round(cor_all_PC3, 2))


# Select phylotypes characteristic of Western and non-Western microbiota ----
bacteroides = d.genus[,"k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides"]
rownames(bacteroides) = NULL

prevotella = d.genus[,"k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella"]
rownames(prevotella) = NULL

ruminococcus = d.genus[,"k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Ruminococcus"]
rownames(ruminococcus) = NULL

treponema = d.genus[,"k__Bacteria;p__Spirochaetes;c__Spirochaetes;o__Spirochaetales;f__Spirochaetaceae;g__Treponema"]
rownames(treponema) = NULL

bifidobacterium = d.genus[,"k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Bifidobacteriales;f__Bifidobacteriaceae;g__Bifidobacterium"]
rownames(bifidobacterium) = NULL

cytophagales = d.order[,"k__Bacteria;p__Bacteroidetes;c__Cytophagia;o__Cytophagales"]
rownames(cytophagales) = NULL

barnesiella = d.genus[,"k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__[Barnesiellaceae];g__Barnesiella"]
rownames(barnesiella) = NULL

select_entero_table = data.frame(ID = rownames(d.genus),
                                 Prevotella = prevotella,
                                 Bacteroides = bacteroides,
                                 Ruminococcus = ruminococcus,
                                 Treponema = treponema,
                                 Bifidobacterium = bifidobacterium,
                                 Cytophagales = cytophagales, 
                                 Barnesiella = barnesiella, 
                                 PC1 = e.pcoa$points[, 1], 
                                 PC2 = e.pcoa$points[, 2], 
                                 PC3 = e.pcoa$points[, 3])

# Test that the mean phylotype abundance is different from that found in the
# meta-analysis of benchmark datasets from curatedMetagenomicData
# (run the accompanying "curated_metagenomes.R" script to get these values)
t.test(prevotella, mu = 0.245)
t.test(bifidobacterium, mu = 0.074)
t.test(bacteroides, mu = 0.230)
t.test(treponema, mu = 0.021)
t.test(cytophagales, mu = 6e-05)
t.test(barnesiella, mu = 0.0126)

# Test that the mean phylotype abundance is different from 0
t.test(bacteroides, mu = 0, alternative = "g")
t.test(prevotella, mu = 0, alternative = "g")
t.test(bifidobacterium, mu = 0, alternative = "g")
t.test(treponema, mu = 0, alternative = "g")
t.test(cytophagales, mu = 0, alternative = "g")
t.test(barnesiella, mu = 0, alternative = "g")

# Western - non-Western histograms (Figure 2)
hist_bacteroides = ggplot(select_entero_table, aes(Bacteroides)) + 
  geom_histogram(boundary = 0, bins = 100) + theme_gray() +
  labs(x = "Bacteroides abundance", y = "Number of subjects")

hist_prevotella = ggplot(select_entero_table, aes(Prevotella)) +
  geom_histogram(boundary = 0, bins = 100) + theme_gray() +
  labs(x = "Prevotella abundance", y = "Number of subjects")

hist_bifidobacterium = ggplot(select_entero_table, aes(Bifidobacterium)) + 
  geom_histogram(boundary = 0, bins = 100) + theme_gray() + 
  labs(x = "Bifidobacterium abundance", y = "Number of subjects")

hist_treponema = ggplot(select_entero_table, aes(Treponema)) + 
  geom_histogram(boundary = 0, bins = 100) + theme_gray() + 
  labs(x = "Treponema abundance", y = "Number of subjects")

hist_cytophagales = ggplot(select_entero_table, aes(Cytophagales)) + 
  geom_histogram(boundary = 0, bins = 100) + theme_gray() + 
  labs(x = "Cytophagales abundance", y = "Number of subjects")

hist_barnesiella = ggplot(select_entero_table, aes(Barnesiella)) + 
  geom_histogram(boundary = 0, bins = 100) + theme_gray() +
  labs(x = "Barnesiella abundance", y = "Number of subjects")

plot_grid(hist_bacteroides, hist_bifidobacterium,
          hist_barnesiella, hist_prevotella, hist_treponema, 
          labels="AUTO", nrow = 2, ncol = 3)

# Prevotella-Bacteroides ratio in the PCoA of weighted UniFrac (Figure SR1)
qplot(PC1, PC2, data = select_entero_table, 
      colour = bacteroides/(bacteroides+prevotella) ) + 
  scale_colour_gradientn(colours=(rainbow(10)), name = "") + 
  labs(x = "PCo1 16.4%", y = "PCo2 13.8%")

cor.test(prevotella, bacteroides, method = "s")


# Coabundance group analyses (CAGs) ----
# Heatmap of OTUs with median abundance >=0.01% OTUs and generation of CAGs
# 1. Spearman's correlations between pairs of OTU are obtained
# 2. A dendrogram is constructed using the Ward's method
# 3. The heatmap is generated.
# In this case, multiple k's were tested for the number of CAGs (not shown)
# Using k = 5

# spearman_tree & CAG_heat_tree are the same but in different order
spearman_matrix = cor(t(abundant_otus), method = "spearman")
spearman_tree = hclust(dist(spearman_matrix), method = "ward")
plot(spearman_tree)
CAGs = factor(cutree(spearman_tree, k = 5))

# Figure S1
CAG_heatmap = aheatmap(spearman_matrix, color = "-Spectral:100", scale = "none",
                       breaks = NA, Rowv=T, Colv=T, width=30, height=8, cexRow=1,
                       hclustfun = "ward", main = "OTU co-ocurrences",
                       distfun = "euclidean", annRow = CAGs)

# Validation of CAGs 
# Randomly split the dataset, compute a correlation matrix and run a Mantel test
subsample = sample(c(1:441), size = 220, replace = F)
train_data = abundant_otus[,subsample]
test_data = abundant_otus[,-subsample]

train_spearman = cor(t(train_data), method = "spearman")
test_spearman = cor(t(test_data), method = "spearman")

mantel(as.dist(train_spearman) ~ as.dist(test_spearman), nperm = 10000, nboot = 10000)

# With 5 CAGs
# 1. The OTUs that belong to each CAG are selected and their abundance is combined.
# In this way, the abundance of each CAG in a given sample is equal to the sum of the abundance of the OTUs that compose it.
# 2. A table is generated with the abundances of each CAG in each sample
CAG_n5 = factor(cutree(tree = spearman_tree, k = 5))
CAG_n5_a = abundant_otus[CAG_n5 == 1,]
CAG_n5_b = abundant_otus[CAG_n5 == 2,]
CAG_n5_c = abundant_otus[CAG_n5 == 3,]
CAG_n5_d = abundant_otus[CAG_n5 == 4,]
CAG_n5_e = abundant_otus[CAG_n5 == 5,]
CAG_n5_matrix = data.frame(a = colSums(CAG_n5_a), b = colSums(CAG_n5_b),
                           c = colSums(CAG_n5_c), d = colSums(CAG_n5_d), 
                           e = colSums(CAG_n5_e))

# To identify the OTUs that belong to each CAG, the names of the OTUs are compared with their respective taxonomy
OTUs_n5_a = row.names(CAG_n5_a)
OTUs_n5_b = row.names(CAG_n5_b)
OTUs_n5_c = row.names(CAG_n5_c)
OTUs_n5_d = row.names(CAG_n5_d)
OTUs_n5_e = row.names(CAG_n5_e)

matrix_n5_a = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n5_a),2], CAG_n5_a)
matrix_n5_b = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n5_b),2], CAG_n5_b)
matrix_n5_c = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n5_c),2], CAG_n5_c)
matrix_n5_d = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n5_d),2], CAG_n5_d)
matrix_n5_e = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% (OTUs_n5_e),2], CAG_n5_e)

# To identify the dominant OTUs in each CAG, the median abundance is calculated
dominant_n5_a = apply(CAG_n5_a, MARGIN = 1, FUN = median)
sort(dominant_n5_a, decreasing = T)
dominant_n5_b = apply(CAG_n5_b, MARGIN = 1, FUN = median)
sort(dominant_n5_b, decreasing = T)
dominant_n5_c = apply(CAG_n5_c, MARGIN = 1, FUN = median)
sort(dominant_n5_c, decreasing = T)
dominant_n5_d = apply(CAG_n5_d, MARGIN = 1, FUN = median)
sort(dominant_n5_d, decreasing = T)
dominant_n5_e = apply(CAG_n5_e, MARGIN = 1, FUN = median)
sort(dominant_n5_e, decreasing = T)

# Mean abundance and standard deviation of each CAG
n5_CAG_mean = apply(CAG_n5_matrix, MARGIN = 2, mean)
n5_CAG_sd = apply(CAG_n5_matrix, MARGIN = 2, sd)
round(n5_CAG_mean*100, 2)
round(n5_CAG_sd*100, 2)

# Graphing the 5 CAGs
# A table is generated with the abundance of each CAG and the values of the PCoA
entero_table = data.frame(PC1 = e.pcoa$points[, 1], PC2 = e.pcoa$points[, 2],
                          PC3 = e.pcoa$points[, 3])
n5_table = data.frame(PC1 = entero_table$PC1, PC2 = entero_table$PC2, 
                      PC3 = entero_table$PC3, a = CAG_n5_matrix$a, 
                      b = CAG_n5_matrix$b, c = CAG_n5_matrix$c, 
                      d = CAG_n5_matrix$d, e = CAG_n5_matrix$e)

# Boxplot of CAG abundance in the whole dataset
melt_n5 = melt(data = CAG_n5_matrix)
CAG = factor(CAG_n5, 
             labels = c("Prevotella", "Lachnospiraceae", "Pathogen",
                        "Akkermansia-Bacteroidales", "Ruminococcaceae"))

# Figure 3A
CAG_boxplot = ggplot(melt_n5, aes(factor(variable), value)) + geom_boxplot() + 
  theme_gray() + theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  labs(x = "CAG", y = "Relative abundance") + 
  scale_fill_manual(breaks = c("Prevotella", "Lachnospiraceae", "Pathogen", "Akkermansia-Bacteroidales", "Ruminococcaceae")) + 
  scale_x_discrete(labels = c("Prevotella", "Lachnospiraceae", "Pathogen", "Akkermansia-Bacteroidales", "Ruminococcaceae"))

# Heatmap with the dendrograms and colors according to the CAGs
heatmap_n5 = aheatmap(spearman_matrix, color = "-Spectral:100", scale = "none", 
                      breaks = NA, Rowv=T, Colv=T, width=30, height=8, cexRow=1,
                      hclustfun = "ward", main = "OTU co-ocurrences", 
                      distfun = "euclidean", annRow = CAG, annCol = CAG, 
                      treeheight = 120)

# Abundance of each CAG along the enterogradient
# Using the weighted UniFrac PCoA
# 2D plots
prevotella.cag = qplot(PC1, PC2, data = n5_table, colour = a) + 
  scale_colour_gradientn(colours=rev(rainbow(10)), name = "Relative\nabundance") + 
  labs(x = "PCo1 16.4%", y = "PCo2 13.8%", title = "Prevotella CAG")

lacho.cag = qplot(PC1, PC2, data = n5_table, colour = b) + 
  scale_colour_gradientn(colours=rev(rainbow(10)), name = "Relative\nabundance") + 
  labs(x = "PCo1 16.4%", y = "PCo2 13.8%", title = "Lachnospiraceae CAG")

pathogen.cag = qplot(PC1, PC2, data = n5_table, colour = c) + 
  scale_colour_gradientn(colours=rev(rainbow(10)), name = "Relative\nabundance") + 
  labs(x = "PCo1 16.4%", y = "PCo2 13.8%", title = "Pathogen CAG")

akk.cag = qplot(PC1, PC2, data = n5_table, colour = d) + 
  scale_colour_gradientn(colours=rev(rainbow(10)), name = "Relative\nabundance") + 
  labs(x = "PCo1 16.4%", y = "PCo2 13.8%", title = "Akkermansia-Bacteroidales CAG")

rumino.cag = qplot(PC1, PC2, data = n5_table, colour = e) + 
  scale_colour_gradientn(colours=rev(rainbow(10)), name = "Relative\nabundance") + 
  labs(x = "PCo1 16.4%", y = "PCo2 13.8%", title = "Ruminococcaceae CAG")

plot_grid(CAG_boxplot, prevotella.cag, lacho.cag, pathogen.cag, akk.cag, rumino.cag, 
          labels=c("A", "B", "C", "D", "E", "F"), ncol = 2, nrow = 3)

# pseudo-3D plots (Figure 3B-F)
scatter3D(n5_table$PC1, n5_table$PC2, n5_table$PC3, colvar = n5_table$a, phi = 40, theta = 40, bty = "g", pch=20, cex=1, clab = "Relative\nabundance", xlab=paste("PCo1",e.PC1, "%"), ylab=paste("PCo2",e.PC2, "%"), zlab=paste("PCo3",e.PC3, "%"), main = "CAG-Prevotella", par(mar=c(2.1,2.1,2.1,2.1)))
scatter3D(n5_table$PC1, n5_table$PC2, n5_table$PC3, colvar = n5_table$b, phi = 40, theta = 40, bty = "g", pch=20, cex=1, clab = "Relative\nabundance", xlab=paste("PCo1",e.PC1, "%"), ylab=paste("PCo2",e.PC2, "%"), zlab=paste("PCo3",e.PC3, "%"), main = "CAG-Lachnospiraceae", par(mar=c(2.1,2.1,2.1,2.1)))
scatter3D(n5_table$PC1, n5_table$PC2, n5_table$PC3, colvar = n5_table$c, phi = 40, theta = 40, bty = "g", pch=20, cex=1, clab = "Relative\nabundance", xlab=paste("PCo1",e.PC1, "%"), ylab=paste("PCo2",e.PC2, "%"), zlab=paste("PCo3",e.PC3, "%"), main = "CAG-Pathogen", par(mar=c(2.1,2.1,2.1,2.1)))
scatter3D(n5_table$PC1, n5_table$PC2, n5_table$PC3, colvar = n5_table$d, phi = 40, theta = 40, bty = "g", pch=20, cex=1, clab = "Relative\nabundance", xlab=paste("PCo1",e.PC1, "%"), ylab=paste("PCo2",e.PC2, "%"), zlab=paste("PCo3",e.PC3, "%"), main = "CAG-Akkermansia-\nBacteroidales", par(mar=c(2.1,2.1,2.1,2.1)))
scatter3D(n5_table$PC1, n5_table$PC2, n5_table$PC3, colvar = n5_table$e, phi = 40, theta = 40, bty = "g", pch=20, cex=1, clab = "Relative\nabundance", xlab=paste("PCo1",e.PC1, "%"), ylab=paste("PCo2",e.PC2, "%"), zlab=paste("PCo3",e.PC3, "%"), main = "CAG-Ruminococcaceae", par(mar=c(2.1,2.1,2.1,2.1)))

# Single-CAG dominated microbiota analyses ----

# Selection of individuals with single-CAG dominated microbiota

# The 95th percentile of the abundance distribution of each CAG is computed
pole_percentile_a = quantile(CAG_n5_matrix$a, probs = 0.95)
pole_percentile_b = quantile(CAG_n5_matrix$b, probs = 0.95)
pole_percentile_c = quantile(CAG_n5_matrix$c, probs = 0.95)
pole_percentile_d = quantile(CAG_n5_matrix$d, probs = 0.95)
pole_percentile_e = quantile(CAG_n5_matrix$e, probs = 0.95)

# Filtering metadata file 
poles_metadata_a = microbio.meta[CAG_n5_matrix$a >= pole_percentile_a,]
poles_metadata_a = poles_metadata_a[-4,] #MI_252_H is classified in 2 groups. It remains in the one with higher abundance
poles_metadata_b = microbio.meta[CAG_n5_matrix$b >= pole_percentile_b,]
poles_metadata_c = microbio.meta[CAG_n5_matrix$c >= pole_percentile_c,]
poles_metadata_d = microbio.meta[CAG_n5_matrix$d >= pole_percentile_d,]
poles_metadata_e = microbio.meta[CAG_n5_matrix$e >= pole_percentile_e,]

# A new metadata table is generated for each CAG and with a CAG variable in a new column
poles_metadata_a = cbind(poles_metadata_a, 
                         pole = rep("Prevotella", nrow(poles_metadata_a)))

poles_metadata_b = cbind(poles_metadata_b, 
                         pole = rep("Lachnospiraceae", nrow(poles_metadata_b)))

poles_metadata_c = cbind(poles_metadata_c, 
                         pole = rep("Pathogen", nrow(poles_metadata_c)))

poles_metadata_d = cbind(poles_metadata_d, 
                         pole = rep("Akkermansia-Bacteroidales", 
                                    nrow(poles_metadata_d)))

poles_metadata_e = cbind(poles_metadata_e, 
                         pole = rep("Ruminococcaceae", nrow(poles_metadata_e)))

poles_matrix_n5 = rbind(poles_metadata_a, poles_metadata_b, poles_metadata_c, 
                        poles_metadata_d, poles_metadata_e)


# Statistical analyses of individuals with single-CAG dominated microbiota in terms of biochemical and anthropometric variables ----

# Continuous variables
# Selection of evaluated variables
metadata_variables = c("bmi", "age","body_fat", "systolic_bp", "diastolic_bp",
                       "LDL", "cholesterol", "triglycerides", "hsCRP", "glucose", 
                       "insulin", "HOMA_IR", "per_protein", "per_animal_protein", 
                       "per_total_fat", "per_saturated_fat", "per_monoinsaturated_fat",
                       "per_polyunsaturated_fat", "per_carbohydrates", "fiber",
                       "glycosylated_hg")

poles_aov = lapply(metadata_variables,
                   FUN = function(x)  Anova(aov(log(poles_matrix_n5[,x]) ~ poles_matrix_n5$pole )))

poles_aov_pvalues = unlist(lapply(poles_aov, 
                                  FUN = function(x) x$`Pr(>F)`[1]))

# Adiponectin cannot be log-transformed. Thus, 4th root transformation is used
poles_adiponectin = aov((poles_matrix_n5$adiponectin)^(1/4) ~ poles_matrix_n5$pole)
adiponectin_aov = Anova(poles_adiponectin)

# HDL and waist circumference are evaluated adjusting by sex
poles_hdl_sex = aov(log(poles_matrix_n5$HDL) ~ poles_matrix_n5$pole + poles_matrix_n5$sex)
hdl_sex_aov = Anova(poles_hdl_sex)

poles_waist_sex = aov(log(poles_matrix_n5$waist) ~ poles_matrix_n5$pole + poles_matrix_n5$sex)
waist_sex_aov = Anova(poles_waist_sex)

# Categorical variables
# City, sex, medicament consumption, stool consistency, occult blood test and age

poles_city = table(poles_matrix_n5$pole, poles_matrix_n5$city)
t(round(prop.table(poles_city, margin = 1)*100, 1))
t(round(prop.table(poles_city, margin = 2)*100, 1))
city_chi = chisq.test(poles_matrix_n5$pole, poles_matrix_n5$city, 
                      simulate.p.value = T)


poles_sex = table(poles_matrix_n5$pole, poles_matrix_n5$sex)
t(round(prop.table(poles_sex, margin = 1)*100, 1))
t(round(prop.table(poles_sex, margin = 2)*100, 1))
sex_chi = chisq.test(poles_matrix_n5$pole, poles_matrix_n5$sex, 
                     simulate.p.value = T)


poles_stool = table(poles_matrix_n5$pole, poles_matrix_n5$stool_consistency)
t(round(prop.table(poles_stool, margin = 1)*100, 1))
t(round(prop.table(poles_stool, margin = 2)*100, 1))
stool_chi = chisq.test(poles_matrix_n5$pole, poles_matrix_n5$stool_consistency, 
                       simulate.p.value = T)


poles_meds = table(poles_matrix_n5$pole, poles_matrix_n5$medicament)
t(round(prop.table(poles_meds, margin = 1)*100, 1))
t(round(prop.table(poles_meds, margin = 2)*100, 1))
meds_chi = chisq.test(poles_matrix_n5$pole, poles_matrix_n5$medicament, 
                      simulate.p.value = T)

poles_s.oculta = table(poles_matrix_n5$pole, poles_matrix_n5$hiden_blood)
t(round(prop.table(poles_s.oculta, margin = 1)*100, 1))
t(round(prop.table(poles_s.oculta, margin = 2)*100, 1))
oculta_chi = chisq.test(poles_matrix_n5$pole, poles_matrix_n5$hiden_blood, 
                        simulate.p.value = T)

tested_variables = c(metadata_variables, "waist", "HDL", "adiponectin", "city", 
                     "sex",  "stool_consistency", "medicament", "hiden_blood") 

complete_pvalues = c(poles_aov_pvalues,
                     waist_sex_aov$`Pr(>F)`[1], hdl_sex_aov$`Pr(>F)`[1], 
                     adiponectin_aov$`Pr(>F)`[1], city_chi$p.value, sex_chi$p.value,
                     stool_chi$p.value, meds_chi$p.value, oculta_chi$p.value)

poles_aov_qvalues = qvalue(complete_pvalues)$qvalues

# Table with mean and standard deviation of each evaluated variable (Table 1)
mean_a = apply(poles_metadata_a[, tested_variables[1:24]], MARGIN = 2,
               FUN = mean, na.rm = T)

sd_a = apply(poles_metadata_a[, tested_variables[1:24]], MARGIN = 2,
             FUN = sd, na.rm = T)

mean_b = apply(poles_metadata_b[, tested_variables[1:24]], MARGIN = 2, 
               FUN = mean, na.rm = T)

sd_b = apply(poles_metadata_b[, tested_variables[1:24]], MARGIN = 2, 
             FUN = sd, na.rm = T)

mean_c = apply(poles_metadata_c[, tested_variables[1:24]], MARGIN = 2, 
               FUN = mean, na.rm = T)

sd_c = apply(poles_metadata_c[, tested_variables[1:24]], MARGIN = 2, 
             FUN = sd, na.rm = T)

mean_d = apply(poles_metadata_d[, tested_variables[1:24]], MARGIN = 2,
               FUN = mean, na.rm = T)

sd_d = apply(poles_metadata_d[, tested_variables[1:24]], MARGIN = 2, 
             FUN = sd, na.rm = T)

mean_e = apply(poles_metadata_e[, tested_variables[1:24]], MARGIN = 2, 
               FUN = mean, na.rm = T)

sd_e = apply(poles_metadata_e[, tested_variables[1:24]], MARGIN = 2, 
             FUN = sd, na.rm = T)

poles_final_table = data.frame(mean_a, sd_a, mean_b, sd_b,
                               mean_c, sd_c, mean_d, sd_d, mean_e, sd_e)
round(poles_final_table, 1)
data.frame(vars = tested_variables, 
           p = round(complete_pvalues, 3), q = round(poles_aov_qvalues, 3))

# Alpha-diversity analyses in individuals with single-CAG dominated microbiota ----

# Filtering the OTU table
count_poles_a  = microbio.rare[rownames(poles_metadata_a),]
count_poles_b  = microbio.rare[rownames(poles_metadata_b),]
count_poles_c  = microbio.rare[rownames(poles_metadata_c),]
count_poles_d  = microbio.rare[rownames(poles_metadata_d),]
count_poles_e  = microbio.rare[rownames(poles_metadata_e),]
count_poles = rbind(count_poles_a, count_poles_b, count_poles_c, count_poles_d, count_poles_e)

# Alpha-diversity metrics calculation
poles_richness = diversityresult(x = count_poles, 
                                 index = "richness", method = "each site")
poles_shannon_i = diversityresult(x = count_poles,
                                  index = "Shannon", method = "each site")
poles_shannon_Je = diversityresult(x = count_poles,
                                   index = "Jevenness", method = "each site")
alpha_div_extremes =cbind(poles_richness, 
                          poles_shannon_i, poles_shannon_Je, 
                          pole = poles_matrix_n5$pole)

# Hypothesis testing of alpha-diversity
# Richness
aov_richness = aov(richness ~ pole, data = alpha_div_extremes)
summary(aov_richness)
TukeyHSD(aov_richness)

# Shannon Index
aov_shannon_i = aov(Shannon ~ pole, data = alpha_div_extremes)
summary(aov_shannon_i)
TukeyHSD(aov_shannon_i)

# Shannon J evenness (Piellou's Index)
aov_shannon_Je = aov(Jevenness ~ pole, data = alpha_div_extremes)
summary(aov_shannon_Je)
TukeyHSD(aov_shannon_Je)

# Alpha-diversity boxplots

diversity_melt = melt(data = alpha_div_extremes)

# Richness
box_richness = ggplot(diversity_melt[diversity_melt$variable == "richness",], 
                      aes(x = pole, y = value))
b.richness = box_richness +
  geom_boxplot(notch = T) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none") +
  labs(x = "CAG", y = "Species richness") + 
  theme_grey()


# Shannon Index
box_shannon_i = ggplot(diversity_melt[diversity_melt$variable == "Shannon",],
                       aes(x = pole, y = value))
b.shannon = box_shannon_i + geom_boxplot(notch = T) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none") +
  labs(x = "CAG", y = "Shannon Index") + 
  theme_grey()


# Shannon J evenness
box_shannon_Je = ggplot(diversity_melt[diversity_melt$variable == "Jevenness",], 
                        aes(x = pole, y = value))
b.shannon_Je = box_shannon_Je + geom_boxplot(notch = T) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none") + 
  labs(x = "CAG", y = "Pielou's J") +
  theme_grey()


# KO richness
data<-read.table("./files/tax4fun.Rdata",header=T)
metagenome_diversity_melt = melt(data = data)
box_ko_richness = ggplot(metagenome_diversity_melt[metagenome_diversity_melt$variable == "N_KO",],
                         aes(x = reorder(pole, -value, FUN=median), y = value))
box_ko_richness + geom_boxplot(notch = T) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none") + 
  labs(x = "CAG", y = "KO richness") + 
  theme_grey()

# Ordering poles in the x-axis in a given order
metagenome_diversity_melt$pole<-factor(metagenome_diversity_melt$pole,
                                       levels=c("Prevotella", "Lachnospiraceae", 
                                                "Pathogen", "Akkermansia-Bacteroidales", 
                                                "Ruminococcaceae"))

box_ko_richness = ggplot(metagenome_diversity_melt[metagenome_diversity_melt$variable == "N_KO",], 
                         aes(x = metagenome_diversity_melt$pole, y = value))

b.ko_richness = box_ko_richness + geom_boxplot(notch = T) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none") + 
  labs(x = "CAG", y = "KO richness") + theme_grey()
b.ko_richness

# Figure S2
plot_grid(b.richness, b.shannon, b.shannon_Je, b.ko_richness,
          labels=c("A", "B", "C", "D"), ncol = 1, nrow = 4)


# Beta-diversity analyses in individuals with single-CAG dominated microbiota ----
# Graphs

# Location of single-CAG dominated microbiota in the complete Weighted UniFrac pointcloud
poles_center = sapply(rownames(microbio.meta),
                      function(x) ifelse(x[1] %in% rownames(poles_matrix_n5), 
                                         yes = as.character(poles_matrix_n5[x, "pole"]),
                                         no = "Center"))

poles_center = factor(poles_center,
                      levels = c("Prevotella", "Lachnospiraceae",
                                 "Pathogen", "Akkermansia-Bacteroidales", 
                                 "Ruminococcaceae"))

graph.poles_center = ggplot(cbind(w, poles_center), aes(V1, V2)) +
  geom_point(aes(color = poles_center)) +
  labs(x = paste("PCo1",e.PC1, "%"), y = paste("PCo2",e.PC2, "%"), title = "PCoA Weighted UniFrac") + 
  scale_color_discrete(name = "CAG pole")
graph.poles_center

# Reduce the dataset to the single-CAG dominated microbiota
samples_poles = rownames(poles_matrix_n5)

# Evaluated variables
BMI_pole = poles_matrix_n5$bmi_class
city_pole = poles_matrix_n5$city
sex_pole = poles_matrix_n5$sex 
age_pole = poles_matrix_n5$age_range
pole = poles_matrix_n5$pole

# weighted UniFrac
dw_poles = dw[samples_poles, samples_poles]
adonis(as.dist(dw_poles) ~ pole)
adonis(as.dist(dw_poles) ~ BMI_pole)
adonis(as.dist(dw_poles) ~ city_pole)
adonis(as.dist(dw_poles) ~ sex_pole)
adonis(as.dist(dw_poles) ~ age_pole)

# unweighted UniFrac
du_poles = du[samples_poles, samples_poles]
adonis(as.dist(du_poles) ~ pole)

# PCoA Graphs
# PCoA weighted
poles.pcoa_weighted = cmdscale(as.dist(dw_poles), k=5, eig = T)
poles.PC1_weighted = round(poles.pcoa_weighted$eig[1]/sum(e.pcoa$eig), 4)* 100
poles.PC2_weighted = round(poles.pcoa_weighted$eig[2]/sum(e.pcoa$eig), 4)* 100
poles.PC3_weighted = round(poles.pcoa_weighted$eig[3]/sum(e.pcoa$eig), 4)* 100

# PCoA unweighted
poles.pcoa_unweighted = cmdscale(as.dist(du_poles), k=5, eig = T)
poles.PC1_unweighted = round(poles.pcoa_unweighted$eig[1]/sum(e.pcoa$eig), 4)* 100
poles.PC2_unweighted = round(poles.pcoa_unweighted$eig[2]/sum(e.pcoa$eig), 4)* 100
poles.PC3_unweighted = round(poles.pcoa_unweighted$eig[3]/sum(e.pcoa$eig), 4)* 100

pole_weighted = as.data.frame(poles.pcoa_weighted$points) 
pole_unweighted = as.data.frame(poles.pcoa_unweighted$points)

# Graphs of weighted and unweighted UniFrac PCoA only with single-CAG dominated microbiota
# Color and ellipses according to different variables

graph.pole_weighted = ggplot(cbind(pole_weighted, pole), aes(V1, V2)) + 
  geom_point(aes(color = pole)) + stat_ellipse(type = "t", mapping = aes(color = pole), level = 0.75) + 
  labs(x = "PCo1 7.1%", y = "PCo2 6.0%") + 
  scale_color_discrete(name = "CAG pole")

graph.pole_unweighted = ggplot(cbind(pole_unweighted, pole), aes(V1, V2)) + 
  geom_point(aes(color = pole)) + 
  stat_ellipse(type = "t", mapping = aes(color = pole), level = 0.75) + 
  labs(x = "PCo1 3.4%", y = "PCo2 2.2%") + scale_color_discrete(name = "CAG pole")

graph.pole_raw = ggplot(cbind(pole_weighted, BMI_pole), aes(V1, V2)) +
  geom_point() +
  labs(x = "PCo1 7.1%", y = "PCo2 6.0%")

graph.pole_city = ggplot(cbind(pole_weighted, city_pole), aes(V1, V2)) +
  geom_point(aes(color = city_pole)) + 
  stat_ellipse(type = "t", mapping = aes(color = city_pole), level = 0.75) + 
  labs(x = "PCo1 7.1%", y = "PCo2 6.0%") +
  scale_color_discrete(name = "City")

graph.pole_bmi = ggplot(cbind(pole_weighted, BMI_pole), aes(V1, V2)) + 
  geom_point(aes(color = BMI_pole)) +
  stat_ellipse(type = "t", mapping = aes(color = BMI_pole), level = 0.75) + 
  labs(x = "PCo1 7.1%", y = "PCo2 6.0%") +
  scale_color_discrete(name = "BMI category", labels = c("Lean", "Obese", "Overweight"))

graph.pole_sex = ggplot(cbind(pole_weighted, sex_pole), aes(V1, V2)) + 
  geom_point(aes(color = sex_pole)) + stat_ellipse(type = "t", mapping = aes(color = sex_pole), level = 0.75) + 
  labs(x = "PCo1 7.1%", y = "PCo2 6.0%") + 
  scale_color_discrete(name = "Sex", labels = c("Male", "Female"))

graph.pole_age = ggplot(cbind(pole_weighted, age_pole), aes(V1, V2)) +
  geom_point(aes(color = age_pole)) + 
  stat_ellipse(type = "t", mapping = aes(color = age_pole), level = 0.75) + 
  labs(x = "PCo1 7.1%", y = "PCo2 6.0%") +
  scale_color_discrete(name = "Age group", labels = c("18 - 40", "41 - 62"))

# Figure 4
plot_grid(graph.pole_weighted,  graph.pole_unweighted, graph.pole_city,
          graph.pole_bmi, graph.pole_sex, graph.pole_age, 
          labels=c("A", "B", "C", "D", "E", "F"), ncol = 2, nrow = 3)

# # Metagenome Prediction ----
# # The OTU table in Qiime format should be imported and then Tax4Fun is run
# # This OTU table needs to be classified with the SILVA database
# poles_biom = importQIIMEData(inputFiles = "./files/microbio_poles.biom")
# # Reference data should be downloaded from http://tax4fun.gobics.de/
# metagenome_prediction = Tax4Fun(Tax4FunInput = poles_biom, folderReferenceData = "SILVA123")
# ko_profile = metagenome_prediction$Tax4FunProfile
# write.table(t(ko_profile), file="ko_profile.txt", sep = "\t", quote = F)
# metadata_poles = data.frame(Samples = rownames(poles_matrix_n5), CAG = poles_matrix_n5$pole)
# 







