# Code used in the meta-analysis of Western and non-Western microbiota presented in the paper:
# de la Cuesta-Zuluaga J, Corrales-Agudelo V, Velásquez-Mejía EP, Carmona JA, Abad JM, Escobar JS. 2018.
# Gut microbiota is associated with obesity and cardiometabolic disease in a population in the midst of Westernization. 
# Scientific Reports. doi: 10.1038/s41598-018-29687-x.
# www.nature.com/articles/s41598-018-29687-x

# Libraries ----
library(curatedMetagenomicData)
library(phyloseq)
library(ggplot2)
library(NMF)
library(plyr)

# Download data
# stool.data <- curatedMetagenomicData("*.metaphlan_bugs_list.stool", dryrun=FALSE)

# Simple phyloseq manipulation ----
#save the tax data as phyloseq object
raw_data = ExpressionSet2phyloseq(HMP_2012.metaphlan_bugs_list.stool()) 
# Subset of healthy males
data_object = subset_samples(raw_data, disease=='healthy' & gender=='male') 
data_object #Summary of the object

# This next part is tricky. If you want the data at a given tax rank,
# you subset by selecting !is.na(rank) & is.na(next_rank)
#
# For example, to create two objects, one with family and the other with genus 
# abundance data
#
# Remember that the curatedMetagenomedata taxa data goes up to strain level
family_data = subset_taxa(data_object, !is.na(Family) & is.na(Genus))
genus_data = subset_taxa(data_object, !is.na(Genus) & is.na(Species))

# To extract data from a given dataset, for example, the genus_data object
genus_data_abundance = otu_table(genus_data) # Phylotype abundance
genus_data_taxonomy = tax_table(genus_data) # Taxonomic classification
genus_data_metadata = sample_data(genus_data) #Metadata

# To extract the abundance of a phylotype at a given taxonomic rank
Prevotella_HMP = otu_table(subset_taxa(genus_data, Genus=="Prevotella"))
# For easier manipulation, it can be transformed
Prevotella_HMP = as.vector(Prevotella_HMP)
summary(Prevotella_HMP)

# Functions ----
extract_taxa = function(dataset){
  # Filter at a given taxonomic rank
  genus_table = subset_taxa(dataset, !is.na(Genus) & is.na(Species))
  # Extract the desired phylotypes
  # Non-western
  Prevotella = otu_table(subset_taxa(genus_table, Genus=="Prevotella"))
  Treponema = otu_table(subset_taxa(genus_table, Genus=="Treponema"))
  # Western
  Bacteroides = otu_table(subset_taxa(genus_table, Genus=="Bacteroides"))
  Bifidobacterium = otu_table(subset_taxa(genus_table, Genus=="Bifidobacterium"))
  Barnesiella = otu_table(subset_taxa(genus_table, Genus=="Barnesiella"))
  Cytophagales = otu_table(subset_taxa(genus_table, Order=="Cytophagales"))
  
  marker_taxa = data.frame(Prevotella = as.vector(Prevotella), 
                           Treponema = as.vector(Treponema),
                           Bacteroides = as.vector(Bacteroides), 
                           Bifidobacterium = as.vector(Bifidobacterium),
                           Barnesiella = as.vector(Barnesiella), 
                           Cytophagales = as.vector(Cytophagales),
                           sample_data(dataset)[,"non_westernized"], 
                           sample_data(dataset)[,"country"])
  
  marker_taxa$non_westernized = as.factor(marker_taxa$non_westernized)
  marker_taxa$country = as.factor(marker_taxa$country)
  marker_taxa
}

col_CAG_markers = function(dataset){
  # Filter at a given taxonomic rank
  genus_table = subset_taxa(dataset, !is.na(Genus) & is.na(Species))
  # Extract the desired phylotypes
  # Non-western
  Prevotella = otu_table(subset_taxa(genus_table, Genus=="Prevotella"))
  Faecalibacterium = otu_table(subset_taxa(genus_table, Genus=="Faecalibacterium"))
  Escherichia = otu_table(subset_taxa(genus_table, Genus=="Escherichia"))
  Akkermansia = otu_table(subset_taxa(genus_table, Genus=="Akkermansia"))
  Methanobrevibacter = otu_table(subset_taxa(genus_table, Genus=="Methanobrevibacter"))
  
  marker_taxa = data.frame(Prevotella = as.vector(Prevotella), 
                           Faecalibacterium = as.vector(Faecalibacterium), 
                           Escherichia = as.vector(Escherichia),
                           Akkermansia = as.vector(Akkermansia), 
                           Methanobrevibacter = as.vector(Methanobrevibacter))
  marker_taxa
}

determine_CAGs = function(x){
  # Reduce the dataset to most abundant phylotypes, compute Spearman's correlation,
  # and cluster the phylotypes using hierarchical clustering with Ward's method.
  # The output can be ploted (dendrogram or heatmap)
  
  # Genera with median abundance >= 0.01%
  g_table = otu_table(subset_taxa(x, Kingdom!="Viruses" & !is.na(Genus) & is.na(Species)))
  genus_median_abund = apply(g_table, MARGIN = 1, FUN = median)
  abundant_genera = g_table[genus_median_abund >= 0.01, ]
  
  # Correlation matrix
  spearman_genus_matrix = cor(t(abundant_genera), method = "spearman")
  spearman_genus_tree = hclust(dist(spearman_genus_matrix), method = "ward")
  results = list(CAG_taxa = rownames(abundant_genera), 
                 Spearman_matrix = spearman_genus_matrix, 
                 Dendrogram = spearman_genus_tree)
  results
}

m_sd = function(x){
  m = round(mean(x, na.rm = T), 5)
  s = round(sd(x, na.rm = T), 5)
  paste(m, "±", s)
}

# Analysis ----

# Convert data to Phyloseq format
Asnicar = ExpressionSet2phyloseq(AsnicarF_2017.metaphlan_bugs_list.stool())
Brito = ExpressionSet2phyloseq(BritoIL_2016.metaphlan_bugs_list.stool())
Feng = ExpressionSet2phyloseq(FengQ_2015.metaphlan_bugs_list.stool())
Heitz_Buschart = ExpressionSet2phyloseq(`Heitz-BuschartA_2016.metaphlan_bugs_list.stool`())
HMP = ExpressionSet2phyloseq(HMP_2012.metaphlan_bugs_list.stool())
Karlsson = ExpressionSet2phyloseq(KarlssonFH_2013.metaphlan_bugs_list.stool())
LeChatelier = ExpressionSet2phyloseq(LeChatelierE_2013.metaphlan_bugs_list.stool())
Liu = ExpressionSet2phyloseq(LiuW_2016.metaphlan_bugs_list.stool())
Loman = ExpressionSet2phyloseq(LomanNJ_2013.metaphlan_bugs_list.stool())
Nielsen = ExpressionSet2phyloseq(NielsenHB_2014.metaphlan_bugs_list.stool())
Obregon_Tito = ExpressionSet2phyloseq(`Obregon-TitoAJ_2015.metaphlan_bugs_list.stool`())
QinJ = ExpressionSet2phyloseq(QinJ_2012.metaphlan_bugs_list.stool())
QinN = ExpressionSet2phyloseq(QinN_2014.metaphlan_bugs_list.stool())
Rampelli = ExpressionSet2phyloseq(RampelliS_2015.metaphlan_bugs_list.stool())
Raymond = ExpressionSet2phyloseq(RaymondF_2016.metaphlan_bugs_list.stool())
Schirmer = ExpressionSet2phyloseq(SchirmerM_2016.metaphlan_bugs_list.stool())
Vatanen = ExpressionSet2phyloseq(VatanenT_2016.metaphlan_bugs_list.stool())
Vincent = ExpressionSet2phyloseq(VincentC_2016.metaphlan_bugs_list.stool())
Vogtmann = ExpressionSet2phyloseq(VogtmannE_2016.metaphlan_bugs_list.stool())
Xie = ExpressionSet2phyloseq(XieH_2016.metaphlan_bugs_list.stool())
Yu = ExpressionSet2phyloseq(YuJ_2015.metaphlan_bugs_list.stool())
Zeller = ExpressionSet2phyloseq(ZellerG_2014.metaphlan_bugs_list.stool())

# Remove diseased individuals, antibiotic consumers and non-adults

Asnicar_clean  = subset_samples(Asnicar,
                                disease == 'healthy' & 
                                  is.na(antibiotics_current_use) & 
                                  age_category=='adult')

Brito_clean  = subset_samples(Brito,
                              is.na(disease) &
                                is.na(antibiotics_current_use) &
                                age_category=='adult')

Feng_clean  = subset_samples(Feng,
                             is.na(disease) & 
                               is.na(antibiotics_current_use) & 
                               age_category=='adult')

Heitz_Buschart_clean  = subset_samples(Heitz_Buschart,
                                       is.na(disease) & 
                                         antibiotics_current_use=='no' & 
                                         age_category=='adult')

HMP_clean  = subset_samples(HMP,
                            disease=='healthy' & 
                              is.na(antibiotics_current_use) & 
                              age_category=='adult')

LeChatelier_clean  = subset_samples(LeChatelier,
                                    is.na(disease) & 
                                      antibiotics_current_use=='no' & 
                                      age_category=='adult')
Liu_nW_clean  = subset_samples(Liu,
                               non_westernized=='yes' & 
                                 is.na(disease) & 
                                 is.na(antibiotics_current_use) &
                                 age_category=='adult')

Liu_W_clean  = subset_samples(Liu,
                              non_westernized=='no' & 
                                is.na(disease) & 
                                is.na(antibiotics_current_use) & 
                                age_category=='adult')

Nielsen_clean  = subset_samples(Nielsen,
                                disease=='healthy' & 
                                  is.na(antibiotics_current_use) &
                                  age_category=='adult')

Obregon_Tito_W_clean  = subset_samples(Obregon_Tito,
                                       non_westernized=='no' & 
                                         disease=='healthy' & 
                                         is.na(antibiotics_current_use) & 
                                         age_category=='adult')

Obregon_Tito_nW_clean  = subset_samples(Obregon_Tito,
                                        non_westernized=='yes' & 
                                          disease=='healthy' & 
                                          is.na(antibiotics_current_use) & 
                                          age_category=='adult')
QinJ_clean  = subset_samples(QinJ,
                             is.na(disease) &
                               is.na(antibiotics_current_use) &
                               age_category=='adult')

QinN_clean  = subset_samples(QinN,
                             is.na(disease) &
                               antibiotics_current_use=='no' &
                               age_category=='adult')

Rampelli_clean  = subset_samples(Rampelli,
                                 is.na(disease) &
                                   is.na(antibiotics_current_use) &
                                   age_category=='adult')

Raymond_clean  = subset_samples(Raymond,
                                is.na(disease) &
                                  antibiotics_current_use=='no' &
                                  age_category=='adult')

Schirmer_clean  = subset_samples(Schirmer,
                                 disease=='healthy' &
                                   antibiotics_current_use=='no' &
                                   age_category=='adult')

Vogtmann_clean  = subset_samples(Vogtmann,
                                 is.na(disease) & 
                                   is.na(antibiotics_current_use) & 
                                   age_category=='adult')

Xie_clean  = subset_samples(Xie,
                            is.na(disease) & 
                              is.na(antibiotics_current_use) &
                              age_category=='adult')

Zeller_clean  = subset_samples(Zeller,
                               is.na(disease) & 
                                 is.na(antibiotics_current_use) & 
                                 age_category=='adult')

# Merge all clean datasets

merged_datasets = merge_phyloseq(Asnicar_clean, Brito_clean, Feng_clean,
                                 Heitz_Buschart_clean, HMP_clean, 
                                 Liu_nW_clean, Liu_W_clean, Nielsen_clean,
                                 Obregon_Tito_W_clean, Obregon_Tito_nW_clean,
                                 QinJ_clean, QinN_clean, Rampelli_clean,
                                 Raymond_clean, Schirmer_clean, Vogtmann_clean,
                                 Xie_clean, Zeller_clean)


# Analyses ----
markers_merged = extract_taxa(merged_datasets)

# Tables
table(markers_merged$country,markers_merged$non_westernized)

# non-Westernized table
table(markers_merged$non_westernized)
# Country Table
table(markers_merged$country)

# Prevotella - Bacteroides coexclusion
# Overall
with(markers_merged, cor.test(Prevotella, Bacteroides, method = "s"))

# by Westernization status
dlply(markers_merged, .variables = "non_westernized", 
      .fun =  function(x) cor.test(x$Prevotella, x$Bacteroides, method = "s"))

# Biomarker abundance
Prevotella_mean_w = aggregate(Prevotella ~ non_westernized, 
                              data = markers_merged, m_sd)
Treponema_mean_w = aggregate(Treponema ~ non_westernized, 
                             data = markers_merged, m_sd)
Bacteroides_mean_w = aggregate(Bacteroides ~ non_westernized, 
                               data = markers_merged, m_sd)
Bifidobacterium_mean_w = aggregate(Bifidobacterium ~ non_westernized, 
                                   data = markers_merged, m_sd)
Barnesiella_mean_w = aggregate(Barnesiella ~ non_westernized, 
                               data = markers_merged, m_sd)
Cytophagales_mean_w = aggregate(Cytophagales ~ non_westernized, 
                                data = markers_merged, m_sd)

# By country
Prevotella_mean_c = aggregate(Prevotella ~ country, 
                              data = markers_merged, m_sd)
Treponema_mean_c = aggregate(Treponema ~ country, 
                             data = markers_merged, m_sd)
Bacteroides_mean_c = aggregate(Bacteroides ~ country,
                               data = markers_merged, m_sd)
Bifidobacterium_mean_c = aggregate(Bifidobacterium ~ country, 
                                   data = markers_merged, m_sd)
Barnesiella_mean_c = aggregate(Barnesiella ~ country, 
                               data = markers_merged, m_sd)
Cytophagales_mean_c = aggregate(Cytophagales ~ country, 
                                data = markers_merged, m_sd)


# By westernization + country 
Prevotella_mean_cw = aggregate(Prevotella ~ non_westernized + country, 
                               data = markers_merged, m_sd)
Treponema_mean_cw = aggregate(Treponema ~ non_westernized + country, 
                              data = markers_merged, m_sd)
Bacteroides_mean_cw = aggregate(Bacteroides ~ non_westernized + country, 
                                data = markers_merged, m_sd)
Bifidobacterium_mean_cw = aggregate(Bifidobacterium ~ non_westernized + country,
                                    data = markers_merged, m_sd)
Barnesiella_mean_cw = aggregate(Barnesiella ~ non_westernized + country, 
                                data = markers_merged, m_sd)
Cytophagales_mean_cw = aggregate(Cytophagales ~ non_westernized + country, 
                                 data = markers_merged, m_sd)

##### Heatmap ----
Prevotella_heatmap_cw = aggregate(Prevotella ~ non_westernized + country, 
                                  data = markers_merged, mean)
Treponema_heatmap_cw = aggregate(Treponema ~ non_westernized + country, 
                                 data = markers_merged, mean)
Bacteroides_heatmap_cw = aggregate(Bacteroides ~ non_westernized + country, 
                                   data = markers_merged, mean)
Bifidobacterium_heatmap_cw = aggregate(Bifidobacterium ~ non_westernized + country, 
                                       data = markers_merged, mean)
Barnesiella_heatmap_cw = aggregate(Barnesiella ~ non_westernized + country, 
                                   data = markers_merged, mean)
Cytophagales_heatmap_cw = aggregate(Cytophagales ~ non_westernized + country, 
                                    data = markers_merged, mean)

heatmap_matrix = data.frame(Prevotella = Prevotella_heatmap_cw$Prevotella,
                            Treponema =  Treponema_heatmap_cw$Treponema,
                            Bacteroides = Bacteroides_heatmap_cw$Bacteroides, 
                            Bifidobacterium = Bifidobacterium_heatmap_cw$Bifidobacterium,
                            Barnesiella = Barnesiella_heatmap_cw$Barnesiella,
                            Cytophagales = Cytophagales_heatmap_cw$Cytophagales)

heatmap(as.matrix(heatmap_matrix), Rowv = NA, Colv = NA, keep.dendro = F)

# Co-abundance groups (CAGs) ----
dataset_list = list(AUT = subset_samples(merged_datasets, country=="AUT"), 
                    CAN = subset_samples(merged_datasets, country=="CAN"), 
                    CHN = subset_samples(merged_datasets, country=="CHN"), 
                    DEU = subset_samples(merged_datasets, country=="DEU"), 
                    DNK = subset_samples(merged_datasets, country=="DNK"), 
                    ESP = subset_samples(merged_datasets, country=="ESP"), 
                    FJI = subset_samples(merged_datasets, country=="FJI"), 
                    FRA = subset_samples(merged_datasets, country=="FRA"), 
                    GBR = subset_samples(merged_datasets, country=="GBR"), 
                    ITA = subset_samples(merged_datasets, country=="ITA"), 
                    LUX = subset_samples(merged_datasets, country=="LUX"), 
                    MNG_w = subset_samples(merged_datasets, country=="MNG" & non_westernized=="no"),
                    MNG_nw = subset_samples(merged_datasets, country=="MNG" & non_westernized=="yes"), 
                    NLD = subset_samples(merged_datasets, country=="NLD"), 
                    PER = subset_samples(merged_datasets, country=="PER"), 
                    TZA = subset_samples(merged_datasets, country=="TZA"), 
                    USA = subset_samples(merged_datasets, country=="USA"))

dataset_list

# CAGs of the different populations
CAG_all_countries = lapply(dataset_list, FUN = determine_CAGs)

plot_CAGs = function(x){
  name = deparse(substitute(x))
  #pdf(paste(name, ".pdf", sep = ""))
  CAG_heatmap = aheatmap(x$Spearman_matrix, color = "-Spectral:100", 
                         scale = "none", breaks = NA, Rowv=T, Colv=T, cexRow=1, 
                         hclustfun = "ward", distfun = "euclidean", main = name)
  #dev.off()
}

# Plotting each heatmap
plot_CAGs(CAG_all_countries$AUT)
plot_CAGs(CAG_all_countries$CAN)
plot_CAGs(CAG_all_countries$CHN)
plot_CAGs(CAG_all_countries$DEU)
plot_CAGs(CAG_all_countries$DNK)
plot_CAGs(CAG_all_countries$ESP)
plot_CAGs(CAG_all_countries$FJI)
plot_CAGs(CAG_all_countries$FRA)
plot_CAGs(CAG_all_countries$GBR)
plot_CAGs(CAG_all_countries$ITA)
plot_CAGs(CAG_all_countries$LUX)
plot_CAGs(CAG_all_countries$MNG_w)
plot_CAGs(CAG_all_countries$MNG_nw)
plot_CAGs(CAG_all_countries$NLD)
plot_CAGs(CAG_all_countries$PER)
plot_CAGs(CAG_all_countries$TZA)
plot_CAGs(CAG_all_countries$USA)

# CAG membership
# Create a plot with presence/absence of all the genera with median abundance >= 0.01% in all datasets
CAG_list = lapply(CAG_all_countries, function(x) x$CAG_taxa)
CAG_vector = unique(as.vector(unlist(CAG_list)))
table(as.vector(unlist(CAG_list)))
membership_table = as.data.frame(lapply(CAG_all_countries, 
                                        FUN = function(y) as.numeric(sapply(CAG_vector, 
                                                                            FUN = function(x) x%in%y$CAG_taxa))))
rownames(membership_table) = unique(as.vector(unlist(CAG_list)))
membership_table = membership_table[order(-rowSums(membership_table)),]
aheatmap(membership_table, scale = "none", border_color = 'white', 
         breaks = NA, Rowv=F, Colv=F, cexRow=1)

# Correlation of Colombian CAG markers in other populations
cor_matrices_country = lapply(dataset_list, FUN = function(x) cor(col_CAG_markers(x)))
