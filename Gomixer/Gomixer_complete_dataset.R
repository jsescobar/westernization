# Piellou's J - All dataset
Piellou = diversityresult(x = microbio.rare, index = "Jevenness", method = "each site")

# Gomixer all dataset.
# Load Gomixer Modules
microbio.gomixer = t(read.table(file = "./files/microbio_selected.gomixer", header = T, sep = "\t", dec = ".", row.names = 1))
microbio.gomixer = microbio.gomixer[order(row.names(microbio.gomixer)),]

# Gomixer modules set A - LPS biosynthesis, mucus degradation and methanogenesis
LPS_biosynthesis = "NA_MF0092"
mucus_degradation = "NA_MF0102"
methanogenesis = c("NA_MF0097", "NA_MF0098", "NA_MF0099")

# Gomixer modules set B - acetate, butyrate and propionate metabolism
acetate_metabolism = c("NA_MF0112", "NA_MF0113")
butyrate_metabolism = c("NA_MF0114", "NA_MF0115", "NA_MF0116", "NA_MF0117")
propionate_metabolism = c("NA_MF0121", "NA_MF0122", "NA_MF0123", "NA_MF0124", "NA_MF0125", "NA_MF0126")

# Gomixer modules set C - lipid, amino acid and carbohydrate degradation
lipid_degradation = c("NA_MF0106", "NA_MF0107", "NA_MF0108", "NA_MF0109", "NA_MF0110", "NA_MF0111")
aminoacid_degradation = c("NA_MF0007", "NA_MF0008", "NA_MF0009", "NA_MF0010", "NA_MF0011", "NA_MF0012", "NA_MF0013", "NA_MF0014", "NA_MF0015", "NA_MF0016", "NA_MF0017", "NA_MF0018", "NA_MF0019", "NA_MF0020", "NA_MF0021", "NA_MF0022", "NA_MF0023", "NA_MF0024", "NA_MF0026", "NA_MF0027", "NA_MF0028", "NA_MF0029", "NA_MF0030", "NA_MF0031", "NA_MF0032", "NA_MF0033", "NA_MF0034", "NA_MF0035", "NA_MF0036", "NA_MF0037", "NA_MF0038", "NA_MF0039", "NA_MF0040", "NA_MF0041", "NA_MF0042", "NA_MF0043")
carbohydrate_degradation = c("NA_MF0045", "NA_MF0046", "NA_MF0047", "NA_MF0048", "NA_MF0049", "NA_MF0050", "NA_MF0051", "NA_MF0052", "NA_MF0053", "NA_MF0054", "NA_MF0055", "NA_MF0056", "NA_MF0057", "NA_MF0058", "NA_MF0059", "NA_MF0060", "NA_MF0061", "NA_MF0062", "NA_MF0063", "NA_MF0064", "NA_MF0065", "NA_MF0066", "NA_MF0067", "NA_MF0068", "NA_MF0069", "NA_MF0070", "NA_MF0071", "NA_MF0072", "NA_MF0073", "NA_MF0074", "NA_MF0075", "NA_MF0076", "NA_MF0077", "NA_MF0078")

# Relative abundance within sets
# Set A
gomx_set_a = data.frame(LPS_biosynthesis = unlist(microbio.gomixer[,LPS_biosynthesis]), mucus_degradation = unlist(microbio.gomixer[,mucus_degradation]), methanogenesis = rowSums(microbio.gomixer[,methanogenesis]), row.names = rownames(microbio.gomixer))
gomx_set_a = gomx_set_a/rowSums(gomx_set_a)


# Set B 
gomx_set_b = data.frame(acetate_metabolism = rowSums(microbio.gomixer[,acetate_metabolism]), butyrate_metabolism = rowSums(microbio.gomixer[,butyrate_metabolism]), propionate_metabolism = rowSums(microbio.gomixer[,propionate_metabolism]))
gomx_set_b = gomx_set_b/rowSums(gomx_set_b)

# Set C
gomx_set_c = data.frame(lipid_degradation = rowSums(microbio.gomixer[,lipid_degradation]), aminoacid_degradation = rowSums(microbio.gomixer[,aminoacid_degradation]), carbohydrate_degradation = rowSums(microbio.gomixer[,carbohydrate_degradation]))
gomx_set_c = gomx_set_c/rowSums(gomx_set_c)



as.vector(a)
# tests
cor_set_a = apply(n5_table[, 4:8], MARGIN = 2, function(x) apply(gomx_set_a, MARGIN = 2, function(y) cor.test(x,y, method = "s")$estimate))

p_set_a = apply(n5_table[, 4:8], MARGIN = 2, function(x) apply(gomx_set_a, MARGIN = 2, function(y) cor.test(x,y, method = "s")$p.value))


cor_set_b = apply(n5_table[, 4:8], MARGIN = 2, function(x) apply(gomx_set_b, MARGIN = 2, function(y) cor.test(x,y, method = "s")$estimate))

p_set_b = apply(n5_table[, 4:8], MARGIN = 2, function(x) apply(gomx_set_b, MARGIN = 2, function(y) cor.test(x,y, method = "s")$p.value))

cor_set_c = apply(n5_table[, 4:8], MARGIN = 2, function(x) apply(gomx_set_c, MARGIN = 2, function(y) cor.test(x,y, method = "s")$estimate))

p_set_c = apply(n5_table[, 4:8], MARGIN = 2, function(x) apply(gomx_set_c, MARGIN = 2, function(y) cor.test(x,y, method = "s")$p.value))
      
fdr_p = qvalue(c(as.vector(p_set_a), as.vector(p_set_b), as.vector(p_set_c)))

qvalue_table = matrix(fdr_p$qvalues, ncol = 5)

