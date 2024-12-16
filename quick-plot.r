library(Rtsne)
source("~/apps/FastPCA.R")
library(gplots)
library(viridis)
library(vegan)
library(readr)


#### fairy coverages
setwd("/data/project_CMO/flye_qc18_dec14/")

# import fairy bin coverages
fairy.bin.in = read.table("metabat_bins/metabat_graph_cov.tsv",head=T)
fairy.bin.cov = fairy.bin.in[,grepl(".gz$",colnames(fairy.bin.in))]
fairy.bin.met = fairy.bin.in[,1:3]
rownames(fairy.bin.cov) = fairy.bin.in[,1]
rownames(fairy.bin.met) = fairy.bin.in[,1]

# import fairy noQC bin coverages
fairy.binq9.in = read.table("metabat_bins/metabat_graph_noqc_cov.tsv",head=T)
fairy.binq9.cov = fairy.binq9.in[,grepl(".gz$",colnames(fairy.binq9.in))]
fairy.binq9.met = fairy.binq9.in[,1:3]
rownames(fairy.binq9.cov) = fairy.binq9.in[,1]
rownames(fairy.binq9.met) = fairy.binq9.in[,1]

# fairy.bin.pca = FastPCA(as.matrix(fairy.bin.cov),top.k = 50,center.use = FALSE, scale.use = FALSE)
# fairy.bin.tsne = Rtsne(fairy.bin.pca$x,dims = 2,check_duplicates = F,pca = F,verbose = T,max_iter = 2500)
# 
# plot(fairy.bin.tsne$Y,pch='.',cex=2)


# import fairy assembly coverages
fairy.asm.in = read.table("assembly_graph_cov.tsv",head=T)
fairy.asm.cov = fairy.asm.in[,grepl(".gz$",colnames(fairy.asm.in))]
fairy.asm.met = fairy.asm.in[,1:3]
rownames(fairy.asm.cov) = fairy.asm.in[,1]
rownames(fairy.asm.met) = fairy.asm.in[,1]

# fairy.asm.pca = FastPCA(as.matrix(fairy.asm.cov),top.k = 50,center.use = FALSE, scale.use = FALSE)
# fairy.asm.tsne = Rtsne(fairy.asm.pca$x,dims = 2,check_duplicates = F,pca = F,verbose = T,max_iter = 2500)
# 
# plot(fairy.asm.tsne$Y,)



# import fairy reads coverages
fairy.read.in = read.table("reads_cov.tsv",head=T)
fairy.read.cov = fairy.read.in[,grepl(".gz$",colnames(fairy.read.in))]
fairy.read.met = fairy.read.in[,1:3]
rownames(fairy.read.cov) = fairy.read.in[,1]
rownames(fairy.read.met) = fairy.read.in[,1]

# fairy.read.pca = FastPCA(as.matrix(fairy.read.cov),top.k = 50,center.use = FALSE, scale.use = FALSE)
# fairy.read.tsne = Rtsne(fairy.read.pca$x,dims = 2,check_duplicates = F,pca = F,verbose = T,max_iter = 2500)
# 
# plot(fairy.read.tsne$Y,)


# import TNFs
tnfs = as.character(read.table("tnfs.txt",stringsAsFactors=F))
tnfs = tnfs[2:length(tnfs)]

lrn.bin = data.frame(read_table("Tetra_metabat_graph_out_all_1500.lrn",comment="%",col_names=F))
lrn.bin$X1 = NULL
lrn.bin.id = read.table("Tetra_metabat_graph_out_all_1500.names",comment="%",row.names=1,stringsAsFactors=F)
rownames(lrn.bin) = lrn.bin.id$V2
colnames(lrn.bin) = tnfs

lrn.asm = data.frame(read_table("Tetra_all_positive_graph_nodes_1500.lrn",comment="%",col_names=F))
lrn.asm$X1 = NULL
lrn.asm.id = read.table("Tetra_all_positive_graph_nodes_1500.names",comment="%",row.names=1,stringsAsFactors=F)
rownames(lrn.asm) = lrn.asm.id$V2
colnames(lrn.asm) = tnfs

lrn.read = data.frame(read_table("/data/scratch/fairy_qc20_out/Tetra_63_1500.lrn",comment="%",col_names=F))
lrn.read$X1 = NULL
lrn.read.id = read.table("/data/scratch/fairy_qc20_out/Tetra_63_1500.names",comment="%",row.names=1,stringsAsFactors=F)
rownames(lrn.read) = lrn.read.id$V2
colnames(lrn.read) = tnfs

### combine 

#bin = MAGs
#asm = full assembly
#read = all reads

all.lrn.cov = rbind(
  cbind(lrn.bin, fairy.bin.cov),
  cbind(lrn.asm[intersect(rownames(lrn.asm),rownames(fairy.asm.cov)),], 
      fairy.asm.cov[intersect(rownames(lrn.asm),rownames(fairy.asm.cov)),]),
  cbind(lrn.read, fairy.read.cov)
  )
# 

#all.lrn.cov = cbind(lrn.asm[intersect(rownames(lrn.asm),rownames(fairy.asm.cov)),], fairy.asm.cov[intersect(rownames(lrn.asm),rownames(fairy.asm.cov)),])
all.lrn.cov = all.lrn.cov[rowSums(all.lrn.cov)>0,colSums(all.lrn.cov)>0]

# run PCA
#all.lrn.cov.pca.bak = all.lrn.cov.pca
all.lrn.cov.pca = FastPCA(all.lrn.cov,top.k = 50,scale.use = T, center.use = T)

#This adds a column of color values
# based on the y values
rbPal <- colorRampPalette(c('red','blue'))
#GCcol <- rbPal(100)[as.numeric(cut(GCpct,breaks = 100))] #GC coloring
ptcol <- rbPal(100)[as.numeric(cut(all.lrn.cov.pca$x[,1],breaks = 100))] #PC1 coloring
names(ptcol) = rownames(all.lrn.cov)
#ptcol <- rbPal(100)[as.numeric(cut(lrn.pca$x[,2],breaks = 100))]
met.all = rbind(fairy.bin.met, fairy.asm.met, fairy.read.met)

plot(all.lrn.cov.pca$x[,1:2],pch='.',col=ptcol, cex=2)
plot(met.all[rownames(all.lrn.cov),"totalAvgDepth"],met.all[rownames(all.lrn.cov),"contigLen"],log='xy',pch='.')

#mysam=sample(1:nrow(all.lrn.cov),20000)
mysam=(1:nrow(all.lrn.cov))[met.all[rownames(all.lrn.cov),"totalAvgDepth"] > 1e-2 & met.all[rownames(all.lrn.cov),"contigLen"] > 5000]
x = all.lrn.cov.pca$x[mysam,]
fairy.bin.tsne = Rtsne(x^.25,dims = 2, check_duplicates = F, pca = F, verbose = T, max_iter = 1000)
plot(fairy.bin.tsne$Y,pch=21,bg=ptcol[mysam],col=NULL, cex=0.5)
plot(fairy.bin.tsne$Y,pch=21,bg=ptcol[mysam],col=NULL, cex=rescale(fairy.asm.met[rownames(all.lrn.cov),"contigLen"]^.25,to=c(0,4)))
library(scales)

rownames(fairy.bin.tsne$Y) = rownames(all.lrn.cov)[mysam]
ptcol2 = ptcol[rownames(x)]
names =sub(".fastq.gz_qc18.fq.gz","",colnames(fairy.asm.cov))

hmmhits = read_table("all_positive_graph_nodes_cant_hyd_cut_tc.tblout",comment = "#",col_names = FALSE)
hmmhits$contigs = sub("(_[^_]+){3}$", "", hmmhits$X1)
intersect(hmmhits$contigs, rownames(x))

pdf(file="sample-tsne.pdf",width=12,height=12)
for(i in 1:ncol(fairy.asm.cov)) {
cexsize = (1:nrow(x))*NA
pick.pts = fairy.asm.cov[rownames(x),i]>0
cexsize[pick.pts] = rescale(fairy.asm.met[rownames(all.lrn.cov),"contigLen"]^.25,to=c(0,4))[pick.pts]

plot(fairy.bin.tsne$Y[,1:2],pch=21,bg="#33333333",col=NULL,cex=0.1, main=nano.meta[names[i],"uniqid"])
points(fairy.bin.tsne$Y[,1:2],pch=21,bg=ptcol2,col=NULL, cex=cexsize)
points(fairy.bin.tsne$Y[hmmhits$contigs,],pch=as.numeric(as.factor(hmmhits$X3)))
}
dev.off()


# all.bin = rbind(cbind(lrn.bin, fairy.bin.cov))
# all.bin.pca = FastPCA(all.bin,top.k = 50, scale.use = T, center.use = T)
 



bin.ids = read.table("metabat_graph_out_contigids.tsv",sep="\t",head=F)
bincol = rainbow(length(levels(as.factor(bin.ids$V1))))
nosplitids = sub("_[^_]*$", "", rownames(fairy.binq9.cov))
nosplitbinids = sub("_[^_]*$", "", rownames(fairy.binq9.cov))

# mycol = match(nosplitids,bin.ids$V2)[mysam]
# plot(fairy.bin.tsne$Y,pch=21,bg=mycol,col=mycol, cex=0.5)


split.agg = aggregate(fairy.binq9.cov,by=list(nosplitbinids), FUN="mean")
rownames(split.agg) = split.agg$Group.1
split.agg$Group.1 = NULL

binlist = bin.ids$V1[match(rownames(split.agg),bin.ids$V2)]

bin.agg = aggregate(split.agg, by=list(binlist), FUN="mean")
rownames(bin.agg) = bin.agg$Group.1
bin.agg$Group.1 = NULL
colnames(bin.agg) = sub(".fastq.gz","",colnames(bin.agg))

nano.meta = read.table("/data/scratch/fairy_qc20_out/nanopore_metadata_dec14.tsv",head=T, sep="\t")
rownames(nano.meta) = make.names(nano.meta$runbarcode)

tophits = read.table("metabat_bins/tophits.sk.tsv",sep="\t",head=F)
rownames(tophits) = sub(".sk.tsv.*","",tophits$V1)

bin.agg = bin.agg[rowSums(bin.agg)>0, colSums(bin.agg)>0]

pdf(file="hm.pdf",width=24,height=24)
hm = heatmap.2(as.matrix(bin.agg)^.5,scale="none",trace="none", mar=c(20,40),col=viridis(100),
          distfun = function(x) vegdist(x, method = "bray"),
          hclustfun = function(x) hclust(x, method = "ward"),
          labCol = paste0(nano.meta[colnames(bin.agg),"uniqid"],"|",nano.meta[colnames(bin.agg),"pcr"]),
#          labRow = tophits[rownames(bin.agg),"V2"],
          cexRow = 0.5, cexCol = 1)
dev.off()

#remove enrichment contam bins
#cut2 = cutree(as.hclust(hm$rowDendrogram),k=2)
bin.nopseudo = bin.agg[! rownames(bin.agg) %in% c("metabat_graph_out.53.fa",
                         "metabat_graph_out.43.fa",
                         "metabat_graph_out.74.fa",
                         "metabat_graph_out.32.fa"),]

pdf(file="hm.pdf",width=32,height=24)
hm = heatmap.2(as.matrix(bin.nopseudo)^.25,scale="none",trace="none", mar=c(20,40),col=viridis(100),
               distfun = function(x) vegdist(x, method = "bray"),
               hclustfun = function(x) hclust(x, method = "ward"),
               labCol = paste0(nano.meta[colnames(bin.agg),"uniqid"],"|",nano.meta[colnames(bin.agg),"pcr"]),
               labRow = tophits[rownames(bin.agg),"V2"],
               cexRow = 0.5, cexCol = 1)
dev.off()

bin.nopseudo.pca = FastPCA(t(bin.nopseudo)^.25,top.k = 20)
bin.nopseudo.tsne = Rtsne(bin.nopseudo.pca$x,perplexity = 10, check_duplicates = F, pca = F, verbose = T, max_iter = 2500)

nclus=8
clus = cutree(as.hclust(hm$colDendrogram),k=nclus)

pdf(file="hm.pdf",width=32,height=24)
hm = heatmap.2(as.matrix(bin.nopseudo)^.25,scale="none",trace="none", mar=c(20,40),col=viridis(100),
               distfun = function(x) vegdist(x, method = "bray"),
               hclustfun = function(x) hclust(x, method = "ward"),
               labCol = paste0(nano.meta[colnames(bin.nopseudo),"uniqid"],"|",nano.meta[colnames(bin.nopseudo),"pcr"]),
               # labRow = tophits[rownames(bin.agg),"V2"],
               cexRow = 0.5, cexCol = 1,
               ColSideColors = rainbow(nclus)[clus])

#plot(bin.nopseudo.tsne$Y,pch=2+as.numeric(as.factor(nano.meta[colnames(bin.nopseudo),"pool"])),col=rainbow(nclus)[clus],cex=6)
plot(bin.nopseudo.tsne$Y,pch=21,bg=rainbow(nclus)[clus],cex=6)
text(bin.nopseudo.tsne$Y, labels=paste0(nano.meta[colnames(bin.nopseudo),"uniqid"],"|",nano.meta[colnames(bin.nopseudo),"pcr"]),cex=1)
plot(bin.nopseudo.tsne$Y,pch=21,bg=rainbow(nclus)[clus],cex=6)
text(bin.nopseudo.tsne$Y, labels=paste0(nano.meta[colnames(bin.nopseudo),"runbarcode"],"|",nano.meta[colnames(bin.nopseudo),"pcr"]),cex=1)
dev.off()


library(tidyr)
library(ggplot2)
library(dplyr)

# Transform the data frame to long format
df = data.frame(t(bin.nopseudo))
#rownames(df) = make.names(nano.meta[rownames(df),"uniqid"], unique = TRUE)
colnames(df) = make.names(tophits[colnames(df),"V2"])

# vec_string <- paste(tophits$V2, collapse = "\n")
# taxa <- data.frame(read_delim(vec_string, delim = ";", col_names = c("do","ph","cl","or","fa","ge","sp")))
# taxa$bin = tophits$V1
# rownames(taxa) = taxa$bin

# Convert the vector into a data frame with one column
# taxa <- data.frame(data = tophits$V2, stringsAsFactors = FALSE)
# 
# # Convert the vector into a data frame
# df <- data.frame(data = tophits$V2, stringsAsFactors = FALSE)
# 
# # Parse the vector into a data frame with missing taxa handled
# df_clean <- df %>%
#   separate_rows(data, sep = ";") %>%                        # Split each entry into key-value pairs
#   separate(data, into = c("key", "value"), sep = ":", fill = "right", extra="merge") %>% # Split keys and values
#   pivot_wider(names_from = key, values_from = value, values_fn = list)       # Reshape into a wide format
# 
df$sample = rownames(df)
# 
# df_long <- df %>%
#   pivot_longer(
#     cols = -sample,            # All columns except 'sample' should be "gathered"
#     names_to = "species",      # New column name for former column headings
#     values_to = "abundance"    # New column name for values
#   )

# # Arrange the order of x-axis (sample)
# neworder = make.names(colnames(bin.nopseudo)[hm$colInd],unique=T)
# df_long$sample <- factor(df_long$sample, levels = neworder)
# 
# ggplot(df_long, aes(x = sample, y = abundance, fill = species)) +
#   geom_bar(stat = "identity", position = "stack") +
#   labs(title = "Stacked Bar Chart of Abundances", x = "Sample", y = "Abundance") +
#   theme_minimal() +
#   guides(fill = guide_legend(ncol = 1)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# 
# ggsave(file="stacked.pdf",width=50,height=50, limitsize=FALSE)
# 
# 
# df_long <- df_long %>%
#   group_by(sample) %>%
#   mutate(percent_abundance = (abundance / sum(abundance)) * 100)
# 
# # Step 3: Plot the data as a percentage stacked bar chart
# ggplot(df_long, aes(x = sample, y = percent_abundance, fill = species)) +
#   geom_bar(stat = "identity", position = "stack") +
#   labs(
#     title = "Percentage Stacked Bar Chart of Abundances",
#     x = "Sample",
#     y = "Percentage Abundance (%)"
#   ) +
#   theme_minimal() +
#   guides(fill = guide_legend(ncol = 1)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# 
# # Use rainbow palette
# n_species <- length(unique(df_long$species)) # Number of species
# rainbow_colors <- rainbow(n_species)         # Generate colors
# 
# # Plot with custom palette
# ggplot(df_long, aes(x = sample, y = percent_abundance, fill = species)) +
#   geom_bar(stat = "identity", position = "stack") +
#   labs(
#     title = "Percentage Stacked Bar Chart of Abundances",
#     x = "Sample",
#     y = "Percentage Abundance (%)"
#   ) +
#   theme_minimal() +
#   guides(fill = guide_legend(ncol = 1)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#   scale_fill_manual(values = rainbow_colors)


df_long <- df %>%
  pivot_longer(cols = -sample, names_to = "species", values_to = "abundance") %>%
  group_by(sample) %>%
  mutate(percent_abundance = (abundance / sum(abundance)) * 100)

# Ensure species order is consistent
df_long$species <- factor(df_long$species)
# df_long$genus = factor(taxa[df_long$taxa,"ge"])

# Arrange the order of x-axis (sample)
neworder = make.names(colnames(bin.nopseudo)[hm$colInd],unique=T)
df_long$sample <- factor(df_long$sample, levels = neworder)


# Create a custom color palette
custom_colors <- setNames(
  rainbow(length(levels(df_long$species))), # Default rainbow colors
  levels(df_long$species)                   # Match species order
)

# Assign black to a specific species (e.g., "speciesB")
custom_colors["NA."] <- "black"

# Plot with the fixed custom palette
ggplot(df_long, aes(x = sample, y = percent_abundance, fill = species)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Percentage Stacked Bar Chart of Abundances",
    x = "Sample",
    y = "Percentage Abundance (%)"
  ) +
  theme_minimal() +
  guides(fill = guide_legend(ncol = 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(values = custom_colors)


ggsave(file="stacked_ra.pdf",width=50,height=50, limitsize=FALSE)



