rm(list = ls())

library(spider)
library(ape)
library(readxl)
library(stringr)
library(tidyr)
library(ggplot2)
library(gghalves)

setwd("/Users/maximilianwagner/Library/Mobile Documents/com~apple~CloudDocs/Gouania/Gouania_phylogenomics/Gouania_Species_delimitation/DNA_Barcoding_gap/")

alignment_raw <- read_xlsx("COI_Gouania_09092020_trimmed.nexus.xlsx")
metadata <- read_xlsx("Gouania_meta_07122021.xlsx")
metadata_spec <- metadata[c("ID2", "Species_short")]

alignment_raw_spec <- merge(metadata_spec, alignment_raw,by = "ID2")
#alignment_raw_spec$ID <- str_c(alignment_raw_spec$Species_short,"_", alignment_raw_spec$ID2)
alignment_raw_spec <- alignment_raw_spec[c("Species_short", "sequence")]

alignment_done <-  t(sapply(strsplit(alignment_raw_spec[,2],""), tolower))

rownames(alignment_done) <- alignment_raw_spec[,1]

alignment_done <- as.DNAbin(alignment_done)


pairwise_dist <- dist.dna(alignment_done)


IDs <-labels(alignment_done)

spec_names <- sapply(strsplit(dimnames(alignment_done)[[1]], split="_"),function(x) paste(x[1],x[2], sep="_"))


inter <- nonConDist(pairwise_dist, IDs); inter
intra <- maxInDist(pairwise_dist, IDs); intra


intra_vs_inter <- as.data.frame(cbind(IDs,inter, intra))

intra_vs_inter$species <- str_split_fixed(intra_vs_inter$IDs, "_", 2)

intra_vs_inter_long <- pivot_longer(intra_vs_inter,
                                    cols=2:3, names_to = "intra_inter", 
                                    values_to = "value")


intra_vs_inter_long <- as.data.frame(as.matrix(intra_vs_inter_long))
intra_vs_inter_long$value <- as.numeric(as.character(intra_vs_inter_long$value))

intra_vs_inter_plot  <- ggplot(intra_vs_inter_long, aes(intra_inter, value, color = species.1)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~species.1, scales = "free"); intra_vs_inter_plot
