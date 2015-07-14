#  interactome_metabolic_proteins_network_analysis.R
#
#  Copyright 2014 Sebastian Kurscheid <sebastian.kurscheid@anu.edu.au>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  

#---------load libraries--------------
library("biomaRt")
library("gdata")
library("GO.db")
library("KEGGREST")
library("ggplot2")
library("gplots")
library("grid")
library("scales")

#---------custom functions------------
keggConv.batch <- function(x, max = 100, org = "mmu", id.type = "ncbi-geneid") {
  if (max > 100){
    on.exit(print("Maximum number of IDs at a given time is 100"))
  } else {
    x <- paste(id.type, x, sep = ":")
    if (length(x > 100)){
      d1 <- split(x, ceiling(seq_along(x)/max))
      s1 <- lapply(d1, function(y){
         keggConv(org, y)
      })
      return(unlist(s1))
    } else {
      d1 <- split(x, ceiling(seq_along(x)/10))
      s1 <- lapply(d1, function(y){
        keggConv(org, y)
      })
      return(unlist(s1))
    }
  }
} # TODO: edit parameters for function call 

#---------global variables----------------------------------------
# for Fisher's Exact test
# can be "greater", "less", or "two.sided"
alternative = "greater"
p.adjust.method = "fdr"

#---------use ENSEMBL biomaRt for annotation data-----------------
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# list available filters
filters <- listFilters(mouse)
# list available attributes
attribs <- listAttributes(mouse)
pages <- attributePages(mouse)
hsap.attribs <- listAttributes(human)


kegg.brite <- read.xls("/Users/u1001407/Dropbox/REM project-Sebastian/Analysis/KEGG_Brite_Hierarchy.xlsx", sheet = 1, as.is = T)
ids <- unlist(lapply(strsplit(kegg.brite$C, " "), function(x) x[1]))
rownames(kegg.brite) <- ids
save(kegg.brite, file = "data/kegg.brite.rda")
total.keggIDs <- keggLink("mmu", "pathway")
total.keggIDs <- unique(total.keggIDs)
length(unique(total.keggIDs))

#-------------whole cell lysate proteome------------------------------------------------------------------------------------------------
wcl <- read.xls("/Users/u1001407/Dropbox/REM project-Sebastian/HL-1 interactome superset.xlsx", sheet = "WCL", as.is = T)
colnames(wcl) <- c("gene_symbol", "ensembl_gene_id")
entrez_ids <- getBM(attributes = c("ensembl_gene_id","entrezgene"), values = wcl[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)
wcl.human_homologs <- getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"), values = wcl[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)

# remove ensembl_gene_ids which have duplicated entrez_ids
entrez_ids <- entrez_ids[-which(duplicated(entrez_ids$ensembl_gene_id)),]
wcl <- merge(wcl, entrez_ids, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T)
wcl.entrezIDs <- unique(wcl[!is.na(wcl$entrezgene),]$entrezgene)
wcl.keggIDs <- keggConv.batch(wcl.entrezIDs)
wcl.keggQ <- lapply(wcl.keggIDs, function(x) keggGet(x))
wcl.pathways <- unique(unlist(lapply(strsplit(names(unlist(lapply(wcl.keggQ, function(x) x[[1]]$"PATHWAY"))), "\\."), function(x) x[3])))
wcl.pathways.genes <- lapply(wcl.pathways, function(x) keggLink("genes", x))
names(wcl.pathways.genes) <- wcl.pathways
wcl.pathways.genes.entrez_ids <- unique(gsub("mmu:", "", as.character(unlist(wcl.pathways.genes))))
wcl.df <- kegg.brite[gsub("mmu", "", wcl.pathways), ]
wcl.df$ID <- rownames(wcl.df)
wcl.df$total <- rep(0, nrow(wcl.df))
wcl.df$total <- sapply(rownames(wcl.df), function(x) length(wcl.pathways.genes[[paste("mmu", x, sep = "")]]))
wcl.df$count <- rep(0, nrow(wcl.df))
wcl.df$frac <- rep(0, nrow(wcl.df))

for (i in rownames(wcl.df)) {
  kL1 <- keggLink("mmu", paste("mmu", i, sep = ""))
  wcl.df[i, ]$count <- length(which(wcl.keggIDs %in% kL1))
  wcl.df[i, ]$frac <- round(length(which(wcl.keggIDs %in% kL1)) / length(kL1) * 100, 2)
}

# extract list of IDs in pathway
wcl.in_path.IDs <- lapply(rownames(wcl.df), function(x) {
  kL1 <- keggLink("mmu", paste("mmu", x, sep = ""))
  in_path <- wcl.keggIDs[which(wcl.keggIDs %in% kL1)]
})

names(wcl.in_path.IDs) <- rownames(wcl.df)

# perform Fisher's Exact Test for each category
bkgd <- length(unique(total.keggIDs))
smpl <- length(wcl.keggIDs)
ftl <- apply(wcl.df[1,], 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = alternative)
})

wcl.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
wcl.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
wcl.df$ft_fdr <- p.adjust(wcl.df$ft_pval, method = "fdr")

#-----------total interactome----------------------------------------------------------------------------------------------------------
interactome <- read.xls("/Users/u1001407/Dropbox/REM project-Sebastian/HL-1 interactome superset.xlsx", sheet = "sheet 1" , as.is = T)
colnames(interactome)[c(1,2)] <- c("ensembl_gene_id", "gene_symbol")
interactome$GO <- as.factor(interactome$GO)
interactome$RBD <- as.factor(interactome$RBD)

interactome.entrez_ids <- getBM(attributes = c("ensembl_gene_id","entrezgene"), values = interactome[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)
interactome.human_homologs <- getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"), values = interactome[,"ensembl_gene_id"], filters = "ensembl_gene_id", mart = mouse)
interactome.mim <- getBM(attributes = c("ensembl_gene_id", "mim_morbid_accession"), values = interactome.human_homologs[,"hsapiens_homolog_ensembl_gene"], filters = "ensembl_gene_id", mart = human)
interactome.WikiGene <- getBM(attributes = c("ensembl_gene_id", "wikigene_description"), values = interactome$ensembl_gene_id, filters = "ensembl_gene_id", mart = mouse)
interactome.InterPro <- getBM(attributes = c("ensembl_gene_id", "interpro_short_description"), values = interactome$ensembl_gene_id, filters = "ensembl_gene_id", mart = mouse)
interactome.GOslim <- getBM(attributes = c("ensembl_gene_id", "goslim_goa_accession", "goslim_goa_description"), values = interactome$ensembl_gene_id, filters = "ensembl_gene_id", mart = mouse)

# remove ensembl_gene_ids which have duplicated entrez_ids
interactome.entrez_ids <- interactome.entrez_ids[-which(duplicated(interactome.entrez_ids$ensembl_gene_id)),]
interactome <- merge(interactome, interactome.entrez_ids, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T)

# better not to merge as 1->many relationships
#interactome <- merge(interactome, interactome.entrez_ids, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T)
#interactome <- merge(interactome, interactome.WikiGene, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T)
#interactome <- merge(interactome, interactome.GOslim, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T)
interactome.entrezIDs <- unique(interactome[!is.na(interactome$entrezgene),]$entrezgene)
interactome.keggIDs <- keggConv.batch(interactome.entrezIDs)
keggQ <- lapply(interactome.keggIDs, function(x) keggGet(x))
interactome.pathways <- unique(unlist(lapply(strsplit(names(unlist(lapply(keggQ, function(x) x[[1]]$"PATHWAY"))), "\\."), function(x) x[3])))
interactome.pathways.genes <- lapply(interactome.pathways, function(x) keggLink("genes", x))
names(interactome.pathways.genes) <- interactome.pathways
interactome.pathways.genes.entrez_ids <- unique(gsub("mmu:", "", as.character(unlist(interactome.pathways.genes))))

# create dataframe for counting hits in pathways
interactome.df <- kegg.brite[gsub("mmu", "", interactome.pathways), ]
interactome.df$source <- rep("Interactome", nrow(interactome.df))
interactome.df$ID <- rownames(interactome.df)
# we are now using WCL as background to test for enrichment
i1 <- intersect(rownames(interactome.df), rownames(wcl.df))
interactome.df$total <- rep(0, nrow(interactome.df))
interactome.df[i1,]$total <- wcl.df[i1,]$count
interactome.df$count <- rep(0, nrow(interactome.df))
interactome.df$frac <- rep(0, nrow(interactome.df))

for (i in rownames(interactome.df)) {
  kL1 <- keggLink("mmu", paste("mmu", i, sep = ""))
  interactome.df[i, ]$count <- length(which(interactome.keggIDs %in% kL1))
  interactome.df[i, ]$frac <- round(length(which(interactome.keggIDs %in% kL1)) / length(kL1) * 100, 2)
}

# extract list of IDs in pathway
interactome.in_path.IDs <- lapply(rownames(interactome.df), function(x) {
  kL1 <- keggLink("mmu", paste("mmu", x, sep = ""))
  in_path <- interactome.keggIDs[which(interactome.keggIDs %in% kL1)]
})

# perform Fisher's Exact Test for each category
bkgd <- length(unique(wcl.keggIDs))
smpl <- length(interactome.keggIDs)

ftl <- apply(interactome.df, 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = alternative)
})

interactome.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
interactome.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
interactome.df$ft_fdr <- p.adjust(interactome.df$ft_pval, method = p.adjust.method, n = nrow(wcl.df))

#---------interactome-summarizing data at "B" level before doing Fisher's Exact test------
interactome.B.df <- data.frame(matrix(ncol = 5, nrow = length(unique(interactome.df$B))))
colnames(interactome.B.df) <- c("B", "A", "total", "count", "source")
interactome.B.df$B <- unique(interactome.df$B)
interactome.B.df$A <- sapply(unique(interactome.df$B), function(x) {A <- unique(interactome.df[which(interactome.df$B %in% x), "A"])})
interactome.B.df$source <- rep("Interactome", nrow(interactome.B.df))
interactome.B.df$total <- sapply(unique(interactome.df$B), function(x) {tot <- sum(interactome.df[which(interactome.df$B %in% x), "total"])})
interactome.B.df$count <- sapply(unique(interactome.df$B), function(x) {count <- sum(interactome.df[which(interactome.df$B %in% x), "count"])})

bkgd <- length(unique(wcl.keggIDs))
smpl <- length(interactome.keggIDs)

ftl <- apply(interactome.B.df, 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = alternative)
})

interactome.B.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
interactome.B.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
interactome.B.df$ft_fdr <- p.adjust(interactome.B.df$ft_pval, method = p.adjust.method, n = nrow(wcl.df))

#-----------GO RNA unrelated-------------------------------------------
# subset the interactome table
interactome.go_rna_unrelated <- interactome[which(interactome$GO == "unrelated"),]
interactome.go_rna_unrelated.entrezIDs <- unique(interactome.go_rna_unrelated[!is.na(interactome.go_rna_unrelated$entrezgene),]$entrezgene)
interactome.go_rna_unrelated.keggIDs <- keggConv.batch(interactome.go_rna_unrelated.entrezIDs)

# dataframe for count data
interactome.go_rna_unrelated.df <- interactome.df
interactome.go_rna_unrelated.df$source <- rep("GO_RNA_unrelated", nrow(interactome.go_rna_unrelated.df))
interactome.go_rna_unrelated.df$ID <- rownames(interactome.go_rna_unrelated.df)
interactome.go_rna_unrelated.df$total <- rep(0, nrow(interactome.go_rna_unrelated.df))

# we are now using WCL as background to test for enrichment
i1 <- intersect(rownames(interactome.go_rna_unrelated.df), rownames(wcl.df))
interactome.go_rna_unrelated.df[i1,]$total <- wcl.df[i1,]$count
interactome.go_rna_unrelated.df$count <- rep(0, nrow(interactome.go_rna_unrelated.df))

for (i in rownames(interactome.go_rna_unrelated.df)) {
  kL1 <- keggLink("mmu", paste("mmu", i, sep = ""))
  interactome.go_rna_unrelated.df[i, ]$count <- length(which(interactome.go_rna_unrelated.keggIDs %in% kL1))
}

# extract list of IDs in pathway
interactome.go_rna_unrelated.in_path.IDs <- lapply(rownames(interactome.go_rna_unrelated.df), function(x) {
  kL1 <- keggLink("mmu", paste("mmu", x, sep = ""))
  in_path <- interactome.go_rna_unrelated.keggIDs[which(interactome.go_rna_unrelated.keggIDs %in% kL1)]
})
names(interactome.go_rna_unrelated.in_path.IDs) <- rownames(interactome.go_rna_unrelated.df)

# perform Fisher's Exact Test for each category
# Using WCL as background
bkgd <- length(unique(wcl.keggIDs))
smpl <- length(interactome.go_rna_unrelated.keggIDs)

ftl <- apply(interactome.go_rna_unrelated.df, 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = alternative)
})

interactome.go_rna_unrelated.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
interactome.go_rna_unrelated.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
interactome.go_rna_unrelated.df$ft_fdr <- p.adjust(interactome.go_rna_unrelated.df$ft_pval, method = p.adjust.method, n = nrow(wcl.df))

# summarizing data at "B" level before doing Fisher's Exact test
interactome.go_rna_unrelated.B.df <- data.frame(matrix(ncol = 5, nrow = length(unique(interactome.go_rna_unrelated.df$B))))
colnames(interactome.go_rna_unrelated.B.df) <- c("B", "A", "total", "count", "source")
interactome.go_rna_unrelated.B.df$B <- unique(interactome.go_rna_unrelated.df$B)
interactome.go_rna_unrelated.B.df$A <- sapply(unique(interactome.go_rna_unrelated.df$B), function(x) {A <- unique(interactome.go_rna_unrelated.df[which(interactome.go_rna_unrelated.df$B %in% x), "A"])})
interactome.go_rna_unrelated.B.df$source <- rep("GO_RNA_unrelated", nrow(interactome.go_rna_unrelated.B.df))
interactome.go_rna_unrelated.B.df$total <- sapply(unique(interactome.go_rna_unrelated.df$B), function(x) {tot <- sum(interactome.go_rna_unrelated.df[which(interactome.go_rna_unrelated.df$B %in% x), "total"])})
interactome.go_rna_unrelated.B.df$count <- sapply(unique(interactome.go_rna_unrelated.df$B), function(x) {count <- sum(interactome.go_rna_unrelated.df[which(interactome.go_rna_unrelated.df$B %in% x), "count"])})

# using WCL as background
ftl <- apply(interactome.go_rna_unrelated.B.df, 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = alternative)
})

interactome.go_rna_unrelated.B.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
interactome.go_rna_unrelated.B.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
interactome.go_rna_unrelated.B.df$ft_fdr <- p.adjust(interactome.go_rna_unrelated.B.df$ft_pval, method = p.adjust.method, n = nrow(wcl.df))

#-----------GO RNA related-------------------------------------------
interactome.go_rna_related <- interactome[-which(interactome$GO == "unrelated"),]

interactome.go_rna_related.entrezIDs <- unique(interactome.go_rna_related[!is.na(interactome.go_rna_related$entrezgene),]$entrezgene)
interactome.go_rna_related.keggIDs <- keggConv.batch(interactome.go_rna_related.entrezIDs)

# we are testing this subset of "interactome", therefore we include all the pathways from "interactome"
interactome.go_rna_related.df <- interactome.df
# TODO - make sure that same background is used in all tests!!!
# we are now using interactome as background to test for enrichment
i1 <- intersect(rownames(interactome.go_rna_related.df), rownames(wcl.df))
interactome.go_rna_related.df$total <- rep(0, nrow(interactome.go_rna_related.df))
interactome.go_rna_related.df[i1,]$total <- wcl.df[i1,]$count
interactome.go_rna_related.df$source <- rep("GO_RNA_related", nrow(interactome.go_rna_related.df))
interactome.go_rna_related.df$ID <- rownames(interactome.go_rna_related.df)
interactome.go_rna_related.df$count <- rep(0, nrow(interactome.go_rna_related.df))
interactome.go_rna_related.df$frac <- rep(0, nrow(interactome.go_rna_related.df))

for (i in rownames(interactome.go_rna_related.df)) {
  kL1 <- keggLink("mmu", paste("mmu", i, sep = ""))
  interactome.go_rna_related.df[i, ]$count <- length(which(interactome.go_rna_related.keggIDs %in% kL1))
}

# extract list of IDs in pathway
interactome.go_rna_related.in_path.IDs <- lapply(rownames(interactome.go_rna_related.df), function(x) {
  kL1 <- keggLink("mmu", paste("mmu", x, sep = ""))
  in_path <- interactome.go_rna_related.keggIDs[which(interactome.go_rna_related.keggIDs %in% kL1)]
})
names(interactome.go_rna_related.in_path.IDs) <- rownames(interactome.go_rna_related.df)

# perform Fisher's Exact Test for each category
# Using WCL as background
bkgd <- length(unique(wcl.keggIDs))
smpl <- length(interactome.go_rna_related.keggIDs)

ftl <- apply(interactome.go_rna_related.df, 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = alternative)
})

interactome.go_rna_related.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
interactome.go_rna_related.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
# changed p.adjust.method to "p.adjust.method" for more conservative control of p values, and set number to number of pathways in WCL
interactome.go_rna_related.df$ft_fdr <- p.adjust(interactome.go_rna_related.df$ft_pval, method = p.adjust.method, n = nrow(wcl.df))

# summarizing data at "B" level before doing Fisher's Exact test
interactome.go_rna_related.B.df <- data.frame(matrix(ncol = 5, nrow = length(unique(interactome.go_rna_related.df$B))))
colnames(interactome.go_rna_related.B.df) <- c("B", "A", "total", "count", "source")
interactome.go_rna_related.B.df$B <- unique(interactome.go_rna_related.df$B)
interactome.go_rna_related.B.df$A <- sapply(unique(interactome.go_rna_related.df$B), function(x) {A <- unique(interactome.go_rna_related.df[which(interactome.go_rna_related.df$B %in% x), "A"])})
interactome.go_rna_related.B.df$source <- rep("GO_RNA_related", nrow(interactome.go_rna_related.B.df))
interactome.go_rna_related.B.df$total <- sapply(unique(interactome.go_rna_related.df$B), function(x) {tot <- sum(interactome.go_rna_related.df[which(interactome.go_rna_related.df$B %in% x), "total"])})
interactome.go_rna_related.B.df$count <- sapply(unique(interactome.go_rna_related.df$B), function(x) {count <- sum(interactome.go_rna_related.df[which(interactome.go_rna_related.df$B %in% x), "count"])})

ftl <- apply(interactome.go_rna_related.B.df, 1, function (x) {
  ct <- as.integer(x["count"])
  tt <- as.integer(x["total"])
  m1 <- matrix(c(ct, tt, smpl - ct, bkgd - tt), 2, 2)
  fisher.test(m1, alternative = alternative)
})

interactome.go_rna_related.B.df$ft_pval <- unlist(lapply(ftl, function(x) {x$p.value}))
interactome.go_rna_related.B.df$ft_OR <- unlist(lapply(ftl, function(x) {x$estimate}))
# changed p.adjust.method to "p.adjust.method" for more conservative control of p values, and set number to number of pathways in WCL
interactome.go_rna_related.B.df$ft_fdr <- p.adjust(interactome.go_rna_related.B.df$ft_pval, method = p.adjust.method, n = nrow(wcl.df))

#-----------plotting of KEGG enrichment analysis results---------------------------------
df1 <- rbind(interactome.B.df[, c("A", "B", "ft_OR", "ft_fdr", "source")],
             interactome.go_rna_related.B.df[, c("A", "B", "ft_OR", "ft_fdr", "source")], 
             interactome.go_rna_unrelated.B.df[, c("A", "B", "ft_OR", "ft_fdr", "source")]
             )

df1$source <- as.factor(df1$source)
df1$source <- factor(df1$source, levels = levels(df1$source)[c(3,1,2)])

df1$ft_OR.cut <- cut(log2(df1$ft_OR), breaks = c(-Inf,-4:4), right = F)

p3 <- ggplot(df1, aes(B, source)) + geom_tile(aes(fill = (df1$ft_OR.cut)))
p3 <- p3 + theme(axis.text.x = element_text(angle = 90))
p3 <- p3 + scale_fill_gradientn(colours= c("red", "green"))
p3

# subsetting for FDR <= 0.025
c1 <- unique(as.character(df1[which(df1$ft_fdr <= 0.025),]$B))
df2 <- df1[which(df1$B %in% c1),]
df2$B <- as.factor(as.character(df2$B))
df2$ft_OR.cut <- cut(log2(df2$ft_OR), breaks = c(-Inf,-4:4), right = F)
df2$B <- factor(df2$B, levels = levels(df2$B)[df2[df2$source == "Interactome", "B"][order(df2[which(df2$source == "Interactome"),]$ft_OR.cut)]])
ggplot(data = df2, aes(x = source, y = B)) + geom_tile(aes(fill = ft_OR.cut), colour = "white") + scale_fill_brewer(palette = "PRGn") + theme(axis.text.x = element_text(angle = 90))

# subsetting Metabolism and Genetic Information Processing
df.metab <- df1[df1$A == "Metabolism",]
df.gip <- df1[df1$A == "Genetic information Processin",]
df3 <- rbind(df.metab[which(df.metab$ft_fdr <= 0.1),], df.gip[df.gip$ft_fdr <= 0.05,])
df3 <- rbind(df.metab, df.gip)
df3$B <- factor(df3$B, levels = levels(df3$B)[df2[df3$source == "Interactome", "B"][order(df2[which(df2$source == "Interactome"),]$ft_OR.cut)]])
p4 <- ggplot(df3, aes(source, B)) + geom_tile(aes(fill = df3$ft_OR.cut))
p4 <- p4 + theme(axis.text.x = element_text(angle = 90))
p4 <- p4 + scale_fill_brewer(palette = "PRGn") 
p4

#---------Plot at KEGG C level-------------------------
dfC <- rbind(interactome.df[, c("A", "B", "C", "ft_OR", "ft_fdr", "source", "count")],
             interactome.go_rna_related.df[, c("A", "B", "C", "ft_OR", "ft_fdr", "source", "count")], 
             interactome.go_rna_unrelated.df[, c("A", "B", "C", "ft_OR", "ft_fdr", "source", "count")]
)

dfC$source <- as.factor(dfC$source)
dfC$source <- factor(dfC$source, levels = levels(dfC$source)[c(3,1,2)])

select1 <- unique(as.character(dfC[which(dfC$ft_fdr <= 0.1 & dfC$ft_OR > 1),]$C))
select1.pathIDs <- paste("mmu", unlist(lapply(strsplit(select1, "\\ "), function(x) x[1])), sep = "")
dfC <- dfC[which(dfC$C %in% select1),]

dfC$C <- as.factor(as.character(dfC$C))
dfC$ft_OR.cut <- cut(log2(dfC$ft_OR), breaks = c(-Inf,-4:4), right = F)
dfC$C <- factor(dfC$C, levels = levels(dfC$C)[dfC[dfC$source == "Interactome", "C"][order(dfC[which(dfC$source == "Interactome"),]$ft_OR.cut, decreasing = T)]])



# formatting labels etc for plotting
l1 <- levels(dfC$ft_OR.cut)
l1 <- gsub("\\[", "", l1)
l1 <- gsub("\\)", "", l1)
levels(dfC$ft_OR.cut) <- l1

l1 <- as.character(levels(dfC$C))
l1 <- unlist(lapply(strsplit(l1, " "), function(x) {
  for (i in 2:length(x)){
    if (i == 2){
      v <- x[i]
    } else {
      v <- paste(v, x[i])
    }
  }
  return(v)
}))
levels(dfC$C) <- l1

levels(dfC$source)[2:3] <- c("RNA-related", "RNA-unrelated")


levels(dfC$C)[3] <- "Ribosome biogenesis"
levels(dfC$C)[6] <- "TCA cycle"
levels(dfC$C)[7] <- "mRNA surveillance"
levels(dfC$C)[11] <- "AA biosynthesis"
levels(dfC$C)[8] <- "H. simplex infection"
levels(dfC$C)[9] <- "Antibiotic biosynthesis"
levels(dfC$C)[12] <- "Glycolysis/Gluconeogenesis"

flevels <- levels(dfC$source)

l1 <- factor(dfC$C, levels = levels(dfC$C)[c(1,2,3,7,4,5,6,12,8,9,10,11)])
dfC$C <- l1
#levels(dfC$C)[1:4] <- c("I", "II", "III", "IV") # Ribosome, RNA transport, Ribosome biogenesis, mRNA surveillance
#levels(dfC$C)[10:12] <- c("V", "VI", "VII") # Antibiotic biosynthesis, Carbon metabolism, AA biosynthesis


p1 <- ggplot(data = dfC, aes(y = source, x = C)) + 
      geom_tile(aes(fill = ft_OR.cut), colour = "white") + 
      scale_fill_manual(values = brewer_pal(pal = "PuOr")(8), labels = levels(dfC$ft_OR.cut)) + #
      theme(axis.text.y = element_text(angle = 0, size = 5), axis.title = element_blank()) +
      guides(fill = guide_legend(label.position = "bottom", direction = "horizontal")) +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 0.8, size = 6.6)) +
      labs(fill = "Log2 OR") +
      scale_y_discrete(limits = rev(flevels)) +
      theme(legend.position = c(0.4,-1.92),
            legend.text = element_text(size = 4),
            legend.text.align = 0.5,
            legend.title = element_text(size = 4, vjust = 5),
            legend.key.size = unit(3.5, "mm"),
            legend.key.width = unit(3.5, "mm"),
            legend.margin = unit(0, "mm"),
            panel.margin = unit(1, "mm"))
p1

ggsave("/Users/u1001407/Dropbox//REM project-Sebastian/Figure_1g_20150603.pdf", plot = p1, scale = 1, height = 50, width = 110, unit = "mm")

# get highlighted pathway maps
for (i in select1.pathIDs) {
  if (length(interactome.go_rna_unrelated.in_path.IDs[[gsub("mmu", "", i)]]) > 0){
    print(i)
    url.rna_unrelated <- mark.pathway.by.objects(i, unlist(interactome.go_rna_unrelated.in_path.IDs[[gsub("mmu", "", i)]]))
    print(url.rna_unrelated)
    download.file(url.rna_unrelated, paste(i, "_RNA_unrelated", ".png", sep = ""))
  }
  if (length(interactome.go_rna_related.in_path.IDs[[gsub("mmu", "", i)]]) > 0){
    url.rna_related <- mark.pathway.by.objects(i, unlist(interactome.go_rna_related.in_path.IDs[[gsub("mmu", "", i)]]))
    download.file(url.rna_related, paste(i, "_RNA_related", ".png", sep = ""))
  }
  if (length(interactome.in_path.IDs[[gsub("mmu", "", i)]]) > 0){
    url.interactome <- mark.pathway.by.objects(i, unlist(interactome.in_path.IDs[[gsub("mmu", "", i)]]))
    download.file(url.interactome, paste(i, "_Interactome", ".png", sep = ""))
  }
}

#---------Venn diagram for metabolic pathways------------------
m1 <- matrix(nrow = length(unique(unlist(interactome.in_path.IDs[c("00020", "00010", "01230", "01200", "01130")]))), ncol = 5)
rownames(m1) <- unique(unlist(interactome.in_path.IDs[c("00020", "00010", "01230", "01200", "01130")]))
colnames(m1) <- c("00020", "00010", "01230", "01200", "01130")

g00020 <- as.character(interactome.in_path.IDs[[c("00020")]])
g00010 <- as.character(interactome.in_path.IDs[[c("00010")]])
g01230 <- as.character(interactome.in_path.IDs[[c("01230")]])
g01200 <- as.character(interactome.in_path.IDs[[c("01200")]])
g01130 <- as.character(interactome.in_path.IDs[[c("01130")]])
universe <- unique(c(g00020, g00010, g01230, g01200, g01130))

g00020.1 <- universe %in% g00020
g00010.1 <- universe %in% g00010
g01230.1 <- universe %in% g01230
g01200.1 <- universe %in% g01200
g01130.1 <- universe %in% g01130

x <- "00020"
m1[which(rownames(m1) %in% (interactome.in_path.IDs[[x]])), x] <- T
x <- "00010"
m1[which(rownames(m1) %in% (interactome.in_path.IDs[[x]])), x] <- T
x <- "01230"
m1[which(rownames(m1) %in% (interactome.in_path.IDs[[x]])), x]
x <- "01200"
m1[which(rownames(m1) %in% (interactome.in_path.IDs[[x]])), x] <- T
x <- "01130"
m1[which(rownames(m1) %in% (interactome.in_path.IDs[[x]])), x] <-T
