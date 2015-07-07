

library(gdata)
library(biomaRt)
library(GO.db)
library(ggplot2)

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
attribs <- listAttributes(mouse)
filters <- listFilters(mouse)
attribs.hsap <- listAttributes(human)

# load IDs
cv.assoc.proteins <- read.xls("/Volumes/MHS//workgroups/jcsmr//PreissLab/Sebastian Kurscheid/Annotations//GO/cardiovascular_associated_proteins.xlsx", sheet = 1, header = T, as.is = T)
interactome <- read.xls("/Users/u1001407/Dropbox/REM project-Sebastian/HL-1 interactome superset.xlsx", sheet = "sheet 1" , as.is = T)
wcl <- read.xls("/Users/u1001407/Dropbox/REM project-Sebastian/HL-1 interactome superset.xlsx", sheet = "WCL GO RNAbind" , as.is = T)
colnames(interactome)[c(1,2)] <- c("ensembl_gene_id", "gene_symbol")
colnames(wcl) <- c("ensembl_gene_id", "gene_symbol", "GO Term", "RBD type")

# biomaRt attribute uniprot_swissprot
mmus.cv.assoc <- getBM(attributes = c("ensembl_gene_id", "uniprot_swissprot"), filters = "uniprot_swissprot", values = cv.assoc.proteins[which(cv.assoc.proteins$Taxon == "10090"), "ID"], mart = mouse)
hsap.cv.assoc <- getBM(attributes = c("ensembl_gene_id", "uniprot_swissprot"), filters = "uniprot_swissprot", values = cv.assoc.proteins[which(cv.assoc.proteins$Taxon == "9606"), "ID"], mart = human)
hsap.cv.assoc.mmus.homologs <- getBM(attributes = c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene"), filters = "uniprot_swissprot", values = cv.assoc.proteins[which(cv.assoc.proteins$Taxon == "9606"), "ID"], mart = human)

i1 <- intersect(unique(mmus.cv.assoc$ensembl_gene_id), interactome$ensembl_gene_id)
i2 <- intersect(unique(hsap.cv.assoc.mmus.homologs$mmusculus_homolog_ensembl_gene), interactome$ensembl_gene_id)
i3 <- c(i1[which(!i1 %in% intersect(i1, i2))], i2[which(!i2 %in% intersect(i1, i2))])

interactome.go_ids <- getBM(attributes = c("ensembl_gene_id", "go_id"), filters = "ensembl_gene_id", values = interactome$ensembl_gene_id, mart = mouse)
# subtract genes which are in common between interactome and WCL
wcl <- wcl[-which(wcl$ensembl_gene_id %in% interactome$ensembl_gene_id),]
wcl.go_ids <- getBM(attributes = c("ensembl_gene_id", "go_id"), filters = "ensembl_gene_id", values =wcl$ensembl_gene_id, mart = mouse)
  
cv.go_terms.bp <- c("GO:0007507", "GO:0048738", "GO:0008015", "GO:0050878", "GO:0001944", "GO:0042060", "GO:0006979", "GO:0016055", "GO:0006520", "GO:0050817", "GO:0006629", "GO:0006936", "GO:0048771", "GO:0051145", "GO:0007517", "GO:0042692", "GO:0048659")
cv.go_terms.cc <- c("GO:0005739", "GO:0005578")

#---------------------some plotting----------------------------------------------------------------------
# from GO.db
xx <- as.list(GOTERM)

go.bp.offspring <- as.list(GOBPOFFSPRING)
go.cc.offspring <- as.list(GOCCOFFSPRING)

interactome.cv.go_bp.offsp <- sapply(cv.go_terms.bp, function(x) length(unique(interactome.go_ids[which(interactome.go_ids$go_id %in% unlist(go.bp.offspring[x])), "ensembl_gene_id"])))
interactome.cv.go_cc.offsp <- sapply(cv.go_terms.cc, function(x) length(unique(interactome.go_ids[which(interactome.go_ids$go_id %in% unlist(go.cc.offspring[x])), "ensembl_gene_id"])))
interactome.cv.go_bp.offsp.IDs <- sapply(cv.go_terms.bp, function(x) unique(interactome.go_ids[which(interactome.go_ids$go_id %in% unlist(go.bp.offspring[x])), "ensembl_gene_id"]))
interactome.cv.go_cc.offsp.IDs <- sapply(cv.go_terms.cc, function(x) unique(interactome.go_ids[which(interactome.go_ids$go_id %in% unlist(go.cc.offspring[x])), "ensembl_gene_id"]))

wcl.cv.go_bp.offsp <- sapply(cv.go_terms.bp, function(x) length(unique(wcl.go_ids[which(wcl.go_ids$go_id %in% unlist(go.bp.offspring[x])), "ensembl_gene_id"])))
wcl.cv.go_cc.offsp <- sapply(cv.go_terms.cc, function(x) length(unique(wcl.go_ids[which(wcl.go_ids$go_id %in% unlist(go.cc.offspring[x])), "ensembl_gene_id"])))
wcl.cv.go_bp.offsp.IDs <- sapply(cv.go_terms.bp, function(x) unique(wcl.go_ids[which(wcl.go_ids$go_id %in% unlist(go.bp.offspring[x])), "ensembl_gene_id"]))
wcl.cv.go_cc.offsp.IDs <- sapply(cv.go_terms.cc, function(x) unique(wcl.go_ids[which(wcl.go_ids$go_id %in% unlist(go.cc.offspring[x])), "ensembl_gene_id"]))

df.go_bp.interactome <- as.data.frame(interactome.cv.go_bp.offsp)
colnames(df.go_bp.interactome) <- "count"
df.go_bp.interactome$group <- rep("interactome", nrow(df.go_bp.interactome))

df.go_bp.wcl <- as.data.frame(wcl.cv.go_bp.offsp)
colnames(df.go_bp.wcl) <- "count"
df.go_bp.wcl$group <- rep("wcl", nrow(df.go_bp.wcl))
df.go_bp.wcl$id <- rownames(df.go_bp.wcl)

df.go_bp <- rbind(df.go_bp.interactome, df.go_bp.wcl)
df.go_bp$id <- c(rownames(df.go_bp.interactome), rownames(df.go_bp.wcl))

df.go_bp$term <- sapply(df.go_bp$id, function(x) xx[x][[1]]@Term)
n1 <- length(unique(interactome.go_ids[which(interactome.go_ids$go_id %in% unlist(go.bp.offspring[cv.go_terms.bp])), "ensembl_gene_id"]))
n2 <- length(unique(wcl.go_ids[which(wcl.go_ids$go_id %in% unlist(go.bp.offspring[cv.go_terms.bp])), "ensembl_gene_id"]))

hist.go_bp <- ggplot(df.go_bp, aes(term, count, group = group, fill = group)) + geom_bar(postion = "dodge", stat = "identity")
hist.go_bp <- hist.go_bp + theme(axis.text.x = element_text(angle = 90))
hist.go_bp <- hist.go_bp + labs(title = paste("CV-associated gene counts\n in GO BP terms for Interactome [N = ", n1, "] and WCL only [N = ", n2, "]", sep = ""))
hist.go_bp

pdf("/Users/u1001407/Dropbox/REM project-Sebastian/WCL_and_Interactome_cardiovascular_assoc_genes_GO_BP_histogram.pdf", paper = "a4r")
hist.go_bp
dev.off()

df.go_cc <- as.data.frame(interactome.cv.go_cc.offsp)
colnames(df.go_cc)[1] <- "count"
df.go_cc$id <- rownames(df.go_cc)
df.go_cc$term <- sapply(rownames(df.go_cc), function(x) xx[x][[1]]@Term)
hist.go_cc <- ggplot(df.go_cc, aes(term, count)) + geom_histogram(stat = "identity", fill = "blue")
n <- length(unique(interactome.go_ids[which(interactome.go_ids$go_id %in% unlist(go.cc.offspring[cv.go_terms.cc])), "ensembl_gene_id"]))
hist.go_cc <- hist.go_cc + theme(axis.text.x = element_text(angle = 0))
hist.go_cc <- hist.go_cc + labs(title = paste("CV-associated Interactome genes\n GO CC [N = ", n, "]", sep = ""))

pdf("/Users/u1001407/Dropbox/REM project-Sebastian/Interactome_cardiovascular_assoc_genes_GO_CC_histogram.pdf", paper = "a4")
hist.go_cc
dev.off()

# make a table with GO terms and associated gene IDs
df1 <- data.frame(matrix(nrow = max(interactome.cv.go_bp.offsp) + 1, ncol = length(interactome.cv.go_bp.offsp)))
colnames(df1) <- names(interactome.cv.go_bp.offsp.IDs)
df1[1,] <- sapply(names(interactome.cv.go_bp.offsp.IDs), function(x) xx[x][[1]]@Term)
for (i in 1:length(interactome.cv.go_bp.offsp.IDs)) {
  df1[1 : length(interactome.cv.go_bp.offsp.IDs[[i]]) + 1 , i] <- as.vector(unlist(interactome.cv.go_bp.offsp.IDs[i]))
}

# make a table with GO terms and associated gene IDs
df1 <- data.frame(matrix(nrow = max(interactome.cv.go_cc.offsp) + 1, ncol = length(interactome.cv.go_cc.offsp)))
colnames(df1) <- names(interactome.cv.go_cc.offsp.IDs)
df1[1,] <- sapply(names(interactome.cv.go_cc.offsp.IDs), function(x) xx[x][[1]]@Term)
for (i in 1:length(interactome.cv.go_cc.offsp.IDs)) {
  df1[1 : length(interactome.cv.go_cc.offsp.IDs[[i]]) + 1 , i] <- as.vector(unlist(interactome.cv.go_cc.offsp.IDs[i]))
}
write.csv(df1, file = "/Users/u1001407/Dropbox/REM project-Sebastian/Interactome_cardiovascular_associated_GO_CC_gene_ensembl_IDs_table.csv")



