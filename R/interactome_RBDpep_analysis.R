#  interactome_RBDpep_analysis.R
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

#---------load libraries---------------------------------------------
library("biomaRt")
library("gdata")
library("ggplot2")
library("Biostrings")
library("bio3d")

#---------create biomaRt connections---------------------------------
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mmus.attribs <- listAttributes(mouse)

mouse#---------load data--------------------------------------------------
RBDpep.hl1 <- read.xls("~/Dropbox/REM project-Sebastian/RBDpep analysis/HL-1 RBDmap peptide.xlsx", sheet = 1, as.is = T)
RBDpep.HeLa <- read.xls("~/Dropbox/REM project-Sebastian/RBDpep analysis/HeLa RBDmap peptides.xlsx", sheet = 1, as.is = T)

# sorting tables by Ensembl Gene ID and start position of fragment to "linearize" data
RBDpep.hl1 <- RBDpep.hl1[order(RBDpep.hl1$ENSMBL.gene.ID, RBDpep.hl1$Start), ]
RBDpep.HeLa <- RBDpep.HeLa[order(RBDpep.HeLa$ENSG, RBDpep.HeLa$Start), ]

#---------get human homologs of mouse [HL-1] proteins----------------
mmus.RBDpep.hsap.homologs <- getBM(attributes = c("ensembl_gene_id", "description", "hsapiens_homolog_ensembl_gene"),
                                   filter = "ensembl_gene_id",
                                   values = RBDpep.hl1$ENSMBL.gene.ID,
                                   mart = mouse)

# have a look at SwissProt/TrEMBL UniProt IDs
mmus.uniprot <- getBM(attributes = c("ensembl_gene_id", "uniprot_sptrembl", "uniprot_swissprot"),
                      filter = "ensembl_gene_id",
                      values = RBDpep.hl1$ENSMBL.gene.ID,
                      mart = mouse)
# create intersect
hsap.i1 <- intersect(RBDpep.HeLa$ENSG, mmus.RBDpep.hsap.homologs$hsapiens_homolog_ensembl_gene)
mmus.i1 <- mmus.RBDpep.hsap.homologs[mmus.RBDpep.hsap.homologs$hsapiens_homolog_ensembl_gene %in% hsap.i1, ]$ensembl_gene_id

length(unique(mmus.i1)) #176/411 - 42%
length(unique(hsap.i1)) #177/368 - 48%

#---------do pairwise alignment of intersect proteins-----------------

aln.blosum62 <- sapply(mmus.i1, function(x) {
  mmus.frag <- RBDpep.hl1[RBDpep.hl1$ENSMBL.gene.ID == x, ]$Fragment.sequence
  hsap.homolog <- mmus.RBDpep.hsap.homologs[mmus.RBDpep.hsap.homologs$ensembl_gene_id == x, ]$hsapiens_homolog_ensembl_gene
  hsap.frag <- RBDpep.HeLa[RBDpep.HeLa$ENSG %in% hsap.homolog, ]$fragmentSequence
  p <- as.character(hsap.frag)
  s <- as.character(mmus.frag)
  pA <- sapply(s, function(sx) {
    p1 <- pairwiseAlignment(pattern = AAStringSet(p),
                            subject = AAString(sx),
                            substitutionMatrix = "BLOSUM62",
                            gapOpening = -12,
                            gapExtension = -5,
                            type = "global-local")
    attr(p1, "pid") <- pid(p1)
    attr(p1, "cStr") <- compareStrings(p1)
    return(p1)
  })
})

# pattern is fragment from human
# subject is fragment from mouse
aln.best <- lapply(aln.blosum62, function(x) lapply(x, function(y) {r1 <- y[which.max(pid(y))]; attr(r1, "pid") <- pid(r1); return(r1)}))

# subject is fragment from human
# pattern is fragment from mouse
aln.best2 <- lapply(aln.blosum62, function(x) lapply(x, function(y) {r1 <- y[which.max(pid(y))]; attr(r1, "pid") <- pid(r1); return(r1)}))
aln.best <- aln.best2

mmus.homolog <- mmus.RBDpep.hsap.homologs[mmus.RBDpep.hsap.homologs$hsapiens_homolog_ensembl_gene == hsap.i1[grep("ENSG00000135316", hsap.i1)], ]$ensembl_gene_id
RBDpep.HeLa[RBDpep.HeLa$ENSG == hsap.i1[grep("ENSG00000135316", hsap.i1)],]
RBDpep.hl1[RBDpep.hl1$ENSMBL.gene.ID == mmus.homolog,]

RBDpep.merge <- RBDpep.hl1[RBDpep.hl1$ENSG %in% mmus.i1,]
RBDpep.merge$mmusHomolog <- as.character(sapply(RBDpep.merge$ENSG, function(x) mmus.RBDpep.hsap.homologs[mmus.RBDpep.hsap.homologs$hsapiens_homolog_ensembl_gene == x,]$ensembl_gene_id))
RBDpep.merge$mmusAlignment <- unlist(lapply(unlist(aln.best[unique(RBDpep.merge$ENSG)]), function(x) compareStrings(x)))
RBDpep.merge$mmusSimilarity <- unlist(lapply(unlist(aln.best[unique(RBDpep.merge$ENSG)]), function(x) attr(x, "pid")))
RBDpep.merge$mmusScore <- unlist(lapply(unlist(aln.best[unique(RBDpep.merge$ENSG)]), function(x) attr(x, "score")))
RBDpep.merge$mmusFragment <- unlist(lapply(unlist(aln.best[unique(RBDpep.merge$ENSG)]), function(x) toString(unaligned(pattern(x)))))
RBDpep.merge$mmusFragmentStart <- unlist(lapply(apply(RBDpep.merge, 1, function(x) RBDpep.hl1[RBDpep.hl1$ENSMBL.gene.ID == x["mmusHomolog"] & RBDpep.hl1$Fragment.sequence == x["mmusFragment"], ]), function(y) unique(y["Fragment.start"])))
RBDpep.merge$mmusFragmentStop <- unlist(lapply(apply(RBDpep.merge, 1, function(x) RBDpep.hl1[RBDpep.hl1$ENSMBL.gene.ID == x["mmusHomolog"] & RBDpep.hl1$Fragment.sequence == x["mmusFragment"], ]), function(y) unique(y["Fragment.stop"])))


# Mouse as pattern
RBDpep.merge <- RBDpep.hl1[RBDpep.hl1$ENSMBL.gene.ID %in% mmus.i1,]
RBDpep.merge$hsapHomolog <- as.character(sapply(RBDpep.merge$ENSMBL.gene.ID, function(x) mmus.RBDpep.hsap.homologs[mmus.RBDpep.hsap.homologs$ensembl_gene_id == x,]$hsapiens_homolog_ensembl_gene))

RBDpep.merge$hsapAlignment <- unlist(lapply(unlist(aln.best[unique(RBDpep.merge$ENSMBL.gene.ID)]), function(x) compareStrings(x)))
RBDpep.merge$hsapSimilarity <- unlist(lapply(unlist(aln.best[unique(RBDpep.merge$ENSMBL.gene.ID)]), function(x) attr(x, "pid")))
RBDpep.merge$hsapScore <- unlist(lapply(unlist(aln.best[unique(RBDpep.merge$ENSMBL.gene.ID)]), function(x) attr(x, "score")))
RBDpep.merge$hsapFragment <- unlist(lapply(unlist(aln.best[unique(RBDpep.merge$ENSMBL.gene.ID)]), function(x) toString(unaligned(pattern(x)))))
RBDpep.merge$hsapFragmentStart <- unlist(lapply(apply(RBDpep.merge, 1, function(x) RBDpep.HeLa[RBDpep.HeLa$fragmentSequence == x["hsapFragment"], ]), function(y) unique(y["fragmentStart"])))
RBDpep.merge$hsapFragmentStop <- unlist(lapply(apply(RBDpep.merge, 1, function(x) RBDpep.HeLa[RBDpep.HeLa$fragmentSequence == x["hsapFragment"], ]), function(y) unique(y["fragmentStop"])))

pdf("~/Dropbox/REM project-Sebastian/RBDpep analysis/MMus_histogram.pdf")
h1 <- hist(RBDpep.merge[!duplicated(RBDpep.merge$hsapFragment),]$hsapSimilarity, labels = T, col = "gray", main = "Frequency of similarity\n HL-1 RBDpep fragments vs. HeLa fragments [N = 237]", xlab = "Similarity [%]")
dev.off()
sum(h1$counts)

write.csv(RBDpep.merge[!duplicated(RBDpep.merge$hsapFragment),], file = "~/Dropbox/REM project-Sebastian/RBDpep analysis/Mmus_nique.csv")

write.csv(RBDpep.merge, file = "~/Dropbox/REM project-Sebastian/RBDpep analysis/Human_mouse_intersect.csv")

