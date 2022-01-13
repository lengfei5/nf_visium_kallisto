##########################################################################
##########################################################################
# Project: axolotl visium analysis
# Script purpose: make Seurat object for visium data and make loope file for 10x software
# Usage example:
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Jan 13 11:24:44 2022
##########################################################################
##########################################################################
library(Seurat)
library(DropletUtils)
library(edgeR)
##########################################
# missed function in DropletUtils
##########################################
estimateAmbience <- function(m, lower=100, by.rank=NULL, round=TRUE, good.turing=TRUE) {
  m <- .rounded_to_integer(m, round)
  totals <- .intColSums(m)
  lower <- .get_lower(totals, lower, by.rank=by.rank)

  if (good.turing) {
    a <- .compute_ambient_stats(m, totals, lower=lower)
    output <- numeric(nrow(m))
    output[!a$discard] <- a$ambient.prop
    names(output) <- rownames(m)
  } else {
    ambient <- totals <= lower
    output <- rowSums(m[,ambient,drop=FALSE])
  }

  output
}

.rounded_to_integer <- function(m, round=TRUE) {
  if (round) {
    cs <- colSums(m)
    rs <- rowSums(m)
    if (!isTRUE(all.equal(cs, round(cs))) ||
        !isTRUE(all.equal(rs, round(rs))))
    {
      m <- round(m)
    }
  }
  m
}

.intColSums <- function(m) {
  # Enforcing discreteness mainly for emptyDrops()'s Monte Carlo step.
  as.integer(round(colSums(m)))
}

#' @importFrom Matrix rowSums
.compute_ambient_stats <- function(m, totals, lower) {
  # This doesn't invalidate 'totals', by definition.
  discard <- rowSums(m) == 0
  if (any(discard)) {
    m <- m[!discard,,drop=FALSE]
  }

  # Computing the average profile from the ambient cells.
  ambient <- totals <= lower # lower => "T" in the text.
  ambient.m <- m[,ambient,drop=FALSE]
  ambient.prof <- rowSums(ambient.m)

  if (sum(ambient.prof)==0) {
    stop("no counts available to estimate the ambient profile")
  }
  ambient.prop <- .safe_good_turing(ambient.prof)

  list(
    m=m, # this MUST have the same number of columns as input.
    discard=discard,
    ambient=ambient,
    ambient.m=ambient.m,
    ambient.prop=ambient.prop
  )
}

.get_lower <- function(totals, lower, by.rank) {
  if (is.null(by.rank)) {
    lower
  } else if (by.rank >= length(totals)) {
    stop("not have enough columns for supplied 'by.rank'")
  } else {
    totals[order(totals, decreasing=TRUE)[by.rank+1]]
  }
}

#' @importFrom edgeR goodTuringProportions
.safe_good_turing <- function(ambient.prof) {
  ambient.prob <- goodTuringProportions(ambient.prof)

  still.zero <- ambient.prob<=0
  if (any(still.zero)) {
    pseudo.prob <- 1/sum(ambient.prof)
    ambient.prob[still.zero] <- pseudo.prob/sum(still.zero)
    ambient.prob[!still.zero] <- ambient.prob[!still.zero] * (1 - pseudo.prob)
  }

  ambient.prob
}


# import read in data
topdir = "./" # source dir
# topdir = '/Volumes/clustertmp/jiwang/visium_axolotl/Amex_d1_183623/'
exp = Matrix::readMM(paste0(topdir, "genecounts.mtx")) #read matrix
bc = read.csv(paste0(topdir, "genecounts.barcodes.txt"), header = F, stringsAsFactors = F)
g = read.csv(paste0(topdir, "genecounts.genes.txt"), header = F, stringsAsFactors = F)
dimnames(exp) = list(paste0(bc$V1,"-1"),g$V1) # number added because of seurat format for barcodes
count.data = Matrix::t(exp)

# plot rankings for number of UMI
br.out <- barcodeRanks(count.data)

pdf("UMIrank.pdf", height = 5, width = 5, useDingbats = F)

plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")

o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")
abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
    legend=c("knee", "inflection"))

dev.off()

# UMI duplication
umi = read.table(paste0(topdir, "umicount.txt"), sep = "\t", header = F, stringsAsFactors = F)

sumUMI = c()
sumi = sum(umi$V4)

for(i in 0:250){
  sumUMI = c(sumUMI, sum(umi$V4[umi$V4>i])/sumi)
}

pdf("UMIduplication.pdf", height = 3.5, width = 7, useDingbats = F)

par(mfrow = c(1,2))

plot(sumUMI, ylim = c(0,1), pch = 20, col = "grey30", ylab = "% of total reads",
     xlab = "More than xx UMI")
diffUMI = sumUMI[-length(sumUMI)] - sumUMI[-1]
plot(diffUMI, ylim = c(0,0.2), pch = 20, col = "grey30", ylab = "Change in % of total reads",
     xlab = "More than xx UMI")

dev.off()

# create Seurat object and add image information
srat <- CreateSeuratObject(counts = count.data, assay = "Spatial") # create object

image <- Read10X_Image(image.dir = paste0(topdir, "mock/outs/spatial/"),
                       filter.matrix = FALSE) # read in the images
image <- image[Cells(x = srat)] # filter image by the spots
DefaultAssay(object = image) <- "Spatial" # set default assay
srat[["slice1"]] <- image # slice name might be changed

# estimate ambient proportion
amb_prop = estimateAmbience(count.data)[rownames(srat@assays$Spatial@meta.features)]

srat@assays$Spatial@meta.features = data.frame(row.names = rownames(srat@assays$Spatial@meta.features),
                                                "ambient_prop" = amb_prop)

# get MT% (genes curated from NCBI chrMT genes)
# mtgenes = c("COX1", "COX2", "COX3", "ATP6", "ND1", "ND5", "CYTB", "ND2", "ND4",
#             "ATP8", "MT-CO1", "COI", "LOC9829747")
# mtgenes = c(mtgenes, paste0("MT", mtgenes), paste0("MT-", mtgenes))
# mtgenes = mtgenes[mtgenes %in% g[,1]]
# srat = PercentageFeatureSet(srat, col.name = "percent.mt", assay = "Spatial",
#                             features = mtgenes)

saveRDS(srat, file = "srat.RDS")
