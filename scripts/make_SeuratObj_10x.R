
#!/usr/local/bin/Rscript --vanilla

library(Seurat)
library(DropletUtils)

# read in data
topdir = "${outbus}" # source dir
exp = Matrix::readMM(paste0(topdir, "/genecounts.mtx")) #read matrix
bc = read.csv(paste0(topdir, "/genecounts.barcodes.txt"), header = F, stringsAsFactors = F)
g = read.csv(paste0(topdir, "/genecounts.genes.txt"), header = F, stringsAsFactors = F)
dimnames(exp) = list(paste0(bc\$V1,"-1"),g\$V1) # number added because of seurat format for barcodes
count.data = Matrix::t(exp)

# get emptyDrops and default cutoff cell estimates
iscell_dd = defaultDrops(count.data, expected = 5000)
eout = emptyDrops(count.data, lower = 200)
eout\$FDR[is.na(eout\$FDR)] = 1
iscell_ed = eout\$FDR<=0.01
meta = data.frame(row.names = paste0(bc\$V1,"-1"),
                  iscell_dd = iscell_dd, iscell_ed = iscell_ed)

# plot rankings for number of UMI
br.out <- barcodeRanks(count.data)
pdf("${params.samplename}_UMIrank.pdf", height = 5, width = 5, useDingbats = F)
plot(br.out\$rank, br.out\$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out\$rank)
lines(br.out\$rank[o], br.out\$fitted[o], col="red")
abline(h=metadata(br.out)\$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)\$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
    legend=c("knee", "inflection"))
dev.off()

# UMI duplication
umi = read.table("${umic}", sep = "\t", header = F, stringsAsFactors = F)
sumUMI = c()
sumi = sum(umi\$V4)
for(i in 0:250){ sumUMI = c(sumUMI, sum(umi\$V4[umi\$V4>i])/sumi) }
pdf("${params.samplename}_UMIduplication.pdf", height = 3.5, width = 7, useDingbats = F)
par(mfrow = c(1,2))
plot(sumUMI, ylim = c(0,1), pch = 20, col = "grey30", ylab = "% of total reads",
     xlab = "More than xx UMI", main = "${params.samplename}")
diffUMI = sumUMI[-length(sumUMI)] - sumUMI[-1]
plot(diffUMI, ylim = c(0,0.2), pch = 20, col = "grey30", ylab = "Change in % of total reads",
     xlab = "More than xx UMI", main = "${params.samplename}")
dev.off()

# create Seurat object
## we're only keeping what might potentially be a cell (by DD or ED)
srat = CreateSeuratObject(counts = count.data[,iscell_dd | iscell_ed],
                          meta.data = meta[iscell_dd | iscell_ed,])
amb_prop = estimateAmbience(count.data)[rownames(srat@assays\$RNA@meta.features)]
srat@assays\$RNA@meta.features = data.frame(row.names = rownames(srat@assays\$RNA@meta.features),
                                            "ambient_prop" = amb_prop)

# get MT% (genes curated from NCBI chrMT genes)
# mtgenes = c("COX1", "COX2", "COX3", "ATP6", "ND1", "ND5", "CYTB", "ND2", "ND4",
#            "ATP8", "MT-CO1", "COI", "LOC9829747")
# mtgenes = c(mtgenes, paste0("MT", mtgenes), paste0("MT-", mtgenes))
# mtgenes = mtgenes[mtgenes %in% g[,1]]
#srat = PercentageFeatureSet(srat, col.name = "percent.mt", assay = "RNA",
#                            features = mtgenes)

# saveRDS(srat, file = "${params.samplename}_srat.RDS")
saveRDS(srat, file = "srat.RDS")
