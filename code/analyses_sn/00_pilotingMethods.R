### LC snRNA-seq analysis
### Build SCE from raw count matrices; perform nuclei calling
### qrsh -l bluejay,mf=32G,h_vmem=34G -pe local 2
#     (which points to this script with `R CMD BATCH`)
### Adopted from https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/build_basic_sce.R
### Initiated MNT 13Dec2021

library(SingleCellExperiment)
library(DropletUtils)
library(rtracklayer)
library(BiocParallel)
library(jaffelab)
library(here)
library(sessioninfo)

here()
# [1] "/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c"


### MNT test 17Dec2021 ==================
load(here("processed_data","SCE","sce_raw_LC.rda"), verbose=T)
# sce.lc, e.out.lc

sample.idx <- splitit(sce.lc$Sample)
i <- "Br6522_LC"


## Create barcode rank plots from `DropletUtils::barcodeRanks()`; get stats ===
BRdf.list <- list()
Sys.time()
# [1] "2021-12-17 12:46:19 EST"
for(i in names(sample.idx)){
  BRdf.list[[i]] <- barcodeRanks(counts(sce.lc[ ,sample.idx[[i]]]))
}
Sys.time()
# 

sapply(BRdf.list, metadata)
#            Br2701_LC Br6522_LC Br8079_LC
# knee       27658     193       196      
# inflection 8999      106       102

# So interestingly this would have detected that 'second knee' for two samples


# Plotting (adapted from the manual under `barcodeRanks()`) ===
for(i in names(BRdf.list)){
png(here("plots","snRNA-seq",paste0("barcodeRankPlot_",i,"_dropletUtils-defaults.png")))
  br.out <- BRdf.list[[i]]
  
  plot(br.out$rank, br.out$total, log="xy", xlab="Barcode Rank", ylab="Total UMIs",
       main=paste0("Barcode rank plot for: ",i), cex.axis=0.8, cex.lab=1.2, las=1)
  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
  abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
         legend=c("knee", "inflection"))
dev.off()
}


e.out.custom <- list()
Sys.time()
# 
#for(i in names(sample.idx)){
cat(paste0("Simulating empty drops for: ", i,"... \n"))

## Run emptyDrops with this data-defined `lower=` param 
set.seed(109)
e.out.br[[i]] <- emptyDrops(counts(sce.lc[ ,sample.idx[[i]]]), niters=20000,
                            lower=customLower,
                            BPPARAM=BiocParallel::MulticoreParam(2))
cat(paste0("\n\t...Simulations complete. \n\t", Sys.time(), "\n\n\n"))
Sys.time()

#}



