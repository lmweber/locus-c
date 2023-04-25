### MNT 29Nov2021 =======
  # Check Read 1 files for discrepant read lengths, as seen before in Tran-Maynard 2021 project:
  # Should be exactly 28bp = 16 [BC] + 12 [UMI]

library(ShortRead)
library(jaffelab)

FASTQ.dir <- "/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/fastq/snRNA-seq/"

### Read in abridged sample info (MNT generated for processing through CR)
samples.lc <- read.table("/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/fastq/snRNA-seq/sample_libs_info.tsv",
                            sep="\t", header=F)$V1

# List files / sample
R1files <- data.frame(
  sampleName = unlist(as.data.frame(
    sapply(samples.lc, function(x){
    rep(x,length(list.files(paste0(FASTQ.dir, x), pattern="R1")))})
    ), use.names=F),
  
  R1 = unlist(as.data.frame(
    sapply(samples.lc,function(x){list.files(paste0(FASTQ.dir,x),
                                                   pattern="R1")})
    ), use.names=F)
)
dim(R1files)  # 12: 4 lanes x 3 sample libs

for(i in 1:nrow(R1files)){
  cat(paste0("Checking R1 length distribution for: ", R1files[i,2], "\n"))
  temp.R1s <- readFastq(paste0(FASTQ.dir, R1files[i,1], "/", R1files[i,2]),
                        withIds=F)
  print(head(sread(temp.R1s), n=4))
  print(table(width(sread(temp.R1s))))
  rm(temp.R1s)
  cat("\n\n")
}
sessionInfo()
