odir <- '/scratch1/xwang505/TWAS/output'
setwd(odir)

library(susieR)
args <- commandArgs(trailingOnly=TRUE)
IDX <- as.numeric(args[1])

# how to distinguish loci{i}
Z <- read.table(paste0("twas_sim_loci", IDX, '.eqtl.genotype.txt.gz'), header = TRUE, sep = "", dec = ".") 
y <- read.table(paste0("twas_sim_loci", IDX, '.eqtl.gexpr.txt.gz'), header = TRUE, sep = "", dec = ".")  

n <- nrow(Z)
p <- ncol(Z)

res <- susie(Z, y, L=10)
coef <- coef(res)

# how to write external_coef{i}.txt.gz
write.table(coef, file = paste0("twas_sim_loci", IDX, ".external_coef.txt"), sep = "\t", row.names = TRUE, col.names = FALSE)

