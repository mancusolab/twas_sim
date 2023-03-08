library(susieR)

args <- commandArgs(trailingOnly=TRUE)
z_path <- args[1]
y_path <- args[2]
out_path <- args[3]

Z <- as.matrix(read.table(z_path, header = FALSE, dec = "."))
y <- as.matrix(read.table(y_path, header = FALSE, dec = "."))

res <- susie(Z, y, L=5)
g_coef <- coef(res)[2:length(coef(res)]

write.table(g_coef, file = paste0("twas_sim_loci", IDX, ".external_coef.txt"), sep = "\t", row.names = TRUE, col.names = FALSE)

