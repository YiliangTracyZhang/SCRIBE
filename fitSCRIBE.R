options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

expression.file <- as.character(args[1])
batch.file <- as.character(args[2])
bio.file <- as.character(args[3])
output.file <- as.character(args[4])

expression.matrix <- as.matrix(read.table(expression.file, head=F))
batch.ind <- read.table(batch.file, head=F)
batch.ind <- batch.ind$V1
bio.ind <- read.table(bio.file, head=F)
bio.ind <- bio.ind$V1

fitted.model <- fitSCRIBE(expression.matrix, batch.ind, bio.ind)
Y.imputed <- fitted.model$Y_impute
write.table(Y.imputed, output.file, quote = F, col.names = F, row.names = F)