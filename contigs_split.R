getwd()
working_path <- getwd()
folder <- "export_splites"
files <- list.files(path = "export_splites", pattern = "*.txt")
data <- read.table(file = paste(working_path,folder,files, sep = "/"), header = TRUE)
abund <- data[,-1]
contig <- as.data.frame(data[,1]); colnames(contig) <- "contig"

for (i in 1:length(abund)){
  num <- paste("col", i, sep = "")
  assign(num, cbind(contig, abund[,i]))
  write.table(cbind(contig, abund[,i]), paste(folder,paste('abund', i, 'txt', sep = '.'), sep = "/"), sep = "\t", quote = F, row.names = F, col.names = F)
}
