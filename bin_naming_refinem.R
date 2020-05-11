getwd()

files <- list.files(pattern = "*.txt")
data2 <- lapply(files, read.table, header=FALSE)
for (i in 1:length(data2)){data2[[i]]<-cbind(data2[[i]],files[i])}
data_rbind <- do.call("rbind", data2) 
colnames(data_rbind)<-c("contigs", "bins")

#write.table(data_rbind, file="bins.txt", sep = "\t", quote = F, row.names = F)
library(tidyverse)
data.temp <- separate(data = data_rbind, col = contigs, into = c("null", "contigs"), sep = ">")
data.temp <- separate(data = data.temp, col = bins, into = c("bin", "num","x","y","z"))
data.temp <- data.temp %>% select(contigs, num)
data.temp$num <- str_c("bin",data.temp$num)

write.table(data.temp, file="bins.txt", sep = "\t", quote = F, row.names = F, col.names = F)
