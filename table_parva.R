setwd("/home/amaniouloux/Documents/data/parva/")

pi_parva_gc<-read.table("par.gc_cons.pn_ps_by_gene.chr_info.txt", header = T)
pi_parva_all<-read.table("par.all_sites.pn_ps_by_gene.chr_info.txt", header = T)

library(dplyr)
library(stringr)

table_parva<-merge(pi_parva_gc, pi_parva_all, by = "gene", all = T)
table_parva<-subset(table_parva, select= -chr.x)
table_parva<-table_parva %>% rename(chr = chr.y)

table_parva <- table_parva %>%
  rename_with(~ gsub("\\.x$", "_gc", .x), ends_with(".x"))
table_parva <- table_parva %>%
  rename_with(~ gsub("\\.y$", "_all", .x), ends_with(".y"))
table_parva$gene <- sapply(str_split(table_parva$gene, "_"), "[[", 1)


setwd("/home/amaniouloux/Documents/data/dataframe")
size<- read.table("chr_category.txt", header = T)
size<-size %>% rename(chr = chromosome)
table_parva<-merge(table_parva, size, by = "chr", all =T)

write.table(table_parva, file="table_parva.txt", sep= "\t", row.names = F)

