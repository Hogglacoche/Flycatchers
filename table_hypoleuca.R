setwd("/home/amaniouloux/Documents/data/hypoleuca/")

pi_hypo_gc<-read.table("pied.gc_cons.pn_ps_by_gene.chr_info.txt", header = T)
pi_hypo_all<-read.table("pied.all_sites.pn_ps_by_gene.chr_info.txt", header = T)

library(dplyr)
library(stringr)

table_hypoleuca<-merge(pi_hypo_gc, pi_hypo_all, by = "gene", all = T)
table_hypoleuca$gene <- sapply(str_split(table_hypoleuca$gene, "_"), "[[", 1)
table_hypoleuca<-subset(table_hypoleuca, select= -chr.x)
table_hypoleuca<-table_hypoleuca %>% rename(chr = chr.y)

table_hypoleuca <- table_hypoleuca %>%
  rename_with(~ gsub("\\.x$", "_gc", .x), ends_with(".x"))
table_hypoleuca <- table_hypoleuca %>%
  rename_with(~ gsub("\\.y$", "_all", .x), ends_with(".y"))

load("~/Documents/data/TPM_values/STAR.RSEM.coll.pied.hybrid.RData")
TPM_brain<-data.coll.pied.hybrid[[1]]
TPM_heart<-data.coll.pied.hybrid[[2]]
TPM_kidney<-data.coll.pied.hybrid[[3]]
TPM_liver<-data.coll.pied.hybrid[[4]]
TPM_testis<-data.coll.pied.hybrid[[5]]

hypobrain<-select(TPM_brain, gene, pied.1, pied.2, pied.3, pied.4, pied.5)
table_hypoleuca<-merge(table_hypoleuca, hypobrain, by = "gene", all = T)
table_hypoleuca<-table_hypoleuca %>% rename(TPM_Brain_pied_1 = pied.1)
table_hypoleuca<-table_hypoleuca%>% rename(TPM_Brain_pied_2 = pied.2)
table_hypoleuca<-table_hypoleuca%>% rename(TPM_Brain_pied_3 = pied.3)
table_hypoleuca<-table_hypoleuca %>% rename(TPM_Brain_pied_4 = pied.4)
table_hypoleuca<-table_hypoleuca %>% rename(TPM_Brain_pied_5 = pied.5)
table_hypoleuca$TPM_Brain_pied_mean <- rowMeans(table_hypoleuca[, c("TPM_Brain_pied_1", "TPM_Brain_pied_2", "TPM_Brain_pied_3", "TPM_Brain_pied_4", "TPM_Brain_pied_5")])

hypoheart<-select(TPM_heart, gene, pied.1, pied.2, pied.3, pied.4, pied.5)
table_hypoleuca<-merge(table_hypoleuca, hypoheart, by = "gene", all = T)
table_hypoleuca<-table_hypoleuca %>% rename(TPM_Heart_pied_1 = pied.1)
table_hypoleuca<-table_hypoleuca%>% rename(TPM_Heart_pied_2 = pied.2)
table_hypoleuca<-table_hypoleuca%>% rename(TPM_Heart_pied_3 = pied.3)
table_hypoleuca<-table_hypoleuca %>% rename(TPM_Heart_pied_4 = pied.4)
table_hypoleuca<-table_hypoleuca %>% rename(TPM_Heart_pied_5 = pied.5)
table_hypoleuca$TPM_Heart_pied_mean <- rowMeans(table_hypoleuca[, c("TPM_Heart_pied_1", "TPM_Heart_pied_2", "TPM_Heart_pied_3", "TPM_Heart_pied_4", "TPM_Heart_pied_5")])

hypokidney<-select(TPM_kidney, gene, pied.1, pied.2, pied.3, pied.4, pied.5)
table_hypoleuca<-merge(table_hypoleuca, hypokidney, by = "gene", all = T)
table_hypoleuca<-table_hypoleuca %>% rename(TPM_Kidney_pied_1 = pied.1)
table_hypoleuca<-table_hypoleuca%>% rename(TPM_Kidney_pied_2 = pied.2)
table_hypoleuca<-table_hypoleuca%>% rename(TPM_Kidney_pied_3 = pied.3)
table_hypoleuca<-table_hypoleuca %>% rename(TPM_Kidney_pied_4 = pied.4)
table_hypoleuca<-table_hypoleuca %>% rename(TPM_Kidney_pied_5 = pied.5)
table_hypoleuca$TPM_Kidney_pied_mean <- rowMeans(table_hypoleuca[, c("TPM_Kidney_pied_1", "TPM_Kidney_pied_2", "TPM_Kidney_pied_3", "TPM_Kidney_pied_4", "TPM_Kidney_pied_5")])

hypoliver<-select(TPM_liver, gene, pied.1, pied.2, pied.3, pied.4, pied.5)
table_hypoleuca<-merge(table_hypoleuca, hypoliver, by = "gene", all = T)
table_hypoleuca<-table_hypoleuca %>% rename(TPM_Liver_pied_1 = pied.1)
table_hypoleuca<-table_hypoleuca%>% rename(TPM_Liver_pied_2 = pied.2)
table_hypoleuca<-table_hypoleuca%>% rename(TPM_Liver_pied_3 = pied.3)
table_hypoleuca<-table_hypoleuca %>% rename(TPM_Liver_pied_4 = pied.4)
table_hypoleuca<-table_hypoleuca %>% rename(TPM_Liver_pied_5 = pied.5)
table_hypoleuca$TPM_Liver_pied_mean <- rowMeans(table_hypoleuca[, c("TPM_Liver_pied_1", "TPM_Liver_pied_2", "TPM_Liver_pied_3", "TPM_Liver_pied_4", "TPM_Liver_pied_5")])

hypotestis<-select(TPM_testis, gene, pied.1, pied.2, pied.3, pied.4, pied.5)
table_hypoleuca<-merge(table_hypoleuca, hypotestis, by = "gene", all = T)
table_hypoleuca<-table_hypoleuca %>% rename(TPM_Testis_pied_1 = pied.1)
table_hypoleuca<-table_hypoleuca %>% rename(TPM_Testis_pied_2 = pied.2)
table_hypoleuca<-table_hypoleuca %>% rename(TPM_Testis_pied_3 = pied.3)
table_hypoleuca<-table_hypoleuca %>% rename(TPM_Testis_pied_4 = pied.4)
table_hypoleuca<-table_hypoleuca %>% rename(TPM_Testis_pied_5 = pied.5)
table_hypoleuca$TPM_Testis_pied_mean <- rowMeans(table_hypoleuca[, c("TPM_Testis_pied_1", "TPM_Testis_pied_2", "TPM_Testis_pied_3", "TPM_Testis_pied_4", "TPM_Testis_pied_5")])

write.table(table_hypoleuca, "table_hypoleuca.txt", sep = "\t", row.names = FALSE)
getwd()
