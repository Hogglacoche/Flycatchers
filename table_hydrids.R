table_hypoleuca<-read.table("table_hypoleuca.txt", header = T)
table_hypoleuca<- table_hypoleuca[, c('gene','chr')]
size<-read.table("chr_category.txt", header = T)
size<-size %>% rename( chr = chromosome)

load("~/Documents/data/TPM_values/STAR.RSEM.coll.pied.hybrid.RData")
TPM_brain<-data.coll.pied.hybrid[[1]]
TPM_heart<-data.coll.pied.hybrid[[2]]
TPM_kidney<-data.coll.pied.hybrid[[3]]
TPM_liver<-data.coll.pied.hybrid[[4]]
TPM_testis<-data.coll.pied.hybrid[[5]]


hybridbrain<-select(TPM_brain, gene, hybrid.1, hybrid.2, hybrid.3)
hybridbrain<-hybridbrain %>% rename(TPM_Brain_hybrid_1 = hybrid.1)
hybridbrain<-hybridbrain%>% rename(TPM_Brain_hybrid_2 = hybrid.2)
hybridbrain<-hybridbrain%>% rename(TPM_Brain_hybrid_3 = hybrid.3)
hybridbrain$TPM_Brain_hybrid_mean <- rowMeans(hybridbrain[, c("TPM_Brain_hybrid_1", "TPM_Brain_hybrid_2", "TPM_Brain_hybrid_3")])

hybridheart<-select(TPM_heart, gene, hybrid.1, hybrid.2, hybrid.3)
table_hybrid<-merge(hybridbrain, hybridheart, by = "gene", all = T)
table_hybrid<-table_hybrid %>% rename(TPM_Heart_hybrid_1 = hybrid.1)
table_hybrid<-table_hybrid%>% rename(TPM_Heart_hybrid_2 = hybrid.2)
table_hybrid<-table_hybrid%>% rename(TPM_Heart_hybrid_3 = hybrid.3)
table_hybrid$TPM_Heart_hybrid_mean <- rowMeans(table_hybrid[, c("TPM_Heart_hybrid_1", "TPM_Heart_hybrid_2", "TPM_Heart_hybrid_3")])

hybridkidney<-select(TPM_kidney, gene, hybrid.1, hybrid.2, hybrid.3)
table_hybrid<-merge(table_hybrid, hybridkidney, by = "gene", all = T)
table_hybrid<-table_hybrid %>% rename(TPM_Kidney_hybrid_1 = hybrid.1)
table_hybrid<-table_hybrid%>% rename(TPM_Kidney_hybrid_2 = hybrid.2)
table_hybrid<-table_hybrid%>% rename(TPM_Kidney_hybrid_3 = hybrid.3)
table_hybrid$TPM_Kidney_hybrid_mean <- rowMeans(table_hybrid[, c("TPM_Kidney_hybrid_1", "TPM_Kidney_hybrid_2", "TPM_Kidney_hybrid_3")])

hybridliver<-select(TPM_liver, gene, hybrid.1, hybrid.2, hybrid.3)
table_hybrid<-merge(table_hybrid, hybridliver, by = "gene", all = T)
table_hybrid<-table_hybrid %>% rename(TPM_Liver_hybrid_1 = hybrid.1)
table_hybrid<-table_hybrid%>% rename(TPM_Liver_hybrid_2 = hybrid.2)
table_hybrid<-table_hybrid%>% rename(TPM_Liver_hybrid_3 = hybrid.3)
table_hybrid$TPM_Liver_hybrid_mean <- rowMeans(table_hybrid[, c("TPM_Liver_hybrid_1", "TPM_Liver_hybrid_2", "TPM_Liver_hybrid_3")])

hybridtestis<-select(TPM_testis, gene, hybrid.1, hybrid.2, hybrid.3)
table_hybrid<-merge(table_hybrid, hybridtestis, by = "gene", all = T)
table_hybrid<-table_hybrid %>% rename(TPM_Testis_hybrid_1 = hybrid.1)
table_hybrid<-table_hybrid %>% rename(TPM_Testis_hybrid_2 = hybrid.2)
table_hybrid<-table_hybrid %>% rename(TPM_Testis_hybrid_3 = hybrid.3)
table_hybrid$TPM_Testis_hybrid_mean <- rowMeans(table_hybrid[, c("TPM_Testis_hybrid_1", "TPM_Testis_hybrid_2", "TPM_Testis_hybrid_3")])


table_hybrid<-merge(table_hybrid, table_hypoleuca, by = "gene", all = T)

table_hybrid<-merge(table_hybrid, size, by ="chr", all = T)

write.table(table_hybrid, "table_hybrids.txt", sep = "\t", row.names = FALSE)
