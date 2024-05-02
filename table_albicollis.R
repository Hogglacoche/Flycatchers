setwd("/home/amaniouloux/Documents/data/albicollis")

bbp_coll_gc<- read.table("bpp_gc_cons.all_genes.coll_results_summ.SW_normed.gene_ids_split.chr_info.chr_genes.txt", header = T)
bbp_coll_all_sites <- read.table("bpp_out.all_sites.all_genes.coll_counts.gene_id_split.chr_info.chr_genes.dnd", header = T)
pi_coll_gc<- read.table("coll.gc_cons.pn_ps_by_gene.chr_info.txt", header = T)
pi_coll_all_sites <- read.table("coll.all_sites.pn_ps_by_gene.chr_info.txt", header = T)
coll_32_gc <- read.table("coll_32.all_genes.coll_only_stats.no_neg_counts.normd.min_ex200bp.chr_info.dnd", header = T)
coll_32_all_sites <- read.table("coll_32.all_genes.all_sites.coll_only_stats.no_neg_counts.min_ex200bp.chr_info.dnd", header = T)

library(dplyr)
library(stringr)


table_albicollis<-merge(bbp_coll_gc, bbp_coll_all_sites, by = "gene", all = T)
table_albicollis$gene <- sapply(str_split(table_albicollis$gene, "_"), "[[", 1)
table_albicollis<-subset(table_albicollis, select= -chr.x)
table_albicollis<-subset(table_albicollis, select= -chr.y)
table_albicollis <- table_albicollis %>%
  rename_with(~ gsub("\\.x$", "_gc_L", .x), ends_with(".x"))
table_albicollis <- table_albicollis %>%
  rename_with(~ gsub("\\.y$", "_all_L", .x), ends_with(".y"))
table_albicollis<-table_albicollis %>% rename(dn_gene_gc_L = dn_gene)
table_albicollis<-table_albicollis %>% rename(ds_gene_gc_L = ds_gene)
table_albicollis<-table_albicollis %>% rename(dn_ds_gene_gc_L = dn_ds_gene)
table_albicollis$dn_gene_all_L<-(table_albicollis$dn_tot_all_L/table_albicollis$dn_norm_tot_all_L)
table_albicollis$ds_gene_all_L<-(table_albicollis$ds_tot_all_L/table_albicollis$ds_norm_tot_all_L)
table_albicollis$dn_ds_gene_all_L<-(table_albicollis$dn_gene_all_L/table_albicollis$ds_gene_all_L)

coll_32_gc$gene <- sapply(str_split(coll_32_gc$gene, "_"), "[[", 1)
table_albicollis<-merge(table_albicollis, coll_32_gc, by = "gene", all=T)
coll_32_all_sites$gene<- sapply(str_split(coll_32_all_sites$gene, "_"), "[[", 1)
table_albicollis<-merge(table_albicollis, coll_32_all_sites, by = "gene", all=T)
table_albicollis <- table_albicollis %>%
  rename_with(~ gsub("\\.x$", "_gc_S", .x), ends_with(".x"))
table_albicollis <- table_albicollis %>%
  rename_with(~ gsub("\\.y$", "_all_S", .x), ends_with(".y"))
table_albicollis$gene <- sapply(str_split(table_albicollis$gene, "_"), "[[", 1)
table_albicollis<-table_albicollis %>% rename(dn_gene_gc_S = dn_gene)
table_albicollis<-table_albicollis %>% rename(ds_gene_gc_S = ds_gene)
table_albicollis<-table_albicollis %>% rename(dn_ds_gene_gc_S = dn_ds_gene)
table_albicollis<-subset(table_albicollis, select= -chr_all_S)
table_albicollis<-subset(table_albicollis, select= -chr_gc_S)
table_albicollis$dn_gene_all_S<-(table_albicollis$dn_tot_all_S/table_albicollis$dn_norm_tot_all_S)
table_albicollis$ds_gene_all_S<-(table_albicollis$ds_tot_all_S/table_albicollis$ds_norm_tot_all_S)
table_albicollis$dn_ds_gene_all_S<-(table_albicollis$dn_gene_all_S/table_albicollis$ds_gene_all_S)

pi_coll_gc$gene <- sapply(str_split(pi_coll_gc$gene, "_"), "[[", 1)
pi_coll_all_sites$gene <- sapply(str_split(pi_coll_all_sites$gene, "_"), "[[", 1)
table_albicollis<-merge(table_albicollis, pi_coll_gc, by= "gene", all = T)
table_albicollis<-merge(table_albicollis, pi_coll_all_sites, by = "gene", all = T)
table_albicollis<-subset(table_albicollis, select= -chr.x)
table_albicollis<-table_albicollis %>% rename(chr = chr.y)
table_albicollis <- table_albicollis %>%
  rename_with(~ gsub("\\.x$", "_gc", .x), ends_with(".x"))
table_albicollis <- table_albicollis %>%
  rename_with(~ gsub("\\.y$", "_all", .x), ends_with(".y"))

load("~/Documents/data/TPM_values/STAR.RSEM.coll.pied.hybrid.RData")
TPM_brain<-data.coll.pied.hybrid[[1]]
TPM_heart<-data.coll.pied.hybrid[[2]]
TPM_kidney<-data.coll.pied.hybrid[[3]]
TPM_liver<-data.coll.pied.hybrid[[4]]
TPM_testis<-data.coll.pied.hybrid[[5]]

colbrain<-select(TPM_brain, gene, coll.1, coll.2, coll.3, coll.4, coll.5)
table_albicollis<-merge(table_albicollis, colbrain, by = "gene", all = T)
table_albicollis<-table_albicollis %>% rename(TPM_Brain_coll_1 = coll.1)
table_albicollis<-table_albicollis %>% rename(TPM_Brain_coll_2 = coll.2)
table_albicollis<-table_albicollis %>% rename(TPM_Brain_coll_3 = coll.3)
table_albicollis<-table_albicollis %>% rename(TPM_Brain_coll_4 = coll.4)
table_albicollis<-table_albicollis %>% rename(TPM_Brain_coll_5 = coll.5)
table_albicollis$TPM_Brain_coll_mean<-rowMeans(table_albicollis[, c("TPM_Brain_coll_1", "TPM_Brain_coll_2", "TPM_Brain_coll_3", "TPM_Brain_coll_4", "TPM_Brain_coll_5")])

colheart<-select(TPM_heart, gene, coll.1, coll.2, coll.3, coll.4, coll.5)
table_albicollis<-merge(table_albicollis, colheart, by = "gene", all = T)
table_albicollis<-table_albicollis %>% rename(TPM_Heart_coll_1 = coll.1)
table_albicollis<-table_albicollis %>% rename(TPM_Heart_coll_2 = coll.2)
table_albicollis<-table_albicollis %>% rename(TPM_Heart_coll_3 = coll.3)
table_albicollis<-table_albicollis %>% rename(TPM_Heart_coll_4 = coll.4)
table_albicollis<-table_albicollis %>% rename(TPM_Heart_coll_5 = coll.5)
table_albicollis$TPM_Heart_coll_mean<-rowMeans(table_albicollis[, c("TPM_Heart_coll_1", "TPM_Heart_coll_2", "TPM_Heart_coll_3", "TPM_Heart_coll_4", "TPM_Heart_coll_5")])

colkidney<-select(TPM_kidney, gene, coll.1, coll.2, coll.3, coll.4, coll.5)
table_albicollis<-merge(table_albicollis, colkidney, by = "gene", all = T)
table_albicollis<-table_albicollis %>% rename(TPM_Kidney_coll_1 = coll.1)
table_albicollis<-table_albicollis %>% rename(TPM_Kidney_coll_2 = coll.2)
table_albicollis<-table_albicollis %>% rename(TPM_Kidney_coll_3 = coll.3)
table_albicollis<-table_albicollis %>% rename(TPM_Kidney_coll_4 = coll.4)
table_albicollis<-table_albicollis %>% rename(TPM_Kidney_coll_5 = coll.5)
table_albicollis$TPM_Kidney_coll_mean <- rowMeans(table_albicollis[, c("TPM_Kidney_coll_1", "TPM_Kidney_coll_2", "TPM_Kidney_coll_3", "TPM_Kidney_coll_4", "TPM_Kidney_coll_5")])

colliver<-select(TPM_liver, gene, coll.1, coll.2, coll.3, coll.4, coll.5)
table_albicollis<-merge(table_albicollis, colliver, by = "gene", all = T)
table_albicollis<-table_albicollis %>% rename(TPM_Liver_coll_1 = coll.1)
table_albicollis<-table_albicollis %>% rename(TPM_Liver_coll_2 = coll.2)
table_albicollis<-table_albicollis %>% rename(TPM_Liver_coll_3 = coll.3)
table_albicollis<-table_albicollis %>% rename(TPM_Liver_coll_4 = coll.4)
table_albicollis<-table_albicollis %>% rename(TPM_Liver_coll_5 = coll.5)
table_albicollis$TPM_Liver_coll_mean <- rowMeans(table_albicollis[, c("TPM_Liver_coll_1", "TPM_Liver_coll_2", "TPM_Liver_coll_3", "TPM_Liver_coll_4", "TPM_Liver_coll_5")])

coltestis<-select(TPM_testis, gene, coll.1, coll.2, coll.3, coll.4, coll.5)
table_albicollis<-merge(table_albicollis, coltestis, by = "gene", all = T)
table_albicollis<-table_albicollis %>% rename(TPM_Testis_coll_1 = coll.1)
table_albicollis<-table_albicollis %>% rename(TPM_Testis_coll_2 = coll.2)
table_albicollis<-table_albicollis %>% rename(TPM_Testis_coll_3 = coll.3)
table_albicollis<-table_albicollis %>% rename(TPM_Testis_coll_4 = coll.4)
table_albicollis<-table_albicollis %>% rename(TPM_Testis_coll_5 = coll.5)
table_albicollis$TPM_Testis_coll_mean <- rowMeans(table_albicollis[, c("TPM_Testis_coll_1", "TPM_Testis_coll_2", "TPM_Testis_coll_3", "TPM_Testis_coll_4", "TPM_Testis_coll_5")])
setwd("/home/amaniouloux/Documents/data/dataframe")
size<- read.table("chr_category.txt", header = T)
size<-size %>% rename(chr = chromosome)
table_albicollis<-merge(table_albicollis, size, by = "chr", all =T)

setwd("/home/amaniouloux/Documents/data/Rstudio/Recombinaison/")
recomb_coll<-read.table("coll.ld_recom_converted.gene_recom_no_cpg_filt.recom_bin_inf.recom_cons_inf.bed", header = T)
recomb_coll<-recomb_coll %>% rename(gene = gene_id)
recomb_coll$gene <- sapply(str_split(recomb_coll$gene, "_"), "[[", 1)
table_albicollis<-merge(table_albicollis, recomb_coll, by = "gene", all = T)

write.table(table_albicollis, file = "table_albicollis.txt", sep = "\t", row.names = FALSE)
getwd()
