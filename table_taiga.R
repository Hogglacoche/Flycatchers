setwd("/home/amaniouloux/Documents/data/taiga/")

pi_taiga_gc<-read.table("taig.gc_cons.pn_ps_by_gene.chr_info.txt", header = T)
pi_taiga_all<-read.table("taig.all_sites.pn_ps_by_gene.chr_info.txt", header = T)
taiga_32_all<-read.table("taig_32.all_genes.all_sites.taig_only_stats.no_neg_counts.min_ex200bp.chr_info.dnd", header = T)
taiga_32_gc<-read.table("taig_32.all_genes.taig_only_stats.no_neg.normd.min_ex200bp.chr_info.dnd", header = T)

library(dplyr)
library(stringr)

taiga_32_gc$gene <- sapply(str_split(taiga_32_gc$gene, "_"), "[[", 1)
taiga_32_all$gene <- sapply(str_split(taiga_32_all$gene, "_"), "[[", 1)
table_taiga<-merge(taiga_32_gc, taiga_32_all, by = "gene", all = T)
table_taiga<-subset(table_taiga, select= -chr.x)
table_taiga<-subset(table_taiga, select= -chr.y)
table_taiga<- table_taiga %>%
  rename_with(~ gsub("\\.x$", "_gc_S", .x), ends_with(".x"))
table_taiga <- table_taiga %>%
  rename_with(~ gsub("\\.y$", "_all_S", .x), ends_with(".y"))

table_taiga$dn_gene_all_S<-(table_taiga$dn_tot_all_S/table_taiga$dn_norm_tot_all_S)
table_taiga$ds_gene_all_S<-(table_taiga$ds_tot_all_S/table_taiga$ds_norm_tot_all_S)
table_taiga$dn_ds_gene_all_S<-(table_taiga$dn_gene_all_S/table_taiga$ds_gene_all_S)

table_taiga<-table_taiga %>% rename(dn_gene_gc_S = dn_gene)
table_taiga<-table_taiga %>% rename(ds_gene_gc_S = ds_gene)
table_taiga<-table_taiga %>% rename(dn_ds_gene_gc_S = dn_ds_gene)

pi_taiga_gc$gene <- sapply(str_split(pi_taiga_gc$gene, "_"), "[[", 1)
pi_taiga_all$gene <- sapply(str_split(pi_taiga_all$gene, "_"), "[[", 1)
table_taiga<-merge(table_taiga, pi_taiga_gc, by= "gene", all = T)
table_taiga<-merge(table_taiga, pi_taiga_all, by = "gene", all = T)
table_taiga<-subset(table_taiga, select= -chr.x)
table_taiga<-table_taiga %>% rename(chr = chr.y)

table_taiga<- table_taiga %>%
  rename_with(~ gsub("\\.x$", "_gc", .x), ends_with(".x"))
table_taiga <- table_taiga %>%
  rename_with(~ gsub("\\.y$", "_all", .x), ends_with(".y"))

setwd("/home/amaniouloux/Documents/data/dataframe")
size<- read.table("chr_category.txt", header = T)
size<-size %>% rename(chr = chromosome)
table_taiga<-merge(table_taiga, size, by = "chr", all =T)

setwd("/home/amaniouloux/Documents/data/Rstudio/Recombinaison/")
recomb_taig<-read.table("taig.ld_recom_converted.gene_recom_no_cpg_filt.recom_bin_inf.recom_cons_inf.bed", header = T)
recomb_taig<-recomb_taig %>% rename(gene = gene_id)
recomb_taig$gene <- sapply(str_split(recomb_taig$gene, "_"), "[[", 1)
table_taiga<-merge(table_taiga, recomb_taig, by = "gene", all = T)

write.table(table_taiga, "table_taiga.txt", sep = "\t", row.names = FALSE)
