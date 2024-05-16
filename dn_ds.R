setwd("/home/amaniouloux/Documents/data/dataframe/")
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(gtools)
library(corrplot)
library(stringr)

calculer_dn_ds <- function(dn_tot, dn_norm_tot, ds_tot, ds_norm_tot) {
  return((dn_tot / dn_norm_tot) / (ds_tot / ds_norm_tot))
}


## coll_dn_ds_gc_L

table_albicollis <- read.table("table_albicollis.txt", header = T) 
count_dnds_gc_L_coll<- aggregate(dn_ds_gene_gc_L ~ chr, data = table_albicollis, FUN = function(x) sum(!is.na(x)))
count_dnds_gc_L_coll <- data.frame(chr = count_dnds_gc_L_coll$chr, count_dnds_gc_L_coll = count_dnds_gc_L_coll$dn_ds_gene_gc_L)
coll_dn_ds_gc_L <- aggregate(cbind(dn_tot_gc_L, dn_norm_tot_gc_L, ds_tot_gc_L, ds_norm_tot_gc_L) ~ chr, table_albicollis, mean)
coll_dn_ds_gc_L$coll_dn_ds_gc_L <- with(coll_dn_ds_gc_L, calculer_dn_ds(dn_tot_gc_L, dn_norm_tot_gc_L, ds_tot_gc_L, ds_norm_tot_gc_L))

## coll_dn_ds_gc_S

table_albicollis <- read.table("table_albicollis.txt", header = T) 
count_dnds_gc_S_coll<- aggregate(dn_ds_gene_gc_S ~ chr, data = table_albicollis, FUN = function(x) sum(!is.na(x)))
count_dnds_gc_S_coll <- data.frame(chr = count_dnds_gc_S_coll$chr, count_dnds_gc_S_coll = count_dnds_gc_S_coll$dn_ds_gene_gc_S)
coll_dn_ds_gc_S <- aggregate(cbind(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S) ~ chr, table_albicollis, mean)
coll_dn_ds_gc_S$coll_dn_ds_gc_S <- with(coll_dn_ds_gc_S, calculer_dn_ds(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S))

## coll_dn_ds_all_L

table_albicollis <- read.table("table_albicollis.txt", header = T) 
count_dnds_all_L_coll<- aggregate(dn_ds_gene_all_L ~ chr, data = table_albicollis, FUN = function(x) sum(!is.na(x)))
count_dnds_all_L_coll <- data.frame(chr = count_dnds_all_L_coll$chr, count_dnds_all_L_coll = count_dnds_all_L_coll$dn_ds_gene_all_L)
coll_dn_ds_all_L <- aggregate(cbind(dn_tot_all_L, dn_norm_tot_all_L, ds_tot_all_L, ds_norm_tot_all_L) ~ chr, table_albicollis, mean)
coll_dn_ds_all_L$coll_dn_ds_all_L <- with(coll_dn_ds_all_L, calculer_dn_ds(dn_tot_all_L, dn_norm_tot_all_L, ds_tot_all_L, ds_norm_tot_all_L))

## coll_dn_ds_all_S

table_albicollis <- read.table("table_albicollis.txt", header = T) 
count_dnds_all_S_coll<- aggregate(dn_ds_gene_all_S ~ chr, data = table_albicollis, FUN = function(x) sum(!is.na(x)))
count_dnds_all_S_coll<- data.frame(chr = count_dnds_all_S_coll$chr, count_dnds_all_S_coll = count_dnds_all_S_coll$dn_ds_gene_all_S)
coll_dn_ds_all_S <- aggregate(cbind(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S) ~ chr, table_albicollis, mean)
coll_dn_ds_all_S$coll_dn_ds_all_S <- with(coll_dn_ds_all_S, calculer_dn_ds(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S))

## taiga_dn_ds_gc_S

table_taiga <- read.table("table_taiga.txt", header = T) 
count_dnds_gc_S_taiga<- aggregate(dn_ds_gene_gc_S ~ chr, data = table_taiga, FUN = function(x) sum(!is.na(x)))
count_dnds_gc_S_taiga<- data.frame(chr = count_dnds_gc_S_taiga$chr, count_dnds_gc_S_taiga = count_dnds_gc_S_taiga$dn_ds_gene_gc_S)
taiga_dn_ds_gc_S <- aggregate(cbind(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S) ~ chr, table_taiga, mean)
taiga_dn_ds_gc_S$taiga_dn_ds_gc_S <- with(taiga_dn_ds_gc_S, calculer_dn_ds(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S))

## taiga_dn_ds_all.S

table_taiga <- read.table("table_taiga.txt", header = T) 
count_dnds_all_S_taiga<- aggregate(dn_ds_gene_all_S ~ chr, data = table_taiga, FUN = function(x) sum(!is.na(x)))
count_dnds_all_S_taiga<- data.frame(chr = count_dnds_all_S_taiga$chr, count_dnds_all_S_taiga = count_dnds_all_S_taiga$dn_ds_gene_all_S)
taiga_dn_ds_all_S <- aggregate(cbind(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S) ~ chr, table_taiga, mean)
taiga_dn_ds_all_S$taiga_dn_ds_all_S<- with(taiga_dn_ds_all_S, calculer_dn_ds(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S))

## merge 
count_dnds <- count_dnds_all_L_coll %>%
  left_join(count_dnds_all_S_coll, by = "chr") %>%
  left_join(count_dnds_gc_L_coll, by = "chr") %>%
  left_join(count_dnds_gc_S_coll, by = "chr") %>%
  left_join(count_dnds_all_S_taiga, by = "chr") %>%
  left_join(count_dnds_gc_S_taiga, by = "chr")

merged_dn_ds <- coll_dn_ds_all_L %>%
  left_join(coll_dn_ds_all_S, by = "chr") %>%
  left_join(coll_dn_ds_gc_L, by = "chr") %>%
  left_join(coll_dn_ds_gc_S, by = "chr") %>%
  left_join(taiga_dn_ds_all_S, by = "chr") %>%
  left_join(taiga_dn_ds_gc_S, by = "chr")

merged_dn_ds <- merged_dn_ds %>%
  rename_with(~ gsub("\\.x$", "_coll", .x), ends_with(".x"))
table_albicollis <- table_albicollis %>%
  rename_with(~ gsub("\\.y$", "_taiga", .x), ends_with(".y"))

merged_dn_ds <- merged_dn_ds[mixedorder(merged_dn_ds$chr), ]
merged_dn_ds <- merged_dn_ds[merged_dn_ds$chr != "Chr22", ]
merged_dn_ds <- merged_dn_ds[merged_dn_ds$chr != "Chr23", ]
merged_dn_ds <- merged_dn_ds[merged_dn_ds$chr != "Chr24", ]
merged_dn_ds <- merged_dn_ds[merged_dn_ds$chr != "Chr25", ]
merged_dn_ds <- merged_dn_ds[merged_dn_ds$chr != "Chr26", ]
merged_dn_ds <- merged_dn_ds[merged_dn_ds$chr != "Chr27", ]
merged_dn_ds <- merged_dn_ds[merged_dn_ds$chr != "ChrLGE22", ]
merged_dn_ds <- merged_dn_ds[merged_dn_ds$chr != "Chr28", ]


################################################################################################

## plot

merged_dn_ds$chr <- factor(merged_dn_ds$chr, levels = mixedsort(merged_dn_ds$chr))

ggplot() +
  geom_point(data = merged_dn_ds, aes(x = chr, y = coll_dn_ds_gc_L, color = "dn/ds GC L"), shape = 16, size = 2) +
  geom_point(data = merged_dn_ds, aes(x = chr, y = coll_dn_ds_gc_S, color = "dn/ds GC S"), shape = 16, size = 2) +
  geom_point(data = merged_dn_ds, aes(x = chr, y = coll_dn_ds_all_L, color = "dn/ds ALL L"), shape = 17, size = 2) +
  geom_point(data = merged_dn_ds, aes(x = chr, y = coll_dn_ds_all_S, color = "dn/ds ALL S"), shape = 17, size = 2) +
  labs(x = "chromosome", y = "dn/ds", color = "Conditions") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

ggplot() +
  geom_point(data = merged_dn_ds, aes(x = chr, y = taiga_dn_ds_gc_S, color = "Taiga dn/ds GC S")) +
  geom_point(data = merged_dn_ds, aes(x = chr, y = taiga_dn_ds_all_S, color = "Taiga dn/ds ALL S")) +
  labs(x = "chromosome", y = "dn/ds", color = "Conditions")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

################################################################################################

# correlation matrix dn/ds chr

matrix_chr <- merged_dn_ds[, c('chr', 'coll_dn_ds_gc_L', 'coll_dn_ds_all_L','coll_dn_ds_gc_S', 'coll_dn_ds_all_S','taiga_dn_ds_gc_S','taiga_dn_ds_all_S')]
matrix_chr <- matrix_chr %>%
  rename(dn_ds_gc_L_coll = coll_dn_ds_gc_L,
         dn_ds_all_L_coll = coll_dn_ds_all_L,
         dn_ds_gc_S_coll = coll_dn_ds_gc_S,
         dn_ds_all_S_coll = coll_dn_ds_all_S,
         dn_ds_gc_S_taiga = taiga_dn_ds_gc_S,
         dn_ds_all_S_taiga = taiga_dn_ds_all_S)

is.numeric(matrix_chr)
matrix_chr<- sapply(matrix_chr, as.numeric)

correlation_matrix <- cor(matrix_chr[, -1], method = "spearman")
corrplot(correlation_matrix, method="circle", type="upper", tl.col="black", tl.srt=45)

###############################################################################################

# correlation matrix gene

setwd("/home/amaniouloux/Documents/data/dataframe/")

table_albicollis <- read.table("table_albicollis.txt", header = T) 
table_taiga <- read.table("table_taiga.txt", header = T)
table_albicollis$gene <- sapply(str_split(table_albicollis$gene, "_"), "[[", 1)
table_taiga$gene <- sapply(str_split(table_taiga$gene, "_"), "[[", 1)

dnds_table_albicollis<- table_albicollis[, c('gene','chr', 'dn_ds_gene_gc_S', 'dn_ds_gene_all_S', 'dn_ds_gene_gc_L', 'dn_ds_gene_all_L')]
dnds_table_albicollis<-na.omit(dnds_table_albicollis)
dnds_table_taiga<- table_taiga[, c('gene','chr', 'dn_ds_gene_all_S', 'dn_ds_gene_gc_S')]
dnds_table_taiga<-na.omit(dnds_table_taiga)

dnds_table<-merge(dnds_table_albicollis, dnds_table_taiga, by = "gene", all = T)
dnds_table <- dnds_table %>%
  rename_with(~ gsub("\\.x$", "_coll", .x), ends_with(".x"))
dnds_table <- dnds_table %>%
  rename_with(~ gsub("\\.y$", "_taiga", .x), ends_with(".y"))
dnds_table<-subset(dnds_table, select= -chr_taiga)
dnds_table<-dnds_table %>% rename(chr = chr_coll)
dnds_table<-dnds_table %>% rename( dn_ds_gene_gc_L_coll = dn_ds_gene_gc_L)
dnds_table<-dnds_table %>% rename(dn_ds_gene_all_L_coll = dn_ds_gene_all_L)
dnds_table <- dnds_table %>% filter(gene != "ENSFALG00000027500")
dnds_table <- dnds_table %>% filter(gene != "ENSFALG00000013927")
dnds_table <- dnds_table %>% filter(gene != "ENSFALG00000012551")
dnds_table <- dnds_table %>% filter(gene != "ENSFALG00000012551")

dnds_table<-subset(dnds_table, select= -gene)
dnds_table<-subset(dnds_table, select= -chr)
dnds_table<-na.omit(dnds_table)
is.numeric(dnds_table)
dnds_table<- sapply(dnds_table, as.numeric)

correlation_matrix<-cor(dnds_table, method="spearman")
correlation_matrix
corrplot(correlation_matrix, method="circle", type="upper", tl.col="black", tl.srt=45)


########################################################################################################

# dn/ds groupe Category macro micro intermediate

table_albicollis <- read.table("table_albicollis.txt", header = T) 
table_taiga <- read.table("table_taiga.txt", header = T)

calculer_dn_ds <- function(dn_tot, dn_norm_tot, ds_tot, ds_norm_tot) {
  return((dn_tot / dn_norm_tot) / (ds_tot / ds_norm_tot))
}

table_albicollis_macro <- subset(table_albicollis, category == "Macro")
table_albicollis_macro<-na.omit(table_albicollis_macro)

table_albicollis_inter <- subset(table_albicollis, category == "Intermediate")
table_albicollis_inter<-na.omit(table_albicollis_inter)

table_albicollis_micro <- subset(table_albicollis, category == "Micro")
table_albicollis_micro<-na.omit(table_albicollis_micro)
table_albicollis_micro <- table_albicollis_micro[table_albicollis_micro$chr != "Chr22", ]
table_albicollis_micro <- table_albicollis_micro[table_albicollis_micro$chr != "Chr25", ]
table_albicollis_micro <- table_albicollis_micro[table_albicollis_micro$chr != "Chr28", ]

table_taiga_macro <- subset(table_taiga, category == "Macro")
table_taiga_macro<-na.omit(table_taiga_macro)

table_taiga_inter <- subset(table_taiga, category == "Intermediate")
table_taiga_inter<-na.omit(table_taiga_inter)

table_taiga_micro <- subset(table_taiga, category == "Micro")
table_taiga_micro<-na.omit(table_taiga_micro)
table_taiga_micro <- table_taiga_micro[table_taiga_micro$chr != "Chr22", ]
table_taiga_micro <- table_taiga_micro[table_taiga_micro$chr != "Chr25", ]
table_taiga_micro <- table_taiga_micro[table_taiga_micro$chr != "Chr28", ]

# coll_dn_ds_gc_L
coll_dn_ds_gc_L_macro <- aggregate(cbind(dn_tot_gc_L, dn_norm_tot_gc_L, ds_tot_gc_L, ds_norm_tot_gc_L) ~ category == "Macro", table_albicollis_macro, mean)
coll_dn_ds_gc_L_macro$coll_dn_ds_gc_L <- with(coll_dn_ds_gc_L_macro, calculer_dn_ds(dn_tot_gc_L, dn_norm_tot_gc_L, ds_tot_gc_L, ds_norm_tot_gc_L))
coll_dn_ds_gc_L_macro<-mean(coll_dn_ds_gc_L_macro$coll_dn_ds_gc_L)

coll_dn_ds_gc_L_inter <- aggregate(cbind(dn_tot_gc_L, dn_norm_tot_gc_L, ds_tot_gc_L, ds_norm_tot_gc_L) ~ category == "Intermediate", table_albicollis_inter, mean)
coll_dn_ds_gc_L_inter$coll_dn_ds_gc_L <- with(coll_dn_ds_gc_L_inter, calculer_dn_ds(dn_tot_gc_L, dn_norm_tot_gc_L, ds_tot_gc_L, ds_norm_tot_gc_L))
coll_dn_ds_gc_L_inter<-mean(coll_dn_ds_gc_L_inter$coll_dn_ds_gc_L)

coll_dn_ds_gc_L_micro <- aggregate(cbind(dn_tot_gc_L, dn_norm_tot_gc_L, ds_tot_gc_L, ds_norm_tot_gc_L) ~ category == "Micro", table_albicollis_micro, mean)
coll_dn_ds_gc_L_micro$coll_dn_ds_gc_L <- with(coll_dn_ds_gc_L_micro, calculer_dn_ds(dn_tot_gc_L, dn_norm_tot_gc_L, ds_tot_gc_L, ds_norm_tot_gc_L))
coll_dn_ds_gc_L_micro<-mean(coll_dn_ds_gc_L_micro$coll_dn_ds_gc_L)

coll_dn_ds_gc_L <- data.frame(coll_dn_ds_gc_L = c(coll_dn_ds_gc_L_macro, coll_dn_ds_gc_L_inter, coll_dn_ds_gc_L_micro),row.names = c("Macro", "Intermediate", "Micro"))

#coll_dn_ds_gc_S

coll_dn_ds_gc_S_macro <- aggregate(cbind(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S) ~ category == "Macro", table_albicollis_macro, mean)
coll_dn_ds_gc_S_macro$coll_dn_ds_gc_S <- with(coll_dn_ds_gc_S_macro, calculer_dn_ds(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S))
coll_dn_ds_gc_S_macro<-mean(coll_dn_ds_gc_S_macro$coll_dn_ds_gc_S)

coll_dn_ds_gc_S_inter <- aggregate(cbind(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S) ~ category == "Intermediate" , table_albicollis_inter, mean)
coll_dn_ds_gc_S_inter$coll_dn_ds_gc_S <- with(coll_dn_ds_gc_S_inter, calculer_dn_ds(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S))
coll_dn_ds_gc_S_inter<-mean(coll_dn_ds_gc_S_inter$coll_dn_ds_gc_S)

coll_dn_ds_gc_S_micro <- aggregate(cbind(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S) ~ category == "Micro", table_albicollis_micro, mean)
coll_dn_ds_gc_S_micro$coll_dn_ds_gc_S <- with(coll_dn_ds_gc_S_micro, calculer_dn_ds(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S))
coll_dn_ds_gc_S_micro<-mean(coll_dn_ds_gc_S_micro$coll_dn_ds_gc_S)

coll_dn_ds_gc_S <- data.frame(coll_dn_ds_gc_S = c(coll_dn_ds_gc_S_macro, coll_dn_ds_gc_S_inter, coll_dn_ds_gc_S_micro),row.names = c("Macro", "Intermediate", "Micro"))

# coll_dn_ds_all_L
coll_dn_ds_all_L_macro <- aggregate(cbind(dn_tot_all_L, dn_norm_tot_all_L, ds_tot_all_L, ds_norm_tot_all_L) ~ category == "Macro", table_albicollis_macro, mean)
coll_dn_ds_all_L_macro$coll_dn_ds_all_L <- with(coll_dn_ds_all_L_macro, calculer_dn_ds(dn_tot_all_L, dn_norm_tot_all_L, ds_tot_all_L, ds_norm_tot_all_L))
coll_dn_ds_all_L_macro<-mean(coll_dn_ds_all_L_macro$coll_dn_ds_all_L)

coll_dn_ds_all_L_inter <- aggregate(cbind(dn_tot_all_L, dn_norm_tot_all_L, ds_tot_all_L, ds_norm_tot_all_L) ~ category == "Intermediate", table_albicollis_inter, mean)
coll_dn_ds_all_L_inter$coll_dn_ds_all_L <- with(coll_dn_ds_all_L_inter, calculer_dn_ds(dn_tot_all_L, dn_norm_tot_all_L, ds_tot_all_L, ds_norm_tot_all_L))
coll_dn_ds_all_L_inter<-mean(coll_dn_ds_all_L_inter$coll_dn_ds_all_L)

coll_dn_ds_all_L_micro <- aggregate(cbind(dn_tot_all_L, dn_norm_tot_all_L, ds_tot_all_L, ds_norm_tot_all_L) ~ category == "Micro", table_albicollis_micro, mean)
coll_dn_ds_all_L_micro$coll_dn_ds_all_L <- with(coll_dn_ds_all_L_micro, calculer_dn_ds(dn_tot_all_L, dn_norm_tot_all_L, ds_tot_all_L, ds_norm_tot_all_L))
coll_dn_ds_all_L_micro<-mean(coll_dn_ds_all_L_micro$coll_dn_ds_all_L)

coll_dn_ds_all_L <- data.frame(coll_dn_ds_all_L = c(coll_dn_ds_all_L_macro, coll_dn_ds_all_L_inter, coll_dn_ds_all_L_micro),row.names = c("Macro", "Intermediate", "Micro"))

#coll_dn_ds_all_S

coll_dn_ds_all_S_macro <- aggregate(cbind(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S) ~ category == "Macro", table_albicollis_macro, mean)
coll_dn_ds_all_S_macro$coll_dn_ds_all_S <- with(coll_dn_ds_all_S_macro, calculer_dn_ds(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S))
coll_dn_ds_all_S_macro<-mean(coll_dn_ds_all_S_macro$coll_dn_ds_all_S)

coll_dn_ds_all_S_inter <- aggregate(cbind(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S) ~ category == "Intermediate" , table_albicollis_inter, mean)
coll_dn_ds_all_S_inter$coll_dn_ds_all_S <- with(coll_dn_ds_all_S_inter, calculer_dn_ds(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S))
coll_dn_ds_all_S_inter<-mean(coll_dn_ds_all_S_inter$coll_dn_ds_all_S)

coll_dn_ds_all_S_micro <- aggregate(cbind(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S) ~ category == "Micro", table_albicollis_micro, mean)
coll_dn_ds_all_S_micro$coll_dn_ds_all_S <- with(coll_dn_ds_all_S_micro, calculer_dn_ds(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S))
coll_dn_ds_all_S_micro<-mean(coll_dn_ds_all_S_micro$coll_dn_ds_all_S)
      
coll_dn_ds_all_S <- data.frame(coll_dn_ds_all_S = c(coll_dn_ds_all_S_macro, coll_dn_ds_all_S_inter, coll_dn_ds_all_S_micro),row.names = c("Macro", "Intermediate", "Micro"))

#taiga_dn_ds_all_S

taiga_dn_ds_all_S_macro <- aggregate(cbind(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S) ~ category == "Macro", table_taiga_macro, mean)
taiga_dn_ds_all_S_macro$taiga_dn_ds_all_S <- with(taiga_dn_ds_all_S_macro, calculer_dn_ds(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S))
taiga_dn_ds_all_S_macro<-mean(taiga_dn_ds_all_S_macro$taiga_dn_ds_all_S)

taiga_dn_ds_all_S_inter <- aggregate(cbind(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S) ~ category == "Intermediate" , table_taiga_inter, mean)
taiga_dn_ds_all_S_inter$taiga_dn_ds_all_S <- with(taiga_dn_ds_all_S_inter, calculer_dn_ds(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S))
taiga_dn_ds_all_S_inter<-mean(taiga_dn_ds_all_S_inter$taiga_dn_ds_all_S)

taiga_dn_ds_all_S_micro <- aggregate(cbind(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S) ~ category == "Micro", table_taiga_micro, mean)
taiga_dn_ds_all_S_micro$taiga_dn_ds_all_S <- with(taiga_dn_ds_all_S_micro, calculer_dn_ds(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S))
taiga_dn_ds_all_S_micro<-mean(taiga_dn_ds_all_S_micro$taiga_dn_ds_all_S)

taiga_dn_ds_all_S <- data.frame(taiga_dn_ds_all_S = c(taiga_dn_ds_all_S_macro, taiga_dn_ds_all_S_inter, taiga_dn_ds_all_S_micro),row.names = c("Macro", "Intermediate", "Micro"))

#taiga_dn_ds_gc_S

taiga_dn_ds_gc_S_macro <- aggregate(cbind(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S) ~ category == "Macro", table_taiga_macro, mean)
taiga_dn_ds_gc_S_macro$taiga_dn_ds_gc_S <- with(taiga_dn_ds_gc_S_macro, calculer_dn_ds(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S))
taiga_dn_ds_gc_S_macro<-mean(taiga_dn_ds_gc_S_macro$taiga_dn_ds_gc_S)

taiga_dn_ds_gc_S_inter <- aggregate(cbind(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S) ~ category == "Intermediate" , table_taiga_inter, mean)
taiga_dn_ds_gc_S_inter$taiga_dn_ds_gc_S <- with(taiga_dn_ds_gc_S_inter, calculer_dn_ds(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S))
taiga_dn_ds_gc_S_inter<-mean(taiga_dn_ds_gc_S_inter$taiga_dn_ds_gc_S)

taiga_dn_ds_gc_S_micro <- aggregate(cbind(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S) ~ category == "Micro", table_taiga_micro, mean)
taiga_dn_ds_gc_S_micro$taiga_dn_ds_gc_S <- with(taiga_dn_ds_gc_S_micro, calculer_dn_ds(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S))
taiga_dn_ds_gc_S_micro<-mean(taiga_dn_ds_gc_S_micro$taiga_dn_ds_gc_S)

taiga_dn_ds_gc_S <- data.frame(taiga_dn_ds_gc_S = c(taiga_dn_ds_gc_S_macro, taiga_dn_ds_gc_S_inter, taiga_dn_ds_gc_S_micro),row.names = c("Macro", "Intermediate", "Micro"))

dndscat<-data.frame(c("Macro", "Intermediate", "Micro"))
names(dndscat)<-c("category")

dn_ds_cat<-cbind(dndscat, coll_dn_ds_gc_L, coll_dn_ds_gc_S, coll_dn_ds_all_L, coll_dn_ds_all_S, taiga_dn_ds_all_S, taiga_dn_ds_gc_S)
dn_ds_cat$category <- factor(dn_ds_cat$category, levels = c("Macro", "Intermediate", "Micro"))

palette_couleurs <- c("coll dn/ds GC L" = "steelblue4", 
                      "coll dn/ds GC S" = "steelblue1",
                      "coll dn/ds ALL L" = "tomato3",
                      "coll dn/ds ALL S" = "tomato",
                      "taiga dn/ds GC S" = "steelblue1",
                      "taiga dn/ds ALL S" = "tomato")

ggplot(dn_ds_cat) +
  geom_point(aes(x = category, y = coll_dn_ds_gc_L, color = "coll dn/ds GC L"), shape = 15, size = 3) +
  geom_point(aes(x = category, y = coll_dn_ds_gc_S, color = "coll dn/ds GC S"), shape = 15, size = 3) +
  geom_point(aes(x = category, y = coll_dn_ds_all_L, color = "coll dn/ds ALL L"), shape = 15, size = 3) +
  geom_point(aes(x = category, y = coll_dn_ds_all_S, color = "coll dn/ds ALL S"), shape = 15, size = 3) +
  geom_point(aes(x = category, y = taiga_dn_ds_gc_S, color = "taiga dn/ds GC S"), shape = 16, size = 3) +
  geom_point(aes(x = category, y = taiga_dn_ds_all_S, color = "taiga dn/ds ALL S"), shape = 16, size = 3) +
  scale_color_manual(values = palette_couleurs) +
  labs(x = "chromosome category", y = "dn/ds", color = "") +
  ggtitle("dn/ds for each category")

######################################################################################################""

# dn/ds groupe new category Macro Micro

table_albicollis <- read.table("table_albicollis.txt", header = T) 
table_taiga <- read.table("table_taiga.txt", header = T)

calculer_dn_ds <- function(dn_tot, dn_norm_tot, ds_tot, ds_norm_tot) {
  return((dn_tot / dn_norm_tot) / (ds_tot / ds_norm_tot))
}

table_albicollis_macro <- subset(table_albicollis, new_category == "Macro")
table_albicollis_macro<-na.omit(table_albicollis_macro)

table_albicollis_micro <- subset(table_albicollis, new_category == "Micro")
table_albicollis_micro<-na.omit(table_albicollis_micro)
table_albicollis_micro <- table_albicollis_micro[table_albicollis_micro$chr != "Chr22", ]
table_albicollis_micro <- table_albicollis_micro[table_albicollis_micro$chr != "Chr25", ]
table_albicollis_micro <- table_albicollis_micro[table_albicollis_micro$chr != "Chr28", ]


table_taiga_macro <- subset(table_taiga, new_category == "Macro")
table_taiga_macro<-na.omit(table_taiga_macro)

table_taiga_micro <- subset(table_taiga, new_category == "Micro")
table_taiga_micro<-na.omit(table_taiga_micro)
table_taiga_micro <- table_taiga_micro[table_taiga_micro$chr != "Chr22", ]
table_taiga_micro <- table_taiga_micro[table_taiga_micro$chr != "Chr25", ]
table_taiga_micro <- table_taiga_micro[table_taiga_micro$chr != "Chr28", ]

# coll_dn_ds_gc_L

coll_dn_ds_gc_L_macro <- aggregate(cbind(dn_tot_gc_L, dn_norm_tot_gc_L, ds_tot_gc_L, ds_norm_tot_gc_L) ~ new_category == "Macro", table_albicollis_macro, mean)
coll_dn_ds_gc_L_macro$coll_dn_ds_gc_L <- with(coll_dn_ds_gc_L_macro, calculer_dn_ds(dn_tot_gc_L, dn_norm_tot_gc_L, ds_tot_gc_L, ds_norm_tot_gc_L))
coll_dn_ds_gc_L_macro<-mean(coll_dn_ds_gc_L_macro$coll_dn_ds_gc_L)

coll_dn_ds_gc_L_micro <- aggregate(cbind(dn_tot_gc_L, dn_norm_tot_gc_L, ds_tot_gc_L, ds_norm_tot_gc_L) ~ new_category == "Micro", table_albicollis_micro, mean)
coll_dn_ds_gc_L_micro$coll_dn_ds_gc_L <- with(coll_dn_ds_gc_L_micro, calculer_dn_ds(dn_tot_gc_L, dn_norm_tot_gc_L, ds_tot_gc_L, ds_norm_tot_gc_L))
coll_dn_ds_gc_L_micro<-mean(coll_dn_ds_gc_L_micro$coll_dn_ds_gc_L)

coll_dn_ds_gc_L <- data.frame(coll_dn_ds_gc_L = c(coll_dn_ds_gc_L_macro, coll_dn_ds_gc_L_micro),row.names = c("Macro", "Micro"))

#coll_dn_ds_gc_S

coll_dn_ds_gc_S_macro <- aggregate(cbind(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S) ~ new_category == "Macro", table_albicollis_macro, mean)
coll_dn_ds_gc_S_macro$coll_dn_ds_gc_S <- with(coll_dn_ds_gc_S_macro, calculer_dn_ds(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S))
coll_dn_ds_gc_S_macro<-mean(coll_dn_ds_gc_S_macro$coll_dn_ds_gc_S)

coll_dn_ds_gc_S_micro <- aggregate(cbind(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S) ~ new_category == "Micro", table_albicollis_micro, mean)
coll_dn_ds_gc_S_micro$coll_dn_ds_gc_S <- with(coll_dn_ds_gc_S_micro, calculer_dn_ds(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S))
coll_dn_ds_gc_S_micro<-mean(coll_dn_ds_gc_S_micro$coll_dn_ds_gc_S)

coll_dn_ds_gc_S <- data.frame(coll_dn_ds_gc_S = c(coll_dn_ds_gc_S_macro, coll_dn_ds_gc_S_micro),row.names = c("Macro", "Micro"))

# coll_dn_ds_all_L

coll_dn_ds_all_L_macro <- aggregate(cbind(dn_tot_all_L, dn_norm_tot_all_L, ds_tot_all_L, ds_norm_tot_all_L) ~ new_category == "Macro", table_albicollis_macro, mean)
coll_dn_ds_all_L_macro$coll_dn_ds_all_L <- with(coll_dn_ds_all_L_macro, calculer_dn_ds(dn_tot_all_L, dn_norm_tot_all_L, ds_tot_all_L, ds_norm_tot_all_L))
coll_dn_ds_all_L_macro<-mean(coll_dn_ds_all_L_macro$coll_dn_ds_all_L)

coll_dn_ds_all_L_micro <- aggregate(cbind(dn_tot_all_L, dn_norm_tot_all_L, ds_tot_all_L, ds_norm_tot_all_L) ~ new_category == "Micro", table_albicollis_micro, mean)
coll_dn_ds_all_L_micro$coll_dn_ds_all_L <- with(coll_dn_ds_all_L_micro, calculer_dn_ds(dn_tot_all_L, dn_norm_tot_all_L, ds_tot_all_L, ds_norm_tot_all_L))
coll_dn_ds_all_L_micro<-mean(coll_dn_ds_all_L_micro$coll_dn_ds_all_L)

coll_dn_ds_all_L <- data.frame(coll_dn_ds_all_L = c(coll_dn_ds_all_L_macro, coll_dn_ds_all_L_micro),row.names = c("Macro", "Micro"))

#coll_dn_ds_all_S

coll_dn_ds_all_S_macro <- aggregate(cbind(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S) ~ new_category == "Macro", table_albicollis_macro, mean)
coll_dn_ds_all_S_macro$coll_dn_ds_all_S <- with(coll_dn_ds_all_S_macro, calculer_dn_ds(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S))
coll_dn_ds_all_S_macro<-mean(coll_dn_ds_all_S_macro$coll_dn_ds_all_S)

coll_dn_ds_all_S_micro <- aggregate(cbind(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S) ~ new_category == "Micro", table_albicollis_micro, mean)
coll_dn_ds_all_S_micro$coll_dn_ds_all_S <- with(coll_dn_ds_all_S_micro, calculer_dn_ds(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S))
coll_dn_ds_all_S_micro<-mean(coll_dn_ds_all_S_micro$coll_dn_ds_all_S)

coll_dn_ds_all_S <- data.frame(coll_dn_ds_all_S = c(coll_dn_ds_all_S_macro, coll_dn_ds_all_S_micro),row.names = c("Macro", "Micro"))

#taiga_dn_ds_all_S

taiga_dn_ds_all_S_macro <- aggregate(cbind(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S) ~ new_category == "Macro", table_taiga_macro, mean)
taiga_dn_ds_all_S_macro$taiga_dn_ds_all_S <- with(taiga_dn_ds_all_S_macro, calculer_dn_ds(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S))
taiga_dn_ds_all_S_macro<-mean(taiga_dn_ds_all_S_macro$taiga_dn_ds_all_S)

taiga_dn_ds_all_S_micro <- aggregate(cbind(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S) ~ new_category == "Micro", table_taiga_micro, mean)
taiga_dn_ds_all_S_micro$taiga_dn_ds_all_S <- with(taiga_dn_ds_all_S_micro, calculer_dn_ds(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S))
taiga_dn_ds_all_S_micro<-mean(taiga_dn_ds_all_S_micro$taiga_dn_ds_all_S)

taiga_dn_ds_all_S <- data.frame(taiga_dn_ds_all_S = c(taiga_dn_ds_all_S_macro, taiga_dn_ds_all_S_micro),row.names = c("Macro", "Micro"))

#taiga_dn_ds_gc_S

taiga_dn_ds_gc_S_macro <- aggregate(cbind(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S) ~ new_category == "Macro", table_taiga_macro, mean)
taiga_dn_ds_gc_S_macro$taiga_dn_ds_gc_S <- with(taiga_dn_ds_gc_S_macro, calculer_dn_ds(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S))
taiga_dn_ds_gc_S_macro<-mean(taiga_dn_ds_gc_S_macro$taiga_dn_ds_gc_S)

taiga_dn_ds_gc_S_micro <- aggregate(cbind(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S) ~ new_category == "Micro", table_taiga_micro, mean)
taiga_dn_ds_gc_S_micro$taiga_dn_ds_gc_S <- with(taiga_dn_ds_gc_S_micro, calculer_dn_ds(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S))
taiga_dn_ds_gc_S_micro<-mean(taiga_dn_ds_gc_S_micro$taiga_dn_ds_gc_S)

taiga_dn_ds_gc_S <- data.frame(taiga_dn_ds_gc_S = c(taiga_dn_ds_gc_S_macro, taiga_dn_ds_gc_S_micro),row.names = c("Macro", "Micro"))

dndscat<-data.frame(c("Macro", "Micro"))
names(dndscat)<-c("new_category")

dn_ds_cat<-cbind(dndscat, coll_dn_ds_gc_L, coll_dn_ds_gc_S, coll_dn_ds_all_L, coll_dn_ds_all_S, taiga_dn_ds_all_S, taiga_dn_ds_gc_S)
dn_ds_cat$category <- factor(dn_ds_cat$new_category, levels = c("Macro", "Micro"))
size<-read.table("chr_category.txt", header = T)
size<-size %>% rename(chr = chromosome)

dn_ds_cat$new_category <- factor(dn_ds_cat$new_category, levels = c("Macro", "Micro"))
library(ggplot2)
palette_couleurs <- c("coll dn/ds GC L" = "steelblue4", 
                      "coll dn/ds GC S" = "steelblue1",
                      "coll dn/ds ALL L" = "tomato3",
                      "coll dn/ds ALL S" = "tomato",
                      "taiga dn/ds GC S" = "steelblue1",
                      "taiga dn/ds ALL S" = "tomato")

ggplot(dn_ds_cat) +
  geom_point(aes(x = new_category, y = coll_dn_ds_gc_L, color = "coll dn/ds GC L"), shape = 15, size = 3) +
  geom_point(aes(x = new_category, y = coll_dn_ds_gc_S, color = "coll dn/ds GC S"), shape = 15, size = 3) +
  geom_point(aes(x = new_category, y = coll_dn_ds_all_L, color = "coll dn/ds ALL L"), shape = 15, size = 3) +
  geom_point(aes(x = new_category, y = coll_dn_ds_all_S, color = "coll dn/ds ALL S"), shape = 15, size = 3) +
  geom_point(aes(x = new_category, y = taiga_dn_ds_gc_S, color = "taiga dn/ds GC S"), shape = 16, size = 3) +
  geom_point(aes(x = new_category, y = taiga_dn_ds_all_S, color = "taiga dn/ds ALL S"), shape = 16, size = 3) +
  geom_line(aes(x = new_category, y = coll_dn_ds_gc_L, group = 1, color = "coll dn/ds GC L")) +
  geom_line(aes(x = new_category, y = coll_dn_ds_gc_S, group = 1, color = "coll dn/ds GC S")) +
  geom_line(aes(x = new_category, y = coll_dn_ds_all_L, group = 1, color = "coll dn/ds ALL L")) +
  geom_line(aes(x = new_category, y = coll_dn_ds_all_S, group = 1, color = "coll dn/ds ALL S")) +
  geom_line(aes(x = new_category, y = taiga_dn_ds_gc_S, group = 1, color = "taiga dn/ds GC S")) +
  geom_line(aes(x = new_category, y = taiga_dn_ds_all_S, group = 1, color = "taiga dn/ds ALL S")) +
  scale_color_manual(values = palette_couleurs) +
  labs(x = "chromosome category", y = "dn/ds", color = "") +
  ggtitle("")

##############################################

# chr ~ size

size<- read.table("chr_category.txt", header = T)
size<-subset(size, select= -category)
size$chromosome <- factor(size$chromosome, levels = mixedsort(size$chromosome))
deux_premieres_lignes <- size[1:2, ]
ligne_fusionnee <- data.frame(chromosome = "chr1", length = sum(deux_premieres_lignes$length))
size <- size[-c(1, 2), ]
size <- rbind(ligne_fusionnee, size)

deux_premieres_lignes1 <- size[c(4,6), ]
ligne_fusionnee <- data.frame(chromosome = "chr4", length = sum(deux_premieres_lignes1$length))
size <- size[-c(4, 6), ]
size <- rbind(ligne_fusionnee, size)

library(ggplot2)
ggplot(size, aes(x = chromosome, y = length)) +
  geom_point() +
  labs(x = "Chromosome", y = "Size in pb", title = "Size of chromosomes of the genus Ficedula")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot(merged_dn_ds$coll_dn_ds_gc_S, merged_dn_ds$taiga_dn_ds_gc_S)

##############################################

# Venn diagram 

gene_coll_dn_ds_gc_L <- aggregate(cbind(dn_tot_gc_L, dn_norm_tot_gc_L, ds_tot_gc_L, ds_norm_tot_gc_L) ~ gene, table_albicollis, mean)
gene_coll_dn_ds_gc_L$coll_dn_ds_gc_L <- with(gene_coll_dn_ds_gc_L, calculer_dn_ds(dn_tot_gc_L, dn_norm_tot_gc_L, ds_tot_gc_L, ds_norm_tot_gc_L))

gene_coll_dn_ds_gc_S <- aggregate(cbind(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S) ~ gene, table_albicollis, mean)
gene_coll_dn_ds_gc_S$coll_dn_ds_gc_S <- with(gene_coll_dn_ds_gc_S, calculer_dn_ds(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S))

gene_coll_dn_ds_all_L <- aggregate(cbind(dn_tot_all_L, dn_norm_tot_all_L, ds_tot_all_L, ds_norm_tot_all_L) ~ gene, table_albicollis, mean)
gene_coll_dn_ds_all_L$coll_dn_ds_all_L <- with(gene_coll_dn_ds_all_L, calculer_dn_ds(dn_tot_all_L, dn_norm_tot_all_L, ds_tot_all_L, ds_norm_tot_all_L))

gene_coll_dn_ds_all_S <- aggregate(cbind(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S) ~ gene, table_albicollis, mean)
gene_coll_dn_ds_all_S$coll_dn_ds_all_S <- with(gene_coll_dn_ds_all_S, calculer_dn_ds(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S))

gene_taiga_dn_ds_gc_S <- aggregate(cbind(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S) ~ gene, table_taiga, mean)
gene_taiga_dn_ds_gc_S$taiga_dn_ds_gc_S <- with(gene_taiga_dn_ds_gc_S, calculer_dn_ds(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S))

gene_taiga_dn_ds_all_S <- aggregate(cbind(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S) ~ gene, table_taiga, mean)
gene_taiga_dn_ds_all_S$taiga_dn_ds_all_S<- with(gene_taiga_dn_ds_all_S, calculer_dn_ds(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S))

## merge 
gene_dnds <- gene_coll_dn_ds_gc_L %>%
  left_join(gene_coll_dn_ds_gc_S, by = "gene") %>%
  left_join(gene_coll_dn_ds_all_L, by = "gene") %>%
  left_join(gene_coll_dn_ds_all_S, by = "gene") %>%
  left_join(gene_taiga_dn_ds_gc_S, by = "gene") %>%
  left_join(gene_taiga_dn_ds_all_S, by = "gene")
  
gene_coll_dn_ds_gc_S<- gene_dnds[, c('gene', 'coll_dn_ds_gc_S')]
gene_coll_dn_ds_gc_S<-na.omit(gene_coll_dn_ds_gc_S)
coll_dn_ds_gc_S <-gene_coll_dn_ds_gc_S$gene 

gene_coll_dn_ds_gc_L<- gene_dnds[, c('gene', 'coll_dn_ds_gc_L')]
gene_coll_dn_ds_gc_L<-na.omit(gene_coll_dn_ds_gc_L)
coll_dn_ds_gc_L <-gene_coll_dn_ds_gc_L$gene 

gene_coll_dn_ds_all_S<- gene_dnds[, c('gene', 'coll_dn_ds_all_S')]
gene_coll_dn_ds_all_S<-na.omit(gene_coll_dn_ds_all_S)
coll_dn_ds_all_S <-gene_coll_dn_ds_all_S$gene 

gene_coll_dn_ds_all_L<- gene_dnds[, c('gene', 'coll_dn_ds_all_L')]
gene_coll_dn_ds_all_L<-na.omit(gene_coll_dn_ds_all_L)
coll_dn_ds_all_L <-gene_coll_dn_ds_all_L$gene 

gene_taiga_dn_ds_gc_S<- gene_dnds[, c('gene', 'taiga_dn_ds_gc_S')]
gene_taiga_dn_ds_gc_S<-na.omit(gene_taiga_dn_ds_gc_S)
taiga_dn_ds_gc_S <-gene_taiga_dn_ds_gc_S$gene

gene_taiga_dn_ds_all_S<- gene_dnds[, c('gene', 'taiga_dn_ds_all_S')]
gene_taiga_dn_ds_all_S<-na.omit(gene_taiga_dn_ds_all_S)
taiga_dn_ds_all_S <-gene_taiga_dn_ds_all_S$gene

venn.plot <- venn.diagram(
  x = list(
    coll_dn_ds_gc_S = coll_dn_ds_gc_S,
    coll_dn_ds_all_S = coll_dn_ds_all_S,
    coll_dn_ds_all_L = coll_dn_ds_all_L,
    coll_dn_ds_gc_L = coll_dn_ds_gc_L,
    taiga_dn_ds_gc_S = taiga_dn_ds_gc_S,
    taiga_dn_ds_all_S = taiga_dn_ds_all_S
  ),
  category.names = c(
    "coll_dn_ds_gc_S", "coll_dn_ds_all_S",
    "coll_dn_ds_all_L", "coll_dn_ds_gc_L",
    "taiga_dn_ds_gc_S", "taiga_dn_ds_all_S"
  ),

  filename = NULL 
)

grid.draw(venn.plot)
 

x = list(
  coll_dn_ds_gc_S = coll_dn_ds_gc_S,
  coll_dn_ds_all_S = coll_dn_ds_all_S,
  coll_dn_ds_all_L = coll_dn_ds_all_L,
  coll_dn_ds_gc_L = coll_dn_ds_gc_L,
  taiga_dn_ds_gc_S = taiga_dn_ds_gc_S,
  taiga_dn_ds_all_S = taiga_dn_ds_all_S
)

library(ggvenn)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#CD531CFF", "#CD834CFF"),
  stroke_size = 0.5, set_name_size = 3
)

library("ggVennDiagram")
venn_plot <- ggVennDiagram(x, label_alpha = 0,set_color = "red")

##################################################################################"""

# common gene

table_albicollis <- read.table("table_albicollis.txt", header = T) 
table_albicollis <- select(table_albicollis, gene, chr,dn_tot_gc_L,dn_norm_tot_gc_L,ds_tot_gc_L,ds_norm_tot_gc_L,dn_tot_gc_S,dn_norm_tot_gc_S,ds_tot_gc_S,ds_norm_tot_gc_S,dn_tot_all_L, dn_norm_tot_all_L, ds_tot_all_L, ds_norm_tot_all_L,dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S)
table_taiga <- read.table("table_taiga.txt", header = T) 
table_taiga<- select(table_taiga,gene,chr, dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S, dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S)
table<-merge(table_albicollis, table_taiga, by="gene", all = T)
table<-na.omit(table)
table<-table %>% rename(dn_tot_gc_S_coll = dn_tot_gc_S.x)
table<-table %>% rename(dn_norm_tot_gc_S_coll = dn_norm_tot_gc_S.x)
table<-table %>% rename(ds_tot_gc_S_coll = ds_tot_gc_S.x)
table<-table %>% rename(ds_norm_tot_gc_S_coll = ds_norm_tot_gc_S.x)
table<-table %>% rename(dn_tot_all_S_coll = dn_tot_all_S.x)
table<-table %>% rename(dn_norm_tot_all_S_coll = dn_norm_tot_all_S.x)
table<-table %>% rename(ds_tot_all_S_coll = ds_tot_all_S.x)
table<-table %>% rename(ds_norm_tot_all_S_coll = ds_norm_tot_all_S.x)
table<-table %>% rename(dn_tot_gc_S_taiga = dn_tot_gc_S.y)
table<-table %>% rename(dn_norm_tot_gc_S_taiga = dn_norm_tot_gc_S.y)
table<-table %>% rename(ds_tot_gc_S_taiga = ds_tot_gc_S.y)
table<-table %>% rename(ds_norm_tot_gc_S_taiga = ds_norm_tot_gc_S.y)
table<-table %>% rename(dn_tot_all_S_taiga = dn_tot_all_S.y)
table<-table %>% rename(dn_norm_tot_all_S_taiga = dn_norm_tot_all_S.y)
table<-table %>% rename(ds_tot_all_S_taiga = ds_tot_all_S.y)
table<-table %>% rename(ds_norm_tot_all_S_taiga = ds_norm_tot_all_S.y)
table<-subset(table, select= -chr.y)
table<-table %>% rename(chr = chr.x)

## coll_dn_ds_gc_L

coll_dn_ds_gc_L <- aggregate(cbind(dn_tot_gc_L, dn_norm_tot_gc_L, ds_tot_gc_L, ds_norm_tot_gc_L) ~ chr, table, mean)
coll_dn_ds_gc_L$coll_dn_ds_gc_L <- with(coll_dn_ds_gc_L, calculer_dn_ds(dn_tot_gc_L, dn_norm_tot_gc_L, ds_tot_gc_L, ds_norm_tot_gc_L))

## coll_dn_ds_gc_S

coll_dn_ds_gc_S <- aggregate(cbind(dn_tot_gc_S_coll, dn_norm_tot_gc_S_coll, ds_tot_gc_S_coll, ds_norm_tot_gc_S_coll) ~ chr, table, mean)
coll_dn_ds_gc_S$coll_dn_ds_gc_S <- with(coll_dn_ds_gc_S, calculer_dn_ds(dn_tot_gc_S_coll, dn_norm_tot_gc_S_coll, ds_tot_gc_S_coll, ds_norm_tot_gc_S_coll))

## coll_dn_ds_all_L

coll_dn_ds_all_L <- aggregate(cbind(dn_tot_all_L, dn_norm_tot_all_L, ds_tot_all_L, ds_norm_tot_all_L) ~ chr, table, mean)
coll_dn_ds_all_L$coll_dn_ds_all_L <- with(coll_dn_ds_all_L, calculer_dn_ds(dn_tot_all_L, dn_norm_tot_all_L, ds_tot_all_L, ds_norm_tot_all_L))

## coll_dn_ds_all_S

coll_dn_ds_all_S <- aggregate(cbind(dn_tot_all_S_coll, dn_norm_tot_all_S_coll, ds_tot_all_S_coll, ds_norm_tot_all_S_coll) ~ chr, table, mean)
coll_dn_ds_all_S$coll_dn_ds_all_S <- with(coll_dn_ds_all_S, calculer_dn_ds(dn_tot_all_S_coll, dn_norm_tot_all_S_coll, ds_tot_all_S_coll, ds_norm_tot_all_S_coll))

## taiga_dn_ds_gc_S

taiga_dn_ds_gc_S <- aggregate(cbind(dn_tot_gc_S_taiga, dn_norm_tot_gc_S_taiga, ds_tot_gc_S_taiga, ds_norm_tot_gc_S_taiga) ~ chr, table, mean)
taiga_dn_ds_gc_S$taiga_dn_ds_gc_S <- with(taiga_dn_ds_gc_S, calculer_dn_ds(dn_tot_gc_S_taiga, dn_norm_tot_gc_S_taiga, ds_tot_gc_S_taiga, ds_norm_tot_gc_S_taiga))

## taiga_dn_ds_all.S

taiga_dn_ds_all_S <- aggregate(cbind(dn_tot_all_S_taiga, dn_norm_tot_all_S_taiga, ds_tot_all_S_taiga, ds_norm_tot_all_S_taiga) ~ chr, table, mean)
taiga_dn_ds_all_S$taiga_dn_ds_all_S<- with(taiga_dn_ds_all_S, calculer_dn_ds(dn_tot_all_S_taiga, dn_norm_tot_all_S_taiga, ds_tot_all_S_taiga, ds_norm_tot_all_S_taiga))

## merge 
merged_dn_ds <- coll_dn_ds_all_L %>%
  left_join(coll_dn_ds_all_S, by = "chr") %>%
  left_join(coll_dn_ds_gc_L, by = "chr") %>%
  left_join(coll_dn_ds_gc_S, by = "chr") %>%
  left_join(taiga_dn_ds_all_S, by = "chr") %>%
  left_join(taiga_dn_ds_gc_S, by = "chr")

merged_dn_ds <- merged_dn_ds %>%
  rename_with(~ gsub("\\.x$", "_coll", .x), ends_with(".x"))
table_albicollis <- table_albicollis %>%
  rename_with(~ gsub("\\.y$", "_taiga", .x), ends_with(".y"))

merged_dn_ds <- merged_dn_ds[merged_dn_ds$chr != "Chr22", ]
merged_dn_ds <- merged_dn_ds[merged_dn_ds$chr != "ChrLGE22", ]
merged_dn_ds <- merged_dn_ds[merged_dn_ds$chr != "Chr25", ]
merged_dn_ds <- merged_dn_ds[merged_dn_ds$chr != "Chr28", ]


matrix_chr <- merged_dn_ds[, c('chr', 'coll_dn_ds_gc_L', 'coll_dn_ds_all_L','coll_dn_ds_gc_S', 'coll_dn_ds_all_S','taiga_dn_ds_gc_S','taiga_dn_ds_all_S')]
matrix_chr <- matrix_chr %>%
  rename(dn_ds_gc_L_coll = coll_dn_ds_gc_L,
         dn_ds_all_L_coll = coll_dn_ds_all_L,
         dn_ds_gc_S_coll = coll_dn_ds_gc_S,
         dn_ds_all_S_coll = coll_dn_ds_all_S,
         dn_ds_gc_S_taiga = taiga_dn_ds_gc_S,
         dn_ds_all_S_taiga = taiga_dn_ds_all_S)

is.numeric(matrix_chr)
matrix_chr<- sapply(matrix_chr, as.numeric)

correlation_matrix <- cor(matrix_chr[, -1], method = "spearman")
corrplot(correlation_matrix, method="circle", type="upper", tl.col="black", tl.srt=45)
correlation_matrix

###################""

