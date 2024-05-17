setwd("/home/amaniouloux/Documents/data/dataframe//")
output_directory <- "/home/amaniouloux/Documents/data/Rstudio/résultats/"
library(dplyr)
library(stringr)
library(corrplot)
library(tidyr)
library(gtools)
library(ggplot2)

size<-read.table("chr_category.txt", header = T)
size<-size %>% rename(chr = chromosome)

table_albicollis <- read.table("table_albicollis.txt", header = T) 
table_hypoleuca <- read.table("table_hypoleuca.txt", header = T)
table_taiga <- read.table("table_taiga.txt", header = T)
#table_taiga <- table_taiga[table_taiga$pn_ps_all <= 1, ]
table_parva <-read.table("table_parva.txt", header = T)
table_albicollis$gene <- sapply(str_split(table_albicollis$gene, "_"), "[[", 1)
table_hypoleuca$gene <- sapply(str_split(table_hypoleuca$gene, "_"), "[[", 1)
table_parva$gene <- sapply(str_split(table_parva$gene, "_"), "[[", 1)
table_taiga$gene <- sapply(str_split(table_taiga$gene, "_"), "[[", 1)

pi_table_albicollis<- table_albicollis[, c('chr', 'category', 'new_category', 'het_zero_gc', 'len_zero_gc', 'het_four_gc', 'len_four_gc', 'het_zero_all','len_zero_all','het_four_all','len_four_all')]
pi_table_albicollis[, -c(1, 2,3)] <- lapply(pi_table_albicollis[, -c(1, 2,3)], as.numeric)
pi_table_albicollis<-na.omit(pi_table_albicollis)

pi_table_taiga<- table_taiga[, c('chr', 'category', 'new_category','het_zero_gc', 'len_zero_gc', 'het_four_gc', 'len_four_gc', 'het_zero_all','len_zero_all','het_four_all','len_four_all')]
pi_table_taiga[, -c(1, 2,3)] <- lapply(pi_table_taiga[, -c(1, 2,3)], as.numeric)
pi_table_taiga<-na.omit(pi_table_taiga)

pi_table_hypoleuca<- table_hypoleuca[, c('chr', 'category', 'new_category', 'het_zero_gc', 'len_zero_gc', 'het_four_gc', 'len_four_gc', 'het_zero_all','len_zero_all','het_four_all','len_four_all')]
pi_table_hypoleuca[, -c(1, 2,3)] <- lapply(pi_table_hypoleuca[, -c(1, 2,3)], as.numeric)
pi_table_hypoleuca<-na.omit(pi_table_hypoleuca)

pi_table_parva<- table_parva[, c('chr', 'category', 'new_category', 'het_zero_gc', 'len_zero_gc', 'het_four_gc', 'len_four_gc', 'het_zero_all','len_zero_all','het_four_all','len_four_all')]
pi_table_parva[, -c(1, 2,3)] <- lapply(pi_table_parva[, -c(1, 2,3)], as.numeric)
pi_table_parva<-na.omit(pi_table_parva)

 calculer_pn_ps <- function(het_zero, len_zero, het_four, len_four) {
  return((het_zero / len_zero) / (het_four / len_four))
}

 calculer_pn <- function(het_zero, len_zero) {
   return(het_zero / len_zero)
 }
 
 calculer_ps <- function(het_four, len_four) {
   return(het_four / len_four)
 }
 
## albicollis_pn_ps_gc
albicollis_pn_ps_gc <- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ chr, pi_table_albicollis, mean)
albicollis_pn_ps_gc$albicollis_pn_gc <- with(albicollis_pn_ps_gc, calculer_pn(het_zero_gc, len_zero_gc))
albicollis_pn_ps_gc$albicollis_ps_gc <- with(albicollis_pn_ps_gc, calculer_ps(het_four_gc, len_four_gc))
albicollis_pn_ps_gc$albicollis_pn_ps_gc <- with(albicollis_pn_ps_gc, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))

## albicollis_pn_ps_all
albicollis_pn_ps_all <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ chr, pi_table_albicollis, mean)
albicollis_pn_ps_all$albicollis_pn_all <- with(albicollis_pn_ps_all, calculer_pn(het_zero_all, len_zero_all))
albicollis_pn_ps_all$albicollis_ps_all <- with(albicollis_pn_ps_all, calculer_ps(het_four_all, len_four_all))
albicollis_pn_ps_all$albicollis_pn_ps_all <- with(albicollis_pn_ps_all, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))

## taiga_pn_ps_gc
taiga_pn_ps_gc <- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ chr, pi_table_taiga, mean)
taiga_pn_ps_gc$taiga_pn_gc <- with(taiga_pn_ps_gc, calculer_pn(het_zero_gc, len_zero_gc))
taiga_pn_ps_gc$taiga_ps_gc <- with(taiga_pn_ps_gc, calculer_ps(het_four_gc, len_four_gc))
taiga_pn_ps_gc$taiga_pn_ps_gc <- with(taiga_pn_ps_gc, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))

## taiga_pn_ps_all
taiga_pn_ps_all <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ chr, pi_table_taiga, mean)
taiga_pn_ps_all$taiga_pn_all <- with(taiga_pn_ps_all, calculer_pn(het_zero_all, len_zero_all))
taiga_pn_ps_all$taiga_ps_all <- with(taiga_pn_ps_all, calculer_ps(het_four_all, len_four_all))
taiga_pn_ps_all$taiga_pn_ps_all <- with(taiga_pn_ps_all, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))

## hypoleuca_pn_ps_gc
hypoleuca_pn_ps_gc <- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ chr, pi_table_hypoleuca, mean)
hypoleuca_pn_ps_gc$hypoleuca_pn_gc <- with(hypoleuca_pn_ps_gc, calculer_pn(het_zero_gc, len_zero_gc))
hypoleuca_pn_ps_gc$hypoleuca_ps_gc <- with(hypoleuca_pn_ps_gc, calculer_ps(het_four_gc, len_four_gc))
hypoleuca_pn_ps_gc$hypoleuca_pn_ps_gc <- with(hypoleuca_pn_ps_gc, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))

## hypoleuca_pn_ps_all
hypoleuca_pn_ps_all <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ chr, pi_table_hypoleuca, mean)
hypoleuca_pn_ps_all$hypoleuca_pn_all <- with(hypoleuca_pn_ps_all, calculer_pn(het_zero_all, len_zero_all))
hypoleuca_pn_ps_all$hypoleuca_ps_all <- with(hypoleuca_pn_ps_all, calculer_ps(het_four_all, len_four_all))
hypoleuca_pn_ps_all$hypoleuca_pn_ps_all <- with(hypoleuca_pn_ps_all, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))

## parva_pn_ps_gc
parva_pn_ps_gc <- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ chr, pi_table_parva, mean)
parva_pn_ps_gc$parva_pn_gc <- with(parva_pn_ps_gc, calculer_pn(het_zero_gc, len_zero_gc))
parva_pn_ps_gc$parva_ps_gc <- with(parva_pn_ps_gc, calculer_ps(het_four_gc, len_four_gc))
parva_pn_ps_gc$parva_pn_ps_gc <- with(parva_pn_ps_gc, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))

## parva_pn_ps_all
parva_pn_ps_all <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ chr, pi_table_parva, mean)
parva_pn_ps_all$parva_pn_all <- with(parva_pn_ps_all, calculer_pn(het_zero_all, len_zero_all))
parva_pn_ps_all$parva_ps_all <- with(parva_pn_ps_all, calculer_ps(het_four_all, len_four_all))
parva_pn_ps_all$parva_pn_ps_all <- with(parva_pn_ps_all, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))

##merge
pi_table <- albicollis_pn_ps_gc %>%
  left_join(hypoleuca_pn_ps_gc, by = "chr") %>%
  left_join(taiga_pn_ps_gc, by = "chr") %>%
  left_join(parva_pn_ps_gc, by = "chr") %>%
  left_join(albicollis_pn_ps_all, by = "chr") %>%
  left_join(hypoleuca_pn_ps_all, by = "chr") %>%
  left_join(taiga_pn_ps_all, by = "chr") %>%
  left_join(parva_pn_ps_all, by = "chr")

##############################################################################
# correlation matrix chr 

pi_table<- pi_table[, c('chr','albicollis_pn_ps_gc', 'hypoleuca_pn_ps_gc', 'taiga_pn_ps_gc', 'parva_pn_ps_gc', 'albicollis_pn_ps_all', 'hypoleuca_pn_ps_all','taiga_pn_ps_all','parva_pn_ps_all')]

pi_table <- pi_table %>% filter(chr != "Chr25")
pi_table <- pi_table %>% filter(chr != "ChrLGE22")

pi_table<-subset(pi_table, select= -chr)
pi_table<- sapply(pi_table, as.numeric)

is.data.frame(pi_table) 
is.numeric(pi_table)                           
correlation_matrix<-cor(pi_table, method = "spearman")
correlation_matrix

corrplot(correlation_matrix, method="circle", type="upper", tl.col="black", tl.srt=45)

col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = correlation_matrix, col = col, symm = TRUE)


#############################################################################

#distribution pN/ps

pi_table<- pi_table[, c('chr','albicollis_pn_ps_gc', 'hypoleuca_pn_ps_gc', 'taiga_pn_ps_gc', 'parva_pn_ps_gc', 'albicollis_pn_ps_all', 'hypoleuca_pn_ps_all','taiga_pn_ps_all','parva_pn_ps_all')]

pi_table$chr <- factor(pi_table$chr, levels = mixedsort(pi_table$chr))
pi_table<-merge(pi_table, size, by ="chr", all =T)
pi_table <- pi_table %>% filter(chr != "Chr24")
pi_table <- pi_table %>% filter(chr != "Chr25")
pi_table <- pi_table %>% filter(chr != "Chr22")
pi_table <- pi_table %>% filter(chr != "Chr26")
pi_table <- pi_table %>% filter(chr != "Chr23")
pi_table <- pi_table %>% filter(chr != "ChrLGE22")
pi_table <- pi_table %>% filter(chr != "ChrFal34")
pi_table <- pi_table %>% filter(chr != "ChrFal36")
pi_table <- pi_table %>% filter(chr != "Chr27")
pi_table <- pi_table %>% filter(chr != "Chr13")
pi_table <- pi_table %>% filter(chr != "Chr28")
pi_table <- pi_table %>% filter(chr != "ChrZ")


palette_couleurs <- c("pn/ps GC coll" = "steelblue1", 
                      "pn/ps GC pied" = "darkorchid",
                      "pn/ps GC parva" = "tomato",
                      "pn/ps GC taiga" = "gold",
                      "pn/ps ALL coll" = "steelblue1",
                      "pn/ps ALL pied" = "darkorchid",
                      "pn/ps ALL parva" = "tomato",
                      "pn/ps ALL taiga" = "gold")

ggplot() +
  geom_point(data = pi_table, aes(x = factor(chr, levels = unique(pi_table$chr[order(pi_table$length)])), 
                                  y = albicollis_pn_ps_gc, color = "pn/ps GC coll"), 
             shape = 16, size = 2) +
  geom_point(data = pi_table, aes(x = factor(chr, levels = unique(pi_table$chr[order(pi_table$length)])), 
                                  y = hypoleuca_pn_ps_gc, color = "pn/ps GC pied"), 
             shape = 16, size = 2) +
  geom_point(data = pi_table, aes(x = factor(chr, levels = unique(pi_table$chr[order(pi_table$length)])), 
                                  y = taiga_pn_ps_gc, color = "pn/ps GC taiga"), 
             shape = 16, size = 2) +
  geom_point(data = pi_table, aes(x = factor(chr, levels = unique(pi_table$chr[order(pi_table$length)])), 
                                  y = parva_pn_ps_gc, color = "pn/ps GC parva"), 
             shape = 16, size = 2) +
  geom_point(data = pi_table, aes(x = factor(chr, levels = unique(pi_table$chr[order(pi_table$length)])), 
                                  y = albicollis_pn_ps_all, color = "pn/ps ALL coll"), 
             shape = 17, size = 2) +
  geom_point(data = pi_table, aes(x = factor(chr, levels = unique(pi_table$chr[order(pi_table$length)])), 
                                  y = hypoleuca_pn_ps_all, color = "pn/ps ALL pied"), 
             shape = 17, size = 2) +
  geom_point(data = pi_table, aes(x = factor(chr, levels = unique(pi_table$chr[order(pi_table$length)])), 
                                  y = taiga_pn_ps_all, color = "pn/ps ALL taiga"), 
             shape = 17, size = 2) +
  geom_point(data = pi_table, aes(x = factor(chr, levels = unique(pi_table$chr[order(pi_table$length)])), 
                                  y = parva_pn_ps_all, color = "pn/ps ALL parva"), 
             shape = 17, size = 2) +
  scale_color_manual(values = palette_couleurs) +
  labs(x = "chromosome", y = "pn/ps", color = "Conditions") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

###################################################################################################

# pn/ps category macro inter micro

# coll_pn_ps_gc

pi_table_albicollis_macro <- subset(pi_table_albicollis, category == "Macro")
pi_table_albicollis_macro<-na.omit(pi_table_albicollis_macro)
coll_pn_ps_gc_macro<- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ category == "Macro", pi_table_albicollis_macro, mean)
coll_pn_ps_gc_macro$coll_pn_ps_gc<- with(coll_pn_ps_gc_macro, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))
coll_pn_ps_gc_macro<-mean(coll_pn_ps_gc_macro$coll_pn_ps_gc)

pi_table_albicollis_inter<- subset(pi_table_albicollis, category == "Intermediate")
pi_table_albicollis_inter<-na.omit(pi_table_albicollis_intermediate)
coll_pn_ps_gc_inter<- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ category == "Intermediate", pi_table_albicollis_inter, mean)
coll_pn_ps_gc_inter$coll_pn_ps_gc<- with(coll_pn_ps_gc_inter, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))
coll_pn_ps_gc_inter<-mean(coll_pn_ps_gc_inter$coll_pn_ps_gc)

pi_table_albicollis_micro <- subset(pi_table_albicollis, category == "Micro")
pi_table_albicollis_micro<-na.omit(pi_table_albicollis_micro)
coll_pn_ps_gc_micro<- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ category == "Micro", pi_table_albicollis_micro, mean)
coll_pn_ps_gc_micro$coll_pn_ps_gc<- with(coll_pn_ps_gc_micro, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))
coll_pn_ps_gc_micro<-mean(coll_pn_ps_gc_micro$coll_pn_ps_gc)

coll_pn_ps_gc<- data.frame(coll_pn_ps_gc= c(coll_pn_ps_gc_macro, coll_pn_ps_gc_inter, coll_pn_ps_gc_micro),row.names = c("Macro", "Intermediate", "Micro"))

# coll_pn_ps_all

pi_table_albicollis_macro <- subset(pi_table_albicollis, category == "Macro")
pi_table_albicollis_macro <- na.omit(pi_table_albicollis_macro)
coll_pn_ps_all_macro <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ category == "Macro", pi_table_albicollis_macro, mean)
coll_pn_ps_all_macro$coll_pn_ps_all <- with(coll_pn_ps_all_macro, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))
coll_pn_ps_all_macro <- mean(coll_pn_ps_all_macro$coll_pn_ps_all)

pi_table_albicollis_intermediate <- subset(pi_table_albicollis, category == "Intermediate")
pi_table_albicollis_intermediate <- na.omit(pi_table_albicollis_intermediate)
coll_pn_ps_all_intermediate <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ category == "Intermediate", pi_table_albicollis_intermediate, mean)
coll_pn_ps_all_intermediate$coll_pn_ps_all <- with(coll_pn_ps_all_intermediate, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))
coll_pn_ps_all_intermediate <- mean(coll_pn_ps_all_intermediate$coll_pn_ps_all)

pi_table_albicollis_micro <- subset(pi_table_albicollis, category == "Micro")
pi_table_albicollis_micro <- na.omit(pi_table_albicollis_micro)
coll_pn_ps_all_micro <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ category == "Micro", pi_table_albicollis_micro, mean)
coll_pn_ps_all_micro$coll_pn_ps_all <- with(coll_pn_ps_all_micro, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))
coll_pn_ps_all_micro <- mean(coll_pn_ps_all_micro$coll_pn_ps_all)

coll_pn_ps_all <- data.frame(coll_pn_ps_all = c(coll_pn_ps_all_macro, coll_pn_ps_all_intermediate, coll_pn_ps_all_micro), row.names = c("Macro", "Intermediate", "Micro"))

# pied_dn_ds_gc

pi_table_hypoleuca_macro <- subset(pi_table_hypoleuca, category == "Macro")
pi_table_hypoleuca_macro<-na.omit(pi_table_hypoleuca_macro)
pied_pn_ps_gc_macro<- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ category == "Macro", pi_table_hypoleuca_macro, mean)
pied_pn_ps_gc_macro$pied_pn_ps_gc<- with(pied_pn_ps_gc_macro, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))
pied_pn_ps_gc_macro<-mean(pied_pn_ps_gc_macro$pied_pn_ps_gc)

pi_table_hypoleuca_inter<- subset(pi_table_hypoleuca, category == "Intermediate")
pi_table_hypoleuca_inter<-na.omit(pi_table_hypoleuca_inter)
pied_pn_ps_gc_inter<- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ category == "Intermediate", pi_table_hypoleuca_inter, mean)
pied_pn_ps_gc_inter$pied_pn_ps_gc<- with(pied_pn_ps_gc_inter, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))
pied_pn_ps_gc_inter<-mean(pied_pn_ps_gc_inter$pied_pn_ps_gc)

pi_table_hypoleuca_micro <- subset(pi_table_hypoleuca, category == "Micro")
pi_table_hypoleuca_micro<-na.omit(pi_table_hypoleuca_micro)
pied_pn_ps_gc_micro<- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ category == "Micro", pi_table_hypoleuca_micro, mean)
pied_pn_ps_gc_micro$pied_pn_ps_gc<- with(pied_pn_ps_gc_micro, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))
pied_pn_ps_gc_micro<-mean(pied_pn_ps_gc_micro$pied_pn_ps_gc)

pied_pn_ps_gc<- data.frame(pied_pn_ps_gc= c(pied_pn_ps_gc_macro, pied_pn_ps_gc_inter, pied_pn_ps_gc_micro),row.names = c("Macro", "Intermediate", "Micro"))

# pied_pn_ps_all

pi_table_hypoleuca_macro <- subset(pi_table_hypoleuca, category == "Macro")
pi_table_hypoleuca_macro <- na.omit(pi_table_hypoleuca_macro)
pied_pn_ps_all_macro <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ category == "Macro", pi_table_hypoleuca_macro, mean)
pied_pn_ps_all_macro$pied_pn_ps_all <- with(pied_pn_ps_all_macro, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))
pied_pn_ps_all_macro <- mean(pied_pn_ps_all_macro$pied_pn_ps_all)

pi_table_hypoleuca_intermediate <- subset(pi_table_hypoleuca, category == "Intermediate")
pi_table_hypoleuca_intermediate <- na.omit(pi_table_hypoleuca_intermediate)
pied_pn_ps_all_intermediate <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ category == "Intermediate", pi_table_hypoleuca_intermediate, mean)
pied_pn_ps_all_intermediate$pied_pn_ps_all <- with(pied_pn_ps_all_intermediate, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))
pied_pn_ps_all_intermediate <- mean(pied_pn_ps_all_intermediate$pied_pn_ps_all)

pi_table_hypoleuca_micro <- subset(pi_table_hypoleuca, category == "Micro")
pi_table_hypoleuca_micro <- na.omit(pi_table_hypoleuca_micro)
pied_pn_ps_all_micro <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ category == "Micro", pi_table_hypoleuca_micro, mean)
pied_pn_ps_all_micro$pied_pn_ps_all <- with(pied_pn_ps_all_micro, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))
pied_pn_ps_all_micro <- mean(pied_pn_ps_all_micro$pied_pn_ps_all)

pied_pn_ps_all <- data.frame(pied_pn_ps_all = c(pied_pn_ps_all_macro, pied_pn_ps_all_intermediate, pied_pn_ps_all_micro), row.names = c("Macro", "Intermediate", "Micro"))

# taiga_pn_ps_gc

pi_table_taiga_macro <- subset(pi_table_taiga, category == "Macro")
pi_table_taiga_macro <- na.omit(pi_table_taiga_macro)
taiga_pn_ps_gc_macro <- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ category == "Macro", pi_table_taiga_macro, mean)
taiga_pn_ps_gc_macro$taiga_pn_ps_gc <- with(taiga_pn_ps_gc_macro, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))
taiga_pn_ps_gc_macro <- mean(taiga_pn_ps_gc_macro$taiga_pn_ps_gc)

pi_table_taiga_intermediate <- subset(pi_table_taiga, category == "Intermediate")
pi_table_taiga_intermediate <- na.omit(pi_table_taiga_intermediate)
taiga_pn_ps_gc_intermediate <- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ category == "Intermediate", pi_table_taiga_intermediate, mean)
taiga_pn_ps_gc_intermediate$taiga_pn_ps_gc <- with(taiga_pn_ps_gc_intermediate, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))
taiga_pn_ps_gc_intermediate <- mean(taiga_pn_ps_gc_intermediate$taiga_pn_ps_gc)

pi_table_taiga_micro <- subset(pi_table_taiga, category == "Micro")
pi_table_taiga_micro <- na.omit(pi_table_taiga_micro)
taiga_pn_ps_gc_micro <- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ category == "Micro", pi_table_taiga_micro, mean)
taiga_pn_ps_gc_micro$taiga_pn_ps_gc <- with(taiga_pn_ps_gc_micro, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))
taiga_pn_ps_gc_micro <- mean(taiga_pn_ps_gc_micro$taiga_pn_ps_gc)

taiga_pn_ps_gc <- data.frame(taiga_pn_ps_gc = c(taiga_pn_ps_gc_macro, taiga_pn_ps_gc_intermediate, taiga_pn_ps_gc_micro), row.names = c("Macro", "Intermediate", "Micro"))

# taiga_pn_ps_all

pi_table_taiga_macro <- subset(pi_table_taiga, category == "Macro")
pi_table_taiga_macro <- na.omit(pi_table_taiga_macro)
taiga_pn_ps_all_macro <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ category == "Macro", pi_table_taiga_macro, mean)
taiga_pn_ps_all_macro$taiga_pn_ps_all <- with(taiga_pn_ps_all_macro, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))
taiga_pn_ps_all_macro <- mean(taiga_pn_ps_all_macro$taiga_pn_ps_all)

pi_table_taiga_intermediate <- subset(pi_table_taiga, category == "Intermediate")
pi_table_taiga_intermediate <- na.omit(pi_table_taiga_intermediate)
taiga_pn_ps_all_intermediate <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ category == "Intermediate", pi_table_taiga_intermediate, mean)
taiga_pn_ps_all_intermediate$taiga_pn_ps_all <- with(taiga_pn_ps_all_intermediate, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))
taiga_pn_ps_all_intermediate <- mean(taiga_pn_ps_all_intermediate$taiga_pn_ps_all)

pi_table_taiga_micro <- subset(pi_table_taiga, category == "Micro")
pi_table_taiga_micro <- na.omit(pi_table_taiga_micro)
taiga_pn_ps_all_micro <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ category == "Micro", pi_table_taiga_micro, mean)
taiga_pn_ps_all_micro$taiga_pn_ps_all <- with(taiga_pn_ps_all_micro, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))
taiga_pn_ps_all_micro <- mean(taiga_pn_ps_all_micro$taiga_pn_ps_all)

taiga_pn_ps_all <- data.frame(taiga_pn_ps_all = c(taiga_pn_ps_all_macro, taiga_pn_ps_all_intermediate, taiga_pn_ps_all_micro), row.names = c("Macro", "Intermediate", "Micro"))

# parva_pn_ps_gc

pi_table_parva_macro <- subset(pi_table_parva, category == "Macro")
pi_table_parva_macro <- na.omit(pi_table_parva_macro)
parva_pn_ps_gc_macro <- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ category == "Macro", pi_table_parva_macro, mean)
parva_pn_ps_gc_macro$parva_pn_ps_gc <- with(parva_pn_ps_gc_macro, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))
parva_pn_ps_gc_macro <- mean(parva_pn_ps_gc_macro$parva_pn_ps_gc)

pi_table_parva_intermediate <- subset(pi_table_parva, category == "Intermediate")
pi_table_parva_intermediate <- na.omit(pi_table_parva_intermediate)
parva_pn_ps_gc_intermediate <- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ category == "Intermediate", pi_table_parva_intermediate, mean)
parva_pn_ps_gc_intermediate$parva_pn_ps_gc <- with(parva_pn_ps_gc_intermediate, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))
parva_pn_ps_gc_intermediate <- mean(parva_pn_ps_gc_intermediate$parva_pn_ps_gc)

pi_table_parva_micro <- subset(pi_table_parva, category == "Micro")
pi_table_parva_micro <- na.omit(pi_table_parva_micro)
parva_pn_ps_gc_micro <- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ category == "Micro", pi_table_parva_micro, mean)
parva_pn_ps_gc_micro$parva_pn_ps_gc <- with(parva_pn_ps_gc_micro, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))
parva_pn_ps_gc_micro <- mean(parva_pn_ps_gc_micro$parva_pn_ps_gc)

parva_pn_ps_gc <- data.frame(parva_pn_ps_gc = c(parva_pn_ps_gc_macro, parva_pn_ps_gc_intermediate, parva_pn_ps_gc_micro), row.names = c("Macro", "Intermediate", "Micro"))

# parva_pn_ps_all

pi_table_parva_macro <- subset(pi_table_parva, category == "Macro")
pi_table_parva_macro <- na.omit(pi_table_parva_macro)
parva_pn_ps_all_macro <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ category == "Macro", pi_table_parva_macro, mean)
parva_pn_ps_all_macro$parva_pn_ps_all <- with(parva_pn_ps_all_macro, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))
parva_pn_ps_all_macro <- mean(parva_pn_ps_all_macro$parva_pn_ps_all)

pi_table_parva_intermediate <- subset(pi_table_parva, category == "Intermediate")
pi_table_parva_intermediate <- na.omit(pi_table_parva_intermediate)
parva_pn_ps_all_intermediate <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ category == "Intermediate", pi_table_parva_intermediate, mean)
parva_pn_ps_all_intermediate$parva_pn_ps_all <- with(parva_pn_ps_all_intermediate, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))
parva_pn_ps_all_intermediate <- mean(parva_pn_ps_all_intermediate$parva_pn_ps_all)

pi_table_parva_micro <- subset(pi_table_parva, category == "Micro")
pi_table_parva_micro <- na.omit(pi_table_parva_micro)
parva_pn_ps_all_micro <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ category == "Micro", pi_table_parva_micro, mean)
parva_pn_ps_all_micro$parva_pn_ps_all <- with(parva_pn_ps_all_micro, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))
parva_pn_ps_all_micro <- mean(parva_pn_ps_all_micro$parva_pn_ps_all)

parva_pn_ps_all <- data.frame(parva_pn_ps_all = c(parva_pn_ps_all_macro, parva_pn_ps_all_intermediate, parva_pn_ps_all_micro), row.names = c("Macro", "Intermediate", "Micro"))


pnpscat<-data.frame(c("Macro", "Intermediate", "Micro"))
names(pnpscat)<-c("category")

pn_ps_cat<-cbind(pnpscat, coll_pn_ps_gc, coll_pn_ps_all, pied_pn_ps_gc, pied_pn_ps_all, taiga_pn_ps_gc, taiga_pn_ps_all, parva_pn_ps_gc, parva_pn_ps_all)
pn_ps_cat$category <- factor(pn_ps_cat$category, levels = c("Macro", "Intermediate", "Micro"))

palette_couleurs <- c("pn/ps GC coll" = "steelblue1", 
                      "pn/ps GC pied" = "darkorchid",
                      "pn/ps GC parva" = "tomato",
                      "pn/ps GC taiga" = "gold",
                      "pn/ps ALL coll" = "steelblue1",
                      "pn/ps ALL pied" = "darkorchid",
                      "pn/ps ALL parva" = "tomato",
                      "pn/ps ALL taiga" = "gold")
ggplot(pn_ps_cat) +
  geom_point(aes(x = category, y = coll_pn_ps_gc, color = "pn/ps GC coll"), shape = 15, size = 3) +
  geom_point(aes(x = category, y = pied_pn_ps_gc, color = "pn/ps GC pied"), shape = 15, size = 3) +
  geom_point(aes(x = category, y = parva_pn_ps_gc, color = "pn/ps GC parva"), shape = 15, size = 3) +
  geom_point(aes(x = category, y = taiga_pn_ps_gc, color = "pn/ps GC taiga"), shape = 15, size = 3) +
  geom_point(aes(x = category, y = coll_pn_ps_all, color = "pn/ps ALL coll"), shape = 16, size = 3) +
  geom_point(aes(x = category, y = pied_pn_ps_all, color = "pn/ps ALL pied"), shape = 16, size = 3) +
  geom_point(aes(x = category, y = parva_pn_ps_all, color = "pn/ps ALL parva"), shape = 16, size = 3) +
  geom_point(aes(x = category, y = taiga_pn_ps_all, color = "pn/ps ALL taiga"), shape = 16, size = 3) +
  scale_color_manual(values = palette_couleurs) +
  labs(x = "chromosome category", y = "pn/ps", color = "") +
  ggtitle("pn/ps (GC) for each category")

##############################################################################################################""

# pn/ps new category macro  micro

# coll_pn_ps_gc

pi_table_albicollis_macro <- subset(pi_table_albicollis, new_category == "Macro")
pi_table_albicollis_macro<-na.omit(pi_table_albicollis_macro)
coll_pn_ps_gc_macro<- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ new_category == "Macro", pi_table_albicollis_macro, mean)
coll_pn_ps_gc_macro$coll_pn_ps_gc<- with(coll_pn_ps_gc_macro, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))
coll_pn_ps_gc_macro<-mean(coll_pn_ps_gc_macro$coll_pn_ps_gc)

pi_table_albicollis_micro <- subset(pi_table_albicollis, new_category == "Micro")
pi_table_albicollis_micro<-na.omit(pi_table_albicollis_micro)
coll_pn_ps_gc_micro<- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ new_category == "Micro", pi_table_albicollis_micro, mean)
coll_pn_ps_gc_micro$coll_pn_ps_gc<- with(coll_pn_ps_gc_micro, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))
coll_pn_ps_gc_micro<-mean(coll_pn_ps_gc_micro$coll_pn_ps_gc)

coll_pn_ps_gc<- data.frame(coll_pn_ps_gc= c(coll_pn_ps_gc_macro, coll_pn_ps_gc_micro),row.names = c("Macro", "Micro"))

# coll_pn_ps_all

pi_table_albicollis_macro <- subset(pi_table_albicollis, new_category == "Macro")
pi_table_albicollis_macro <- na.omit(pi_table_albicollis_macro)
coll_pn_ps_all_macro <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ new_category == "Macro", pi_table_albicollis_macro, mean)
coll_pn_ps_all_macro$coll_pn_ps_all <- with(coll_pn_ps_all_macro, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))
coll_pn_ps_all_macro <- mean(coll_pn_ps_all_macro$coll_pn_ps_all)


pi_table_albicollis_micro <- subset(pi_table_albicollis, new_category == "Micro")
pi_table_albicollis_micro <- na.omit(pi_table_albicollis_micro)
coll_pn_ps_all_micro <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ new_category == "Micro", pi_table_albicollis_micro, mean)
coll_pn_ps_all_micro$coll_pn_ps_all <- with(coll_pn_ps_all_micro, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))
coll_pn_ps_all_micro <- mean(coll_pn_ps_all_micro$coll_pn_ps_all)

coll_pn_ps_all <- data.frame(coll_pn_ps_all = c(coll_pn_ps_all_macro, coll_pn_ps_all_micro), row.names = c("Macro", "Micro"))

# pied_dn_ds_gc

pi_table_hypoleuca_macro <- subset(pi_table_hypoleuca, new_category == "Macro")
pi_table_hypoleuca_macro<-na.omit(pi_table_hypoleuca_macro)
pied_pn_ps_gc_macro<- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ new_category == "Macro", pi_table_hypoleuca_macro, mean)
pied_pn_ps_gc_macro$pied_pn_ps_gc<- with(pied_pn_ps_gc_macro, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))
pied_pn_ps_gc_macro<-mean(pied_pn_ps_gc_macro$pied_pn_ps_gc)


pi_table_hypoleuca_micro <- subset(pi_table_hypoleuca, new_category == "Micro")
pi_table_hypoleuca_micro<-na.omit(pi_table_hypoleuca_micro)
pied_pn_ps_gc_micro<- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ new_category == "Micro", pi_table_hypoleuca_micro, mean)
pied_pn_ps_gc_micro$pied_pn_ps_gc<- with(pied_pn_ps_gc_micro, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))
pied_pn_ps_gc_micro<-mean(pied_pn_ps_gc_micro$pied_pn_ps_gc)

pied_pn_ps_gc<- data.frame(pied_pn_ps_gc= c(pied_pn_ps_gc_macro, pied_pn_ps_gc_micro),row.names = c("Macro", "Micro"))

# pied_pn_ps_all

pi_table_hypoleuca_macro <- subset(pi_table_hypoleuca, new_category == "Macro")
pi_table_hypoleuca_macro <- na.omit(pi_table_hypoleuca_macro)
pied_pn_ps_all_macro <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ new_category == "Macro", pi_table_hypoleuca_macro, mean)
pied_pn_ps_all_macro$pied_pn_ps_all <- with(pied_pn_ps_all_macro, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))
pied_pn_ps_all_macro <- mean(pied_pn_ps_all_macro$pied_pn_ps_all)

pi_table_hypoleuca_micro <- subset(pi_table_hypoleuca, new_category == "Micro")
pi_table_hypoleuca_micro <- na.omit(pi_table_hypoleuca_micro)
pied_pn_ps_all_micro <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ new_category == "Micro", pi_table_hypoleuca_micro, mean)
pied_pn_ps_all_micro$pied_pn_ps_all <- with(pied_pn_ps_all_micro, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))
pied_pn_ps_all_micro <- mean(pied_pn_ps_all_micro$pied_pn_ps_all)

pied_pn_ps_all <- data.frame(pied_pn_ps_all = c(pied_pn_ps_all_macro, pied_pn_ps_all_micro), row.names = c("Macro", "Micro"))

# taiga_pn_ps_gc

pi_table_taiga_macro <- subset(pi_table_taiga, new_category == "Macro")
pi_table_taiga_macro <- na.omit(pi_table_taiga_macro)
taiga_pn_ps_gc_macro <- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ new_category == "Macro", pi_table_taiga_macro, mean)
taiga_pn_ps_gc_macro$taiga_pn_ps_gc <- with(taiga_pn_ps_gc_macro, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))
taiga_pn_ps_gc_macro <- mean(taiga_pn_ps_gc_macro$taiga_pn_ps_gc)

pi_table_taiga_micro <- subset(pi_table_taiga, category == "Micro")
pi_table_taiga_micro <- na.omit(pi_table_taiga_micro)
taiga_pn_ps_gc_micro <- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ category == "Micro", pi_table_taiga_micro, mean)
taiga_pn_ps_gc_micro$taiga_pn_ps_gc <- with(taiga_pn_ps_gc_micro, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))
taiga_pn_ps_gc_micro <- mean(taiga_pn_ps_gc_micro$taiga_pn_ps_gc)

taiga_pn_ps_gc <- data.frame(taiga_pn_ps_gc = c(taiga_pn_ps_gc_macro, taiga_pn_ps_gc_micro), row.names = c("Macro", "Micro"))

# taiga_pn_ps_all

pi_table_taiga_macro <- subset(pi_table_taiga, new_category == "Macro")
pi_table_taiga_macro <- na.omit(pi_table_taiga_macro)
taiga_pn_ps_all_macro <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ new_category == "Macro", pi_table_taiga_macro, mean)
taiga_pn_ps_all_macro$taiga_pn_ps_all <- with(taiga_pn_ps_all_macro, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))
taiga_pn_ps_all_macro <- mean(taiga_pn_ps_all_macro$taiga_pn_ps_all)

pi_table_taiga_micro <- subset(pi_table_taiga, new_category == "Micro")
pi_table_taiga_micro <- na.omit(pi_table_taiga_micro)
taiga_pn_ps_all_micro <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ new_category == "Micro", pi_table_taiga_micro, mean)
taiga_pn_ps_all_micro$taiga_pn_ps_all <- with(taiga_pn_ps_all_micro, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))
taiga_pn_ps_all_micro <- mean(taiga_pn_ps_all_micro$taiga_pn_ps_all)

taiga_pn_ps_all <- data.frame(taiga_pn_ps_all = c(taiga_pn_ps_all_macro, taiga_pn_ps_all_micro), row.names = c("Macro", "Micro"))

# parva_pn_ps_gc

pi_table_parva_macro <- subset(pi_table_parva, new_category == "Macro")
pi_table_parva_macro <- na.omit(pi_table_parva_macro)
parva_pn_ps_gc_macro <- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ new_category == "Macro", pi_table_parva_macro, mean)
parva_pn_ps_gc_macro$parva_pn_ps_gc <- with(parva_pn_ps_gc_macro, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))
parva_pn_ps_gc_macro <- mean(parva_pn_ps_gc_macro$parva_pn_ps_gc)

pi_table_parva_micro <- subset(pi_table_parva, new_category == "Micro")
pi_table_parva_micro <- na.omit(pi_table_parva_micro)
parva_pn_ps_gc_micro <- aggregate(cbind(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc) ~ new_category == "Micro", pi_table_parva_micro, mean)
parva_pn_ps_gc_micro$parva_pn_ps_gc <- with(parva_pn_ps_gc_micro, calculer_pn_ps(het_zero_gc, len_zero_gc, het_four_gc, len_four_gc))
parva_pn_ps_gc_micro <- mean(parva_pn_ps_gc_micro$parva_pn_ps_gc)

parva_pn_ps_gc <- data.frame(parva_pn_ps_gc = c(parva_pn_ps_gc_macro, parva_pn_ps_gc_micro), row.names = c("Macro", "Micro"))

# parva_pn_ps_all

pi_table_parva_macro <- subset(pi_table_parva, new_category == "Macro")
pi_table_parva_macro <- na.omit(pi_table_parva_macro)
parva_pn_ps_all_macro <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ new_category == "Macro", pi_table_parva_macro, mean)
parva_pn_ps_all_macro$parva_pn_ps_all <- with(parva_pn_ps_all_macro, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))
parva_pn_ps_all_macro <- mean(parva_pn_ps_all_macro$parva_pn_ps_all)

pi_table_parva_micro <- subset(pi_table_parva, new_category == "Micro")
pi_table_parva_micro <- na.omit(pi_table_parva_micro)
parva_pn_ps_all_micro <- aggregate(cbind(het_zero_all, len_zero_all, het_four_all, len_four_all) ~ new_category == "Micro", pi_table_parva_micro, mean)
parva_pn_ps_all_micro$parva_pn_ps_all <- with(parva_pn_ps_all_micro, calculer_pn_ps(het_zero_all, len_zero_all, het_four_all, len_four_all))
parva_pn_ps_all_micro <- mean(parva_pn_ps_all_micro$parva_pn_ps_all)

parva_pn_ps_all <- data.frame(parva_pn_ps_all = c(parva_pn_ps_all_macro, parva_pn_ps_all_micro), row.names = c("Macro", "Micro"))


pnpscat<-data.frame(c("Macro", "Micro"))
names(pnpscat)<-c("new_category")

pn_ps_cat<-cbind(pnpscat, coll_pn_ps_gc, coll_pn_ps_all, pied_pn_ps_gc, pied_pn_ps_all, taiga_pn_ps_gc, taiga_pn_ps_all, parva_pn_ps_gc, parva_pn_ps_all)
pn_ps_cat$new_category <- factor(pn_ps_cat$new_category, levels = c("Macro", "Micro"))

palette_couleurs <- c("pn/ps GC coll" = "steelblue1", 
                      "pn/ps GC pied" = "darkorchid",
                      "pn/ps GC parva" = "tomato",
                      "pn/ps GC taiga" = "gold",
                      "pn/ps ALL coll" = "steelblue1",
                      "pn/ps ALL pied" = "darkorchid",
                      "pn/ps ALL parva" = "tomato",
                      "pn/ps ALL taiga" = "gold")

linetype_palette <- c("pn/ps GC coll" = "dashed",
                      "pn/ps GC pied" = "dashed",
                      "pn/ps GC parva" = "dashed",
                      "pn/ps GC taiga" = "dashed",
                      "pn/ps ALL coll" = "solid",
                      "pn/ps ALL pied" = "solid",
                      "pn/ps ALL parva" = "solid",
                      "pn/ps ALL taiga" = "solid")

ggplot(pn_ps_cat) +
  geom_point(aes(x = new_category, y = coll_pn_ps_all, color = "pn/ps ALL coll"), shape = 16, size = 3) +
  geom_point(aes(x = new_category, y = pied_pn_ps_all, color = "pn/ps ALL pied"), shape = 16, size = 3) +
  geom_point(aes(x = new_category, y = parva_pn_ps_all, color = "pn/ps ALL parva"), shape = 16, size = 3) +
  geom_point(aes(x = new_category, y = taiga_pn_ps_all, color = "pn/ps ALL taiga"), shape = 16, size = 3) +
  geom_point(aes(x = new_category, y = coll_pn_ps_gc, color = "pn/ps GC coll"), shape = 17, size = 3) +
  geom_point(aes(x = new_category, y = pied_pn_ps_gc, color = "pn/ps GC pied"), shape = 17, size = 3) +
  geom_point(aes(x = new_category, y = parva_pn_ps_gc, color = "pn/ps GC parva"), shape = 17, size = 3) +
  geom_point(aes(x = new_category, y = taiga_pn_ps_gc, color = "pn/ps GC taiga"), shape = 17, size = 3) +
  geom_line(aes(x = new_category, y = coll_pn_ps_all, group = 1, color = "pn/ps ALL coll", linetype ="pn/ps ALL coll"), size = 1)+
  geom_line(aes(x = new_category, y = pied_pn_ps_all, group = 1, color = "pn/ps ALL pied", linetype ="pn/ps ALL pied"), size = 1)+
  geom_line(aes(x = new_category, y = parva_pn_ps_all, group = 1, color = "pn/ps ALL parva", linetype ="pn/ps ALL parva"), size = 1)+
  geom_line(aes(x = new_category, y = taiga_pn_ps_all, group = 1, color = "pn/ps ALL taiga", linetype ="pn/ps ALL taiga"), size = 1)+
  geom_line(aes(x = new_category, y = coll_pn_ps_gc, group = 1, color = "pn/ps GC coll", linetype ="pn/ps GC coll"), size = 1)+
  geom_line(aes(x = new_category, y = pied_pn_ps_gc, group = 1, color = "pn/ps GC pied", linetype ="pn/ps GC pied"), size = 1)+
  geom_line(aes(x = new_category, y = parva_pn_ps_gc, group = 1, color = "pn/ps GC parva", linetype ="pn/ps GC parva"), size = 1)+
  geom_line(aes(x = new_category, y = taiga_pn_ps_gc, group = 1, color = "pn/ps GC taiga", linetype ="pn/ps GC taiga"), size = 1)+
  scale_color_manual(values = palette_couleurs) +
  scale_linetype_manual(values = linetype_palette) +
  labs(x = "chromosome category", y = "pn/ps", color = "", linetype = "") +
  ggtitle("") +
  theme(
    legend.title = element_text(size = 15),          
    legend.text = element_text(size = 14),          
    axis.title.x = element_text(size = 16),          
    axis.title.y = element_text(size = 16),        
    axis.text.x = element_text(size = 16),           
    axis.text.y = element_text(size = 12)            
  )


y###################################################################################################

# distribution piN

setwd("/home/amaniouloux/Documents/data/dataframe//")
output_directory <- "/home/amaniouloux/Documents/data/Rstudio/résultats/"

size<-read.table("chr_category.txt", header = T)
size<-size %>% rename(chr = chromosome)

pi_table<- pi_table[, c('chr','albicollis_pn_gc', 'hypoleuca_pn_gc', 'taiga_pn_gc', 'parva_pn_gc', 'albicollis_pn_all', 'hypoleuca_pn_all','taiga_pn_all','parva_pn_all')]

pi_table$chr <- factor(pi_table$chr, levels = mixedsort(pi_table$chr))
pi_table<-merge(pi_table, size, by ="chr", all =T)
pi_table <- pi_table %>% filter(chr != "ChrFal34")
pi_table <- pi_table %>% filter(chr != "ChrFal36")
pi_table <- pi_table %>% filter(chr != "ChrZ")

palette_couleurs <- c("pn GC coll" = "steelblue1", 
                      "pn GC pied" = "darkorchid",
                      "pn GC parva" = "tomato",
                      "pn GC taiga" = "gold",
                      "pn ALL coll" = "steelblue1",
                      "pn ALL pied" = "darkorchid",
                      "pn ALL parva" = "tomato",
                      "pn ALL taiga" = "gold")
ggplot() +
  geom_point(data = pi_table, aes(x = reorder(chr, length), y = albicollis_pn_gc, color = "pn GC coll"), shape = 16, size = 2) +
  geom_point(data = pi_table, aes(x = reorder(chr, length), y = hypoleuca_pn_gc, color = "pn GC pied"), shape = 16, size = 2) +
  geom_point(data = pi_table, aes(x = reorder(chr, length), y = parva_pn_gc, color = "pn GC parva"), shape = 16, size = 2) +
  geom_point(data = pi_table, aes(x = reorder(chr, length), y = taiga_pn_gc, color = "pn GC taiga"), shape = 16, size = 2) +
  geom_point(data = pi_table, aes(x = reorder(chr, length), y = albicollis_pn_all, color = "pn ALL coll"), shape = 17, size = 2) +
  geom_point(data = pi_table, aes(x = reorder(chr, length), y = hypoleuca_pn_all, color = "pn ALL pied"), shape = 17, size = 2) +
  geom_point(data = pi_table, aes(x = reorder(chr, length), y = parva_pn_all, color = "pn ALL parva"), shape = 17, size = 2) +
  geom_point(data = pi_table, aes(x = reorder(chr, length), y = taiga_pn_all, color = "pn ALL taiga"), shape = 17, size = 2) +
  scale_color_manual(values = palette_couleurs) +
  labs(x = "chromosome", y = "ps", color = "Conditions") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  

###########################################################################################################

# distribution piN

ssetwd("/home/amaniouloux/Documents/data/dataframe//")
output_directory <- "/home/amaniouloux/Documents/data/Rstudio/résultats/"

size<-read.table("chr_category.txt", header = T)
size<-size %>% rename(chr = chromosome)

pi_table<- pi_table[, c('chr','albicollis_ps_gc', 'hypoleuca_ps_gc', 'taiga_ps_gc', 'parva_ps_gc', 'albicollis_ps_all', 'hypoleuca_ps_all','taiga_ps_all','parva_ps_all')]

pi_table$chr <- factor(pi_table$chr, levels = mixedsort(pi_table$chr))
pi_table<-merge(pi_table, size, by ="chr", all =T)
pi_table <- pi_table %>% filter(chr != "ChrFal34")
pi_table <- pi_table %>% filter(chr != "ChrFal36")
pi_table <- pi_table %>% filter(chr != "ChrZ")

palette_couleurs <- c("ps GC coll" = "steelblue1", 
                      "ps GC pied" = "darkorchid",
                      "ps GC parva" = "tomato",
                      "ps GC taiga" = "gold",
                      "ps ALL coll" = "steelblue1",
                      "ps ALL pied" = "darkorchid",
                      "ps ALL parva" = "tomato",
                      "ps ALL taiga" = "gold")
ggplot() +
  geom_point(data = pi_table, aes(x = reorder(chr, length), y = albicollis_ps_gc, color = "ps GC coll"), shape = 16, size = 2) +
  geom_point(data = pi_table, aes(x = reorder(chr, length), y = hypoleuca_ps_gc, color = "ps GC pied"), shape = 16, size = 2) +
  geom_point(data = pi_table, aes(x = reorder(chr, length), y = parva_ps_gc, color = "ps GC parva"), shape = 16, size = 2) +
  geom_point(data = pi_table, aes(x = reorder(chr, length), y = taiga_ps_gc, color = "ps GC taiga"), shape = 16, size = 2) +
  geom_point(data = pi_table, aes(x = reorder(chr, length), y = albicollis_ps_all, color = "ps ALL coll"), shape = 17, size = 2) +
  geom_point(data = pi_table, aes(x = reorder(chr, length), y = hypoleuca_ps_all, color = "ps ALL pied"), shape = 17, size = 2) +
  geom_point(data = pi_table, aes(x = reorder(chr, length), y = parva_ps_all, color = "ps ALL parva"), shape = 17, size = 2) +
  geom_point(data = pi_table, aes(x = reorder(chr, length), y = taiga_ps_all, color = "ps ALL taiga"), shape = 17, size = 2) +
  scale_color_manual(values = palette_couleurs) +
  labs(x = "chromosome", y = "ps", color = "Conditions") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  

#################

# gene count

taiga_pnps_gc<- aggregate(pn_ps_gc ~ chr, data = table_taiga, FUN = function(x) sum(!is.na(x)))
taiga_pnps_all<- aggregate(pn_ps_all ~ chr, data = table_taiga, FUN = function(x) sum(!is.na(x)))
parva_pnps_gc<-aggregate(pn_ps_gc ~ chr, data = table_parva, FUN = function(x) sum(!is.na(x)))
parva_pnps_all<- aggregate(pn_ps_all ~ chr, data = table_parva, FUN = function(x) sum(!is.na(x)))
pied_pnps_gc<- aggregate(pn_ps_gc ~ chr, data = table_hypoleuca, FUN = function(x) sum(!is.na(x)))
pied_pnps_all<- aggregate(pn_ps_all ~ chr, data = table_hypoleuca, FUN = function(x) sum(!is.na(x)))
coll_pnps_gc<- aggregate(pn_ps_gc ~ chr, data = table_albicollis, FUN = function(x) sum(!is.na(x)))
coll_pnps_all<- aggregate(pn_ps_all ~ chr, data = table_albicollis, FUN = function(x) sum(!is.na(x)))

count_pnps <- taiga_pnps_gc %>%
  left_join(taiga_pnps_all, by = "chr") %>%
  left_join(parva_pnps_gc, by = "chr") %>%
  left_join(parva_pnps_all, by = "chr") %>%
  left_join(pied_pnps_gc, by = "chr") %>%
  left_join(pied_pnps_all, by = "chr") %>%
  left_join(coll_pnps_gc, by = "chr") %>%
  left_join(coll_pnps_all, by = "chr")
  


