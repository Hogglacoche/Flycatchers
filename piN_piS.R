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
table_parva <-read.table("table_parva.txt", header = T)


pi_table_albicollis<- table_albicollis[, c('gene', 'chr', 'category', 'new_category', 'het_zero_gc', 'len_zero_gc', 'het_four_gc', 'len_four_gc', 'het_zero_all','len_zero_all','het_four_all','len_four_all')]
pi_table_albicollis[, -c(1, 2,3,4)] <- lapply(pi_table_albicollis[, -c(1, 2,3,4)], as.numeric)
pi_table_albicollis<-na.omit(pi_table_albicollis)

pi_table_taiga<- table_taiga[, c('gene','het_zero_gc', 'len_zero_gc', 'het_four_gc', 'len_four_gc', 'het_zero_all','len_zero_all','het_four_all','len_four_all')]
pi_table_taiga[, -c(1)] <- lapply(pi_table_taiga[, -c(1)], as.numeric)
pi_table_taiga<-na.omit(pi_table_taiga)

pi_table_hypoleuca<- table_hypoleuca[, c('gene', 'het_zero_gc', 'len_zero_gc', 'het_four_gc', 'len_four_gc', 'het_zero_all','len_zero_all','het_four_all','len_four_all')]
pi_table_hypoleuca[, -c(1)] <- lapply(pi_table_hypoleuca[, -c(1)], as.numeric)
pi_table_hypoleuca<-na.omit(pi_table_hypoleuca)

pi_table_parva<- table_parva[, c('gene', 'het_zero_gc', 'len_zero_gc', 'het_four_gc', 'len_four_gc', 'het_zero_all','len_zero_all','het_four_all','len_four_all')]
pi_table_parva[, -c(1)] <- lapply(pi_table_parva[, -c(1)], as.numeric)
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
 
pi_table<-merge(pi_table_albicollis, pi_table_hypoleuca, by = 'gene', all = T)
pi_table <- pi_table %>%
  rename_with(~ gsub("\\.x$", "_coll", .x), ends_with(".x"))
pi_table <- pi_table %>%
  rename_with(~ gsub("\\.y$", "_pied", .x), ends_with(".y"))
pi_table<-merge(pi_table, pi_table_taiga, by = 'gene', all = T)
pi_table<-merge(pi_table, pi_table_parva, by = 'gene', all = T)
pi_table <- pi_table %>%
  rename_with(~ gsub("\\.x$", "_taiga", .x), ends_with(".x"))
pi_table <- pi_table %>%
  rename_with(~ gsub("\\.y$", "_parva", .x), ends_with(".y"))
pi_table<-na.omit(pi_table)


### chr

## albicollis_pn_ps_gc
albicollis_pn_ps_gc <- aggregate(cbind(het_zero_gc_coll, len_zero_gc_coll, het_four_gc_coll, len_four_gc_coll) ~ chr, pi_table, mean)
albicollis_pn_ps_gc$albicollis_pn_gc <- with(albicollis_pn_ps_gc, calculer_pn(het_zero_gc_coll, len_zero_gc_coll))
albicollis_pn_ps_gc$albicollis_ps_gc <- with(albicollis_pn_ps_gc, calculer_ps(het_four_gc_coll, len_four_gc_coll))
albicollis_pn_ps_gc$albicollis_pn_ps_gc <- with(albicollis_pn_ps_gc, calculer_pn_ps(het_zero_gc_coll, len_zero_gc_coll, het_four_gc_coll, len_four_gc_coll))

## albicollis_pn_ps_all
albicollis_pn_ps_all <- aggregate(cbind(het_zero_all_coll, len_zero_all_coll, het_four_all_coll, len_four_all_coll) ~ chr, pi_table, mean)
albicollis_pn_ps_all$albicollis_pn_all <- with(albicollis_pn_ps_all, calculer_pn(het_zero_all_coll, len_zero_all_coll))
albicollis_pn_ps_all$albicollis_ps_all <- with(albicollis_pn_ps_all, calculer_ps(het_four_all_coll, len_four_all_coll))
albicollis_pn_ps_all$albicollis_pn_ps_all <- with(albicollis_pn_ps_all, calculer_pn_ps(het_zero_all_coll, len_zero_all_coll, het_four_all_coll, len_four_all_coll))

## taiga_pn_ps_gc
taiga_pn_ps_gc <- aggregate(cbind(het_zero_gc_taiga, len_zero_gc_taiga, het_four_gc_taiga, len_four_gc_taiga) ~ chr, pi_table, mean)
taiga_pn_ps_gc$taiga_pn_gc <- with(taiga_pn_ps_gc, calculer_pn(het_zero_gc_taiga, len_zero_gc_taiga))
taiga_pn_ps_gc$taiga_ps_gc <- with(taiga_pn_ps_gc, calculer_ps(het_four_gc_taiga, len_four_gc_taiga))
taiga_pn_ps_gc$taiga_pn_ps_gc <- with(taiga_pn_ps_gc, calculer_pn_ps(het_zero_gc_taiga, len_zero_gc_taiga, het_four_gc_taiga, len_four_gc_taiga))

## taiga_pn_ps_all
taiga_pn_ps_all <- aggregate(cbind(het_zero_all_taiga, len_zero_all_taiga, het_four_all_taiga, len_four_all_taiga) ~ chr, pi_table, mean)
taiga_pn_ps_all$taiga_pn_all <- with(taiga_pn_ps_all, calculer_pn(het_zero_all_taiga, len_zero_all_taiga))
taiga_pn_ps_all$taiga_ps_all <- with(taiga_pn_ps_all, calculer_ps(het_four_all_taiga, len_four_all_taiga))
taiga_pn_ps_all$taiga_pn_ps_all <- with(taiga_pn_ps_all, calculer_pn_ps(het_zero_all_taiga, len_zero_all_taiga, het_four_all_taiga, len_four_all_taiga))

## hypoleuca_pn_ps_gc
hypoleuca_pn_ps_gc <- aggregate(cbind(het_zero_gc_pied, len_zero_gc_pied, het_four_gc_pied, len_four_gc_pied) ~ chr, pi_table, mean)
hypoleuca_pn_ps_gc$hypoleuca_pn_gc <- with(hypoleuca_pn_ps_gc, calculer_pn(het_zero_gc_pied, len_zero_gc_pied))
hypoleuca_pn_ps_gc$hypoleuca_ps_gc <- with(hypoleuca_pn_ps_gc, calculer_ps(het_four_gc_pied, len_four_gc_pied))
hypoleuca_pn_ps_gc$hypoleuca_pn_ps_gc <- with(hypoleuca_pn_ps_gc, calculer_pn_ps(het_zero_gc_pied, len_zero_gc_pied, het_four_gc_pied, len_four_gc_pied))

## hypoleuca_pn_ps_all
hypoleuca_pn_ps_all <- aggregate(cbind(het_zero_all_pied, len_zero_all_pied, het_four_all_pied, len_four_all_pied) ~ chr, pi_table, mean)
hypoleuca_pn_ps_all$hypoleuca_pn_all <- with(hypoleuca_pn_ps_all, calculer_pn(het_zero_all_pied, len_zero_all_pied))
hypoleuca_pn_ps_all$hypoleuca_ps_all <- with(hypoleuca_pn_ps_all, calculer_ps(het_four_all_pied, len_four_all_pied))
hypoleuca_pn_ps_all$hypoleuca_pn_ps_all <- with(hypoleuca_pn_ps_all, calculer_pn_ps(het_zero_all_pied, len_zero_all_pied, het_four_all_pied, len_four_all_pied))

## parva_pn_ps_gc
parva_pn_ps_gc <- aggregate(cbind(het_zero_gc_parva, len_zero_gc_parva, het_four_gc_parva, len_four_gc_parva) ~ chr, pi_table, mean)
parva_pn_ps_gc$parva_pn_gc <- with(parva_pn_ps_gc, calculer_pn(het_zero_gc_parva, len_zero_gc_parva))
parva_pn_ps_gc$parva_ps_gc <- with(parva_pn_ps_gc, calculer_ps(het_four_gc_parva, len_four_gc_parva))
parva_pn_ps_gc$parva_pn_ps_gc <- with(parva_pn_ps_gc, calculer_pn_ps(het_zero_gc_parva, len_zero_gc_parva, het_four_gc_parva, len_four_gc_parva))

## parva_pn_ps_all
parva_pn_ps_all <- aggregate(cbind(het_zero_all_parva, len_zero_all_parva, het_four_all_parva, len_four_all_parva) ~ chr, pi_table, mean)
parva_pn_ps_all$parva_pn_all <- with(parva_pn_ps_all, calculer_pn(het_zero_all_parva, len_zero_all_parva))
parva_pn_ps_all$parva_ps_all <- with(parva_pn_ps_all, calculer_ps(het_four_all_parva, len_four_all_parva))
parva_pn_ps_all$parva_pn_ps_all <- with(parva_pn_ps_all, calculer_pn_ps(het_zero_all_parva, len_zero_all_parva, het_four_all_parva, len_four_all_parva))

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

pi_table<- pi_table[, c('chr','albicollis_pn_ps_all','albicollis_pn_ps_gc', 'hypoleuca_pn_ps_all', 'hypoleuca_pn_ps_gc','taiga_pn_ps_all','taiga_pn_ps_gc', 'parva_pn_ps_all','parva_pn_ps_gc')]

pi_table <- pi_table %>% filter(chr != "Chr25")
pi_table <- pi_table %>% filter(chr != "Chr28")
pi_table <- pi_table %>% filter(chr != "ChrLGE22")

pi_table<-subset(pi_table, select= -chr)
pi_table<- sapply(pi_table, as.numeric)

is.data.frame(pi_table) 
is.numeric(pi_table)                           
correlation_matrix<-cor(pi_table, method = "spearman")
correlation_matrix

new_labels <- c("COLL ALL", "COLL GC", "PIED ALL", "PIED GC", "TAIGA ALL", "TAIGA GC", "PARVA ALL", "PARVA GC")

colnames(correlation_matrix) <- new_labels
rownames(correlation_matrix) <- new_labels

corrplot(correlation_matrix, method = "circle", type = "upper", 
				 tl.col = "black", tl.srt = 45, addCoef.col = "black")

col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = correlation_matrix, col = col, symm = TRUE)

cor.test(pi_table$albicollis_pn_ps_gc, pi_table$taiga_pn_ps_gc, method = "spearman")

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


##############################################################################################################""

# pn/ps new category macro  micro

pi_table_macro <- subset(pi_table, new_category == "Macro")
pi_table_macro<-na.omit(pi_table_macro)

pi_table_micro <- subset(pi_table, new_category == "Micro")
pi_table_micro<-na.omit(pi_table_micro)

# coll_pn_ps_gc

coll_pn_ps_gc_macro<- aggregate(cbind(het_zero_gc_coll, len_zero_gc_coll, het_four_gc_coll, len_four_gc_coll) ~ new_category == "Macro", pi_table_macro, mean)
coll_pn_ps_gc_macro$coll_pn_ps_gc<- with(coll_pn_ps_gc_macro, calculer_pn_ps(het_zero_gc_coll, len_zero_gc_coll, het_four_gc_coll, len_four_gc_coll))
coll_pn_ps_gc_macro<-mean(coll_pn_ps_gc_macro$coll_pn_ps_gc)

coll_pn_ps_gc_micro<- aggregate(cbind(het_zero_gc_coll, len_zero_gc_coll, het_four_gc_coll, len_four_gc_coll) ~ new_category == "Micro", pi_table_micro, mean)
coll_pn_ps_gc_micro$coll_pn_ps_gc<- with(coll_pn_ps_gc_micro, calculer_pn_ps(het_zero_gc_coll, len_zero_gc_coll, het_four_gc_coll, len_four_gc_coll))
coll_pn_ps_gc_micro<-mean(coll_pn_ps_gc_micro$coll_pn_ps_gc)

coll_pn_ps_gc<- data.frame(coll_pn_ps_gc= c(coll_pn_ps_gc_macro, coll_pn_ps_gc_micro),row.names = c("Macro", "Micro"))

# coll_pn_ps_all

coll_pn_ps_all_macro <- aggregate(cbind(het_zero_all_coll, len_zero_all_coll, het_four_all_coll, len_four_all_coll) ~ new_category == "Macro", pi_table_macro, mean)
coll_pn_ps_all_macro$coll_pn_ps_all <- with(coll_pn_ps_all_macro, calculer_pn_ps(het_zero_all_coll, len_zero_all_coll, het_four_all_coll, len_four_all_coll))
coll_pn_ps_all_macro <- mean(coll_pn_ps_all_macro$coll_pn_ps_all)

coll_pn_ps_all_micro <- aggregate(cbind(het_zero_all_coll, len_zero_all_coll, het_four_all_coll, len_four_all_coll) ~ new_category == "Micro", pi_table_micro, mean)
coll_pn_ps_all_micro$coll_pn_ps_all <- with(coll_pn_ps_all_micro, calculer_pn_ps(het_zero_all_coll, len_zero_all_coll, het_four_all_coll, len_four_all_coll))
coll_pn_ps_all_micro <- mean(coll_pn_ps_all_micro$coll_pn_ps_all)

coll_pn_ps_all <- data.frame(coll_pn_ps_all = c(coll_pn_ps_all_macro, coll_pn_ps_all_micro), row.names = c("Macro", "Micro"))

# pied_dn_ds_gc

pied_pn_ps_gc_macro<- aggregate(cbind(het_zero_gc_pied, len_zero_gc_pied, het_four_gc_pied, len_four_gc_pied) ~ new_category == "Macro", pi_table_macro, mean)
pied_pn_ps_gc_macro$pied_pn_ps_gc<- with(pied_pn_ps_gc_macro, calculer_pn_ps(het_zero_gc_pied, len_zero_gc_pied, het_four_gc_pied, len_four_gc_pied))
pied_pn_ps_gc_macro<-mean(pied_pn_ps_gc_macro$pied_pn_ps_gc)

pied_pn_ps_gc_micro<- aggregate(cbind(het_zero_gc_pied, len_zero_gc_pied, het_four_gc_pied, len_four_gc_pied) ~ new_category == "Micro", pi_table_micro, mean)
pied_pn_ps_gc_micro$pied_pn_ps_gc<- with(pied_pn_ps_gc_micro, calculer_pn_ps(het_zero_gc_pied, len_zero_gc_pied, het_four_gc_pied, len_four_gc_pied))
pied_pn_ps_gc_micro<-mean(pied_pn_ps_gc_micro$pied_pn_ps_gc)

pied_pn_ps_gc<- data.frame(pied_pn_ps_gc= c(pied_pn_ps_gc_macro, pied_pn_ps_gc_micro),row.names = c("Macro", "Micro"))

# pied_pn_ps_all

pied_pn_ps_all_macro <- aggregate(cbind(het_zero_all_pied, len_zero_all_pied, het_four_all_pied, len_four_all_pied) ~ new_category == "Macro", pi_table_macro, mean)
pied_pn_ps_all_macro$pied_pn_ps_all <- with(pied_pn_ps_all_macro, calculer_pn_ps(het_zero_all_pied, len_zero_all_pied, het_four_all_pied, len_four_all_pied))
pied_pn_ps_all_macro <- mean(pied_pn_ps_all_macro$pied_pn_ps_all)

pied_pn_ps_all_micro <- aggregate(cbind(het_zero_all_pied, len_zero_all_pied, het_four_all_pied, len_four_all_pied) ~ new_category == "Micro", pi_table_micro, mean)
pied_pn_ps_all_micro$pied_pn_ps_all <- with(pied_pn_ps_all_micro, calculer_pn_ps(het_zero_all_pied, len_zero_all_pied, het_four_all_pied, len_four_all_pied))
pied_pn_ps_all_micro <- mean(pied_pn_ps_all_micro$pied_pn_ps_all)

pied_pn_ps_all <- data.frame(pied_pn_ps_all = c(pied_pn_ps_all_macro, pied_pn_ps_all_micro), row.names = c("Macro", "Micro"))

# taiga_pn_ps_gc

taiga_pn_ps_gc_macro <- aggregate(cbind(het_zero_gc_taiga, len_zero_gc_taiga, het_four_gc_taiga, len_four_gc_taiga) ~ new_category == "Macro", pi_table_macro, mean)
taiga_pn_ps_gc_macro$taiga_pn_ps_gc <- with(taiga_pn_ps_gc_macro, calculer_pn_ps(het_zero_gc_taiga, len_zero_gc_taiga, het_four_gc_taiga, len_four_gc_taiga))
taiga_pn_ps_gc_macro <- mean(taiga_pn_ps_gc_macro$taiga_pn_ps_gc)

taiga_pn_ps_gc_micro <- aggregate(cbind(het_zero_gc_taiga, len_zero_gc_taiga, het_four_gc_taiga, len_four_gc_taiga) ~ new_category == "Micro", pi_table_micro, mean)
taiga_pn_ps_gc_micro$taiga_pn_ps_gc <- with(taiga_pn_ps_gc_micro, calculer_pn_ps(het_zero_gc_taiga, len_zero_gc_taiga, het_four_gc_taiga, len_four_gc_taiga))
taiga_pn_ps_gc_micro <- mean(taiga_pn_ps_gc_micro$taiga_pn_ps_gc)

taiga_pn_ps_gc <- data.frame(taiga_pn_ps_gc = c(taiga_pn_ps_gc_macro, taiga_pn_ps_gc_micro), row.names = c("Macro", "Micro"))

# taiga_pn_ps_all

taiga_pn_ps_all_macro <- aggregate(cbind(het_zero_all_taiga, len_zero_all_taiga, het_four_all_taiga, len_four_all_taiga) ~ new_category == "Macro", pi_table_macro, mean)
taiga_pn_ps_all_macro$taiga_pn_ps_all <- with(taiga_pn_ps_all_macro, calculer_pn_ps(het_zero_all_taiga, len_zero_all_taiga, het_four_all_taiga, len_four_all_taiga))
taiga_pn_ps_all_macro <- mean(taiga_pn_ps_all_macro$taiga_pn_ps_all)

taiga_pn_ps_all_micro <- aggregate(cbind(het_zero_all_taiga, len_zero_all_taiga, het_four_all_taiga, len_four_all_taiga) ~ new_category == "Micro", pi_table_micro, mean)
taiga_pn_ps_all_micro$taiga_pn_ps_all <- with(taiga_pn_ps_all_micro, calculer_pn_ps(het_zero_all_taiga, len_zero_all_taiga, het_four_all_taiga, len_four_all_taiga))
taiga_pn_ps_all_micro <- mean(taiga_pn_ps_all_micro$taiga_pn_ps_all)

taiga_pn_ps_all <- data.frame(taiga_pn_ps_all = c(taiga_pn_ps_all_macro, taiga_pn_ps_all_micro), row.names = c("Macro", "Micro"))

# parva_pn_ps_gc

parva_pn_ps_gc_macro <- aggregate(cbind(het_zero_gc_parva, len_zero_gc_parva, het_four_gc_parva, len_four_gc_parva) ~ new_category == "Macro", pi_table_macro, mean)
parva_pn_ps_gc_macro$parva_pn_ps_gc <- with(parva_pn_ps_gc_macro, calculer_pn_ps(het_zero_gc_parva, len_zero_gc_parva, het_four_gc_parva, len_four_gc_parva))
parva_pn_ps_gc_macro <- mean(parva_pn_ps_gc_macro$parva_pn_ps_gc)

parva_pn_ps_gc_micro <- aggregate(cbind(het_zero_gc_parva, len_zero_gc_parva, het_four_gc_parva, len_four_gc_parva) ~ new_category == "Micro", pi_table_micro, mean)
parva_pn_ps_gc_micro$parva_pn_ps_gc <- with(parva_pn_ps_gc_micro, calculer_pn_ps(het_zero_gc_parva, len_zero_gc_parva, het_four_gc_parva, len_four_gc_parva))
parva_pn_ps_gc_micro <- mean(parva_pn_ps_gc_micro$parva_pn_ps_gc)

parva_pn_ps_gc <- data.frame(parva_pn_ps_gc = c(parva_pn_ps_gc_macro, parva_pn_ps_gc_micro), row.names = c("Macro", "Micro"))

# parva_pn_ps_all

parva_pn_ps_all_macro <- aggregate(cbind(het_zero_all_parva, len_zero_all_parva, het_four_all_parva, len_four_all_parva) ~ new_category == "Macro", pi_table_macro, mean)
parva_pn_ps_all_macro$parva_pn_ps_all <- with(parva_pn_ps_all_macro, calculer_pn_ps(het_zero_all_parva, len_zero_all_parva, het_four_all_parva, len_four_all_parva))
parva_pn_ps_all_macro <- mean(parva_pn_ps_all_macro$parva_pn_ps_all)

parva_pn_ps_all_micro <- aggregate(cbind(het_zero_all_parva, len_zero_all_parva, het_four_all_parva, len_four_all_parva) ~ new_category == "Micro", pi_table_micro, mean)
parva_pn_ps_all_micro$parva_pn_ps_all <- with(parva_pn_ps_all_micro, calculer_pn_ps(het_zero_all_parva, len_zero_all_parva, het_four_all_parva, len_four_all_parva))
parva_pn_ps_all_micro <- mean(parva_pn_ps_all_micro$parva_pn_ps_all)

parva_pn_ps_all <- data.frame(parva_pn_ps_all = c(parva_pn_ps_all_macro, parva_pn_ps_all_micro), row.names = c("Macro", "Micro"))


pnpscat<-data.frame(c("Macro", "Micro"))
names(pnpscat)<-c("new_category")

pn_ps_cat<-cbind(pnpscat, coll_pn_ps_gc, coll_pn_ps_all, pied_pn_ps_gc, pied_pn_ps_all, taiga_pn_ps_gc, taiga_pn_ps_all, parva_pn_ps_gc, parva_pn_ps_all)
pn_ps_cat$new_category <- factor(pn_ps_cat$new_category, levels = c("Macro", "Micro"))

palette_couleurs <- c("pn/ps GC coll" = "tomato", 
                      "pn/ps GC pied" = "darkorchid",
                      "pn/ps GC parva" = "tomato",
                      "pn/ps GC taiga" = "gold",
                      "pn/ps ALL coll" = "steelblue4",
                      "pn/ps ALL pied" = "darkorchid",
                      "pn/ps ALL parva" = "tomato",
                      "pn/ps ALL taiga" = "gold")

linetype_palette <- c("pn/ps GC coll" = "solid",
                      "pn/ps GC pied" = "dashed",
                      "pn/ps GC parva" = "dashed",
                      "pn/ps GC taiga" = "dashed",
                      "pn/ps ALL coll" = "solid",
                      "pn/ps ALL pied" = "solid",
                      "pn/ps ALL parva" = "solid",
                      "pn/ps ALL taiga" = "solid")

ggplot(pn_ps_cat) +
  geom_point(aes(x = new_category, y = coll_pn_ps_all, color = "pn/ps ALL coll"), shape = 16, size = 3) +
  geom_point(aes(x = new_category, y = coll_pn_ps_gc, color = "pn/ps GC coll"), shape = 17, size = 3) +
  geom_line(aes(x = new_category, y = coll_pn_ps_all, group = 1, color = "pn/ps ALL coll", linetype ="pn/ps ALL coll"), size = 2)+
  geom_line(aes(x = new_category, y = coll_pn_ps_gc, group = 1, color = "pn/ps GC coll", linetype ="pn/ps GC coll"), size = 2)+
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
  

#############################################################################

# ANOVA


pi_table_albicollis<- table_albicollis[, c('gene', 'chr', 'new_category', 'het_zero_gc', 'len_zero_gc', 'het_four_gc', 'len_four_gc', 'het_zero_all','len_zero_all','het_four_all','len_four_all')]
pi_table_albicollis[, -c(1, 2,3,4)] <- lapply(pi_table_albicollis[, -c(1, 2,3,4)], as.numeric)
pi_table_albicollis<-na.omit(pi_table_albicollis)

pi_table_taiga<- table_taiga[, c('gene','het_zero_gc', 'len_zero_gc', 'het_four_gc', 'len_four_gc', 'het_zero_all','len_zero_all','het_four_all','len_four_all')]
pi_table_taiga[, -c(1)] <- lapply(pi_table_taiga[, -c(1)], as.numeric)
pi_table_taiga<-na.omit(pi_table_taiga)

pi_table_hypoleuca<- table_hypoleuca[, c('gene', 'het_zero_gc', 'len_zero_gc', 'het_four_gc', 'len_four_gc', 'het_zero_all','len_zero_all','het_four_all','len_four_all')]
pi_table_hypoleuca[, -c(1)] <- lapply(pi_table_hypoleuca[, -c(1)], as.numeric)
pi_table_hypoleuca<-na.omit(pi_table_hypoleuca)

pi_table_parva<- table_parva[, c('gene', 'het_zero_gc', 'len_zero_gc', 'het_four_gc', 'len_four_gc', 'het_zero_all','len_zero_all','het_four_all','len_four_all')]
pi_table_parva[, -c(1)] <- lapply(pi_table_parva[, -c(1)], as.numeric)
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

pi_table<-merge(pi_table_albicollis, pi_table_hypoleuca, by = 'gene', all = T)
pi_table <- pi_table %>%
	rename_with(~ gsub("\\.x$", "_coll", .x), ends_with(".x"))
pi_table <- pi_table %>%
	rename_with(~ gsub("\\.y$", "_pied", .x), ends_with(".y"))
pi_table<-merge(pi_table, pi_table_taiga, by = 'gene', all = T)
pi_table<-merge(pi_table, pi_table_parva, by = 'gene', all = T)
pi_table <- pi_table %>%
	rename_with(~ gsub("\\.x$", "_taiga", .x), ends_with(".x"))
pi_table <- pi_table %>%
	rename_with(~ gsub("\\.y$", "_parva", .x), ends_with(".y"))
pi_table<-na.omit(pi_table)

# coll_pn_ps_gc

coll_pn_ps_gc<- aggregate(cbind(het_zero_gc_coll, len_zero_gc_coll, het_four_gc_coll, len_four_gc_coll) ~ gene, pi_table, mean)
coll_pn_ps_gc$coll_pn_ps_gc<- with(coll_pn_ps_gc, calculer_pn_ps(het_zero_gc_coll, len_zero_gc_coll, het_four_gc_coll, len_four_gc_coll))

coll_pn_ps_gc<-merge(coll_pn_ps_gc, table_albicollis, by = "gene", all = T)
coll_pn_ps_gc<-na.omit(coll_pn_ps_gc)
coll_pn_ps_gc<- coll_pn_ps_gc[coll_pn_ps_gc$coll_pn_ps_gc >= 0 & coll_pn_ps_gc$coll_pn_ps_gc <= 100, ]

# coll_pn_ps_all

coll_pn_ps_all <- aggregate(cbind(het_zero_all_coll, len_zero_all_coll, het_four_all_coll, len_four_all_coll) ~ gene, pi_table, mean)
coll_pn_ps_all$coll_pn_ps_all <- with(coll_pn_ps_all, calculer_pn_ps(het_zero_all_coll, len_zero_all_coll, het_four_all_coll, len_four_all_coll))

coll_pn_ps_all<-merge(coll_pn_ps_all, table_albicollis, by = "gene", all = T)
coll_pn_ps_all<-na.omit(coll_pn_ps_all)
coll_pn_ps_all<- coll_pn_ps_all[coll_pn_ps_all$coll_pn_ps_all >= 0 & coll_pn_ps_all$coll_pn_ps_all <= 100, ]


# pied_dn_ds_gc

pied_pn_ps_gc<- aggregate(cbind(het_zero_gc_pied, len_zero_gc_pied, het_four_gc_pied, len_four_gc_pied) ~ gene, pi_table, mean)
pied_pn_ps_gc$pied_pn_ps_gc<- with(pied_pn_ps_gc, calculer_pn_ps(het_zero_gc_pied, len_zero_gc_pied, het_four_gc_pied, len_four_gc_pied))

pied_pn_ps_gc<-merge(pied_pn_ps_gc, table_hypoleuca, by = "gene", all = T)
pied_pn_ps_gc<-na.omit(pied_pn_ps_gc)
pied_pn_ps_gc<- pied_pn_ps_gc[pied_pn_ps_gc$pied_pn_ps_gc >= 0 & pied_pn_ps_gc$pied_pn_ps_gc <= 100, ]


# pied_pn_ps_all

pied_pn_ps_all<- aggregate(cbind(het_zero_all_pied, len_zero_all_pied, het_four_all_pied, len_four_all_pied) ~ gene, pi_table, mean)
pied_pn_ps_all$pied_pn_ps_all <- with(pied_pn_ps_all, calculer_pn_ps(het_zero_all_pied, len_zero_all_pied, het_four_all_pied, len_four_all_pied))

pied_pn_ps_all<-merge(pied_pn_ps_all, table_hypoleuca, by = "gene", all = T)
pied_pn_ps_all<-na.omit(pied_pn_ps_all)
pied_pn_ps_all<- pied_pn_ps_all[pied_pn_ps_all$pied_pn_ps_all >= 0 & pied_pn_ps_all$pied_pn_ps_all <= 100, ]

# taiga_pn_ps_gc

taiga_pn_ps_gc<- aggregate(cbind(het_zero_gc_taiga, len_zero_gc_taiga, het_four_gc_taiga, len_four_gc_taiga) ~ gene, pi_table, mean)
taiga_pn_ps_gc$taiga_pn_ps_gc <- with(taiga_pn_ps_gc, calculer_pn_ps(het_zero_gc_taiga, len_zero_gc_taiga, het_four_gc_taiga, len_four_gc_taiga))

taiga_pn_ps_gc<-merge(taiga_pn_ps_gc, table_taiga, by = "gene", all = T)
taiga_pn_ps_gc<-na.omit(taiga_pn_ps_gc)
taiga_pn_ps_gc<- taiga_pn_ps_gc[taiga_pn_ps_gc$taiga_pn_ps_gc >= 0 & taiga_pn_ps_gc$taiga_pn_ps_gc <= 100, ]


# taiga_pn_ps_all

taiga_pn_ps_all <- aggregate(cbind(het_zero_all_taiga, len_zero_all_taiga, het_four_all_taiga, len_four_all_taiga) ~ gene, pi_table, mean)
taiga_pn_ps_all$taiga_pn_ps_all <- with(taiga_pn_ps_all, calculer_pn_ps(het_zero_all_taiga, len_zero_all_taiga, het_four_all_taiga, len_four_all_taiga))

taiga_pn_ps_all<-merge(taiga_pn_ps_all, table_taiga, by = "gene", all = T)
taiga_pn_ps_all<-na.omit(taiga_pn_ps_all)
taiga_pn_ps_all<- taiga_pn_ps_all[taiga_pn_ps_all$taiga_pn_ps_all >= 0 & taiga_pn_ps_all$taiga_pn_ps_all <= 100, ]


# parva_pn_ps_gc

parva_pn_ps_gc <- aggregate(cbind(het_zero_gc_parva, len_zero_gc_parva, het_four_gc_parva, len_four_gc_parva) ~ gene, pi_table, mean)
parva_pn_ps_gc$parva_pn_ps_gc <- with(parva_pn_ps_gc, calculer_pn_ps(het_zero_gc_parva, len_zero_gc_parva, het_four_gc_parva, len_four_gc_parva))

parva_pn_ps_gc<-merge(parva_pn_ps_gc, table_parva, by = "gene", all = T)
parva_pn_ps_gc<-na.omit(parva_pn_ps_gc)
parva_pn_ps_gc<- parva_pn_ps_gc[parva_pn_ps_gc$parva_pn_ps_gc >= 0 & parva_pn_ps_gc$parva_pn_ps_gc <= 100, ]


# parva_pn_ps_all

parva_pn_ps_all<- aggregate(cbind(het_zero_all_parva, len_zero_all_parva, het_four_all_parva, len_four_all_parva) ~ gene, pi_table, mean)
parva_pn_ps_all$parva_pn_ps_all <- with(parva_pn_ps_all, calculer_pn_ps(het_zero_all_parva, len_zero_all_parva, het_four_all_parva, len_four_all_parva))

parva_pn_ps_all<-merge(parva_pn_ps_all, table_parva, by = "gene", all = T)
parva_pn_ps_all<-na.omit(parva_pn_ps_all)
parva_pn_ps_all<- parva_pn_ps_all[parva_pn_ps_all$parva_pn_ps_all >= 0 & parva_pn_ps_all$parva_pn_ps_all <= 100, ]


pi_anova <- lm(coll_pn_ps_all ~ new_category, data = coll_pn_ps_all)
anova_result <- anova(pi_anova)
print(anova_result)

pi_anova <- lm(pied_pn_ps_all ~ new_category, data = pied_pn_ps_all)
anova_result <- anova(pi_anova)
print(anova_result)

pi_anova <- lm(parva_pn_ps_all ~ new_category, data = parva_pn_ps_all)
anova_result <- anova(pi_anova)
print(anova_result)

pi_anova <- lm(taiga_pn_ps_all ~ new_category, data = taiga_pn_ps_all)
anova_result <- anova(pi_anova)
print(anova_result)


pi_anova <- lm(coll_pn_ps_gc ~ new_category, data = coll_pn_ps_gc)
anova_result <- anova(pi_anova)
print(anova_result)

pi_anova <- lm(pied_pn_ps_gc ~ new_category, data = pied_pn_ps_gc)
anova_result <- anova(pi_anova)
print(anova_result)

pi_anova <- lm(parva_pn_ps_gc ~ new_category, data = parva_pn_ps_gc)
anova_result <- anova(pi_anova)
print(anova_result)

pi_anova <- lm(taiga_pn_ps_gc ~ new_category, data = taiga_pn_ps_gc)
anova_result <- anova(pi_anova)
print(anova_result)


merged_pn_ps <- coll_pn_ps_all %>%
	left_join(coll_pn_ps_gc, by = "gene")

table_albicollis <- read.table("table_albicollis.txt", header = T) 
table_albicollis <- table_albicollis[table_albicollis$chr != "Chr22", ]
table_albicollis <- table_albicollis[table_albicollis$chr != "Chr25", ]
table_albicollis <- table_albicollis[table_albicollis$chr != "Chr28", ]
setwd("/home/amaniouloux/Documents/data/Rstudio/Recombinaison/")
recomb<- read.table("coll.ld_recom_converted.gene_recom_no_cpg_filt.recom_bin_inf.recom_cons_inf.bed", header = T) 
recomb$gene_id <- sapply(str_split(recomb$gene_id, "_"), "[[", 1)
recomb<-recomb %>% rename(gene = gene_id)

table_anova <- merge(merged_pn_ps, table_albicollis, by = "gene", all = T)
table_anova <- table_anova[, c('gene', 'chr','density', 'new_category', 'density', 'length','coll_pn_ps_all', 'coll_pn_ps_gc', 'TPM_Brain_coll_mean', 'TPM_Testis_coll_mean', 'TPM_Heart_coll_mean', 'TPM_Liver_coll_mean', 'TPM_Kidney_coll_mean')]
table_anova <- merge(table_anova, recomb, by = "gene", all = T)
tau<-tau_results_coll_gene %>% rename(gene = Gene)
table_anova <- merge(table_anova, tau, by = "gene", all = T)

table_anova<-na.omit(table_anova)

cat<-as.factor(table_anova$new_category)
is.numeric(table_anova$gene_recom_coll)

####
library(car)
vif(lm(coll_pn_ps_gc ~ density+gene_recom_coll + cat   + Tau+TPM_Brain_coll_mean+TPM_Testis_coll_mean+TPM_Heart_coll_mean+TPM_Liver_coll_mean+TPM_Kidney_coll_mean, data = table_anova))
coll_pn_ps_gc <- lm(coll_pn_ps_gc ~ density+gene_recom_coll + density +Tau+cat+TPM_Brain_coll_mean+TPM_Testis_coll_mean+TPM_Heart_coll_mean+TPM_Liver_coll_mean+TPM_Kidney_coll_mean, data = table_anova)
anova(coll_pn_ps_gc)
Anova(coll_pn_ps_gc, type = "II")

coll_pn_ps_all <- lm(coll_pn_ps_all ~ gene_recom_coll + density+cat + Tau , data = table_anova)
anova(coll_pn_ps_all)
Anova(coll_pn_ps_all, type = "II")



taiga_dn_ds_gc_S <- aggregate(cbind(dn_tot_gc_S_taiga, dn_norm_tot_gc_S_taiga, ds_tot_gc_S_taiga, ds_norm_tot_gc_S_taiga) ~ gene, table, mean)
taiga_dn_ds_gc_S$taiga_dn_ds_gc_S <- with(taiga_dn_ds_gc_S, calculer_dn_ds(dn_tot_gc_S_taiga, dn_norm_tot_gc_S_taiga, ds_tot_gc_S_taiga, ds_norm_tot_gc_S_taiga))

## taiga_dn_ds_all.S

taiga_dn_ds_all_S <- aggregate(cbind(dn_tot_all_S_taiga, dn_norm_tot_all_S_taiga, ds_tot_all_S_taiga, ds_norm_tot_all_S_taiga) ~ gene, table, mean)
taiga_dn_ds_all_S$taiga_dn_ds_all_S<- with(taiga_dn_ds_all_S, calculer_dn_ds(dn_tot_all_S_taiga, dn_norm_tot_all_S_taiga, ds_tot_all_S_taiga, ds_norm_tot_all_S_taiga))

## merge 
merged_dn_ds <- merge(taiga_dn_ds_all_S, taiga_dn_ds_gc_S, by = "gene")
table_taiga <- read.table("table_taiga.txt", header = T) 
table_taiga <- table_taiga[table_taiga$chr != "Chr22", ]
table_taiga <- table_taiga[table_taiga$chr != "Chr25", ]
table_taiga <- table_taiga[table_taiga$chr != "Chr28", ]
setwd("/home/amaniouloux/Documents/data/Rstudio/Recombinaison/")
recomb<- read.table("taig.ld_recom_converted.gene_recom_no_cpg_filt.recom_bin_inf.recom_cons_inf.bed", header = T) 
recomb$gene_id <- sapply(str_split(recomb$gene_id, "_"), "[[", 1)
recomb<-recomb %>% rename(gene = gene_id)

table_anova <- select(merged_dn_ds, gene, taiga_dn_ds_gc_S, taiga_dn_ds_all_S)
table_anova <- merge(table_anova, table_taiga, by = "gene", all = T)
table_anova <- table_anova[, c('gene', 'chr', 'new_category', 'length','taiga_dn_ds_gc_S', 'taiga_dn_ds_all_S')]
table_anova <- merge(table_anova, recomb, by = "gene", all = T)

cat<-as.factor(table_anova$new_category)
recom_bin_taiga<-as.factor(table_anova$recom_bin_taiga)
as.numeric(table_anova$gene_recom_taiga)
chr<-unique(table_anova$chr)

##
taiga_dn_ds_all_S <- lm(taiga_dn_ds_all_S ~ length + gene_recom_taig + recom_bin_taig+ cat + chr + chr:gene_recom_taig, data = table_anova)
anova(taiga_dn_ds_all_S)

taiga_dn_ds_gc_S <- lm(taiga_dn_ds_gc_S ~ length + gene_recom_taig + recom_bin_taig+ cat + chr + chr : gene_recom_taig, data = table_anova)
anova(taiga_dn_ds_gc_S)
##




