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
table_albicollis$gene <- sapply(str_split(table_albicollis$gene, "_"), "[[", 1)
table_hypoleuca$gene <- sapply(str_split(table_hypoleuca$gene, "_"), "[[", 1)
table_parva$gene <- sapply(str_split(table_parva$gene, "_"), "[[", 1)
table_taiga$gene <- sapply(str_split(table_taiga$gene, "_"), "[[", 1)

pi_table_albicollis<- table_albicollis[, c('chr', 'het_zero_gc', 'len_zero_gc', 'het_four_gc', 'len_four_gc', 'het_zero_all','len_zero_all','het_four_all','len_four_all')]
pi_table_albicollis[, -1] <- sapply(pi_table_albicollis[, -1], as.numeric)
pi_table_albicollis<-na.omit(pi_table_albicollis)

pi_table_taiga<- table_taiga[, c('chr', 'het_zero_gc', 'len_zero_gc', 'het_four_gc', 'len_four_gc', 'het_zero_all','len_zero_all','het_four_all','len_four_all')]
pi_table_taiga[, -1] <- sapply(pi_table_taiga[, -1], as.numeric)
pi_table_taiga<-na.omit(pi_table_taiga)

pi_table_hypoleuca<- table_hypoleuca[, c('chr', 'het_zero_gc', 'len_zero_gc', 'het_four_gc', 'len_four_gc', 'het_zero_all','len_zero_all','het_four_all','len_four_all')]
pi_table_hypoleuca[, -1] <- sapply(pi_table_hypoleuca[, -1], as.numeric)
pi_table_hypoleuca<-na.omit(pi_table_hypoleuca)

pi_table_parva<- table_parva[, c('chr', 'het_zero_gc', 'len_zero_gc', 'het_four_gc', 'len_four_gc', 'het_zero_all','len_zero_all','het_four_all','len_four_all')]
pi_table_parva[, -1] <- sapply(pi_table_parva[, -1], as.numeric)
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

## hypoleuca_pn_ps_all
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
# correlation matrix gene 

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
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Orientation verticale
#################################################################################

# albicollis GC pn~pS

pi_table_albicollis <- read.table("table_albicollis.txt", header = T)
pi_table_albicollis<- pi_table_albicollis[, c('chr', 'pn_gc', 'pn_all', 'ps_gc', 'ps_all', 'pn_ps_gc', 'pn_ps_all')]
pi_table_albicollis<-na.omit(pi_table_albicollis)
pi_table_albicollis[pi_table_albicollis == "NULL"] <- 0
pi_table_albicollis$pn_gc<-as.numeric(pi_table_albicollis$pn_gc)
pi_table_albicollis$ps_gc<-as.numeric(pi_table_albicollis$ps_gc)

pi_table_albicollis <- pi_table_albicollis[mixedorder(pi_table_albicollis$chr), ]
groupes <- split(pi_table_albicollis, pi_table_albicollis$chr)

piN_gc <- groupes$Chr1[2] 
piN_gc<- sapply(piN_gc, as.numeric)
piS_gc<- groupes$Chr1[4] 
piS_gc<- sapply(piS_gc, as.numeric)
cor(piN_gc, piS_gc)

piN_all <- groupes$Chr1[3] 
piN_all<- sapply(piN_all, as.numeric)
piS_all<- groupes$Chr1[5] 
piS_all<- sapply(piS_all, as.numeric)
cor(piN_all, piS_all)

resultats_corr <- list()

for (chr in names(groupes)) {
  piN <- as.numeric(groupes[[chr]][[2]]) 
  piS <- as.numeric(groupes[[chr]][[4]])  
  cor_chr <- cor(piN, piS)
    resultats_corr[[chr]] <- cor_chr
}
chromosomes <- names(resultats_corr)
extract_numbers <- function(x) {
  as.numeric(gsub("[^0-9]", "", x))
}

chromosomes <- setdiff(chromosomes, c("Chr22", "ChrLGE22","Chr23"))

index_a_retirer1 <- which(names(resultats_corr) == "Chr23")
resultats_corr <- resultats_corr[-index_a_retirer1]
index_a_retirer <- which(is.na(resultats_corr))
resultats_corr <- resultats_corr[-index_a_retirer]


index_ord <- order(sapply(chromosomes, extract_numbers))
resultats_corr_ord <- resultats_corr[index_ord]

resultats_df <- do.call(rbind, resultats_corr_ord)
resultats_df <- data.frame(Chromosome = names(resultats_corr_ord), Correlation = unlist(resultats_corr_ord))
resultats_df$Chromosome <- factor(resultats_df$Chromosome, levels = unique(resultats_df$Chromosome))

chromosome_numbers <- sapply(chromosomes, extract_numbers)
sorted_chromosomes <- chromosomes[order(chromosome_numbers)]
resultats_df$Chromosome <- factor(resultats_df$Chromosome, levels = sorted_chromosomes)

ggplot(resultats_df, aes(x = Chromosome, y = Correlation, color = Color)) +
  geom_point(size = 3) +  # Ajouter les points
  geom_smooth(method = "lm", se = FALSE) + 
  labs(title = "Correlation between pN and pS for each chromosome (Ficedula albicollis / GC) ",
       x = "Chromosome",
       y = "Correlation") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +

########################################
# albicollis ALL pn~pS

pi_table_albicollis <- read.table("table_albicollis.txt", header = T)
pi_table_albicollis<- pi_table_albicollis[, c('chr', 'pn_gc', 'pn_all', 'ps_gc', 'ps_all', 'pn_ps_gc', 'pn_ps_all')]
pi_table_albicollis<-na.omit(pi_table_albicollis)
pi_table_albicollis[pi_table_albicollis == "NULL"] <- 0
pi_table_albicollis$pn_gc<-as.numeric(pi_table_albicollis$pn_gc)
pi_table_albicollis$ps_gc<-as.numeric(pi_table_albicollis$ps_gc)

pi_table_albicollis <- pi_table_albicollis[mixedorder(pi_table_albicollis$chr), ]
groupes <- split(pi_table_albicollis, pi_table_albicollis$chr)

piN_gc <- groupes$Chr1[2] 
piN_gc<- sapply(piN_gc, as.numeric)
piS_gc<- groupes$Chr1[4] 
piS_gc<- sapply(piS_gc, as.numeric)
cor(piN_gc, piS_gc)

piN_all <- groupes$Chr1[3] 
piN_all<- sapply(piN_all, as.numeric)
piS_all<- groupes$Chr1[5] 
piS_all<- sapply(piS_all, as.numeric)
cor(piN_all, piS_all)

resultats_corr <- list()

for (chr in names(groupes)) {
  piN <- as.numeric(groupes[[chr]][[3]]) 
  piS <- as.numeric(groupes[[chr]][[5]])  
  cor_chr <- cor(piN, piS)
  resultats_corr[[chr]] <- cor_chr
}
chromosomes <- names(resultats_corr)
extract_numbers <- function(x) {
  as.numeric(gsub("[^0-9]", "", x))
}

chromosomes <- setdiff(chromosomes, c("ChrLGE22"))

index_a_retirer1 <- which(names(resultats_corr) == "ChLGE22")
resultats_corr <- resultats_corr[-index_a_retirer1]


index_ord <- order(sapply(chromosomes, extract_numbers))
resultats_corr_ord <- resultats_corr[index_ord]

resultats_df <- do.call(rbind, resultats_corr_ord)
resultats_df <- data.frame(Chromosome = names(resultats_corr_ord), Correlation = unlist(resultats_corr_ord))
resultats_df$Chromosome <- factor(resultats_df$Chromosome, levels = unique(resultats_df$Chromosome))

chromosome_numbers <- sapply(chromosomes, extract_numbers)
sorted_chromosomes <- chromosomes[order(chromosome_numbers)]
resultats_df$Chromosome <- factor(resultats_df$Chromosome, levels = sorted_chromosomes)

ggplot(resultats_df, aes(x = Chromosome, y = Correlation, color = Color)) +
  geom_point(size = 3) +  # Ajouter les points
  geom_smooth(method = "lm", se = FALSE) + 
  labs(title = "Correlation between pN and pS for each chromosome (Ficedula albicollis / GC) ",
       x = "Chromosome",
       y = "Correlation") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values = c("red", "blue"))


########################################
# albicollis pn/ps
setwd("/home/amaniouloux/Documents/data/dataframe//")

pi_table_albicollis <- read.table("table_albicollis.txt", header = T)
pi_table_albicollis<- pi_table_albicollis[, c('chr', 'pn_ps_gc', 'pn_ps_all')]
pi_table_albicollis<-na.omit(pi_table_albicollis)
pi_table_albicollis[pi_table_albicollis == "NULL"] <- 0
pi_table_albicollis$pn_ps_gc<- sapply(pi_table_albicollis$pn_ps_gc, as.numeric)
pi_table_albicollis$pn_ps_gc-as.numeric(pi_table_albicollis$pn_ps_gc)
pi_table_albicollis$pn_ps_all<- sapply(pi_table_albicollis$pn_ps_all, as.numeric)
pi_table_albicollis$pn_ps_all<-as.numeric(pi_table_albicollis$pn_ps_all)

pi_table_albicollis <- pi_table_albicollis[mixedorder(pi_table_albicollis$chr), ]
groupes <- split(pi_table_albicollis, pi_table_albicollis$chr)

resultats_corr <- list()

for (chr in names(groupes)) {
  pn_ps_gc <- as.numeric(groupes[[chr]][[2]]) 
  pn_ps_all <- as.numeric(groupes[[chr]][[3]])  
  cor_chr <- cor(pn_ps_gc, pn_ps_all, method = "spearman")
  resultats_corr[[chr]] <- cor_chr
}
chromosomes <- names(resultats_corr)
extract_numbers <- function(x) {
  as.numeric(gsub("[^0-9]", "", x))
}

chromosomes <- setdiff(chromosomes, c("Chr22", "ChrLGE22","Chr24", "Chr25"))

index_a_retirer <- which(is.na(resultats_corr))
resultats_corr <- resultats_corr[-index_a_retirer]


index_ord <- order(sapply(chromosomes, extract_numbers))
resultats_corr_ord <- resultats_corr[index_ord]

resultats_df <- do.call(rbind, resultats_corr_ord)
resultats_df <- data.frame(Chromosome = names(resultats_corr_ord), Correlation = unlist(resultats_corr_ord))
resultats_df$Chromosome <- factor(resultats_df$Chromosome, levels = unique(resultats_df$Chromosome))

chromosome_numbers <- sapply(chromosomes, extract_numbers)
sorted_chromosomes <- chromosomes[order(chromosome_numbers)]
resultats_df$Chromosome <- factor(resultats_df$Chromosome, levels = sorted_chromosomes)

ggplot(resultats_df, aes(x = Chromosome, y = Correlation)) +
  geom_point(color = "black", size = 3) +  # Ajouter les points
  geom_smooth(method = "lm", se = FALSE) + 
  labs(title = "Correlation between pN/pS GC and pN/pS all sites for each chromosome (Ficedula albicollis) ",
       x = "Chromosome",
       y = "Correlation") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 


#############################################################################
# hypoleuca pn/ps

setwd("/home/amaniouloux/Documents/data/dataframe//")

pi_table_hypoleuca <- read.table("table_hypoleuca.txt", header = T)
pi_table_hypoleuca<- pi_table_hypoleuca[, c('chr', 'pn_ps_gc', 'pn_ps_all')]
pi_table_hypoleuca<-na.omit(pi_table_hypoleuca)
pi_table_hypoleuca[pi_table_hypoleuca == "NULL"] <- 0
pi_table_hypoleuca$pn_ps_gc<- sapply(pi_table_hypoleuca$pn_ps_gc, as.numeric)
pi_table_hypoleuca$pn_ps_gc-as.numeric(pi_table_hypoleuca$pn_ps_gc)
pi_table_hypoleuca$pn_ps_all<- sapply(pi_table_hypoleuca$pn_ps_all, as.numeric)
pi_table_hypoleuca$pn_ps_all<-as.numeric(pi_table_hypoleuca$pn_ps_all)

pi_table_hypoleuca <- pi_table_hypoleuca[mixedorder(pi_table_hypoleuca$chr), ]
groupes <- split(pi_table_hypoleuca, pi_table_hypoleuca$chr)

resultats_corr <- list()

for (chr in names(groupes)) {
  pn_ps_gc <- as.numeric(groupes[[chr]][[2]]) 
  pn_ps_all <- as.numeric(groupes[[chr]][[3]])  
  cor_chr <- cor(pn_ps_gc, pn_ps_all, method = "spearman")
  resultats_corr[[chr]] <- cor_chr
}
chromosomes <- names(resultats_corr)
extract_numbers <- function(x) {
  as.numeric(gsub("[^0-9]", "", x))
}

chromosomes <- setdiff(chromosomes, c("Chr22", "ChrLGE22","Chr23", "Chr24"))

index_a_retirer <- which(is.na(resultats_corr))
resultats_corr <- resultats_corr[-index_a_retirer]

index_ord <- order(sapply(chromosomes, extract_numbers))
resultats_corr_ord <- resultats_corr[index_ord]

resultats_df <- do.call(rbind, resultats_corr_ord)
resultats_df <- data.frame(Chromosome = names(resultats_corr_ord), Correlation = unlist(resultats_corr_ord))
resultats_df$Chromosome <- factor(resultats_df$Chromosome, levels = unique(resultats_df$Chromosome))

chromosome_numbers <- sapply(chromosomes, extract_numbers)
sorted_chromosomes <- chromosomes[order(chromosome_numbers)]
resultats_df$Chromosome <- factor(resultats_df$Chromosome, levels = sorted_chromosomes)

library(ggplot2)
ggplot(resultats_df, aes(x = Chromosome, y = Correlation)) +
  geom_point(color = "black", size = 3) +  # Ajouter les points
  geom_smooth(method = "lm", se = FALSE) + 
  labs(title = "Correlation between pN/pS GC and pN/pS all sites for each chromosome (Ficedula hypoleuca) ",
       x = "Chromosome",
       y = "Correlation") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

#############################################################################
# taiga pn/ps

setwd("/home/amaniouloux/Documents/data/dataframe//")

pi_table_taiga <- read.table("table_taiga.txt", header = T)
pi_table_taiga<- pi_table_taiga[, c('chr', 'pn_ps_gc', 'pn_ps_all')]
pi_table_taiga<-na.omit(pi_table_taiga)
pi_table_taiga[pi_table_taiga == "NULL"] <- 0
pi_table_taiga$pn_ps_gc<- sapply(pi_table_taiga$pn_ps_gc, as.numeric)
pi_table_taiga$pn_ps_gc-as.numeric(pi_table_taiga$pn_ps_gc)
pi_table_taiga$pn_ps_all<- sapply(pi_table_taiga$pn_ps_all, as.numeric)
pi_table_taiga$pn_ps_all<-as.numeric(pi_table_taiga$pn_ps_all)

pi_table_taiga <- pi_table_taiga[mixedorder(pi_table_taiga$chr), ]
groupes <- split(pi_table_taiga, pi_table_taiga$chr)

resultats_corr <- list()

for (chr in names(groupes)) {
  pn_ps_gc <- as.numeric(groupes[[chr]][[2]]) 
  pn_ps_all <- as.numeric(groupes[[chr]][[3]])  
  cor_chr <- cor(pn_ps_gc, pn_ps_all, method = "spearman")
  resultats_corr[[chr]] <- cor_chr
}
chromosomes <- names(resultats_corr)
extract_numbers <- function(x) {
  as.numeric(gsub("[^0-9]", "", x))
}

chromosomes <- setdiff(chromosomes, c("ChrLGE22"))


index_a_retirer <- which(is.na(resultats_corr))
resultats_corr <- resultats_corr[-index_a_retirer]

index_ord <- order(sapply(chromosomes, extract_numbers))
resultats_corr_ord <- resultats_corr[index_ord]

resultats_df <- do.call(rbind, resultats_corr_ord)
resultats_df <- data.frame(Chromosome = names(resultats_corr_ord), Correlation = unlist(resultats_corr_ord))
resultats_df$Chromosome <- factor(resultats_df$Chromosome, levels = unique(resultats_df$Chromosome))

chromosome_numbers <- sapply(chromosomes, extract_numbers)
sorted_chromosomes <- chromosomes[order(chromosome_numbers)]
resultats_df$Chromosome <- factor(resultats_df$Chromosome, levels = sorted_chromosomes)

ggplot(resultats_df, aes(x = Chromosome, y = Correlation)) +
  geom_point(color = "black", size = 3) +  # Ajouter les points
  geom_smooth(method = "lm", se = FALSE) + 
  labs(title = "Correlation between pN/pS GC and pN/pS all sites for each chromosome (Ficedula taiga) ",
       x = "Chromosome",
       y = "Correlation") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

#############################################################################
# parva pn/ps

setwd("/home/amaniouloux/Documents/data/dataframe//")

library(tidyr)
library(gtools)

pi_table_parva <- read.table("table_parva.txt", header = T)
pi_table_parva<- pi_table_parva[, c('chr', 'pn_ps_gc', 'pn_ps_all')]
pi_table_parva<-na.omit(pi_table_parva)
pi_table_parva[pi_table_parva == "NULL"] <- 0
pi_table_parva$pn_ps_gc<- sapply(pi_table_parva$pn_ps_gc, as.numeric)
pi_table_parva$pn_ps_gc-as.numeric(pi_table_parva$pn_ps_gc)
pi_table_parva$pn_ps_all<- sapply(pi_table_parva$pn_ps_all, as.numeric)
pi_table_parva$pn_ps_all<-as.numeric(pi_table_parva$pn_ps_all)

pi_table_parva <- pi_table_parva[mixedorder(pi_table_parva$chr), ]
groupes <- split(pi_table_parva, pi_table_parva$chr)

resultats_corr <- list()

for (chr in names(groupes)) {
  pn_ps_gc <- as.numeric(groupes[[chr]][[2]]) 
  pn_ps_all <- as.numeric(groupes[[chr]][[3]])  
  cor_chr <- cor(pn_ps_gc, pn_ps_all, method = "spearman")
  resultats_corr[[chr]] <- cor_chr
}
chromosomes <- names(resultats_corr)
extract_numbers <- function(x) {
  as.numeric(gsub("[^0-9]", "", x))
}

chromosomes <- setdiff(chromosomes, c("ChrLGE22", "Chr22", "Chr24", "chr25"))

index_a_retirer <- which(is.na(resultats_corr))
resultats_corr2 <- resultats_corr[-index_a_retirer]

index_ord <- order(sapply(chromosomes, extract_numbers))
resultats_corr_ord <- resultats_corr[index_ord]

resultats_df <- do.call(rbind, resultats_corr_ord)
resultats_df <- data.frame(Chromosome = names(resultats_corr_ord), Correlation = unlist(resultats_corr_ord))
resultats_df$Chromosome <- factor(resultats_df$Chromosome, levels = unique(resultats_df$Chromosome))

chromosome_numbers <- sapply(chromosomes, extract_numbers)
sorted_chromosomes <- chromosomes[order(chromosome_numbers)]
resultats_df$Chromosome <- factor(resultats_df$Chromosome, levels = sorted_chromosomes)

ggplot(resultats_df, aes(x = Chromosome, y = Correlation)) +
  geom_point(color = "black", size = 3) +  # Ajouter les points
  geom_smooth(method = "lm", se = FALSE) + 
  labs(title = "Correlation between pN/pS GC and pN/pS all sites for each chromosome (Ficedula parva) ",
       x = "Chromosome",
       y = "Correlation") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 


#########################################################################################

# cor matrix chr

setwd("/home/amaniouloux/Documents/data/dataframe//")
output_directory <- "/home/amaniouloux/Documents/data/Rstudio/résultats/"

table_albicollis <- read.table("table_albicollis.txt", header = T) 
table_hypoleuca <- read.table("table_hypoleuca.txt", header = T)
table_taiga <- read.table("table_taiga.txt", header = T)
table_parva <-read.table("table_parva.txt", header = T)
table_albicollis$gene <- sapply(str_split(table_albicollis$gene, "_"), "[[", 1)
table_hypoleuca$gene <- sapply(str_split(table_hypoleuca$gene, "_"), "[[", 1)
table_parva$gene <- sapply(str_split(table_parva$gene, "_"), "[[", 1)
table_taiga$gene <- sapply(str_split(table_taiga$gene, "_"), "[[", 1)

pi_table_albicollis<- table_albicollis[, c('gene','chr', 'pn_ps_gc', 'pn_ps_all')]
pi_table_albicollis<-na.omit(pi_table_albicollis)
pi_table_hypoleuca<- table_hypoleuca[, c('gene','chr','pn_ps_gc', 'pn_ps_all')]
pi_table_hypoleuca<-na.omit(pi_table_hypoleuca)
pi_table_taiga<- table_taiga[, c('gene','chr', 'pn_ps_gc', 'pn_ps_all')]
pi_table_taiga<-na.omit(pi_table_taiga)
pi_table_parva<- table_parva[, c('gene','chr','pn_ps_gc', 'pn_ps_all')]
pi_table_parva<- na.omit(pi_table_parva)

pi_table<-merge(pi_table_albicollis, pi_table_hypoleuca, by = "gene", all = T)
pi_table <- pi_table %>%
  rename_with(~ gsub("\\.x$", "_coll", .x), ends_with(".x"))
pi_table <- pi_table %>%
  rename_with(~ gsub("\\.y$", "_pied", .x), ends_with(".y"))
pi_table<-merge(pi_table, pi_table_parva, by = "gene", all = T)
pi_table<-merge(pi_table, pi_table_taiga, by = "gene", all = T)
pi_table <- pi_table %>%
  rename_with(~ gsub("\\.x$", "_parva", .x), ends_with(".x"))
pi_table <- pi_table %>%
  rename_with(~ gsub("\\.y$", "_taiga", .x), ends_with(".y"))
pi_table<-subset(pi_table, select= -chr_pied)
pi_table<-subset(pi_table, select= -chr_taiga)
pi_table<-subset(pi_table, select= -chr_parva)
pi_table<-pi_table %>% rename(chr = chr_coll)
pi_table[pi_table == "NULL"] <- NA
pi_table<-subset(pi_table, select= -gene)

pi_table <- pi_table %>% filter(chr != "ChrLGE22")

cols_to_convert <- names(pi_table)[-1]
pi_table[cols_to_convert] <- lapply(pi_table[cols_to_convert], as.numeric)
is.numeric(pi_table$pn_ps_all_coll)
as.data.frame(pi_table)
is.data.frame(pi_table)
as.character(pi_table$chr)

pi_table_mean_chr <- aggregate(cbind(pn_ps_gc_coll, pn_ps_all_coll, pn_ps_gc_pied, pn_ps_all_pied, pn_ps_gc_taiga, pn_ps_all_taiga, pn_ps_gc_parva, pn_ps_all_parva) ~ chr, pi_table, sum)
pi_table_mean_chr <- aggregate(cbind(pn_ps_gc_coll, pn_ps_all_coll, pn_ps_gc_pied, pn_ps_all_pied, pn_ps_gc_taiga, pn_ps_all_taiga, pn_ps_gc_parva, pn_ps_all_parva) ~ chr, pi_table, mean)
correlation_matrix <- cor(pi_table_mean_chr[, -1], method ="spearman")
corrplot(correlation_matrix, method="circle", type="upper", tl.col="black", tl.srt=45)

plot(pi_table_mean_chr$pn_ps_all_taiga, pi_table_mean_chr$pn_ps_all_coll)

###############################################################################

# pn/ps category

size<- read.table("chr_category.txt", header = T)
size<-size %>% rename(chr = chromosome)
pi_table_cat<-merge(pi_table, size, by = "chr", all = T)

categories <- c("Macro", "Intermediate", "Micro")

pi_table_cat<-na.omit(pi_table_cat)

pi_table_cat <- pi_table_cat %>% filter(chr != "Chr24")
pi_table_cat <- pi_table_cat %>% filter(chr != "Chr25")
pi_table_cat <- pi_table_cat %>% filter(chr != "Chr22")
pi_table_cat <- pi_table_cat %>% filter(chr != "Chr26")
pi_table_cat <- pi_table_cat %>% filter(chr != "Chr23")
pi_table_cat <- pi_table_cat %>% filter(chr != "ChrLGE22")
pi_table_cat <- pi_table_cat %>% filter(chr != "ChrFal34")
pi_table_cat <- pi_table_cat %>% filter(chr != "ChrFal36")
pi_table_cat <- pi_table_cat %>% filter(chr != "Chr27")
pi_table_cat <- pi_table_cat %>% filter(chr != "Chr13")
pi_table_cat <- pi_table_cat %>% filter(chr != "Chr28")

pi_table_cat_mean <- aggregate(cbind(albicollis_pn_ps_gc, albicollis_pn_ps_all, hypoleuca_pn_ps_gc, hypoleuca_pn_ps_all, taiga_pn_ps_gc, taiga_pn_ps_all, parva_pn_ps_gc, parva_pn_ps_all) ~ category, pi_table_cat, mean)

pi_table_cat_mean$category <- factor(pi_table_cat_mean$category, levels = c("Macro", "Intermediate", "Micro"))

palette_couleurs <- c("pn/ps GC coll" = "steelblue1", 
                      "pn/ps GC pied" = "darkorchid",
                      "pn/ps GC parva" = "tomato",
                      "pn/ps GC taiga" = "gold",
                      "pn/ps ALL coll" = "steelblue1",
                      "pn/ps ALL pied" = "darkorchid",
                      "pn/ps ALL parva" = "tomato",
                      "pn/ps ALL taiga" = "gold")
ggplot() +

  geom_point(data = pi_table_cat_mean, aes(x = category, y = albicollis_pn_ps_all, color = "pn/ps ALL coll"), shape = 16, size = 2) +
  geom_point(data = pi_table_cat_mean, aes(x = category, y = hypoleuca_pn_ps_all, color = "pn/ps ALL pied"), shape = 16, size = 2) +
  geom_point(data = pi_table_cat_mean, aes(x = category, y = parva_pn_ps_all, color = "pn/ps ALL parva"), shape = 16, size = 2) +
  geom_point(data = pi_table_cat_mean, aes(x = category, y = taiga_pn_ps_all, color = "pn/ps ALL taiga"), shape = 16, size = 2) +
  
  scale_color_manual(values = palette_couleurs) +
  labs(x = "chromosome", y = "pn/ps", color = "Conditions") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Orientation verticale

###################################################################################################

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
  


