setwd("/home/amaniouloux/Documents/data/dataframe")
table_albicollis<-read.table("table_albicollis.txt", header = T)
table_taiga<-read.table("table_taiga.txt", header = T)
size<-read.table("chr_category.txt", header = T)
size<-size %>% rename(chr = chromosome)

library(ggplot2)

table_albicollis_recomb <- table_albicollis[, c("gene", "chr","length","category","gene_recom_coll","recom_bin_coll","win_cons_pc")]
table_albicollis_recomb<-na.omit(table_albicollis_recomb)
table_albicollis_recomb <- subset(table_albicollis_recomb, chr != "ChrLGE22")
table_albicollis_recomb <- subset(table_albicollis_recomb, chr != "Chr25")

table_taiga_recomb <- table_taiga[, c("gene", "chr","length","category","gene_recom_taig","recom_bin_taig","win_cons_pc")]
table_taiga_recomb<-na.omit(table_taiga_recomb)
table_taiga_recomb <- subset(table_taiga_recomb, chr != "ChrLGE22")
table_taiga_recomb <- subset(table_taiga_recomb, chr != "Chr25")

table_recomb<-merge(table_albicollis_recomb, table_taiga_recomb, by = "gene", all = T)
table_recomb<-subset(table_recomb, select= -chr.y)
table_recomb<-subset(table_recomb, select= -length.y)
table_recomb<-subset(table_recomb, select= -category.y)
table_recomb<-table_recomb %>% rename(chr = chr.x)
table_recomb<-table_recomb %>% rename(length = length.x)
table_recomb<-table_recomb %>% rename(category = category.x)
table_recomb<-na.omit(table_recomb)

#################################################################################################################################

table_taiga_recomb$category <- factor(table_taiga_recomb$category,
                                      levels = c("Macro", "Intermediate", "Micro"))
ggplot(table_taiga_recomb, aes(x = category, y = gene_recom_taig, fill = category)) +
  geom_boxplot() +
  labs(x = "Category", y = "Recombination rate",
       title = "Gene recombination rate for each category (F.albicilla)") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5)) 

table_albicollis_recomb$category <- factor(table_albicollis_recomb$category,
                                           levels = c("Macro", "Intermediate", "Micro"))
ggplot(table_albicollis_recomb, aes(x = category, y = gene_recom_coll, fill = category)) +
  geom_boxplot() +
  labs(x = "Category", y = "Recombination rate",
       title = "Gene recombination rate for each category (F.albicollis)") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5)) 


table_taiga_recomb$chr <- reorder(table_taiga_recomb$chr, -table_taiga_recomb$length)

ggplot(table_taiga_recomb, aes(x = chr, y = gene_recom_taig)) +
  geom_boxplot() +
  labs(x = "Chromosome", y = "Recombination rate",
       title = "Gene recombination rate for each chromosome (F.albicilla)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Orientation verticale

table_albicollis_recomb$chr <- reorder(table_albicollis_recomb$chr, -table_albicollis_recomb$length)
ggplot(table_albicollis_recomb, aes(x = chr, y = gene_recom_coll)) +
  geom_boxplot() +
  labs(x = "Chromosome", y = "Recombination rate",
       title = "Gene recombination rate for each chromosome (F.albicollis)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Orientation verticale

########################################################################################################################################

table_taiga_recomb$recom_bin_taig <- factor(table_taiga_recomb$recom_bin_taig, levels = c("low", "med", "high"))
ggplot(table_taiga_recomb, aes(x = chr, y = gene_recom_taig, color = recom_bin_taig)) +
  geom_point() +
  labs(title = "Recombination for each chromosome category (F.albicilla)",
       x = "Chromosome category",
       y = "Recombination rate",
       color = "Recombination
    rate level") +
  scale_color_manual(values = c("low" = "green", "med" = "blue", "high" = "red"),
                     labels = c("low", "med", "high")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  


table_albicollis_recomb$recom_bin_coll <- factor(table_albicollis_recomb$recom_bin_coll, levels = c("low", "med", "high"))
ggplot(table_albicollis_recomb, aes(x = category, y = gene_recom_coll, color = recom_bin_coll)) +
  geom_point() +
  labs(title = "Recombination for each chromosome category (F.albicollis)",
       x = "Chromosome category",
       y = "Recombination rate",
       color = "Recombination
    rate level") +
  scale_color_manual(values = c("low" = "green", "med" = "blue", "high" = "red"),
                     labels = c("low", "med", "high")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  

#############################################################################################################################
nombre_high <- sum(table_taiga_recomb$recom_bin_taig == "high")
nombre_med <- sum(table_taiga_recomb$recom_bin_taig == "med")
nombre_low <- sum(table_taiga_recomb$recom_bin_taig == "low")
variance_by_group <- aggregate(gene_recom_taig ~ recom_bin_taig, data = table_taiga_recomb, FUN = var)
variance_by_chr <- aggregate(gene_recom_taig ~ chr, data = table_taiga_recomb, FUN = var)

nombre_high <- sum(table_albicollis_recomb$recom_bin_coll == "high")
nombre_med <- sum(table_albicollis_recomb$recom_bin_coll == "med")
nombre_low <- sum(table_albicollis_recomb$recom_bin_coll == "low")
variance_by_group <- aggregate(gene_recom_coll ~ recom_bin_coll, data = table_albicollis_recomb, FUN = var)
variance_by_chr <- aggregate(gene_recom_coll ~ chr, data = table_albicollis_recomb, FUN = var)

#################################################################################################################"

### cor dn/ds 
table_recomb <- table_recomb %>% filter(gene != "ENSFALG00000027500")
table_recomb <- table_recomb %>% filter(gene != "ENSFALG00000013927")
table_recomb <- table_recomb %>% filter(gene != "ENSFALG00000012551")
table_recomb <- table_recomb %>% filter(gene != "ENSFALG00000012551")
plot(table_recomb$gene_recom_coll, table_recomb$dn_ds_gene_gc_S_coll)
cor(table_recomb$gene_recom_coll, table_recomb$dn_ds_gene_gc_S_coll, method = "spearman")
cor(table_recomb$gene_recom_coll, table_recomb$dn_ds_gene_gc_L, method = "spearman")
cor(table_recomb$gene_recom_coll, table_recomb$dn_ds_gene_all_S_coll, method = "spearman")
cor(table_recomb$gene_recom_coll, table_recomb$dn_ds_gene_all_L, method = "spearman")

cor(table_recomb$gene_recom_taig, table_recomb$dn_ds_gene_gc_S_taiga, method = "spearman")
cor(table_recomb$gene_recom_taig, table_recomb$dn_ds_gene_all_S_taiga, method = "spearman")

# no correlation between dn/ds and recombination gene 


### cor dn/ds chr

table_recomb_chr <- aggregate(cbind(gene_recom_coll, gene_recom_taig) ~ chr, table_recomb, mean)

calculer_dn_ds <- function(dn_tot, dn_norm_tot, ds_tot, ds_norm_tot) {
  return((dn_tot / dn_norm_tot) / (ds_tot / ds_norm_tot))
}

## coll_dn_ds_gc_L
table_albicollis <- read.table("table_albicollis.txt", header = T) 
coll_dn_ds_gc_L <- aggregate(cbind(dn_tot_gc_L, dn_norm_tot_gc_L, ds_tot_gc_L, ds_norm_tot_gc_L) ~ chr, table_albicollis, mean)
coll_dn_ds_gc_L$coll_dn_ds_gc_L <- with(coll_dn_ds_gc_L, calculer_dn_ds(dn_tot_gc_L, dn_norm_tot_gc_L, ds_tot_gc_L, ds_norm_tot_gc_L))

## coll_dn_ds_gc_S
table_albicollis <- read.table("table_albicollis.txt", header = T) 
coll_dn_ds_gc_S <- aggregate(cbind(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S) ~ chr, table_albicollis, mean)
coll_dn_ds_gc_S$coll_dn_ds_gc_S <- with(coll_dn_ds_gc_S, calculer_dn_ds(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S))

## coll_dn_ds_all_L
table_albicollis <- read.table("table_albicollis.txt", header = T) 
coll_dn_ds_all_L <- aggregate(cbind(dn_tot_all_L, dn_norm_tot_all_L, ds_tot_all_L, ds_norm_tot_all_L) ~ chr, table_albicollis, mean)
coll_dn_ds_all_L$coll_dn_ds_all_L <- with(coll_dn_ds_all_L, calculer_dn_ds(dn_tot_all_L, dn_norm_tot_all_L, ds_tot_all_L, ds_norm_tot_all_L))

## coll_dn_ds_all_S
table_albicollis <- read.table("table_albicollis.txt", header = T) 
coll_dn_ds_all_S <- aggregate(cbind(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S) ~ chr, table_albicollis, mean)
coll_dn_ds_all_S$coll_dn_ds_all_S <- with(coll_dn_ds_all_S, calculer_dn_ds(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S))

## taiga_dn_ds_gc_S
table_taiga <- read.table("table_taiga.txt", header = T) 
taiga_dn_ds_gc_S <- aggregate(cbind(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S) ~ chr, table_taiga, mean)
taiga_dn_ds_gc_S$taiga_dn_ds_gc_S <- with(taiga_dn_ds_gc_S, calculer_dn_ds(dn_tot_gc_S, dn_norm_tot_gc_S, ds_tot_gc_S, ds_norm_tot_gc_S))

## taiga_dn_ds_all.S
table_taiga <- read.table("table_taiga.txt", header = T) 
taiga_dn_ds_all_S <- aggregate(cbind(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S) ~ chr, table_taiga, mean)
taiga_dn_ds_all_S$taiga_dn_ds_all_S<- with(taiga_dn_ds_all_S, calculer_dn_ds(dn_tot_all_S, dn_norm_tot_all_S, ds_tot_all_S, ds_norm_tot_all_S))

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

table_recomb_dnds<-merge(table_recomb_chr, merged_dn_ds, by = 'chr', all = T)
table_recomb_dnds<-na.omit(table_recomb_dnds)

plot(table_recomb_dnds$gene_recom_coll, table_recomb_dnds$coll_dn_ds_all_L)
cor(table_recomb_dnds$gene_recom_coll, table_recomb_dnds$coll_dn_ds_all_L)

plot(table_recomb_dnds$gene_recom_coll, table_recomb_dnds$coll_dn_ds_all_S)
cor(table_recomb_dnds$gene_recom_coll, table_recomb_dnds$coll_dn_ds_all_S)

plot(table_recomb_dnds$gene_recom_coll, table_recomb_dnds$coll_dn_ds_gc_L)
cor(table_recomb_dnds$gene_recom_coll, table_recomb_dnds$coll_dn_ds_gc_L)

plot(table_recomb_dnds$gene_recom_coll, table_recomb_dnds$coll_dn_ds_gc_S)
cor(table_recomb_dnds$gene_recom_coll, table_recomb_dnds$coll_dn_ds_gc_S)

plot(table_recomb_dnds$gene_recom_taig, table_recomb_dnds$taiga_dn_ds_all_S)
cor(table_recomb_dnds$gene_recom_taig, table_recomb_dnds$taiga_dn_ds_all_S)

plot(table_recomb_dnds$gene_recom_taig, table_recomb_dnds$taiga_dn_ds_gc_S)
cor(table_recomb_dnds$gene_recom_taig, table_recomb_dnds$taiga_dn_ds_gc_S)

## no correlation chr ~ recombination + dn/ds short branch
## negative correlation (low) chr ~ recombination + dn/ds long branch

### cor dn/ds chr category

table_recomb_cat <- aggregate(cbind(gene_recom_coll, gene_recom_taig) ~ category, table_recomb, mean)
merged_dn_ds_size<-merge(merged_dn_ds, size, by="chr", all =T )
merged_dn_ds_cat<- aggregate(cbind(coll_dn_ds_all_L, coll_dn_ds_all_S, coll_dn_ds_gc_L, coll_dn_ds_gc_S, taiga_dn_ds_all_S, taiga_dn_ds_gc_S) ~ category, merged_dn_ds_size, mean)

table_recomb_cat<-merge(table_recomb_cat, merged_dn_ds_cat, by="category", all = T)

data <- data.frame(
  gene_recom_coll = table_recomb_cat$gene_recom_coll,
  coll_dn_ds_all_L = table_recomb_cat$coll_dn_ds_all_L,
  coll_dn_ds_all_S = table_recomb_cat$coll_dn_ds_all_S,
  coll_dn_ds_gc_L = table_recomb_cat$coll_dn_ds_gc_L,
  coll_dn_ds_gc_S = table_recomb_cat$coll_dn_ds_gc_S,
  gene_recom_taig = table_recomb_cat$gene_recom_taig,
  taiga_dn_ds_all_S = table_recomb_cat$taiga_dn_ds_all_S,
  taiga_dn_ds_gc_S = table_recomb_cat$taiga_dn_ds_gc_S,
  category = table_recomb_cat$category
)

palette_couleurs <- c("coll dn/ds GC L" = "steelblue4", 
                      "coll dn/ds GC S" = "steelblue1",
                      "coll dn/ds ALL L" = "tomato3",
                      "coll dn/ds ALL S" = "tomato",
                      "taiga dn/ds GC S" = "steelblue1",
                      "taiga dn/ds ALL S" = "tomato")

ggplot(table_recomb_cat) +
  geom_point(aes(x = gene_recom_coll, y = coll_dn_ds_gc_L, color = "coll dn/ds GC L", shape = category), size = 3) +
  geom_point(aes(x = gene_recom_coll, y = coll_dn_ds_gc_S, color = "coll dn/ds GC S", shape = category), size = 3) +
  geom_point(aes(x = gene_recom_coll, y = coll_dn_ds_all_L, color = "coll dn/ds ALL L", shape = category), size = 3) +
  geom_point(aes(x = gene_recom_coll, y = coll_dn_ds_all_S, color = "coll dn/ds ALL S", shape = category), size = 3) +
  
  scale_color_manual(values = palette_couleurs) +
  labs(x = "recombination rate", y = "dn/ds", color = "") +
  ggtitle("dn/ds ~ recombination for each category")


