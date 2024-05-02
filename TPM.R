setwd("/home/amaniouloux/Documents/data/dataframe")
table_albicollis<-read.table("table_albicollis.txt", header =T)
table_hypoleuca<-read.table("table_hypoleuca.txt", header = T)
table_hybrid<-read.table("table_hybrids.txt", header = T)
size<-read.table("chr_category.txt", header = T)

library(ggplot2)
library(gtools)
library(dplyr)
library(corrplot)

# albicollis gene
calculate_tau <- function(table_albicollis) {

  tau_results <- data.frame(Gene = table_albicollis$gene)
  
  tpm_columns <- grep("^TPM_\\w+_coll_mean", names(table_albicollis))
  
  num_tissues <- length(tpm_columns)
  
  max_values <- apply(table_albicollis[tpm_columns], 1, max)
  normalized_table <- table_albicollis[tpm_columns] / max_values
  
  max_col_index <- max.col(normalized_table)
  
  max_tissue <- names(table_albicollis)[tpm_columns][max_col_index]
  
  tau_results$Max_Tissue <- max_tissue
  
  tau_values <- apply(normalized_table, 1, function(row) {
    sum(1 - row) / 4
  })
  
  tau_results$Tau <- tau_values
  
  return(tau_results)
}

tau_results_coll_gene <- calculate_tau(table_albicollis)
tau_results_coll_gene <- na.omit(tau_results_coll_gene)


# hypoleuca gene

calculate_tau <- function(table_hypoleuca) {

  tau_results <- data.frame(Gene = table_hypoleuca$gene)
  
  tpm_columns <- grep("^TPM_\\w+_pied_mean", names(table_hypoleuca))
  
  num_tissues <- length(tpm_columns)
  
  max_values <- apply(table_hypoleuca[tpm_columns], 1, max)
  normalized_table <- table_hypoleuca[tpm_columns] / max_values
  
  max_col_index <- max.col(normalized_table)
  
  max_tissue <- names(table_hypoleuca)[tpm_columns][max_col_index]
  
  tau_results$Max_Tissue <- max_tissue
  
  tau_values <- apply(normalized_table, 1, function(row) {
    sum(1 - row) / 4
  })
  
  tau_results$Tau <- tau_values
  
  return(tau_results)
}

tau_results_pied_gene <- calculate_tau(table_hypoleuca)
tau_results_pied_gene <- na.omit(tau_results_pied_gene)

# hybrids gene

calculate_tau <- function(table_hybrid) {
  
  tau_results <- data.frame(Gene = table_hybrid$gene)
  
  tpm_columns <- grep("^TPM_\\w+_hybrid_mean", names(table_hybrid))
  
  num_tissues <- length(tpm_columns)
  
  max_values <- apply(table_hybrid[tpm_columns], 1, max)
  normalized_table <- table_hybrid[tpm_columns] / max_values
  
  max_col_index <- max.col(normalized_table)
  
  max_tissue <- names(table_hybrid)[tpm_columns][max_col_index]
  
  tau_results$Max_Tissue <- max_tissue
  
  tau_values <- apply(normalized_table, 1, function(row) {
    sum(1 - row) / 4
  })
  
  tau_results$Tau <- tau_values
  
  return(tau_results)
}

tau_results_hybrid_gene <- calculate_tau(table_hybrid)
tau_results_hybrid_gene <- na.omit(tau_results_hybrid_gene)
################################################################################################################################################################

# albicollis chromosome
calculate_tau <- function(table_albicollis) {
  
  tau_results <- data.frame(chr = table_albicollis$chr)
  
  tpm_columns <- grep("^TPM_\\w+_coll_mean", names(table_albicollis))
  
  num_tissues <- length(tpm_columns)
  
  max_values <- apply(table_albicollis[tpm_columns], 1, max)
  normalized_table <- table_albicollis[tpm_columns] / max_values
  
  max_col_index <- max.col(normalized_table)
  
  max_tissue <- names(table_albicollis)[tpm_columns][max_col_index]
  
  tau_results$Max_Tissue <- max_tissue
  
  tau_values <- apply(normalized_table, 1, function(row) {
    sum(1 - row) / 4
  })
  
  tau_results$Tau <- tau_values
  
  return(tau_results)
}

tau_results_coll_chr <- calculate_tau(table_albicollis)
tau_results_coll_chr<-na.omit(tau_results_coll_chr)
mean_by_chromosome_coll <- aggregate(Tau ~ chr, data = tau_results_coll_chr, sum)
mean_by_chromosome_coll <- aggregate(Tau ~ chr, data = tau_results_coll_chr, mean)

tau_results_coll_chr$chr <- factor(tau_results_coll_chr$chr, levels = mixedsort(unique(tau_results_coll_chr$chr)))

 ggplot(tau_results_coll_chr, aes(x = chr, y = Tau)) +
  geom_point() +
  labs(x = "Chromosome", y = "Mean Tau Index", title = " Distribution Mean Tau Index by Chromosome") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Tourner les étiquettes de l'axe x

size<-size %>% rename( chr = chromosome)
mean_by_chromosome_coll<-merge(mean_by_chromosome_coll, size, by = "chr", all = T)
mean_by_chromosome_coll <- na.omit(mean_by_chromosome_coll)

mean_by_chromosome_coll$chr <- factor(mean_by_chromosome_coll$chr, levels = mixedsort(unique(mean_by_chromosome_coll$chr)))
mean_by_chromosome_coll <- mean_by_chromosome_coll[order(-mean_by_chromosome_coll$length), ]
mean_by_chromosome_coll$chr <- factor(mean_by_chromosome_coll$chr, levels = mean_by_chromosome_coll$chr)

plot <- ggplot(mean_by_chromosome_coll, aes(x = chr, y = Tau)) +
  geom_point() +
  labs(x = "Chromosome", y = "Mean Tau Index", title = "Mean Tau Index by Chromosome") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Tourner les étiquettes de l'axe x
print(plot)


#hypoleuca chromosome

calculate_tau <- function(table_hypoleuca) {
  
  tau_results <- data.frame(chr = table_hypoleuca$chr)
  
  tpm_columns <- grep("^TPM_\\w+_pied_mean", names(table_hypoleuca))
  
  num_tissues <- length(tpm_columns)
  
  max_values <- apply(table_hypoleuca[tpm_columns], 1, max)
  normalized_table <- table_hypoleuca[tpm_columns] / max_values
  
  max_col_index <- max.col(normalized_table)
  
  max_tissue <- names(table_hypoleuca)[tpm_columns][max_col_index]
  
  tau_results$Max_Tissue <- max_tissue
  
  tau_values <- apply(normalized_table, 1, function(row) {
    sum(1 - row) / 4
  })
  
  tau_results$Tau <- tau_values
  
  return(tau_results)
}

tau_results_pied_chr <- calculate_tau(table_hypoleuca)
tau_results_pied_chr<-na.omit(tau_results_pied_chr)
mean_by_chromosome_pied <- aggregate(Tau ~ chr, data = tau_results_pied_chr, sum)
mean_by_chromosome_pied<- aggregate(Tau ~ chr, data = tau_results_pied_chr, mean)

tau_results_pied_chr$chr <- factor(tau_results_pied_chr$chr, levels = mixedsort(unique(tau_results_pied_chr$chr)))

plot <- ggplot(tau_results_pied_chr, aes(x = chr, y = Tau)) +
  geom_point() +
  labs(x = "Chromosome", y = "Mean Tau Index", title = " Distribution Mean Tau Index by Chromosome") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Tourner les étiquettes de l'axe x
print(plot)

mean_by_chromosome_pied<-merge(mean_by_chromosome_pied, size, by = "chr", all = T)
mean_by_chromosome_pied <- na.omit(mean_by_chromosome_pied)
mean_by_chromosome_pied$chr <- factor(mean_by_chromosome_pied$chr, levels = mixedsort(unique(mean_by_chromosome_pied$chr)))
mean_by_chromosome_pied <- mean_by_chromosome_pied[order(-mean_by_chromosome_pied$length), ]
mean_by_chromosome_pied$chr <- factor(mean_by_chromosome_pied$chr, levels = mean_by_chromosome_pied$chr)

plot <- ggplot(mean_by_chromosome_pied, aes(x = chr, y = Tau)) +
  geom_point() +
  labs(x = "Chromosome", y = "Mean Tau Index", title = "Mean Tau Index by Chromosome") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Tourner les étiquettes de l'axe x
print(plot)

#hybrid chromosome

  calculate_tau <- function(table_hybrid) {
  
  tau_results <- data.frame(chr = table_hybrid$chr)
  
  tpm_columns <- grep("^TPM_\\w+_hybrid_mean", names(table_hybrid))
  
  num_tissues <- length(tpm_columns)
  
  max_values <- apply(table_hybrid[tpm_columns], 1, max)
  normalized_table <- table_hybrid[tpm_columns] / max_values
  
  max_col_index <- max.col(normalized_table)
  
  max_tissue <- names(table_hybrid)[tpm_columns][max_col_index]
  
  tau_results$Max_Tissue <- max_tissue
  
  tau_values <- apply(normalized_table, 1, function(row) {
    sum(1 - row) / 4
  })
  
  tau_results$Tau <- tau_values
  
  return(tau_results)
}

tau_results_hybrid_chr <- calculate_tau(table_hybrid)
tau_results_hybrid_chr<-na.omit(tau_results_hybrid_chr)
mean_by_chromosome_hybrid <- aggregate(Tau ~ chr, data = tau_results_hybrid_chr, sum)
mean_by_chromosome_hybrid<- aggregate(Tau ~ chr, data = tau_results_hybrid_chr, mean)

tau_results_hybrid_chr$chr <- factor(tau_results_hybrid_chr$chr, levels = mixedsort(unique(tau_results_hybrid_chr$chr)))

plot <- ggplot(tau_results_hybrid_chr, aes(x = chr, y = Tau)) +
  geom_point() +
  labs(x = "Chromosome", y = "Mean Tau Index", title = " Distribution Mean Tau Index by Chromosome") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Tourner les étiquettes de l'axe x
print(plot)

mean_by_chromosome_hybrid<-merge(mean_by_chromosome_hybrid, size, by = "chr", all = T)
mean_by_chromosome_hybrid <- na.omit(mean_by_chromosome_hybrid)
mean_by_chromosome_hybrid$chr <- factor(mean_by_chromosome_hybrid$chr, levels = mixedsort(unique(mean_by_chromosome_hybrid$chr)))
mean_by_chromosome_hybrid <- mean_by_chromosome_hybrid[order(-mean_by_chromosome_hybrid$length), ]
mean_by_chromosome_hybrid$chr <- factor(mean_by_chromosome_hybrid$chr, levels = mean_by_chromosome_hybrid$chr)

plot <- ggplot(mean_by_chromosome_pied, aes(x = chr, y = Tau)) +
  geom_point() +
  labs(x = "Chromosome", y = "Mean Tau Index", title = "Mean Tau Index by Chromosome") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Tourner les étiquettes de l'axe x
print(plot)

#######################################################################################################################################################

#correlation gene 

tau_results<-merge(tau_results_coll_gene, tau_results_pied_gene, by ="Gene", all =T)
tau_results <- tau_results %>%
  rename_with(~ gsub("\\.x$", "_coll", .x), ends_with(".x"))
tau_results <- tau_results %>%
  rename_with(~ gsub("\\.y$", "_pied", .x), ends_with(".y"))
tau_results_clean<- tau_results[complete.cases(tau_results$Tau_coll, tau_results$Tau_pied), ]

correlation_df <- data.frame(
  Gene = tau_results$Gene,  # Noms des gènes
  Correlation = sapply(tau_results$Gene, function(gene) {
    cor(tau_results$Tau_coll[tau_results$Gene == gene], tau_results$Tau_pied[tau_results$Gene == gene])
  })  # Calcul de la corrélation pour chaque gène
)

correlation_df <- tau_results %>%
  group_by(Gene) %>%
  summarise(Correlation = cor(Tau_coll, Tau_pied))

cor(tau_results$Tau_coll, tau_results$Tau_pied)
cor.test(tau_results$Tau_coll, tau_results$Tau_pied)

ggplot(tau_results_clean, aes(x = Tau_coll, y = Tau_pied)) +
  geom_point(aes(color = ifelse(Tau_coll > 0.85 & Tau_pied > 0.85, "Both High", 
                                ifelse(Tau_coll > 0.85, "Tau_coll High",
                                       ifelse(Tau_pied > 0.85, "Tau_pied High", "None"))))) +
  scale_color_manual(name = "",
                     values = c("Both High" = "purple", "Tau_coll High" = "red", "Tau_pied High" = "blue", "none" = "grey"),
                     labels = c("Both specific expression", "F. albicollis specific expression", "F. hypoleuca specific expression", "Wide-spread expression")) +
  labs(x = "tau in 5 Ficedula albicollis tissues", y = "tau in 5 Ficedula hypoleuca tissues", title = "Scatter Plot of gene specific expression between F. albicollis et F. hypoleuca (n=14607 genes)") +
  theme_minimal()
 
coll_0.85 <- sum(tau_results_clean$Tau_coll > 0.85)
pied_0.85 <- sum(tau_results_clean$Tau_pied > 0.85)
pied_less_0.85 <- sum(tau_results_clean$Tau_pied < 0.85)
coll_less_0.85 <- sum(tau_results_clean$Tau_coll < 0.85)
both_above_0.85 <- sum(tau_results_clean$Tau_coll > 0.85 & tau_results_clean$Tau_pied > 0.85)

#####################################################################################################################################################################################
# chromosome expression specific

count_0.85_coll <- tau_results_coll_chr %>%
  group_by(chr) %>%
  summarize(Count_Above_0.85 = sum(Tau > 0.85),
            Count_Below_0.85 = sum(Tau <= 0.85),
            Ratio = (Count_Above_0.85 / (Count_Above_0.85 + Count_Below_0.85)) * 100)

count_0.85_pied <- tau_results_pied_chr %>%
  group_by(chr) %>%
  summarize(Count_Above_0.85 = sum(Tau > 0.85),
            Count_Below_0.85 = sum(Tau <= 0.85),
            Ratio = (Count_Above_0.85 / (Count_Above_0.85 + Count_Below_0.85)) * 100)


#############################################################################################################################################################################################
# boxplot

TPM_coll <- select(table_albicollis, gene, TPM_Brain_coll_mean, TPM_Heart_coll_mean, TPM_Kidney_coll_mean, TPM_Liver_coll_mean, TPM_Testis_coll_mean) 
TPM_pied <- select(table_hypoleuca, gene, TPM_Brain_pied_mean, TPM_Heart_pied_mean, TPM_Kidney_pied_mean, TPM_Liver_pied_mean, TPM_Testis_pied_mean)
TPM<- merge(TPM_coll, TPM_pied, by = "gene", all = T)
TPM<-subset(TPM, select= -gene)

pied<-subset(tau_results_pied_gene, Tau >= 0.85 & Tau <= 1)

boxplot(tau_results,
        main = "Tissue specificity (tau) distribution for each species",
        names = c("albicollis", "hypoleuca"),  # Noms des boîtes
        ylab = "tau",
        col = c("skyblue", "red"))


#albicollis
filtered_data_coll<- subset(tau_results_coll_gene, Tau >= 0.85 & Tau <= 1)
ggplot(tau_results_coll_gene, aes(x = Max_Tissue, y = Tau, fill = Max_Tissue)) +
  geom_violin(color = "black", alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values = c("TPM_Testis_coll_mean" = "lightblue",
                               "TPM_Brain_coll_mean" = "lightpink",
                               "TPM_Heart_coll_mean" = "tomato",
                               "TPM_Liver_coll_mean" = "lightyellow",
                               "TPM_Kidney_coll_mean" = "lightgreen"),
                    guide = FALSE) +
  scale_x_discrete(labels = c("TPM_Testis_coll_mean" = "Testis (n=4799)",
                              "TPM_Brain_coll_mean" = "Brain (n=5321)",
                              "TPM_Heart_coll_mean" = "Heart (n=1858)",
                              "TPM_Liver_coll_mean" = "Liver (n=1143)",
                              "TPM_Kidney_coll_mean" = "Kidney (n=1734)")) +
  labs(title = "Tau for each tissue in F. albicollis  (n= 14 855 genes)",
       x = "Tissue",
       y = "Tau")

#hypoleuca
filtered_data_pied<- subset(tau_results_pied_gene, Tau >= 0.85 & Tau <= 1)
ggplot(filtered_data_pied, aes(x = Max_Tissue, y = Tau, fill = Max_Tissue)) +
  geom_violin(color = "black", alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values = c("TPM_Testis_pied_mean" = "lightblue",
                               "TPM_Brain_pied_mean" = "lightpink",
                               "TPM_Heart_pied_mean" = "tomato",
                               "TPM_Liver_pied_mean" = "lightyellow",
                               "TPM_Kidney_pied_mean" = "lightgreen"),
                    guide = FALSE) +
  scale_x_discrete(labels = c("TPM_Testis_pied_mean" = "Testis (n=4995)",
                              "TPM_Brain_pied_mean" = "Brain (n=5438)",
                              "TPM_Heart_pied_mean" = "Heart (n=2026)",
                              "TPM_Liver_pied_mean" = "Liver (n=1125)",
                              "TPM_Kidney_pied_mean" = "Kidney (n=1095)")) +
  labs(title = "Tau for each tissue in F. hypoleuca (n=14679 genes)",
       x = "Tissue",
       y = "Tau")

#hybrid
filtered_data_hybrid<- subset(tau_results_hybrid_gene, Tau >= 0.85 & Tau <= 1)
ggplot(filtered_data_hybrid, aes(x = Max_Tissue, y = Tau, fill = Max_Tissue)) +
  geom_violin(color = "black", alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values = c("TPM_Testis_hybrid_mean" = "lightblue",
                               "TPM_Brain_hybrid_mean" = "lightpink",
                               "TPM_Heart_hybrid_mean" = "tomato",
                               "TPM_Liver_hybrid_mean" = "lightyellow",
                               "TPM_Kidney_hybrid_mean" = "lightgreen"),
                    guide = FALSE) +
  scale_x_discrete(labels = c("TPM_Testis_hybrid_mean" = "Testis (n=)",
                              "TPM_Brain_hybrid_mean" = "Brain (n=)",
                              "TPM_Heart_hybrid_mean" = "Heart (n=)",
                              "TPM_Liver_hybrid_mean" = "Liver (n=)",
                              "TPM_Kidney_hybrid_mean" = "Kidney (n=)")) +
  labs(title = "Significant specific expression for each tissue in hybrids  (n= 14618 genes)",
       x = "Tissue",
       y = "Tau")

#########################################################################################
tau_results_chr<-merge(mean_by_chromosome_coll, mean_by_chromosome_pied, by="chr", all = T)
tau_results_chr<-subset(tau_results_chr, select= -length)
tau_results_chr<-subset(tau_results_chr, select= -category)
tau_results_chr<-tau_results_chr %>% rename(Mean_coll = Tau.x)
tau_results_chr<-tau_results_chr %>% rename(Mean_taiga = Tau.y)

cor(mean_by_chromosome_coll$Tau, mean_by_chromosome_pied$Tau)

tau_results_gene<-merge(tau_results_pied_gene, tau_results_coll_gene, by="Gene", all =T)
tau_results_gene<-na.omit(tau_results_gene)
cor(tau_results_gene$Tau.x, tau_results_gene$Tau.y)

tau_results_gene<-merge(tau_results_hybrid_gene, tau_results_coll_gene, by ="Gene", all = T)
tau_results_gene<-na.omit(tau_results_gene)
cor(tau_results_gene$Tau.x, tau_results_gene$Tau.y)

tau_results_gene<-merge(tau_results_hybrid_gene, tau_results_pied_gene, by ="Gene", all = T)
tau_results_gene<-na.omit(tau_results_gene)
cor(tau_results_gene$Tau.x, tau_results_gene$Tau.y)

########################################################################################

tau_results$heart_coll <- ifelse(tau_results$Max_Tissue_coll == "TPM_Heart_coll_mean", tau_results$Tau_coll, NA)
tau_results$heart_pied <- ifelse(tau_results$Max_Tissue_pied == "TPM_Heart_pied_mean", tau_results$Tau_pied, NA)
tau_results$brain_coll <- ifelse(tau_results$Max_Tissue_coll == "TPM_Brain_coll_mean", tau_results$Tau_coll, NA)
tau_results$brain_pied <- ifelse(tau_results$Max_Tissue_pied == "TPM_Brain_pied_mean", tau_results$Tau_pied, NA)
tau_results$testis_coll <- ifelse(tau_results$Max_Tissue_coll == "TPM_Testis_coll_mean", tau_results$Tau_coll, NA)
tau_results$testis_pied <- ifelse(tau_results$Max_Tissue_pied == "TPM_Testis_pied_mean", tau_results$Tau_pied, NA)
tau_results$kidney_coll <- ifelse(tau_results$Max_Tissue_coll == "TPM_Kidney_coll_mean", tau_results$Tau_coll, NA)
tau_results$kidney_pied <- ifelse(tau_results$Max_Tissue_pied == "TPM_Kidney_pied_mean", tau_results$Tau_pied, NA)
tau_results$liver_coll <- ifelse(tau_results$Max_Tissue_coll == "TPM_Liver_coll_mean", tau_results$Tau_coll, NA)
tau_results$liver_pied <- ifelse(tau_results$Max_Tissue_pied == "TPM_Liver_pied_mean", tau_results$Tau_pied, NA)

cor(tau_results$heart_coll, tau_results$heart_pied, use = "pairwise.complete.obs")
cor(tau_results$brain_coll, tau_results$brain_pied, use = "pairwise.complete.obs")
cor(tau_results$testis_coll, tau_results$testis_pied, use = "pairwise.complete.obs")
cor(tau_results$kidney_coll, tau_results$kidney_pied, use = "pairwise.complete.obs")
cor(tau_results$liver_coll, tau_results$liver_pied, use = "pairwise.complete.obs")

matrice_correlation <- cor(tau_results[,6:15], use = "pairwise.complete.obs")
corrplot(matrice_correlation, method="circle", type="upper", tl.col="black", tl.srt=45)

tau_results$heart_coll <- ifelse(tau_results$Max_Tissue_coll == "TPM_Heart_coll_mean" & tau_results$Tau_coll >= 0.85, tau_results$Tau_coll, NA)
tau_results$heart_pied <- ifelse(tau_results$Max_Tissue_pied == "TPM_Heart_pied_mean" & tau_results$Tau_pied >= 0.85, tau_results$Tau_pied, NA)
tau_results$brain_coll <- ifelse(tau_results$Max_Tissue_coll == "TPM_Brain_coll_mean" & tau_results$Tau_coll >= 0.85, tau_results$Tau_coll, NA)
tau_results$brain_pied <- ifelse(tau_results$Max_Tissue_pied == "TPM_Brain_pied_mean" & tau_results$Tau_pied >= 0.85, tau_results$Tau_pied, NA)
tau_results$testis_coll <- ifelse(tau_results$Max_Tissue_coll == "TPM_Testis_coll_mean" & tau_results$Tau_coll >= 0.85, tau_results$Tau_coll, NA)
tau_results$testis_pied <- ifelse(tau_results$Max_Tissue_pied == "TPM_Testis_pied_mean" & tau_results$Tau_pied >= 0.85, tau_results$Tau_pied, NA)
tau_results$kidney_coll <- ifelse(tau_results$Max_Tissue_coll == "TPM_Kidney_coll_mean" & tau_results$Tau_coll >= 0.85, tau_results$Tau_coll, NA)
tau_results$kidney_pied <- ifelse(tau_results$Max_Tissue_pied == "TPM_Kidney_pied_mean" & tau_results$Tau_pied >= 0.85, tau_results$Tau_pied, NA)
tau_results$liver_coll <- ifelse(tau_results$Max_Tissue_coll == "TPM_Liver_coll_mean" & tau_results$Tau_coll >= 0.85, tau_results$Tau_coll, NA)
tau_results$liver_pied <- ifelse(tau_results$Max_Tissue_pied == "TPM_Liver_pied_mean" & tau_results$Tau_pied >= 0.85, tau_results$Tau_pied, NA)

matrice_correlation_significant <- cor(tau_results[,6:13], use = "pairwise.complete.obs")
corrplot(matrice_correlation_significant)

#############################################################################################"

# TAU by chr

size<-size %>% rename(chr = chromosome)

# croissaant
#tau_results_coll_chr$chr <- reorder(tau_results_coll_chr$chr, tau_results_coll_chr$Tau, FUN = median)
#tau_results_pied_chr$chr <- reorder(tau_results_pied_chr$chr, tau_results_pied_chr$Tau, FUN = median)
#tau_results_hybrid_chr$chr <- reorder(tau_results_hybrid_chr$chr, tau_results_hybrid_chr$Tau, FUN = median)

# significant gene 
#tau_results_coll_chr<- subset(tau_results_coll_chr, Tau >= 0.85 & Tau <= 1)
#tau_results_pied_chr<- subset(tau_results_pied_chr, Tau >= 0.85 & Tau <= 1)
#tau_results_hybrid_chr<- subset(tau_results_hybrid_chr, Tau >= 0.85 & Tau <= 1)

tau_results_coll_chr<-merge(tau_results_coll_chr, size, by = "chr", all = T)
tau_results_coll_chr <- subset(tau_results_coll_chr, chr != "ChrFal34")
tau_results_coll_chr <- subset(tau_results_coll_chr, chr != "ChrFal36")
ggplot(tau_results_coll_chr, aes(x = chr, y = Tau, fill = category)) +
  geom_boxplot() +
  xlab("Chromosome") +
  ylab("Tau") +
  ggtitle("Significant tau for each chromosome (F.albicollis)") +
  scale_fill_discrete(name = "Category") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


tau_results_pied_chr<-merge(tau_results_pied_chr, size, by = "chr", all = T)
tau_results_pied_chr <- subset(tau_results_pied_chr, chr != "ChrFal34")
tau_results_pied_chr <- subset(tau_results_pied_chr, chr != "ChrFal36")
ggplot(tau_results_pied_chr, aes(x = chr, y = Tau, fill = category)) +
  geom_boxplot() +
  xlab("Chromosome") +
  ylab("Tau") +
  ggtitle("Signficant tau for each chromosome (F.hypoleuca)") +
  scale_fill_discrete(name = "Category") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 


tau_results_hybrid_chr<-merge(tau_results_hybrid_chr, size, by = "chr", all = T)
tau_results_hybrid_chr <- subset(tau_results_hybrid_chr, chr != "ChrFal34")
tau_results_hybrid_chr <- subset(tau_results_hybrid_chr, chr != "ChrFal36")
ggplot(tau_results_hybrid_chr, aes(x = chr, y = Tau, fill = category)) +
  geom_boxplot() +
  xlab("Chromosome ordered by length") +
  ylab("Tau") +
  ggtitle("Significant tau for each chromosome (Hybrid)") +
  scale_fill_discrete(name = "Category") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 


category_order <- c("Macro", "Intermediate", "Micro")
tau_results_coll_chr$Max_Tissue <- gsub("TPM_Kidney_coll_mean", "Kidney", tau_results_coll_chr$Max_Tissue)
tau_results_coll_chr$Max_Tissue <- gsub("TPM_Liver_coll_mean", "Liver", tau_results_coll_chr$Max_Tissue)
tau_results_coll_chr$Max_Tissue <- gsub("TPM_Heart_coll_mean", "Heart", tau_results_coll_chr$Max_Tissue)
tau_results_coll_chr$Max_Tissue <- gsub("TPM_Brain_coll_mean", "Brain", tau_results_coll_chr$Max_Tissue)
tau_results_coll_chr$Max_Tissue <- gsub("TPM_Testis_coll_mean", "Testis", tau_results_coll_chr$Max_Tissue)
ggplot(tau_results_coll_chr, aes(x = Max_Tissue, y = Tau, fill = category)) +
  geom_boxplot() +
  xlab("Tissue") +  
  ylab("Tau") +
  ggtitle("Significant tau for each tissue (F.albicollis)") +
  scale_fill_manual(name = "Category", values = c("Macro" = "blue", "Intermediate" = "green", "Micro" = "red"), labels = category_order) +
  theme() 

category_order <- c("Macro", "Intermediate", "Micro")
tau_results_pied_chr$Max_Tissue <- gsub("TPM_Kidney_pied_mean", "Kidney", tau_results_pied_chr$Max_Tissue)
tau_results_pied_chr$Max_Tissue <- gsub("TPM_Liver_pied_mean", "Liver", tau_results_pied_chr$Max_Tissue)
tau_results_pied_chr$Max_Tissue <- gsub("TPM_Heart_pied_mean", "Heart", tau_results_pied_chr$Max_Tissue)
tau_results_pied_chr$Max_Tissue <- gsub("TPM_Brain_pied_mean", "Brain", tau_results_pied_chr$Max_Tissue)
tau_results_pied_chr$Max_Tissue <- gsub("TPM_Testis_pied_mean", "Testis", tau_results_pied_chr$Max_Tissue)
ggplot(tau_results_pied_chr, aes(x = Max_Tissue, y = Tau, fill = category)) +
  geom_boxplot() +
  xlab("Tissue") +  
  ylab("Tau") +
  ggtitle("Tau for each tissue (F.hypoleuca)") +
  scale_fill_manual(name = "Category", values = c("Macro" = "blue", "Intermediate" = "green", "Micro" = "red"), labels = category_order) +
  theme() 

category_order <- c("Macro", "Intermediate", "Micro")
tau_results_hybrid_chr$Max_Tissue <- gsub("TPM_Kidney_hybrid_mean", "Kidney", tau_results_hybrid_chr$Max_Tissue)
tau_results_hybrid_chr$Max_Tissue <- gsub("TPM_Liver_hybrid_mean", "Liver", tau_results_hybrid_chr$Max_Tissue)
tau_results_hybrid_chr$Max_Tissue <- gsub("TPM_Heart_hybrid_mean", "Heart", tau_results_hybrid_chr$Max_Tissue)
tau_results_hybrid_chr$Max_Tissue <- gsub("TPM_Brain_hybrid_mean", "Brain", tau_results_hybrid_chr$Max_Tissue)
tau_results_hybrid_chr$Max_Tissue <- gsub("TPM_Testis_hybrid_mean", "Testis", tau_results_hybrid_chr$Max_Tissue)
ggplot(tau_results_hybrid_chr, aes(x = Max_Tissue, y = Tau, fill = category)) +
  geom_boxplot() +
  xlab("Tissue") +  
  ylab("Tau") +
  ggtitle("Tau for each tissue (Hybrid)") +
  scale_fill_manual(name = "Category", values = c("Macro" = "blue", "Intermediate" = "green", "Micro" = "red"), labels = category_order) +
  theme() 

