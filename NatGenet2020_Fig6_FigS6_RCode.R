### Finalised code for PBM-based RBC SNPs analysis, Nature Genetics, 2020
### Kian Hong Kock (Martha Bulyk Lab), kianhongkock@g.harvard.edu, mlbulyk@genetics.med.harvard.edu

### Set working directory as the folder containing this piece of code and the required files
### Read in original PBM dataset files

setwd("PBM_datasets/")
list_files_PBM <- list.files()
list_files_PBM <- list_files_PBM[-1]

for (i in 1:length(list_files_PBM)){
  name_PBM_file <- strsplit(list_files_PBM[i], split = ".txt")[[1]]
  assign(name_PBM_file, read.table(file = list_files_PBM[i], header = FALSE))
  if (dim(get(name_PBM_file))[1] == 32897){
    assign(name_PBM_file, read.table(file = list_files_PBM[i], header = TRUE))
  }
}

### All the indices should be the same for each row in the PBM dataset files - confirm this is the case. 
### If so, vector of indices can be precomputed, and we can populate the new PBM E-score files with
### the master vector of indices.

for (i in 1:length(list_files_PBM)){
  name_PBM_file <- strsplit(list_files_PBM[i], split = ".txt")[[1]]
  print(setdiff(KLF1_REF_R2_8mers$X8.mer, get(name_PBM_file)[,1])  )
}
### Output for this was all character(0) - indicates that all the 8-mers were in the same order in all PBM files

for (i in 1:length(list_files_PBM)){
  name_PBM_file <- strsplit(list_files_PBM[i], split = ".txt")[[1]]
  print(name_PBM_file)  
  df <- get(name_PBM_file)
  for (j in 1:2){
    print(paste("Column ", as.character(j), sep = ""))
    print((df[c(1:5), j]))
  }
}

### Output for this was all 8-mers - indicates that first two columns of all PBM files were all 8-mers

for (i in 1:length(list_files_PBM)){
  name_PBM_file <- strsplit(list_files_PBM[i], split = ".txt")[[1]]
  df <- get(name_PBM_file)
  if (dim(df)[2] == 20){
    print(name_PBM_file)  
    print(range(df[, 4]))
    print(range(df[, 7]))
  }
}   

### For PBM files with 5 columns, Column #3 is the E-score (checked using controlloop + print for dim(df)[2] == 5)
### For PBM files with 4 columns, Column #4 is the E-score (checked $V4)
### For PBM files with 7 columns, Column #3 is the E-score (checked colnames)
### For PBM files with 20 columns, Columns #4 and #7 are the E-scores (checked using controlloop + print for dim(df)[2] == 20)
### Use column #4 (AMADID #015681) as E-score for PBM files with 20 columns

### Obtain only E-score columns for PBM files

for (i in 1:length(list_files_PBM)){
  name_PBM_file <- strsplit(list_files_PBM[i], split = ".txt")[[1]]
  number_of_columns <- dim(get(name_PBM_file))[2]
  PBM_df <- get(name_PBM_file)
  if (number_of_columns == 5){
    Escore_column <- 3
  }
  if (number_of_columns == 4){
    Escore_column <- 4
  }
  if (number_of_columns == 7){
    Escore_column <- 3
  }
  if (number_of_columns == 20){
    Escore_column <- 4
  }
  vec_eightmers <- c(as.character(PBM_df[,1]), as.character(PBM_df[,2]))
  vec_Escore <- c(PBM_df[,Escore_column], PBM_df[,Escore_column])
  newdf <- cbind(vec_eightmers, vec_Escore)
  write.table(newdf, file = paste(name_PBM_file, "Escore.txt", sep = "_"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

### Getting TCF4 and HES1/5/7 PBMs from CIS-BP mouse PBM file - note that these only contain one 8-mer orientation;  
### the accompanying Jupyter notebook fixes this by assigning the reverse complement the same E-score.

cisbp_mouse <- read.table(file = "CIS-BP/musmusculus_PBM_cisbp_Escores.txt", header = TRUE)

grep("M0191_1.02", colnames(cisbp_mouse))
# [1] 605 606
PBM_HES1 <- cisbp_mouse[, c(1, 605, 606)]
write.table(PBM_HES1[, c(1:2)], file = "HES1_R1_Escore.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

grep("M0215_1.02", colnames(cisbp_mouse))
# [1] 671 672
PBM_HES5 <- cisbp_mouse[, c(1, 671, 672)]
write.table(PBM_HES5[, c(1:2)], file = "HES5_R1_Escore.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

grep("M0192_1.02", colnames(cisbp_mouse))
# [1] 663 664
PBM_HES7 <- cisbp_mouse[, c(1, 663, 664)]
write.table(PBM_HES7[, c(1:2)], file = "HES7_R1_Escore.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

grep("M0219_1.02", colnames(cisbp_mouse))
# [1] 637 638
PBM_TCF4 <- cisbp_mouse[, c(1, 637, 638)]
write.table(PBM_TCF4[, c(1:2)], file = "TCF4_R1_Escore.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

### Identifying perturbed TF binding, and bootstrapping for 3263 RBC SNPs using 15 TFs + GATA average
### Criteria: E-score >0.35 for one allele, E-score < 0.3 for other allele)

setwd("..")
setwd("empirical_background/")

vec_background_files <- list.files()

### Reading in background distribution for each PBM dataset

for (i in 1:length(vec_background_files)){
  name <- strsplit(vec_background_files[i], split = ".txt")[[1]][1]
  vec_background_from_file <- as.numeric(readLines(con = vec_background_files[i]))
  assign(name, vec_background_from_file)
}

setwd("..")

foreground_3263SNPs_035030 <- read.table(file = "foreground_3263SNPs.txt", header = FALSE, sep = ",")

SigCtr_vec_boostrap_pval <- c()
SigCtr_vec_bootstrap_mean <- c()
SigCtr_vec_obs_over_exp <- c()
# eightmer_SigCtr_vec_boostrap_pval <- c()
# eightmer_SigCtr_vec_bootstrap_mean <- c()
# eightmer_SigCtr_vec_obs_over_exp <- c()

### 100,000 samples from each background (of all common SNPs > 10% allele frequency in dbGaP)

for (module in 1:length(foreground_3263SNPs_035030$V1)){
  char_module <- as.character(foreground_3263SNPs_035030[module, 1])
  vec_background_to_sample <- get(char_module)
  vec_test_sample <- c()
#  vec_cases_sample <- c()
  for (i in 1:100000){
    int_test_sample <- sample(vec_background_to_sample, size = 3263, replace = TRUE)
    vec_test_sample <- c(vec_test_sample, sum(int_test_sample))
    # vec_cases_sample <- c(vec_cases_sample, sum(int_test_sample>0))
  }
  observed <- foreground_3263SNPs_035030[module, 3] + foreground_3263SNPs_035030[module, 4]  
  SigCtr_vec_boostrap_pval[module] <- length(which(vec_test_sample > observed))/length(vec_test_sample)
  SigCtr_vec_bootstrap_mean[module] <- mean(vec_test_sample)
  SigCtr_vec_obs_over_exp[module] <- observed/mean(vec_test_sample)
  # observed_8mer <- foreground_3263SNPs_035030[module, 5]
  # eightmer_SigCtr_vec_boostrap_pval[module] <- length(which(vec_cases_sample > observed_8mer))/length(vec_cases_sample)
  # eightmer_SigCtr_vec_bootstrap_mean[module] <- mean(vec_cases_sample)
  # eightmer_SigCtr_vec_obs_over_exp[module] <- observed_8mer/mean(vec_cases_sample)
}

bootstrap_df_foreground_3263SNPs_035030 <- data.frame(list((foreground_3263SNPs_035030$V1), 
                                                           SigCtr_vec_bootstrap_mean,
                                                           SigCtr_vec_obs_over_exp,
                                                           SigCtr_vec_boostrap_pval))###,
                                                           # eightmer_SigCtr_vec_bootstrap_mean,
                                                           # eightmer_SigCtr_vec_obs_over_exp,
                                                           # eightmer_SigCtr_vec_boostrap_pval))

colnames(bootstrap_df_foreground_3263SNPs_035030) <- c("Individual_PBM", "Bootstrap_mean", "Obs_over_Exp", "Bootstrap_pval")###, "Eightmer_Bootstrap_mean", "Eightmer_Obs_over_Exp", "Eightmer_Bootstrap_pval")

### Table of observed/expected ratios and q-values to report

datatable_NatGenet2020_Fig6C <- bootstrap_df_foreground_3263SNPs_035030[,c(1,3,4)]

### Obtain Benjamini-Hochberg-corrected q-values for bootstrapping analysis

datatable_NatGenet2020_Fig6C$Bootstrap_qval <- p.adjust(datatable_NatGenet2020_Fig6C$Bootstrap_pval, method = "fdr")
  
write.table(datatable_NatGenet2020_Fig6C, file = "NatGenet2020_Fig6C_values.txt", quote = FALSE, sep = "\t", row.names = FALSE)

### Plotting barcharts for Figure 6C

library(ggplot2)
ordered_barchart_general <- bootstrap_df_foreground_3263SNPs_035030[order(-bootstrap_df_foreground_3263SNPs_035030$Obs_over_Exp),]

rownames <- c()
for (i in 1:length(ordered_barchart_general$Individual_PBM)){
  rownames <- c(rownames, strsplit(as.character(ordered_barchart_general$Individual_PBM[i]), split = "_")[[1]][1])
}
ordered_barchart_general$Individual_PBM <- rownames

ggplot(data=ordered_barchart_general, aes(x= reorder(Individual_PBM, -Obs_over_Exp), y=Obs_over_Exp)) +
  geom_bar(stat="identity") +
  xlab("Individual PBM datasets") + ylab("Ratio of Observed / Expected of binary cases") +
  ggtitle("TFs (from individual PBM datasets) with greater # of perturbed binding events for all 3263 RBC SNPs than expected (background: all common SNPs in dbSNP") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_flip()

ggsave("PerturbedTFbinding_PBM_3263.pdf", device = "pdf")

### Getting boxplots and barcharts of TF PBM E-scores-rsID pairs - set your working directory to one containing the data files of interest

list_input_files <- list.files()
library(ggplot2)

for (i in 1:length(list_input_files)){
  name <- strsplit(list_input_files[i], split = ".txt")[[1]][1]
  assign(name, read.table(file = list_input_files[i], header = FALSE))
  TF_name <- strsplit(name, split = "_")[[1]][1]
  ### May need to adjust the index for following command, depending on naming format
  rsID_name <- strsplit(name, split = "_")[[1]][4]
  boxplot_title <- paste(TF_name, " PBM E-scores for ", rsID_name, " 8-mers", sep = "")
  boxplot_filename <- paste("Boxplot_", name, ".pdf", sep = "")
  allele1_character <- strsplit(as.character(get(name)$V1[1]), split = "")[[1]][8]
  allele2_character <- strsplit(as.character(get(name)$V2[1]), split = "")[[1]][8]
  allele1_label <- paste(allele1_character, "allele", sep = " ")
  allele2_label <- paste(allele2_character, "allele", sep = " ")
  vec_allele_labels <- c(rep(allele1_label, 8), rep(allele2_label, 8))
  vec_allele_Escores <- c(get(name)$V3, get(name)$V4)
  df_boxplot <- data.frame(list(vec_allele_Escores, vec_allele_labels))
  colnames(df_boxplot) <- c("E-score", "Allele")
  p <- ggplot(df_boxplot, aes(Allele, `E-score`))
  p + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) + labs(title = boxplot_title) + theme(plot.title = element_text(hjust = 0.5))
  ggsave(boxplot_filename, device = "pdf", width = 15, units = "cm")
  
  vec_Escore_plot <- c()
  vec_index_plot <- c()
  vec_allele_plot <- c()
  for (j in 1:8){
    vec_Escore_plot <- c(vec_Escore_plot, as.numeric(get(name)[j, 3]), as.numeric(get(name)[j, 4]))
    vec_index_plot <- c(vec_index_plot, paste("Position", as.character(j), sep = "_"), paste("Position", as.character(j), sep = "_"))
    vec_allele_plot <- c(vec_allele_plot, allele1_character, allele2_character)
  }
  Escore_plot <- data.frame(vec_Escore_plot, vec_index_plot, vec_allele_plot)
  ggplot(Escore_plot, aes(fill = vec_allele_plot, y = vec_Escore_plot, x = vec_index_plot)) +
    geom_bar(position = "dodge", stat = "identity") +
    labs(title = name) + xlab("Position in 15-mer") + ylab("E-score") + 
    theme(plot.title = element_text(hjust = 0.5)) + labs(title = boxplot_title) + 
    geom_hline(yintercept = 0.35, linetype = "dashed", color = "red", size = 1) + guides(fill=guide_legend(title=NULL))
  barcharts_filename <- paste("Barcharts_", name, ".pdf", sep = "")
  ggsave(barcharts_filename, device = "pdf")
}

### Statistical tests for significance: Supplementary Data Table 3, Figure 6, and Supplementary Figure 6
### Fisher's exact test (two-sided) for Supplementary Data Tables 3A and 3B

supptable3a_PU1_Gata2only <- matrix(data = c(5079, 2424, (6208), (19640)), nrow = 2, ncol = 2)
fisher.test(supptable3a_PU1_Gata2only, alternative = "two.sided")

supptable3a_KLF1_Gata1only <- matrix(data = c(1355, 1019, (1362), (7118)), nrow = 2, ncol = 2)
fisher.test(supptable3a_KLF1_Gata1only, alternative = "two.sided")

supptable3b_Smad2only <- matrix(data = c(4549, 3881, (11189), (48940)), nrow = 2, ncol = 2)
fisher.test(supptable3b_Smad2only, alternative = "two.sided")

supptable3b_Tcf7l2only <- matrix(data = c(4549, 3881, (175), (3095)), nrow = 2, ncol = 2)
fisher.test(supptable3b_Tcf7l2only, alternative = "two.sided")

### Calculating Figures 6 and S6 Wilcoxon signed-rank test (two-sided) p-values 
### from comparisons of PBM E-scores for allele 1 versus allele 2

setwd(..)
setwd("rsID_Escores")
vec_pval <- c()
name_pval <- c()
list_PBM_files <- list.files()

for (i in 1:length(list_PBM_files)){
  df_pval_calc <- read.table(file = list_PBM_files[i], header = FALSE, sep = "\t")
  wilcox_pval <- wilcox.test(df_pval_calc$V3, df_pval_calc$V4, paired = TRUE)$p.val
  vec_pval <- c(vec_pval, wilcox_pval)
  name_pval <- c(name_pval, list_PBM_files[i])  
}

df_PBM_pval <- as.data.frame(cbind(name_pval, vec_pval))
colnames(df_PBM_pval) <- c("PBM-rsID combination", "Wilcoxon signed-rank test p-val")

setwd(..)
write.table(df_PBM_pval, file = "Fig6_S6_pval.txt", sep = "\t", quote = FALSE)

### Examining perturbed TF binding in fine-mapped SNPs from Ulirsch Nat Genet 2019

STF_PP0001 <- matrix(c(911, 646, 1000, 706), ncol = 2, nrow = 2)
STF_PP001 <- matrix(c(542, 1015, 522, 1184), ncol = 2, nrow = 2)
GATA_PP0001 <- matrix(c(276, 172, 1635, 1180), ncol = 2, nrow = 2)
GATA_PP001 <- matrix(c(158, 290, 906, 1909), ncol = 2, nrow = 2)
rownames(STF_PP0001) <- c("Fine-mapped", "Not in fine-mapped list")
colnames(STF_PP0001) <- c("Perturbed in STF binding", "Not perturbed in STF binding")
rownames(STF_PP001) <- c("Fine-mapped", "Not in fine-mapped list")
colnames(STF_PP001) <- c("Perturbed in STF binding", "Not perturbed in STF binding")

rownames(GATA_PP0001) <- c("Fine-mapped", "Not in fine-mapped list")
colnames(GATA_PP0001) <- c("Perturbed in GATA binding", "Not perturbed in GATA binding")
rownames(GATA_PP001) <- c("Fine-mapped", "Not in fine-mapped list")
colnames(GATA_PP001) <- c("Perturbed in GATA binding", "Not perturbed in GATA binding")

STF_PP001
STF_PP0001
GATA_PP001
GATA_PP0001

fisher.test(STF_PP001, alternative = "two.sided")
fisher.test(STF_PP0001, alternative = "two.sided")
fisher.test(GATA_PP001, alternative = "two.sided")
fisher.test(GATA_PP0001, alternative = "two.sided")
