library(ggplot2)
library(extrafont)

###Pausing index (PI) calculation###
###Importing data from bwtool into R###

den_26_rep1 <- read.table("den_26x_rep1_new.txt", header = FALSE)
num_26_rep1 <- read.table("num_26x_rep1_new.txt", header = FALSE)

den_26_rep2 <- read.table("den_26x_rep2_new.txt", header = FALSE)
num_26_rep2 <- read.table("num_26x_rep2_new.txt", header = FALSE)

den_52_rep1 <- read.table("den_52x_rep1_new.txt", header = FALSE)
num_52_rep1 <- read.table("num_52x_rep1_new.txt", header = FALSE)

den_52_rep2 <- read.table("den_52x_rep2_new.txt", header = FALSE)
num_52_rep2 <- read.table("num_52x_rep2_new.txt", header = FALSE)

gene_names <- read.table("genes.txt", header = FALSE)

###Combining data into single dataframe and using PI formula###
Pausing_df <- data.frame(Genes = gene_names$V1, 
                         CTD_26_rep2 = (num_26_rep2$V10/350)/(den_26_rep2$V10/den_26_rep2$V4), 
                         CTD_52_rep2 = (num_52_rep2$V10/350)/(den_52_rep2$V10/den_52_rep2$V4))

###cleaning###
good <- complete.cases(Pausing_df)
Pausing_df_clean <- Pausing_df[good, ]
table(is.na(Pausing_df_clean))
dim(Pausing_df_clean)
table(duplicated(Pausing_df_clean$Genes))

Pausing_df_clean2 <- Pausing_df_clean[!duplicated(Pausing_df_clean[, "Genes"]), ]
rownames(Pausing_df_clean2) <- Pausing_df_clean2$Genes
Pausing_df_clean2 <- Pausing_df_clean2[, -1]

###Stats for all the genes###
boxplot(Pausing_df_clean2$CTD_26_rep2, Pausing_df_clean2$CTD_52_rep2,
         outline = FALSE,
         xlab = "Samples",
         ylab = "Pausing Index",
         names = c("CTD_26x", "CTD_52x"),
         main="Protein-coding genes",
         col = c("orange", "aquamarine"))
 
wilcox.test(Pausing_df_clean2$CTD_26_rep2, Pausing_df_clean2$CTD_52_rep2, paired = TRUE, alternative = "two.sided")
median(Pausing_df_clean2$CTD_26_rep2)
median(Pausing_df_clean2$CTD_52_rep2)
###How many genes are paused /PI >2?###
sum(Pausing_df_clean2[, 1] >= 2) 
sum(Pausing_df_clean2[, 2] >= 2) 

###Clustering genes based on their lengths###
Pausing_df_clean2$means <- rowMeans(Pausing_df_clean2)
PI_genes_ranked <- Pausing_df_clean2[order(Pausing_df_clean2$means),]

###remove Inf values###
PI_genes_ranked <- PI_genes_ranked[!is.infinite(PI_genes_ranked[, "means"]), ]

###G0,G1,G2,G3 clusters where G0 genes with PI =0, and G3 are most paused###
G0 <- subset(PI_genes_ranked, means == 0)
dim(G0) 

above_0 <- subset(PI_genes_ranked, means != 0)
quantiles <- quantile(above_0$means, probs = seq(0, 1, 1/4))
quantiles

G1 <- subset(above_0, means < quantiles[2]) 
dim(G1)
G2 <- subset(above_0, means >= quantiles[2] & means < quantiles[4]) 
dim(G2)
G3 <- subset(above_0, means >= quantiles[4]) 
dim(G3)

###plotting###
x <-  c(nrow(G0), nrow(G1), nrow(G2), nrow(G3))
labels <-  c("G0","G1","G2","G3")
png(file = "genes_cluster2.jpg")
pie(x, labels = x, main = "Gene clusters based on PI",col = rainbow(length(x)))
legend("topright", c("G0","G1","G2","G3"), cex = 1,
       fill = rainbow(length(x)))
dev.off()

###Comparing PI in cluster G3###
boxplot(G3$CTD_26_rep2, G3$CTD_52_rep2,
        outline = FALSE,
        xlab = "Samples",
        ylab = "Pausing Index",
        names = c("CTD_26_rep2", "CTD_52_rep2"),
        main="G3",
        col = c("orange", "aquamarine"))

wilcox.test(G3$CTD_26_rep2, G3$CTD_52_rep2, paired = TRUE, alternative = "two.sided")
median(G3$CTD_26_rep2)
median(G3$CTD_52_rep2)
