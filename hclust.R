library(dplyr)
library(ggplot2)

#### RNAseq clustering
Hct_RNAseq_WTDelR_pc <- read.csv("./Hct116_RNAseq_WTDelR_pcgenes.csv")
Hct_RNAseq_WTDelR_pc_sig <- Hct_RNAseq_WTDelR_pc %>% filter(Del_vs_WT.is.DE == 1 | R_vs_WT.is.DE == 1 | Del_vs_WT.is.DE == -1 | R_vs_WT.is.DE == -1)

Hct_RNAseq_WTDelR_pcsigCPM <- Hct_RNAseq_WTDelR_pc_sig %>% select(1, 10:18)

write.csv(Hct_RNAseq_WTDelR_pcsigCPM, file = "./Hct_RNAseq_WTDelR_pcsigCPM.csv")
#Remove the first column on csv file and re-read
Hct_RNAseq_WTDelR_pcsigCPM_edit <- read.csv("./Hct_RNAseq_WTDelR_pcsigCPM.csv")

#write.csv(Hct_RNAseq_WTDelR, file="./Hct_RNAsew_WTDelR.csv")
#Hct_RNAseq_WTDelR_rmdup <- read.csv("./Hct_RNAsew_WTDelR_edited.csv")
#Hct_RNAseq_WTDelR_edit1 <- Hct_RNAseq_WTDelR_rmdup[, -(9:11)]
#Hct_RNAseq_WTDelR_edit <- Hct_RNAseq_WTDelR_edit1[, -1]

#Reading in data
mydata <- data.frame(Hct_RNAseq_WTDelR_pcsigCPM_edit)
#Convert dataframe to a matrix
mydata_matrix <- mydata %>% 
  dplyr::select(-Gene_name) %>%
  as.matrix()

#assign rownames
rownames(mydata_matrix) <- mydata$Gene_name

# Scale the data (to obtain z-scores). 
mydata_matrix <- mydata_matrix %>%
  # transpose the matrix so genes are as columns
  t() %>%
  # apply scalling to each column of the matrix (genes)
  scale() %>%
  # transpose back so genes are as rows again
  t()

hclustfunc <- function(x) hclust(x, method = "complete")
distfunc <- function(x) dist(x, method = "euclidean")
d <- distfunc(mydata_matrix)
fit <- hclustfunc(d)

#plot dendogram only
plot(fit)
groups <- cutree(fit, k=6)
write.csv(groups, file = "./Hct116_RNAseq_WTDelRpcsig_clusters.csv" )

#Add a rectangle to identify cluster
#rect.hclust(fit, k=4, border = "red")

#Plot heatmap with dendogram 

my_pallete <- colorRampPalette(c("blue", "white", "red"))(n=299)

heatmap.2(as.matrix(mydata_matrix), dendrogram ="both",trace="none",margin=c(8,9), hclust=hclustfunc, distfun = distfunc, RowSideColors =as.character(groups), col = my_pallete, 
          key = TRUE, density.info = "none", Colv = FALSE)

#Export heatmap as PDF to preserve resolution. 

dev.off()
