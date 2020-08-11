library(tidyverse)
require(data.table)
library("survival")
library("survminer")

DMPs<-read_tsv("/Users/fazal2/Desktop/tidy_practice/B3_Case1_Genes_Beta0.4_0.05FDR.txt")
head(DMPs)
clinical<-read_tsv("/Users/fazal2/Desktop/tidy_practice/tgct_tcga_clinical_data.tsv")
head(clinical)
colnames(clinical)<-c("ID","PID","Sample_ID","CancerType","DSF_Time","DFS_Status")
#reform_clin=reshape2::melt(clinical)
#colnames(reform_clin)<-c("ID","PID","Sample_ID","CancerType","DFS_Status","variable","Time")
#head(reform_clin)

mRNA<-read_tsv("/Users/fazal2/Desktop/tidy_practice/data_RNA_Seq_v2_expression_median.txt")
reform_mRNA=reshape2::melt(mRNA)
head(reform_mRNA)
colnames(reform_mRNA)<-c("Gene","SampleID","expression")
dim(reform_mRNA)

#create a list of the files from your target directory
  test=reform_mRNA %>% 
    filter(Gene %in% DMPs$Gene) %>% 
    group_by(SampleID) %>% 
    pivot_wider(names_from = SampleID, values_from = expression) %>% 
    drop_na()
  
  dim(test)

#calculate mean of each sample
Samples_Mean=test %>% 
  summarise(across(where(is.numeric), mean))

head(Samples_Mean)

reform_means=melt(Samples_Mean)
head(reform_means)
colnames(reform_means)<-c("Sample_ID","Geneexpression")

#mapp clinical data
expression_Clinical=clinical %>% 
  filter(Sample_ID%in%reform_means$Sample_ID) %>% 
  full_join(reform_means,by="Sample_ID") %>% 
  drop_na() %>% 
  select(c("Sample_ID","CancerType","DFS_Status","DSF_Time","Geneexpression"))

dim(expression_Clinical)

#spliting data by median Beta value
Median_expression=median(expression_Clinical$Geneexpression)

expression_Clinical$level = "NA" # creates a new variable filled with NAs
high = expression_Clinical$Geneexpression>=Median_expression
low =  expression_Clinical$Geneexpression<Median_expression
expression_Clinical$level[high]="High"
expression_Clinical$level[low]="Low"

expression_Clinical$event = "NA" # creates a new variable filled with NAs
event_occured = expression_Clinical$DFS_Status=="1:Recurred/Progressed"
no_event =  expression_Clinical$DFS_Status=="0:DiseaseFree"
expression_Clinical$event[event_occured]="1"
expression_Clinical$event[no_event]="0"


#expression_Clinical$CancerType <- as.character(expression_Clinical$CancerType)
expression_Clinical$CancerType[expression_Clinical$CancerType == "Non-Seminomatous Germ Cell Tumor" | expression_Clinical$CancerType == "Non-Seminomatous Germ Cell Tumo" | expression_Clinical$CancerType == "Embryonal Carcinoma" ] <- "Non-Seminoma"

write_tsv(expression_Clinical, path = "/Users/fazal2/Desktop/tidy_practice/B3_Case1_Beta0.4_0.05FDR_RNAseq_Clinical.tsv")

#run KM
fit <- survfit(Surv(as.numeric(expression_Clinical$DSF_Time), as.numeric(expression_Clinical$event)) ~ expression_Clinical$level, data = expression_Clinical)
print(fit)
summary(fit)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


surv_diff <- survdiff(Surv(as.numeric(expression_Clinical$DSF_Time), as.numeric(expression_Clinical$event)) ~ expression_Clinical$level, data = expression_Clinical)
surv_diff

