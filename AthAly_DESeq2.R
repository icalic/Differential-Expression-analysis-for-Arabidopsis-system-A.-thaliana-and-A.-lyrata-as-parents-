library(DESeq2)
directory<-"/Users/irinacalic/Desktop/transcriptome/"
sampleFiles<-grep("txt", list.files(directory), value=TRUE)
sampleCondition<-vector(length=length(sampleFiles))
sampleGenotype<-vector(length=length(sampleFiles))
for (i in 1:length(sampleFiles)){
  a<-unlist(strsplit(sampleFiles[i],"_"))
  sampleGenotype[i]<-a[1]          
  sampleCondition[i]<-a[2]
}
sampleGenotype<-factor(sampleGenotype)
sampleCondition<-factor(sampleCondition)

sampleTable<-data.frame(sampleName = sampleFiles, fileName=sampleFiles, genotype=sampleGenotype, condition=sampleCondition)
sampleTable

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design= ~ condition + genotype + condition:genotype)

#The design formula used here for setting up the ddsHTSeq object for DESeq2 is following:
#design = ~ condition + genotype + genotype:condition

sampleTable
# sampleName                    fileName genotype  condition
# 1    alymn_flowering_93320.txt   alymn_flowering_93320.txt alymn  flowering
# 2    alymn_flowering_93332.txt   alymn_flowering_93332.txt alymn  flowering
# 3    alymn_flowering_93376.txt   alymn_flowering_93376.txt alymn  flowering
# 4     alymn_seedling_93334.txt    alymn_seedling_93334.txt alymn   seedling
# 5     alymn_seedling_93336.txt    alymn_seedling_93336.txt alymn   seedling
# 6     alymn_seedling_93340.txt    alymn_seedling_93340.txt alymn   seedling
# 7     alymn_seedling_93380.txt    alymn_seedling_93380.txt alymn   seedling
# 8   alymn_vegetative_93358.txt  alymn_vegetative_93358.txt alymn vegetative
# 9   alymn_vegetative_93360.txt  alymn_vegetative_93360.txt alymn vegetative
# 10  alymn_vegetative_93362.txt  alymn_vegetative_93362.txt alymn vegetative
# 11  alymn_vegetative_93382.txt  alymn_vegetative_93382.txt alymn vegetative
# 12  colfri_flowering_93366.txt  colfri_flowering_93366.txt colfri  flowering
# 13  colfri_flowering_93368.txt  colfri_flowering_93368.txt colfri  flowering
# 14  colfri_flowering_93370.txt  colfri_flowering_93370.txt colfri  flowering
# 15  colfri_flowering_93384.txt  colfri_flowering_93384.txt colfri  flowering
# 16   colfri_seedling_93318.txt   colfri_seedling_93318.txt colfri   seedling
# 17   colfri_seedling_93322.txt   colfri_seedling_93322.txt colfri   seedling
# 18   colfri_seedling_93324.txt   colfri_seedling_93324.txt colfri   seedling
# 19 colfri_vegetative_93342.txt colfri_vegetative_93342.txt colfri vegetative
# 20 colfri_vegetative_93344.txt colfri_vegetative_93344.txt colfri vegetative
# 21 colfri_vegetative_93348.txt colfri_vegetative_93348.txt colfri vegetative

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design= ~ condition + genotype + condition:genotype)


ddsHTSeq
# class: DESeqDataSet
# dim: 17482 21
# metadata(1): version
# assays(1): counts
# rownames(17482): AT1G01010 AT1G01020 ... AT5G67630 AT5G67640
# rowData names(0):
# colnames(21): alymn_flowering_93320.txt alymn_flowering_93332.txt ... colfri_vegetative_93344.txt colfri_vegetative_93348.txt
# colData names(2): genotype condition

dds_irina<-DESeq(ddsHTSeq)
dds_irina

#to filter for low expressed genes with a mean final read count <10 across samples that passed RNA integrity test

dds_filterlowcount<-estimateSizeFactors(dds_irina)
keep<-rowSums(counts(dds_filterlowcount, normalized=TRUE)>=10)>=21
dds_keep<-dds_filterlowcount[keep,]
dds_keep

#out of 17482 ortholog genes between Athaliana and Alyrata, 12689 genes were yielded after filter for low expressed genes

# likelihood ratio test will determine when adding an interaction genotype:interaction improves the model
ddsLRT_irina <- DESeq(dds_irina, test="LRT", reduced= ~ genotype+condition)##tests whether the interaction brings something
res<-results(ddsLRT_irina)
summary(res)
head(res)# 
sig<-res[which(res$padj<0.05),]###only the gene that have significant genotype:condition interaction
nrow(sig)##2899 genes in this case =>  create boxplot to verify that. 
head(sig)

dds_irina$group <- factor(paste0(dds_irina$genotype, dds_irina$condition))
design(dds_irina) <- ~ group
levels(dds_irina$group)
dds_irina$group<- relevel(dds_irina$group, ref = "alymnseedling")###repeat for all levels to have all pairwise foldchange estimates.
levels(dds_irina$group)
dds_test <- DESeq(dds_irina)
# using pre-existing size factors
# estimating dispersions
# found already estimated dispersions, replacing these
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

resultsNames(dds_test)
# [1] "Intercept" "group_alymnseedling_vs_alymnflowering" "group_alymnvegetative_vs_alymnflowering"  "group_colfriflowering_vs_alymnflowering"
# [5] "group_colfriseedling_vs_alymnflowering" "group_colfrivegetative_vs_alymnflowering"

#write results
colfri_seedling_alymn_seedling<-as.matrix(results(dds_test, contrast=c("group", "colfriseedling", "alymnseedling")))
write.csv(colfri_seedling_alymn_seedling, "/Users/irinacalic/Desktop/colfri_seedling_alymn_seedling.csv")
colfri_vegetative_alymn_vegetative<-as.matrix(results(dds_test, contrast=c("group", "colfrivegetative", "alymnvegetative")))
write.csv(colfri_vegetative_alymn_vegetative, "/Users/irinacalic/Desktop/colfri_vegetative_alymn_vegetative.csv")
colfri_flowering_alymn_flowering<-as.matrix(results(dds_test, contrast=c("group", "colfriflowering", "alymnflowering")))
write.csv(colfri_flowering_alymn_flowering, "/Users/irinacalic/Desktop/colfri_flowering_alymn_flowering.csv")


