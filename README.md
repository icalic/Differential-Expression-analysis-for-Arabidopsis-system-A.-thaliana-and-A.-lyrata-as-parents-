# Differential-Expression-analysis-for-Arabidopsis-system-A.-thaliana-and-A.-lyrata-as-parents-
# DESeq analysis in R for parental transcriptome data (A. thaliana and A. lyrata) 

#Initial step is to generate count reads associated with expression with HTSeq. 
#to quantify gen expression, it can be used following command within HTSeq 

htseq-count [options] <alignment.sam> <annotation.gtf> > /path/file_count.count

#note: alignment files are unique sam files which are derived within mapping, where gff file is annotation for corresponding reference genotype

#DESeq2 provides analysis of count data from RNA_Seq to detect differentially expressed genes. DESeq works based on negative binomial generalized linear model.
#As input files is necessary to gather all count file that passed RNA Integrity check (for RNA Integrity QC, please refer to [the](https://github.com/icalic/RNA-Integrity-)).
#All count files should be in txt format placed at the working directory.
#If working with two or more different species (e.g., close relatives such as Arabidopsis thaliana, Arabidopsis lyrata or Arabidopsis halleri), it is necessary to obtain ortholog genes within species. 

#DESeq2 is run in R using Bioconductor available: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#citation:
Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi:10.1186/s13059-014-0550-8.
