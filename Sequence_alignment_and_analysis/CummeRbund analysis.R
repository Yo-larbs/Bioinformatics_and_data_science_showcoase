BiocManager::install("cummeRbund")
install.packages('dplyr')
install.packages("RColorBrewer") #this is a colour palette
install.packages("ggVennDiagram") # to make venn diagram
install.packages("tidyr")
install.packages("ggplot2")
install.packages("ggvennDiagram")
library(cummeRbund) # cummeRbund is the R package designed to read cuffdiff output, allows you to visualise and manipulate data
library(tidyr)
library(dplyr) # data manipulation
library(RColorBrewer) 
library(ggVennDiagram) 
library(ggplot2)

options(max.print = 1000000)

cuff_data <- readCufflinks('outputs/cuffdiff') #cummeRbund directly reads cuffdiff output, so have to give cuffdigg folder

pdf(file="scatterplot_6hr_14hr.pdf") #specifying output file
csScatter(genes(cuff_data), 'X6hr', 'X14hr') #cummeRbund has built in scattergraph function, X6hr and x14hr refer to the 2 time conditions
dev.off()

pdf(file="scatterplot_6hr_26hr.pdf")
csScatter(genes(cuff_data), 'X6hr', 'X26hr')
dev.off()

pdf(file="scatterplot_6hr_26hr.pdf")
csScatter(genes(cuff_data), 'X14hr', 'X26hr')
dev.off()

pdf(file="scatterplot_matrix.pdf")
csScatterMatrix(genes(cuff_data)) # scatter matrix plots scattergraph for all samples
dev.off() #dev.off() stops writing to my the pdf
                
myGeneid=getSig(cuff_data,level="genes",alpha=0.05) #sig gene function allows you to get list of item (in this case genes, could be transcripts etc..) that meet the threshold (in this case 0.05),)
myGene=getGenes(cuff_data,myGeneid) #get gene actually creates object for the gene list you give  it from the cuffset specified, it contains all the information for those genes
                 
                
#pdf(file="heatmap.pdf")
csHeatmap(myGene,cluster='both',labRow=F) # this gives heatmap across 3 samples for the geneset you specify (in this case mygene)
#dev.off()
                
#pdf(file="heatmap_similarity.pdf")
#csDistHeat(
#dev.off
                  

                  
test_1.all_genes = getSig(cuff_data,x='X6hr',y='X14hr',level="genes",alpha=0.99) # here we gather all genes that were tested for differential expression between 6hr and 14hr, and get the information for all of them)
test_1.genes = getGenes(cuff_data,test_1.all_genes)
test_1.upreg=subset(diffData(test_1.genes),sample_1=='X6hr' & sample_2=='X14hr' & q_value<=0.05 & log2_fold_change>2,select= c(gene_id,log2_fold_change,q_value))
#here i create a subset using the object 'test_1.genes'
#from this object i extract all the differential test data using diffData(). since the object holds all the information about the list of genes i gave it (including every test those genes were involved in) i specify that i only want the tests between 6 and 14 hr
# with a siginificant q value (not p value since it runs multiple tests) and a log fold over 2 (big increase in expression), i then create a dataframe holding this information

test_1.upreg_volc=subset(diffData(test_1.genes),sample_1=='X6hr' & sample_2=='X14hr',select= c(gene_id,log2_fold_change,q_value))
#  do a similar thing as above but since this is used for volcano plots i want all tests for those times not just the signifcant ones
test_1.upreg_genes = as.vector(test_1.upreg$gene_id)
#  save the gene id from my dataframe as a vector for further manipulation
test_1.upreg_genes_names=getGenes(cuff_data,test_1.upreg_genes)
#  take this vector for the significantly expressed high log2fold and create an object giving me the information on these genes
test_1.upreg_genes_names=featureNames(test_1.upreg_genes_names)
# cummeRbund saves genes using arbitrary geneID from cufflinks (in this case using XLOC...) featureNames()gets the short gene names which allow further analysis (such as gene ontology)

test_1.upreg_volc$col=ifelse(test_1.upreg_volc$log2_fold_change > 2 & test_1.upreg$q_value<=0.05, "blue","black")
test_1.upreg_volc$col[test_1.upreg_volc$log2_fold_change<=-2 & test_1.upreg$q_value<=0.05]="red"
# here create columns for my volcano plot which assigns colours depending on increased or decreased expression 
                 
test_1.upreg_volc$label=ifelse(test_1.upreg_volc$col == "blue","up","down")
test_1.upreg_volc$label[test_1.upreg_volc$col=="red"]="down"
test_1.upreg_volc$label[test_1.upreg_volc$col=="black"]="no"
#  create another column for labelling in my legend depending on level of differential expression                  

pdf(file="volcano_6hr_14hr.pdf")
ggplot(data=test_1.upreg_volc,aes(x=log2_fold_change,y=-log10(q_value)))+
  geom_point(aes(colour=label),size=0.5)+
  labs(colour="gene expression")+
  theme_minimal()+
  ggtitle("Saccharomyces cerevisiae gene expression 6hrs vs 14hrs")+
  scale_color_manual(values=c("red","black","blue"))+
  xlim(-10, 15)
dev.off()
# here construct volcano plot using ggplot, x is log2fold change and y is -log10(qvalue),  colour the points using colour label and aes
# set the limits of x to make graph more legible since there is an outlier which squishes the rest of the graph


# do the exact same things for test 2 and test 3 (6hr vs 26hr and 14hr v s 26hr)             
test_2.all_genes = getSig(cuff_data,x='X6hr',y='X26hr',level="genes",alpha=0.99)
test_2.genes = getGenes(cuff_data,test_2.all_genes)
test_2.upreg=subset(diffData(test_2.genes),sample_1=='X6hr' & sample_2=='X26hr' & q_value<0.05 & log2_fold_change>2,select= c(gene_id,log2_fold_change,q_value))
test_2.upreg_volc=subset(diffData(test_2.genes),sample_1=='X6hr' & sample_2=='X26hr',select= c(gene_id,log2_fold_change,q_value))
test_2.upreg_genes = as.vector(test_2.upreg$gene_id)
test_2.upreg_genes_names=getGenes(cuff_data,test_2.upreg_genes)
test_2.upreg_genes_names=featureNames(test_2.upreg_genes_names)
test_2.upreg_volc$col=ifelse(test_2.upreg_volc$log2_fold_change > 2 & test_2.upreg$q_value<=0.05, "blue","black")
test_2.upreg_volc$col[test_2.upreg_volc$log2_fold_change<=-2 &test_1.upreg$q_value<=0.05]="red"

test_2.upreg_volc$label=ifelse(test_2.upreg_volc$col == "blue","up","down")
test_2.upreg_volc$label[test_2.upreg_volc$col=="red"]="down"
test_2.upreg_volc$label[test_2.upreg_volc$col=="black"]="no"

                  
pdf(file="volcano_6hr_26hr.pdf")
ggplot(data=test_2.upreg_volc,aes(x=log2_fold_change,y=-log10(q_value),colour=col))+
  geom_point(aes(colour=label),size=0.5)+
  labs(colour="gene expression")+
  theme_minimal()+
  ggtitle("Saccharomyces cerevisiae gene expression 6hrs vs 26hrs")
  scale_color_manual(values=c("red","black","blue"))+
  xlim(-30, 30)
  dev.off()
                  
test_3.all_genes = getSig(cuff_data,x='X14hr',y='X26hr',level="genes",alpha=0.99)
test_3.genes = getGenes(cuff_data,test_3.all_genes)
test_3.upreg=subset(diffData(test_3.genes),sample_1=='X14hr' & sample_2=='X26hr' & q_value<0.05 & log2_fold_change>2,select= c(gene_id,log2_fold_change,q_value))
test_3.upreg_volc=subset(diffData(test_3.genes),sample_1=='X14hr' & sample_2=='X26hr',select= c(gene_id,log2_fold_change,q_value))
test_3.upreg_genes = as.vector(test_3.upreg$gene_id)
test_3.upreg_genes_names=getGenes(cuff_data,test_3.upreg_genes)
test_3.upreg_genes_names=featureNames(test_3.upreg_genes_names)
test_3.upreg_genes_names = cbind(test_3.upreg_genes_names,test_3.upreg)
test_3.upreg_genes_names= test_3.upreg_genes_names[,c(2,4,5)]
test_3.upreg_volc$col=ifelse(test_3.upreg_volc$log2_fold_change > 2 & test_3.upreg$q_value<=0.05, "blue","black")
test_3.upreg_volc$col[test_3.upreg_volc$log2_fold_change<=-2 &test_1.upreg$q_value<=0.05]="red"
test_3.upreg_volc$label=ifelse(test_3.upreg_volc$col == "blue","up","down")
test_3.upreg_volc$label[test_3.upreg_volc$col=="red"]="down"
test_3.upreg_volc$label[test_3.upreg_volc$col=="black"]="no"


pdf(file="volcano_14hr_26hr.pdf")
ggplot(data=test_3.upreg_volc,aes(x=log2_fold_change,y=-log10(q_value),colour=col))+
  geom_point(aes(colour=label),size=0.5)+
  labs(colour="gene expression")+
  theme_minimal()+
  ggtitle("Saccharomyces cerevisiae gene expression 14hrs vs 26hrs")+
  scale_color_manual(values=c("red","black","blue"))+
  xlim(-10, 20)
dev.off()
                  
all_test.sig_genes = getSig(cuff_data,level="genes",alpha=0.99)
all_test.genes = getGenes(cuff_data,all_test.sig_genes)
all_test.upreg=subset(diffData(all_test.genes), q_value<0.05, select= c(gene_id,log2_fold_change,q_value))
all_test.upreg_genes = as.vector(all_test.upreg$gene_id)
all_test.upreg_genes_names=getGenes(cuff_data,all_test.upreg_genes)
all_test.upreg_genes_names_list=featureNames(all_test.upreg_genes_names)  
n_distinct(all_test.upreg_genes_names_list$gene_short_name)

#here i try to recreate the heatmap from paper using our own results. they compared expression for set of genes associated with gluceogenesis across samples
#wanted to see if our results gave similar pattern of expression
                  
gene_comp = c('TDH3','FBA1','GPM1','TPI1','PGI1','ENO1','PGK1','TDH2','FBP1','ERT1','SDL1','GPM3','GPM2','TDH1','MDH2','PYC1','PCK1')
gene_comp_list=getGenes(cuff_data,gene_comp)
gene_matrix = repFpkm(gene_comp_list)
#first of all  create a list containing the genes they used in heatmap and get the Fpkm for those genes across all replicates (3 replicates per time)
# in the heatmap they used Z score

gene_comp_feature = featureNames(gene_comp_list)
# get the actual gene short names for the genes as a dataframe with GeneID
gene_matrix_pivot=gene_matrix %>% pivot_wider(names_from = rep_name, values_from = fpkm)
# pivot wider before to help manipulate data, get list of genes with fpkm for each sample as columns
gene_matrix_pivot=gene_matrix_pivot[,c(1,9:17)]
#  keep only the geneID [1], the and the fpkm for each replciate [9:17]
gene_matrix_pivot=gene_matrix_pivot %>%group_by(gene_id)%>%summarise_all(na.omit)
# because of the pivot wider data frame had multiple rows of the same geneID, each row only having one entry for fpkm across the columns]
#summarise by geneID , and get rid of NA, meaning all values for that geneID combined into one row giving desired dataframe shape

gene_matrix_zscore <- t(scale(t(gene_matrix_pivot[,2:10])))
# apply z score but first  transpose it because scale applies columnwise ,  then transpose it back
gene_matrix_zscore=gene_matrix_zscore[,c(2,1,3,5,4,6,7,8,9)] 
#  reorder the columns to get replicates in order e.g 6hr_0, 6hr_1, 6hr_2, 14hr_0 etc
gene_matrix_zscore=cbind(gene_comp_feature$gene_short_name,gene_matrix_zscore)
#  attach the gene short names back so each gene is recognisable by short name 
colnames(gene_matrix_zscore)[1]='Gene'
#  rename the column to gene (was given internal name)
gene_matrix_zscore=as_data_frame(gene_matrix_zscore)
# turn it back into dataframe (was matrix before hand)
gene_heatmap=gene_matrix_zscore %>% pivot_longer(!Gene,names_to = "Sample",values_to = "z score")
# for the heatmap need to get it back to longer format (e.g each row having gene Id, which replicate and zscore)
#used pivot longer for this (pretty much reversed pivot wider but with z score instead of fpkm)
gene_heatmap$`z score`=as.numeric(gene_heatmap$`z score`)
# had to convert z score back to numerical value
gene_heatmap$Sample <- factor(gene_heatmap$Sample, levels=c("X6hr_0", "X6hr_1", "X6hr_2","X14hr_0","X14hr_1","X14hr_2","X26hr_0","X26hr_1","X26hr_2"))
gene_heatmap$Gene <- factor(gene_heatmap$Gene, levels=gene_comp[17:1])
# 2 lines above force the heatmap to display in specified order, done to make it visually comparable to original papers

pdf(file="heat map z score.pdf")
ggplot(gene_heatmap, aes(x=Sample,y=Gene, fill= `z score`)) + 
  geom_tile()+
  theme(legend.position = "top")+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(legend.text=element_text(size=10),
  legend.title=element_text(size=15))+
  scale_fill_distiller(palette = "YlOrRd",direction = 1)
#construct heatmap using ggplot
dev.off()
candidate_genes=subset(test_3.upreg_genes_names,!(test_3.upreg_genes_names$gene_short_name%in%test_1.upreg_genes_names$gene_short_name))

# now trying to get a venn diagram containing list of genes, wanted to calculate it mathematically first before using graphing to check 
test_1.sig_genes=getSig(cuff_data,level="genes",alpha=0.05) # get list of significant genes
test_1.sig_genes_list=getGenes(cuff_data,test_1.sig_genes) #get all data on those genes
test_1.sig_list=subset(diffData(test_1.sig_genes_list),sample_1=='X6hr' & sample_2=='X14hr' & q_value<0.05,select= c(gene_id,log2_fold_change,q_value))
#get dataframe specifying specific test we want (test1 6 vs 14) with specified significance
test_1.sig_genes_names=getGenes(cuff_data,test_1.sig_list$gene_id)
# get list of significantly differential genes
test_1.sig_genes_names=featureNames(test_1.sig_genes_names)
# get the list of short gene names 

# do the same for the other 2 tests

test_2.sig_genes=getSig(cuff_data,level="genes",alpha=0.05)
test_2.sig_genes_list=getGenes(cuff_data,test_2.sig_genes)
test_2.sig_list=subset(diffData(test_2.sig_genes_list),sample_1=='X6hr' & sample_2=='X26hr' & q_value<0.05,select= c(gene_id,log2_fold_change,q_value))
test_2.sig_genes_names=getGenes(cuff_data,test_2.sig_list$gene_id)
test_2.sig_genes_names=featureNames(test_2.sig_genes_names)

test_3.sig_genes=getSig(cuff_data,level="genes",alpha=0.99)
test_3.sig_genes_list=getGenes(cuff_data,test_3.sig_genes)
test_3.sig_list=subset(diffData(test_3.sig_genes_list),sample_1=='X14hr' & sample_2=='X26hr' & q_value<0.05,select= c(gene_id,log2_fold_change,q_value))
test_3.sig_genes_names=getGenes(cuff_data,test_3.sig_list$gene_id)
test_3.sig_genes_names=featureNames(test_3.sig_genes_names)



x14v26=intersect(test_3.sig_genes_names$gene_short_name,test_2.sig_genes_names$gene_short_name)
#here i use the intersect function to find how many genes intersect between the two sets (for test 3(14vs26)and test 2(6vs26))
x6v14=intersect(test_1.sig_genes_names$gene_short_name,test_3.sig_genes_names$gene_short_name)
# i do the exact same things for test 1 and test 3 (6hr vs 14hr and 14hr v s 26hr)
x6v26=intersect(test_2.sig_genes_names$gene_short_name, test_1.sig_genes_names$gene_short_name)
# i do the exact same things for test 2 and test 1 (6hr vs 26hr and 6hr vs 14hr)

x6x14x26=(intersect(intersect(test_2.sig_genes_names$gene_short_name,test_1.sig_genes_names$gene_short_name),test_3.sig_genes_names$gene_short_name))
length(intersect(x6x14x26,x6v14))      
#get intersect of all 3 lists, if you imagine a 3 circle venn diagram it would be the part in middle

cuff_data
# now use ggvennDiagram to plot out venn diagram
x = list(test_1.sig_genes_names$gene_short_name,test_2.sig_genes_names$gene_short_name,test_3.sig_genes_names$gene_short_name)

png(file="genes venn diagram.png")
ggVennDiagram(x,label = "count",c("6hr vs 14hr ","6hr vs 26hr","14hr vs 26hr"),label_size = 9,label_alpha = 0,set_size = 6)+scale_fill_gradient(low = "#FFFFFF", high = "#FFFFFF") + theme(legend.position = "none")# set background to white FFFFFFF
dev.off()

sink("list for gene ontology.txt")
noquote(test_3.upreg_genes_names$gene_short_name)
sink()

test_3.sig_list_GO = subset(test_3.sig_list)
test_3.sig_list_GO_names= getGenes(cuff_data,test_3.sig_list_GO$gene_id)
test_3.sig_list_GO_names=featureNames(test_3.sig_list_GO_names)


sink("list 2 for gene ontology.txt")
noquote(test_3.sig_list_GO_names$gene_short_name)
sink()

top17=top_n(test_3.upreg,17,log2_fold_change)
top17_gene=getGenes(cuff_data,top17$)
top17_gene_matrix = repFpkm(top17_gene)
#first of all  create a list containing the genes they used in heatmap and get the Fpkm for those genes across all replicates (3 replicates per time)
# in the heatmap they used Z score

MAplot(genes(cuff_data),"X6hr","X14hr",logMode = T,pseudocount = 5,useCount=T,smooth=FALSE)
MAplot(genes(cuff_data),"X6hr","X26hr",logMode = T,pseudocount = 5,useCount=T,smooth=FALSE)
MAplot(genes(cuff_data),"X14hr","X26hr",logMode = T,pseudocount = 5,useCount=T,smooth=FALSE)

features(getGene(cuff_data,'FBP1'))
