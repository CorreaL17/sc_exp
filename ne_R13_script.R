#Load the required Libraries to run the experiment 
library(readr)
library(dplyr)
library(ggplot2)
library(Seurat)



 
setwd("/Users/luis/Desktop/LC_Analysis/sc_from_toxpaper/Unzipped_files")


#This list all the directories in my folder. Each sample for the single cells is stored as a folder here
dirs<-list.dirs(path=".")

#Print to see that we actually have the folder we want. I see theres some folders I do not want in here. So I will truncate to pic the single cell ones
print(dirs)
sc_dirs<-dirs[2:9]

#I will generate a list of sample names by selecting the information from the folder name that comes after scRNA
sampleName=gsub("./scRNA_","",sc_dirs)


#I will create an empty list to put in all the samples in 
list_cds<-list()

#the following loop will go through the samples an initialize the 10X data to an item in spot i for the length of samples i have. This will also crease seurat object
# A gene willl be kept if it is expressed at least within 3 cells 
# A cell will be kept if it has 200 genes expressed 
for(i in 1:(length(sc_dirs))){
  cds<-Read10X(sc_dirs[i])
  CD8 <- CreateSeuratObject(counts = cds, project = sampleName[i], min.cells = 3, min.features = 200)
  list_cds[[i]]<-CD8
}

#Individually normalizing, finding variable features and removing cells with too many mitochondiral genes 
for (i in 1:length(list_cds)) {
  list_cds[[i]] <- NormalizeData(list_cds[[i]], verbose = FALSE)
  list_cds[[i]] <- FindVariableFeatures(list_cds[[i]], selection.method = "vst", 
                                        nfeatures = 2000, verbose = FALSE)
  list_cds[[i]][["percent.mt"]]<- PercentageFeatureSet(list_cds[[i]], pattern = "^mt-")
}


#Here is a violing plot of 1 sample for Fature number, counts per cells, percent of counts that are mitochondrial 
VlnPlot(list_cds[[1]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



for (i in 1:length(list_cds)) {
  list_cds[[i]] <- subset(list_cds[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
}






#Seurat objects arre then merged into 
big_cds<-merge(list_cds[[1]],
               c(list_cds[[2]],
                 list_cds[[3]],
                 list_cds[[4]],
                 list_cds[[5]],
                 list_cds[[6]],
                 list_cds[[7]],
                 list_cds[[8]]),
              project="TEX")


big_cds
big_cds <- NormalizeData(big_cds)
big_cds <- FindVariableFeatures(big_cds, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(big_cds)
big_cds <- ScaleData(big_cds, features = all.genes)
big_cds <- RunPCA(big_cds, features = VariableFeatures(object = big_cds))
big_cds<-RunTSNE(big_cds,dims = 1:5)

ElbowPlot(big_cds)

#Elbow plot reveals that most of the variance is explained withing the first 5 PCs 
big_cds <- RunUMAP(big_cds, dims = 1:5, group.by="active")






#generating plot for DANG
DimPlot(big_cds, reduction = "umpa",label=TRUE, group.by = "orig.ident")
ggsave("UMAP-batching.jpg",dpi = 400)



#lookin thoruhg genes in the principle components 
VizDimLoadings(big_cds, dims =2, reduction = "pca")



#Attempt at clustering 
big_cds <- FindNeighbors(big_cds, reduction = "umap", dims = 1:5)
big_cds <- FindClusters(big_cds, resolution = 0.3)



FeaturePlot(big_cds, features = c("Me1","Me2","Tox","Lipa"))


write_rds(big_cds,"Tox_Seurat.rds")







