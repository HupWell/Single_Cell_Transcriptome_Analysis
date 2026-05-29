#wget -r -np -k -p -e robots=off https://ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/GSE176078/suppl/ #Terminal

setwd("~/ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/GSE240112/suppl")  
untar("GSE240112_RAW.tar",exdir = "GSE240112_RAW")
dir='GSE240112_RAW/'
## 整理文件夹
library(stringr)
# 列出所有文件
files <- list.files("~/ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/GSE240112/suppl/GSE240112_RAW/",full.names = T)
files 

# 遍历文件并移动到对应的样本文件夹
for (file in files) {
  # 提取样本名
  #file <- files[1]
  sample_name <- gsub("_barcodes.tsv.gz|_features.tsv.gz|_matrix.mtx.gz","", basename(file))
  # 创建样本文件夹（如果不存在）
  sample_folder <- file.path("~/ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/GSE240112/suppl/GSE240112_RAW", sample_name)
  if (!dir.exists(sample_folder)) {
    dir.create(sample_folder, showWarnings = FALSE)
  }
  # 移动文件到对应的样本文件夹
  target_file <- file.path(sample_folder, gsub(paste0(sample_name,"_"),"", basename(file)) )
  file.rename(file, target_file)
}
library(fs)
dir_tree("~/ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/GSE240112/suppl/GSE240112_RAW")
###### step1: 导入数据 ######   
samples <- list.dirs("~/ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/GSE240112/suppl/GSE240112_RAW", recursive = F, full.names = F)
samples
scRNAlist <- lapply(samples, function(pro){
  #pro <- samples[1]
  print(pro)
  folder <- file.path("~/ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/GSE240112/suppl/GSE240112_RAW", pro)
  folder
  counts <- Read10X(folder, gene.column = 2)
  sce <- CreateSeuratObject(counts, project=pro, min.cells = 3)
  return(sce)
})
names(scRNAlist) <-  samples
scRNAlist
PT1 <- scRNAlist[1]
View(PT1)
do.call(rbind,lapply(scRNAlist, dim))
sce.all=merge(x=scRNAlist[[1]],
              y=scRNAlist[ -1 ],
              add.cell.ids = samples  ) 
names(sce.all@assays$RNA@layers)
sce.all[["RNA"]]$counts 
LayerData(sce.all, assay = "RNA", layer = "counts")
sce.all
sce.all <- JoinLayers(sce.all)
sce.all
dim(sce.all[["RNA"]]$counts )
as.data.frame(sce.all@assays$RNA$counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all$orig.ident) 
length(sce.all$orig.ident)
library(stringr)
phe = sce.all@meta.data
table(phe$orig.ident)
phe$group = str_split(phe$orig.ident,'[_]',simplify = T)[,2] 
phe$tissue <- ifelse(phe$orig.ident %in% c("GSM7681687_PT1","GSM7681688_PT2","GSM7681689_PT5"),"PT",
                     ifelse(phe$orig.ident %in% c("GSM7681685_NT7","GSM7681686_NT8"),"NT","RT"))
table(phe$tissue)
sce.all@meta.data = phe
phe$tissue <- factor(phe$tissue, levels = c("PT", "RT","NT"))
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(data.table)
library(dplyr)
#计算线粒体基因比例
mito_genes=rownames(sce.all)[grep("^MT-", rownames(sce.all),ignore.case = T)] 
print(mito_genes) #可能是13个线粒体基因，小鼠数据基因名为小写"^mt-"
#sce.all=PercentageFeatureSet(sce.all, "^MT-", col.name = "percent_mito")
sce.all=PercentageFeatureSet(sce.all, features = mito_genes, col.name = "percent_mito")
fivenum(sce.all@meta.data$percent_mito)
#计算核糖体基因比例
ribo_genes=rownames(sce.all)[grep("^Rp[sl]", rownames(sce.all),ignore.case = T)]
print(ribo_genes)
sce.all=PercentageFeatureSet(sce.all,  features = ribo_genes, col.name = "percent_ribo")
fivenum(sce.all@meta.data$percent_ribo)
#计算红血细胞基因比例
Hb_genes=rownames(sce.all)[grep("^Hb[^(p)]", rownames(sce.all),ignore.case = T)]
print(Hb_genes)
sce.all=PercentageFeatureSet(sce.all,  features = Hb_genes,col.name = "percent_hb")
fivenum(sce.all@meta.data$percent_hb)
head(sce.all@meta.data)
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito",
           "percent_ribo", "percent_hb")
feats <- c("nFeature_RNA", "nCount_RNA")
p1=VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
  NoLegend()
p1 
w=length(unique(sce.all$orig.ident))/3+5;w
ggsave(filename="Vlnplot1.pdf",plot=p1,width = w,height = 5)
feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2=VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3, same.y.lims=T) + 
  scale_y_continuous(breaks=seq(0, 100, 5)) +
  NoLegend()
p2 
w=length(unique(sce.all$orig.ident))/2+5;w
ggsave(filename="Vlnplot2.pdf",plot=p2,width = w,height = 5)

p3=FeatureScatter(sce.all, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
p3
ggsave(filename="Scatterplot.pdf",plot=p3)
if(F){
  selected_c <- WhichCells(sce.all, expression = nFeature_RNA > 500)
  selected_f <- rownames(sce.all)[Matrix::rowSums(sce.all@assays$RNA$counts > 0 ) > 3]
  sce.all.filt <- subset(sce.all, features = selected_f, cells = selected_c)
  dim(sce.all) 
  dim(sce.all.filt) 
}
sce.all.filt =  sce.all
# par(mar = c(4, 8, 2, 1))
# 这里的C 这个矩阵，有一点大，可以考虑随抽样 
C=subset(sce.all.filt,downsample=100)@assays$RNA$counts
dim(C)
C=Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100

most_expressed <- order(apply(C, 1, median), decreasing = T)[50:1]

pdf("TOP50_most_expressed_gene.pdf",width=14)
boxplot(as.matrix(Matrix::t(C[most_expressed, ])),
        cex = 0.1, las = 1, 
        xlab = "% total count per cell", 
        col = (scales::hue_pal())(50)[50:1], 
        horizontal = TRUE)
dev.off()
rm(C)
selected_mito <- WhichCells(sce.all.filt, expression = percent_mito < 25)
selected_ribo <- WhichCells(sce.all.filt, expression = percent_ribo > 3)
selected_hb <- WhichCells(sce.all.filt, expression = percent_hb < 1 )
length(selected_hb)
length(selected_ribo)
length(selected_mito)
sce.all.filt <- subset(sce.all.filt, cells = selected_mito)
sce.all.filt <- subset(sce.all.filt, cells = selected_ribo)
sce.all.filt <- subset(sce.all.filt, cells = selected_hb)
dim(sce.all.filt)
table(sce.all.filt$orig.ident)
length(sce.all.filt$orig.ident)
feats <- c("nFeature_RNA", "nCount_RNA")
p1_filtered=VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
  NoLegend()
w=length(unique(sce.all.filt$orig.ident))/3+5;w 
ggsave(filename="Vlnplot1_filtered.pdf",plot=p1_filtered,width = w,height = 5)

feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2_filtered=VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + 
  NoLegend()
w=length(unique(sce.all.filt$orig.ident))/2+5;w 
ggsave(filename="Vlnplot2_filtered.pdf",plot=p2_filtered,width = w,height = 5) 
sce.all.filt <- NormalizeData(sce.all.filt, 
                              normalization.method = "LogNormalize",
                              scale.factor = 1e4) 
sce.all.filt <- FindVariableFeatures(sce.all.filt)
p4 <- VariableFeaturePlot(sce.all.filt) 
p4
sce.all.filt <- ScaleData(sce.all.filt)
sce.all.filt <- RunPCA(sce.all.filt, features = VariableFeatures(object = sce.all.filt))
##可视化PCA结果
VizDimLoadings(sce.all.filt, dims = 1:2, reduction = "pca")
DimPlot(sce.all.filt, reduction = "pca") + NoLegend()
DimHeatmap(sce.all.filt, dims = 1:12, cells = 500, balanced = TRUE)
seuratObj <- RunHarmony(sce.all.filt, "orig.ident")
names(seuratObj@reductions)
seuratObj <- RunUMAP(seuratObj,  dims = 1:15, 
                     reduction = "harmony")
DimPlot(seuratObj,reduction = "umap",label=F ) 
seuratObj <- RunTSNE(seuratObj, dims = 1:15, 
                     reduction = "harmony")

DimPlot(seuratObj,reduction = "tsne",label=F ) 
sce.all.filt=seuratObj

sce.all.filt <- FindNeighbors(sce.all.filt, reduction = "harmony",
                              dims = 1:15) 

sce.all.filt.all=sce.all.filt
#设置不同的分辨率，观察分群效果(选择哪一个？)
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  sce.all.filt.all=FindClusters(sce.all.filt.all, #graph.name = "CCA_snn", 
                                resolution = res, algorithm = 1)
}
colnames(sce.all.filt.all@meta.data)
apply(sce.all.filt.all@meta.data[,grep("RNA_snn",colnames(sce.all.filt.all@meta.data))],2,table)

p1_dim=plot_grid(ncol = 3, DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.01") + 
                   ggtitle("louvain_0.01"), DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.1") + 
                   ggtitle("louvain_0.1"), DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.2") + 
                   ggtitle("louvain_0.2"))
ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_low.pdf",width = 14)

p1_dim=plot_grid(ncol = 3, DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.8") + 
                   ggtitle("louvain_0.8"), DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.1") + 
                   ggtitle("louvain_1"), DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.3") + 
                   ggtitle("louvain_0.3"))
ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_high.pdf",width = 18)

p2_tree=clustree(sce.all.filt.all@meta.data, prefix = "RNA_snn_res.")
ggsave(plot=p2_tree, filename="Tree_diff_resolution.pdf")
table(sce.all.filt.all@active.ident) 

sel.clust = "RNA_snn_res.0.8"
sce.all.int <- SetIdent(sce.all.filt.all, value = sel.clust)
table(sce.all.int@active.ident) 
colnames(sce.all.int@meta.data) 
dir.create("./3-Celltype")
setwd("./3-Celltype")
scRNA=sce.all.int
# 加载参考数据集
library(SingleR)
library(celldex)
refRNA <- celldex::HumanPrimaryCellAtlasData()

# 提取表达矩阵
data_matrix <- GetAssayData(scRNA, slot = "data")

# 使用 SingleR 进行细胞类型注释
annotations <- SingleR(test = data_matrix, ref = refRNA, labels = refRNA$label.main)

# 将注释结果添加到 Seurat 对象
scRNA$labels <- annotations$labels

# 可视化注释结果
DimPlot(scRNA, group.by = "labels", reduction = "umap", label = TRUE)
# DEG analysis
table(scRNA@meta.data$tissue) 

Idents(scRNA) <- scRNA$tissue
Idents(scRNA) <- scRNA@meta.data$tissue

deg_results <- FindMarkers(scRNA, ident.1 = "PT", ident.2 = "NT",
                           logfc.threshold = 0.25, min.pct = 0.25)
deg_results$gene <- rownames(deg_results)
write.csv(deg_results, "deg_normal_vs_tumor.csv", row.names = FALSE)
library(sceasy)
scRNA[["RNA"]] <- as(scRNA[["RNA"]], Class = "Assay")

sceasy::convertFormat(scRNA, from="seurat", to="anndata",
                      outFile='scRNA.h5ad')