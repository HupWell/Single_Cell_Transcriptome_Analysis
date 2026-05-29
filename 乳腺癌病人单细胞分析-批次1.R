rm(list = ls()) # 清空所有变量
# 下载 GSE176078
#wget -r -np -k -p -e robots=off https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/ #Terminal
if ("package:Seurat" %in% search()) {
  detach("package:Seurat", unload = TRUE)
}
if ("package:SeuratObject" %in% search()) {
  detach("package:SeuratObject", unload = TRUE)
}
gz_files <- list.files(pattern = "^GSM\\d+_CID\\d+\\.tar\\.gz$")
i <- 2
sample_name <- gsub("\\.tar\\.gz$", "", gz_files[i])
if (!dir.exists(sample_name)) {
  dir.create(sample_name)
}
untar(gz_files[i], exdir = sample_name)
library(fs)
# 2. 找到barcodes.tsv|genes.tsv|matrix.mtx
all_files <- list.files(path = sample_name,all.files = T,recursive = T,full.names = T)
all_files
f <- paste0(sample_name,"/",basename(all_files))
f
file.rename(all_files[1], f[1])
file.rename(all_files[2], f[2])
file.rename(all_files[3], f[3])

# 删除剩余的空文件夹
dir_path <- list.dirs(path = sample_name,recursive = F)
dir_path
unlink(dir_path, recursive = TRUE)
for(i in 1:length(gz_files)) {
  # 这里i=1注释掉，试试看一个样本的结果，理解，然后再走for循环
  # i <- 2
  # 提取样本名称：去除.raw_gene_bc_matrices.tar.gz后缀
  # 例如：GSM4546310_Pat06-A.raw_gene_bc_matrices.tar.gz → GSM4546310_Pat06-A
  sample_name <- gsub("\\.tar\\.gz$", "", gz_files[i])
  sample_name
  print(sample_name)
  
  # 创建目标目录
  if (!dir.exists(sample_name)) {
    dir.create(sample_name)
  }
  
  untar(gz_files[i], exdir = sample_name)
  
  # 2. 找到barcodes.tsv|genes.tsv|matrix.mtx
  all_files <- list.files(path = sample_name,all.files = T,recursive = T,full.names = T)
  all_files
  f <- paste0(sample_name,"/",basename(all_files))
  f
  file.rename(all_files[1], f[1])
  file.rename(all_files[2], f[2])
  file.rename(all_files[3], f[3])
  
  # 删除剩余的空文件夹
  dir_path <- list.dirs(path = sample_name,recursive = F)
  dir_path
  unlink(dir_path, recursive = TRUE)
}
dir_tree("./")
rename_10x_files_uncompressed <- function(files) {
  renamed_count <- 0
  for (file in files) {
    base_name <- basename(file)
    dir_name <- dirname(file)
    
    # 智能识别文件类型，并指定新的标准名称
    if (grepl("barcode", base_name, ignore.case = TRUE)) {
      new_name <- file.path(dir_name, "barcodes.tsv")
    } else if (grepl("gene", base_name, ignore.case = TRUE)) {
      # 注意：10X的 "genes" 文件对应 Seurat 的 "features"
      new_name <- file.path(dir_name, "features.tsv") 
    } else if (grepl("matrix", base_name, ignore.case = TRUE)) {
      new_name <- file.path(dir_name, "matrix.mtx")
    } else {
      next # 如果不是这三个文件，跳过
    }
    
    # 只有当文件名需要改变时才执行重命名
    if (file != new_name) {
      tryCatch({
        file.rename(from = file, to = new_name)
        cat("成功重命名:", base_name, "->", basename(new_name), "\n")
        renamed_count <- renamed_count + 1
      }, error = function(e) {
        cat("错误：无法重命名", base_name, ":", conditionMessage(e), "\n")
      })
    }
  }
  return(renamed_count)
}
rename_10x_files <- function(files) {
  renamed <- FALSE
  for (file in files) {
    base_name <- basename(file)  # 获取文件名（不含路径）
    
    # 智能识别文件类型
    if (grepl("barcode", base_name, ignore.case = TRUE)) {
      new_name <- file.path(dirname(file), "barcodes.tsv.gz")
    } else if (grepl("feature|gene", base_name, ignore.case = TRUE)) {
      new_name <- file.path(dirname(file), "features.tsv.gz")
    } else if (grepl("matrix", base_name, ignore.case = TRUE)) {
      new_name <- file.path(dirname(file), "matrix.mtx.gz")
    } else {
      next  # 如果不是这三个文件，跳过
    }
    
    # 只有当文件名不同时才执行重命名
    if (file != new_name) {
      file.rename(file, new_name)
      cat("重命名成功:", base_name, "->", basename(new_name), "\n")
      renamed <- TRUE
    }
  }
  return(renamed)
}
sample_dirs <- list.dirs(path = "./", full.names = TRUE, recursive = FALSE)
files_to_rename <- list.files(path = sample_dirs, 
                              pattern = "\\.tsv$|\\.mtx$", # 匹配以 .tsv 或 .mtx 结尾的文件
                              full.names = TRUE)
for (s_dir in sample_dirs) {
  cat("\n--- 正在检查文件夹:", basename(s_dir), "---\n")
  
  # a. 找到文件夹内所有的 .tsv 和 .mtx 文件
  files_to_rename <- list.files(path = s_dir, 
                                pattern = "\\.tsv$|\\.mtx$", # 匹配以 .tsv 或 .mtx 结尾的文件
                                full.names = TRUE)
  
  if (length(files_to_rename) > 0) {
    # b. 调用我们修改后的函数来重命名文件
    count <- rename_10x_files_uncompressed(files_to_rename)
    if (count == 0) {
      cat("文件已是标准格式，无需重命名。\n")
    }
  } else {
    cat("未找到需要处理的文件。\n")
  }
}
samples <- list.dirs("./", recursive = F, full.names = F)
samples
seurat_objects_list <- list()
# -------------------- 步骤 1: 加载 Seurat --------------------
library(Seurat)
library(harmony)
cat("Seurat 已加载。当前版本: ", as.character(packageVersion("Seurat")), "\n")


# -------------------- 步骤 2: 读取所有 10X 数据并创建 Seurat 对象列表 --------------------
# 定义数据所在的主目录
untar("176078_RAW.tar",exdir = "GSE176078_RAW")
main_data_path <- "~/ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE240112_RAW" 
setwd("~/ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_RAW" )

# 获取所有样本文件夹的路径 (以 "GSM" 开头)
sample_dirs <- list.dirs(path = main_data_path, full.names = TRUE, recursive = FALSE)
sample_dirs <- sample_dirs[grepl("GSM", basename(sample_dirs))]

# 创建一个空列表来存放 Seurat 对象
seurat_objects_list <- list()

# 循环读取每个样本文件夹
for (folder in sample_dirs) {
  sample_name <- basename(folder)
  cat("\n--- 正在处理样本:", sample_name, "---\n")
  
  # 使用 Read10X 读取数据
  # 关键点1: gene.column=1 
  counts <- Read10X(data.dir = folder, gene.column = 1)
  
  # 使用 CreateSeuratObject 创建 Seurat 对象
  # 关键点2: 确保 s_obj 是一个 Seurat 对象
  s_obj <- CreateSeuratObject(counts = counts, 
                              project = sample_name,
                              min.cells = 3,
                              min.features = 200)
  
  # 关键点3: 将创建好的 Seurat 对象存入列表
  seurat_objects_list[[sample_name]] <- s_obj
  
  cat("样本", sample_name, "已创建为 Seurat 对象。\n")
  # 我们可以立即验证一下
  cat("刚刚创建的对象的类型是: ", class(seurat_objects_list[[sample_name]]), "\n")
}


# -------------------- 步骤 3: 合并 Seurat 对象 --------------------
# 现在 seurat_objects_list 绝对是一个包含 Seurat 对象的列表了

# 准备样本名用于添加细胞ID前缀
sample_names_for_merge <- names(seurat_objects_list)

# 判断是否需要合并
if (length(seurat_objects_list) > 1) {
  cat("\n--- 正在合并", length(seurat_objects_list), "个 Seurat 对象... ---\n")
  
  # 直接调用 merge()，让 Seurat v5 的 S4 方法派发机制工作
  # 这是 Seurat v5 的标准合并方法
  merged_seurat <- merge(x = seurat_objects_list[[1]], 
                         y = seurat_objects_list[2:length(seurat_objects_list)], 
                         add.cell.ids = sample_names_for_merge)
  
} else if (length(seurat_objects_list) == 1) {
  cat("\n--- 只有一个样本，无需合并。正在准备对象... ---\n")
  merged_seurat <- seurat_objects_list[[1]]
  # 即使只有一个样本，最好也用 RenameCells 统一细胞名称格式，便于后续处理
  # merged_seurat <- RenameCells(merged_seurat, add.cell.id = names(seurat_objects_list)[1])
  
} else {
  stop("错误：没有找到任何可以处理的样本！")
}


# -------------------- 步骤 4: (仅限多样本合并后) 整合 Layers --------------------
# 这一步只有在进行了 merge 操作后才需要
if (length(seurat_objects_list) > 1) {
  cat("\n--- 正在使用 JoinLayers 整合数据层... ---\n")
  merged_seurat <- JoinLayers(merged_seurat)
}


# -------------------- 最终检查 --------------------
cat("\n--- 处理完成！ ---\n")

# 检查最终对象的类型
cat("最终对象的类型是: ", class(merged_seurat), "\n")

# 打印最终对象，查看摘要信息
print(merged_seurat)

# 数据预处理与分析
merged_seurat <- merged_seurat %>%
  Seurat::NormalizeData(verbose = FALSE) %>%  # 归一化数据
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%  # 选择2000个变异基因
  ScaleData(verbose = FALSE) %>%  # 数据缩放
  RunPCA(npcs = 50, verbose = FALSE)  # 主成分分析（PCA）

# Harmony批次效应校正
merged_seurat <- merged_seurat %>% RunHarmony(group.by.vars = "orig.ident")

# 获取Harmony嵌入数据
harmony_embeddings <- Embeddings(merged_seurat, 'harmony')

# 降维和可视化
dims = 1:30
merged_seurat<- merged_seurat %>%
  RunUMAP(reduction = "harmony", dims = dims) %>%  # 使用Harmony嵌入做UMAP
  RunTSNE(reduction = "harmony", dims = dims)  %>% # 使用Harmony嵌入做t-SNE
  FindNeighbors(reduction = "harmony", dims = dims) 
merged_seurat <- FindClusters(merged_seurat, 
                              reduction = "harmony", 
                              dims = dims, 
                              resolution = 0.5) 
# 可视化结果
DimPlot(merged_seurat, reduction = "umap")  # UMAP可视化
DimPlot(merged_seurat, reduction = "tsne")  # t-SNE可视化


library(celldex)
library(assertthat)
library(monocle)
library(Seurat)
library(tidyverse)
library(Matrix)
library(stringr)
library(dplyr)
library(tricycle)   
library(scattermore)
library(scater)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleR)
library(CCA)
library(clustree)
library(cowplot)
library(monocle)
library(tidyverse)
library(SCpubr)
library(UCell)
library(irGSEA)
library(GSVA)
library(GSEABase)
library(harmony)
library(plyr)
library(randomcoloR)
library(CellChat)
library(future)
library(ggplot2)
library(ggforce)
library(ggsci)
genes <- list("Breast_Cancer_Cells"=c('EPCAM','KRT19','KRT8','ERBB2','ESR1'),
              "B_cells"=c("CD79A",  "MS4A1"),
              "Myeloid"=c("CSF1R", "CSF3R", "CD68"),
              "T_cells"=c('CD3D','CD3E','CD8A','CD4'),
              "NK_cells"=c('NCAM1','GNLY','NKG7','KLRD1'),
              "Macrophages"=c('MRC1','MSR1','APOE'),
              "Endothelial"=c('PECAM1','VWF','CDH5','CLDN5'),
              "Plasma"=c("MZB1", "IGHG1"),
              "Fibroblasts"=c("PDGFRB" ,"LUM")
)

DotPlot(merged_seurat, features = genes,cols = "RdYlBu") +
  RotatedAxis()
ann.ids <- c("T_cells", #cluster0
             "NK_cells",#cluster1
             "Breast_Cancer_Cells",#cluster2
             "Myeloid",#cluster3
             "Fibroblasts",#cluster4
             "Endothelial",#cluster5
             "Fibroblasts",#cluster6
             "Breast_Cancer_Cells",#cluster7
             "T_cells",#cluster8
             "Breast_Cancer_Cells",#cluster9
             "B_cells",#cluster10
             "B_cells",#cluster11
             "Breast_Cancer_Cells",#cluster12
             "T_cells",#cluster13
             "Breast_Cancer_Cells",#cluster14
             "T_cells",#cluster15
             "Myeloid",#cluster16
             "Plasma",#cluster17
             "Endothelial",#cluster18
             "Breast_Cancer_Cells",#cluster19
             "Breast_Cancer_Cells",#cluster20
             "Breast_Cancer_Cells",#cluster21
             "Breast_Cancer_Cells",#cluster22
             "Macrophages",#cluster23
             "Breast_Cancer_Cells",#cluster24
             "T_cells",#cluster25
             "Breast_Cancer_Cells",#cluster26
             "Breast_Cancer_Cells",#cluster27
             "T_cells"#cluster28
)
seuratidens=mapvalues(Idents(merged_seurat), from = levels(Idents(merged_seurat)), to = ann.ids)
Idents(merged_seurat)=seuratidens
merged_seurat$cellType=Idents(merged_seurat)
DimPlot(merged_seurat, reduction = "umap", label = T, label.size = 3.5,pt.size = 2.5)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")

saveRDS(merged_seurat,"adata_annotated.rds")
rds <- readRDS("adata_annotated.rds")
#if(!require(sceasy))devtools::install_github("cellgeni/sceasy")

library(sceasy)


#if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
#remotes::install_github("mojaveazure/seurat-disk",force = TRUE)
#if(!reticulate::py_module_available('anndata'))reticulate::py_install('anndata')
#if(!reticulate::py_module_available('scanpy'))reticulate::py_install('scanpy')
dir.create('seurat2scanpy')
merged_seurat[["RNA"]] <- as(merged_seurat[["RNA"]], Class = "Assay")
sceasy::convertFormat(merged_seurat, from="seurat", to="anndata",
                      outFile='rds.h5ad')
library(stringr)
phe = merged_seurat@meta.data
table(phe$orig.ident)
