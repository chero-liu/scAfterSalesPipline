#must source /home/liuxuan/miniconda3/bin/activate /gpfs/oe-software/conda_envs/scrna_envs/cytoTRACE
#Author :lip , liuxuan
library("optparse")
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/gpfs/oe-software/conda_envs/scrna_envs/cytoTRACE/bin/python")
Sys.setenv(PYTHONPATH = "/gpfs/oe-software/conda_envs/scrna_envs/cytoTRACE/bin/python/python3.8/site-packages")
library(CytoTRACE)
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(SeuratDisk)
source("/gpfs/oe-software/conda_envs/scrna_envs/cytoTRACE/plotCytoTRACE.R")
#=command line parameters setting=============================
option_list = list(
    make_option( c( "--input","-i"), type = "character",
        help = "the input rds file,can be multiple files"),
    make_option( c("--output","-o"),type="character", default = "./",
        help="the output directory of CytoTRACE results." ),
    make_option( c("--assay","-a"),type="character", default = "RNA",
        help="the assays info,default is RNA." ),
    make_option( c("--groupby", "-g" ), type = "character", default = "new_celltype",
        help = "column in metadat ,new_celltype ,group or clusters"),
    make_option( c("--pthread","-p"), type="integer", default = 6,
                 help = "the thread of the analysis CPU."),
    make_option( c("--reduct","-r"), type="character", default = NULL,
                 help = "[OPTIONAL]the previous calculated reduction result used in the featureplot. If this parameter is not provided, CytoTRACE will perform its own dimensionality reduction and the CytoTRACE scores will not be plotted based on the original UMAP coordinates."),
    make_option( c("--batch","-b"), type="character", default = "TRUE",
                 help = "[OPTIONAL] generates single-cell predictions of differentiation status across multiple, heterogeneous scRNA-seq batches/datasets."),
    make_option( c("--subnew_celltype"),type = "character", default = "all"),
    make_option( c("--subsampleid"),type = "character", default = "all"),
    make_option( c("--subgroup"),type = "character", default = "all"),
    make_option( c("--subclusters"),type = "character", default = "all"),
    make_option( c("--prefix"),type = "character", default = "prefix")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#####################function
CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}


#1.check file and param
if ( is.null(opt$input) ){
    stop("the seurat object is NOT AVAILABLE!")
}else{
    filetype = tools::file_ext(opt$input)
    if(filetype =="rds"){
        seurat_ob = readRDS(opt$input)
        if ( seurat_ob@version < 3 ){
            seurat_ob = UpdateSeuratObject(seurat_ob)
        }
    } else if(filetype == "h5seurat"){
        seurat_ob = SeuratDisk::LoadH5Seurat(opt$input)
    }
    Seurat::DefaultAssay(seurat_ob) <- opt$assay
}


if (opt$subnew_celltype != "all"){
    subnew_celltype = str_split(opt$subnew_celltype,",")[[1]]
    print(subnew_celltype)
    seurat_ob = seurat_ob[,seurat_ob@meta.data$new_celltype %in% subnew_celltype]
    print(unique(seurat_ob@meta.data$new_celltype))
}
if (opt$subsampleid != "all"){
    subsampleid = str_split(opt$subsampleid,",")[[1]]
    seurat_ob = seurat_ob[,seurat_ob@meta.data$sampleid %in% subsampleid]
    print(unique(seurat_ob@meta.data$sampleid))
}
if (opt$subgroup != "all"){
    subgroup = str_split(opt$subgroup,",")[[1]]
    seurat_ob = seurat_ob[,seurat_ob@meta.data$group %in% subgroup]
    print(unique(seurat_ob@meta.data$group))
}
if (opt$subclusters != "all"){
    subclusters = str_split(opt$subclusters,",")[[1]]
    seurat_ob = seurat_ob[,seurat_ob@meta.data$clusters %in% subclusters]
    print(unique(seurat_ob@meta.data$clusters))
}


if ( is.null(opt$output) ){
    print("NO output directory specified,the current directory will be used!")
    output_dir = getwd()
}else{
    if ( file.exists(opt$output)){
        output_dir = opt$output
    }else{
        output_dir = opt$output
        dir.create(output_dir, recursive = T)
    }
}
setwd(output_dir)


if ( is.null(opt$reduct) ){
    print("No embedding information provided, CytoTRACE will perform its own dimensionality reduction and the CytoTRACE scores will not be plotted based on the original UMAP coordinates.")
    reduct = NULL
}else{
    reduct = opt$reduct
}




if ( as.logical(opt$batch) ){
    print("Function iCytoTRACE will be used to generates single-cell predictions of differentiation status across multiple, heterogeneous scRNA-seq batches/datasets.")
    split <- SplitObject(seurat_ob,split.by="sampleid")
    names(split) = NULL
    for (i in 1:length(split)){
    split[[i]] = as.data.frame(split[[i]][["RNA"]]@counts)
    }

    expr_results <- iCytoTRACE(split,ncores= 8,enableFast = TRUE)

    pheno_mat = as.character(seurat_ob@meta.data[,opt$groupby])
    names(pheno_mat) = rownames(seurat_ob@meta.data)

}else{
    print("CytoTRACE function will be used.")
    expr_mat = as.data.frame(seurat_ob[[opt$assay]]@counts)
    pheno_mat = as.character(seurat_ob@meta.data[,opt$groupby])
    names(pheno_mat) = rownames(seurat_ob@meta.data)

    expr_results <- CytoTRACE(mat = expr_mat,ncores= opt$pthread,enableFast = TRUE)
}



#1.gene_barplot
plotCytoGenes(expr_results, numOfGenes = 15)
corPearson_df = tibble::rownames_to_column(as.data.frame(expr_results$cytoGenes),var = "geneID")
names(corPearson_df) = c("geneID","corPearson")
write.table(corPearson_df,file.path("corPearson_gene.xls"),sep="\t",quote=F,row.names=F,col.names=T)


if ( is.null(reduct) ){
    options(repr.plot.width = 16, repr.plot.height=6)
    plotCytoTRACE(cyto_obj = expr_results)
    plotCytoTRACE(cyto_obj = expr_results, phenotype = pheno_mat,outputDir = 'phenotype_')
}else{
    options(repr.plot.width = 16, repr.plot.height=6)
    ## 获取原始的数据的umap坐标
    umap <- as.data.frame(seurat_ob@reductions[[reduct]]@cell.embeddings)
    ## 基于原始umap坐标绘制CytoTRACE分数
    plotCytoTRACE(cyto_obj = expr_results, emb = umap)
    plotCytoTRACE(cyto_obj = expr_results, phenotype = pheno_mat, emb = umap,outputDir = 'phenotype_')
}

#删除中间文件
pdf_file_path <- file.path("Rplots.pdf")
# 检查文件是否存在
if (file.exists(pdf_file_path)) {
  # 如果存在，则删除文件
  file.remove(pdf_file_path)
  cat("Rplots.pdf 文件已删除。\n")
} else {
  cat("Rplots.pdf 文件不存在。\n")
}

if (!file.exists(file.path("CytoTRACE分析说明.docx"))) {
        file.copy(
            "/public/scRNA_works/Documents/CytoTRACE分析说明.docx",
            file.path("CytoTRACE分析说明.docx")
        )
    }

print("Convert pdf to png...")
system("ls *.pdf | xargs -I {} -P 8 bash -c '/data/software/ImageMagick/ImageMagick-v7.0.8-14/bin/convert  -verbose -density 200 -trim  ${0}  -quality 100  -flatten  ${0%.pdf}.png' {}")
