source("/home/liuchenglong/script/lclFunc.r")
source("/home/liuchenglong/script/Get_colors.R")
library(dplyr)
library(ggplot2)
library("argparse")

parser <- ArgumentParser('Manualanno for Seurat object')
parser$add_argument('--manualfile', default = "auto", help='Path to manual annotation file')
parser$add_argument('--barcodefile', default = "auto", help='Path to barcode file')
parser$add_argument('--input', help='Path to Seurat input file')
parser$add_argument('--outdir', help='Output directory')
parser$add_argument('--palette', default = "customecol2", help='Color palette to use')
parser$add_argument('--reduct', default = "umap", help='Reduction method to use')
parser$add_argument('--orderbyFreq', default = "yes", help='Order by frequency')
args <- parser$parse_args()

manualfile <- args$manualfile
barcodefile <- args$barcodefile
input_path <- args$input
outdir_path <- args$outdir
palette <- args$palette
reduct <- args$reduct
orderbyFreq <- args$orderbyFreq
input <- input_path
outdir <- outdir_path
createDir(outdir)
rds = readH5seurat(input)

if (manualfile != "auto" && barcodefile != "auto"){
    manualfile <- read.csv(manualfile, sep="\t")
    barcodefile <- read.csv(barcodefile)
    if(colnames(manualfile)[2] != colnames(barcodefile)[2]){
        stop("Please provide the same column name for manualfile and barcodefile")
    }
    useCol = colnames(manualfile)[2]

    refs = unique(rds@meta.data[,colnames(manualfile)[1]])
    rds@meta.data[,colnames(manualfile)[2]] = ""
    for(i in 1:length(refs)){
        ref = refs[i]
        rds@meta.data[,colnames(manualfile)[2]][which(rds@meta.data[,colnames(manualfile)[1]] == ref)] = manualfile[,2][which(manualfile[,1] == ref)]
    }

    for (rb in rds@meta.data$rawbc){
        if (rb %in% barcodefile[,1]){
            rds@meta.data[,colnames(barcodefile)[2]][which(rds@meta.data$rawbc == rb)] = barcodefile[,2][which(barcodefile[,1] == rb)]
        }
    }
}else if (manualfile != "auto" && barcodefile == "auto") {
    manualfile <- read.csv(manualfile, sep="\t")
    useCol = colnames(manualfile)[2]
    refs = unique(rds@meta.data[,colnames(manualfile)[1]])
    rds@meta.data[,colnames(manualfile)[2]] = ""
    for(i in 1:length(refs)){
        ref = refs[i]
        rds@meta.data[,colnames(manualfile)[2]][which(rds@meta.data[,colnames(manualfile)[1]] == ref)] = manualfile[,2][which(manualfile[,1] == ref)]
    }
}else if (manualfile == "auto" && barcodefile != "auto") {
    barcodefile <- read.csv(barcodefile)
    useCol = colnames(barcodefile)[2]
    rds@meta.data[,colnames(barcodefile)[2]] = ""
    for (rb in rds@meta.data$rawbc){
        if (rb %in% barcodefile[,1]){
            rds@meta.data[,colnames(barcodefile)[2]][which(rds@meta.data$rawbc == rb)] = barcodefile[,2][which(barcodefile[,1] == rb)]
        }
    }
}else {
    stop("Please provide manualfile or barcodefile")
}


rds=rds[,rds@meta.data[,useCol] != ""]
rds=rds[,rds@meta.data[,useCol] != "Delete"]
rds=rds[,rds@meta.data[,useCol] != "delete"]
if(paste0(useCol,"_col") %in% colnames(rds@meta.data)){
    rds@meta.data <- rds@meta.data[ , !(names(rds@meta.data) %in% c(paste0(useCol,"_col")))]
}

getManualannoResult(rds,outdir,useCol="new_celltype",addFreqCount="no", saveRDS = 'yes',label = FALSE, legend = 'yes',orderbyFreq = orderbyFreq,reduct = reduct,palette=palette,pointsize=0.5)
