# 2024_0705
# jhyu
suppressPackageStartupMessages(library(loupeR, lib='/data/software/conda_envs/scrna_envs/loupeR'))
suppressPackageStartupMessages( library("optparse") )
suppressPackageStartupMessages( library("Seurat") )
suppressPackageStartupMessages(library("dplyr"))
option_list = list(
    make_option( c("--input", "-i"), type = "character", default = "TRUE",
                 help = "the seurat object or metadata object."),
    make_option( c("--format", "-f"), type = "character", default = "seurat",
                 help = "the input object format, rds or metadata (csv)."),
    make_option( c("--outdir","-o"),type="character", default = "./",
                help="the outdir directory of results.", metavar="character"),
    make_option( c("--project","-p"),type="character", default = NULL,
                help="the outdir directory of results.", metavar="character"),
    make_option( c("--groupby", "-g"), type = "character", default = NULL,
                help = "[OPTIONAL]visualize cells in seperate plot split by this groupping variable."));
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (opt$format == "seurat"){
    seurat_obj = OESingleCell::ReadX(input = opt$input, informat = 'h5seurat', verbose = F)
}else{
    seurat_obj = readRDS(opt$input)
}

# 兼容老版本脚本未存储rawbc情况
if (!"rawbc" %in% colnames(seurat_obj@meta.data)) { seurat_obj[["rawbc"]] <- seurat_obj[["orig.ident"]] }

if ( is.null(opt$outdir) ){
    output_dir = getwd()
} else {
    if ( file.exists(opt$outdir) ){
        output_dir = opt$outdir
    } else {
        output_dir = opt$outdir
        dir.create(output_dir,recursive=T)
    }
}

# 正式运行

if (!is.null(opt$project)){
    if (opt$project == "Mobi" || opt$project == "BD") {
        # barcode长度不一样，墨卓（20个碱基）10x（16个碱基）
        # 用该表替换 RenameCells
        data <- read.table(gzfile("/data/software/cellranger/cellranger-7.0.1/lib/python/cellranger/barcodes/3M-february-2018.txt.gz"), header = F)
        # 获取 rawbc
        seurat_barcodes <- seurat_obj@meta.data$rawbc
        
        # 确保抽取数量一致
        if (length(seurat_barcodes) != dim(seurat_obj)[2]) {
            stop("Seurat object does not have the correct number of rows.")
        }
        
        # 随机抽取
        random_barcodes <- sample(data$V1, dim(seurat_obj)[2], replace = FALSE)
        
        if (opt$project == "Mobi") {
            print("执行：Mobi")
            # 提取 seurat_obj 行名的后缀部分（- 及其之后的部分）
            suffixes <- sub("^[^-]+", "", seurat_barcodes)
            # 替换 -1 之前的部分
            final_barcodes <- paste0(random_barcodes, suffixes)
        } else if (opt$project == "BD") {
            print("BD")
            # barcode长度不一样，BD（22个碱基）10x（16个碱基）
            # 获取样本数
            unique_sampleids <- unique(seurat_obj@meta.data$sampleid)
            
            # 创建从 sampleid 到 suffix 的映射
            suffix_mapping <- setNames(paste0("-", seq_along(unique_sampleids)), unique_sampleids)
            
            # 提取样本 ID 前缀
            sampleid_prefix <- sub("_.*", "", seurat_barcodes)
            
            # 用随机条形码替换中间序列
            new_barcodes <- random_barcodes
            
            # 使用映射根据 sampleid_prefix 添加后缀
            suffix <- suffix_mapping[seurat_obj@meta.data$sampleid]
            final_barcodes <- paste0(new_barcodes, suffix)
        }
        # 更新 Seurat 对象的行名
        seurat_obj@meta.data$new_rawbc <- final_barcodes
    }
    simplified_meta = seurat_obj@meta.data %>%
                                dplyr::select(rawbc,new_rawbc)
    write.table(simplified_meta, quote = F,sep =",",row.names = F,
                file.path(output_dir,paste0("rawbc_metadata.csv",collapse = "")))
    # 修改Cells信息
    cc = paste0(seurat_obj@meta.data$new_rawbc)
    seurat_obj <- RenameCells(seurat_obj,new.names=cc)
}else{
    cc = paste0(seurat_obj@meta.data$rawbc)
    seurat_obj <- RenameCells(seurat_obj,new.names=cc)
}
if (!is.null(opt$groupby)){
        seurat_obj@meta.data = seurat_obj@meta.data[,c("sampleid","group","clusters",opt$groupby)]
        }else{
        seurat_obj@meta.data = seurat_obj@meta.data[,c("sampleid","group","clusters")]
    }

# 获取所有需要删除的文件名
files_to_delete <- list.files(
path = output_dir, 
pattern = "converted.cloupe", 
full.names = TRUE, 
recursive = TRUE
)

# 如果找到了需要删除的文件，则删除它们
if (length(files_to_delete) > 0) {
    file.remove(files_to_delete)
}
setup()
setwd(output_dir)
create_loupe_from_seurat(seurat_obj)
