suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(OESingleCell))
suppressPackageStartupMessages( library(dplyr))
suppressPackageStartupMessages( library(optparse))
suppressPackageStartupMessages( library(ggplot2))
suppressPackageStartupMessages(library("stringr"))
#=command line parameters setting=============================
option_list = list(
    make_option( c("--input", "-i" ), type = "character",
         help = "[REQUIRED]the seurat object saved as RDS format."),
    make_option( c("--reduct", "-r" ), type = "character", default = "tsne",
         help = "[OPTIONAL]the previous calculated reduction result used in the Dimplot."),
    make_option( c("--pointsize ", "-s" ), type = "double", default = 0.5,
         help = "[OPTIONAL]the point size in the plot."),
    make_option( c("--output","-o"),type="character", default = "./",
        help="the output directory of results." ),
    make_option( c("--palette" ), type = "character",  default = NULL,
                help = "[Optional]选填，根据需求指定 Get_colors.R 中的离散型色板名."),
    make_option( c("--cellcycle_color" ), type = "character",  default = NULL,
                help = "[Optional]选填，逗号分隔，用十六进制颜色码指定每个细胞阶段的颜色"),
    make_option( c("--subnew_celltype"),type = "character", default = "all"),
    make_option( c("--subsampleid"),type = "character", default = "all"),
    make_option( c("--subgroup"),type = "character", default = "all"),
    make_option( c("--subclusters"),type = "character", default = "all"),
    make_option( c("--prefix"),type = "character", default = "prefix")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
source("/public/scRNA_works/pipeline/scRNA-seq_further_analysis/Get_colors.R")
#=================================================================================
#parse the command line parameters
#=================================================================================
if ( is.null(opt$input) ){
    stop("the seurat object is NOT AVAILABLE!")
}else{
    seurat_ob = readRDSMC(opt$input, cores = 10)
    if ( seurat_ob@version < 3){
        seurat_ob = UpdateSeuratObject(seurat_ob) #make sure the seurat object match with the latest seurat package
    }
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

if ( is.null(opt$palette ) ){
    print("没有指定色板，将采用rds中注释的颜色或者默认色板.")
    palette = "customecol2"
}else{
    palette = opt$palette
}

if ( is.null(opt$output) ){
    print("NO output directory specified,the current directory will be used!")
    output_dir = getwd()
}else{
    if ( file.exists(opt$output) ){
        output_dir = opt$output
    }else{
        output_dir = opt$output
        dir.create(output_dir, recursive = T)
    }
}

if ( is.null(opt$reduct) ){
    reduct = "tsne"
}else{
    reduct = opt$reduct
}

if ( is.null(opt$pointsize) ){
    pointsize = 0.5
}else{
    pointsize = opt$pointsize
}

output_dir = normalizePath(output_dir)

#=================================================================================
# the main workflow
#=================================================================================
metadata = seurat_ob@meta.data
if ( !"Phase" %in% colnames(seurat_ob@meta.data) ){
	genes.inuse = rownames(GetAssayData(seurat_ob, slot="counts"))
	s.genes = CaseMatch(search = cc.genes$s.genes, match = genes.inuse)
	g2m.genes = CaseMatch(search = cc.genes$g2m.genes, match = genes.inuse)
	seurat_ob <- CellCycleScoring(object = seurat_ob, s.features = s.genes, g2m.features = g2m.genes,set.ident = F) 
}

#get colors
if ( !is.null(opt$cellcycle_color)){
        color_list = unlist(strsplit( opt$cellcycle_color,",",perl = T))
        color_file = data.frame(Phase=c("G1","G2M","S"),Phase_col=color_list)
        meta_anno = color_anno(seurat_ob@meta.data, color_file)
    } else {
        meta_anno = seurat_ob@meta.data
    }

color_use = get_colors(meta_anno, "Phase", palette)
seurat_ob = AddMetaData( seurat_ob, metadata = color_use[["object_meta"]])
# user_color_pal = color_use[["user_color_pal"]]
new_celltype_pal = color_use[["new_celltype_pal"]]
new_celltype_pal = na.omit(new_celltype_pal)
print("color.use :")
print(new_celltype_pal)
print(table(seurat_ob@meta.data[, "Phase_col" ]))

seurat_ob = SetIdent(seurat_ob, value = "Phase")
nlevel = length(unique(seurat_ob@meta.data[,"Phase"]))
ggtsne = DimPlot(object = seurat_ob, reduction = reduct , pt.size = pointsize ) + theme( plot.title = element_text(hjust = 0.5)) +scale_colour_manual( values = new_celltype_pal)
ggsave(file.path(output_dir,"Cellcycle.pdf"),bg="white")
ggsave(file.path(output_dir,"Cellcycle.png"), dpi = 1000 ,limitsize = F,bg="white")

cellcycle_results = seurat_ob@meta.data %>% dplyr::rename( "Barcode" = "rawbc") %>% select( Barcode, sampleid, clusters, group, S.Score, G2M.Score, Phase,Phase_col)
write.table(cellcycle_results, quote = F,  file.path(output_dir,"Cellcycle_results.xls"), sep ="\t",row.names =F)
