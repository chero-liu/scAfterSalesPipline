suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}

option_list = list(
    make_option(c("-i", "--input"), type="character",
          help="The input exprssion matrix in RDS format."),
    make_option(c("-g", "--genelist"),type="character",
           help="The genelist used to generate the heatmap with header."),
    make_option(c("-c", "--collapseby"),type="character",
           help="The variable level of heatmap. e.g. clusters or sampleid or celltype"),
    make_option(c("-r", "--rowcluster"),type="logical", default = "F",
           help="boolean values determining if rows should be clustered or hclust object,"),
    make_option(c("-l", "--colcluster"),type="logical", default = "F",
           help="boolean values determining if columns should be clustered or hclust object"),
    make_option(c("-n", "--showname"),type="logical", default = NULL,
           help="Whether to display the row name"),
    make_option(c("-w", "--gaps_row"),type = "character", default = NULL,
           help="Whether to display gaps by row"),
    make_option(c("-L", "--gaps_col"),type = "character", default = NULL,
           help="Whether to display gaps by col"),
    make_option( c("--outdir","-o"),type="character", default = "./",
        help="the output directory of Clustering results." ),
    make_option( c("--assay","-a"),type="character", default = "RNA",
        help="the assay to use in case of multimodal data." ),
    make_option( c("--topn"), type="integer", default = 25,
                 help = "[OPTIONAL] the number of top DEGs to visualizse."),
    make_option( c("--topby"), type = "character", default = "avg_logFC",
                 help="the column used to pick top n marker gene to visulize.The
                 option can be one of the column in the input marker genes table."),
    make_option( c("--order"), type = "character", default = NULL,
                help = "[OPTIONAL] specify the order of colnames to present from left to right.
                        OR can also be used to show subset of levels for factors specifyed by --collapedby"),
    make_option( c("--var2use", "-q" ), type = "character", default = NULL,
                help = "[OPTIONAL]The column name in cell metadata used as identity
                        of each cell combined with levels4var."),
    make_option( c("--levels4var", "-u" ), type = "character", default = NULL,
                help = "[OPTIONAL] subset of factor levels for the specified factor by --var2use."),
    make_option(c("-x", "--rowanno"), type="character", default = NULL,
                help="A table with row annotation."),
    make_option(c("-y", "--colanno"), type="character", default = NULL,
                help="A table with col annotation."),
    make_option(c( "--palette"), type="character", default = NULL,
                help="默认色板"),
    make_option( c("--use_color_anno" ), type = "logical",  default = TRUE,
                help = "[Optional]是否采用rds中注释的颜色信息，默认采用，若无则自动重新注释颜色。"),
    make_option( c("--color_file" ), type = "character",  default = NULL,
                help = "[Optional]选填，输入以tab分隔的文件，第一列的列名为metadata列名，第一列为该列元素，第二列为对应的颜色."),
    make_option( c("--subnew_celltype"),type = "character", default = "all"),
    make_option( c("--subsampleid"),type = "character", default = "all"),
    make_option( c("--subgroup"),type = "character", default = "all"),
    make_option( c("--subclusters"),type = "character", default = "all"),
    make_option( c("--prefix"),type = "character", default = "prefix")

    # make_option(c("-s", "--sign"), type="logical", default = "F",
                # help="Boolean values determine whether or not to add significance markers.")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#=================================================================================
# parse the command line parameters
#=================================================================================
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

if ( is.null(opt$assay) ){
    assay = "RNA"
} else {
    assay = opt$assay
}

if ( is.null( opt$topn )){
    topn = 25
} else {
    topn = opt$topn
}

if ( !is.null( opt$gaps_row )){
    gaps_row = as.numeric(unlist(strsplit( opt$gaps_row,",",perl = T)))
} else {
    gaps_row = opt$gaps_row
}

if ( !is.null( opt$gaps_col )){
    gaps_col = as.numeric(unlist(strsplit( opt$gaps_col,",",perl = T)))
} else {
    gaps_col = opt$gaps_col
}

#####################################################################
#read in the single cell expression matrix and the genelist
#####################################################################
seurat_ob = readRDS(opt$input)

if ( seurat_ob@version < 3 ){
    seurat_ob = UpdateSeuratObject(seurat_ob)
}

source("/public/scRNA_works/pipeline/scRNA-seq_further_analysis/Get_colors.R")
if ( is.null(opt$palette ) ){
    print("没有指定色板，将采用rds中注释的颜色或者默认色板.")
    palette = "customecol2"
}else{
    palette = opt$palette
}
split.bys = c("group","sampleid","new_celltype")[which(c("group","sampleid","new_celltype") %in% colnames(seurat_ob@meta.data))]
## 设置group为factor，指定顺序
if (class(seurat_ob@meta.data$group) == "character") {
    seurat_ob@meta.data$group = factor(seurat_ob@meta.data$group, levels = unique(seurat_ob@meta.data$group))
}
new_celltype_pal_list <- list()
for (split.by in split.bys){
    if ( ! opt$use_color_anno ){
        seurat_ob@meta.data = seurat_ob@meta.data[ ,!grepl(paste0("^",split.by,"_col$" ), colnames(seurat_ob@meta.data))]
    }
    if ( !is.null(opt$color_file)){
        color_file = read.delim(opt$color_file, sep="\t", header = T)
        meta_anno = color_anno(seurat_ob@meta.data, color_file)
    } else {
        meta_anno = seurat_ob@meta.data
    }
    color_use = get_colors(meta_anno, split.by, palette)
    # seurat_ob = Seurat::AddMetaData( seurat_ob, metadata = color_use[["object_meta"]])
    # user_color_pal = color_use[["user_color_pal"]]
    new_celltype_pal = color_use[["new_celltype_pal"]]
    new_celltype_pal = na.omit(new_celltype_pal)
    new_celltype_pal_list[[split.by]] <- new_celltype_pal
}



#get the subset of cells used for visualization if necessay
if ( !is.null(opt$levels4var)){
    if ( is.null(opt$var2use ) ){
        print("NO cell identity column name AVAILABLE! The clusters column will be used as default.")
        ident2use = "clusters"
    }else{
        ident2use = opt$var2use
    }
    cluster_list = unlist(strsplit( opt$levels4var,",",perl = T))
    seurat_ob = SubsetData( seurat_ob, subset.name = ident2use, accept.value = cluster_list)
    seurat_ob@meta.data[[ident2use]]=factor(seurat_ob@meta.data[[ident2use]],levels = sort(unique(seurat_ob@meta.data[[ident2use]])))
}


if (opt$subnew_celltype != "all"){
    subnew_celltype = str_split(opt$subnew_celltype,",")[[1]]
    seurat_ob = seurat_ob[,seurat_ob@meta.data$new_celltype %in% subnew_celltype]
    seurat_ob@meta.data$new_celltype = factor(as.character(seurat_ob@meta.data$new_celltype), levels = subnew_celltype)
    print(unique(seurat_ob@meta.data$new_celltype))
}
if (opt$subsampleid != "all"){
    subsampleid = str_split(opt$subsampleid,",")[[1]]
    seurat_ob = seurat_ob[,seurat_ob@meta.data$sampleid %in% subsampleid]
    seurat_ob@meta.data$sampleid = factor(as.character(seurat_ob@meta.data$sampleid), levels = subsampleid)
    print(unique(seurat_ob@meta.data$sampleid))
}
if (opt$subgroup != "all"){
    subgroup = str_split(opt$subgroup,",")[[1]]
    seurat_ob = seurat_ob[,seurat_ob@meta.data$group %in% subgroup]
    seurat_ob@meta.data$group = factor(as.character(seurat_ob@meta.data$group), levels = subgroup)
    print(unique(seurat_ob@meta.data$group))
}
if (opt$subclusters != "all"){
    subclusters = str_split(opt$subclusters,",")[[1]]
    seurat_ob = seurat_ob[,seurat_ob@meta.data$clusters %in% subclusters]
    seurat_ob@meta.data$clusters = factor(as.character(seurat_ob@meta.data$clusters), levels = subclusters)

    print(unique(seurat_ob@meta.data$clusters))
}


DefaultAssay(seurat_ob) = assay

genelist = read.delim(opt$genelist, sep="\t", header = T)
filtered_gene = data.frame()
if (dim(genelist)[2]>1){
    if ( "cluster" %in% colnames(genelist) ) {
        markers2vis = genelist
        markers2vis[["cluster"]] = factor(markers2vis[["cluster"]] , levels = sort(unique(seurat_ob@meta.data[[opt$collapseby]])))
        topn_markers  = markers2vis %>% group_by(cluster) %>%
            arrange(p_val,desc(avg_logFC),desc(gene_diff)) %>%
            top_n(opt$topn,.data[[opt$topby]]) %>% arrange(cluster)
        write.table( as.data.frame(topn_markers), file.path(output_dir,paste0("order_",basename(opt$genelist)) ),quote = F,row.names=F,sep="\t")
        topn_markers = topn_markers %>% mutate(folder_suffix = paste0("cluster",cluster)) %>% select(cluster,gene,folder_suffix)
        markers2vis = as.vector(topn_markers$gene)
    }else if ( "gaps_row" %in% colnames(genelist) ) {
        # gaps_row = NULL
        markers2vis = CaseMatch(search = as.vector(as.vector(genelist[,1])),match = rownames(seurat_ob))
        # new_genelist = genelist
        filtered_gene = genelist[! genelist[,1] %in% names(markers2vis ),1]
        if(length(filtered_gene)!=0){
            filtered_gene = as.data.frame(filtered_gene)
            colnames(filtered_gene) = "Gene"
            write.table(filtered_gene,file.path(output_dir,"filtered_gene.xls"),quote = F,row.names=F)
            print("There are some mismatched gene symbol, Please check filtered_gene.xls for the genename.")
        }
        # new_genelist$gaps_row = factor(new_genelist$gaps_row,levels = unique(new_genelist$gaps_row))
        # row_num = c(as.data.frame(table(new_genelist$gaps_row))[,2])
        # gap_row_ind<-c() 
        # m<-0
        # for (i in row_num){
        #   gap_row_ind <- c(gap_row_ind,(i+m))
        #   m <- (i+m)
        # }
        # #markers2vis = as.vector(new_genelist[,1])
        # gaps_row = gap_row_ind[1:(length(gap_row_ind)-1)]
        # print(gaps_row)
    }else{
        up = filter(genelist,FoldChange > 1) %>% arrange(desc(log2FoldChange ))  %>% top_n(topn,log2FoldChange ) %>% select(gene)
        down = filter(genelist,FoldChange < 1) %>% arrange(log2FoldChange )  %>% top_n(as.numeric(paste0("-",topn)),log2FoldChange ) %>% select(gene)
        genelist = rbind(up, down)
        genelist[,1]=factor(as.character(genelist[,1]))
        markers2vis = genelist[,1]
    }
}else{
    markers2vis = CaseMatch(search = as.vector(as.vector(genelist[,1])),match = rownames(seurat_ob))
    filtered_gene = genelist[! genelist[,1] %in% names(markers2vis ),1]
    if(length(filtered_gene)!=0){
        gaps_row_new = list()
        for (i in gaps_row){
            gap_temp=i
            for (gene in filtered_gene){
                if ( which(genelist[,1] == gene) > i){
                    gap_temp = gap_temp
                }else{
                    gap_temp = gap_temp - 1
                }
            }
            gaps_row_new = c(gaps_row_new,gap_temp)
        }
        gaps_row_new = unlist(gaps_row_new)
        filtered_gene = as.data.frame(filtered_gene)
        colnames(filtered_gene) = "Gene"
        write.table(filtered_gene,file.path(output_dir,"filtered_gene.xls"),quote = F,row.names=F)
        print("There are some mismatched gene symbol, Please check filtered_gene.xls for the genename.")
    }else{
        gaps_row_new = gaps_row
    }
}


if (length(markers2vis) <= 1) stop("The number of matching gene is smaller than 2, please check the input gene list.")


# count = as.matrix(seurat_ob@data[markers2vis,])
count = as.matrix(GetAssayData(seurat_ob, slot = "data")[markers2vis,])
meta.data = seurat_ob@meta.data
collapseby = opt$collapseby

meta.data$id = rownames(meta.data)
collapsed_count = vector()
if ( !collapseby %in% colnames(meta.data) ){
    stop("NO specified column found!")
}

collapsed_group = meta.data %>% group_by(.dots = collapseby) %>% do(data.frame(cellid = paste0(.$id,collapse = ",")))
if (collapseby == "clusters")  collapsed_group$clusters = paste("cluster",collapsed_group$clusters,sep="_")

for ( cells in collapsed_group$cellid ){
    samplex = unlist(strsplit(cells, ",", perl =T))
    collapsed_count= cbind(collapsed_count,rowMeans( count[,samplex,drop=F] ))
}
collapsed_count = as.matrix( collapsed_count )
collapsed_group = as.data.frame(collapsed_group)
colnames(collapsed_count) = as.matrix(collapsed_group[,1])
if (!is.null(opt$order)){
    order = unlist(strsplit( opt$order,",",perl = T))
    if (collapseby == "clusters") order = paste("cluster",order,sep="_")
    collapsed_count=collapsed_count[,order]
}

data = tibble::rownames_to_column(as.data.frame(collapsed_count),"GeneID")
write.table(data, file.path(output_dir,"heatmap_count.xls"),quote = F, row.names = F, sep = "\t")
if (dim(collapsed_count)[2]>2) {
    data_scaled = tibble::rownames_to_column(as.data.frame(pheatmap:::scale_rows(log2(collapsed_count+0.0001))),"GeneID")
    write.table(data_scaled, file.path(output_dir,"heatmap_count_scaled.xls"),quote = F, row.names = F, sep = "\t")
}
#=================================================================================
# data transforamtion & heatmap plotting 
#=================================================================================
ind <- apply(collapsed_count, 1, mean) > 0
collapsed_count_filter <- collapsed_count[ind, ]

heatmap_colors <- colorRampPalette(c("#406AA8", "white", "#D91216"))(n=299)
# palette <-colorRampPalette(c("navy", "white", "firebrick3"))(299)

annotation_colors_row=NA
annotation_colors_col=NA
annotation_colors = NA
if ( is.null(opt$rowanno) ){
    annotation_row = NA
    annotation_colors_row = NA
} else {
    annotation_row = read.delim(opt$rowanno, row.names=1)
    annotation_row[,1] = as.vector(annotation_row[,1])
    annotation_colors_row = CustomCol2(1:length(unique(annotation_row[,1])))
    names(annotation_colors_row) = unique(annotation_row[,1])
    annotation_colors_row = list(annotation_colors_row)
    names(annotation_colors_row) = colnames(annotation_row)
    # annotation_colors = annotation_colors_row
    if( names(annotation_colors_row) %in% names(new_celltype_pal_list) ) {
        col_use = names(annotation_colors_row)
        annotation_colors_row[[col_use]] = new_celltype_pal_list[[col_use]][match( names(annotation_colors_row[[col_use]]),names(new_celltype_pal_list[[col_use]]))]
    }
}

if( is.null(opt$colanno) ){
    annotation_col = NA
    annotation_colors_col = NA
}else{
    annotation_col = read.delim(opt$colanno, row.names=1)
    annotation_col[,1] = factor(annotation_col[,1], levels= unique(annotation_col[,1]))
    annotation_colors_col = CustomCol2(1:length(unique(annotation_col[,1])))
    names(annotation_colors_col) = unique(annotation_col[,1])
    annotation_colors_col = list(annotation_colors_col)
    names(annotation_colors_col) = colnames(annotation_col)
    if( names(annotation_colors_col) %in% names(new_celltype_pal_list) ) {
        col_use = names(annotation_colors_col)
        annotation_colors_col[[col_use]] = new_celltype_pal_list[[col_use]][match( names(annotation_colors_col[[col_use]]),names(new_celltype_pal_list[[col_use]]))]
    }
}

annotation_colors = ifelse( is.null(opt$rowanno) ,annotation_colors_col , c( annotation_colors_row, annotation_colors_col))
names(annotation_colors) = c(names(annotation_colors_row), names(annotation_colors_col))
# print(annotation_colors)
# if ( opt$sign && opt$rowanno == "up_down" ){
  # pmt <- exps$pvalue  #提取p值
  # #判断显著性
  # if (!is.null(pmt)){
    # ssmt <- pmt< 0.01
    # pmt[ssmt] <-'**'
    # smt <- pmt >0.01& pmt <0.05
    # pmt[smt] <- '*'
    # pmt[!ssmt&!smt]<- ''
  # } else {
    # pmt <- F
  # }
  # #构建显著性注释矩阵
  # display_number <- matrix(0,nrow(collapsed_count_filter),ncol(collapsed_count_filter))
  # for(i in 1:dim(collapsed_count_filter)[1]){
    # display_number[i,] <- pmt[i]
  # }
  
# }else{
  # display_number <- opt$sign
# }
if ( "gaps_row" %in% colnames(genelist) ) {
    df = collapsed_count
    df = df[which(rowSums(df) > 0),]
    new_genelist = genelist[ toupper(genelist[,1]) %in% toupper(rownames(df )),]
    new_genelist$gaps_row = factor(new_genelist$gaps_row,levels = unique(new_genelist$gaps_row))
    row_num = c(as.data.frame(table(new_genelist$gaps_row))[,2])
    gap_row_ind<-c() 
    m<-0
    for (i in row_num){
        gap_row_ind <- c(gap_row_ind,(i+m))
        m <- (i+m)
    }
    gaps_row_new = gap_row_ind[1:(length(gap_row_ind)-1)]
    print(gaps_row_new)
}
## parsing show rowname
if (is.null(opt$showname)) {
    showname=ifelse(dim(collapsed_count_filter)[1]>100,FALSE,TRUE) 
} else {
    showname = opt$showname}
cellwidth=36
cellheight=ifelse(showname,12,8) # if showname =F, set the panel height to 576/72point = 8 inches
p = pheatmap(log2(collapsed_count_filter+0.0001),
        color=heatmap_colors,
        cex=1,
        border=F,
        angle_col=45,
        treeheight_row=36, treeheight_col=36,
        lwd=1.1,
        cellheight=cellheight, cellwidth=cellwidth,
        scale=ifelse(dim(collapsed_count_filter)[2]==2,"none","row"),
        show_rownames=showname,
        gaps_row = gaps_row_new,
        gaps_col = gaps_col,
        annotation_col = annotation_col,
        annotation_row = annotation_row,
        annotation_colors = annotation_colors,
        # display_numbers = display_number,
        cluster_rows=opt$rowcluster ,
        cluster_cols=opt$colcluster)
ggsave(file.path(output_dir,"heatmap.pdf"), plot= p, 
    width=ifelse((is.null(opt$rowanno) & is.null(opt$colanno)),((36*dim(collapsed_count_filter)[2]+144)/50), ((36*dim(collapsed_count_filter)[2]+144)/50)+3), 
    height = ifelse(showname,(12*dim(collapsed_count_filter)[1]+108)/50,(8*dim(collapsed_count_filter)[1]+108)/50),bg="white",limitsize = FALSE)
ggsave(file.path(output_dir,"heatmap.png"), plot= p, dpi=1000,
    width=ifelse((is.null(opt$rowanno) & is.null(opt$colanno)),((36*dim(collapsed_count_filter)[2]+144)/50), ((36*dim(collapsed_count_filter)[2]+144)/50)+3), 
    height = ifelse(showname,(12*dim(collapsed_count_filter)[1]+108)/50,(8*dim(collapsed_count_filter)[1]+108)/50),bg="white",limitsize = FALSE)

