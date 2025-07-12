# module load OESingleCell/3.0.d
.libPaths("/public/scRNA_works/works/guokaiqi/software/R/x86_64-conda-linux-gnu-library/4.0")
source("/public/scRNA_works/pipeline/scRNA-seq_further_analysis/CellChat_functions.R")
source("/home/liuchenglong/script/Get_colors.R")

suppressPackageStartupMessages({
    library(Seurat)
    library(CellChat)
    library(patchwork)
    library(ggalluvial)
    library(NMF)
    library(ggplot2)
    library(optparse)
    library(ComplexHeatmap)
    library(future)
    library(tidyverse)
    library(corrplot)
    library(circlize)
    library(pheatmap)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(Cairo)
    library(stringr)
})

options(future.globals.maxSize= Inf ) # setting the maxumium mermory usage much bigger in case of big data
 
option_list = list(
    make_option( c("--input", "-i"), type = "character", default = NULL,
                help = "the seurat object saved as R object in RDS format."),
    make_option( c("--informat", "-f" ), type = "character", default = "rds",
                 help = "seurat:读取h5,rds:读取rds(jhy)."),
    make_option(c("--assay"),type="character", default = "RNA",
                help="[OPTIONAL] the assay used to calulation in case of multimodal data."),
    make_option( c("--species", "-s"), type="character", default="human",
                help="the species:human, mouse."),
    make_option( c("--column4cell", "-c"), type="character", default= NULL,
                help="the cell type annotation column in the cell annotation meta data."),
    make_option( c("--celltype", "-l"), type = "character", default = NULL,
                help = "[OPTIONAL] the comma seperated list for desired cell types,to be used with the parameter --column4cell. If NULL, all cell type will be used."),
    make_option( c("--groupby", "-g"), type="character", default=NULL,
                help="the groupping variable in the metadata for cellchat."),
    make_option( c("--contrast", "-d"), type="character", default=NULL,
                help="the string used to compare. e.g. CASE1:CON+CASE2:CON"),
    make_option( c("--output", "-o"),type="character", default = "./",
                help="the output directory of results.", metavar="character"),
    make_option( c("--threads", "-p"), type = "integer", default = "10 ",
                 help="[OPTIONAL] the threads use to run cellchat, 10 as default."),
    make_option( c("--subsetby", "-q" ), type = "character", default = NULL,
                help = "[OPTIONAL]The column name in cell metadata used as identity
                        of each cell combined with which_cells."),
    make_option( c("--which_cells", "-u" ), type = "character", default = NULL,
                help = "[OPTIONAL] subset of factor levels for the specified factor by --ident2use."),
    make_option( c("--rds" ), type = "character",  default = NULL,
                help = "[Optional]The cellchat_list.rds or cellchat_results.rds input to replot."),
    make_option( c("--use_color_anno" ), type = "logical",  default = TRUE,
                help = "[Optional]是否采用rds中注释的颜色信息，默认采用，若无则自动重新注释颜色。"),
    make_option( c("--color_file" ), type = "character",  default = NULL,
                help = "[Optional]选填，输入以tab分隔的文件，第一列的列名为metadata列名，第一列为该列元素，第二列为对应的颜色."),
    make_option( c("--palette" ), type = "character",  default = NULL,
                help = "[Optional]选填，根据需求指定 Get_colors.R 中的离散型色板名."),
    make_option( c("--split", "-y" ), type = "character",  default = NULL,
                help = "[Optional]针对cellchat输入可视化list文件，进行指定样本或者组别拆分，从而进行分析(jhy)."),
    make_option( c("--extraSignaling", "-x" ), type = "character", default = NULL,
                help = "[OPTIONAL]The extra signaling of interest to visualize specified by the user." ),
    make_option( c("--subnew_celltype"),type = "character", default = "all"),
    make_option( c("--subsampleid"),type = "character", default = "all"),
    make_option( c("--subgroup"),type = "character", default = "all"),
    make_option( c("--subclusters"),type = "character", default = "all"),
    make_option( c("--prefix"),type = "character", default = "prefix")
    );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# CustomCol2 <- function(n){
#   my_palette=c(
#     "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
#     "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
#     "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
#     "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
#     "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
#   return(my_palette[n])
# }

CustomCol2 <- function(n){
  my_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
        "#F29403", "#F781BF", "#BC9DCC", "#A65628", "#54B0E4", 
        "#222F75", "#1B9E77", "#B2DF8A", "#E3BE00", "#FB9A99", 
        "#E7298A", "#910241", "#00CDD1", "#A6CEE3", "#CE1261", 
        "#5E4FA2", "#8CA77B", "#00441B", "#DEDC00", "#B3DE69", 
        "#8DD3C7", "#999999")
    if (length(my_palette[n]) <= length(my_palette)) {
        colors <- my_palette[n]
    }
    else {
        colors <- (grDevices::colorRampPalette(my_palette))(length(my_palette[n]))
    }
    return(colors)
}


# =================================================================================
# parse the command line parameters
# =================================================================================
if ( is.null(opt$column4cell ) ){
    print("NO cell type column name AVAILABLE! The celltype annotation will be used as default.")
    ident2use = "celltype"
}else{
    ident2use = opt$column4cell
}

if ( is.null(opt$palette ) ){
    print("没有指定色板，将采用rds中注释的颜色或者默认色板.")
    palette = "cellchat"
}else{
    palette = opt$palette
}

# change the default assay for reduction if necessary
if (!is.null(opt$assay)) {
    assay = opt$assay
} else {
    assay = "RNA"
}

if (opt$species == "human") {
    # Set the ligand-receptor interaction database
    CellChatDB <- CellChatDB.human
    PPI.species <- PPI.human
} else if (opt$species == "mouse"){
    CellChatDB <- CellChatDB.mouse
    PPI.species <- PPI.mouse
}

if (!is.null(opt$groupby)) {
    groupby <- opt$groupby
} else {
    groupby <- NULL
}

if (!is.null(opt$contrast)) {
    contrasts = unlist(strsplit(opt$contrast, "\\+", perl = T))
}

if (is.null(opt$threads)) {
    threads <- 10
} else {
    threads <- opt$threads
}

if (is.null(opt$output)) {
    print("NO output directory specified,the current directory will be used!")
    output_dir <- getwd()
} else {
    if (file.exists(opt$output)) {
        output_dir <- opt$output
    } else {
        output_dir <- opt$output
        dir.create(output_dir, recursive = T)
    }
}
output_dir <- normalizePath(output_dir)
###########################################
if ( is.null(opt$input) ){
    if ( is.null(opt$rds) ){
        stop("The seurat object is NOT AVAILABLE!")
    }else{
        print("no input rds, read cellchat result rds")
    }
} else {
    if ( opt$informat == "seurat" ){
        print("读取h5")
        seurat_ob = OESingleCell::ReadX(input = opt$input, informat = 'h5seurat', verbose = F)
    }
    else if ( opt$informat == "rds" ){
        print("读取rds（默认）")
        seurat_ob = readRDS(opt$input)
    }else{
        print("请确保 -f 参数使用正确")
    }

    if ( ! opt$use_color_anno ){
        seurat_ob@meta.data = seurat_ob@meta.data[ ,!grepl(paste0("^",ident2use,"_col$" ), colnames(seurat_ob@meta.data))]
    }
    if ( !is.null(opt$color_file)){
        color_file = read.delim(opt$color_file, sep="\t", header = T)
        meta_anno = color_anno(seurat_ob@meta.data, color_file)
    } else {
        meta_anno = seurat_ob@meta.data
    }
    color_use = get_colors(meta_anno, ident2use, palette)
    seurat_ob = AddMetaData( seurat_ob, metadata = color_use[["object_meta"]])
    # user_color_pal = color_use[["user_color_pal"]]
    new_celltype_pal = color_use[["new_celltype_pal"]]
    new_celltype_pal = na.omit(new_celltype_pal)
    print("color.use :")
    print(new_celltype_pal)
    print(table(seurat_ob@meta.data[, paste0( ident2use,"_col" )]))
    if ('rawbc' %in% colnames(seurat_ob@meta.data)) {
        Barcode_content = 'rawbc'
    }else{
        Barcode_content = 'orig.ident'
    }
    simplified_meta = seurat_ob@meta.data %>%
                              dplyr::rename( "Barcode" = Barcode_content) %>%
                              dplyr::select( Barcode, sampleid, group,!!ident2use, !!paste0( ident2use,"_col" ))
    write.table(simplified_meta, quote = F,sep =",",row.names = F,
              file.path(output_dir,paste0(ident2use,"_col.metadata.csv",collapse = "")))

#get the subset of cells used for visualization if necessay
    if ( !is.null(opt$which_cells)){
        if ( is.null(opt$subsetby ) ){
            print("NO cell identity column name AVAILABLE! The clusters column will be used as default.")
            subsetby = "clusters"
        }else{
            subsetby = opt$subsetby
        }
        cluster_list = unlist(strsplit( opt$which_cells,",",perl = T))
        # seurat_ob = SubsetData( seurat_ob, subset.name = subsetby, accept.value = cluster_list)
        desired_cells= seurat_ob@meta.data[which(seurat_ob@meta.data[,subsetby] %in% cluster_list),]
        seurat_ob = seurat_ob[, rownames(desired_cells)]
    }
    if ( !is.null(opt$celltype)){
        cluster_list = unlist(strsplit( opt$celltype,",",perl = T))
        # seurat_ob = SubsetData(seurat_ob,cells =
                               # OldWhichCells( seurat_ob, subset.name= ident2use, accept.value = cluster_list)) 
        desired_cells= seurat_ob@meta.data[which(seurat_ob@meta.data[,ident2use] %in% cluster_list),]
        seurat_ob = seurat_ob[, rownames(desired_cells)]
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
# =================================================================================
# Load data
# =================================================================================
if ( is.null(opt$rds) ){
    seurat_ob = SetIdent(seurat_ob, value = ident2use)
    if (!is.null(groupby)) {
        sp <- SplitObject(seurat_ob, split.by = groupby)
    } else {
        sp <- list()
        sp[['all']] <- seurat_ob
    }
    if(class(seurat_ob@meta.data[,ident2use]) == "factor"){sp = lapply(sp, function(x){ x@meta.data[,ident2use] = droplevels(x@meta.data[,ident2use] )
                                                       return(x)}) }
    cellchat_list <- list()
    for(samples in names(sp)){
        # data.input <- GetAssayData(sp[[samples]], assay = assay, slot = "data")  # normalized data matrix
        #label <- Idents(sp[[samples]])
        # meta <- data.frame(group = label, row.names = names(label))
          # create a dataframe of the cell labels
        # cellchat_list[[samples]] <- createCellChat(object = data.input, meta = meta, group.by = "group")  # Create a CellChat object
        cellchat_list[[samples]] <- createCellChat(object = sp[[samples]], group.by = ident2use)  # Create a CellChat object
        CellChatDB.use <- CellChatDB
        cellchat_list[[samples]]@DB <- CellChatDB.use
        # subset the expression data of signaling genes for saving computation cost
        cellchat_list[[samples]] <- subsetData(cellchat_list[[samples]]) # This step is necessary even if using the whole database
        future::plan("multiprocess", workers = 10)  # do parallel
        cellchat_list[[samples]] <- identifyOverExpressedGenes(cellchat_list[[samples]])
        cellchat_list[[samples]] <- identifyOverExpressedInteractions(cellchat_list[[samples]])
        cellchat_list[[samples]] <- projectData(cellchat_list[[samples]], PPI.species)
        cellchat_list[[samples]] <- computeCommunProb(cellchat_list[[samples]])
        cellchat_list[[samples]] <- filterCommunication(cellchat_list[[samples]], min.cells = 10) # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
        cellchat_list[[samples]] <- computeCommunProbPathway(cellchat_list[[samples]]) # Infer the cell-cell communication at a signaling pathway level
        cellchat_list[[samples]] <- aggregateNet(cellchat_list[[samples]]) # Calculate the aggregated cell-cell communication network        
    }
    saveRDS(cellchat_list, file.path(output_dir, "cellchat_list.rds"))
}else{
    cellchat_list = readRDS(opt$rds)
    if(!is.null(opt$split)){
    print("开始执行特定数据的可视化展示(jhy)")
    split_list <- unlist(strsplit(opt$split, ",", perl = TRUE))
    cellchat_list <- setNames(lapply(split_list, function(name) {
    return(cellchat_list[[name]])
    }), split_list)
    }
    cellchat_list = lapply(cellchat_list, function(x){
                           x = updateCellChat(x)
                            if ( ! opt$use_color_anno ){
                                x@meta = x@meta[ ,!grepl(paste0("^",ident2use,"_col$" ), colnames(x@meta))]
                            }
                            if ( !is.null(opt$color_file)){
                                color_file = read.delim(opt$color_file, sep="\t", header = T)
                                meta_anno = color_anno(x@meta, color_file)
                            } else { meta_anno = x@meta }
                           color_use = get_colors(meta_anno, ident2use, palette)
                           x = addMeta(x, meta = color_use[["object_meta"]])
                           return(x)
                           })

    if ( !is.null(opt$which_cells)){
        if ( is.null(opt$subsetby ) ){
            print("NO cell identity column name AVAILABLE! The clusters column will be used as default.")
            subsetby = "clusters"
        }else{
            subsetby = opt$subsetby
        }
        cluster_list = unlist(strsplit( opt$which_cells,",",perl = T))
        cellchat_list = lapply(cellchat_list, function(x){
                           x = subsetCellChat_comma(x, group.by=subsetby, idents.use = cluster_list)
                           if(class(x@meta[,ident2use]) == "factor"){x@meta[,ident2use] = droplevels(x@meta[,ident2use] )}
                           return(x) })
    }
}
if( length(cellchat_list)==1 ){
    cellchat <- cellchat_list[[1]]
    color_use = get_colors(cellchat@meta, ident2use, palette)
    new_celltype_pal = color_use[["new_celltype_pal"]]
    groupSize <- as.numeric(table(cellchat@idents))
    sub_col = new_celltype_pal
    if ( length(rownames(cellchat@net$count)) != length(new_celltype_pal)){
        sub_col = subset( new_celltype_pal, names(new_celltype_pal) %in% rownames(cellchat@net$count) )
    }
    pdf(file.path(output_dir, 'interaction_number_network.pdf'))
    par(mfrow = c(1,1), xpd=T)
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, arrow.size=1, title.name = "Number of interactions", color.use = sub_col)
    dev.off()
    interaction_number_data = cbind(rownames(cellchat@net$count), cellchat@net$count)
    colnames(interaction_number_data) = c("ligand_cell/receptor_cell",colnames(interaction_number_data)[-1])
    write.table(interaction_number_data , file.path(output_dir,  "interaction_number_network.xls") , sep="\t", quote=F, row.names=F)

    pdf(file.path(output_dir, 'interaction_strength_network.pdf'))
    par(mfrow = c(1,1), xpd=T)
    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, arrow.size=1, title.name = "Interaction weights/strength", color.use = sub_col)
    dev.off()
    interaction_strength_data = cbind(rownames(cellchat@net$weight), cellchat@net$weight)
    colnames(interaction_strength_data) = c("ligand_cell/receptor_cell",colnames(interaction_strength_data)[-1])
    write.table(interaction_strength_data , file.path(output_dir, "interaction_strength_network.xls") , sep="\t", quote=F, row.names=F)

    pdf(file.path(output_dir, "interaction_number_heatmap.pdf"))
    par(mfrow = c(1, 1), xpd = TRUE)
    ht = netVisual_heatmap_c(cellchat, measure = "count", color.heatmap = "Reds", 
                        color.use = sub_col, title.name = "Number of interactions")
    ComplexHeatmap::draw(ht, ht_gap = unit(0.5, "cm"))
    dev.off()
    pdf(file.path(output_dir, "interaction_strength_heatmap.pdf"))
    par(mfrow = c(1, 1), xpd = TRUE)
    ht = netVisual_heatmap_c(cellchat, measure = "weight", color.heatmap = "Reds", 
                        color.use = sub_col, title.name = "Interaction strength")
    ComplexHeatmap::draw(ht, ht_gap = unit(0.5, "cm"))
    dev.off()

    nrow = ceiling(length(unique(new_celltype_pal))/2)
    mat <- cellchat@net$count
    pdf(file.path(output_dir, "interaction_number_network_groupby.pdf"), height = 7 * (nrow-1))
    par(mfrow = c(nrow, 2), xpd = TRUE)
    for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, label.edge= F, edge.weight.max = max(mat), title.name = rownames(mat)[i], color.use = sub_col)
    }
    dev.off()

    mat <- cellchat@net$weight
    pdf(file.path(output_dir, "interaction_strength_network_groupby.pdf"), height = 7 * (nrow-1))
    par(mfrow = c(nrow, 2), xpd = TRUE)
    for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, label.edge= F, edge.weight.max = max(mat), title.name = rownames(mat)[i], color.use = sub_col)
    }
    dev.off()

    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

    outdir <- file.path(output_dir, "signaling_pathway_visualize")
    dir.create(outdir, recursive = T)
    pathways = c()
    for(i in cellchat@netP$pathways){ if(sum(cellchat@netP$prob[,,i]) !=0 ){pathways = c(pathways, i)} }
    for(pathway in pathways){
        pathways.show <- pathway

        pdf(file.path(outdir, paste0(pathway,'_Circle_plot.pdf')))
        par(mfrow=c(1,1), xpd=TRUE)
        netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", color.use = sub_col)
        dev.off()
        # Chord diagram
        pdf(file.path(outdir, paste0(pathway,'_Chord_diagram.pdf')),width = length(sub_col)*0.2+6 , height = length(sub_col)*0.2+6)
        par(mfrow=c(1,1))
        netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", color.use = sub_col)
        dev.off()
        # Heatmap
        pdf(file.path(outdir, paste0(pathway,'_Heatmap.pdf')))
        par(mfrow=c(1,1))
        print(netVisual_heatmap_c(cellchat, signaling = pathways.show, color.heatmap = "Reds", color.use = sub_col))
        dev.off()
        # L-R_contribution
        p = netAnalysis_contribution(cellchat, signaling = pathways.show )
        ggsave(file.path(outdir, paste0(pathway,'_L-R_contribution.pdf')), plot=p,bg="white" )
        # GeneExpression
        p = plotGeneExpression(cellchat, signaling = pathways.show , color.use = sub_col)
        ggsave(file.path(outdir, paste0(pathway,'_GeneExpression.pdf')), plot=p,bg="white" ) 
        # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
        pdf(file.path(outdir, paste0(pathway,'_signalingRole.pdf')))
        p = netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, font.size = 4, color.use = sub_col)
        dev.off()
    }
    p <- netVisual_bubble(cellchat, remove.isolate = FALSE, return.data = T)
    p$communication$pval <- gsub("1", "p > 0.05", p$communication$pval)
    p$communication$pval <- gsub("2", "0.01 < p < 0.05", p$communication$pval)
    p$communication$pval <- gsub("3", "p < 0.01", p$communication$pval)
    write.table(p$communication, file.path(output_dir, "communication.xls"), sep = "\t", row.names = F, quote = F)

    pp <- netVisual_bubble(cellchat, remove.isolate = FALSE)
    ggsave(file.path(output_dir, "significant_interactions_bubble_plot.pdf"), plot = pp, height = length(unique(p$communication$interaction_name)) / 7 + 1, width = length(unique(p$communication$source.target)) / 5 + 1, limitsize = FALSE,bg="white")


    if ( !is.null(opt$extraSignaling) ){
      extra_signaling = read.delim(opt$extraSignaling, sep=",", header = T)
      if (length(extra_signaling$interaction_name) == 0){
      filtered_pairedLR = extra_signaling$interaction_name
     } else {
      filtered_pairedLR = subset(extra_signaling, interaction_name %in% p$communication$interaction_name)
     }      
      filtered_source = Reduce(intersect,list(extra_signaling$source,p$communication$source))
      filtered_target = Reduce(intersect,list(extra_signaling$target,p$communication$target))
      if (length(filtered_source) == 0 && length(filtered_target) == 0 && length(filtered_pairedLR) >0) {
        filtered_pp <- netVisual_bubble(cellchat, pairLR.use = filtered_pairedLR['interaction_name'], remove.isolate = FALSE)
          } else if (length(filtered_source) > 0 && length(filtered_target) > 0 && length(filtered_pairedLR) == 0) {
          filtered_pp <- netVisual_bubble(cellchat, sources.use = filtered_source, targets.use = filtered_target, remove.isolate = FALSE)
        } else if (length(filtered_source) > 0 && length(filtered_target) > 0 && length(filtered_pairedLR) > 0) {
          filtered_pp <- netVisual_bubble(cellchat, pairLR.use = filtered_pairedLR['interaction_name'],sources.use = filtered_source, targets.use = filtered_target, remove.isolate = FALSE)
        }
        ggsave(file.path(output_dir, "assigned_interactions_bubble_plot.pdf"), plot = filtered_pp, height = max(20,length(unique(extra_signaling$interaction_name))), width = length(unique(extra_signaling$source.target))/2 + 8, limitsize = FALSE,bg="white")
    }

    gg1 <- netAnalysis_signalingRole_scatter(cellchat, color.use = sub_col)
    ggsave(file.path(output_dir, 'Signaling_role.pdf'),plot=gg1,bg="white") 

    if (!file.exists(file.path(output_dir, "CellChat细胞通讯分析说明.docx"))) {
        file.copy(
            "/public/scRNA_works/Documents/CellChat细胞通讯分析说明.docx",
            file.path(output_dir, "CellChat细胞通讯分析说明.docx")
        )
    }
} else {

    cellchat <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))
    color_use = get_colors(cellchat@meta, ident2use, palette)
    new_celltype_pal = color_use[["new_celltype_pal"]]

    gg1 <- compareInteractions(cellchat, show.legend = F, group = names(cellchat_list), 
                               group.levels = names(cellchat_list), color.use = CustomCol2(1:length(names(cellchat_list))))
    gg2 <- compareInteractions(cellchat, show.legend = F, group = names(cellchat_list), 
                               group.levels = names(cellchat_list), color.use = CustomCol2(1:length(names(cellchat_list))), measure = "weight")
    gg <- gg1 + gg2
    ggsave(file.path(output_dir, "compare_interaction_number_strength.pdf"), plot = gg, height = 5,bg="white")
    for(contrast in contrasts){
        contrast = unlist(strsplit(contrast, ":", perl = T))
        case = contrast[1]
        control = contrast[2]

        pdf(file.path(output_dir, paste0("diff_",case,"-vs-",control,"_interaction_number_network.pdf")))
        par(mfrow = c(1, 1), xpd = T)
        netVisual_diffInteraction_c(cellchat, comparison = c(control, case), weight.scale = T, color.use = new_celltype_pal)
        dev.off()
        pdf(file.path(output_dir, paste0("diff_",case,"-vs-",control,"_interaction_strength_network.pdf")))
        par(mfrow = c(1, 1), xpd = T)
        netVisual_diffInteraction_c(cellchat, comparison = c(control, case), measure = "weight", weight.scale = T, color.use = new_celltype_pal)
        dev.off()

        pdf(file.path(output_dir, paste0("diff_",case,"-vs-",control,"_interaction_number_heatmap.pdf")))
        par(mfrow = c(1, 1), xpd = T)
        print(netVisual_heatmap_c(cellchat, comparison = c(control, case), color.use = new_celltype_pal))
        dev.off()
        pdf(file.path(output_dir, paste0("diff_",case,"-vs-",control,"_interaction_strength_heatmap.pdf")))
        par(mfrow = c(1, 1), xpd = T)
        print(netVisual_heatmap_c(cellchat, measure = "weight", comparison = c(control, case), color.use = new_celltype_pal))
        dev.off()

        case_index = which(levels(cellchat@meta$datasets) == case)
        control_index = which(levels(cellchat@meta$datasets) == control)
        gg1 <- netVisual_bubble_oe(cellchat, comparison = c(control_index, case_index), max.dataset = case_index, title.name = paste0("Increased signaling in ", case), angle.x = 45, remove.isolate = F)+ guides( size = guide_legend(order = 1 ), fill = guide_legend(order = 2 ))
        gg2 <- netVisual_bubble_oe(cellchat, comparison = c(control_index, case_index), max.dataset = control_index, title.name = paste0("Decreased signaling in ", case), angle.x = 45, remove.isolate = F)+ guides( size = guide_legend(order = 1 ), fill = guide_legend(order = 2 ))
        gg <- gg1 + gg2
        ggsave(file.path(output_dir, paste0("diff_", case, "-vs-", control, "_bubble_plot.pdf")), plot = gg, width = max(length(unique(gg1$data$group.names))/3*2.5,length(unique(gg2$data$group.names))/3*2.5), height = max(length(unique(gg1$data$interaction_name))/5+1,length(unique(gg2$data$interaction_name))/5+1), limitsize = FALSE,bg="white")
    }
    count.max <- getMaxWeight(cellchat_list, attribute = c("idents", "count"))
    weight.max <- getMaxWeight(cellchat_list, attribute = c("idents", "weight"))
    for (i in 1:length(cellchat_list)) {
        group_name = names(cellchat_list)[i]
        groupSize <- as.numeric(table(cellchat_list[[i]]@idents))
        sub_col = get_colors(cellchat_list[[i]]@meta, ident2use, palette)[["new_celltype_pal"]]
        pdf(file.path(output_dir, paste0(group_name, "_interaction_number_network.pdf")))
        par(mfrow = c(1, 1), xpd = TRUE)
        netVisual_circle(cellchat_list[[i]]@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, arrow.size = 0.5, edge.weight.max = count.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", group_name), color.use = sub_col)
        dev.off()
        pdf(file.path(output_dir, paste0(group_name, "_interaction_strength_network.pdf")))
        par(mfrow = c(1, 1), xpd = TRUE)
        netVisual_circle(cellchat_list[[i]]@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, arrow.size = 0.5, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction strength - ", group_name), color.use = sub_col)
        dev.off()

        value_max_row = signif(max(unlist(sapply(cellchat_list, function (x) {rowSums(x@net$count) }))), 1)
        value_max_col = signif(max(unlist(sapply(cellchat_list, function (x) {colSums(x@net$count) }))), 1)
        value_max = signif(max(unlist(sapply(cellchat_list, function (x) {x@net$count }))), 1)
        pdf(file.path(output_dir, paste0(group_name, "_interaction_number_heatmap.pdf")))
        par(mfrow = c(1, 1), xpd = TRUE)
        ht = netVisual_heatmap_c(cellchat_list[[i]], measure = "count", color.heatmap = "Reds", 
                            value_max = value_max, value_max_row = value_max_row, value_max_col = value_max_col,
                            color.use = sub_col, title.name = paste0("Number of interactions - ", group_name))
        ComplexHeatmap::draw(ht, ht_gap = unit(0.5, "cm"))
        dev.off()
        value_max_row = signif(max(unlist(sapply(cellchat_list, function (x) {rowSums(x@net$weight) }))), 1)
        value_max_col = signif(max(unlist(sapply(cellchat_list, function (x) {colSums(x@net$weight) }))), 1)
        value_max = signif(max(unlist(sapply(cellchat_list, function (x) {x@net$weight }))), 1)
        pdf(file.path(output_dir, paste0(group_name, "_interaction_strength_heatmap.pdf")))
        par(mfrow = c(1, 1), xpd = TRUE)
        ht = netVisual_heatmap_c(cellchat_list[[i]], measure = "weight", color.heatmap = "Reds", 
                            value_max = value_max, value_max_row = value_max_row, value_max_col = value_max_col,
                            color.use = sub_col, title.name = paste0("Interaction strength - ", group_name))
        ComplexHeatmap::draw(ht, ht_gap = unit(0.5, "cm"))
        dev.off()

        interaction_number_data <- cbind(rownames(cellchat_list[[i]]@net$count), cellchat_list[[i]]@net$count)
        colnames(interaction_number_data) <- c("ligand_cell/receptor_cell", colnames(interaction_number_data)[-1])
        write.table(interaction_number_data, file.path(output_dir, paste0(group_name, "_interaction_number_network.xls")), sep = "\t", quote = F, row.names = F, col.names = T)

        interaction_strength_data <- cbind(rownames(cellchat_list[[i]]@net$weight), cellchat_list[[i]]@net$weight)
        colnames(interaction_strength_data) <- c("ligand_cell/receptor_cell", colnames(interaction_strength_data)[-1])
        write.table(interaction_strength_data, file.path(output_dir, paste0(group_name, "_interaction_strength_network.xls")), sep = "\t", quote = F, row.names = F, col.names = T)
    }
    cellchat_list_length = sapply(cellchat_list,function(x){ return(length(levels(x@idents))) })
    do.stat = ifelse( length(unique(cellchat_list_length)) == 1, TRUE, FALSE)
    gg <- rankNet(cellchat, mode = "comparison", stacked = F, comparison = c(1:length(cellchat_list)), do.stat = do.stat, color.use = CustomCol2(1:length(cellchat_list)))
    ggsave(file.path(output_dir, "information_flow_of_each_signaling_pathway.pdf"), plot = gg, height = max(7, length(unique(gg$data$name))/10), width = 5,bg="white") # TODO

    p <- netVisual_bubble_oe(cellchat, sources.use = names(new_celltype_pal), targets.use = names(new_celltype_pal), comparison = c(1:length(cellchat_list)), angle.x = 45)
    ggsave(file.path(output_dir, "significant_interactions_bubble_plot.pdf"), plot = p, height = length(unique(p$data$interaction_name))/6, width = max(7, length(unique(p$data$group.names))/2 + 8 ), limitsize = FALSE,bg="white") 
    p <- netVisual_bubble_oe(cellchat, sources.use = names(new_celltype_pal), targets.use = names(new_celltype_pal), comparison = c(1:length(cellchat_list)), angle.x = 45, return.data = T)
    p$communication$pval <- gsub("1", "p > 0.05", p$communication$pval)
    p$communication$pval <- gsub("2", "0.01 < p < 0.05", p$communication$pval)
    p$communication$pval <- gsub("3", "p < 0.01", p$communication$pval)
    write.table(p$communication, file.path(output_dir, "communication.xls"), sep = "\t", row.names = F, quote = F)

    if ( !is.null(opt$extraSignaling) ){
      extra_signaling = read.delim(opt$extraSignaling, sep=",", header = T)
      if (length(extra_signaling$interaction_name) == 0){
        filtered_pairedLR = extra_signaling$interaction_name
       } else {
        filtered_pairedLR = subset(extra_signaling, interaction_name %in% p$communication$interaction_name)
       }   
      filtered_source = Reduce(intersect,list(extra_signaling$source,p$communication$source))
      filtered_target = Reduce(intersect,list(extra_signaling$target,p$communication$target))
      #nlevels = length(unique(p$communication$dataset))
      if (length(filtered_source) == 0 && length(filtered_target) == 0 && length(filtered_pairedLR) >0) {
          filtered_pp <- netVisual_bubble(cellchat, pairLR.use = filtered_pairedLR['interaction_name'], comparison = c(1:length(cellchat_list)), remove.isolate = FALSE)
        } else if (length(filtered_source) > 0 && length(filtered_target) > 0 && length(filtered_pairedLR) == 0) {
          filtered_pp <- netVisual_bubble(cellchat, sources.use = filtered_source, targets.use = filtered_target, comparison = c(1:length(cellchat_list)),remove.isolate = FALSE)
        } else if (length(filtered_source) > 0 && length(filtered_target) > 0 && length(filtered_pairedLR) > 0) {
          filtered_pp <- netVisual_bubble(cellchat, pairLR.use = filtered_pairedLR['interaction_name'],sources.use = filtered_source, targets.use = filtered_target,comparison = c(1:length(cellchat_list)), remove.isolate = FALSE)
        }
        ggsave(file.path(output_dir, "assigned_interactions_bubble_plot.pdf"), plot = filtered_pp, height = max(20,length(unique(extra_signaling$interaction_name))), width = length(unique(extra_signaling$source.target))/2 + 8, limitsize = FALSE,bg="white")
    }

    outdir <- file.path(output_dir, "signaling_pathway_visualize")
    dir.create(outdir, recursive = T)
    nrow = ceiling(length(cellchat_list) / 2)
    li = list()
    for(i in 1:length(cellchat_list)){li[[i]] = cellchat_list[[i]]@netP$pathways}
    pathways = Reduce(intersect, li)
    for (pathway in pathways) {
        pathways.show <- pathway
        weight.max <- getMaxWeight(cellchat_list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
        value_max_row = signif(max(unlist(sapply(cellchat_list, function (x) {rowSums(x@netP$prob[, , pathways.show]) }))), 1)
        value_max_col = signif(max(unlist(sapply(cellchat_list, function (x) {colSums(x@netP$prob[, , pathways.show]) }))), 1)
        value_max = signif(max(unlist(sapply(cellchat_list, function (x) {x@netP$prob[, , pathways.show] }))), 1)
        # Circle plot
        pdf(file.path(outdir, paste0(pathway, "_Circle_plot.pdf")), width = 10)
        par(mfcol = c(nrow, 2), xpd = TRUE)
        for (i in 1:length(cellchat_list)) {
            netVisual_aggregate(cellchat_list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, 
                                signaling.name = paste(pathways.show, names(cellchat_list)[i]), color.use = get_colors(cellchat_list[[i]]@meta, ident2use, palette)[["new_celltype_pal"]] )
        }
        dev.off()
        # Heatmap
        for (i in 1:length(cellchat_list)) {
            pdf(file.path(outdir, paste0(pathway, "_Heatmap-", names(cellchat_list)[i], ".pdf")), width = 7)
            par(mfrow = c(1, 1), xpd = TRUE)
            ht <- netVisual_heatmap_c(cellchat_list[[i]], signaling = pathways.show, color.heatmap = "Reds", 
                  value_max = value_max, value_max_row = value_max_row, value_max_col = value_max_col, 
                  title.name = paste(pathways.show, "signaling ", names(cellchat_list)[i]), color.use = get_colors(cellchat_list[[i]]@meta, ident2use, palette)[["new_celltype_pal"]] )
            ComplexHeatmap::draw(ht, ht_gap = unit(0.5, "cm"))
            dev.off()
        }
        # L-R_contribution
        for (i in 1:length(cellchat_list)) {
            p = netAnalysis_contribution(cellchat_list[[i]], signaling = pathways.show )
            ggsave(file.path(outdir, paste0(pathway,'_L-R_contribution-', names(cellchat_list)[i], '.pdf')), plot=p,bg="white" )
        }
        # Chord diagram
        pdf(file.path(outdir, paste0(pathway, "_Chord_diagram.pdf")), width = nrow*7+20 , height = nrow*7+5 )
        par(mfrow = c(nrow, 2), xpd = TRUE)
        for (i in 1:length(cellchat_list)) {
            netVisual_aggregate(cellchat_list[[i]], signaling = pathways.show, vertex.label.cex= 1.3,layout = "chord", 
                                signaling.name = paste(pathways.show, names(cellchat_list)[i]), color.use = get_colors(cellchat_list[[i]]@meta, ident2use, palette)[["new_celltype_pal"]])
        }
        dev.off()
        # Gene Expression
        cellchat@meta$datasets <- factor(cellchat@meta$datasets, levels = names(cellchat_list)) # set factor level
        p <- plotGeneExpression(cellchat, signaling = pathways.show, split.by = "datasets", colors.ggplot = T, color.use = CustomCol2(1:length(cellchat_list)))
        ggsave(file.path(outdir, paste0(pathway, "_GeneExpression.pdf")), plot = p, height = 10, width = 12,bg="white")
    }
    if (!file.exists(file.path(output_dir, "CellChat细胞通讯分析说明-多组比较.docx"))) {
        file.copy(
            "/public/scRNA_works/Documents/CellChat细胞通讯分析说明-多组比较.docx",
            file.path(output_dir, "CellChat细胞通讯分析说明-多组比较.docx")
        )
    }
}

if ( is.null(opt$rds) ){
    saveRDS(cellchat, file.path(output_dir, "cellchat_results.rds"))
}

# convert pdf to png
setwd(output_dir)
print("Convert pdf to png...")
system("for strength in `ls *_strength.pdf`; do /home/dongjiaoyang/miniconda3/envs/OESingleCell/bin/convert  -verbose -density 500 -trim $strength -quality 100 -flatten ${strength/.pdf/.png} & done")
system("for heatmap in `ls *_heatmap.pdf`; do /home/dongjiaoyang/miniconda3/envs/OESingleCell/bin/convert  -verbose -density 500 -trim $heatmap -quality 100 -flatten ${heatmap/.pdf/.png} & done")
system("for network in `ls *_network.pdf`; do /home/dongjiaoyang/miniconda3/envs/OESingleCell/bin/convert  -verbose -density 500 -trim $network -quality 100 -flatten ${network/.pdf/.png} & done")
if (is.null(opt$groupby)) {
    system("/home/dongjiaoyang/miniconda3/envs/OESingleCell/bin/convert  -verbose -density 500 -trim Signaling_role.pdf -quality 100 -flatten Signaling_role.png &")
    system("for groupby in `ls *_network_groupby.pdf`; do /home/dongjiaoyang/miniconda3/envs/OESingleCell/bin/convert  -verbose -density 500 -trim $groupby -quality 100 -flatten ${groupby/.pdf/.png} & done &")
}
system("for i in  `ls *_signaling_pathway.pdf`; do /home/dongjiaoyang/miniconda3/envs/OESingleCell/bin/convert  -verbose -density 500 -trim $i -quality 100 -flatten ${i/.pdf/.png} ; done; wait")

setwd(outdir)
system("for diagram in `ls *_diagram.pdf`; do /home/dongjiaoyang/miniconda3/envs/OESingleCell/bin/convert  -verbose -density 500 -trim $diagram -quality 100 -flatten ${diagram/.pdf/.png} & done")
system("for Circle in `ls *_Circle_plot.pdf`; do /home/dongjiaoyang/miniconda3/envs/OESingleCell/bin/convert  -verbose -density 500 -trim $Circle -quality 100 -flatten ${Circle/.pdf/.png} & done")
system("for GeneExpression in `ls *_GeneExpression.pdf`; do /home/dongjiaoyang/miniconda3/envs/OESingleCell/bin/convert  -verbose -density 500 -trim $GeneExpression -quality 100 -flatten ${GeneExpression/.pdf/.png} & done")
system("for Heatmap in  `ls *_Heatmap*.pdf`; do /home/dongjiaoyang/miniconda3/envs/OESingleCell/bin/convert  -verbose -density 500 -trim $Heatmap -quality 100 -flatten ${Heatmap/.pdf/.png} & done")
if (is.null(opt$groupby)) {
    system("for i in `ls *_signalingRole.pdf`; do /home/dongjiaoyang/miniconda3/envs/OESingleCell/bin/convert  -verbose -density 500 -trim $i -quality 100 -flatten ${i/.pdf/.png} & done")
}
system("for i in  `ls *_L-R_contribution*.pdf`; do /home/dongjiaoyang/miniconda3/envs/OESingleCell/bin/convert  -verbose -density 500 -trim $i -quality 100 -flatten ${i/.pdf/.png} ; done; wait")

