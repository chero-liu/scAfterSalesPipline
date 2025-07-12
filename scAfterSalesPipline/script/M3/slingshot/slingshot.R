#!/usr/bin/env Rscript
# Created by: guokaiqi
# Created on: 2022-05-12 20:50:19
# module load OESingleCell/3.0.d
### https://bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html#using-slingshot

.libPaths("/home/guokaiqi/R/x86_64-conda-linux-gnu-library/4.0")
source("/public/scRNA_works/pipeline/scRNA-seq_further_analysis/Get_colors.R")
# Normalization
FQnorm <- function(counts) {
    rk <- apply(counts, 2, rank, ties.method = "min")
    counts.sort <- apply(counts, 2, sort)
    refdist <- apply(counts.sort, 1, median)
    norm <- apply(rk, 2, function(r) {
        refdist[r]
    })
    rownames(norm) <- rownames(counts)
    return(norm)
}

#=================================================================================
# customized function definition
#=================================================================================

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("slingshot"))
suppressPackageStartupMessages(library("OESingleCell"))
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("uwot"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("grDevices"))
suppressPackageStartupMessages(library("stringr"))
# ================command line parameters setting============================= #
option_list = list(
    make_option( c("--RDS", "-i"), type = "character",
                 help = "the seurat object saved as R object in RDS format."),
    make_option( c("--informat", "-f"), type = "character",default = "rds",
                 help = "--RDS文件的文件类型，h5seurat:读取h5seurat,rds:读取rds."),
    make_option( c("--SCE", "-v"), type = "character",
                 help = "the SingleCellExperiment(sce) object saved as R object in RDS format."),
    make_option( c("--gamlist", "-e"), type = "character", default = NULL,
                help = "[OPTIONAL]SingleCellExperiment of GAM models after fitGAM, used to demonstrate the various tests."),
    make_option( c("--seurat","-b"), type = "logical", default=FALSE,
                help="Whether to perform slingshot operations based on the seurat object."),
    make_option( c("--output","-o"),type="character", default = "./",
                help="the output directory of results.", metavar="character"),
    make_option( c("--topn", "-n"), type="integer", default = 25,
                 help = "the number of top markers for each cluster to visualizse."),
    make_option(c("--toptype"), type = "character", default = "both",
                 help = "choose from up,down,or both DEGs to visualizse."),
    make_option( c("--curve", "-c"), type = "double", default = NULL,
                 help="Select the curve used to draw the heatmap."),
    make_option( c("--groupby", "-g"), type = "character", default = NULL,
                help = "[OPTIONAL]The grouppinig variable in the metadata for
                        separate the cells to visulize marker genes."),
    make_option( c("--start", "-a"), type = "character", default = NULL,
                 help = "[OPTIONAL]The starting group specified in groupby."),
    make_option( c("--end", "-z"), type = "character", default = NULL,
                 help = "[OPTIONAL]The ending group specified in groupby."),
    make_option( c("--pointsize", "-s"), type = "double", default = NULL,
                help = "[OPTIONAL]the point size in the plot."),
    make_option( c("--reduct"), type = "character", default = "UMAP",
                help = "the previous calculated reduction result used in the slingshot. You can use either PCA/UMAP/TSNE for reduction."),
    make_option( c("--alpha2use"), type = "double", default = 0.5,
                help = "[OPTIONAL]the opacity of the points on the Smoother plot."),
    make_option( c("--use_color_anno" ), type = "logical",  default = TRUE,
                help = "[Optional]是否采用rds中注释的颜色信息，默认采用，若无则自动重新注释颜色。"),
    make_option( c("--color_file" ), type = "character",  default = NULL,
                help = "[Optional]选填，输入以tab分隔的文件，第一列的列名为metadata列名，第一列为该列元素，第二列为对应的颜色."),
    make_option( c("--palette" ), type = "character",  default = NULL,
                help = "[Optional]选填，根据需求指定离散型色板名."),
    make_option( c("--test"), type = "character", default = "startVsEndTest",
                help = "[OPTIONAL]Test methods between lineages for the Smoother plot. Could be startVsEndTest, patternTest or diffEndTest."),
    make_option(c("--vismethod","-m"), type= "character",default=NULL,
                help = "the visulization methods for the marker genes of each cell cluster.
                the methods can be heatmap,Smoothersplot."),
    make_option( c("--extraGene", "-x"), type = "character", default = NULL,
                 help = "[OPTIONAL]The extra gene list of interest to visualize specified by the user."),
    make_option( c("--aimCurve", "-d"), type = "character", default = NULL,
                 help = "[OPTIONAL]The custom curve to show."),
    make_option( c("--predicate"), type = "character", default = NULL,
                help = "The conditional expression to subset cells used for subtyping.[default: %(default)s]"),
    make_option( c("--shownames"), type = "character", default = NULL,
                 help = "[OPTIONAL]热图是否展示基因名."),
    make_option( c("--clusters_num"), type = "integer", default = 4,
                help = "热图聚类分的模块数量.[default: %(default)s]"),
    make_option( c("--annotation"), type = "character", default = NULL,
                help = "用于热图模块各module注释.[default: %(default)s]"),
    make_option( c("--subnew_celltype"),type = "character", default = "all"),
    make_option( c("--subsampleid"),type = "character", default = "all"),
    make_option( c("--subgroup"),type = "character", default = "all"),
    make_option( c("--subclusters"),type = "character", default = "all"),
    make_option( c("--prefix"),type = "character", default = "prefix")
    );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
# ======================================================= #
# ========== parse the command line parameters ========== #
# ======================================================= #
if (is.null(opt$RDS) & is.null(opt$SCE)){
    stop("There is no seurat object or SCE object AVAILABLE!")
}

if ( !is.null(opt$RDS) ){
    if (opt$informat == "rds"){
        seurat_ob <- readRDS(opt$RDS)
    } else {
        seurat_ob <- OESingleCell::ReadX(input = opt$RDS, informat = opt$informat)
    }
    if ( seurat_ob@version < 3){
        seurat_ob <- UpdateSeuratObject(seurat_ob) #make sure the seurat object match with the latest seurat package
    }
   if ( !is.null(opt$predicate) ){
      desired_cells= subset(seurat_ob@meta.data, eval( parse(text=opt$predicate)))
      seurat_ob = seurat_ob[, rownames(desired_cells)]
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

 }


if ( is.null( opt$groupby ) ){
    print( "NO groupping variable AVAILABLE for cell groupping! The default cell clusters id will be used!")
    groupby <- "clusters"
}else{
    groupby <- unlist(strsplit(opt$groupby,","))
}

if ( is.null(opt$pointsize) ){
    pointsize <- 1
} else {
    pointsize <- opt$pointsize
}

if ( is.null(opt$reduct) ){
    print( "NO specified reduction AVAILABLE! The default reduction will be used!")
    reduct <- "UMAP"
}else{
    reduct <- toupper(opt$reduct)
}

if ( is.null(opt$output) ){
    print("NO output directory specified,the current directory will be used!")
    output_dir = getwd()
}else{
    if ( file.exists(opt$output) ){
        output_dir = opt$output
    }else{
        output_dir = opt$output
        dir.create(output_dir,recursive = T)
    }
}
output_dir = normalizePath(output_dir )

# my_palette=c(
    # "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    # "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    # "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    # "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    # "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
if ( is.null(opt$palette ) ){
    print("没有指定色板，将采用rds中注释的颜色或者默认色板.")
    palette = "customecol2"
}else{
    palette = opt$palette
}


#======================== Prepare the SingleCellExperiment object ========================#
if ( is.null(opt$SCE) ){
    print("the SCE object is NOT AVAILABLE and will be built from the seurat object.")
    if ( opt$seurat == F ){
        sce <- SingleCellExperiment(assays = list(counts = seurat_ob$RNA@counts), colData = seurat_ob@meta.data)
        # Gene Filtering
        geneFilter <- apply(assays(sce)$counts, 1, function(x) {
            sum(x >= 3) >= 10
        })
        sce <- sce[geneFilter, ]
        assays(sce)$norm <- FQnorm(assays(sce)$counts)
        ## Dimensionality Reduction
        ## pca
        pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
        rd1 <- pca$x[, 1:2]
        # pdf(file.path(output_dir, 'pca.pdf') )
        # plot(rd1, col = rgb(0, 0, 0, .5), pch = 16, asp = 1)
        # dev.off()
        ## UMAP
        reduct <- "UMAP"
        rd2 <- uwot::umap(t(log1p(assays(sce)$norm)) )
        colnames(rd2) <- c("UMAP1", "UMAP2")
        # pdf('umap-slingshot2.pdf')
        # plot(rd2, col = rgb(0, 0, 0, .5), pch = 16, asp = 1)
        # dev.off()
        reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)
        saveRDS(rd2, file.path(output_dir, "rd2.rds"))
    } else {
        sce <- as.SingleCellExperiment(seurat_ob)
    }
    #========================== Using Slingshot =============================#
    if ( !is.null(opt$start) ){
        start_clus = unlist(strsplit(opt$start,","))
    } else { start_clus = opt$start }
    if ( !is.null(opt$end) ){
        end_clus = unlist(strsplit(opt$end,","))
    } else { end_clus = opt$end }
    sce <- slingshot(sce, clusterLabels = groupby, reducedDim = reduct, start.clus= start_clus, end.clus = end_clus)
    saveRDS(sce, file.path(output_dir, 'sce.rds') )
} else {
    sce <- readRDS(opt$SCE)
}

if ( ! opt$use_color_anno ){
    colData(sce) = colData(sce)[ ,!grepl(paste0("^",groupby,"_col$" ), colnames(colData(sce)))]
}
if ( !is.null(opt$color_file)){
    color_file = read.delim(opt$color_file, sep="\t", header = T)
    meta_anno = color_anno(colData(sce), color_file)
} else {
    meta_anno = colData(sce)
}
color_use = get_colors(meta_anno, groupby, palette)
colData(sce) = color_use[["object_meta"]]
# user_color_pal = color_use[["user_color_pal"]]
new_celltype_pal = color_use[["new_celltype_pal"]]
new_celltype_pal = na.omit(new_celltype_pal)
print("color.use :")
print(new_celltype_pal)
print(table(colData(sce)[, paste0( groupby,"_col" )]))
if ('rawbc' %in% colnames(colData(sce))) {
    Barcode_content = 'rawbc'
}else{
    Barcode_content = 'orig.ident'
}

colData(sce)[[groupby]] <- factor( colData(sce)[[groupby]] )
new_celltype_pal = new_celltype_pal[levels(colData(sce)[[groupby]])]
if ( is.null(opt$aimCurve) ){
    ### Clustering plot 
    # groupby_clust <- colData(sce)[[paste0( groupby,"_col" )]]
    # levels(groupby_clust) <- c(1:length(levels(colData(sce)[[groupby]])))
    pdf(file.path(output_dir, 'Clustering_slingshot.pdf'),width = 8)
    par(mai=c(1,1,1,2))
    plot(reducedDims(sce)[[reduct]], col = colData(sce)[[paste0( groupby,"_col" )]], pch = 16, asp = 1, cex = pointsize)
    xy=par("usr")
    lines(SlingshotDataSet(sce), lwd = 1, col = "black")
    legend(x=xy[2]+xinch(0.2), y=xy[4], xpd = TRUE,                                    #图例位置
        legend = names(new_celltype_pal),        #图例内容
        col = new_celltype_pal,                 #图例颜色
        pch = 19 )
    dev.off()

    ### Pseudotime plot 
    colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
    pt <- slingPseudotime(sce)
    nms <- colnames(pt)
    nc <- ifelse(length(nms)>3,3,length(nms))
    nr <- ceiling(length(nms)/nc)
    pdf(file.path(output_dir, 'Slingshot_Pseudotime.pdf'), width = nc*5 , height = nr*5)
    par(mfrow = c(nr, nc))
    for (i in nms) {
      colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
      plotcol <- colors[cut(pt[,i], breaks=100)]
      plotcol[is.na(plotcol)] = "lightgray"
      plot(reducedDims(sce)[[reduct]], col = plotcol, pch = 16, cex = pointsize, main = i)
      lines(SlingshotDataSet(sce)@curves[[i]], lwd = 2, col = 'black') # also could get the lineages by "type = 'lineages'"
    }
    dev.off()
} else {
    aimCurve = paste0("curve",opt$aimCurve)
    groupby_clust <- colData(sce)[[groupby]]
    levels(groupby_clust) <- c(1:length(levels(colData(sce)[[groupby]])))
    pdf(file.path(output_dir, 'Clustering_slingshot_with_1_curve.pdf'),width = 8)
    par(mai=c(1,1,1,2))
    plot(reducedDims(sce)[[reduct]], col = colData(sce)[[paste0( groupby,"_col" )]], pch = 16, asp = 1, cex = pointsize)
    xy=par("usr")
    lines(SlingshotDataSet(sce)@curves[[aimCurve]], lwd = 1, col = "black")
    legend(x=xy[2]+xinch(0.2), y=xy[4], xpd = TRUE,                                    #图例位置
        legend = names(new_celltype_pal),        #图例内容
        col = new_celltype_pal,                 #图例颜色
        pch = 19 )    
    dev.off()
}    

# =========================================================================== #
# =========================== Downstream Analysis =========================== #
# =========================================================================== #

if ( !is.null(opt$vismethod) ){
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("pheatmap"))
    suppressPackageStartupMessages(library("tibble"))
    suppressPackageStartupMessages(library("gridExtra"))
    suppressPackageStartupMessages(library("tradeSeq"))
    source("/public/scRNA_works/pipeline/scRNA-seq_further_analysis/plotsmoothers_new.r")

    pt <- slingPseudotime(sce)
    vismethods <- unlist(strsplit(opt$vismethod,","))
    if ( is.null(opt$extraGene) ){ 
        if ( is.null(opt$gamlist )){ 
            print("Start to fitGAM!")
            set.seed(123)
            sce0 <- fitGAM(counts=as.matrix(assays(sce)$counts), sds=SlingshotDataSet(sce), pseudotime = pt, cellWeights = slingCurveWeights(sce), verbose = T ) # take some time
            saveRDS(sce0, file.path(output_dir, "gamlist_sce.rds"))
        } else { 
            sce0 <- readRDS( opt$gamlist ) }
    } else {
        gene <- read.delim(opt$extraGene, sep = "\t")
        if (dim(gene)[2] > 1) {
            up <- filter(gene, FoldChange > 1) %>%
                arrange(desc(log2FoldChange)) %>%
                top_n(opt$topn, log2FoldChange) %>%
                select(gene)
            down <- filter(gene, FoldChange < 1) %>%
                arrange(log2FoldChange) %>%
                top_n(as.numeric(paste0("-", opt$topn)), log2FoldChange) %>%
                select(gene)
            if (opt$toptype == "up") {
                gene <- up
            } else if (opt$toptype == "down") {
                gene <- down
            } else {
                gene <- rbind(up, down)
            }
        }
        filtered <- gene[!gene[, 1] %in% rownames(sce),1 ]
        if(length(filtered)!=0){
        filtered = as.data.frame(filtered)
        write.table(filtered,file.path(output_dir,"genes_not_matched.xls"),quote = F,row.names=F)
        print("There are some mismatched gene symbol, Please check genes_not_matched.xls for the genename.")
        }
        gene[, 1] <- factor(as.character(gene[, 1]))
        topgenes <- Seurat::CaseMatch(search = levels(gene[, 1]), match = rownames(sce))
        print("Start to fitGAM!")
        set.seed(123)
        sce0 = fitGAM(counts=as.matrix(assays(sce)$counts), sds=SlingshotDataSet(sce), pseudotime = pt, cellWeights = slingCurveWeights(sce), genes = topgenes, verbose = T )
    }

    for ( vismethod in vismethods ){
        if ( vismethod == "heatmap" ){
            if ( !is.null(opt$curve) ){
                pt = subset(pt, select= paste0("curve", opt$curve))
            }
            ATres <- associationTest(sce0, lineages= T, global = T )
            ATres <- as.data.frame(ATres)
            write.table(rownames_to_column(as.data.frame(ATres),var="GeneID"),file.path(output_dir, "Dynamic_genes.xls"), sep="\t",quote=F, row.names = F)
            for (i in colnames(pt)){
                curve_number <- gsub("curve", "", i)
                if ( is.null(opt$extraGene)){
                    topgenes <- rownames(ATres[order(-ATres[,paste0("waldStat_",curve_number)], ATres[,paste0("pvalue_",curve_number)]), ])[1:opt$topn]
                }
                if(is.null(opt$shownames)){
                    shownames = ifelse(length(topgenes) > 100, FALSE, TRUE)
                } else {
                    shownames = opt$shownames
                }
                heatmap = Dynamic_genes_heatmap(sce0, topgenes, nPoints = 100, lineage = curve_number, num_clusters = opt$clusters_num, show_rownames = as.logical(shownames))
                # ggsave(file.path(output_dir, paste0('Dynamic_genes_heatmap_',i,'.png')),plot=p, width = 8, height = 9)
                ggsave(file.path(output_dir, paste0('Dynamic_genes_heatmap_',i,'.pdf')),plot = heatmap, width = 8, height = 9,bg="white")
                gene_clusters <- cutree(heatmap$tree_row, k = opt$clusters_num)
                gene_clustering <- data.frame(gene_clusters) %>% tibble::rownames_to_column(var = "gene")
                gene_clustering[, 2] <- as.character(gene_clustering[, 2])
                colnames(gene_clustering) <- c("gene","gene_module")
                ATres_sub =  ATres[topgenes, c(paste0("waldStat_",curve_number), paste0("pvalue_",curve_number), paste0("df_",curve_number)) ] %>% rownames_to_column(var="gene")
                gene_clustering <- left_join(gene_clustering, ATres_sub, by = "gene")
                if(!is.null(opt$annotation)){
                    anno <- read.delim(opt$annotation, stringsAsFactors = F, sep = "\t", quote = "")
                    gene_clustering <- left_join(gene_clustering, anno, by = c("gene" = "id"))
                    gene_clustering[is.na(gene_clustering)] <- "--"
                }
                write.table(gene_clustering, file.path(output_dir, paste0("Dynamic_genes_heatmap_",i,"_module_anno.xls")), sep = "\t", quote = F, row.names = F)
            } 
        } 
        if ( vismethod == "Smoothersplot" ){
            if ( is.null(opt$extraGene)){
                if (opt$test == "patternTest" ){
                    Res <- patternTest(sce0, global =T )
                } else if (opt$test == "diffEndTest"){
                    Res <- diffEndTest(sce0, global =T , pairwise=TRUE)
                } else {
                    Res <- startVsEndTest(sce0, lineages=T, global =T )
                }
                Res <- as.data.frame(Res)
                curves <- NULL
                if ( !is.null(opt$curve) ){
                    curves <- paste0("_lineage",opt$curve)
                }
            topgenes <- rownames(Res[order(-Res[,paste0("waldStat",curves)], Res[,paste0("df",curves)]),])[1:opt$topn]
            write.table(rownames_to_column(Res[topgenes,],var="GeneID"),file.path(output_dir, paste0(opt$test,"_Top",opt$topn,"_genes.xls")), sep="\t",quote=F, row.names = F)
            }
            #source("/public/scRNA_works/pipeline/scRNA-seq_further_analysis/plotsmoothers.r")
            # for (sigGeneStart in topgenes){
                # colData(sce0)=cbind(colData(sce0),colData(sce))
                # p = plotSmoothers(sce0, assays(sce)$counts, gene = sigGeneStart, alpha = opt$alpha2use, pointCol = groupby )
                # ggsave(file.path(output_dir, paste0('pseudotime_genes_Smoothers_plot_',sigGeneStart,'.png')),plot=p, width = 6, height = 5)
                # ggsave(file.path(output_dir, paste0('pseudotime_genes_Smoothers_plot_',sigGeneStart,'.pdf')),plot=p, width = 6, height = 5)
            # } 
            colData(sce0)=cbind(colData(sce0),colData(sce))
            gs <- lapply(topgenes, function(x) plotSmoothers_sce(sce0[topgenes,], assays(sce)$counts[topgenes,], gene = x, lwd = 1, aim_curve=opt$aimCurve,
                                        alpha = opt$alpha2use, pointCol = groupby , point_colors = new_celltype_pal, border = TRUE) + labs(title = x) )
            i <- 1
            while (i <=  length(gs)){
              gs_sub <- gs[i:(i+9)]
              gs_sub <- gs_sub[which(!sapply(gs_sub, is.null))]

              pdf(file.path(output_dir,paste0("pseudotime_genes_Smoothers_plot",ifelse(i==1,"",i%/%10+1),".pdf", collapse = "_")),
                              width = 10, height = length(gs_sub)*1.5)
              grid.arrange(grobs = gs_sub, ncol=2) 
              dev.off()
              # png(file.path(output_dir,paste0("pseudotime_genes_Smoothers_plot",ifelse(i==1,"",i%/%10+1),".png",collapse = "_")),
                              # width = 10, height = length(gs_sub)*1.5 , res = 96, units = "in")
              # grid.arrange(grobs = gs_sub, ncol=2) 
              # dev.off()
              i=i+10
            }
        }
        if ( vismethod == "featureplot" ){
            if ( is.null(opt$extraGene)){
                if (opt$test == "patternTest" ){
                    Res <- patternTest(sce0, global =T )
                } else if (opt$test == "diffEndTest"){
                    Res <- diffEndTest(sce0, global =T , pairwise=TRUE)
                } else {
                    Res <- startVsEndTest(sce0, lineages=T, global =T )
                }
                Res <- as.data.frame(Res)
                curves <- NULL
                if ( !is.null(opt$curve) ){
                    curves <- paste0("_lineage",opt$curve)
                }
            topgenes <- rownames(Res[order(-Res[,paste0("waldStat",curves)], Res[,paste0("df",curves)]),])[1:opt$topn]
            write.table(rownames_to_column(Res[topgenes,],var="GeneID"),file.path(output_dir, paste0(opt$test,"_Top",opt$topn,"_genes.xls")), sep="\t",quote=F, row.names = F)
            }
            # colors <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(99)
            colors <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)
            for (gene in topgenes){
                gene_expr <- log1p(assays(sce)$counts[gene, ])
                gene_expr_data = as.data.frame(gene_expr)
                gene_expr_data$plotcol <- colors[cut(gene_expr_data$gene_expr, breaks=100)]
                gene_expr_data = gene_expr_data[order(gene_expr_data$gene_expr),]
                data_label <- round(seq(min(gene_expr), max(gene_expr), length.out = 5),2)
                atLims = seq(0.5,0.7, length.out = 5)
                pdf(file.path(output_dir,paste0("genes_",gene,"_Featureplot.pdf")))
                par(fig=c(0, 0.9, 0, 0.9))
                plot(reducedDims(sce)[[reduct]][rownames(gene_expr_data),], main = gene, pch = 16, cex = pointsize,
                     axes = F, xlab = "", ylab = "", col = gene_expr_data$plotcol)
                lines(SlingshotDataSet(sce), lwd = 2, col = "black")
                par(fig=c(0.5, 1, 0, 0.9),new=TRUE,adj = 0)
                image(x=c(0.92,0.97),y=seq(0.5,0.7,length.out=100),z=matrix(seq(0,1,length.out=100),nrow=1), 
                      col = colors, axes=F, xlim=c(0,1), ylim=c(0,1), xlab=NA, ylab=NA)
                text(x = 1, y = atLims, labels = data_label, srt = 0, cex = 0.8, xpd = TRUE)
                dev.off()
            }
        }
    }
}

if(!file.exists(file.path(output_dir, "Slingshot分析说明.docx"))){
  file.copy("/public/scRNA_works/Documents/Slingshot分析说明.docx",
  file.path(output_dir, "Slingshot分析说明.docx"))
}
setwd(output_dir)
print("Convert pdf to png...")
system("for i in  `ls *.pdf`;do /data/software/ImageMagick/ImageMagick-v7.0.8-14/bin/convert  -verbose -density 500 -trim  $i  -quality 100  -flatten  ${i/.pdf/.png}  ;done")

