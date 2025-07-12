#!/home/lipeng/miniconda3/envs/Seurat3.1/bin/Rscript
.libPaths("/home/lipeng/miniconda3/envs/Seurat3.1/lib/R/library")
##fix table information and change default color of heatmap color
##add scFEA type 20220815 by lip
##add ggstatsplot in vlnplot 20221009 by lip

print("must activate envs && source /home/lipeng/miniconda3/bin/activate Seurat3.1")

library("optparse")

option_list = list(
    make_option( c("--RDS", "-v"), type = "character", default = NULL,
                 help = "the seurat object saved as R object in RDS format."),
    make_option( c("--output","-o"),type="character", default = "./",
                help="the output directory of results.", metavar="character"),
    make_option(c("--collapseby","-q"),type="character",
           help="The variable level of heatmap. e.g. clusters or sampleid or celltype"),
    make_option( c("--colors", "-c"),  type = "character", default = NULL,
                help = "colors choise for Heatmap picture : redwhiteblue, redblackgreen ,yellowblackblue "),
    make_option( c("--topn", "-n"), type="integer", default = 10,
                 help = "the number of top Term for heatmap"),
    make_option( c("--pvalue", "-p"), type="double", default = 0.05,
                 help = "the pvalue "),
    make_option( c("--type", "-t"), type = "character", default = "REACTOME ",
                 help="KEGG or REACTOME(only human) ,or scFEA(human & mouse,so slowly)"),
    make_option( c("--gmt", "-g"), type = "character", default = NULL,
                 help="[OPTION:]input your KEGG or REACTOME gmtfile"),
    make_option( c("--species" ,"-s"), type="character", default = "human",
                help="species of you input file,human or mouse"),
    make_option( c("--contrast", "-d"),type = "character",default = "group:all:all",
            help = "[Required]levels of a factor used to compare with for final differenetial results.
                    The format is Factor:interesting_level:reference_level."),
    make_option( c("--cpu"), type="integer", default = 2,
                 help = "the number of cpu."),
    make_option( c("--pair"), type = "character", default = "False",
        help = "Whether to display p value ,True or False"),
    make_option( c("--method", "-m"), type = "character", default = "VISION",
                 help = "method supports VISION, AUCell, ssGSEA"),
    make_option( c("--predicate" ), type = "character", default = NULL, 
             help = "[OPTIONAL]The column name in cell metadata used as identity of each cell combined with which_cell."),
    make_option( c("--subnew_celltype"),type = "character", default = "all"),
    make_option( c("--subsampleid"),type = "character", default = "all"),
    make_option( c("--subgroup"),type = "character", default = "all"),
    make_option( c("--subclusters"),type = "character", default = "all"),
    make_option( c("--prefix"),type = "character", default = "prefix")
    );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

suppressWarnings({
    suppressPackageStartupMessages(library("VISION"))
    suppressPackageStartupMessages( library("Seurat"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("gridExtra"))
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("tibble"))
    suppressPackageStartupMessages(library("RColorBrewer"))
    suppressPackageStartupMessages(library("wesanderson"))
    suppressPackageStartupMessages(library(pheatmap))
    suppressPackageStartupMessages(library(scMetabolism))
    suppressPackageStartupMessages(library(rsvd))
    library(stringr)
})

#1.check file and param
if ( is.null(opt$RDS) ){
    stop("the seurat object is NOT AVAILABLE!")
}else{
    seurat_ob = readRDS(opt$RDS)
    dim(seurat_ob)
}
if ( file.exists(opt$output)){
    output_dir = opt$output
}else{
    output_dir = opt$output
    dir.create(output_dir,recursive = TRUE)
}
#if ( seurat_ob@version < 4){
#    my_test = Seurat::UpdateSeuratObject(object=seurat_ob) #make sure the seurat object match with the latest seurat package
#} else {
    my_test = seurat_ob
#}
if ( !is.null(opt$colors) ){
    if ( opt$colors == "redwhiteblue" ){
        palette <- colorRampPalette(c("RoyalBlue2", "White", "Red2"))(n=256)
    }else if ( opt$colors == "redblackgreen" ){
        palette <- colorRampPalette(c("Green", "Black", "Red"))(n=256)
    }else if ( opt$colors == "yellowblackblue" ){
        palette <- colorRampPalette(c("Blue", "Black", "Yellow"))(n=256)
    }
}else {
    palette <- colorRampPalette(c("steelblue", "White", "firebrick2"))(n=256)
}

if (!is.null(opt$predicate)) {
    df <- my_test@meta.data
    desired_cells <- subset(df, eval(parse(text = opt$predicate)))
    my_test <- my_test[, rownames(desired_cells)]
}

if (opt$subnew_celltype != "all"){
    subnew_celltype = str_split(opt$subnew_celltype,",")[[1]]
    print(subnew_celltype)
    my_test = my_test[,my_test@meta.data$new_celltype %in% subnew_celltype]
    print(unique(my_test@meta.data$new_celltype))
}
if (opt$subsampleid != "all"){
    subsampleid = str_split(opt$subsampleid,",")[[1]]
    my_test = my_test[,my_test@meta.data$sampleid %in% subsampleid]
    print(unique(my_test@meta.data$sampleid))
}
if (opt$subgroup != "all"){
    subgroup = str_split(opt$subgroup,",")[[1]]
    my_test = my_test[,my_test@meta.data$group %in% subgroup]
    print(unique(my_test@meta.data$group))
}
if (opt$subclusters != "all"){
    subclusters = str_split(opt$subclusters,",")[[1]]
    my_test = my_test[,my_test@meta.data$clusters %in% subclusters]
    print(unique(my_test@meta.data$clusters))
}

print(dim(my_test))
print(table(my_test$group))

#2.some default module
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100))

CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}

my_sc.metabolism.Seurat <-function (obj, gmtFile,method, imputation=F, ncores,metabolism.type){
    countexp <- obj@assays$RNA@counts
    countexp <- data.frame(as.matrix(countexp))
    if (imputation == F) {
        countexp2 <- countexp
    }else if (imputation == T) {
        cat("Start imputation...\n")
        cat("Citation: George C. Linderman, Jun Zhao, Yuval Kluger. Zero-preserving imputation of scRNA-seq data using low-rank approximation. bioRxiv. doi: https://doi.org/10.1101/397588 \n")
        result.completed <- alra(as.matrix(countexp))
        countexp2 <- result.completed[[3]]
        row.names(countexp2) <- row.names(countexp)
    }
    cat("Start quantify the metabolism activity...\n")
    if (method == "VISION") {
        library(VISION)
        n.umi <- colSums(countexp2)
        scaled_counts <- t(t(countexp2)/n.umi) * median(n.umi)
        vis <- Vision(scaled_counts, signatures = gmtFile)
        options(mc.cores = ncores)
        vis <- analyze(vis)
        signature_exp <- data.frame(t(vis@SigScores))
    }
    if (method == "AUCell") {
        library(AUCell)
        library(GSEABase)
        cells_rankings <- AUCell_buildRankings(as.matrix(countexp2),
            nCores = ncores, plotStats = F)
        geneSets <- getGmt(gmtFile)
        cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
        signature_exp <- data.frame(getAUC(cells_AUC))
    }
    if (method == "ssGSEA") {
        library(GSVA)
        library(GSEABase)
        geneSets <- getGmt(gmtFile)
        gsva_es <- gsva(as.matrix(countexp2), geneSets, method = c("ssgsea"),
            kcdf = c("Poisson"), parallel.sz = ncores)
        signature_exp <- data.frame(gsva_es)
    }
    if (method == "ssGSVA") {
        library(GSVA)
        library(GSEABase)
        geneSets <- getGmt(gmtFile)
        gsva_es <- gsva(as.matrix(countexp2), geneSets, method = c("gsva"),
            kcdf = c("Poisson"), parallel.sz = ncores)
        signature_exp <- data.frame(gsva_es)
    }
    obj@assays$METABOLISM$score <- signature_exp
    obj
}

center_scale <- function(data){
	data = data
	center_data = scale(data,center=T,scale=F)
	return(center_data)
}

#3.check gtmFile
if(opt$type %in% names(my_test@assays)){
    print("The score already calculated and written to the rds")
    metabolism.matrix  <- as.data.frame(my_test[[opt$type]]@counts)
}else {
    if(opt$type == "scFEA"){
        scFEA_py = "/home/lipeng/miniconda3/envs/Metabolism/share/scFEA/src/scFEA.py"
        moduleGene_dir = "/home/lipeng/miniconda3/envs/Metabolism/share/scFEA/data/"
        if ( opt$species =="human") {
            module_file1 = "module_gene_human_m168.csv"
            module_file2 = "cmMat_human_c70_m168.csv"
            module_file3 = "cName_human_c70_m168.csv"
        }else if ( opt$species =="mouse") {
            module_file1 = "module_gene_mouse_m168.csv"
            module_file2 = "cmMat_mouse_c70_m168.csv"
            module_file3 = "cName_mouse_c70_m168.csv"
        }
        in_dir = file.path(output_dir,"input")
        out_dir = file.path(output_dir,"output")
        if (!file.exists(in_dir)){
            dir.create(in_dir, recursive = T)
        }
        if (!file.exists(out_dir)){
            dir.create(out_dir, recursive = T)
        }
        mtx = as.matrix(GetAssayData(my_test))
        mtx = as.data.frame(mtx) %>% tibble::rownames_to_column(var="genes")
        mtx_file = file.path(in_dir, "/exp_mat.csv")
        write.table(mtx, mtx_file, sep = ",", quote = FALSE, row.names=F)
        if("scFEA" %in% names(my_test@assays)){
            print("scFEA score has added in this obj,skip calculated")
            metabolism.matrix <- as.data.frame(my_test[['scFEA']]@counts)
        } else {
            system(glue::glue("module purge && source /home/lipeng/miniconda3/bin/activate /home/lipeng/miniconda3/envs/Metabolism &&
            python {scFEA_py}   --input_dir {in_dir}  --test_file exp_mat.csv  --sc_imputation True  --data_dir {moduleGene_dir}  --moduleGene_file {module_file1} --stoichiometry_matrix {module_file2}  --cName_file {module_file3} --res_dir {out_dir} && source /home/lipeng/miniconda3/bin/deactivate && module load OESingleCell/2.0.0 ")) # #&& source ~/.bashrc
            path_score = paste0(out_dir,"/balance_result.csv")
            score_met  = read.delim(path_score,sep=",",row.names=1,check.names=F)
            metabolism.matrix  = as.data.frame(t(score_met))
        }
    } else {
        if( is.null(opt$gmt)){
          if ( opt$species =="mouse") { 
            print("Loading default mouse gmtFile")
            opt$gmt = glue::glue("/public/scRNA_works/works/donghongjie/project/test/test_scmetabolism/data/{opt$type}_metabolism.gmt")
        	  countexp.Seurat<-my_sc.metabolism.Seurat(obj = my_test,method = opt$method, imputation = F, ncores = opt$cpu, metabolism.type = opt$type,gmtFile = opt$gmt)
          }else{ 
            print("Loading default Human gmtFile")
            countexp.Seurat<-sc.metabolism.Seurat(obj = my_test,method = opt$method, imputation = F, ncores = opt$cpu, metabolism.type = opt$type)
          }
        } else {
        	print("Loading your gmtFile")
        	countexp.Seurat<-my_sc.metabolism.Seurat(obj = my_test,method = opt$method, imputation = F, ncores = opt$cpu, metabolism.type = opt$type,gmtFile = opt$gmt)
        }
    	metabolism.matrix <- as.data.frame(countexp.Seurat@assays$METABOLISM$score)
    }
}
metabolism_score = tibble::rownames_to_column(metabolism.matrix,var = opt$type)
write.table(metabolism_score,quote = F,sep ="\t",row.names = F,col.names=T,
    file.path(output_dir,paste0("metabolism_",opt$type,"_score.xls",collapse = "")))

#4.save score

if ( opt$type == "KEGG" ){
    my_test[['KEGG']] = CreateAssayObject(counts = metabolism.matrix)
    count = as.data.frame(my_test@assays$KEGG@counts)
}else if ( opt$type == "REACTOME" ){
    my_test[['REACTOME']] = CreateAssayObject(counts = metabolism.matrix)
    count = as.data.frame(my_test@assays$REACTOME@counts)
}else if(opt$type == "scFEA"){
    my_test[['scFEA']] = CreateAssayObject(counts = metabolism.matrix )
    count = as.data.frame(my_test@assays$scFEA@counts)
}

filename= as.character(lapply(strsplit(basename(opt$RDS),split=".rds",fixed=TRUE),head,1))
if(!opt$type %in% names(seurat_ob@assays)){
	saveRDS(my_test,file.path(output_dir,paste0(filename,".metabolism_",opt$type,"_score.rds",collapse = "")))
}
#saveRDS(my_test,file.path(output_dir,paste0("metabolism_",opt$type,"_score.rds",collapse = "")))

heatmax = max(metabolism.matrix)
heatmin = min(metabolism.matrix)

meta.data = my_test@meta.data
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

ave_score = tibble::rownames_to_column(as.data.frame(collapsed_count),var = opt$type)
write.table(ave_score,quote = F,sep ="\t",row.names = F,col.names=T,
    file.path(output_dir,paste0("average_",opt$type,"_score.xls",collapse = "")))

#5.plot
#prepare data
if ( dim(collapsed_count)[2]==2){
    scale_type  ="none"
    heatmap_data = as.data.frame(t(apply(collapsed_count,1,center_scale)))
    names(heatmap_data) = colnames(collapsed_count)
}else{
    scale_type="row"
    heatmap_data = collapsed_count
}
#remove rows which sd=0
heatmap_data <- heatmap_data[apply(heatmap_data, 1, function(x) sd(x)!=0),]

title = paste0("metabolism_",opt$type)
heatmap_plot = pheatmap(heatmap_data,main = title, treeheight_row=7, 
        lwd=1,scale=scale_type,border=F,border_color = "white",
        color = palette , show_rownames = T,show_colnames =T,,fontsize_col = 7,fontsize_row = 7, 
        cluster_rows = T , cluster_cols = F , angle_col = 45,
		cellheight = 8, cellwidth = 9
		)
ggsave(paste0(output_dir, "/average_", opt$type, "_heatmap.pdf"), plot = heatmap_plot, height = nrow(heatmap_data)/7+1, width=9,bg="white")
ggsave(paste0(output_dir, "/average_", opt$type, "_heatmap.png"), plot = heatmap_plot, height = nrow(heatmap_data)/7+1, width = 9, dpi = 600, limitsize = F,bg="white")
#ggsave(paste0(output_dir,"/average_",opt$type,"_heatmap.pdf"),plot=heatmap_plot,width=(dim(collapsed_count)[2]))
#ggsave(paste0(output_dir,"/average_",opt$type,"_heatmap.png"),plot=heatmap_plot,width=(dim(collapsed_count)[2]),dpi = 1000 ,limitsize = F)

##6.Find_diff
assay_metadata = my_test@meta.data
if ( is.null(opt$contrast ) ){ #no contrast is provided
    stop("NO contrast string is AVAILABLE!")
}else{
    contrast = opt$contrast
}
#contrast = paste0(collapseby,":all:all")
contrasts = unlist( strsplit(contrast,":",perl = T) )
all_levels = as.vector(unique(assay_metadata[,contrasts[1]]))
if ( contrasts[2] == "all" & contrasts[3] != "all" ){
    case_levels = all_levels[-which(all_levels==contrasts[3])] #delete the reference level
    all_comparisions = paste(contrasts[1],case_levels,contrasts[3],sep = ":")
}else if ( contrasts[2] != "all" & contrasts[3] == "all" ){
    ref_levels = all_levels[-which(all_levels==contrasts[2])] #delete the interested level
    all_comparisions = paste(contrasts[1],contrasts[2],ref_levels,sep = ":")
}else if ( contrasts[2] == "all" & contrasts[3] == "all" ){
    all_comparisions = lapply(all_levels,
                function(x) paste(contrasts[1],x,paste0(all_levels[-which(all_levels==x)],collapse = ","),sep = ":"))
    all_comparisions = unlist(all_comparisions)
}else{
    all_comparisions = contrast
}
diff_outdir = paste0(output_dir,"/1.Diff_Term",sep="")
dir.create(diff_outdir,recursive = TRUE)

if ( contrasts[2] == "all" & contrasts[3] == "all" ){
    # find the significant differential pathway for each Interested_group against all other Interested_group
    Diff_results = c()
    for ( contrastx in all_comparisions ){
		contrastsx = unlist(strsplit(contrastx,':',perl =T))
        DEG_META_tmp = FindMarkers( my_test, 
                                ident.1 = as.character(unlist( strsplit(contrastsx[2],",",perl = T) )), 
								ident.2 = as.character(unlist( strsplit(contrastsx[3],",",perl = T) )), 
								test.use = "bimod", 
								min.pct = 0, 
								group.by = contrastsx[1], 
								logfc.threshold = -Inf, 
								assay = opt$type, 
								slot = "data" )
        DEG_META_tmp$cluster = as.character(contrastsx[2])
        DEG_META = DEG_META_tmp %>% rownames_to_column(var = "geneset")
        Diff_results = rbind(Diff_results, DEG_META) 
    }
#	Diff_results = Diff_results %>% select(geneset, avg_log2FC, pct.1 ,pct.2, p_val_adj,cluster)
    colnames(Diff_results) = c("geneset", "pval", "avg_log2FC","pct_1", "pct_2", "padj", "Interested_group")
    #Diff_results = Diff_results[, c("geneset", "pct_1", "pct_2","pval",  "padj","avg_log2FC", "Interested_group")]
    Diff_results = Diff_results[, c("geneset","pval",  "padj","avg_log2FC", "Interested_group")]
    Diff_results = Diff_results%>% select( geneset, everything())
    topn_markers  = Diff_results %>% group_by(Interested_group) %>% 
              filter(pval < opt$pvalue ) %>%
              arrange(pval,desc(avg_log2FC)) %>%
              # filter(gene_diff > pct_fold_cutoff)  %>% 
              top_n(opt$topn,avg_log2FC)
			  
    dim(Diff_results)
	dim(topn_markers)
    write.table(Diff_results, file = file.path(diff_outdir,paste0(opt$type,"_all_results.xls", collapse = "")), quote=F, sep="\t", col.names=T, row.names=F)
    write.table(topn_markers, file = file.path(diff_outdir,paste0("top", opt$topn,"_", opt$type,"_Significant_results.xls", collapse = "")), quote=F, sep="\t", col.names=T, row.names=F)

#####plot
    sort_topn = topn_markers[order(topn_markers$Interested_group),]
    sort_topn$geneset = factor(sort_topn$geneset,levels = unique(sort_topn$geneset))
    plot_list = as.matrix(unique(sort_topn$geneset))
    plot_data = heatmap_data[as.vector(plot_list),]

    topn_score = tibble::rownames_to_column(as.data.frame(plot_data),var = opt$type)
    write.table(topn_score,quote = F,sep ="\t",row.names = F,col.names=T,
    file.path(diff_outdir,paste0("Significant_",opt$type,"_score.xls",collapse = "")))

    heatmap_plot2 = pheatmap(plot_data,main = "Significant_Term", treeheight_row=7, 
            lwd=1,scale=scale_type,
            border_color = "white",color = palette , show_rownames = T,show_colnames =T, 
            cluster_rows = F , cluster_cols = F , angle_col = 45 ,fontsize_col = 9,
            fontsize_row = 8, cellheight = 12, cellwidth = 12
        )
    ggsave(paste0(diff_outdir,"/Significant_",opt$type,"_heatmap.pdf"),plot=heatmap_plot2,width=(dim(plot_data)[2]+3.5),height = (dim(plot_data)[2]+3),bg="white")
    ggsave(paste0(diff_outdir,"/Significant_",opt$type,"_heatmap.png"),plot=heatmap_plot2,width=(dim(plot_data)[2]+3.5),height = (dim(plot_data)[2]+3),dpi = 600 ,limitsize = F,bg="white")
}else {
    for ( contrastx in all_comparisions ){
		contrastsx = unlist(strsplit(contrastx,':',perl =T))
        Diff_results = FindMarkers( my_test, 
                                ident.1 = as.character(unlist( strsplit(contrastsx[2],",",perl = T) )), 
								ident.2 = as.character(unlist( strsplit(contrastsx[3],",",perl = T) )), 
								test.use = "bimod", 
								min.pct = 0, 
								group.by = contrastsx[1], 
								logfc.threshold = -Inf, 
								assay = opt$type, 
								slot = "data" )

        Diff_results = Diff_results %>% rownames_to_column(var = "geneset")
        colnames(Diff_results) = c("geneset", "pval", "avg_log2FC","pct_1", "pct_2", "padj")
        Diff_results = Diff_results%>% select( geneset, pval,  padj,avg_log2FC )

        write.table(Diff_results,file.path(diff_outdir,
                        paste0("diffexp_genesets_", opt$type, "_score4",contrastsx[1],"_",contrastsx[2],
                        "-vs-",contrastsx[3],".xls")), quote = F,sep = '\t',row.names = F)


        Diff_results$just = ifelse( Diff_results$avg_log2FC<0,0,1)
        Diff_results$is.just = Diff_results$just==1 
        pval_cutoff=opt$pvalue
        Diff_results$is.sig = Diff_results$pval<pval_cutoff


        pathway2vis = Diff_results %>% filter( pval < pval_cutoff ) %>% group_by( is.just ) %>% arrange(abs(avg_log2FC)) %>% top_n(opt$topn,abs(avg_log2FC))

            pp = ggplot(pathway2vis, aes(reorder(geneset, avg_log2FC), avg_log2FC)) +
                    geom_col(aes(fill=is.just)) +
                    scale_fill_manual(values=c( "#6CC570","#2A5078")) +
                    coord_flip() +
                    labs(x="Pathway", y=paste0("avg_log2FC value of ", opt$type ) ) +
                    theme_minimal() +
                    geom_text( aes(x= geneset, y=0, label = geneset), hjust = pathway2vis$just, size = 3.5 )+
                    theme(axis.text.y=element_blank()) +
                    theme(panel.grid =element_blank())
            pp = pp + labs( fill = paste0("avg_log2FC value > 0") )

            ggsave(file.path(diff_outdir,
                paste0("top", opt$topn,"_diffexp_genesets_",opt$type,"_score4",contrastsx[1],"_",contrastsx[2],"-vs-",contrastsx[3],"_barplot.pdf")),plot = pp, height = 7,width = max(nchar(pathway2vis$geneset))*0.12+4, bg="white")
            ggsave(file.path(diff_outdir,
                paste0("top", opt$topn,"_diffexp_genesets_",opt$type,"_score4",contrastsx[1],"_",contrastsx[2],"-vs-",contrastsx[3],"_barplot.png")),plot = pp, height = 7,width = max(nchar(pathway2vis$geneset))*0.12+4,bg="white",dpi=300)

            write.table(pathway2vis[,1:4],file.path(diff_outdir,
                        paste0("top", opt$topn,"_diffexp_genesets_",opt$type,"_score4",contrastsx[1],"_",contrastsx[2],
                        "-vs-",contrastsx[3],".xls")), quote = F,sep = '\t',row.names = F)
            
            plot_list <- pathway2vis$geneset
    }
}

###vlnplot
if(collapseby =="clusters"){
        nlevels = sort(unique(my_test@meta.data$clusters))
        user_color_pal = CustomCol2(nlevels)
        nlevel = length(unique(my_test@meta.data[,collapseby]))
    }else if(collapseby =="new_celltype" & "new_celltype_col" %in% colnames(my_test@meta.data)){
            nlevel = length(unique(my_test@meta.data[,collapseby]))
            nlevel_list = sort(as.character(unique(my_test@meta.data[,collapseby])))
            tmp_df <- unique(my_test@meta.data[c("new_celltype","new_celltype_col")])
            new_celltype_pal <- as.vector(tmp_df$new_celltype_col)
            names(new_celltype_pal) <-  as.vector(tmp_df$new_celltype)
            new_celltype_pal = as.list(new_celltype_pal)
            user_color_pal = new_celltype_pal[nlevel_list]
    }else {
        nlevel = length(unique(my_test@meta.data[,collapseby]))
        user_color_pal = CustomCol2(1:nlevel)
    }

if( is.null(opt$gmt)){
    for(pathway in plot_list){
        pathway_name = gsub('_-_|_/_|__|,_','_',gsub("[[:space:]]", "_",gsub("[(.*)]","", pathway)))
        p1 = VlnPlot(my_test,assay = opt$type,features =pathway,group.by= collapseby,pt.size=0,cols = user_color_pal) + 
                NoLegend() +  geom_boxplot(width=.2)
		p2 = ggstatsplot::ggbetweenstats(as.data.frame(p1$data),y=!!sym(pathway),x=ident,centrality.plotting=F,bf.message =F,results.subtitle =F ,pairwise.comparisons = as.logical(opt$pair))+
        ggplot2::scale_colour_manual( values = user_color_pal) + #ggtitle(label = paste0("N = " ,dim(p1$data)[1])) +
        theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(),axis.line = element_line(colour = "black"),
			  axis.text.x = element_text(angle = 45,hjust = 1)) + 
        xlab(collapseby)
		#两组间比较加P值，annotation 
		if(length(table(p1$data$ident)) ==2){
            library(ggpubr)
			p_test = ggstatsplot::ggbetweenstats(as.data.frame(p1$data),y=!!sym(pathway),
                                 x=ident,centrality.plotting=F,bf.message =F,
                                 results.subtitle =T ,pairwise.comparisons =TRUE)
            p2 = p2 +  geom_signif(comparisons = list(names(table(p1$data$ident))),test = "wilcox.test",
									map_signif_level = F,na.rm = T,y_position = max(p2$data[pathway])+0.1,
									annotation = paste0("p = ",(as.list(p_test$labels$subtitle[[3]])[[3]]))
                  )
        }
        ggsave(paste0(diff_outdir,"/",opt$type,"_",pathway_name,"_violin_plot.png"),plot=p2,dpi = 300 ,limitsize = F,bg="white")
        ggsave(paste0(diff_outdir,"/",opt$type,"_",pathway_name,"_violin_plot.pdf"),plot=p2,limitsize = F,bg="white")
    }
} else {
    for(pathway in plot_list){
        pathway_name = gsub('_-_|_/_|__|,_','_',gsub("[[:space:]]", "_", pathway))
        p1 = VlnPlot(my_test,assay = opt$type,features =pathway,group.by= collapseby,pt.size=0,cols = user_color_pal) + 
                NoLegend() +  geom_boxplot(width=.2)
		p2 = ggstatsplot::ggbetweenstats(as.data.frame(p1$data),y=!!sym(pathway),x=ident,centrality.plotting=F,bf.message =F,results.subtitle =F ,pairwise.comparisons = as.logical(opt$pair))+
        ggplot2::scale_colour_manual( values = user_color_pal) + #ggtitle(label = paste0("N = " ,dim(p1$data)[1])) +
        theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(),axis.line = element_line(colour = "black"),
			  axis.text.x = element_text(angle = 45,hjust = 1)) + 
        xlab(collapseby)
		#两组间比较加P值，annotation 
		if(length(table(p1$data$ident)) ==2){
            library(ggpubr)
			p_test = ggstatsplot::ggbetweenstats(as.data.frame(p1$data),y=!!sym(pathway),
                                 x=ident,centrality.plotting=F,bf.message =F,
                                 results.subtitle =T ,pairwise.comparisons =TRUE)
            p2 = p2 +  geom_signif(comparisons = list(names(table(p1$data$ident))),test = "wilcox.test",
									map_signif_level = F,na.rm = T,y_position = max(p2$data[pathway])+0.1,
									annotation = paste0("p = ",(as.list(p_test$labels$subtitle[[3]])[[3]]))
                  )
        }
        ggsave(paste0(diff_outdir,"/",opt$type,"_",pathway_name,"_violin_plot.png"),plot=p2,dpi = 300 ,limitsize = F,bg="white")
        ggsave(paste0(diff_outdir,"/",opt$type,"_",pathway_name,"_violin_plot.pdf"),plot=p2,limitsize = F,bg="white")
    }
}

if(opt$type == "scFEA"){
    if(!file.exists(file.path(output_dir, "scFEA代谢分析说明.docx"))){
        file.copy("/public/scRNA_works/Documents/scFEA代谢分析说明.docx",
        file.path(output_dir, "scFEA代谢分析说明.docx"))}
} else {
    if(!file.exists(file.path(output_dir, "scMetabolism代谢分析说明.docx"))){
        file.copy("/public/scRNA_works/Documents/scMetabolism代谢分析说明.docx",
        file.path(output_dir, "scMetabolism代谢分析说明.docx"))}
}
