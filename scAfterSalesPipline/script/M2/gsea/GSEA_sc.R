#!/usr/bin/env Rscript

suppressWarnings({
    suppressPackageStartupMessages(library("Seurat"))
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("tibble"))
    suppressPackageStartupMessages(library("future"))
    suppressPackageStartupMessages(library("data.table"))
    suppressPackageStartupMessages(library("OESingleCell"))
    library(stringr)
})

#=================================================================================
# command line parameters setting
#=================================================================================
option_list = list(
    make_option( c("--input", "-i" ), type = "character",
                 help = "[Required] The input expression matrix in seurat format (rds file)."),
    make_option( c("--contrast", "-c"), type = "character",
                 help = "[Required] levels of a factor used to compare with for final differenetial results.
                    The format is Factor:interesting_level:reference_level."),
    make_option( c("--gmt", "-g" ), type = "character",
                 help = "[Required] the gene sets in gmt format from the MsigDB/KEGG/GO database etc."),
    make_option( c("--splitby", "-s"), type = "character", default = NULL,
                 help = "[OPTIONAL] The variable in the metadata used to run separately."),
    make_option( c("--ident2use", "-q" ), type = "character", default = NULL,
                 help = "[OPTIONAL] The column name in cell metadata used as identity of each cell combined with which_cell."),
    make_option( c("--which_cells", "-u" ), type = "character", default = NULL,
                 help = "[OPTIONAL] The subset of cluster ids used for analysis."),
    make_option( c("--minsize", "-m" ), type = "numeric", default = 15,
                 help = "[OPTIONAL] The minsize gene numbers of term for filtering gmt file ."),
    make_option( c("--maxsize"), type = "numeric", default = 500,
                 help = "[OPTIONAL] The maxsize gene numbers of term for filtering gmt file ."),
    make_option( c("--process", "-p"), type = "numeric", default = 8,
                 help = "[OPTIONAL] 多线程运行,默认线程数8."),
    make_option( c("--outdir","-o"), type="character", default = "./",
                 help="the output directory of GSEA results." ),
    make_option( c("--subgene"),type = "character", default = "all"),
    make_option( c("--subnew_celltype"),type = "character", default = "all"),
    make_option( c("--subsampleid"),type = "character", default = "all"),
    make_option( c("--subgroup"),type = "character", default = "all"),
    make_option( c("--subclusters"),type = "character", default = "all"),
    make_option( c("--prefix"),type = "character", default = "prefix"));
opt_parser = OptionParser(option_list=option_list);
opts = parse_args(opt_parser);

#=================================================================================
# parse the command line parameters
#=================================================================================
if ( !is.null(opts$input )){
    seurat_ob = readRDSMC( opts$input, cores = availableCores())
    # if the input seurat object version is less than 3, upgrade it to version 3
    if ( seurat_ob@version < 3 ){
        seurat_ob = UpdateSeuratObject(seurat_ob) # make sure the seurat object match with the latest seurat package
    }
}else{
    stop("NO seurat object is AVAILABLE!")
}
if ( !is.null(opts$which_cells)){
    cluster_list = unlist(strsplit( opts$which_cells,",",perl = T))
    seurat_ob = SubsetData(seurat_ob,  subset.name = opts$ident2use, accept.value = cluster_list)
}

if (opts$subnew_celltype != "all"){
    subnew_celltype = str_split(opts$subnew_celltype,",")[[1]]
    seurat_ob = seurat_ob[,seurat_ob@meta.data$new_celltype %in% subnew_celltype]
    print(unique(seurat_ob@meta.data$new_celltype))
}
if (opts$subsampleid != "all"){
    subsampleid = str_split(opts$subsampleid,",")[[1]]
    seurat_ob = seurat_ob[,seurat_ob@meta.data$sampleid %in% subsampleid]
    print(unique(seurat_ob@meta.data$sampleid))
}
if (opts$subgroup != "all"){
    subgroup = str_split(opts$subgroup,",")[[1]]
    seurat_ob = seurat_ob[,seurat_ob@meta.data$group %in% subgroup]
    print(unique(seurat_ob@meta.data$group))
}
if (opts$subclusters != "all"){
    subclusters = str_split(opts$subclusters,",")[[1]]
    seurat_ob = seurat_ob[,seurat_ob@meta.data$clusters %in% subclusters]
    print(unique(seurat_ob@meta.data$clusters))
}

if (opts$subgene != "all"){
    print(seurat_ob)
    data=read.table(opts$subgene,sep= '\t',header=T)
    seurat_ob = seurat_ob[rownames(seurat_ob) %in% data[,1],]
    print(seurat_ob)
}


if ( is.null(opts$contrast ) ){
    stop("NO contrast string is AVAILABLE!")
}else{
    contrast = opts$contrast
}
contrasts = unlist( strsplit(contrast,":",perl = T) )
col = contrasts[1]
case_group = contrasts[2]
control_group = contrasts[3]

if ( is.null(opts$gmt ) ){
    stop("NO gmt file is AVAILABLE!")
}else{
    gmt = opts$gmt
}
gmts = unlist( strsplit(gmt,",",perl = T) )

if ( is.null(opts$splitby ) ){
    split.by = NULL
}else{
    split.by = opts$splitby
}

if ( is.null(opts$outdir) ){
    print("NO output directory specified, the current directory will be used!")
    output_dir = getwd()
}else{
    if ( file.exists(opts$outdir) ){
        output_dir = opts$outdir
    }else{
        output_dir = opts$outdir
        dir.create(output_dir, recursive = T)
    }
}
output_dir = normalizePath(output_dir)


if (!is.null(split.by)){
    seurat_ob_split = SplitObject(seurat_ob, split.by = split.by)
    Idents(seurat_ob) = split.by
    mtx_list = lapply(seurat_ob_split, function(i){
      mtx = cbind(as.matrix(GetAssayData(i)[,i@meta.data[,col]==case_group]),
                    as.matrix(GetAssayData(i)[,i@meta.data[,col]==control_group]));mtx = as.data.frame(mtx) %>% rownames_to_column(var="genes");mtx})

    lapply(seq_along(mtx_list), function(i){
        out_dir = output_dir
        if (!file.exists(out_dir)){
            dir.create(out_dir, recursive = T)
        }
        mtx_file = file.path(out_dir,paste(names(mtx_list)[i],"_",case_group,"-vs-",control_group,".txt",sep=""))
        print(mtx_file)
        # write.table(mtx_list[[i]], mtx_file, sep = "\t", quote = FALSE, row.names=F)
        fwrite(mtx_list[[i]], file = mtx_file, sep = "\t", quote = FALSE, row.names=F)
        cls3 = paste0(sub("_[ATCG]{16,}.*","",colnames(mtx_list[[i]])[2:ncol(mtx_list[[i]])]),collapse=" ")
        cls1 = paste(ncol(mtx_list[[i]])-1,"2","1",sep=" ")
        cls2 = paste("#",case_group,control_group,sep=" ")
        cls = c(cls1,cls2,cls3)
        cls_file = file.path(out_dir,paste(names(mtx_list)[i],"_",case_group,"-vs-",control_group,".cls",sep=""))
        print(cls_file)
        write.table(cls, cls_file, row.names=F, quote =F, col.names=F)
        for( gmt_file in gmts ){
            setwd(out_dir)
            system(glue::glue("module purge && /gpfs/oe-scrna/guopengyu/GSEA/oebio2 gsea gsea -c {cls_file} -g {gmt_file} {mtx_file} --min_size {opts$minsize} --max_size {opts$maxsize} -p {opts$process} "))
        }
    })
}else{
    mtx = cbind(as.matrix(GetAssayData(seurat_ob)[,seurat_ob@meta.data[,col]==case_group]),
                as.matrix(GetAssayData(seurat_ob)[,seurat_ob@meta.data[,col]==control_group]))
    mtx = as.data.frame(mtx) %>% rownames_to_column(var="genes")
    out_dir = output_dir
    if (!file.exists(out_dir)){
        dir.create(out_dir, recursive = T)
    }
    mtx_file = file.path(out_dir, paste(case_group,"-vs-",control_group,".txt",sep=""))
    # write.table(mtx, mtx_file, sep = "\t", quote = FALSE, row.names=F)
    fwrite(mtx, file = mtx_file, sep = "\t", quote = FALSE, row.names=F)
    cls3_case_group = rep(case_group,times=length(seurat_ob@meta.data[which(seurat_ob@meta.data[,col]==case_group),][,col]))
    cls3_control_group = rep(control_group,times=length(seurat_ob@meta.data[which(seurat_ob@meta.data[,col]==control_group),][,col]))
    cls3 = c(cls3_case_group,cls3_control_group)
    cls3 =  paste0(cls3,collapse=" ")
    cls1 = paste(ncol(mtx)-1,"2","1",sep=" ")
    cls2 = paste("#", case_group, control_group, sep=" ")
    cls = c(cls1,cls2,cls3)
    cls_file = file.path(out_dir,paste(case_group,"-vs-",control_group,".cls",sep=""))
    write.table(cls, cls_file, row.names=F,quote =F,col.names=F)

    # 导出表达量不全为0的基因列表
    # gene_exp = FetchData(seurat_ob,vars=rownames(seurat_ob),slot="data")
    # gene_list = vector()
    # for(i in colnames(gene_exp)) {
    #     if(max(gene_exp[[i]]) > 0) {
    #         gene_list = c(gene_list, i)
    #     }
    # }
    gene_list = rownames(seurat_ob)[which(rowSums(seurat_ob) >0)]
    gene_df = data.frame(gene= gene_list)
    write.csv(gene_df, paste0(out_dir, "/genelist.txt"), quote=F, row.names=F)

    for( gmt_file in gmts ){
        gmt_file2 = normalizePath(gmt_file)
        setwd(out_dir)
        system(glue::glue("module purge && /gpfs/oe-scrna/guopengyu/GSEA/oebio2 gsea gsea -c {cls_file} -g {gmt_file2} {mtx_file} --min_size {opts$minsize} --max_size {opts$maxsize} -p {opts$process} "))

        ## 过滤gmt文件，剔除表达量均为0的gene
        filterpy <-"/gpfs/oe-scrna/guopengyu/script/filter_gmt.py"
        out_file =  ifelse (length(grep("kegg", basename(gmt_file))) > 0, "gene_gsea_kegg_background.xls", "gene_gsea_gobp_background.xls")
        print("开始过滤gmt文件，仅保留表达量不全为0的基因")
        cmd <- glue::glue('/gpfs/oe-scrna/guopengyu/conda/miniforge3/bin/python {filterpy} -gmt {gmt_file2} ',
                                    '-o {paste0(out_dir,"/", out_file)} -gene {paste0(out_dir, "/genelist.txt")}')
        futile.logger::flog.info(glue::glue("running:{cmd}"))
        system(cmd)
        }
    
    file.remove(paste0(out_dir, "/genelist.txt"))
}


