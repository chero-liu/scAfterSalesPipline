createDir <- function(dir){
    if(!dir.exists(dir)){
        dir.create(dir, recursive=T)
    }
}

#--------------------------------------------------------------------------------------------------------------------------
reorderLevelsByFrequency <- function(factor_var) {
  if(!is.factor(factor_var)) {
    factor_var = factor(factor_var)
  }
  freq <- table(factor_var)
  sorted_levels <- names(sort(freq, decreasing = TRUE))
  factor_var <- factor(factor_var, levels = sorted_levels)
  return(factor_var)
}

reorderClusters <- function(rds, colname) {
  rds@meta.data[[colname]] = reorderLevelsByFrequency(rds@meta.data[[colname]])
  old_clusters <- as.integer(rds@meta.data[[colname]])
  unique_clusters <- sort(unique(old_clusters))
  new_cluster_mapping <- setNames(seq_along(unique_clusters), unique_clusters)
  rds@meta.data[[colname]] <- new_cluster_mapping[as.character(old_clusters)]
  rds@meta.data[[colname]] <- factor(rds@meta.data[[colname]], levels = seq_along(unique_clusters))
  return(rds)
}

dirnameScript <- function(){
    cmd <- commandArgs(trailingOnly = FALSE)
    scriptName <- sub('--file=', "", cmd[grep('^--file=', cmd)])
    path <- system(paste0("realpath ", scriptName), intern=T)
    scriptDirname <- dirname(path)
    return(path)
}

readH5seurat <- function(input,informat = "h5seurat"){
    #env /gpfs/oe-scrna/liuchenglong/envs/h5seurat/bin/R
    rds = OESingleCell::ReadX(input = input, informat =  informat, verbose = F)
    return(rds)
}

saveH5seurat <- function(rds, outdir,update = FALSE,informat = "h5seurat"){
    #env /gpfs/oe-scrna/liuchenglong/envs/h5seurat/bin/R
    rds = OESingleCell::SaveX(rds, output = outdir, update = update, outformat = informat)
    return(rds)
}

getColorOrder = function(rds,useCol){
    colors_data <- rds@meta.data[!duplicated(rds@meta.data[,useCol]), ][,c(useCol, paste0(useCol,"_col"))]
    colors_data[,useCol] = factor(colors_data[,useCol],levels = levels(rds@meta.data[,useCol]))
    colors_data = colors_data %>% arrange(!!sym(useCol))
    colors2use = colors_data[,paste0(useCol,"_col")]
    return(colors2use)
}

homologene_seurat_replace <- function(data_ob,assay2use,inTaxidgene,outTaxidgene){
    counts=data_ob@assays[[assay2use]]@counts
    counts_filted = counts[as.vector(inTaxidgene),]
    rownames(counts_filted) = as.character(outTaxidgene)
    data_ob@assays[[assay2use]]@counts=counts_filted
    expr=data_ob@assays[[assay2use]]@data
    gene=rownames(expr)
    expr_filtered = expr[as.vector(inTaxidgene),]
    rownames(expr_filtered) = as.character(outTaxidgene)
    data_ob@assays[[assay2use]]@data=expr_filtered

    print("Since the change is changed to scale only for high-variable genes, scale.data is discarded from the homologous converted RDS")
    features=data_ob@assays[[assay2use]]@meta.features
    features_filtered = features[as.vector(inTaxidgene),]
    features_filtered = as.matrix(features_filtered)
    rownames(features_filtered) = as.character(outTaxidgene)
    features_filtered = as.data.frame(features_filtered)
    data_ob@assays[[assay2use]]@meta.features=features_filtered

    return(data_ob)
}

homologene_transformed = function(data_ob,inTaxid,outTaxid,output,assay="RNA"){
    suppressPackageStartupMessages(library("Seurat"))
    suppressPackageStartupMessages(library("homologene"))
    suppressPackageStartupMessages(library("OESingleCell"))
    assay2use = assay
    counts=data_ob@assays[[assay2use]]@counts
    gene=row.names(counts)
    in2out = homologene(gene,inTax=inTaxid,outTax=outTaxid)
    index <- duplicated(in2out[,1])
    in2out <-in2out[!index,]
    index<-duplicated(in2out[,2])
    in2out <- in2out[!index,]
    homologene <- in2out[,1:2]
    write.table(homologene,file.path(output,paste0("homologene.",inTaxid,"_2_",outTaxid,".xls",collapse=".")),sep="\t",row.names=F,quote=F)
    data_ob <- homologene_seurat_replace(data_ob,assay2use,in2out[,1],in2out[,2])
    
    return(data_ob)
}

create_group_columns <- function(data, cluster_col_name) {
  unique_clusters <- sort(unique(data[[cluster_col_name]]))
  for (cluster in unique_clusters) {
    group_col_name <- paste0("group_", cluster)
    data[[group_col_name]] <- ifelse(data[[cluster_col_name]] == cluster,
                                     as.character(cluster),
                                    #  paste0(setdiff(unique_clusters, cluster), collapse = "_"),
                                     "other")
  }
  return(data)
}