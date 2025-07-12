#=================================================================================
# customized function definition
#=================================================================================
filter_switchgenes = function (sce, allgenes = FALSE, pathway_name = NULL, genelists = GeneSwitches:::gs_genelists,
    genetype = c("Surface proteins", "TFs"), zero_pct = 0.9,
    r2cutoff = 0.03, direction = c("up", "down"), topnum = 1e+05)
{
    if (allgenes == TRUE) {
        toplotgl <- rowData(sce)
        toplotgl$feature_type <- "All genes"
    }
    else if (!is.null(pathway_name)) {
        gl <- c()
        for (pn in pathway_name) {
            pgl <- data.frame(msigdb_h_c2_c5[pn], pn, stringsAsFactors = FALSE)
            colnames(pgl) <- c("feature_name", "feature_type")
            gl <- rbind(gl, pgl)
        }
        multi <- gl$feature_name[duplicated(gl$feature_name)]
        if (length(multi) > 0) {
            gl <- ddply(gl, .(feature_name), paste)[, c(1, 3)]
            rownames(gl) <- gl$feature_name
            colnames(gl)[2] <- "types"
            gl$feature_type <- gl$types
            gl[multi, ]$feature_type <- "Multiple"
            genestoplot <- intersect(rownames(sce), gl$feature_name)
            toplotgl <- rowData(sce)[genestoplot, ]
                        toplotgl$feature_type <- gl[genestoplot, ]$feature_type
            toplotgl$types <- gl[genestoplot, ]$types
        }
        else {
            rownames(gl) <- gl$feature_name
            genestoplot <- intersect(rownames(sce), gl$feature_name)
            toplotgl <- rowData(sce)[genestoplot, ]
            toplotgl$feature_type <- gl[genestoplot, ]$feature_type
        }
    }
    else {
        genelists_sub <- genelists[genelists$genetypes %in% genetype,
            ]
        genelists_sub <- genelists_sub[!duplicated(genelists_sub$genenames),
            ]
        rownames(genelists_sub) <- genelists_sub$genenames
        genestoplot <- intersect(rownames(sce), genelists_sub$genenames)
        toplotgl <- rowData(sce)[genestoplot, ]
        toplotgl$feature_type <- genelists_sub[genestoplot, ]$genetypes
    }
    toplotgl_sub <- toplotgl[toplotgl$zerop_gene < zero_pct &
        toplotgl$prd_quality == 1 & toplotgl$pseudoR2s > r2cutoff &
        toplotgl$direction %in% direction ,]
    toplotgl_sub <- toplotgl_sub[grep("^(RPS|RPL)",toplotgl_sub$geneID,invert=T) , ]
    if (nrow(toplotgl_sub) > topnum) {
        toplotgl_sub <- toplotgl_sub[order(toplotgl_sub$pseudoR2s,
            decreasing = TRUE), ]
        toplotgl_sub <- toplotgl_sub[1:topnum, ]
    }
    return(toplotgl_sub)
}

plot_pathway_density=function (switch_pw_re, toplotgl_sig, pw_direction = c("up","down"), orderbytime = TRUE)
{
    switch_pw_re <- switch_pw_re[switch_pw_re$direction %in%                                          
        pw_direction, ]                                 
    if (orderbytime == TRUE) {                                                                         
        switch_pw_re <- switch_pw_re[order(switch_pw_re$switch_at_time,
            decreasing = TRUE), ]
    }
    else {
        switch_pw_re <- switch_pw_re[order(switch_pw_re$FDR,
            decreasing = TRUE), ]
    }
    gl <- c()
    for (i in 1:nrow(switch_pw_re)) {
        pn <- rownames(switch_pw_re)[i]
        pgl <- data.frame(gset_list[pn], pn, stringsAsFactors = FALSE)
        colnames(pgl) <- c("Genes", "Pathways")
        tgn <- length(unique(pgl$Genes))
        genestoplot <- intersect(rownames(toplotgl_sig), pgl$Genes)
        pgl <- pgl[pgl$Genes %in% genestoplot, ]
        toplotgl <- toplotgl_sig[pgl$Genes, ]
        pgl <- as.data.frame(cbind(pgl, toplotgl))
        pgl <- pgl[pgl$direction == switch_pw_re[i, ]$direction,
            ]
        pgl$Pathways <- paste0(pgl$Pathways, "(", nrow(pgl),
            "/", tgn, ")")
        gl <- rbind(gl, pgl)
    }
    gl$Pathways <- factor(gl$Pathways, levels = unique(gl$Pathways))
    p <- ggplot(gl, aes(x = switch_at_time, y = Pathways, fill = direction,
        col = direction))
    p <- p + theme_classic()
    p <- p + xlab("Pseudo-timeline") + geom_density_ridges(alpha = 0.5) +
        labs(fill = "Regulation", color = "Regulation")
    p <- p + theme(text = element_text(size = 12, family = "Helvetica"),
        panel.background = element_rect(fill = "white", colour = NA),
        axis.line = element_line(colour = "black"), legend.key.size = unit(10,
            "pt"), legend.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 11, colour = "black")) +
        scale_color_manual(values = c("forestgreen", "chocolate2")) +    
        scale_fill_manual(values = c("forestgreen", "chocolate2"))
    return(p)                                                           
}  

distinct_genes=
function (toplotgl_Rsub1, toplotgl_Rsub2, path1name = "Path1Genes",                                   
    path2name = "Path2Genes", r2cutoff = 0.05, scale_timeline = FALSE,                                
    path1time = NULL, path2time = NULL, bin = 100, topn=NULL)                                             
{
    toplotgl_Rsub1$genenames <- rownames(toplotgl_Rsub1)                                              
    toplotgl_Rsub2$genenames <- rownames(toplotgl_Rsub2)                                              
    gl1 <- toplotgl_Rsub1$genenames
    gl2 <- toplotgl_Rsub2$genenames
    glin1 <- setdiff(gl1, gl2)
    glin2 <- setdiff(gl2, gl1)
    gs_p1 <- toplotgl_Rsub1[glin1, ]
    gs_p2 <- toplotgl_Rsub2[glin2, ]
    gs_p1$Paths <- path1name
    gs_p2$Paths <- path2name
    if (scale_timeline == TRUE) {        
steptime1 <- (max(path1time) - min(path1time))/bin                                            
gs_p1$switch_at_time <- round((gs_p1$switch_at_time - min(path1time))/steptime1)      
steptime2 <- (max(path2time) - min(path2time))/bin                                                    
gs_p2$switch_at_time <- round((gs_p2$switch_at_time - min(path2time))/steptime2)    }
    toplotgl <- rbind(gs_p1[, c("geneID", "zerop_gene", "switch_at_time",  "pvalues", "FDR", "pseudoR2s", "estimates", "prd_quality","direction", "switch_at_timeidx", "genenames", "Paths")],gs_p2[, c("geneID", "zerop_gene", "switch_at_time", "pvalues","FDR", "pseudoR2s", "estimates", "prd_quality", "direction","switch_at_timeidx", "genenames", "Paths")])
    if (!is.null(topn)) {
       toplotgl = toplotgl[head(order(toplotgl $pseudoR2s,decreasing=T),topn),]
       print(paste0("plotting  top " ,topn," switching genes."))
    } else {
        toplotgl <- toplotgl[toplotgl$pseudoR2s > r2cutoff, ]
        print(paste0("plotting  switching genes with r2 > " ,r2cutoff, " (",nrow(toplotgl ),"genes)"))
    }
    return(toplotgl)
}

#=================================================================================
# command line parameters setting
#=================================================================================
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(monocle))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(OESingleCell))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(GeneSwitches,lib.loc="/home/dongjiaoyang/bin/Rpackages/"))

option_list = list(
    make_option( c("--seurat", "-i"), type = "character", 
                 help = "the seurat object saved as R object in RDS format."),
    make_option( c("--upstream_obj", "-m"), type = "character", 
                 help = "the monocle or slingshot object saved as R object in RDS format."),
    make_option( c("--upstream", "-f"), type = "character",  default = "monocle",
                 help = "上游依赖任务，monocle或者slingshot."),
    make_option( c("--output","-o"),type="character", default = "./",
                help="the output directory of results."),
    make_option( c("--states","-s"),type="character", default = NULL,
                help="states for the analysis. Default: automatically taking all states. To specify two branches, use colon to seperate two branches, prebranch is needed."),
    make_option( c("--bn_cutoff","-b"),type="double", default = 0.05,                                 
                help="binarize cut-offs. choose appropriate cutoff based on histogram if needed."),
    make_option( c("--reference","-r"),type="character",                                 
                help="directory containing TFs.txt or enrichment gmt files. OR mouse human."),
    make_option( c("--showtf"),type="logical", default = TRUE,                                 
                help="whether to show TFs on the timeline."),
    make_option( c("--topn","-n"),type="integer", default = 50,                           
                help="Number of top genes to show on the timeline.")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#=================================================================================
#parse the command line parameters
#=================================================================================
if ( is.null(opt$output) ){
    print("NO output directory specified,the current directory will be used!")
    output_dir = getwd()
}else{
    if ( file.exists(opt$output)){
        output_dir = opt$output
    }else{
        output_dir = opt$output
        dir.create(output_dir, recursive = T)
        dir.create(paste0(output_dir,"/enrichment"), recursive = T)
    }
}
upstream_obj=readRDS(opt$upstream_obj)
if (opt$upstream == "slingshot"){
    if (is.null(opt$states)){
        branch = as.list(names(upstream_obj@int_metadata$slingshot@curves))
        num=length(branch)
    } else {
        states = unlist( strsplit(opt$states,":",perl = T) )
        num = length(states)
        branch = list()
        for ( i in 1:num) {
            branch[[i]]= names(upstream_obj@int_metadata$slingshot@curves)[as.integer(states[i])]
        }
    }
} else {
    seurat_ob=readRDSMC(opt$seurat,20)
    if ( seurat_ob@version < 3){
            seurat_ob = UpdateSeuratObject(seurat_ob) #make sure the seurat object match with the latest seurat package
    }
    logexpdata=as.matrix(GetAssayData(seurat_ob))

    if (is.null(opt$states)){
        num=1
        branch=list(as.integer(unique(phenoData(upstream_obj)@data[order(phenoData(upstream_obj)$Pseudotime),"State"])))
    } else {
        states= unlist( strsplit(opt$states,":",perl = T) )
        num=length(states)
        branch=list()
        for ( i in 1:num) {
            branch[[i]]= as.integer(unlist( strsplit(states[i],",",perl = T) ))
        }
    }
}
# new_cds <- buildBranchCellDataSet(cds_subset, branch_point = branchpoint, progenitor_method = "duplicate")
# cell_fate1 = unique(pData(new_cds)[which(pData(new_cds)$Branch == unique(pData(new_cds)$Branch)[1]),]$State)#cell_fate2
# cell_fate2 = unique(pData(new_cds)[which(pData(new_cds)$Branch == unique(pData(new_cds)$Branch)[2]),]$State)
# branch_labels = c( paste("State",paste(sort(setdiff(cell_fate1,cell_fate2)),collapse="-")) ,  paste("State",paste(sort(setdiff(cell_fate2,cell_fate1)),collapse="-")) )

########## refs
if (tolower(opt$reference)=="human"){
    ref_dir= "/data/database/geneswiches_database/human"
} else if (tolower(opt$reference)=="mouse"){
    ref_dir= "/data/database/geneswiches_database/mouse"
} else {
    ref_dir= opt$reference
}

#1.TF
if (opt$showtf==T){
    TFlist=read.table(file.path(ref_dir,"TFs.xls"),sep="\t",header= T)
    rownames(TFlist)=TFlist$genenames
}
#2. enrich
gmt_file=c(
kegg=file.path(ref_dir,"kegg.gmt"),
go_bp=file.path(ref_dir,"go_bp.gmt"),
go_cc=file.path(ref_dir,"go_cc.gmt"),
go_mf=file.path(ref_dir,"go_mf.gmt"),
go=file.path(ref_dir,"go.gmt"))


#=================================================================================
# main
#=================================================================================
sce=list()
sg_pw=list()
for ( i in 1:num){
    ###### prepare
    if (opt$upstream == "slingshot"){
        curve_num = gsub("curve","",branch[[i]])
        sce[[i]] = convert_slingshot(upstream_obj, paste0("slingPseudotime_",curve_num), assayname = names(assays(upstream_obj))[2])
        branch_name = branch[[i]]
    } else {
        sce[[i]] <- convert_monocle2(monocle2_obj = upstream_obj, 
                               states = branch[[i]] , expdata = logexpdata)
        branch_name = paste0("branch_State",paste0(branch[[i]],collapse="-"))
    }
    pdf(file.path(output_dir,paste0(branch_name,".hist.pdf")))
    h <- hist(as.matrix(assays(sce[[i]])$expdata), breaks = 200, plot = FALSE)
    {plot(h, freq = FALSE, xlim = c(0,2), ylim = c(0,1), main = "Histogram of gene expression",
    xlab = "Gene expression", col = "darkgoldenrod2", border = "grey")
    abline(v=0.05, col="blue")}
    dev.off()
    ###### binarize & fit glm & save raw result
    sce[[i]] <- binarize_exp(sce[[i]], fix_cutoff = TRUE, binarize_cutoff = opt$bn_cutoff)
    sce[[i]] <- find_switch_logistic_fastglm(sce[[i]], downsample = TRUE, show_warning = FALSE)
    saveRDS(sce,file.path(output_dir,"sce_result.rds"))
    
    
    all_data=rowData(sce[[i]]) %>% as.data.frame %>% arrange(-pseudoR2s) %>% filter(prd_quality==1) %>%  dplyr::rename(GeneID = geneID) %>% select(GeneID,zerop_gene,switch_at_time,pseudoR2s,direction,FDR,switch_at_timeidx)
    write.table(all_data ,file.path(output_dir,paste0("all_switching_genes_for_",branch_name,".xls")),quote=F, sep="\t",row.names=F)       
    sig_data=all_data %>% filter(zerop_gene < 0.9 ,pseudoR2s > 0.2)
    write.table(sig_data ,file.path(output_dir,paste0("sig_switching_genes_for_",branch_name,"_R2-0.2.xls")),quote=F, sep="\t",row.names=F)       
    
    
    ###### visualization
    sg_allgenes <- filter_switchgenes(sce[[i]], allgenes = TRUE, topnum = opt$topn)
    
    if (opt$showtf==T){
        ## top 15 surface proteins and TFs only
        sg_gtypes <- filter_switchgenes(sce[[i]], allgenes = FALSE, topnum = 15,
                                        genelists = TFlist , genetype = c("Surface proteins", "TFs")) 
        
        ## combine switching genes and remove duplicated genes from sg_allgenes
        sg_vis <- rbind(sg_gtypes, sg_allgenes[setdiff(rownames(sg_allgenes), rownames(sg_gtypes)),])
    } else {
        sg_vis <-  sg_allgenes
    }
    timeline = plot_timeline_ggplot(sg_vis, timedata = sce[[i]]$Pseudotime, txtsize = 3)
    ggsave(file.path(output_dir,paste0("timeline_for_",branch_name,".pdf")),timeline,bg="white")
    ggsave(file.path(output_dir,paste0("timeline_for_",branch_name,".png")),timeline,dpi=1000,bg="white")
    
    ###### enrichment
    # prep: filter sig [feature meta]
    sg_pw[[i]] <- filter_switchgenes(sce[[i]],   allgenes = TRUE, r2cutoff = 0) # use all switching genes
    
    for ( gmt  in names(gmt_file)){
        if(!file.exists(gmt_file[[gmt]])) next 
        gset_list =GSEABase::geneIds(GSEABase::getGmt(con=gmt_file[gmt]))
        # enrich
        switch_pw <- find_switch_pathway(rowData(sce[[i]]), sig_FDR = 0.05, pathways = gset_list, toplotgl_sig=sg_pw[[i]])
        # filter
        skip_to_next <- FALSE
        switch_pw_reduce=tryCatch({
            reduce_pathways(switch_pw, pathways = gset_list, redundant_pw_rate = 0.8)
        }, error = function(e) {
            message(paste0("An error occurred in ",gmt,".gmt analysis process ,but it has been caught and ignored"))
            skip_to_next <<- TRUE
        })
        if(skip_to_next) next

        if (nrow(switch_pw_reduce) > 0) {
            pathway = plot_pathway_density(switch_pw_reduce[1:10,], sg_pw[[i]], orderbytime = TRUE)
            ggsave(file.path(output_dir,paste0("enrichment/",gmt,"_top10_significant_pathways_for_",branch_name,".pdf")),pathway,width=14,bg="white")
            ggsave(file.path(output_dir,paste0("enrichment/",gmt,"_top10_significant_pathways_for_",branch_name,".png")),pathway,width=14,bg="white")
        } else {
            warning(paste0("No significnat pathway for ",gmt," in branch",branch[[i]], "\nNum of sig switching genes:",nrow(sg_pw[[i]])))
        }
    }
    sg_pw[[i]] <- filter_switchgenes(sce[[i]], allgenes = TRUE, r2cutoff = 0.03) # for next section
}

# if two branches ======================================================================
if (num == 2 ) {
    sg_disgs <- distinct_genes(sg_pw[[1]], sg_pw[[2]], 
      #path1name = paste0("S",paste0(branch[[1]],collapse="-")),path2name =  paste0("S",paste0(branch[[2]],collapse="-")),
      path1time = sce[[1]]$Pseudotime, path2time = sce[[2]]$Pseudotime,
      topn=30)
    
    diff= plot_timeline_ggplot(sg_disgs, timedata = sce[[1]]$Pseudotime, color_by = "Paths", 
                         iffulltml = FALSE, txtsize = 3)
    ggsave(file.path(output_dir,"distinct_switching_genes_between_two_branches.pdf"),diff,bg="white")
    ggsave(file.path(output_dir,"distinct_switching_genes_between_two_branches.png"),diff,bg="white")
    
    
    ## scales
    sg_disgs_scaled <- distinct_genes(sg_pw[[1]], sg_pw[[2]], 
      #path1name = paste0("S",paste0(branch[[1]],collapse="-")),path2name =  paste0("S",paste0(branch[[2]],collapse="-")),
      path1time = sce[[1]]$Pseudotime, path2time = sce[[2]]$Pseudotime,
      topn=30,
    #  r2cutoff =0.15,
      scale_timeline = T, bin = 100) 
    
    diff_scaled= plot_timeline_ggplot(sg_disgs_scaled, timedata = 1:100, color_by = "Paths", 
                         iffulltml = FALSE, txtsize = 3)
    ggsave(file.path(output_dir,"distinct_switching_genes_between_two_branches_scaled.pdf"),diff_scaled,bg="white")
    ggsave(file.path(output_dir,"distinct_switching_genes_between_two_branches_scaled.png"),diff_scaled,bg="white")
}

if (!file.exists(file.path(output_dir, "GeneSwitches基因开关分析说明.docx"))) {
	file.copy( "/public/scRNA_works/Documents/GeneSwitches基因开关分析说明.docx", file.path(output_dir, "GeneSwitches基因开关分析说明.docx"))
}

