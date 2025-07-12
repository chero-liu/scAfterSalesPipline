### 调色板
CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}
# source /home/liuhongyan/miniconda3/bin/activate scvdj
### https://github.com/IndrajeetPatil/ggstatsplot 
suppressWarnings({
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("ggstatsplot"))
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("optparse",lib.loc="/home/dongjiaoyang/miniconda3/envs/OESingleCell/lib/R/library/"))
    suppressPackageStartupMessages(library("stringr"))
})
option_list = list(
    make_option( c("--input","-i"), type ="character",
                help="the list file of geneset ."),
    make_option( c("--output","-o"),type="character", default = "./",
                help="the output directory of results.", metavar="character"),
    make_option( c("--groupby", "-g"), type = "character", default = NULL,
                help = "[OPTIONAL]The grouppinig variable in the metadata for
                        separate the cells to visulize marker genes."),
    make_option( c("--centrality", "-c"),type= "logical",default = FALSE,
                help = "[OPTIONAL]Whether to show the centrality and the number of cells."),
    make_option( c("--vismethod","-m"), type= "character",default="geneset",
                help = "the visulization methods for the marker genes of each cell cluster, can only used geneset or ggstatsplot"),
    make_option( c("--pvalue", "-p"),type= "character",default = NULL,
                help = "[OPTIONAL]use like GROUP1:GROUP2+GROUP2:GROUP3 to add pvalue on vlnplot or boxplot.")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
options(stringsAsFactors = F)
data = read.delim( opt$input, sep="\t", header = T,quote="")
#保持groupby原有levels
data = as.data.frame(data)
groupby= as.character(opt$groupby)
#data[,groupby]=factor(data[,groupby],levels = c(unique(data[,groupby])))


if ( is.null(opt$vismethod) ){
    print("The vlnplot and featureplot will be used 'geneset'")
    vismethods = "geneset"
}else if ( !is.null(opt$vismethod)){
    print("The vlnplot and featureplot will be used 'ggstatsplot'")
    vismethods = opt$vismethod
}else{
  print("check!")
}


if ( is.null(opt$output) ){
    print("NO output directory specified,the current directory will be used!")
    output_dir = getwd()
}else{
    if ( file.exists(opt$output)){
        output_dir = opt$output
    }else{
        output_dir = opt$output
        dir.create(output_dir, recursive = T)
    }
}
groupby= as.character(opt$groupby)
comp = list()
if ( !is.null(opt$pvalue) ){
  all_comparisions = list()
  contrasts_list = unlist(strsplit(opt$pvalue, "\\+", perl = T))
  # contrasts_list = unlist(strsplit(pvalue, "\\+", perl = T)) ##delete
  for ( contrast in contrasts_list ){
    contrast = paste0(groupby,":",contrast)
    contrasts = unlist( strsplit(contrast,":",perl = T) )
    assay_metadata=data
    all_levels = as.vector(unique(assay_metadata[,contrasts[1]]))
    if ( contrasts[2] == "all" & contrasts[3] != "all" ){
        all_levels = all_levels[-which(all_levels==contrasts[3])] #delete the reference level
        all_comparisions = paste(all_levels,contrasts[3],sep = ":")
        break
    }else if( contrasts[2] == "all" & contrasts[3] == "all" ){
        all_comparisions = "all"
        break
    }else if ( contrasts[2] != "all" & contrasts[3] == "all" ){
        all_levels = all_levels[-which(all_levels==contrasts[2])]
        all_comparisions = paste(contrasts[2],all_levels,sep = ":")
        break
    }else{
        if ( !contrasts[2] %in% all_levels | !contrasts[3] %in% all_levels){
          print(paste0(contrasts[2],":",contrasts[3],"所选分组中细胞数为0,请检查分组比较信息。已跳过该分组。"))
        }else if ( table(assay_metadata[,groupby])[contrasts[2]]<=1 | table(assay_metadata[,groupby])[contrasts[3]]<=1){
          print(paste0(contrasts[2],":",contrasts[3],"所选分组中细胞数小于2,请检查分组比较信息。已跳过该分组。"))
        }else{
          all_comparisions = c(all_comparisions , paste0(contrasts[2],":",contrasts[3]))
        }
    }
  }
  my_comparisons = sort(unique((data[,groupby]))) ##找到分组信息
  comp=list()
  ##若为all:all,则通过循环生成各个分组两两结合比对的list
  if ( all_comparisions == "all" ){
    for(a in 1:(length(my_comparisons)-1)){
      for(b in 1:(length(my_comparisons)-a)){
      list=c(as.character(my_comparisons[a]),as.character(my_comparisons[a+b]))
      comp=c(comp,list(list))
      }
    }
  }else{
    for( i in all_comparisions ){
      list=strsplit(i,":",perl = T)
      comp=c(comp,list)
    }
  }
}
colors2use = CustomCol2(1:length(unique(data[,groupby])))


#增加判断，如果是ggstatsplot,保存结果
if ( vismethods == "ggstatsplot" ){
  for (i in colnames(data)[-(1:2)]){
      p = ggstatsplot::ggbetweenstats(data,x = !!sym(groupby),y = !!sym(i) ,
                                          plot.type = "boxviolin",
                                          results.subtitle =FALSE,
                                          messages = FALSE,
                                          pairwise.comparisons =FALSE, 
                                          mean.label.size = 0,
                                          centrality.plotting = opt$centrality,
                                          ylab = paste(i)) + 
              scale_color_manual(values= colors2use) +
              theme(axis.text.x = element_text(size=8,colour="black",angle = 30,vjust = 0.85,hjust = 0.75),
                    axis.text.y = element_text(size=8,colour="black"),
                    panel.grid.major =element_blank(), 
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(colour = "black"))
      savenames = paste0(i,"_ggstatsplot_violin_boxplot")
      if ( !is.null(opt$pvalue) & length(comp) != 0 ){
          suppressPackageStartupMessages( library("ggsignif") )
          p = p + geom_signif(comparisons = comp, map_signif_level = F, test ="wilcox.test", step_increase = 0.1)
          savenames = paste0(i,"_ggstatsplot_violin_boxplot_P_value")
      }
      ggsave(file.path(output_dir,paste0(savenames,".pdf")),plot = p,bg="white")
      ggsave(file.path(output_dir,paste0(savenames,".png")),plot = p,bg="white")
  }
}

#增加判断，如果是geneset,保存结果
if ( vismethods == "geneset" ){
    for (i in colnames(data)[-(1:2)]){
    p = ggstatsplot::ggbetweenstats(data,x = !!sym(groupby),y = !!sym(i) ,
                                        plot.type = "boxviolin",
                                        results.subtitle =FALSE,
                                        messages = FALSE,
                                        pairwise.comparisons =FALSE, 
                                        mean.label.size = 0,
                                        centrality.plotting = opt$centrality,
                                        ylab = paste(i,"_Score")) + 
            scale_color_manual(values= colors2use) +
            theme(axis.text.x = element_text(size=8,colour="black",angle = 30,vjust = 0.85,hjust = 0.75),
                  axis.text.y = element_text(size=8,colour="black"),
                  panel.grid.major =element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"))
    savenames = paste0(i,"_score_violin_boxplot")
    if ( !is.null(opt$pvalue) & length(comp) != 0 ){
        suppressPackageStartupMessages( library("ggsignif") )
        p = p + geom_signif(comparisons = comp, map_signif_level = F, test ="wilcox.test", step_increase = 0.1)
        savenames = paste0(i,"_score_violin_boxplot_P_value")
    }
    ggsave(file.path(output_dir,paste0(savenames,".pdf")),plot = p,bg="white")
    ggsave(file.path(output_dir,paste0(savenames,".png")),plot = p,bg="white")
  }
}


