###使用方式,数据输入有两种，1.上传的输入数据列名为gene,p_val,FC,cluster列的tab分隔文件；2.上游依赖任务模块为findallmarkers的结果文件all_markers_for_each_cluster_anno.tsv
#Rscirpt multiple_volcano.R -i xxx.file -o ./ -n NULL -p 0.05 -f 1.5
#
library(tidyverse)
library(ggrepel) 
library(dplyr)
suppressWarnings({
  suppressPackageStartupMessages( library("optparse") )
  suppressPackageStartupMessages(library("tidyverse"))
  suppressPackageStartupMessages(library("ggplot2"))
  suppressPackageStartupMessages(library("ggrepel"))
})
CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[1:n])
}
#=command line parameters setting=============================
option_list = list(
    make_option( c("--input","-i"), type ="character",
               help="the list files of Differential gene to be visulized."),  
  make_option( c("--genelist","-u"), type ="character",
               help="the upload files of Differential gene to be visulized."),
  make_option( c("--output","-o"),type="character", default = "./",
               help="the output directory of results.", metavar="character"),
  make_option( c("--topn", "-n"), type="character", default = NULL,
               help = "the number of top Differential gene for each cluster to visualizse."),
  make_option( c("-p", "--pvalue"), type = "double", default = 0.05,
               help = "the Padj of the gene differential expression."),
    make_option( c("-f", "--FC"), type = "double", default = 1.5,
               help = "the Foldchange of the gene differential expression.")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(is.null(opt$input) & !is.null(opt$genelist)){
    df=read.csv(opt$genelist,sep="\t",header = T)
}else if (!is.null(opt$input)){
    df=read.csv(opt$input,sep="\t",header = T)
    df=df %>% mutate(FC=2**(avg_log2FC)) %>%mutate(log2FC=log2(FC))
}else{
    print("输入文件有误！！！")
}


# 自定义一个函数，之后我们直接调用它
mutiVolcano = function(df,         # 绘图数据
                       P = 0.05,   # P值卡值
                       FC = 1.5,   # FC卡值
                       topn=NULL,
                       GroupName = c("Sig","Not Sig"),      # 分组标签
                       pointColor = c("#CC3333","#0099CC"), # 分组散点的颜色
                       barFill = "#efefef",  # 柱子的颜色
                       pointSize = 0.9,      # 散点的大小
                       labeltype = "1",      # 标记差异基因的选项，标记类型有"1"和"2"两种选项
                       labelNum = 5,         # 当标记类型为1时，待标记的散点个数
                       labelName =NULL,      # 当标记类型为2时，待标记的散点名称
                       tileLabel =  "Label", # 标记比较对的选项，选项有“Label”和“Num”，Label时显示分组名称，Num时显示数字，防止因为标签太长导致的不美观
                       tileColor = NULL      # 比较对的颜色
                       ){
  # 数据分组 根据p的卡值分组
  dfSig = df %>% 
    mutate(log2FC = log2(FC)) %>%
    filter(FC>{{FC}} | FC <(1/{{FC}})) %>%
    mutate(Group = ifelse(p_val<0.05,GroupName[[1]],GroupName[[2]])) %>%
    mutate(Group = factor(Group,levels=GroupName)) %>%
    mutate(cluster = factor(cluster,levels=unique(cluster)))   # cluster的顺序是文件中出现的顺序
  
  # 柱形图数据整理
  dfBar = dfSig %>%
    group_by(cluster) %>%
    summarise(min = min(log2FC,na.rm = T),
              max = max(log2FC,na.rm = T)
              )
  # 散点图数据整理
  dfJitter = dfSig %>%
    mutate(jitter = jitter(as.numeric(cluster),factor = 2))
  dfJitter
  
  # 整理标记差异基因的数据
  if(labeltype == "1"){
    if(!is.null(topn)){
    # 标记一
    # 每组P值最小的几个
    dfLabel = dfJitter %>%
      group_by(cluster) %>%
      slice_min(p_val,n=labelNum,with_ties = F) %>%slice_head(n = topn)%>%
      ungroup()
    }else if(is.null(topn)){
      dfLabel = dfJitter %>%
      group_by(cluster) %>%
      slice_min(p_val,n=labelNum,with_ties = F) %>%
      ungroup()
    }
  }else if(labeltype == "2"){
    # 标记二
    # 指定标记
    dfLabel = dfJitter %>%
      filter(gene %in% labelName)
  }else{
    dfLabel = dfJitter %>% slice()
  }
    
  # 绘图
  p = ggplot()+
    # 绘制柱形图
    geom_col(data = dfBar,aes(x=cluster,y=max),fill = barFill)+
    geom_col(data = dfBar,aes(x=cluster,y=min),fill = barFill)+
    # 绘制散点图
    geom_point(data = dfJitter,
               aes(x = jitter, y = log2FC, color = Group),
               size = pointSize,
               show.legend = NA
               )+
    # 绘制中间的标签方块
    ggplot2::geom_tile(data = dfSig,
                       ggplot2::aes(x = cluster, y = 0, fill = cluster), 
                       color = "black",
                       height = log2(FC) * 1.5,
                       # alpha = 0.3,
                       show.legend = NA
                       ) + 
    # 标记差异基因
    ggrepel::geom_text_repel(
      data = dfLabel,
      aes(x = jitter,                   # geom_text_repel 标记函数
          y = log2FC,          
          label=gene),        
      min.segment.length = 0.1,
      max.overlaps = 10000,                    # 最大覆盖率，当点很多时，有些标记会被覆盖，调大该值则不被覆盖，反之。
      size=3,                                  # 字体大小
      box.padding=unit(0.5,'lines'),           # 标记的边距
      point.padding=unit(0.1, 'lines'), 
      segment.color='black',                   # 标记线条的颜色
      show.legend=F)#+
  
  if(tileLabel=="Label"){
     p =
      p +
      geom_text(data = dfSig,aes(x = cluster,y = 0,label = cluster))+
      ggplot2::scale_fill_manual(values = tileColor,
                                 guide = NULL # 不显示该图例
                                 )
  }else if(tileLabel=="Num"){
    # 如果比较对的名字太长，可以改成数字标签
    p =
      p +
      geom_text(data = dfSig,aes(x = cluster,y = 0,label = as.numeric(cluster)),show.legend = NA)+
      ggplot2::scale_fill_manual(values = tileColor,
                                 labels = c(paste0(1:length(unique(dfSig$cluster)),": ",unique(dfSig$cluster))))
  }

  
    
  # 修改主题
  p = p+ggplot2::scale_color_manual(values = pointColor)+
    theme_classic()+
    ggplot2::scale_y_continuous(n.breaks = 5) + 
    ggplot2::theme(
                   legend.position = "right", 
                   legend.title = ggplot2::element_blank(), 
                   legend.background = ggplot2::element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   # axis.title.x = element_blank(),
                   axis.line.x = element_blank()
                   ) + 
    ggplot2::xlab("Clusters") + ggplot2::ylab("log2FC") + 
    # ggplot2::guides(fill = ggplot2::guide_legend())
    guides(color = guide_legend(override.aes = list(size = 3)))
    
    return(p)
}

if(is.null(opt$topn)){
    topn=NULL
}else{
    topn=as.numeric(opt$topn)
}
# 调用函数画图
plot=mutiVolcano(
  df = df,    # 绘图数据
  P = opt$pvalue,   # P值卡值
  FC = opt$FC,   # FC卡值
  topn=topn, #是否只高亮topN
  GroupName = c(paste0("Pvalue<",opt$pvalue),paste0("Pvalue>",opt$pvalue)),      # 分组标签
  pointColor = c("#CC3333","grey"), # 分组散点的颜色
  barFill = "#efefef",   # 柱子的颜色
  pointSize = 0.9,       # 散点的大小
  labeltype = "1",       # 标记差异基因的选项，标记类型有"1"和"2"两种选项
  labelNum = 5,          # 当标记类型为1时，待标记的散点个数
  labelName =c("ID1","ID2029"),   # 当标记类型为2时，待标记的散点名称
  tileLabel =  "Label",           # 标记比较对的选项，选项有“Label”和“Num”，Label时显示分组名称，Num时显示数字，防止因为标签太长导致的不美观
  #tileColor = RColorBrewer::brewer.pal(length(unique(df$cluster)),"Set3")
  tileColor = CustomCol2(length(unique(df$cluster)))# 比较对的颜色
)

if (!dir.exists(opt$output)){
  dir.create(opt$output,recursive = T)
}

ggsave(plot,file=paste0(opt$output,"/multiple_volcano_plot.pdf"),width=30)
ggsave(plot,file=paste0(opt$output,"/multiple_volcano_plot.png"),width=30)
write.table(df,file=paste0(opt$output,"/multiple_volcano_plot_diffexp_gene_result.tsv"),sep="\t",row.names=F,quote=F)
