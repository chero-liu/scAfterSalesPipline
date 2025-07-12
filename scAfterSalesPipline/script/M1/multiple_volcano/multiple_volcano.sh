module purge && module load OESingleCell/3.0.d


# library(stringr)
# cts = c("Macro1,Macro2,Macro3,Macro4,Macro5,Macro6,Macro7,Macro8,Macro9,Mono_macro1,Mono_macro2")
# cts=str_split(cts,",")[[1]]
# dflist = list()
# for (ct in cts){
#     input = paste0("/gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2023110250/result/M1/Mono_Macrophage_DC/SubDiff/",ct,"/STB-vs-BS/group_STB-vs-BS-all_diffexp_genes.xls")
#     data = read.csv(input,sep = "\t",header = T)
#     data$cluster = ct
#     data = data[,c("gene","FoldChange","p.value","cluster")]
#     colnames(data) = c("gene","FC","p_val","cluster")
#     dflist[[ct]] = data
# }
# df = Reduce(rbind,dflist)
# write.table(df,"/gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2023110250/result/M1/Mono_Macrophage_DC/SubDiff/all_diffexp_genes.file",row.names = F,quote = F,sep = "\t")



Rscript /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/script/M1/multiple_volcano/multiple_volcano.r \
    -u /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2023110250/result/M1/CD4T/multiple_volcano/all_diffexp_genes.file \
    -o /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2023110250/result/M1/CD4T/multiple_volcano \
    -n 10 \
    -p 0.05 \
    -f 1.5

