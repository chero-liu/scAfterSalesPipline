source /home/lipeng/miniconda3/bin/activate Seurat3.1 #激活环境
# Rscript scmetabolism_scFEA.v1.9.R -h
#     -v: ***.rds        #rds文件
#     -q: new_celltype   #分组信息
#     -c: redwhiteblue   #热图色板,默认海蓝-白-砖红
#     -n: 10             #差异筛选展示的top
#     -s:species         #human or mouse
#     -t: KEGG           #or REACTOME，or scFEA，对应不同数据集(默认提供KEGG和REACTOME结果)
#     -p: 0.05           #pvalue
#     -m: VISION         #其他算法依赖的环境待测试；
#     -d: contrast       #差异分组 不是all:all的数据，可配合predicate抽取部分数据绘图
#     --pair:            #显示P值
#     --predicate        #细胞筛选条件，表达式规则和R语法规则一致，如 --predicate "group %in% c('control')"  

Rscript /public/scRNA_works/pipeline/scRNA-seq_further_analysis/scmetabolism_scFEA.v1.9.R \
    -v /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2023113282/result/M1/B_Plasma/Clustering/data_ob_v3.rds  \
    -o /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2023113282/result/M2/scMetabolism/B_Plasma/clusters \
    -n 10 \
    -q clusters \
    -t KEGG \
    -p 0.05 \
    -m VISION \
    --cpu 4 \
    -s human \
    --pair True


Rscript /public/scRNA_works/pipeline/scRNA-seq_further_analysis/scmetabolism_scFEA.v1.9.R \
    -v /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2023113282/result/M1/B_Plasma/Clustering/data_ob_v3.rds  \
    -o /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2023113282/result/M2/scMetabolism/B_Plasma/clusters \
    -n 10 \
    -t REACTOME \
    -q clusters \
    -p 0.05 \
    -m VISION \
    --cpu 4 \
    -s human \
    --pair True
