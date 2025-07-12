
# Description: draw venn graph
# Options:
#         -i     <infile>    The input files                                                       [Required]
#         -f     <str>       Inputfile format: matrix or sets. (default: matrix)                   [Optional]
#         -l     <str>       The label names is shown on the venn graph when '-f' is sets          [Optional]
#                            default: labels is inputfile's names when there is no option '-l'
#         -c     <str>       The inputfile specified columns is used draw venn                     [Optional]
#                            default: 1(for sets), all colnum(for matrix)
#         -o     <outdir>    The output directory of result. (default: ./)                         [Optional]
#         -m     <str>       Draw venn graph method: VennDiagram, UpSetR, Petals                   [Optional]
#                            default: VennDiagram(<=5), UpSetR(>5 && <=15), Petals(>15)
#         -s     <infile>    color config file for 'Petals'                                        [Optional]
#                            Without this option, the color is the default color
#                            2 colnum(labels[tab]color):
#                            Group2-vs-Group1     #8470FF
#                            Group3-vs-Group2     #20B2AA
#         -h|help            print help info
# Example:
#         perl venn_graph.pl -i file1.xls file2.xls file3.xls -f sets -l name1 name2 name3 -c 1 1 1 -o outdir/
#         perl venn_graph.pl -i matrix.xls -c 2 3 4 5 -o outdir/

/home/lipeng/miniconda3/bin/perl /public/scRNA_works/pipeline/scRNA-seq_further_analysis/venn_graph.pl \
    -i /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024051785/result/M1/Astrocyte_res0.05/Diffexp/Narcissin-vs-Model/group_Narcissin-vs-Model-all_diffexp_genes.xls /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024051785/result/M1/Endothelial_cell_res0.09/Diffexp/Narcissin-vs-Model/group_Narcissin-vs-Model-all_diffexp_genes.xls /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024051785/result/M1/Microglia_cell_res0.09/Diffexp/Narcissin-vs-Model/group_Narcissin-vs-Model-all_diffexp_genes.xls \
    -f sets \
    -l Astrocyte Endothelial_cell Microglia_cell \
    -c 1 1 1 \
    -o /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024051785/20241012/result/DZOE2024051785-唐增颖老师单细胞外泌体小鼠脑20241012/venn
