#! /bin/bash

workdir=$1
infile=$2
TF=$3
DirForScriptSelf=$(cd "$(dirname "$0")";pwd)

# n=${infile%-diff-*-*-FC-[0-9].[0-9].ppi_network.tsv}
# n=${n#group_}
p=${infile%.tsv}

mkdir -p $workdir/$p
awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $2, $4, $5}' $workdir/$infile > $workdir/$p/gene2gene_network_tmp.txt
infile=$p/gene2gene_network_tmp.txt
#如果差异没有生成空文件
if [ $(cat $workdir/$infile | wc -l) = 1 ];then
  touch $workdir/$p/$p.html
  exit
fi


awk 'BEGIN{FS=OFS="\t"}NR>1{print $1}' $workdir/$infile >> $workdir/$p/tmp1
awk 'BEGIN{FS=OFS="\t"}NR>1{print $3}' $workdir/$infile >> $workdir/$p/tmp1
sort $workdir/$p/tmp1 |uniq -c|sed 's/ /\t/g'| awk 'BEGIN{FS=OFS="\t"}{print $(NF-1),$NF}' |awk '{print $0"\t"(FNR-1)}'|awk 'BEGIN{FS=OFS="\t"}{print $2,$3,$1}'  > $workdir/$p/network.tmp && rm $workdir/$p/tmp1
awk 'BEGIN{FS=OFS="\t"}NR>1{print $1,$2}' $workdir/$infile >> $workdir/$p/tmp2
awk 'BEGIN{FS=OFS="\t"}NR>1{print $3,$4}' $workdir/$infile >> $workdir/$p/tmp2
###inside color
sort -u $workdir/$p/tmp2 | awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=$1;b[$1]=$2;c[$1]=$3;next}{if ($1==a[$1]) print $0"_TFs";else print $0}' $TF - |awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=$1;b[$1]=$2;next}{if (a[$1]==$1) print $0,b[$1]}' - $workdir/$p/network.tmp |sed 1'i\name\tid\tsize\tgroup' > $workdir/$p/col  && rm $workdir/$p/tmp2
###relationship
awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=$1;b[$1]=$2;c[$1]=$3;next}{if (a[$1]==$1) print $3,b[$1]}' $workdir/$p/network.tmp $workdir/$infile |awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=$1;b[$1]=$2;c[$1]=$3;next}{if (a[$1]==$1) print $2,b[$1],"1","black"}' $workdir/$p/network.tmp - | sed 1'i\prot1\tprot2\tscore\tcol'> $workdir/$p/line && rm $workdir/$p/network.tmp
###networkD3
/data/software/conda_envs/scrna_envs/PPI_loc/bin/Rscript ${DirForScriptSelf}/network_3d.r -l $workdir/$p/line -n $workdir/$p/col -o $workdir/$p.html
# rm $workdir/$p/line $workdir/$p/col
rm -r $workdir/$p\_files
# rm $workdir/$infile
rm -rf $workdir/$p
