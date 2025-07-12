suppressPackageStartupMessages(library("optparse"))

homologene_seurat_replace <- function(seurat_ob,assay2use,inTaxidgene,outTaxidgene){
   
    # for seurat object replacing homologene in type example counts, data 
    ##  counts gene replace 
    counts=seurat_ob@assays[[assay2use]]@counts
    counts_filted = counts[as.vector(inTaxidgene),]
    rownames(counts_filted) = as.character(outTaxidgene)
    seurat_ob@assays[[assay2use]]@counts=counts_filted

    ## data gene replace
    expr=seurat_ob@assays[[assay2use]]@data
    gene=rownames(expr)
    expr_filtered = expr[as.vector(inTaxidgene),]
    rownames(expr_filtered) = as.character(outTaxidgene)
    seurat_ob@assays[[assay2use]]@data=expr_filtered

    print("由于更改为只针对高变基因进行scale,因此同源转换后的rds舍弃scale.data")
    ## scale.data （scale是所有基因的，恢复这几行代码即可）
    # scale=seurat_ob@assays[[assay2use]]@scale.data
    # gene=rownames(scale)
    # scale_filtered = scale[as.vector(inTaxidgene),]
    # rownames(scale_filtered) = as.character(outTaxidgene)
    # seurat_ob@assays[[assay2use]]@scale.data=scale_filtered 
    
    ## scale.data （scale是高变基因的，恢复这几行代码即可）
    # scale=seurat_ob@assays[[assay2use]]@scale.data
    # gene=rownames(scale)
    # ## 因此我们的scale只有2000，因此我们这里取交集
    # print("我们的scale只有2000，因此我们这里取交集,请注意后续分析，替换后scale基因过少")
    # intersection <- intersect(inTaxidgene, gene)
    # # 获取交集基因位置
    # common_indices <- match(intersection, in2out[,1])
    # # 提取转化基因名
    # extracted_values <- in2out[,2][common_indices]
    # scale_filtered = scale[as.vector(intersection),]
    # rownames(scale_filtered) = as.character(extracted_values)
    # seurat_ob@assays[[assay2use]]@scale.data=scale_filtered
    
    ## meta.features
    features=seurat_ob@assays[[assay2use]]@meta.features
    features_filtered = features[as.vector(inTaxidgene),]
    features_filtered = as.matrix(features_filtered)
    rownames(features_filtered) = as.character(outTaxidgene)
    features_filtered = as.data.frame(features_filtered)
    seurat_ob@assays[[assay2use]]@meta.features=features_filtered

    return(seurat_ob)
}


blast_homo <- function(output_dir, in_gtf, in_fa, in_features, out_gtf, out_fa, out_features ){
    ######  blast  function ####################
    system(paste0("less -S ", in_gtf, " | awk -F\'\t\' -v OFS=\'\t\' \'$3==\"exon\"\' - |awk -F\'\t|; \' -v OFS=\'\t\' \'{{for (f=1; f <= NF; f+=1) {if ($f ~ /gene_id /) {print $1,$2,$3,$4,$5,$6,$7,$8,$f}}}}\' |sed \'s/gene_id //\' |awk -F \"\t\" -v OFS=\"\t\" \'{print $1,$2,$3,$4,$5,$6,$7,$8,\"transcript_id \"$9\"; gene_id \"$9\"; gene_name \"$9\";\"}\' > ", output_dir, "/in_gene_test.gtf" ) ) 

    system(paste0("less -S ", out_gtf,  " | awk -F\'\t\' -v OFS=\'\t\' \'$3==\"exon\"\' - |awk -F\'\t|; \' -v OFS=\'\t\' \'{{for (f=1; f <= NF; f+=1) {if ($f ~ /gene_id /) {print $1,$2,$3,$4,$5,$6,$7,$8,$f}}}}\' |sed \'s/gene_id //\' |awk -F \"\t\" -v OFS=\"\t\" \'{print $1,$2,$3,$4,$5,$6,$7,$8,\"transcript_id \"$9\"; gene_id \"$9\"; gene_name \"$9\";\"}\' > ", output_dir, "/out_gene_test.gtf " ) ) 

    system(glue::glue("module pure && module load gffread && gffread -w {output_dir}/in_gene.fa -g {in_fa} {output_dir}/in_gene_test.gtf && gffread -w {output_dir}/out_gene.fa -g {out_fa} {output_dir}/out_gene_test.gtf" ) )

    system(glue::glue("module pure && module load blast+/2.7.1 && makeblastdb -in {output_dir}/out_gene.fa -dbtype nucl  && blastn -db {output_dir}/out_gene.fa -query {output_dir}/in_gene.fa  -out {output_dir}/blast_results.xls  -evalue 1e-5 -num_threads 10 -outfmt \"6 qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs\" -max_target_seqs=1 " ) )

    suppressPackageStartupMessages( library("dplyr") )
    infeatures <- read.table(in_features, sep="\t",stringsAsFactors=F, quote=""  )
    hits <- read.table(file.path(output_dir,"blast_results.xls"), sep="\t",stringsAsFactors=F,quote="")
    ## filter and unique 
    hits <- hits %>% dplyr::filter(V3 > 70) %>% dplyr::select(V1,V2) %>% dplyr::rename( in_geneid = V1, out_geneid = V2) %>%  dplyr::distinct(in_geneid, .keep_all=T) %>% dplyr::distinct(out_geneid, .keep_all=T)
    if ( nrow(hits)==0 ) {
        stop("two genomes don't have homologene" )
    }
    ## process feature.csv and seurat genename not match problem
    infeatures[,2] <- gsub("_","-",infeatures[,2])
    ## blast result geneid  replaceing gene names.
    infeatures <- infeatures %>% dplyr::select(V1,V2) %>% dplyr::rename( in_geneid = V1, in_genenames = V2)

    if( !is.null(out_features) ){
        outfeatures <- read.table(out_features, sep="\t",stringsAsFactors=F ,quote="" )
        ## process feature.csv and seurat genename not match problem
        outfeatures[,2] <- gsub("_","-",outfeatures[,2])
        ## blast result geneid  replaceing gene names.
        outfeatures <- outfeatures %>% dplyr::select(V1,V2) %>% dplyr::rename( out_geneid = V1, out_genenames = V2)
        hits <- hits %>% dplyr::left_join(y=infeatures,by="in_geneid") %>% dplyr::left_join(y=outfeatures,by="out_geneid") %>% dplyr::select(in_genenames, out_genenames) 
    }else{
        ## 读取 gtf ## 
        suppressPackageStartupMessages( library("rtracklayer") )
        gtf <- rtracklayer::import(out_gtf)
        gtf_df <- as.data.frame(gtf)
        if ("gene_name" %in% colnames(gtf_df) ) {
            outfeatures <- gtf_df %>% dplyr::filter(type=="exon") %>% dplyr::select(gene_id, gene_name) %>% dplyr::distinct(gene_id, gene_name, .keep_all=T) %>% dplyr::rename(out_geneid = gene_id, out_genenames= gene_name)
            hits <- hits %>% dplyr::left_join(y=infeatures,by="in_geneid") %>% dplyr::left_join(y=outfeatures,by="out_geneid") %>% dplyr::select(in_genenames, out_genenames)
        }else{
            hits <- hits %>% dplyr::left_join(y=infeatures,by="in_geneid")  %>% dplyr::select(in_genenames, out_geneid) 
        }
    }

    #unique
    index <- duplicated(hits[,1])
    hits <-hits[!index,]
    index<-duplicated(hits[,2])
    hits <- hits[!index,]
    write.table(hits,file.path(output_dir,"homologene_in_2_out.xls"),sep="\t",quote = F,row.names=F)

    ### remove files 
    file.remove(file.path(output_dir,list.files(path =output_dir, pattern ="^out")), file.path(output_dir,list.files(path =output_dir, pattern ="^in")) )
    return(hits)
}

#=command line parameters setting=============================
option_list = list(
    make_option( c("--RDS", "-v"), type = "character", default = NULL,
        help = "the seurat object saved as R object in RDS format."),
    make_option( c("--output","-o"),type="character", default = "./",
        help="the output directory of results.", metavar="character"),
    make_option( c("--assay"), type = "character", default = "RNA",
        help = "[OPTIONAL]the default assay to use  for the this run.[default:RNA]"),
    make_option( c("--inTaxid","-i"),type="character", default = NULL,
        help="The taxonomy ID consistenting with input RDS corresponding species" ),
    make_option( c("--outTaxid","-t"),type="character", default = NULL,
        help="The taxonomy ID of the species that needs to be converted into homologous genes" ),
    make_option( c("--genelist","-g"),type="character", default = NULL,
        help="two species's genelist matched, first colum genes is consistenting with input RDS corresponding species, second colum genes are other species' homologous genes"),
    make_option( c("--blast"),type = "logical", default = FALSE,
        help="whether doing blast combined with ingenome and outgenome"), 
    make_option( c("--ingenome"),type="character", default = NULL,
        help="genome fasta, gtf and feature.tsv consistenting with knowns RDS corresponding species, the three are separated by commas." ),
    make_option( c("--outgenome"),type="character", default = NULL,
        help="genome fasta, gtf and feature.tsv of the species that needs to be converted into homologous genes, the three are separated by commas.")
    );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



#=================================================================================
# parse the command line parameters
#=================================================================================
# output directory setting
if ( is.null(opt$output) ){
	print("NO output directory specified,the current directory will be used!")
	output_dir = getwd()
}else{
	if ( file.exists(opt$output)){
		output_dir = opt$output
	}else{
		output_dir = opt$output
		dir.create(output_dir,recursive=T)
	}
}

if ( is.null(opt$assay) ){
    assay2use = "RNA"
}else{
    assay2use = opt$assay
}

if ( opt$blast ) {
	print("begining blasting")
    ingenome = unlist( strsplit(opt$ingenome,",",perl = T) )
    outgenome = unlist( strsplit(opt$outgenome,",",perl = T) )
    ## Check whether the input parameters meet requirements
    if ( length(ingenome) == 3 & length(outgenome) %in% c(2,3)  ) {
        if ( ! basename(ingenome[1]) %in% "genome.fa" | ! basename(outgenome[1]) %in% "genome.fa" ) {
            stop("ingenome and outgenome info' first file must be genome fasta." )
        }
        if ( ! ( basename(ingenome[2]) %in% "genes.gtf" | basename(ingenome[2]) %in% "genes.gtf.gz" ) | ! ( basename(outgenome[2]) %in% "genes.gtf" | basename(outgenome[2]) %in% "genes.gtf.gz" ) ) {
            stop("ingenome and outgenome info' second file must be genes.gtf." )
        }
        in_fa <- ingenome[1]
        out_fa <- outgenome[1]
        in_gtf <- ingenome[2]
        out_gtf <- outgenome[2]
        if (! ( basename(ingenome[3]) %in% "features.tsv" | basename(ingenome[3]) %in% "features.tsv.gz" )  ) { 
            stop("ingenome info' third file must be features.tsv." )
        }else{
            if ( length(unlist(strsplit(basename(ingenome[3]),"[.]", perl=T)) ) == 3 ) {
                system(glue::glue("gunzip -c {ingenome[3]} > {output_dir}/in_features.tsv"))
                in_features <- file.path(output_dir,"in_features.tsv")
            }else{
                in_features <- ingenome[3]
            }
        }
        if( length(outgenome) == 3 )
            if (! ( basename(outgenome[3]) %in% "features.tsv" | basename(outgenome[3]) %in% "features.tsv.gz" )  ) { 
                stop("ingenome info' third file must be features.tsv." )
            }else{
                if ( length(unlist(strsplit(basename(outgenome[3]),"[.]", perl=T)) ) == 3 ) {
                    system(glue::glue("gunzip -c {outgenome[3]} > {output_dir}/out_features.tsv"))
                    out_features <- file.path(output_dir,"out_features.tsv")
                }else{
                    out_features <- outgenome[3]
                }
        }else{
            out_features <- NULL 
        }
    }else{
        stop("ingenome info must be tree files and outgenome info must be have fasta and gtf file " )
    }
}


if ( is.null( opt$RDS ) ) {
	if ( opt$blast) {
		in2out <- blast_homo(output_dir, in_gtf, in_fa, in_features, out_gtf, out_fa, out_features )
    }else{
		stop("rds is NOT AVAILABLE and blast is not true." ) 
	}    
}else{
	suppressPackageStartupMessages(library("Seurat"))
	suppressPackageStartupMessages(library("homologene"))
	suppressPackageStartupMessages(library("OESingleCell"))
	print("Reading rds")
	seurat_ob = readRDSMC(opt$RDS)
	if ( seurat_ob@version < 3){
		seurat_ob = UpdateSeuratObject(seurat_ob) #make sure the seurat object match with the latest seurat package
	}
    if ( opt$blast ) {
        ## find homologene list by blast 
        in2out <- blast_homo(output_dir, in_gtf, in_fa, in_features, out_gtf, out_fa, out_features )
        ##  counts, data gene replace 
        seurat_ob <- homologene_seurat_replace(seurat_ob,assay2use,in2out[,1],in2out[,2])
        ## save rds 
		print("saving rds") 
		saveRDSMC(seurat_ob,file.path(output_dir,paste0("in_2_out.",basename(opt$RDS))))
    }else{
        if ( is.null(opt$inTaxid) & is.null(opt$genelist) ){
            stop("the inTax id is NOT AVAILABLE!")
        }else{
            inTaxid=opt$inTaxid
        }

        if ( is.null(opt$outTaxid) & is.null(opt$genelist) ){
            stop("the outTax id is NOT AVAILABLE!")
        }else{
            outTaxid=opt$outTaxid
        }   
        ## homologene and counts, data gene replace ## 
        if ( !is.null(opt$genelist) ){
            in2out <- read.delim(opt$genelist,sep="\t",header=T)
            inTaxid <- colnames(in2out)[1]
            outTaxid <- colnames(in2out)[2]
            ## unique
            index <- duplicated(in2out[,1])
            in2out <-in2out[!index,]
            index<-duplicated(in2out[,2])
            in2out <- in2out[!index,]
            ## confirm gene in seurat object 
            counts=seurat_ob@assays[[assay2use]]@counts
            filtered_gene = in2out[,1][!in2out[,1]  %in% row.names(counts)]
            if(length(filtered_gene)!=0){
                filtered_gene = as.data.frame(filtered_gene)
                colnames(filtered_gene) = "Gene"
                write.table(filtered_gene,file.path(output_dir,"genes_not_matched_seurat_genes.xls"),quote = F,row.names=F)
                print("There are some mismatched gene symbol, Please check genes_not_matched_seurat_genes.xls for the genename.")
            }
            final_gene = in2out[,1][in2out[,1]  %in% row.names(counts)]
            rownames(in2out)=as.vector(in2out[,1])
            in2out <- in2out[as.vector(final_gene),]
            write.table(in2out,file.path(output_dir,paste0("homologene.",inTaxid,"_2_",outTaxid,".xls",collapse=".")),sep="\t",row.names=F,quote=F)
            ##  counts, data gene replace 
            seurat_ob <- homologene_seurat_replace(seurat_ob,assay2use,in2out[,1],in2out[,2])
        }else{
            ## find homologene list by homologene package  
            counts=seurat_ob@assays[[assay2use]]@counts
            gene=row.names(counts)
            # inTax is input species tax id , outTax is the transformed species tax id
            in2out = homologene(gene,inTax=inTaxid,outTax=outTaxid)
            ## unique
            index <- duplicated(in2out[,1])
            in2out <-in2out[!index,]
            index<-duplicated(in2out[,2])
            in2out <- in2out[!index,]
            homologene <- in2out[,1:2]
            write.table(homologene,file.path(output_dir,paste0("homologene.",inTaxid,"_2_",outTaxid,".xls",collapse=".")),sep="\t",row.names=F,quote=F)
            ##  counts, data gene replace 
            seurat_ob <- homologene_seurat_replace(seurat_ob,assay2use,in2out[,1],in2out[,2])
        }
        ## save RDS ## 
        print("saving rds") 
		saveRDSMC(seurat_ob,file.path(output_dir,paste0(inTaxid,"_2_",outTaxid,".",basename(opt$RDS))))
	}
}
