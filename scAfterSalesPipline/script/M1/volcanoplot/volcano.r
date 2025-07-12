#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
# suppressPackageStartupMessages(library(oeRtools))
library(ggplot2)
library(ggrepel)
library(stringr)

option_list <- list(
    make_option(c("-i", "--diff"),
        type = "character",
        help = "diff-all.gene.xls", metavar = "character"
    ),
    make_option(c("-P", "--prefix"),
        type = "character", default = "auto",
        help = "", metavar = "character"
    ),
    make_option(c("-p", "--pvalue"),
        type = "double",
        help = paste0(
            "pvalue ratio threshold.Filtering can be performed ",
            "using any one of (-p), (-f) at a time"
        ), metavar = "double"
    ),
    make_option(c("-q", "--fdr"),
        type = "double",
        help = paste0(
            "fdr ratio threshold.Filtering can be performed ",
            "using any one of (-p), (-f) at a time"
        ), metavar = "double"
    ),
    make_option(c("-f", "--foldchange"),
        type = "double", default = 1.5,
        help = "foldchange threshold [default %default]", metavar = "double"
    ),
    make_option(c("--symbol_gene"),
        type = "character", default = NULL,
        help = "label gene list",
        metavar = "character"
    ),
    make_option(c("-t", "--type"),
        type = "character", default = "gene",
        help = "transcript or gene type:mRNA,circRNA,lncRNA,miRNA,gene",
        metavar = "character"
    ),
    make_option(c("-o", "--outputdir"),
        type = "character", default = "./",
        help = "output directory for results", metavar = "character"
    ),
    make_option(c("--xlim"),
        type = "character", default = "auto",
        help = "set x-axis limit, e.g., '10' or '-10,7'", metavar = "character"
    ),
    make_option(c("--xcoord"),
        type = "character", default = "auto",
        help = "Adjusted the position of Conte Info, e.g., '-10,7'", metavar = "character"
    ),
    make_option(c("--ycoord"),
        type = "character", default = "auto",
        help = "Adjusted the position of Conte Info, e.g., '5,7'", metavar = "character"
    ),
    make_option(c("--colors"),
        type = "character", default = "auto",
        help = " 'Significant Up', 'Significant Down', 'Up', 'Down', 'Non-significant', auto=c('#ff4d40', '#4682b4', '#febfca', '#86cdf9', '#d2dae2')",
        metavar = "character"
    ),
    make_option(c("--geneTextColors"),
        type = "character", default = "auto",
        help = " 'Significant Up', 'Significant Down', 'Up', 'Down', 'Non-significant', auto: Gene text color follows the color of the dots",
        metavar = "character"
    ),
    make_option(c("--symbol_topn"),
        type = "double", default = NULL,
        help = "Displays the topn gene name", metavar = "double"
    ),
    make_option(c("--symbol_fc"),
        type = "double", default = 2,
        help = "Displays the topn gene name", metavar = "double"
    )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

identifySymbolGeneType <- function(x) {
  is_file_path <- function(x) {
    grepl("\\.xls$|\\.xlsx$|\\.csv$|\\.txt$", x)
  }
  
  is_gene_list <- function(x) {
    !is_file_path(x) && grepl(",", x)
  }
  
  is_single_gene <- function(x) {
    !is_file_path(x) && !grepl(",", x)
  }
  
  if (is_file_path(x)) {
    return("filepath")
  } else if (is_gene_list(x)) {
    return("genelist")
  } else if (is_single_gene(x)) {
    return("singlegene")
  } else {
    return("stop")
  }
}

getXlimOptimized <- function(xlim, log2FC) {
  if (xlim == "auto") {
    return(c(-1, 1) * round(max(abs(log2FC)) * 1.1 + 0.5))
  }
  xlimValues <- as.numeric(unlist(str_split(xlim, ",")))
  if (length(xlimValues) == 1) {
    return(c(-xlimValues, xlimValues))
  } else if (length(xlimValues) == 2) {
    return(xlimValues)
  } else {
    stop("Invalid xlim specification. Please provide a single number or a comma-separated pair, e.g., '10' or '-10,10'.")
  }
}

convertToNumeric <- function(x) {
  result <- tryCatch({
    numeric_result <- as.numeric(x) 
    if(any(is.na(numeric_result))) { 
      stop("The input cannot be converted to numeric.") 
    }
    return(numeric_result) 
  }, error=function(e){ 
    stop("Error: ", e$message) 
  })
  return(result)
}

volcano_plot <- function(prefix, FC, Pvalue, FDR) {
    data <- read.delim(
        opt$diff,
        header = T, row.names = 1, sep = "\t", quote = "",
        comment.char = ""
    )

    if (!is.null(FDR)) {
        data <- data[which(data["q.value"] != ""), ]
        P <- data["q.value"]
        picname <- c(
            paste0(
                opt$outputdir, "/", prefix, "-volcano", "-q-val-", FDR,
                "-FC-", FC, ".", opt$type
            )
        )
        type <- c("q-value")
        thresh <- FDR
    }else if(!is.null(Pvalue)) {
        data <- data[which(data["p.value"] != ""), ]
        P <- data["p.value"]
        picname <- c(
            paste0(
                opt$outputdir, "/", prefix, "-volcano", "-p-val-",
                Pvalue, "-FC-", FC, ".", opt$type
            )
        )
        type <- c("p-value")
        thresh <- Pvalue
    }else{
        stop("Please provide either pvalue or fdr")
    }

    log2FC <- data$log2FoldChange
    df <- data.frame(P, log2FC)
    colnames(df)[1] <- "P"
    df$P[which(df$P < 5E-300 )] = 5E-300
    tmp <- df[which(df$log2FC != "Inf" & df$log2FC != "-Inf"), ]
    max <- max(tmp$log2FC)
    min <- min(tmp$log2FC)
    df$log2FC <- as.numeric(sub("-Inf", min, df$log2FC))
    df$log2FC <- as.numeric(sub("Inf", max, df$log2FC))

    df$group <- "Non-significant"
    df$group[df$log2FC < 0 & df$P < thresh] <- "Down"
    df$group[df$log2FC > 0 & df$P < thresh] <- "Up"
    df$group[df$log2FC < -log2(FC) & df$P < thresh] <- "Significant Down"
    df$group[df$log2FC > log2(FC) & df$P < thresh] <- "Significant Up"
    df$group <- factor(
        df$group,
        levels = c(
            "Significant Up",
            "Significant Down",
            "Up",
            "Down",
            "Non-significant"
        )
    )

##########label
    if(!is.null(opt$symbol_gene)){
       symbol_gene_type = identifySymbolGeneType(opt$symbol_gene)
       if(symbol_gene_type == "filepath"){
       show <- read.delim(opt$symbol_gene,header = T,  sep = "\t", quote = "",comment.char = "")
            df_label<-df[row.names(df) %in% show[,1],]
       }else if(symbol_gene_type == "genelist"){
            df_label<-df[row.names(df) %in% unlist(str_split(opt$symbol_gene, "\\s*,\\s*")),]
       }else if (symbol_gene_type == "singlegene"){
            df_label<-df[row.names(df) %in% opt$symbol_gene,]
       }else{
            stop("please check your --symbol_gene!")
       }
       print("genes will be labeled: ")
       print(rownames(df_label))
       df_label$label <-as.character(rownames(df_label))
    }else if(!is.null(opt$symbol_topn)){
    df_up <- subset(
        df, P < ifelse(!is.null(Pvalue), Pvalue, FDR) & df$log2FC > log2(FC)
    )
    df_up <- head(df_up[order(df_up$log2FC, decreasing = T), ], as.numeric(opt$symbol_topn))

    df_down <- subset(
        df, P < ifelse(!is.null(Pvalue), Pvalue, FDR) & df$log2FC < -log2(FC)
    )
    df_down <- head(df_down[order(df_down$log2FC, decreasing = F), ], as.numeric(opt$symbol_topn))
       df_label <- rbind(df_up, df_down)
       df_label$label <- rownames(df_label)
    }
    else{
       show <- (df$P< opt$pvalue & abs(df$log2FC) > log(opt$symbol_fc, 2) )  # 基因名展示阈值

       df_label<-df[show,]
       df_label$label <-as.character(rownames(df_label))
    }
#################################

    if (opt$colors == "auto") {
        colors <- c("#ff4d40", "#4682b4", "#febfca", "#86cdf9", "#d2dae2")
    } else {
        colors <- unlist(str_split(opt$colors, "\\s*,\\s*"))
        if (length(colors) < 5) {
            stop("Invalid colors specification. Please provide at least 5 colors separated by commas.")
        }
    }

    p <- ggplot(
        df, aes(x = log2FC, y = -log10(P), colour = group)
    ) +
        xlim(getXlimOptimized(opt$xlim, df$log2FC)) +
        geom_point(alpha = 0.8, size = 1) +
        scale_color_manual(
            values = colors[levels(df$group) %in% unique(df$group)]
        ) +
        geom_vline(
            xintercept = c(-log2(FC), log2(FC)), lty = 2, col = "grey",
            lwd = 0.4
        ) +
        geom_hline(
            yintercept = -log10(ifelse(!is.null(Pvalue), Pvalue, FDR)),
            lty = 2, col = "grey", lwd = 0.4
        ) +
        labs(
            x = expression(paste(log[2], " FoldChange")),
            y = ifelse(
                !is.null(Pvalue),
                expression(paste("-", log[10], "p-value")),
                expression(paste("-", log[10], "q-value"))
            )
        ) +
        labs(
            title = substitute(
                A*log[2]~B,
                list(
                    A = paste0(
                        prefix, ": ",
                        type, "<", thresh, " && |"
                    ),
                    B = paste0("FC|>", round(log2(FC), 2))
                )
            )
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 12),
            legend.position = "right",
            legend.title = element_blank(),
            panel.grid = element_blank(),
            legend.key.height = unit(0.5, "cm"),
            legend.text = element_text(size = 10),
            panel.border = element_rect(
                color = "black", size = 1, fill = NA
            ),
            axis.text = element_text(size = 9),
            axis.title = element_text(size = 13)
        )

    if (opt$xcoord == "auto") {
        xcoord = getXlimOptimized(opt$xlim, df$log2FC) * c(1, 0.7)
    } else {
        xcoord <- convertToNumeric(unlist(str_split(opt$xcoord, ",")))
        if (length(xcoord) != 2) {
            stop("Invalid xcoord specification. Please provide a comma-separated pair, e.g., '-10,7'.")
        }
    }

    if (opt$ycoord == "auto") {
        ycoord <- min(-log10(df$P)) +
            diff(
                range(
                    -log10(df$P)[!is.infinite(-log10(df$P))]
                )
            ) * .97
        ycoord = c(ycoord,ycoord)
    } else {
        ycoord <- convertToNumeric(unlist(str_split(opt$ycoord, ",")))
        if (length(ycoord) != 2) {
            stop("Invalid ycoord specification. Please provide a comma-separated pair, e.g., '6,7'.")
        }
    }

    count_info <- data.frame(
        x = xcoord,
        y = ycoord,
        label = c(
            paste0(
                "Down\n",
                ifelse(!is.null(Pvalue), "p", "q"),
                " < ", ifelse(!is.null(Pvalue), Pvalue, FDR), ": ",
                nrow(
                    subset(
                        df,
                        P < ifelse(
                            !is.null(Pvalue), Pvalue, FDR
                        ) & df$log2FC < 0
                    )
                ),
                "\nlog2FC < ", round(-log2(FC),2), ": ",
                nrow(subset(df, df$log2FC < -log2(FC))),
                "\nSig. : ",
                                nrow(
                    subset(
                        df,
                        P < ifelse(
                            !is.null(Pvalue), Pvalue, FDR
                        ) & df$log2FC < -log2(FC)
                    )
                )
            ),
            paste0(
                "Up\n",
                ifelse(!is.null(Pvalue), "p", "q"),
                " < ", ifelse(!is.null(Pvalue), Pvalue, FDR), ": ",
                nrow(
                    subset(
                        df,
                        P < ifelse(
                            !is.null(Pvalue), Pvalue, FDR
                        ) & df$log2FC > 0
                    )
                ),
                "\nlog2FC > ", round(log2(FC),2), ": ",
                nrow(subset(df, df$log2FC > log2(FC))),
                "\nSig. : ",
                nrow(
                    subset(
                        df,
                        P < ifelse(
                            !is.null(Pvalue), Pvalue, FDR
                        ) & df$log2FC > log2(FC)
                    )
                )
            )
        )
    )

    p <- p + geom_text(
        parse = F,
        data = count_info,
        color = "black",
        size = 3,
        hjust = 0,
        aes(x = x, y = y, label = label)
    )

    if (opt$geneTextColor != "auto"){
        geneTextColors = str_split(opt$geneTextColors, "\\s*,\\s*")[[1]]
        geneTextColorsIdx = data.frame(sig = levels(df$group), gtc = geneTextColors)
        df_label <- merge(df_label, geneTextColorsIdx, by.x = "group", by.y = "sig", all.x = T)
        if (nrow(df_label) > 0) {
            p <- p + geom_text_repel(
                data = df_label,
                aes(label = label),
                size = 3,
                color = df_label$gtc,
                box.padding = .6,
                max.overlaps = Inf,
                min.segment.length = 0.25,
                segment.alpha = .4,
                show.legend = FALSE
            )
        }
    }else{
        if (nrow(df_label) > 0) {
            p <- p + geom_text_repel(
                data = df_label,
                aes(label = label),
                size = 3,
                box.padding = .6,
                max.overlaps = Inf,
                min.segment.length = 0.25,
                segment.alpha = .4,
                show.legend = FALSE
            )
        }
    }


    picname <- paste0(
        prefix,"-volcano", ifelse(!is.null(Pvalue), "-p-val-", "-q-val-"),
        ifelse(!is.null(Pvalue), Pvalue, FDR), "-FC-", FC, ".", opt$type
    )

    dir.create(file.path(opt$outputdir), showWarnings = FALSE)
    ggsave(file.path(opt$outputdir, paste0(picname,".pdf")),p, width = 10, height = 8)
    ggsave(file.path(opt$outputdir,paste0(picname,".png")), p, dpi=300,height = 8, width = 10,bg = "white")

    colnames(df)[1] <- type
    colnames(df)[3] <- "regulation"
    df$genes <- rownames(df)
    write.csv(df, file.path(opt$outputdir, paste0(picname,".csv")), row.names = F, quote = F)
}

generate_prefix <- function(filepath, type = "auto", delimiter = "-vs-", split_str = "-all") {
  if (is.null(filepath) || filepath == "") {
    stop("File path cannot be empty")
  }
  
  filename <- basename(filepath)
  
  if (is.null(type) || type == "") {
    stop("Type cannot be empty")
  }
  
  if (type == "auto") {
    if (grepl(delimiter, filename)) {
      prefix <- strsplit(filename, split = split_str)[[1]][1]
    } else {
      prefix <- "case-vs-control"
    }
  } else {
    prefix <- type
  }
  
  return(prefix)
}

prefix <- generate_prefix(opt$diff, opt$prefix)

volcano_plot(prefix, opt$foldchange, opt$pvalue, opt$fdr)
