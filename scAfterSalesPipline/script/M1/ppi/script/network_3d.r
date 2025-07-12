library(networkD3)
suppressPackageStartupMessages(library("optparse"))
option_list = list(
        make_option(c("-l", "--links"), type="character", default=NULL, help="top 300 ppi network results file name", metavar="character"),
	make_option(c("-n", "--nodes"), type="character", default=NULL, help="top 300 ppi network nodes", metavar="character"),
        make_option(c("-o", "--picname"), type="character", default=NULL, help="output picture name", metavar="character")
       );
opt_parser = OptionParser(option_list=option_list,epilogue = "Rscript networkD3.r -l top_300_diff-Isc-vs-A_Con_protein2protein_network.xls_res.xls -n top_300_diff-Isc-vs-A_Con_protein2protein_network.xls_tmp -o top_300_diff-Isc-vs-A_Con_protein2protein_network.xls_tmp.html");
opt = parse_args(opt_parser);
if (is.null(opt$links) | is.null(opt$nodes) | is.null(opt$picname)){
	print_help(opt_parser)
	stop("--links --nodes --picname must be supplied", call.=FALSE)
}

MisLinks<-read.table(opt$links,header=T,sep="\t")
MisNodes<-read.table(opt$nodes,header=T,sep="\t")
ColourScale <- 'd3.scaleOrdinal()
            .domain(["Up", "Down","Up_TFs","Down_TFs"])
           .range(["#FF0000","#00FF00","#A020F0","#0000FF"]);'
html<-forceNetwork(Links = MisLinks, Nodes = MisNodes, Source = "prot1", Target = "prot2", Value = "score", NodeID = "name", Group = "group", opacity = 0.8, linkColour=MisLinks$col,arrows=TRUE,linkDistance = networkD3::JS("function(d) { return d.value*120; }"),Nodesize = "size" ,zoom = TRUE,opacityNoHover=TRUE,legend=TRUE,colourScale = JS(ColourScale))
saveNetwork(html,opt$picname,selfcontained=TRUE)
