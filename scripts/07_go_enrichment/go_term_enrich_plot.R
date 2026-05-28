#----------------------------------run the GO analysis----------------------------------#
library(topGO)
setwd('/Users/zhonh0b/Desktop/epigenetic_nanopore/GO_analysis/')

mult_files = list.files("./",pattern = "diffexpr-results-*.")

# bp,cc,mf mean the three components: biological process, molecular function, cellular component
# for each one, there would be a file, so totally 3 files output
setwd('./GO_analysis/')
mult_files=list.files("./",pattern = "*.txt")

for (go_category in c('bp','cc','mf')) {
  annot_filename = '/Users/zhonh0b/Desktop/epigenetic_nanopore/aip_updated_transcript_level_annotation_merged.tsv'
  gene_id_to_go = readMappings(file=annot_filename,sep = '\t')
  gene_id_to_go = gene_id_to_go[gene_id_to_go != 'no_hit']
  gene_names = names(gene_id_to_go)
  
  for (m in mult_files) {
    print(paste("Current file:", m))
    genes_of_interest_filename = paste0("./", m)
    genes_of_interest = scan(genes_of_interest_filename, character(0), sep="\n")
    
    # hack: spis/smic genes have spaces in them, but gene_id_to_go has had
    # spaces stripped from them - so we've got to strip spaces from 
    # genes_of_interest as well.
    genes_of_interest = gsub(" ", "", genes_of_interest) # this one is just get the right format of the gene, if yours are already correct, then no need to do this step
    
    genelist = factor(as.integer(gene_names %in% genes_of_interest))
    names(genelist) = gene_names
    
    GOdata = try(new("topGOdata", ontology=toupper(go_category), allGenes=genelist, gene2GO=gene_id_to_go, annotationFun=annFUN.gene2GO))
    
    # handle error
    if (class(GOdata) == "try-error") {
      print (paste0("Error for file", m, "!"))
      next
    }
    
    # weight01 is the default algorithm used in Alexa et al. (2006)
    weight01.fisher <- runTest(GOdata, statistic = "fisher")
    
    # generate a results table (for only the top 1000 GO terms)
    #   topNodes: highest 1000 GO terms shown
    #   numChar: truncates GO term descriptions at 1000 chars (basically, disables truncation)
    results_table = GenTable(GOdata, P_value=weight01.fisher, orderBy="P_value", topNodes=1000, numChar=1000)
    
    # write it out into a file for python post-processing
    output_filename = paste0("./topGO_output/", go_category, "_", m)
    write.table(results_table, file=output_filename, quote=FALSE, sep='\t')
  }
}


# plot
library(ggplot2)
library(scales)
setwd("./topGO_output")

ntop <- 30 # select the top 30 go term for visualization

# plot bp
for (m in mult_files) {
  print(m)
  m1<-gsub(".txt","",m)
  goEnrichment<-read.table(paste0("bp_",m),sep='\t',header = T,row.names = 1)
  ggdata <- goEnrichment[1:ntop,]
  ggdata$Term <- str_wrap(ggdata$Term, width = 70)
  ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order
  ggplot(ggdata,
         aes(x = Term, y = -log10(P_value), size = -log10(P_value), fill = -log10(P_value))) +
    
    expand_limits(y = 1) +
    geom_point(shape = 21) +
    scale_size(range = c(2.5,12.5)) +
    scale_fill_continuous(low = 'royalblue', high = 'red4') +
    
    xlab('') + ylab('Enrichment score') +
    labs(
      title = 'GO Biological processes',
      subtitle = 'Top 30 terms ordered by Kolmogorov-Smirnov p-value',
      caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001') +
    
    geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
               linetype = c("dotted", "longdash", "solid"),
               colour = c("black", "black", "black"),
               size = c(0.5, 1.5, 3)) +
    
    theme_bw(base_size = 24) +
    theme(
      legend.position = 'right',
      legend.background = element_rect(),
      plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
      plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
      plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
      
      axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10,colour = 'black'),
      axis.text.y = element_text(angle = 0, size = 11, face = 'bold', colour = 'black',lineheight = 0.85,),
      axis.title = element_text(size = 12, face = 'bold'),
      axis.title.x = element_text(size = 12, face = 'bold'),
      axis.title.y = element_text(size = 12, face = 'bold'),
      axis.line = element_line(colour = 'black'),
      
      #Legend
      legend.key = element_blank(), # removes the border
      legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
      legend.text = element_text(size = 14, face = "bold"), # Text size
      title = element_text(size = 14, face = "bold")) +
    theme(plot.margin = margin(10, 10, 10, 20))+
    coord_flip()
  ggsave(paste0("../BP_plot_top30/bp_",m1,".png"),  
         dpi = 600,
         width = 11,
         height = 9,
         bg = "white")
}

# plot mf
for (m in mult_files) {
  print(m)
  m1<-gsub(".txt","",m)
  goEnrichment<-read.table(paste0("mf_",m),sep='\t',header = T,row.names = 1)
  ggdata <- goEnrichment[1:ntop,]
  ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order
  ggplot(ggdata,
         aes(x = Term, y = -log10(P_value), size = -log10(P_value), fill = -log10(P_value))) +
    
    expand_limits(y = 1) +
    geom_point(shape = 21) +
    scale_size(range = c(2.5,12.5)) +
    scale_fill_continuous(low = 'royalblue', high = 'red4') +
    
    xlab('') + ylab('Enrichment score') +
    labs(
      title = 'GO Molecular Funcation',
      subtitle = 'Top 30 terms ordered by Kolmogorov-Smirnov p-value',
      caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001') +
    
    geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
               linetype = c("dotted", "longdash", "solid"),
               colour = c("black", "black", "black"),
               size = c(0.5, 1.5, 3)) +
    
    theme_bw(base_size = 24) +
    theme(
      legend.position = 'right',
      legend.background = element_rect(),
      plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
      plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
      plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
      
      axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
      axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
      axis.title = element_text(size = 12, face = 'bold'),
      axis.title.x = element_text(size = 12, face = 'bold'),
      axis.title.y = element_text(size = 12, face = 'bold'),
      axis.line = element_line(colour = 'black'),
      
      #Legend
      legend.key = element_blank(), # removes the border
      legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
      legend.text = element_text(size = 14, face = "bold"), # Text size
      title = element_text(size = 14, face = "bold")) +
    
    coord_flip()
  ggsave(paste0("../MF_plot_top30/mf_",m1,".png"))
}

# plot cc
for (m in mult_files) {
  print(m)
  m1<-gsub(".txt","",m)
  goEnrichment<-read.table(paste0("cc_",m),sep='\t',header = T,row.names = 1)
  ggdata <- goEnrichment[1:ntop,]
  ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order
  ggplot(ggdata,
         aes(x = Term, y = -log10(P_value), size = -log10(P_value), fill = -log10(P_value))) +
    
    expand_limits(y = 1) +
    geom_point(shape = 21) +
    scale_size(range = c(2.5,12.5)) +
    scale_fill_continuous(low = 'royalblue', high = 'red4') +
    
    xlab('') + ylab('Enrichment score') +
    labs(
      title = 'GO Cellular Compoment',
      subtitle = 'Top 30 terms ordered by Kolmogorov-Smirnov p-value',
      caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001') +
    
    geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
               linetype = c("dotted", "longdash", "solid"),
               colour = c("black", "black", "black"),
               size = c(0.5, 1.5, 3)) +
    
    theme_bw(base_size = 24) +
    theme(
      legend.position = 'right',
      legend.background = element_rect(),
      plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
      plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
      plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
      
      axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
      axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
      axis.title = element_text(size = 12, face = 'bold'),
      axis.title.x = element_text(size = 12, face = 'bold'),
      axis.title.y = element_text(size = 12, face = 'bold'),
      axis.line = element_line(colour = 'black'),
      
      #Legend
      legend.key = element_blank(), # removes the border
      legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
      legend.text = element_text(size = 14, face = "bold"), # Text size
      title = element_text(size = 14, face = "bold")) +
    
    coord_flip()
  ggsave(paste0("../CC_plot_top30/cc_",m1,".png"))
}


