library(IsoformSwitchAnalyzeR)
rm(list=ls())
setwd('/Users/zhonh0b/Desktop/epigenetic_nanopore/isoformSwitchAnalyzeR/only_3_replicates_default/')

### Import quantifications
count=read.csv('../input_data/OUT.transcript_model_grouped_counts.tsv',header = T,sep = '\t')

tpm_data=read.csv('../input_data/OUT.transcript_model_grouped_tpm.tsv',header = T,sep = '\t')

colnames(tpm_data)[1]='isoform_id'
colnames(count)[1]='isoform_id'

count=count[,c(1,2,3,4,6,7,8)]
tpm_data=tpm_data[,c(1,2,3,4,6,7,8)]

### Make design matrix
myDesign <- data.frame(
  sampleID = c('A1','A2','A3','H1','H2','H3'),
  condition = c(rep('Aposymbiotic',3),rep('Symbiotic',3)),
  batch=factor(c(1,2,1,2,3,2)) # rna_batch + date_batch
)

### Create switchAnalyzeRlist
aSwitchList <- importRdata(
  isoformCountMatrix   = count,
  isoformRepExpression = tpm_data,
  designMatrix         = myDesign,
  isoformExonAnnoation = "../input_data/OUT.transcript_models.gtf",
  isoformNtFasta       = "../input_data/aip_expressed_isoform.fa",
  #onlyConsiderFullORF = T, # this is to make sure we only add the orf that is annotated with both a start and stop codon
)


summary(aSwitchList)

## step1: filter the non-expressed gene/isoforms, identify isoform switches, annotates open reading frames, switches with and extracts both the 
## nucleotide and peptide sequences and output them as two seperate fasta files
### Create SwitchAnalyzeRlist

### Filter
mySwitchList <- preFilter(aSwitchList,
                          removeSingleIsoformGenes = T)   

### Test for isoform switches
mySwitchList <- isoformSwitchTestDEXSeq( mySwitchList,
                                         reduceToSwitchingGenes = T,
                                         reduceFurtherToGenesWithConsequencePotential = F,
                                         alpha = 0.05,
                                         dIFcutoff = 0.1,
                                         onlySigIsoforms = F)   

extractSwitchSummary(mySwitchList)
write.csv(mySwitchList$isoformFeatures, "./isoformfeatures_combined_batch.tr.csv")

### If analysing (some) novel isoforms (else use CDS from ORF as explained in importRdata() )
mySwitchList=analyzeORF(mySwitchList,orfMethod = 'longest')
#mySwitchList <- addORFfromGTF( mySwitchList ,pathToGTF = '../input_data/aiptasia_genome.dups_removed.gtf')
#mySwitchList <- analyzeNovelIsoformORF( mySwitchList,analysisAllIsoformsWithoutORF = TRUE)

### Extract Sequences
mySwitchList <- extractSequence( mySwitchList, pathToOutput = './')

### Summary
extractSwitchSummary(mySwitchList)


### step2: plot all the isoform switches and their annotation
### need to perform the external analysis first, such as cpc2, pfam, signalip, iupread2a, analyzedeeptmhmm, analyzedeeploc2
### Add annotation
SwitchListAnalyzed <- analyzeCPC2(switchAnalyzeRlist = mySwitchList,
                                  pathToCPC2resultFile = './external_output/cpc2_output.txt',
                                  removeNoncodinORFs = TRUE)

SwitchListAnalyzed <- analyzePFAM(switchAnalyzeRlist = SwitchListAnalyzed,
                                  pathToPFAMresultFile = './external_output/pfam_output.txt')

SwitchListAnalyzed <- analyzeSignalP(switchAnalyzeRlist = SwitchListAnalyzed,
                                     pathToSignalPresultFile = './external_output/singnalp_output.txt')

SwitchListAnalyzed <- analyzeIUPred2A(switchAnalyzeRlist = SwitchListAnalyzed,
                                      pathToIUPred2AresultFile = './external_output/isoformSwitchAnalyzeR_isoform_aa_IUPred2A.result') # OR

SwitchListAnalyzed <- analyzeDeepLoc2(switchAnalyzeRlist = SwitchListAnalyzed,
                                      pathToDeepLoc2resultFile = './external_output/analyzeDeepLoc2_result.csv')

SwitchListAnalyzed <- analyzeDeepTMHMM(switchAnalyzeRlist = SwitchListAnalyzed,
                                       pathToDeepTMHMMresultFile = './external_output/deeptmhmm_output.gff3')

SwitchListAnalyzed

## predict the alternative splicing
SwitchListAnalyzed <- analyzeAlternativeSplicing( SwitchListAnalyzed )

table(SwitchListAnalyzed$AlternativeSplicingAnalysis$IR) # ES,MEE,MES,A3,A5,atss,ATTS,

### plot blow the bar plot for all the alternative splicing detected events
# Load necessary library
library(ggplot2)
library(reshape2)

# Create the data frame
data <- data.frame(
  Event = c("IR", "ES", "MEE", "MES", "A3", "A5", "ATSS", "ATTS"),
  `0` = c(288,327,382,332,320,312,196,175),
  `1` = c(93,58,10,44,69,77,197,218),
  `2` = c(9,8,1,15,4,4,0,0),
  `3` = c(3, 0, 0, 2,0,0,0,0)
)

# Melt the data frame for ggplot2
data_melted <- melt(data, id.vars = "Event", variable.name = "Count", value.name = "Frequency")

# Create the bar plot
p <- ggplot(data_melted, aes(x = Event, y = Frequency, fill = Count)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  #labs(title = "Alternative Splicing Events", x = "Splicing Event Type", y = "Frequency") +
  scale_fill_brewer(palette = "Set1",labels = c("0", "1", "2", "3"))+
  ylab(label = '')+
  xlab(label = '')
p+
  theme(
    # Remove background grids
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Add axis lines
    axis.line = element_line(color = "black"),
    # Remove panel border
    panel.border = element_blank(),
    # Change legend labels
    legend.text = element_text(size = 14),
    axis.text= element_text(size = 16),
    legend.title = element_text(size = 16)
    
  )
# Save the plot
ggsave("./plots/alternative_splicing_events.png")


data <- data.frame(
  Event = c("IR", "ES", "MEE", "MES", "A3", "A5", "ATSS", "ATTS"),
  `0` = c(288,327,382,332,320,312,196,175),
  `1` = c(93,58,10,44,69,77,197,218),
  `2` = c(9,8,1,15,4,4,0,0),
  `3` = c(3,0,0,2,0,0,0,0),
  check.names = FALSE
)

# Melt
data_melted <- melt(data, id.vars = "Event", variable.name = "Count", value.name = "Frequency")

plot_df <- data_melted %>%
  filter(Count != "0", Frequency > 0) %>%
  mutate(
    Event = factor(Event, levels = c("IR", "ES", "MEE", "MES", "A3", "A5", "ATSS", "ATTS")),
    Count = factor(Count, levels = c("1", "2", "3"))
  )

# Total positive counts per event
totals <- plot_df %>%
  group_by(Event) %>%
  summarise(Total = sum(Frequency), .groups = "drop")

# Plot
p <- ggplot(plot_df, aes(x = Event, y = Frequency, fill = Count)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  #geom_text(aes(label = Frequency),
  #          position = position_stack(vjust = 0.5),
  #          size = 4,
  #          color = "white",
  #          fontface = "bold") +
  geom_text(data = totals,
            aes(x = Event, y = Total, label = Total),
            inherit.aes = FALSE,
            vjust = -0.5,
            size = 4.5,
            ) +
  scale_fill_manual(
    values = c("1" = "#8ECFC9", "2" = "#FFBE7A", "3" = "#FA7F6F"),
    name = "Number of\nAS events"
  ) +
  labs(
    x = NULL,
    y = "Number of genes"
  ) +
  expand_limits(y = max(totals$Total) * 1.12) +
  theme_classic(base_size = 24) +
  theme(
    axis.text.x = element_text(size = 14,  color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 16, ),
    legend.title = element_text(size = 14, ),
    legend.text = element_text(size = 12),
    legend.position = "right"
  )

p

ggsave(
  "./plots/alternative_splicing_events_publication.png",
  plot = p,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

consequencesOfInterest <- c(
  'intron_retention',
  'coding_potential',
  'ORF_seq_similarity',
  'NMD_status', #Nonsense-Mediated Decay status,is a cellular process that degrades mRNAs containing premature stop codons, thereby preventing the production of truncated proteins.
  'domains_identified',
  'IDR_identified', #IDR stands for Intrinsically Disordered Region. IDR_identified likely refers to the identification or characterization of regions within proteins that lack a stable tertiary structure under physiological conditions.
  'IDR_type', # IDR_type could refer to the classification or categorization of intrinsically disordered regions based on their characteristics or functional roles.
  'signal_peptide_identified' #A signal peptide is a short peptide sequence found at the N-terminus of newly synthesized proteins in the cell. It guides the protein to its correct location within the cell or outside the cell (secretion).
)

SwitchListAnalyzed <- analyzeSwitchConsequences(
  SwitchListAnalyzed,
  consequencesToAnalyze = consequencesOfInterest,
  showProgress=FALSE
)

extractSwitchSummary(SwitchListAnalyzed, filterForConsequences = FALSE)


## modify the gene id and gene name
SwitchListAnalyzed$isoformFeatures$gene_name=gsub('_gene$','',SwitchListAnalyzed$isoformFeatures$gene_id)
SwitchListAnalyzed$isoformFeatures$gene_ID=SwitchListAnalyzed$isoformFeatures$gene_id
SwitchListAnalyzed$isoformFeatures$gene_id=SwitchListAnalyzed$isoformFeatures$gene_name

gene_symbol=read.csv('../CC7_genesymbol.csv')
head(gene_symbol)
gene_symbol$gene_symbol=gsub(' ','',gene_symbol$gene_symbol)
for (i in 1:393){
  if (SwitchListAnalyzed$isoformFeatures$gene_id[i] %in% gene_symbol$Query) {
    print('haha')
    gene_name_new=gene_symbol[gene_symbol$Query==SwitchListAnalyzed$isoformFeatures$gene_id[i],]$gene_symbol
    SwitchListAnalyzed$isoformFeatures$gene_name[i]=gene_name_new
  }
}

all_top_switches=extractTopSwitches(
  SwitchListAnalyzed, 
  extractGenes = TRUE,
  filterForConsequences = FALSE, 
  n = Inf, 
  sortByQvals = TRUE
)


switchPlotTopSwitches(
  switchAnalyzeRlist = SwitchListAnalyzed, 
  n = Inf,                                          
  filterForConsequences = FALSE,
  fileType = "pdf",                                 
  pathToOutput = "./dtu-plots"
)

head(all_top_switches)
write.csv(all_top_switches,'./all_top_sig_switch_genes_based_q_values.csv')


isoquant_deg=read.csv('../../DEG_analysis/diffexpr-results-isoquant-genes-6samples_update_from_transcript_model_gene_level.csv',header = T,row.names = 1)
head(isoquant_deg)
isoquant_deg=isoquant_deg[isoquant_deg$padj<0.05,]
dim(isoquant_deg)

isoquant_det=read.csv('../../DEG_analysis/diffexpr-results-isoquant-genes-6samples_update_from_transcript_model.csv',header = T,row.names = 1)
head(isoquant_det)
isoquant_det=isoquant_det[isoquant_det$padj<0.05,]
dim(isoquant_det)

all_top_switches_consequence=all_top_switches[all_top_switches$switchConsequencesGene,]

isoquant_deg$Gene=gsub('_gene$','',isoquant_deg$Gene)
isoquant_det$Gene=gsub('_gene$','',isoquant_det$Gene)

# draw the upset plot between the deg, det, deu
library(UpSetR)
gene_list <- list(
  DEG = isoquant_deg$Gene,
  DET = isoquant_det$Gene,
  DTU = all_top_switches$gene_id
)

png("upset_deg_det_dtu_genes.png", width = 24, height = 14, units = "in", res = 600)
upset(fromList(gene_list), order.by = "freq",sets.bar.color = "black",set_size.numbers_size=7,
      text.scale=6,att.color='black',point.size=10,mainbar.y.label='Gene Intersection',)
dev.off()


isoquant_deg[isoquant_deg$Gene %in% all_top_switches_consequence$gene_id,]
all_top_switches_consequence[all_top_switches_consequence$gene_id  %in% isoquant_deg$Gene,]

write.csv(isoquant_deg[isoquant_deg$Gene %in% intersect(isoquant_deg$Gene,all_top_switches$gene_id),]
          ,'isoquant_degs_all_top_sig_switches.csv')
write.table(isoquant_deg[isoquant_deg$Gene %in% intersect(isoquant_deg$Gene,all_top_switches_consequence$gene_id),][,1]
            ,'isoquant_degs_all_top_sig_switches_genes_names.txt',row.names = F,quote = F,col.names = F)

write.csv(isoquant_det[isoquant_det$Gene %in% intersect(isoquant_det$Gene,all_top_switches$gene_id),]
          ,'isoquant_dets_all_top_sig_switches.csv')
write.table(isoquant_det[isoquant_det$Gene %in% intersect(isoquant_det$Gene,all_top_switches_consequence$gene_id),][,1]
            ,'isoquant_dets_all_top_sig_switches_genes_names.txt',row.names = F,quote = F,col.names = F)

write.csv(isoquant_det[isoquant_det$Gene %in% Reduce(intersect, list(isoquant_deg$Gene,isoquant_det$Gene,all_top_switches$gene_id)),]
          ,'isoquant_degs_dets_all_top_sig_switches.csv')
write.table(isoquant_det[isoquant_det$Gene %in% Reduce(intersect, list(isoquant_deg$Gene,isoquant_det$Gene,all_top_switches_consequence$gene_id)),][,1]
            ,'isoquant_degs_dets_all_top_sig_switches_genes_names.txt',row.names = F,quote = F,col.names = F)


## splicing
# Splicing
extractSplicingSummary( SwitchListAnalyzed, 
                        asFractionTotal = FALSE, 
                        plotGenes =FALSE,
                        alpha = 0.05,
                        dIFcutoff = 0.1,returnResult = T)

# only consider significant isoforms, meaning only analyzing genes where at least two isoforms which both have significant usage changes in opposite direction (quite strict)
extractSplicingSummary( SwitchListAnalyzed, 
                        asFractionTotal = FALSE, 
                        plotGenes =FALSE,
                        alpha = 0.05,
                        dIFcutoff = 0.1,returnResult = T,onlySigIsoforms = T)

splicingEnrichment <- extractSplicingEnrichment(
  SwitchListAnalyzed,
  splicingToAnalyze='all',
  returnResult=TRUE,
  returnSummary=TRUE
)

pdf("./plots/splicing-genome-wide.pdf", 15,10)
extractSplicingGenomeWide( SwitchListAnalyzed )
dev.off()

splicingEnrichment <- extractSplicingEnrichment(
  SwitchListAnalyzed,
  splicingToAnalyze='all',
  returnResult=TRUE,
  returnSummary=TRUE,
  onlySigIsoforms = T
)

pdf("./plots/splicing-genome-wide_strict.pdf", 15,10)
extractSplicingGenomeWide( SwitchListAnalyzed )
dev.off()

bioMechanismeAnalysis <- analyzeSwitchConsequences(
  SwitchListAnalyzed, 
  consequencesToAnalyze = 'all',#c('tss','tts','intron_structure'),
  showProgress = FALSE
)$switchConsequence

bioMechanismeAnalysis <- bioMechanismeAnalysis[which(bioMechanismeAnalysis$isoformsDifferent),]

myConsequences <- SwitchListAnalyzed$switchConsequence
myConsequences <- myConsequences[which(myConsequences$isoformsDifferent),]
myConsequences$isoPair <- paste(myConsequences$isoformUpregulated, myConsequences$isoformDownregulated)

write.csv(myConsequences, "./myconsequences.csv")

bioMechanismeAnalysis$isoPair <- paste(bioMechanismeAnalysis$isoformUpregulated, bioMechanismeAnalysis$isoformDownregulated)
bioMechanismeAnalysis <- bioMechanismeAnalysis[which(bioMechanismeAnalysis$isoPair %in% myConsequences$isoPair),]  

write.csv(bioMechanismeAnalysis, "./biomechanalysis.csv")

extractSplicingSummary( SwitchListAnalyzed )
extractSplicingEnrichment( SwitchListAnalyzed )
#extractSplicingEnrichmentComparison( SwitchListAnalyzed )
extractSplicingGenomeWide( SwitchListAnalyzed )

### Analyse (predicting) switch consequences
SwitchListAnalyzed <- analyzeSwitchConsequences( SwitchListAnalyzed,
                                                 consequencesToAnalyze = 'all')
extractSwitchSummary(SwitchListAnalyzed,filterForConsequences = F)
extractSwitchSummary(SwitchListAnalyzed,filterForConsequences = T) # the difference mean how many genes with switch have significant isoform switches functional consequences 


### post analysis of isoform switches with consequences
## analysis of individual isoforem switching
## we can extract the top switching genes by q-values
extractTopSwitches(
  SwitchListAnalyzed, 
  filterForConsequences = TRUE, 
  n = NA, 
  sortByQvals = TRUE)

## or we can extract the top switching genes by dIF values
extractTopSwitches(
  SwitchListAnalyzed, 
  filterForConsequences = TRUE, 
  n = NA, 
  sortByQvals = F)

switchingIso <- extractTopSwitches( 
  SwitchListAnalyzed, 
  filterForConsequences = TRUE, 
  n = NA,                  # n=NA: all features are returned
  extractGenes = FALSE,    # when FALSE isoforms are returned
  sortByQvals = TRUE
)

### plot the overlap genes with the det and deg
##"AIPGENE10029" "AIPGENE29161" "AIPGENE4696"  "AIPGENE15748" "AIPGENE21105" "AIPGENE22997" "AIPGENE9034"  "AIPGENE6206" 
# "AIPGENE4696"  "AIPGENE15748" "AIPGENE21105" "AIPGENE22997" "AIPGENE9034"  "AIPGENE6206"
switchPlot(SwitchListAnalyzed,gene = 'AIPGENE6206') # mix, GRINA, maybe differrnt network here # intersting 
switchPlot(SwitchListAnalyzed,gene = 'AIPGENE9034') # down, srsf4, mmd sensitivy # interesting
switchPlot(SwitchListAnalyzed,gene = 'AIPGENE22997')  # up, CTHRC1, different coding structure  # may report 
switchPlot(SwitchListAnalyzed,gene = 'AIPGENE21105') # up, v1g163628, turn into the non-coding thing.. # interesting 
switchPlot(SwitchListAnalyzed,gene = 'AIPGENE15748')  # up, MIP, both are coding
switchPlot(SwitchListAnalyzed,gene = 'AIPGENE4696') # up, totally different structure, transcript630 increase use in symbiotic, and upregulated 

switchPlot(SwitchListAnalyzed,gene = 'AIPGENE29161') 
switchPlot(SwitchListAnalyzed,gene = 'AIPGENE10029') 

write.csv(c("AIPGENE4696_gene","AIPGENE15748_gene","AIPGENE21105_gene","AIPGENE22997_gene","AIPGENE9034_gene","AIPGENE6206_gene"),'../../GO_analysis/overlap_genes_with_det_deg_dtu_with_consequence.txt',quote = F,row.names = F)


## genome wide summaries
extractConsequenceSummary(
  SwitchListAnalyzed,
  consequencesToAnalyze='all',
  plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
  asFractionTotal = FALSE      # enables analysis of fraction of significant features
)


### extract genome-wide analysis of isoform switching
extractConsequenceSummary(SwitchListAnalyzed)
extractSplicingSummary(SwitchListAnalyzed)

# check the splicing/consequence enrichment
extractConsequenceEnrichment(SwitchListAnalyzed) 
extractSplicingEnrichment(SwitchListAnalyzed)

# comparison of enrichment
extractConsequenceEnrichmentComparison(SwitchListAnalyzed)
extractSplicingEnrichmentComparison(SwitchListAnalyzed)

# analyze of genome-wide changes in isoform usage
extractConsequenceGenomeWide(SwitchListAnalyzed)
extractSplicingGenomeWide(SwitchListAnalyzed)


##### augmenting small upstreamn ORFs
# run ORF analysis on longest ORF
mySwitchList2 = SwitchListAnalyzed
mySwitchList2=extractSequence(mySwitchList2,genomeObject = './external_output/isoformSwitchAnalyzeR_isoform_nt.fasta')
mySwitchList2<- analyzeORF(mySwitchList2, orfMethod='longest')
mySwitchList2$orfAnalysis$orfTransciptLength[is.na(mySwitchList2$orfAnalysis$orfTransciptLength)] <- 0
mean(mySwitchList2$orfAnalysis$orfTransciptLength)
# run ORF analysis on most upstream ORF
mySwitchList2 <- analyzeORF(mySwitchList2, orfMethod = 'mostUpstream', minORFlength = 50)
mySwitchList2$orfAnalysis$orfTransciptLength[is.na(mySwitchList2$orfAnalysis$orfTransciptLength)] <- 0

mean(mySwitchList2$orfAnalysis$orfTransciptLength)

# calculate pairwise difference
summary(
  mySwitchList2$orfAnalysis$orfTransciptLength -
    mySwitchList2$orfAnalysis$orfTransciptLength[
      match(
        mySwitchList2$orfAnalysis$isoform_id,
        mySwitchList2$orfAnalysis$isoform_id
      )
    ]
)

table(mySwitchList2$isoformFeatures$codingPotential, exclude = NULL)
mySwitchList2$isoformFeatures$codingPotential <- NA
mySwitchList2$isoformFeatures$codingPotential[which(mySwitchList2$isoformFeatures$codingPotentialValue > 0.75)] <- TRUE
mySwitchList2$isoformFeatures$codingPotential[which(mySwitchList2$isoformFeatures$codingPotentialValue < 0.25)] <- FALSE

table(mySwitchList2$isoformFeatures$codingPotential, exclude = NULL)


### analyzing the biological mechanisms behind isoform switching
# there are three mechanisms: alternative transcription start site (aTSS); alternative splicing (AS); alternative transcription termination site (aTSS)
### analyze the biological mechanisms
bioMechanismeAnalysis <- analyzeSwitchConsequences(
  mySwitchList2, 
  consequencesToAnalyze = c('tss','tts','intron_structure'),
  showProgress = FALSE
)$switchConsequence # only the consequences are interesting here

### subset to those with differences
bioMechanismeAnalysis <- bioMechanismeAnalysis[which(bioMechanismeAnalysis$isoformsDifferent),]

### extract the consequences of interest already stored in the switchAnalyzeRlist
myConsequences <- mySwitchList2$switchConsequence
myConsequences <- myConsequences[which(myConsequences$isoformsDifferent),]
myConsequences$isoPair <- paste(myConsequences$isoformUpregulated, myConsequences$isoformDownregulated) # id for specific iso comparison

### Obtain the mechanisms of the isoform switches with consequences
bioMechanismeAnalysis$isoPair <- paste(bioMechanismeAnalysis$isoformUpregulated, bioMechanismeAnalysis$isoformDownregulated)
bioMechanismeAnalysis <- bioMechanismeAnalysis[which(bioMechanismeAnalysis$isoPair %in% myConsequences$isoPair),]  # id for specific iso comparison

### Create list with the isoPair ids for each consequence
AS   <- bioMechanismeAnalysis$isoPair[ which( bioMechanismeAnalysis$featureCompared == 'intron_structure')]
aTSS <- bioMechanismeAnalysis$isoPair[ which( bioMechanismeAnalysis$featureCompared == 'tss'             )]
aTTS <- bioMechanismeAnalysis$isoPair[ which( bioMechanismeAnalysis$featureCompared == 'tts'             )]

mechList <- list(
  AS=AS,
  aTSS=aTSS,
  aTTS=aTTS
)

### Create Venn diagram
library(VennDiagram)
#> Loading required package: grid
#> Loading required package: futile.logger
#> 
#> Attaching package: 'futile.logger'
#> The following object is masked from 'package:mgcv':
#> 
#>     scat
myVenn <- venn.diagram(
  x = mechList,
  col='transparent',
  alpha=0.4,
  fill=RColorBrewer::brewer.pal(n=3,name='Dark2'),
  filename=NULL
)

### Plot the venn diagram
grid.newpage() ; grid.draw(myVenn)

write.csv(bioMechanismeAnalysis, "biomechanalysis_subset.csv")
write.csv(myConsequences, "myconsequences_subset.csv")

#### overview plot
ggplot(data=mySwitchList2$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
    size=1
  ) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
  geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
  facet_wrap( ~ condition_2) +
  #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
  scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
  labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
  theme_bw()


### Switch vs Gene changes:
ggplot(data=mySwitchList2$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
    size=1
  ) + 
  facet_wrap(~ condition_2) +
  #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
  geom_hline(yintercept = 0, linetype='dashed') +
  geom_vline(xintercept = 0, linetype='dashed') +
  scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
  labs(x='Gene log2 fold change', y='dIF') +
  theme_bw()


### compare the condition difference behind the biological consequence of the three mechanism
library(dplyr)

# 1) get all significant isoform pairs (with consequences)
pairs <- mySwitchList2$switchConsequence %>%
  filter(isoformsDifferent) %>%
  dplyr::select(isoformUpregulated, isoformDownregulated) %>%
  mutate(
    isoPair = paste(isoformUpregulated, isoformDownregulated),
    # pick a neutral "anchor" isoform per pair (lexicographically smallest)
    anchor_iso = if_else(isoformUpregulated < isoformDownregulated,
                         isoformUpregulated, isoformDownregulated)
  ) %>%
  distinct()

# 2) bring in dIF and condition labels for that anchor isoform
anchor_if <- mySwitchList2$isoformFeatures %>%
  dplyr::select(isoform_id, dIF, condition_1, condition_2) %>%
  rename(anchor_iso = isoform_id)

dir_tbl <- pairs %>%
  left_join(anchor_if, by = "anchor_iso") %>%
  mutate(
    favored_condition = if_else(dIF >= 0, condition_2, condition_1)
  ) %>%
  dplyr::select(isoPair, favored_condition, dIF, condition_1, condition_2) %>%
  distinct()

# sanity check: should now have both conditions
table(dir_tbl$favored_condition)
summary(dir_tbl$dIF)




##############################################################
#### compare the gene level and transcriptome level
library(dplyr)
library(tibble)

tx_res=read.csv('../../DEG_analysis/diffexpr-results-isoquant-genes-6samples_update_from_transcript_model_transcript_level_raw.csv',header=T,row.names = 1)
gene_res=read.csv('../../DEG_analysis/diffexpr-results-isoquant-genes-6samples_update_from_transcript_model_gene_level_raw.csv',header=T,row.names = 1)


# --- Keep only valid padj rows, then filter significant sets ---
tx_res   <- tx_res   %>% filter(!is.na(padj))
gene_res <- gene_res %>% filter(!is.na(padj))

sig_txs   <- tx_res   %>% filter(padj < 0.05) %>% filter(abs(log2FoldChange)>1)
sig_genes <- gene_res %>% filter(padj < 0.05) %>% filter(abs(log2FoldChange)>1)

# --- tx2gene mapping ---
tx2gene <- read.table("../input_data/tx2gene.tsv", header = FALSE, sep = "\t",
                      col.names = c("transcript_id", "gene_id")) %>%
  distinct(transcript_id, gene_id)
#tx2gene$gene_id=gsub('_gene$','',tx2gene$gene_id)

# --- Genes that have >=1 significant transcript ---
genes_with_sig_tx <- tx2gene %>%
  semi_join(sig_txs, by = "transcript_id") %>%
  distinct(gene_id)

# --- Overlap metrics ---
n_DE_genes             <- nrow(sig_genes)
n_DE_tx                <- nrow(sig_txs)
colnames(sig_genes)[1]='gene_id'
n_DE_genes_with_sig_tx <- sig_genes %>% semi_join(genes_with_sig_tx, by = "gene_id") %>% nrow()
n_DE_genes_no_sig_tx   <- n_DE_genes - n_DE_genes_with_sig_tx

DE_genes_with_sig_tx <- sig_genes %>% semi_join(genes_with_sig_tx, by = "gene_id")
DE_genes_with_sig_tx <- DE_genes_with_sig_tx$gene_id
DE_genes_no_sig_tx   <- sig_genes[!(sig_genes$gene_id %in% DE_genes_with_sig_tx),]$gene_id

write.csv(DE_genes_no_sig_tx,'../../GO_analysis/DE_genes_no_sig_tx.txt',row.names = F,quote = F)
write.csv(DE_genes_with_sig_tx,'../../GO_analysis/DE_genes_with_sig_tx.txt',row.names = F,quote = F)

# Fraction of significant transcripts that belong to DE genes
sig_tx_in_DE_genes <- sig_txs %>%
  semi_join(tx2gene %>% semi_join(sig_genes, by = "gene_id"),
            by = "transcript_id")
n_sig_tx_in_DE_genes <- nrow(sig_tx_in_DE_genes)

# Optional: distribution of how many sig transcripts per DE gene
sig_tx_counts_per_gene <- tx2gene %>%
  semi_join(sig_txs, by = "transcript_id") %>%
  count(gene_id, name = "n_sig_tx")

cat("\n=== Overlap summary ===\n")
cat("Significant genes (FDR<0.05 & abs(logFC)>1):", n_DE_genes, "\n")
cat("Significant transcripts (FDR<0.05 & abs(logFC)>1):", n_DE_tx, "\n")
cat("DE genes with >=1 significant transcript:", n_DE_genes_with_sig_tx, "\n")
cat("DE genes with NO significant transcript:", n_DE_genes_no_sig_tx, "\n")
cat(sprintf("%% of DE genes with >=1 sig transcript: %.1f%%\n",
            100 * n_DE_genes_with_sig_tx / max(1, n_DE_genes)))
cat(sprintf("Sig transcripts that belong to DE genes: %d (%.1f%%)\n",
            n_sig_tx_in_DE_genes,
            100 * n_sig_tx_in_DE_genes / max(1, n_DE_tx)))

# Optional: peek at genes with >=2 sig transcripts
sig_tx_counts_per_gene %>% filter(n_sig_tx >= 2) %>% arrange(desc(n_sig_tx)) %>% head()


library(VennDiagram)
venn.plot <- draw.pairwise.venn(
  area1 = n_DE_genes,
  area2 = n_DE_tx,
  cross.area = n_DE_genes_with_sig_tx,
  category = c("DE Genes", "DE Transcripts"),
  fill = c("skyblue", "orange"), alpha = 0.5
)
grid.draw(venn.plot)
ggsave('../only_3_replicates_default/plots/venn_plot_deg_vs_de_transcript.png')

library(VennDiagram)
library(grid)

grid.newpage()
venn.plot <- draw.pairwise.venn(
  area1      = n_DE_genes,
  area2      = n_DE_tx,
  cross.area = n_DE_genes_with_sig_tx,
  category   = c("DE Genes", "DE Transcripts"),
  fill       = c("#3B82F6", "#F59E0B"),  # blue & amber
  alpha      = 0.5,
  lty        = "blank",
  cex        = 2,
  cat.cex    = 1.5,
  cat.pos    = c(-20, 20),
  cat.dist   = c(0.05, 0.05)
)
grid.draw(venn.plot)

head(gene_res)
head(tx_res)
head(tx2gene)

################ calculate the correlations between the DEG logFC and the DE transcript logFC ###################
#########################################################
tx_fc_summary <- tx_res %>%
  inner_join(tx2gene, by = "transcript_id") %>%
  group_by(gene_id) %>%
  summarize(mean_tx_LFC = mean(log2FoldChange, na.rm = TRUE),
            n_tx = n())

# Merge with gene-level table
colnames(gene_res)[1]='gene_id'
gene_tx_fc <- gene_res %>%
  inner_join(tx_fc_summary, by = "gene_id")

cor_test <- cor.test(gene_tx_fc$log2FoldChange, gene_tx_fc$mean_tx_LFC, method = "pearson")
cor_test # p-value < 2.2e-16, r:0.979


library(ggplot2)

r_val <- unname(cor_test$estimate)
p_val <- cor_test$p.value

p <- ggplot(gene_tx_fc, aes(x = mean_tx_LFC, y = log2FoldChange)) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", linewidth = 0.7, color = "grey50") +
  geom_point(size = 2.2, alpha = 0.6, shape = 16) +
  geom_smooth(method = "lm", color = "#D62728", se = TRUE,
              linewidth = 1) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = sprintf("Pearson r = %.2f\nP = %.2e", r_val, p_val),
    hjust = 1.1, vjust = -0.5,
    size = 6
  ) +
  labs(
    x = "Mean transcript-level log2 fold change",
    y = "Gene-level log2 fold change",
    title = "Correlation between gene- and transcript-level changes"
  ) +
  coord_equal() +
  theme_classic(base_size = 24) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.7)
  )

p

ggsave(
  "./plots/correlation_of_mean_logFC_between_gene_transcript.png",
  plot = p,
  width = 6,
  height = 6,
  dpi = 600,
  bg = "white"
)

## find genes with opposite-direction isoforms
tx_fc_signs <- tx_res %>%
  inner_join(tx2gene, by = "transcript_id") %>%
  group_by(gene_id) %>%
  summarize(
    n_isoforms = n(),
    n_pos = sum(log2FoldChange > 0),
    n_neg = sum(log2FoldChange < 0),
    has_both_directions = (n_pos > 0 & n_neg > 0)
  )

isoform_switch_genes <- tx_fc_signs %>% filter(has_both_directions)
nrow(isoform_switch_genes)
head(isoform_switch_genes)

isoform_switch_annot <- isoform_switch_genes %>%
  inner_join(gene_res, by = "gene_id") %>%
  arrange(desc(abs(log2FoldChange))) %>%
  dplyr::select(gene_id, log2FoldChange, n_isoforms, n_pos, n_neg)
write_csv(isoform_switch_annot,'isoform_switch_anno_from_count.csv')

# plot the example gene
example_gene <- "AIPGENE6206_gene"  # replace with a real one
tx_res %>%
  inner_join(tx2gene, by = "transcript_id") %>%
  filter(gene_id == example_gene) %>%
  ggplot(aes(x = transcript_id, y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_col() +
  scale_fill_manual(values = c("red", "blue")) +
  labs(title = paste("Isoform expression changes in", example_gene),
       x = "Transcript ID", y = "log2FoldChange") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
