###############################################################################
#####################  Importing and filtering data  ##########################
###############################################################################
# We need to import the functions in the asynt.R script.
# You need to have the Intervals package installed on your system for this to work
setwd("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\Paper_GRC_transcription\\GENETICS_submission\\Revisions\\04_HGT\\asynt")
home <- getwd()
home

###############################################################################
########################  single scaffold alignment plot  #####################
###############################################################################
# We need to import the functions in the asynt.R script.
# You need to have the Intervals package installed on your system for this to work
source("asynt.R")

#if we have a single genomic region we are interested in, we can visualise the alignments diretcly
# First import the alignment data
alignments <- import.paf("grc2_vs_rickettsia_1to1.1aln.paf.gz")
alignments
# note that you could also have used an alignment generated using mummer
# with the nucmer and show-coords tools), using import.nucmer()

# Now we subset by scaffold to just the  reference and query sequences we are interested in.
alignments <- subset(alignments, Rlen >= 100 & query=="SUPER_GRC2" & reference == "ptg000159c")
alignments
ref_data <- import.genome(fai_file="Rickettsiaceae.finalassembly.fa.fai")

alignments <- reverse.references(alignments, reference_lens=ref_data$seq_len, ref_names="ptg000159c")

#and make the plot
par(mar = c(2,0,1,0))
plot.alignments(alignments, sigmoid=TRUE, tick_spacing=1e4, las=1)

#focus in on a specific region by setting the first and last base in the reference and query
plot.alignments(alignments, sigmoid=TRUE, Rfirst=0, Rlast =1753321, Qfirst=100002455, Qlast=101755776, tick_spacing=100000, las=1)

#show the entire alignemnt 
#plot.alignments(alignments, sigmoid=TRUE, Rfirst=0, Rlast =1753321, Qfirst=0, Qlast=101162943, tick_spacing=10000, las=2)
#plot.alignments(alignments, sigmoid=TRUE, Rfirst=0, Rlast =1753321, Qfirst=0, Qlast=111327891, tick_spacing=10000, las=2)

# Total bacteria length 
bac_length <- 1753321+98675   
GRC2_length <- 111327891
bac_length
#show the entire alignemnt 
plot.alignments(alignments, sigmoid=TRUE, Rfirst=0, Rlast =bac_length, Qfirst=0, Qlast=GRC2_length, tick_spacing=1000000, las=2)

801942*0.06575
801942*0.10959

###############################################################################
######################  multiple scaffold alignment plot  #####################
###############################################################################
# We need to import the functions in the asynt.R script.
# You need to have the Intervals package installed on your system for this to work
source("asynt.R")

# If we have multiple scaffolds making up a chromosome,
# we can string them together, either automatically or manually

#import alignments
alignments <- import.paf("grc2_vs_rickettsia_1to1.1aln.paf.gz")

#Next we import load scaffold length data which is necessary to plot the scaffolds
ref_data <- import.genome(fai_file="Rickettsiaceae.finalassembly.fa.fai")
query_data <- import.genome(fai_file="idBraCopr2.1.primary.masked.fa.fai")

#now define the scaffolds we're interested in
query_scafs <- "SUPER_GRC2"
reference_scafs <- c("ptg000159c", "ptg000322c")

#keep only alignments involving these scaffolds
alignments <- subset(alignments, query %in% query_scafs & reference %in% reference_scafs)

#plot alinments
par(mar = c(4,4,4,0), xpd=NA)
plot.alignments.multi(alignments, reference_lens=ref_data$seq_len, query_lens=query_data$seq_len, sigmoid=T)
