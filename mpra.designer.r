#!/usr/bin/env Rscript


### To do
# - Unify mask format into a comma-separated range not data.frame
# - add scramble
# - drop anchors with "None" motifs
# - design parameter sheet in the xlsx file
# - frequency of each TF combination
#	before & after random selection of PRDM1-only 
# - freaeuency of tf combination in conjunction with motif
# IMPORTANT: if downsampling is to be done (per TF combinations), motif scanning better be done before the downssampling
#	because of any missing motif instance due to higher marginal score
# - validation of meme motif names using read_meme in universalmotif


library(tidyverse)
library(Biostrings)
library(openxlsx)
library(ggplot2)
#library(PWMEnrich)
library(dplyr)
library(universalmotif)
library(Biostrings)
library(data.table)
library(ggplot2)
library(reshape2)
library('yaml', quiet=TRUE)

source(sprintf("%s/commonR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/motifR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/basicR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/genomeR.r", Sys.getenv("COMMON_LIB_BASE")))
source("../MPRA.common.r")

config = yaml.load_file("./config.MPRA.yml")


## Required config elements
configElemL=c(
	"main_title",
	"src_anchor",
	"src_cobind",
	"src_motif",
	"motifNameL",
	"adapter5",
	"adapter3",
	"scanWindow",
	"genome",
	"mutationMode",
	"pv1",
	"pv2"
)

for( elem in configElemL ) stopifnot(elem %in% names(config))

## Config elements and validation
main_title 	= config[["main_title"]]
src.anchor 	= config[["src_anchor"]]
src.cobind	= config[["src_cobind"]]
src.motif 	= config[["src_motif"]]
motifNameL	= config[["motifNameL"]]
downsampleL = config[["downsample"]]
adapter5	= config[["adapter5"]]
adapter3	= config[["adapter3"]]
scanWindow	= config[["scanWindow"]]
genome		= config[["genome"]]
mutationMode = config[["mutationMode"]]
pv1 		= config[["pv1"]]
pv2 		= config[["pv2"]]

stopifnot(pv2 >= pv1)
if(!is.null(downsampleL)) stopifnot(all(unlist(downsampleL) > 0))

##############################
## Destination files
des.log = "log.txt"

des.scan.fa = "anchor.scan.fa"
des.scan.bed = "anchor.scan.bed"
des.fimo = "anchor.fimo.txt"
des.motif = "motif.meme"

des.comboBarplot = "barplot.TF_combo.pdf"
desPrefix.maskOverlap = "maskOverlap"

desPrefix = "mpra"


# initial counts
N.motif = length(motifNameL)



## Start logging
writeLog("=====================================
MPRA-ChIP:Designer
=====================================", des.log, append=FALSE)

writeLog(sprintf("- Title: %s
- Anchor bed file: %s
- TF overlap table: %s",
	main_title,
	src.anchor,
	src.cobind),
	des.log, append=TRUE) 

writeLog(sprintf("- Motif file: %s
- Motif names:
%s
- Motif scan window: %d bp
- Genome: %s
- P-value
\tStringent: %.e
\tMarginal: %.e",
	src.motif,
	paste(sprintf("\t%s", motifNameL), collapse="\n"),
	scanWindow,
	genome,
	pv1,
	pv2),
	des.log, append=TRUE) 


writeLog(sprintf("- Downsampling
%s",
	paste(sprintf("\t%s: %d", names(downsampleL), unlist(downsampleL)), collapse="\n")),
	des.log, append=TRUE) 

writeLog(sprintf("- Mutation mode: %s
- Adapter
\t5': %s
\t3': %s,",
	mutationMode,
	pv1,
	pv2),
	des.log, append=TRUE)
writeLog("", des.log, append=TRUE)




###############################
## Step 0: TF overlap Boolean table & MEME motif validation
# Note: need to reimplemented in a generalized form

# overlap table
writeLog(sprintf("%s: Reading TF overlap table - %s", Sys.time(), src.cobind), des.log, append=TRUE)
data.cobind = read.delim(src.cobind, header=TRUE, row.names=4, stringsAsFactors=FALSE)
df.cobind = data.frame(data.cobind[,6:ncol(data.cobind)] > 0)
rownames(df.cobind) = rownames(data.cobind)
#anchorNameL = rownames(data.cobind)
stopifnot(all(motifNameL %in% colnames(df.cobind)))
df.cobind = df.cobind[,motifNameL]

# MEME moti names & filtering
writeLog(sprintf("%s: Filtering motifs", Sys.time()), des.log, append=TRUE)
motifL = read_meme(src.motif)
stopifnot(all(motifNameL %in% sapply(motifL, function(x) x@altname)))
writeLog(sprintf("%s: Saving %s", Sys.time(), des.motif), des.log, append=TRUE)
motifL = filter_motifs(motifL, altname=motifNameL)
write_meme(motifL, "motif.meme", overwrite=TRUE)


#writeLog(sprintf("%s: Writing TF overlap stat. - %s", Sys.time(), "factorCobindStat.txt"), des.log, append=TRUE)
#df.tfStat = apply(df.cobind, 2, sum)
#write.table( data.frame(df.tfStat), "factorCobindStat.txt", row.names=TRUE, col.names=FALSE, quote = FALSE, sep="\t" )



## anchor bed and sequence
# load / resize / reorder (to match cobind)
writeLog(sprintf("%s: Resizing scan window and extracting DNA sequences - %s", Sys.time(), des.scan.fa), des.log, append=TRUE)
scan.bed = resizeBed(readBedFile(src.anchor), scanWindow)
colnames(scan.bed) = c("Chr","Start","End","Name","Null","Direc")
stopifnot(setequal(rownames(df.cobind), scan.bed[,4]))
rownames(scan.bed) = scan.bed$Name
anchorNameL = scan.bed$Name
# reordering of cobind table to match with anchor
df.cobind = df.cobind[anchorNameL,]
writeBedFile(scan.bed, des.scan.bed)

if(file.exists(des.scan.fa)){
	scan.fa = readDNAStringSet(des.scan.fa)
	stopifnot(all(names(scan.fa) == anchorNameL))
}else{
	scan.fa = extractDNA(scan.bed, genome)
	scan.fa = scan.fa[anchorNameL]
	writeXStringSet(scan.fa, des.scan.fa)
}




###############################
## Step 2: Initial motif scan

# Initial FIMO Output:
#	a list of 
#	- instance: data.frame of identified motif instance in a format
#		Name: sequence name to scan
#		Offset: starting location (1-base) within the sequence
#		Score: motif score
#		Direc: motif direction +/-
#		Seq: identified motif instance sequence
#	- status: vector of scan status, same length with the number of input sequence
#		the status vector has one of Stringent / Marginal / None
#		each element denote the status of the motif scan
#			Stringent: stringent instance found with score > cutoff1, can be more than one instances
#			Marginal: single instance of motif found with maximum score > cutoff2
#			None: no instance found i.e. score is < cutoff2 in all position / strand
# Action
# - drop marginal instances when an anchor has stringent ones
# - select best marginal one(s) if no stringent instances per anchor
writeLog(sprintf("%s: Fimo motif scan - %s", Sys.time(), des.fimo), des.log, append=TRUE)
fimoL = scanMemeMPRA(des.scan.fa, des.motif, pv1, pv2, motifNameL, des.fimo=des.fimo, verbose=FALSE)

## reformatting of FIMO results
# List of motif --> Data.frame
#	- instance: N x M data.frame where N is the number of anchors and M is the number of motifs
#	- status: data.frame of the same dimension denoting the status of the motif, i.e. stringent / marginal / none
writeLog(sprintf("%s: Converting FIMO results to data.frame", Sys.time()), des.log, append=TRUE)
df.fimo = makeFimoDataFrame(fimoL, anchorNameL, motifNameL)



###############################
## Step 3: Filtering motif instances by TF overlap
writeLog(sprintf("%s: Intersecting cobinding x motif scan", Sys.time()), des.log, append=TRUE)
## Intersection of cobind x motif scan
# - Ignore motif instances without cobinding
# - drop anchors without desired motif but only with cobinding
stopifnot(all(rownames(df.cobind) == rownames(df.fimo[["status"]])))
result.intersect = intersectCobindMotif(df.cobind, df.fimo)
df.cobind.filtered = result.intersect[["cobind"]]
df.fimo.filtered = result.intersect[["motifScan"]]
comboNameL = boolTableToVector(df.cobind.filtered)

## Frequency of combinations before downsampling
freq.cobind = getComboFrequency(df.cobind)
freq.intersect = getComboFrequency(df.cobind.filtered)



## Optional downsampling
## Note: PERMANENT CHAGE
#	- df.cobind.filered
#	- df.fimo.filtered
# Optional subsequent sorting to preserve the original order, but not implemented yet

# anchor name vector to use down the road
anchorNameL.filtered = rownames(df.cobind.filtered)
if( ! is.null(downsampleL) ){
	writeLog(sprintf("%s: Downsampling\n%s",
				Sys.time(),
				paste(sprintf("\t- %s: %d", names(downsampleL), unlist(downsampleL)), collapse="\n")),
			des.log, append=TRUE)
	df.cobind.filtered = downsampleCombo(df.cobind.filtered, downsampleL)
	anchorNameL.filtered = rownames(df.cobind)[rownames(df.cobind) %in% rownames(df.cobind.filtered)]
	df.cobind.filtered = df.cobind.filtered[anchorNameL.filtered,]
	df.fimo.filtered[["instance"]] = df.fimo.filtered[["instance"]][anchorNameL.filtered,]
	df.fimo.filtered[["status"]] = df.fimo.filtered[["status"]][anchorNameL.filtered,]
}
# DNA sequence after downampling
scan.fa.filtered = scan.fa[rownames(df.cobind.filtered)]


writeLog(sprintf("%s: Visualizing TF combination counts - %s", Sys.time(), des.comboBarplot), des.log, append=TRUE)
# Frequencies after downsampling
freq.down = getComboFrequency(df.cobind.filtered)
stopifnot(all(names(freq.cobind) == names(freq.intersect)))
stopifnot(all(names(freq.cobind) == names(freq.down)))
df.freq = data.frame(Cobind = freq.cobind, Intersect = freq.intersect, Downsample = freq.down)

## combo frequency barplot
N.crs = calcCrsCount(setNames(df.freq[,"Downsample"], rownames(df.freq)))
g = drawComboFreqPlot(df.freq) +
		labs(title = main_title, 
			subtitle=sprintf("Frequency of combinations with TF\nNumber of expected CRS = %d", N.crs))
ggsave(des.comboBarplot, g, width=8, height=8)




###############################
## Step 4: Mask generation and mutation

# Convert fimo motif instances to a list of mask data.frames
writeLog(sprintf("%s: Converting motif instances to mask", Sys.time()), des.log, append=TRUE)
maskL = instanceFimoToMask(df.fimo.filtered[["instance"]])

## Assess overlap between motif mask sets of two TFs
checkMaskOverlap(maskL, desPrefix.maskOverlap)



################################
# Step 5: Combinatorial mutation
writeLog(sprintf("%s: Generating mutations", Sys.time()), des.log, append=TRUE)
# prepare combination Boolean table for various numbers of TFs
comboBoolTableL=list()
for( i in 1:N.motif ) comboBoolTableL[[i]] = makeComboBoolTable(i)

# Group masks by anchor names
maskL.groupByAnchor = groupMaskByAnchor(maskL)

# Generate mutated sequence
df.crs = generateMutatedSeqSet(scan.fa.filtered, maskL.groupByAnchor, comboBoolTableL, mutationMode, maxTryDiff=10, validate=TRUE)
# Assitn unique id
rownames(df.crs) = sprintf("CRS%05d", 1:nrow(df.crs))

if(FALSE){
	# test code
	wt = 1
	mut = 4
	compareTwoStrings(df.crs[wt,"Sequence"], df.crs[mut,"Sequence"], df.crs[mut,"Mask"], width=100)
}



#####################################
## Step 6: Motif score compare before and after mutation
writeLog(sprintf("%s: Motif scan in designed sequences", Sys.time()), des.log, append=TRUE)

## Mutated max motif score, including WT
maxMotifScore.crs = getMaxMemeScore(setNames(df.crs$Sequence, rownames(df.crs)), des.motif, motifNameL, pv=0.01, des="crs.fimo.txt")
stopifnot(all(rownames(maxMotifScore.crs) == rownames(df.crs)))

writeLog(sprintf("%s: Comparing motif scores before vs after mutation", Sys.time()), des.log, append=TRUE)
# Note:The visualization functions below was from previous implementation of the motif score data.frame
# 	where the rownames were the combiniation of anchor name and motif list.
# 	Thus, the motif score data.frame is temporarily assigned with rownames to reuse the functions.
tmp.score = maxMotifScore.crs
rownames(tmp.score) = sprintf("%s:%s", df.crs$AnchorName, df.crs$MotifMutated)
# motif score: WT vs Mutant for all
compareMotifScoreAll("motifScoreCompareAll", tmp.score,
		main="Motif Score WT vs Mut: All", mode="all")
compareMotifScoreAll("motifScoreCompareAll.Target", tmp.score,
		main="Motif Score WT vs Mut: Targeted Motifs", mode="target")
compareMotifScoreAll("motifScoreCompareAll.Nontarget", tmp.score,
		main="Motif Score WT vs Mut: Nontarget Motifs", mode="nontarget")
compareMotifScoreAll("motifScoreCompareAll.Scramble", tmp.score,
		main="Motif Score WT vs Scramble", mode="scramble")

# motif score compare for each TF before vs after single TF motif mutation
compareMotifScoreSingleMut("motifScoreComparePair.singleOnly", tmp.score,
		main = "Motif Score WT vs Mut: Single Mutation Only", mode = "exclusive")
compareMotifScoreSingleMut("motifScoreComparePair.all", tmp.score,
		main = "Motif Score WT vs Mut: Including Multiple Mutations", mode="all")
compareMotifScoreScramble("motifScoreComparePair.Scramble", tmp.score,
		main="Motif Score WT vs Scramble")



#################################
## Step 7: output prep
writeLog(sprintf("%s: Exporting results", Sys.time()), des.log, append=TRUE)

exportDesign(desPrefix, config,
	anchor = scan.bed[anchorNameL.filtered,],
	crs = df.crs,
	cobind = df.cobind.filtered[anchorNameL.filtered,],
	motifStatus.raw = df.fimo[["status"]][anchorNameL.filtered,],
	motifStatus.filtered = df.fimo.filtered[["status"]][anchorNameL.filtered,],
	maxMotifScore = maxMotifScore.crs)




q()










########################################################
## Debug/test code
if(FALSE){
	tmp1 = boolTableToVector(df.mask!="NA")
	tmp2 = df.cobind.select
	colnames(tmp2) = sub("Cobind.","",names(tmp2))
	tmp2 = boolTableToVector(tmp2)

	tmp1 = names(tmp1)[tmp1=="PRDM1"]
	tmp2 = names(tmp2)[tmp2=="PRDM1"]
	name.diff = setdiff(tmp1, tmp2)

	df.cobind[name.diff,]
	df.motifStatus.raw[name.diff,]
	df.mask[name.diff,]
	# df.cobind
	# df.cobind.select
	# df.motifStatus.raw
	# #df.motifStatus.filtered  
	# df.mask
	#stopifnot(all( (df.motifStatus.filtered != "None") == (df.mask != "NA") ))
	stopifnot( all(rownames(df.cobind.select) == rownames(df.motifStatus.raw)) )
	stopifnot( all(rownames(df.cobind.select) == rownames(df.mask)) )
}


### adapter sequence motif score
if(FALSE){
	seq1="AGGACCGGATCAACT"
	seq2="CATTGCGTGAACCGA"

	tmp = getMaxMotifScore(c(seq1, seq2), srcL.motif)
}