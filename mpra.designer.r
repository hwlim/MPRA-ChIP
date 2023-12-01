#!/usr/bin/env Rscript


### To do
# - Unify mask format into a comma-separated range not data.frame

library(tidyverse)
library(Biostrings)
library(openxlsx)

source(sprintf("%s/commonR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/basicR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/motifR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/genomeR.r", Sys.getenv("COMMON_LIB_BASE")))
source("MPRA.common.r")


#source(sprintf("%s/scRNAseq/seurat.common.r", Sys.getenv("LIMLAB_BASE")))

src.anchor = "Data.Peak/DE_Ctrl_PRDM1.bed"
src.overlap = "Data.Peak/Overlap_PF.txt"
srcL.motif = list()
srcL.motif[["PRDM1"]] = "Data.Motif/prdm1.motif"
srcL.motif[["EOMES"]] = "Data.Motif/eomes.motif"
srcL.motif[["FOXA"]] = "Data.Motif/foxa1_1_FG.motif"
srcL.motif[["GATA"]] = "Data.Motif/gata4.motif"
srcL.motif[["SOX17"]] = "Data.Motif/sox17.motif"


# motif scan window size and genome
scanWindow = 270
genome = "hg38"
mutationMode = "shuffle"
# fraction of PRDM1-only peaks to retain
anchorOnlyFrac = 0.1

src.anchor.fa = sprintf("DE_Ctrl_PRDM1.%sbp.fa", scanWindow)
src.anchor.bed = sprintf("DE_Ctrl_PRDM1.%sbp.bed", scanWindow)
des.log = sprintf("log.txt")



## Start logging
writeLog("=====================================
MPRA-ChIP:Designer
=====================================", des.log, append=FALSE)

writeLog(sprintf("
- anchor bed file: %s
- TF overlap table: %s
- Motif list:
%s
- Motif scan window: %d bp
- Genome: %s
- Mutation mode: %s
- Anchor TF only fracdtion: %.02f",
	src.anchor, src.overlap,
	paste(sprintf("\t%s : %s", names(srcL.motif), unlist(srcL.motif)), collapse="\n"),
	scanWindow, genome, mutationMode, anchorOnlyFrac), 
	des.log, append=TRUE)
writeLog("", des.log, append=TRUE)




###############################
## Step 0: TF overlap Boolean table
# Note: need to reimplemented in a generalized form
writeLog(sprintf("%s: Reading TF overlap table - %s", Sys.time(), src.overlap), des.log, append=TRUE)
data.overlap = read.delim(src.overlap, header=TRUE, row.names=4, stringsAsFactors=FALSE)
df.cobind = data.frame(PRDM1 = TRUE,
						EOMES = data.overlap$EOMES > 0,
						FOXA = data.overlap$FOXA2_maki > 0,
						GATA = data.overlap$GATA4_DE_WT | data.overlap$GATA6_DE_WT,
						SOX17 = data.overlap$SOX17 > 0)
rownames(df.cobind) = rownames(data.overlap)
writeLog(sprintf("%s: Writing TF overlap stat. - %s", Sys.time(), "factorCobindStat.txt"), des.log, append=TRUE)
df.tfStat = apply(df.cobind, 2, sum)
write.table( data.frame(df.tfStat), "factorCobindStat.txt", row.names=TRUE, col.names=FALSE, quote = FALSE, sep="\t" )


## Random selection of anchor-only sites
if(anchorOnlyFrac < 1){
	writeLog(sprintf("%s: Selecting party of anchor-only sites randomly", Sys.time()), des.log, append=TRUE)
	anchorOnly.all = names(which(apply(df.cobind, 1, function(x) all(x[-1]==FALSE))))
	anchorOnly.select = sample(anchorOnly.all, length(anchorOnly.all)*anchorOnlyFrac)
	anchorOnly.discard = setdiff(anchorOnly.all, anchorOnly.select)
}else{
	writeLog(sprintf("%s: Using all anchor-only sites ", Sys.time()), des.log, append=TRUE)
	anchorOnly.all = rownames(df.cobind)
	anchorOnly.select = anchorOnly.all
	anchorOnly.discard = NULL
}
anchorName.select = setdiff(rownames(df.cobind), anchorOnly.discard)
writeLog(sprintf("%s: Total %d anchor selected", Sys.time(), length(anchorName.select)), des.log, append=TRUE)
df.cobind.select = df.cobind[anchorName.select,]



###############################
## Step 1: 
# motif scan window
writeLog(sprintf("%s: Reading anchor file - %s", Sys.time(), src.anchor), des.log, append=TRUE)
bed.anchor= readBedFile(src.anchor)
rownames(bed.anchor) = bed.anchor[,4]
colnames(bed.anchor) = c("Chr","Start","End","AnchorName","Score","Direc")
stopifnot(setequal(bed.anchor[,4], rownames(df.cobind)))
bed.anchor.select = bed.anchor[anchorName.select,]
bed.anchor.select = resizeBed(bed.anchor.select, scanWindow)
# 4th column must be unique value
writeLog(sprintf("%s: Total number of anchor regions - %d", Sys.time(), nrow(bed.anchor)), des.log, append=TRUE)
writeBedFile(bed.anchor.select, src.anchor.bed)



# DNA sequence to scan the motif
if(file.exists(src.anchor.fa)){
	writeLog(sprintf("%s: Reading existing anchor fasta file - %s", Sys.time(), src.anchor.fa), des.log, append=TRUE)
	data.fa = readDNAStringSet(src.anchor.fa)
	if(!setequal(names(data.fa), bed.anchor.select[,4])){
		writeLog(sprintf("%s: Exsiting fasta file doesn't match with anchor bed file; re-generating", Sys.time(), src.anchor.fa), des.log, append=TRUE)
		data.fa = extractDNA(bed.anchor.select, genome, des=src.anchor.fa)
	}
}else{
	writeLog(sprintf("%s: Extracting genomic sequence and creating %s", Sys.time(), src.anchor.fa), des.log, append=TRUE)
	data.fa = extractDNA(bed.anchor.select, genome, des=src.anchor.fa)
}
# reorder to match bed.anchor.select
data.fa=data.fa[anchorName.select]

# tmp1 = table(boolTableToVector(df.cobind))
# tmp2 = table(boolTableToVector(df.cobind.select))

# tmp2 = sort(table(boolTableToVector(df.cobind.select )), decreasing=TRUE)
# writeLog(sprintf("Count of TF combinations
# %s", paste(sprintf("\t%s - %d", names(tmp2), tmp2), collapse="\n")),
# 	des.log, append=TRUE)

# barplot(matrix(cbind(as.numeric(tmp1), as.numeric(tmp2))))

# counts
N.motif = length(srcL.motif)
motifNameL = names(srcL.motif)
N.anchor = nrow(bed.anchor.select)





###############################
## Step 2: Initial motif scan
# Output:
#	resultL[[motifName]][[instance]]:	data.frame of identified motif instance
#	resultL[[motifName]][[status]]:		vector of motif scan status in all achor windows
motifScanL=list()
for( motifName in names(srcL.motif) ){
	# motifName = "PRDM1"
	write(sprintf("Scanning %s: %s", motifName, srcL.motif[[motifName]]), stderr())
	motif = readMotif.Homer(srcL.motif[[motifName]])
	#if(!is.null(motifScanL[[motifName]])) next
	motifScanL[[motifName]] = scanMotifMPRA(data.fa, motif$PWM, motif$threshold, 0, parallel=4, verbose=TRUE)
}

## motif scan stat
# Write a motif scan statistics table
df.stat = motifScanListToCntTable(motifScanL)
write.table(df.stat, sprintf("motifStat.0.prelim.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

# debug code
if(FALSE){
	# number of identified motif instances stratified by types
	sapply(motifScanL, function(x) table(x$status))
	# number of anchors that has at least one motif either stringent or marginal
	sapply(motifScanL, function(x) length(unique(x$instance$Name)))
}


###############################
## Step 2: Filtering motif instances by TF overlap
# Select motif instances corresponding to TF overlap only
motifScanL.filtered = filterMotifByTF(motifScanL, df.cobind.select)
df.stat.filtered = motifScanListToCntTable(motifScanL.filtered)
write.table(df.stat.filtered, sprintf("motifStat.1.withTF.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")




## data.frame
df.motifStatus.raw = motifScanListToStatusTable(motifScanL)
df.motifStatus.filtered = motifScanListToStatusTable(motifScanL.filtered)




###############################
# Motif mask generation in 1-base cooridnate
maskL.stringent = list()
maskL = list()
for( name in names(motifScanL.filtered) ){
	# name = "EOMES"
	instance = motifScanL.filtered[[name]]$instance
	status = motifScanL.filtered[[name]]$status
	instance.stringent = instance[status[instance$Name]=="Stringent",]

	maskL.stringent[[name]] = instanceToMask(instance.stringent)
	maskL[[name]] = instanceToMask(instance)
}

###############################
## Assess overlap between motif mask sets of two TFs
# 1. between stringent only
checkMaskOverlap(maskL.stringent, "maskOverlap.stringent")
# 2. between stringent + marginal
checkMaskOverlap(maskL, "maskOverlap.both")


# ###############################
# # Combinatorial mutation

# prepare combination Boolean table for various numbers of TFs
comboBoolTableL=list()
for( i in 1:N.motif ) comboBoolTableL[[i]] = makeComboBoolTable(i)

# Group masks by anchor names
maskL.groupByAnchor = groupMaskByAnchor(maskL)
df.mask = maskListToTable(maskL.groupByAnchor)[anchorName.select,]

# Generate mutated sequence
df.mut = generateMutatedSeqSet(data.fa, maskL.groupByAnchor, comboBoolTableL, mutationMode, maxTryDiff=5, validate=TRUE)
#compareTwoStrings(df.mut[1,"Sequence"], df.mut[6,"Sequence"], df.mut[6,"Mask"], width=100)


# temporary test code comparing sequences before and after mutation
if(FALSE){
	compareTwoSequence(as.character(data.fa[1]), df.mut[1,"Sequence"])
	compareTwoSequence(as.character(data.fa[1]), data.seq[1])
}

## Motif score


# ## Original max motif score
# maxMotifScore.raw = NULL
# for( motifName in names(motifScanL) ){
# 	if(is.null(maxMotifScore.raw)){
# 		maxMotifScore.raw = motifScanL[[motifName]]$maxScore0
# 	}else{
# 		maxMotifScore.raw = cbind(maxMotifScore.raw, motifScanL[[motifName]]$maxScore)
# 	}
# }
# colnames(maxMotifScore.raw) = names(motifScanL)

## Mutated max motif score
maxMotifScore.mut = getMaxMotifScore(df.mut$Sequence, srcL.motif, parallel=4)
rownames(maxMotifScore.mut) = rownames(df.mut)

# head(maxMotifScore.raw)
head(maxMotifScore.mut,20)


# motif score: WT vs Mutant for all
compareMotifScoreAll("motifScoreCompareAll", maxMotifScore.mut, main="Motif Score WT vs Mut: All", mode="all")
compareMotifScoreAll("motifScoreCompareAll.Target", maxMotifScore.mut, main="Motif Score WT vs Mut: Targeted Motifs", mode="target")
compareMotifScoreAll("motifScoreCompareAll.Nontarget", maxMotifScore.mut, main="Motif Score WT vs Mut: Nontarget Motifs", mode="nontarget")

# motif score compare for each TF before vs after single TF motif mutation
compareMotifScoreSingleMut("motifScoreComparePair.singleOnly",
		maxMotifScore.mut,
		main = "Motif Score WT vs Mut: Single Mutation Only",
		mode = "exclusive")
compareMotifScoreSingleMut("motifScoreComparePair.all",
		maxMotifScore.mut,
		main = "Motif Score WT vs Mut: Including Multiple Mutations",
		mode="all")





data.mut = data.frame(Id = rownames(df.mut), df.mut, maxMotifScore.mut)

names(df.cobind.select) = sprintf("Cobind.%s", names(df.cobind.select))
names(df.motifStatus.raw) = sprintf("Motif.%s", names(df.motifStatus.raw))
names(df.motifStatus.filtered) = sprintf("Motif.Cobind.%s", names(df.motifStatus.filtered))
df.cobind.select = df.cobind.select[anchorName.select,]
df.motifStatus.raw = df.motifStatus.raw[anchorName.select,]
df.motifStatus.filtered = df.motifStatus.filtered[anchorName.select,]
data.anchor = data.frame(
	bed.anchor.select,
	df.cobind.select,
	df.motifStatus.raw,
	df.motifStatus.filtered
)



######################################################
## Output


wb = createWorkbook("MPRA-ChIP")
addWorksheet(wb, "Anchor")
addWorksheet(wb, "Design")
writeData(wb, sheet = 1, data.anchor)
writeData(wb, sheet = 2, data.mut)
saveWorkbook(wb, "mpra.design.xlsx", overwrite = TRUE)

# plan
# - design parameter sheet in the xlsx file
# - frequency of each TF combination
#	before & after random selection of PRDM1-only 
# - freaeuency of tf combination in conjunction with motif