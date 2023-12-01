

##################################################################################################################
## Utility functions for MPRA-ChIP:design
##################################################################################################################

## Convert a given mask data.frame into a comma-separated range
#	e.g.
#	Input:
#		10	20
#		30	40
#	Output:
#		10-20,30-40
maskToRange=function(mask){
	# mask=maskL.groupByAnchor[[1]][[2]]
	range = apply(mask, 1, function(x) sprintf("%d-%d", x[1], x[2]))
	range = paste(range, collapse=",")
	return(range)
}

## Convert a given comma-separated range to mask data.frame
#	e.g.
#	Input:
#		10-20,30-40
#	Output:
#		10	20
#		30	40
rangeToMask=function(range){
	mask = do.call(rbind, strsplit(strsplit(range, ",")[[1]],"-"))
	mask = sapply(data.frame(mask), as.numeric)
	colnames(mask) = c("Start","End")
	return(mask)
}


## Convert mask list groupped by anchor into a data.frame
maskListToTable=function(maskList){
	# maskList = maskL.groupByAnchor
	motifNameL = unique(unlist(sapply(maskList, function(x) names(x))))
	df = data.frame(matrix("NA", length(maskList), length(motifNameL)))
	rownames(df) = names(maskList)
	colnames(df) = motifNameL
	for( anchorName in names(maskList) ){
		range = sapply(maskList[[anchorName]], maskToRange)
		df[anchorName, names(range)] = range
	}

	return(df)
}

compareMotifScoreAll=function(desPrefix, motifScore, main=NULL, mode="all"){

	stopifnot(mode %in% c("all","nontarget","target"))
	width=600
	height=600

	axisMin = min(motifScore)
	axisMax = max(motifScore)
	margin = (axisMax - axisMin) * 0.05
	axisMin = axisMin - margin
	axisMax = axisMax + margin
	axisLim = c(axisMin, axisMax)

	# split WT and Mutant score
	index.wt = which(sapply(rownames(motifScore), function(x) strsplit(x, ":")[[1]][2]=="WT"))
	score.wt = motifScore[index.wt,]
	score.mut = motifScore[-index.wt,]
	rownames(score.wt) = sub(":WT$","", rownames(score.wt))
	# resize WT score to match mutant score to compare
	score.wt = score.wt[sapply(rownames(score.mut), function(x) strsplit(x, ":")[[1]][1]),]

	if(mode=="all"){
		# All WT vs Mut scores
		df = data.frame(unlist(score.wt), unlist(score.mut))
	}else if(mode=="target"){
		# WT vs Mut for mutated target motifs only
		tmp = sapply(rownames(score.mut), function(x) strsplit(strsplit(x, ":")[[1]][2],",")[[1]])
		df = data.frame(matrix(0, length(unlist(tmp)), 2))
		index.row = unlist(sapply(1:length(tmp), function(x) rep(x, length(tmp[[x]]))))
		index.col = unlist(tmp)
		for( i in 1:length(index.row) ){
			df[i,1] = score.wt[index.row[i], index.col[i]]
			df[i,2] = score.mut[index.row[i], index.col[i]]
		}
	}else{
		# WT vs Mut for unmutated target motifs only
		motifNameL=colnames(motifScore)
		tmp = sapply(rownames(score.mut), function(x) setdiff(motifNameL, strsplit(strsplit(x, ":")[[1]][2],",")[[1]]))
		df = data.frame(matrix(0, length(unlist(tmp)), 2))
		index.row = unlist(sapply(1:length(tmp), function(x) rep(x, length(tmp[[x]]))))
		index.col = unlist(tmp)
		for( i in 1:length(index.row) ){
			df[i,1] = score.wt[index.row[i], index.col[i]]
			df[i,2] = score.mut[index.row[i], index.col[i]]
		}

	}
	png(sprintf("%s.png", desPrefix), width=width, height=height)
	par(las=1)
	tmp = drawDensity2D(df, xlim=axisLim, ylim=axisLim, transferFun=function(x) x^0.25)
	abline(a=0, b=1, col="purple")
	title(main)
	dev.off()
}




# draw pairwise densities of maximum motif scores comparing WT and single-motif mutant sequences
# Input:
#	desPrefix: destination prefix including path
#	motifScore: data.frame of motif score for the WT / mutated sequences
#		Row names are optional, but should be in the same order with the annotation data.frame
#		Column names are motif names
#	mode:
#		- all: all mutation case considered including additiona motif mutation
#			e.g. for PPARG, double mutation of PPARG and CEBPB also visualized
#		- exclusive: only exact single motif mutation cases are considred
# Output:
#	Scatter plot matrix
#		Row: single mutation
#		Column: motif to compare before vs after mutation
compareMotifScoreSingleMut=function(desPrefix, motifScore, main="NULL", mode="all"){

	stopifnot(mode %in% c("all","exclusive"))
	anchorNameL = sapply(rownames(motifScore), function(x) strsplit(x, ":")[[1]][1])
	motifNameL = colnames(motifScore)
	N.motif = length(motifNameL)

	width=200*(N.motif+1)
	height=200*(N.motif+1)

	axisLimL=list()	
	for( motifName in motifNameL ){
		axisMin = min(motifScore[,motifName])
		axisMax = max(motifScore[,motifName])
		margin = (axisMax - axisMin) * 0.05
		axisMin = axisMin - margin
		axisMax = axisMax + margin
		axisLim = c(axisMin, axisMax)
		axisLimL[[motifName]] = axisLim
	}

	## score data split
	# split WT and Mutant score
	index.wt = which(sapply(rownames(motifScore), function(x) strsplit(x, ":")[[1]][2]=="WT"))
	score.wt = motifScore[index.wt,]
	score.mut = motifScore[-index.wt,]
	rownames(score.wt) = sub(":WT$","", rownames(score.wt))
	# resize WT score to match mutant score to compare
	score.wt = score.wt[sapply(rownames(score.mut), function(x) strsplit(x, ":")[[1]][1]),]

	motifMutatedL = sapply(rownames(score.mut), function(x) strsplit(x, ":")[[1]][2])
	motifMutatedL.split = sapply(motifMutatedL, function(x) strsplit(x, ",")[[1]])

	png(sprintf("%s.png", desPrefix), width=width, height=height)
	par(mfrow=c(N.motif+1,N.motif+1), mar=c(1.5,2,1.5,1.5), oma=c(2,0,0,0), las=1)
	for( i in 0:N.motif ){
		# motifName="PRDM1"
		for( j in 0:N.motif ){
			if(i==0){
				if(j==0){
					plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
				}else{
					motifToCompare = motifNameL[j]
					plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
					text(x = 0.5, y = 0.2, sprintf("%s", motifToCompare), cex=2)
				}
			}else{
				motifMut = motifNameL[i]
				if(mode=="exclusive"){
					indexSelect = which(motifMutatedL==motifMut)
				}else{
					indexSelect = which(sapply(motifMutatedL.split, function(x) motifMut %in% x))
				}
				#anchorSelect = anchorNameL[indexSelect]
				if(j==0){
					plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
					text(x = 0.5, y = 0.5, sprintf("Mutated\n%s", motifMut), cex=2)
				}else{
					motifToCompare = motifNameL[j]
					#score.wt = motifScore[sprintf("%s:WT", anchorSelect), motifToCompare]
					#score.mut = motifScore[indexSelect, motifToCompare]
					s.wt = score.wt[indexSelect, motifToCompare]
					s.mut = score.mut[indexSelect, motifToCompare]
					tmp.df = data.frame(WT=s.wt, Mut=s.mut)
					axisLim = axisLimL[[motifToCompare]]
					tmp = drawDensity2D(tmp.df, transferFun=function(x) x^0.25, xlim=axisLim, ylim=axisLim)
					abline(a=0, b=1, col="purple")
				}
			}
		}
	}
	mtext(main, side=3, outer=TRUE, line=-4, cex=2)
	dev.off()
}


## Utility function to validate mutation results to see if they mutated location is within the mask range
# Input:
#	wt: WT sequence
#	mut: mutant sequence
#	mask: a list of mask for mutation
#	MotifMutated: comma-separted mutated motif name list
# Output:
#	Valid, Invalid, Same: if they are valid or not
validateMutant=function(wt, mut, mask, motifMutated){
	motifMutated = strsplit(motifMutated,",")[[1]]
	mask = do.call(rbind, mask[motifMutated])
	# location of mismtach
	ind = compareTwoSequence(wt,mut)
	if(length(ind)==0) return("Same")
	for( i in ind ){
		isValid=FALSE
		for( j in 1:nrow(mask) ){
			if( i >= mask[j,1] && i <= mask[j,2] ){
				isValid=TRUE
				break
			}
		}
		if(!isValid) return("Invalid")
	}
	return("Valid")
}





## Compare two sequence and print out the location of mismatch locations
compareTwoSequence=function(s1, s2){
	s1 = strsplit(s1, "")[[1]]
	s2 = strsplit(s2, "")[[1]]
	return(which(s1!=s2))

	# test code
	if(FALSE){
		i=3; j=11
		compareTwoSequence(df.mut[i,2], df.mut[j,2])
		maskL.groupByAnchor[[df.mut[i,1]]]
	}
}


## Perform motif scan and return maximum score for a list of motifs
# Input: 
#	fa: string vector or Biostring object of DNA sequences
#		each element must have name, which will become rownames for the result
#	motifL: a list of Homer motif files with name
#		e.g. motifL[["motifName"]] = "path/to/homer.motif"
# Output:
#	N x M data.frame of maximum motif score per sequence where
#		N: number of sequences
#		M: number of motifs
getMaxMotifScore=function(fa, motifL, parallel=4){
	if(FALSE){
		fa = df.mut$Sequence
		motifL = srcL.motif
		parallel=4
	}
	require(Biostrings)
	# convert to Biostring object if needed
	if(class(fa) == "character") fa = DNAStringSet(fa)
	df.maxScore = NULL
	for( motifName in names(motifL) ){
		# motifName = "PRDM1"
		write(sprintf("Scanning %s: %s", motifName, motifL[[motifName]]), stderr())
		motif = readMotif.Homer(motifL[[motifName]])
		tmp = getMotifScoreFa(fa, motif$PWM, parallel=parallel)
		maxScore = apply(cbind(tmp[[1]], tmp[[2]]), 1, max)
		if(is.null(df.maxScore)){
			df.maxScore = maxScore
		}else{
			df.maxScore = cbind(df.maxScore, maxScore)
		}
	}
	df.maxScore  = data.frame(df.maxScore)
	colnames(df.maxScore) = names(motifL)
	rownames(df.maxScore) = names(fa)
	return(df.maxScore)
}


## Generate combinatorial mutation for a give list of input sequence
# Input:
#	faL: DNA sequence list in Biostring format with unique anchor name
#	maskL: mutation mask list groupped by anchor name
#	comboBoolTableL: a list of preset combinatorial true table in various size
#	maxTryDiff: maximum number of mutation to avoid making same sequence after mutation
#	validate: Validate if all the mismatch after mutation fall into the mask region
# Output:
#	df.mut: a data.frame mutated sequence (plus original sequence)
generateMutatedSeqSet=function(faL, maskL, comboBoolTableL, mutationMode, maxTryDiff=1, validate=FALSE){

	# input validation
	stopifnot(setequal(names(faL), names(maskL)))

	# convert Biostring to character string
	data.seq = sapply(faL, toString)

	# Initializing output data.frame
	# df.mut will have two columns: Sequence & MotifMutated
	#	- AnchorName: anchor name
	#	- Sequence: mutated (or WT) sequence
	#	- MotifMutated: comma-separated list of motifs mutated
	df.mut = NULL
	# progress bar
	pbar = txtProgressBar(min = 1, max = length(data.seq), initial = 1, style=3) 
	for( i in 1:length(data.seq) ){
		# i=1
		anchorName = names(data.seq)[i]
		s = data.seq[i]

		mask = maskL.groupByAnchor[[anchorName]]
		motiflist = names(mask)
		n = length(motiflist)

		# generate all combinatorial mutations including WT
		comboTable = comboBoolTableL[[n]]
		colnames(comboTable) = motiflist
		s.mut = generateMutatedSeqCombo(s, mask, comboTable, mutationMode, maxTryDiff=maxTryDiff)

		if(validate){
			# check if the mismatch locations are within the mask range
			for( j in 2:nrow(s.mut) ){
				tmp = validateMutant(s.mut[1,1], s.mut[j,1], mask, s.mut[j,2])
				if(tmp != "Valid") write(tmp, stderr())
			}
		}

		df.mut = rbind(df.mut, data.frame(AnchorName=anchorName, s.mut))
		setTxtProgressBar(pbar, i)
	}
	close(pbar)
	# assign unique name for each mutated sequence by combining anchor name and mutated motifs
	rownames(df.mut) = sprintf("%s:%s", df.mut$AnchorName, df.mut$MotifMutated)

	return(df.mut)
}




# if(FALSE){
# 	# debug code
# 	df.mut = generateMutatedSeqSet(data.fa, maskL.groupByAnchor, comboBoolTableL, mutationMode, maxTryDiff=3, validate=TRUE)

# 	anchoName="chr20-1"
# 	s = toString(data.fa[1])
# 	mask = maskL.groupByAnchor[[anchoName]]
# 	comboTable = comboBoolTableL[[length(mask)]]
# 	colnames(comboTable) = names(mask)
# 	mode="shuffle"
# 	maxTryDiff=3
# }


## For a give sequence and mask list, generate a set of sequence by combinatorial mutations
# Input:
#	s: DNA sequence in R character string
#	mask: a list of masks
#		e.g. mask[["TF1"]] = <N x 2 data.frame> of masks for start and end 1-base coordinate)
#	comboTable: N x M Boolean table of the combinations of motifs to mutate
#		N: number of combinations
#		M: number of motifs
#		Column names must match with motif list in mask
# Output: a list containing two vectors
#	- Sequence: a vector of mutated sequence potentially including the original one
#	- MotifMutated: a vector of comma-separated motif names mutated
generateMutatedSeqCombo=function(s, mask, comboTable, mode="shuffle", maxTryDiff=1){
	# comboTable=comboBoolTableL[[n]]
	if(FALSE){
		mode=mutationMode
	}
	motiflist=names(mask)
	stopifnot(length(motiflist) == ncol(comboTable))
	stopifnot(all(motiflist == colnames(comboTable)))

	# a vector of mutated sequence (including WT)
	seqL = NULL
	# a comma-separated mutated motif list
	MotifMutatedL = NULL
	# mask string
	maskStrL = NULL
	# mutation
	for( i in 1:nrow(comboTable) ){
		motif.select = names(which(comboTable[i,]))
		if(length(motif.select)==0){
			s.mut = s
			MotifMutated = "WT"
			maskStr = "NA"
		}else{
			df.mask = do.call(rbind, mask[motif.select])
			df.mask = df.mask[order(df.mask[,1], df.mask[,2]),]
			s.mut = mutateSeq(s, df.mask, mode, maxTryDiff)
			MotifMutated = paste(motif.select, collapse=",")
			maskStr = maskToRange(df.mask)
		}
		seqL = c(seqL, s.mut)
		MotifMutatedL = c(MotifMutatedL, MotifMutated)
		maskStrL = c(maskStrL, maskStr)
	}

	return(data.frame(Sequence = seqL, MotifMutated = MotifMutatedL, Mask=maskStrL))
}



## make random DNA sequence of lenght n
# input: length of sequence
# output:
#	- string if asVecor=FALSE
#	- character vector if asVector=TRUE
makeRandomDNA=function(n, asVector=FALSE){
	s = c("A","C","G","T")[sample(1:4, n, replace=TRUE)]
	if(! asVector) s = paste(s, collapse="")
	return(s)
}


## shuffle given DNA sequence
# input:
#	- string if asVecor=FALSE
#	- character vector if asVector=TRUE
# output:
#	- string if asVecor=FALSE
#	- character vector if asVector=TRUE
shuffleDNA=function(s, asVector=FALSE){
	if(! asVector){
		stopifnot(length(s)==1)
		s = paste(strsplit(s, "")[[1]][sample(1:nchar(s))], collapse="")
	}else{
		s = s[sample(1:length(s))]
	}
	return(s)
}


## Mutate an input sequence within the given mask regions
# Input:
#	- s: character string of DNA
#	- mask: n x 2 data.frame of start/end ranges to mutate
#	- mode: mutation mode. {random, shuffle}
#		random: random sequence
#		shuffle: bases are shuffled in each mask region
# Output:
#	- Mutated character string of DNA
mutateSeq = function(s, mask, mode, maxTryDiff=1){

	if(FALSE){
		# test code
		mutateSeq("ACGTACGTACGT", data.frame(c(1,7),c(4,8)), mode="random")
		mutateSeq("ACGTACGTACGT", data.frame(c(1,7),c(4,8)), mode="shuffle")
	}
	stopifnot(mode %in% c("random", "shuffle"))

	# string to vector of individual character
	s = strsplit(s, "")[[1]]
	# mutate
	if(mode == "random"){
		for( i in 1:nrow(mask) ){
			tryCnt=1
			while(TRUE){
				s.mut = makeRandomDNA(mask[i,2] - mask[i,1] + 1, asVector=TRUE)
				if(any(s.mut != s[mask[i,1]:mask[i,2]]) || tryCnt>=maxTryDiff) break
				tryCnt = tryCnt + 1
			}
			if(tryCnt > 1) write(sprintf("Tried %d", tryCnt), stderr())
			s[mask[i,1]:mask[i,2]] = s.mut
		}
	}else{
		for( i in 1:nrow(mask) ){
			tryCnt=1
			while(TRUE){
				s.mut = shuffleDNA(s[mask[i,1]:mask[i,2]], asVector=TRUE)
				if(any(s.mut != s[mask[i,1]:mask[i,2]]) || tryCnt>=maxTryDiff) break
				tryCnt = tryCnt + 1
			}
			if(tryCnt > 1) write(sprintf("Tried %d", tryCnt), stderr())
			s[mask[i,1]:mask[i,2]] = s.mut
		}
	}
	# convert character vector to string
	s = paste(s, collapse="")
	return(s)
}



## Group given list of masks for multiple TFs by anchor names
# Input:
#	maskL: a list of mask for multiple TFs groupped by TF name
# Oupput:
#	maskL.byAnchor: a list of mask for multiple TFs groupped by anchor name
#		i.e. maskL.byAnchor[[anchorName]] = list( TF1=<data.frame of mask>, TF2=<data.frame of mask> ... )
groupMaskByAnchor=function(maskL){
	anchorNameL = unique(unlist(sapply(maskL, function(x) x$Name)))
	maskL.byAnchor = list()
	for( anchorName in anchorNameL ) maskL.byAnchor[[anchorName]] = list()
	for( factorName in names(maskL) ){
		# factorName="PRDM1"
		mask = maskL[[factorName]] %>% split(f = as.factor(.$Name))
		mask = lapply(mask, function(x) x[,-1])
		for( anchorName in names(mask) ){
			if(is.null(maskL.byAnchor[[anchorName]])) maskL.byAnchor[[anchorName]]=list()
			maskL.byAnchor[[anchorName]][[factorName]] = mask[[anchorName]]
		}
	}

	return(maskL.byAnchor)

	# debug code
	if(FALSE){
		head(df.motifStatus.filtered)
		maskL.byAnchor[["chr20-1"]]
		table(sapply(maskL.byAnchor, length))
		table(apply(df.motifStatus.filtered != "None", 1, sum))
	}
}


## Convert Boolean table to a vector of comma-separted column names of TRUE
boolTableToVector=function(tf){
	v = apply(tf, 1, function(x) paste(colnames(tf)[unlist(x)], collapse=","))
	return(v)
}

## Convert motifScanList object into a data.frame of the number of anchors having each category of motif instances
# Input: a list of motif scan results
# Output: a data.frame of N x 5 where N is the number of factors
#	Columns are:
#		TF: name of transcription factor (or motif)
#		Stringent / Marginal: the number of anchors having Stringent / Marginal motif instances
#		None: the number of anchors without even marginal instance
#		StrOrMarginal: the number of anchors having either Stringent or Marginal instance
motifScanListToCntTable=function(motifScanL){
	df.cnt = data.frame(TF=names(motifScanL), Stringent=0, Marginal=0, None=0)
	rownames(df.cnt) = names(motifScanL)
	tmp.cnt = sapply(motifScanL, function(x) table(x$status))
	for( name in names(tmp.cnt) ){
		for( type in names(tmp.cnt[[name]]) ) df.cnt[name, type] = tmp.cnt[[name]][type]
	}
	df.cnt$StrOrMarginal = df.cnt$Stringent + df.cnt$Marginal

	return(df.cnt)
}


## Convert a list of motif scan results to status table 
# Input:
#	- motifScanL:
#		a list sof motif scan result
# Output:
#	- data.frame with N x M
#		N: number of anchors
#		M: number of motifs
#		each elemenet is one of Stringent / Marginal / None
motifScanListToStatusTable=function(motifScanL){
	anchorNameL = names(motifScanL[[1]]$status)
	motifNameL = names(motifScanL)
	N.anchor = length(anchorNameL)
	N.motif = length(motifScanL)

	# initialize output data.frame
	df.status = data.frame(matrix("", N.anchor, N.motif))
	rownames(df.status) = anchorNameL
	colnames(df.status) = motifNameL

	for( motifName in motifNameL ){
		stopifnot(all(rownames(df.status) == names(motifScanL[[motifName]]$status)))
		df.status[,motifName] = motifScanL[[motifName]]$status
	}

	return(df.status)
}


## Combine motif masks and TF overlap table to generate final motif combinations to mutate
# Some scan window overlapping TF1 may not contain stringent or marginal motif instances of TF1
# In this case, TF1 is ignored in this region and not mutated for TF1
# Input:
#	- motifScanL:
#		a list sof motif scan result
#	- df.cobind: Boolean data.frame of TF overlap
#		Row: motif scan window, Rownames are the scan window name (4th column of the scan window bed file)
#		Col: TF or motif name
# Output:
#	- a list of filtered motif scan results
# Note:
#	"Stringent" or "Marginal" is set to "None" if no factor overlap at the anchor
#	Thus, length of the status is still the same with the original
filterMotifByTF=function(motifScanL, df.cobind){
	motifScanL.filtered = list()
	for( motifName in names(motifScanL) ){
		# scan results for a motif
		scan.select = motifScanL[[motifName]]
		# anchor names that has TF overlap
		anchorName.withTF = rownames(df.cobind)[df.cobind[,motifName]]

		# intersection with TF overlap
		scan.select$instance = scan.select$instance[scan.select$instance$Name %in% anchorName.withTF,]
		scan.select$status[! names(scan.select$status) %in% anchorName.withTF] = "None"
		motifScanL.filtered[[motifName]] = scan.select
	}
	return(motifScanL.filtered)
} 



## Make a Boolean matrix of TF combination
# Input: A vector of TF names
# Output: Boolean data frame of all possible combinations of TF selection
#	Column: TF names
#	Row: TRUE/FALSE of TF selection, 2^length(nameL)
makeComboBoolTable = function(n, nameL=NULL){
	library(data.table)
	boolTable = do.call(CJ, replicate(n, 0:1, FALSE)) > 0
	if(!is.null(nameL)) colnames(boolTable) = nameL
	return(boolTable)
}


## For a given list of motif masks, check their pairwise overlap between two TFs
# Input:
#	- list of motif mask generated from instanceToMask
# Output:
#	- A pairwise plot matrkx of
#		Upper half: histogram of overlapping length between two TF motif masks
#		Lower half: number of motif scan windows where at least 1bp overlap happens between to TFs masks
checkMaskOverlap=function(maskL, desPrefix){
	n = length(maskL)
	nameL = names(maskL)

	# plotting layout
	mat = matrix(0, nrow=n, ncol=n)
	ind = 1
	for( i in 1:n ){
		for( j in i:n ){
			mat[i,j] = ind
			ind = ind + 1
			if(i!=j){
				mat[j,i] = ind
				ind = ind + 1
			}
		}
	}

	# convert to 0-base coordinate to use bedtools
	for( i in 1:length(maskL) ) maskL[[i]][,2] = maskL[[i]][,2] - 1

	# check overlap and pairwise draw plot matrix
	pdf(sprintf("%s.pdf", desPrefix), width=1.5*n, height=1.5*n)
	layout(mat)
	write(sprintf("Checking overlap between"), stderr())
	for( i in 1:n ){
		# i=1
		for( j in i:n ){
			# j=2

			if(i==j){
				# diagonal element: TF name
				name = nameL[i]
				par(mar = c(0,0,0,0))
				plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
				text(x = 0.5, y = 0.5, name, cex=2)
			}else{
				name1 = nameL[i]
				name2 = nameL[j]
				mask1 = maskL[[i]]
				mask2 = maskL[[j]]
				write(sprintf("  - %s vs %s", name1,  name2), stderr())

				# check overlap
				tmp = bedTools.2in("intersectBed", maskL[[name1]], maskL[[name2]], "-wo", direct=FALSE, debug=FALSE)
				n.overlap = length(unique(tmp[,1]))

				# upper triangle: histogram of overlapping bp
				par(mar=c(2,2,2,0.5), las=1)
				hist(tmp[,7], plot=TRUE, main=NULL)
				# lower triangle: number of scanwindows that has overlapping mask
				par(mar = c(0,0,0,0))
				plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
				text(x = 0.5, y = 0.5, n.overlap, cex=2)
			}
		}
	}
	dev.off()
	system(sprintf("convert -density 300 %s.pdf %s.png", desPrefix, desPrefix))
}






## convert motif instance data.frame to mask in bed format by merging overlapping motif instances
# mask range is 1-based cooridnate
instanceToMask=function(instance){
	stopifnot(all(c("Name","Offset","Score","Direc","Seq") %in% names(instance)))
	motifLen = nchar(instance[1,"Seq"])
	# convert instance to bed format
	bed.instance = data.frame(Name = instance$Name,
				Start = instance$Offset - 1,
				End = instance$Offset + motifLen - 1,
				None = "NULL",
				Score = instance$Score,
				Direc=instance$Direc)
	# order by coordinate
	bed.instance = bed.instance[ order(bed.instance$Name, bed.instance$Start, bed.instance$End), ]
	# merge to make mask
	mask = bedTools.1in("mergeBed", bed.instance, direct=FALSE, debug=FALSE)
	mask[,2] = mask[,2] + 1
	colnames(mask) = c("Name","Start","End")
	return(mask)
}

## apply motif score threshold for a given score matrix
#	i.e. if stringent match exist (>cutoff1), select all of them
# 	if no stringent match, find single maximum value position above cutoff 2
# Return:
#	list of
#	- Boolean matrix of the same size denoting selection of the motif instances
#	- String vector of motif scan status of each scan window, Stringent / Marginal / None
#		in the same order with the input scan window names
applyCutoffMotif=function(score, cutoff1, cutoff2){

	# bed window X location Boolean matrix of motif instance
	# initial stringent filtering
	scanBoolMatrix = score > cutoff1

	tf.strong = apply(scanBoolMatrix, 1, sum) > 0
	tf.marginal = !tf.strong

	# marginal filtering
	# apply marginal cutoff for regions without stringent match
	# by finding single maximum score position > cutoff2
	findMarginalMax = function(scoreVec, cutoff2){
		scanBoolVec = rep(FALSE, length(scoreVec))
		ind.max = which.max(scoreVec)
		if(scoreVec[ind.max] > cutoff2) scanBoolVec[ind.max] = TRUE
		return(scanBoolVec)
	}

	scanBoolMatrix.marginal = t(apply(score[tf.marginal,],1, function(x) findMarginalMax(x, cutoff2)))
	scanBoolMatrix[tf.marginal] = scanBoolMatrix.marginal

	# status vector: stringent / marginal / none
	status=rep("None", nrow(score))
	names(status) = rownames(score)
	status[tf.strong] = "Stringent"
	status[tf.marginal][apply(scanBoolMatrix.marginal, 1, any)] = "Marginal"

	return(list(scanBoolMatrix=scanBoolMatrix, status=status))
}

## Motif scan for MPRA sequence design
# Input:
#	- fa: R string or fasta object Biostring
#	- pwm: PWM object from PWMEnrich package
#	- cutoff1: stringent cutoff
#	- cutoff2: marginal cutoff
#	- parallel: number of cores to use
# Output:
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
scanMotifMPRA=function(fa, pwm, cutoff1, cutoff2, parallel=NULL, verbose=FALSE){
	
	if(class(fa) == "character") fa = DNAStringSet(fa)

	len.motif = length(pwm)

	if(verbose) write("Scanning ...", stderr())
	scoreL = getMotifScoreFa(fa, pwm, parallel)
	scoreAll = cbind(scoreL[[1]], scoreL[[2]])
	maxScore = apply(scoreAll, 1, max)

	## apply stringent / marginal threshold
	## i.e. if stringent match exist (>cutoff1), select them
	# 	if no stringent match, find single maximum value position above cutoff 2
	if(verbose) write("Applying threshold ...", stderr())
	result.scan = applyCutoffMotif(scoreAll, cutoff1, cutoff2)
	tf = result.scan[["scanBoolMatrix"]]
	status = result.scan[["status"]]

	# split forward / reverse strand portions
	col.half = ncol(tf)/2
	tf.f = tf[,1:col.half]
	tf.r = tf[,(col.half+1):ncol(tf)]

	# total number of hits
	N.hit = sum(tf)

	# select window names with hits only
	ind = apply( tf, 1, any )
	nameL = rownames(tf)[ind]

	# Initialize empty data.frame for motif instances found
	if(verbose) write("Making result data.frame ...", stderr())
	result = as.data.frame(matrix(0, ncol=5, nrow=N.hit), stringsAsFactor=FALSE)
	colnames(result) = c("Name","Offset","Score","Direc","Seq")
	index=1

	# iterate over motif scan window
	# note: subseq function from Biostring is too slow, hence sequences are processed in character string using a native R function
	for( i in 1:length(nameL) ){
	# i=1
		name = nameL[i]
		w = width(fa[name])
		
		fa.str = as.character(fa[name])
		fa.rc = as.character(reverseComplement(fa[name]))
		# motif instance locations in the current forward scan window
		offsetL = which(tf.f[name,])
		for( offset in offsetL ){
			s = offset
			e = offset + len.motif - 1
			seq = substr(fa.str, s, e)
			score = scoreL[["forward"]][name, offset]
			result[index,] = data.frame(name, offset, score, "+", seq, stringsAsFactors=FALSE)
			index = index + 1
		}
		
		# motif instance locations in the current forward scan window
		offsetL = which(tf.r[name,])
		for( offset in offsetL ){
			e = offset + len.motif - 1
			s.rc = w - e + 1
			e.rc = s.rc + len.motif - 1
			seq = substr(fa.rc, s.rc, e.rc)
			score = scoreL[["reverse"]][name,offset]
			result[index,] = data.frame(name, offset, score, "-", seq, stringsAsFactors=FALSE)
			index = index + 1
		}
	}

	return(list(instance=result, status=status, maxScore=maxScore))
}



## Compare mismatches between two sequences and print results in BLAST-like format
# Written by Christopher Ahn; Revised by Hee Woong Lim
# Input:
#	- sequence1 (e.g. "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTTAAAAAAAAAACCCCCCCCCC")
#	- sequence2 (e.g. "AAAAAAAAAACCCCCACCCCGGGGGAAGGGTTTATTTTTTAAAAAAAAAACCCCCCCCCCGGGGTGGGGGTTTTTTTTTTAAAAAAAAAACCCCGCCCCC")
#	- character vector of masked regions; e.g. maskVector = c(mask1="10-20,1-3,35-40,55-70", mask2="10-20,30-40,90-120")
#		Each range is comma separated in each mask region.
#		Masked region coordinates should be relative to the input sequence length, and region can be larger than the sequence length.
#		Masked regions do not have to be in order, but regions coming later will overwrite any overlapping previous regions.
#	- width of each line of sequence output. If width = 50, 50 characters of each sequence will be printed in each line. Default = 100
compareTwoStrings = function(seq1, seq2, mask, width=100) {
	# Check if input sequences are longer than 0, and if both sequences have the same length
	if (nchar(seq1) < 1 || nchar(seq2) < 1) stop("Sequence is not valid.")
	if (nchar(seq1) != nchar(seq2)) stop("Sequences must have the same length.")
	# empty mask name handling
	if (is.null(names(mask)) || any(is.na(names(mask)))) names(mask) = sprintf("Mask%d", 1:length(mask))

	# Check if width is valid, and truncate if longer than sequence length
	if (!is.numeric(width)) {
		stop("Width parameter is not numeric.")
	} else if ( width%%1!=0 || width < 1 ) {
		stop("Width parameter needs to be an integer larger than zero.")
	} else if ( width > nchar(seq1) ) width = nchar(seq1)

	# Initialize mask sequence vector
	namedMaskSeq=NULL
	# Generate mask string
	for (name in names(mask)) {

		## Initialize mask sequence with empty spaces
		maskSeq = rep(" ", nchar(seq1))

		## Parse the mask regions
		ranges = strsplit(strsplit(mask[name], ",")[[1]], "-")

		## Loop through each range and update mask sequence
		for (range in ranges) {
			start = as.numeric(range[1])
			end = as.numeric(range[2])

			## Check if ranges are valid
			if ( is.na(start) || is.na(end) ) {
				stop(sprintf('Incorrect range: "start = %d" and "end = %d". Make sure there are no negative numbers."', start, end))
			} else if ( start > end ) {
				stop(sprintf('Incorrect range: "start = %d" and "end = %d". Make sure start >= end."', start, end))
			}

			## mask the regions with "|" as inclusive boundaries and "-" for all other masked regions
			maskSeq[start:end] = "-"
			maskSeq[start] = "|"
			maskSeq[end] = "|"
		}

		namedMaskSeq = c(namedMaskSeq, paste(maskSeq,collapse=""))
	}
	# assign names to each mask sequence and truncate name if longer than 9 chars
	names(namedMaskSeq) = substr(names(mask), 1, 9)

	# Generate mismatch string
	mismatches = paste(ifelse(strsplit(seq1, "")[[1]] == strsplit(seq2, "")[[1]], " ", "*"), collapse = "")

	# Print
	for (i in 1:ceiling((nchar(seq1))/width) ) {

		# get start and end coordinates
		currentStart = 1 + (width * (i - 1))
		currentEnd = currentStart + width - 1

		# if at the last window, handle cases where output line is shorter than previous lines
		if ( currentEnd > nchar(seq1) ) {
			currentEnd = nchar(seq1)
			tmpWidth = currentEnd - width

			if ( tmpWidth < 1 || tmpWidth == width ) {
				break
			} else {
				width = currentEnd - currentStart
			}
		}

		# Select segment to print
		seq1_window = substr(seq1, start = currentStart, stop = currentEnd)
		seq2_window = substr(seq2, start = currentStart, stop = currentEnd)
		mismatch_window = substr(mismatches, start = currentStart, stop = currentEnd)

		# print sequence and mismatches; left-align first column
		write(sprintf("%-9s%5d %*s%5d", "Sequence1", currentStart, width, seq1_window, currentEnd), stdout())
		write(sprintf("%-9s%5s %s", "Mismatch", " ", mismatch_window), stdout())
		write(sprintf("%-9s%5d %*s%5d", "Sequence2", currentStart, width, seq2_window, currentEnd), stdout())

		# truncate masked regions by window size and print
		for (name in names(namedMaskSeq)) {
			mask_window = substr(namedMaskSeq[name], start = currentStart, stop = currentEnd)
			write(sprintf("%-9s%5s %s", name, " ", mask_window), stdout())
		}
		
		## add empty line after current row block
		write("", stdout())
	}
}

compareTwoStrings(df.mut[1,"Sequence"], df.mut[6,"Sequence"], df.mut[6,"Mask"], width=100)
