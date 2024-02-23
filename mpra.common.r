

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
## Note: this doesn't handle single range case; need update
rangeToMask=function(range){
	if(FALSE){
		range = df.fimo.filtered[[1]][,1][1]
		range = df.fimo.filtered[[1]][,1][4]
	}
	mask = do.call(rbind, strsplit(strsplit(range, ",")[[1]],"-"))
	mask = as.data.frame(t(apply(mask, 1, as.numeric)))
	#mask = sapply(data.frame(mask), as.numeric)
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
	if(FALSE){
		compareMotifScoreAll("motifScoreCompareAll.Target", maxMotifScore.mut, main="Motif Score WT vs Mut: Targeted Motifs", mode="target")
		desPrefix="motifScoreCompareAll.Target"
		motifScore=maxMotifScore.mut
		main="Motif Score WT vs Mut: Targeted Motifs"
		mode="target"
	}
	stopifnot(mode %in% c("all","nontarget","target","scramble"))
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
	index.scram = which(sapply(rownames(motifScore), function(x) strsplit(x, ":")[[1]][2]=="Scramble"))
	score.wt = motifScore[index.wt,]
	score.scram = motifScore[index.scram,]
	score.mut = motifScore[-c(index.wt, index.scram),]
	rownames(score.wt) = sub(":WT$","", rownames(score.wt))
	rownames(score.scram) = sub(":Scramble$","", rownames(score.scram))
	# resize WT score to match mutant score to compare
	score.wt.tiled = score.wt[sapply(rownames(score.mut), function(x) strsplit(x, ":")[[1]][1]),]

	if(mode=="all"){
		# All WT vs Mut scores
		df = data.frame(unlist(score.wt.tiled), unlist(score.mut))
	}else if(mode=="target"){
		# WT vs Mut for mutated target motifs only
		tmp = sapply(rownames(score.mut), function(x) strsplit(strsplit(x, ":")[[1]][2],",")[[1]])
		df = data.frame(matrix(0, length(unlist(tmp)), 2))
		index.row = unlist(sapply(1:length(tmp), function(x) rep(x, length(tmp[[x]]))))
		index.col = unlist(tmp)
		for( i in 1:length(index.row) ){
			df[i,1] = score.wt.tiled[index.row[i], index.col[i]]
			df[i,2] = score.mut[index.row[i], index.col[i]]
		}
	}else if(mode=="nontarget"){
		# WT vs Mut for unmutated target motifs only
		motifNameL=colnames(motifScore)
		tmp = sapply(rownames(score.mut), function(x) setdiff(motifNameL, strsplit(strsplit(x, ":")[[1]][2],",")[[1]]))
		df = data.frame(matrix(0, length(unlist(tmp)), 2))
		index.row = unlist(sapply(1:length(tmp), function(x) rep(x, length(tmp[[x]]))))
		index.col = unlist(tmp)
		for( i in 1:length(index.row) ){
			df[i,1] = score.wt.tiled[index.row[i], index.col[i]]
			df[i,2] = score.mut[index.row[i], index.col[i]]
		}
	}else{
		# WT vs Scramble
		score.scram = score.scram[rownames(score.wt),]
		df = data.frame(unlist(score.wt), unlist(score.scram))
	}
	png(sprintf("%s.png", desPrefix), width=width, height=height)
	par(las=1)
	tmp = drawDensity2D(df, xlim=axisLim, ylim=axisLim, transferFun=function(x) x^0.25)
	abline(a=0, b=1, col="purple")
	title(main, xlab="WT", ylab="Mutant (or Scramble)")
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


compareMotifScoreScramble=function(desPrefix, motifScore, main="NULL"){
	if(FALSE){
		desPrefix = "motifScoreComparePair.Scramble"
		motifScore=maxMotifScore.mut
		main="Motif Score WT vs Scramble"

	}
	anchorNameL = sapply(rownames(motifScore), function(x) strsplit(x, ":")[[1]][1])
	motifNameL = colnames(motifScore)
	N.motif = length(motifNameL)

	width=200*(N.motif+1)
	height=300

	motifScore = motifScore[grep(":WT$|:Scramble$", rownames(motifScore)),]
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
	score.scram = motifScore[-index.wt,]
	rownames(score.wt) = sub(":WT$","", rownames(score.wt))
	rownames(score.scram) = sub(":Scramble$","", rownames(score.scram))
	score.scram = score.scram[rownames(score.wt),]

	#motifMutatedL = sapply(rownames(score.mut), function(x) strsplit(x, ":")[[1]][2])
	#motifMutatedL.split = sapply(motifMutatedL, function(x) strsplit(x, ",")[[1]])

	png(sprintf("%s.png", desPrefix), width=width, height=height)
	par(mfrow=c(1,N.motif+1), mar=c(3,2,2,1.5), oma=c(3,0,4,0), las=1)
	for( i in 0:N.motif ){
		# motifName="PRDM1"
		if(i==0){
				plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
				text(x = 0.5, y = 0.5, "Scramble", cex=2)
		}else{
			motifToCompare = motifNameL[i]
			s.wt = score.wt[, motifToCompare]
			s.scram = score.scram[, motifToCompare]
			tmp.df = data.frame(WT=s.wt, Scramble=s.scram)
			axisLim = axisLimL[[motifToCompare]]
			tmp = drawDensity2D(tmp.df, transferFun=function(x) x^0.25, xlim=axisLim, ylim=axisLim)
			title(motifToCompare, xlab="WT")
			abline(a=0, b=1, col="purple")
		}
	}
	mtext(main, side=3, outer=TRUE, line=1, cex=2)
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


# Get the maximum motif score for a given list of motifs after FIMO
#
# Input:
#	fa: character vector containing DNA sequence with names of unique id
#	meme: MEME motif file
#	motifNameL: a list of motif names (alt_name) in MEME file
#	pv: p-value threshold for motif scan
#		IF not motif instance below this threshold, zero score is assigned to the sequence
#	des: destination file to save the fimo results
# Output:
#	N x M data.frame containing the maximum motif score per sequence per motif
#		N: number of sequences
#		M: number of motifs
getMaxMemeScore=function(fa, meme, motifNameL, pv=0.01, des="crs.fimo.txt"){
	if(FALSE){
		fa = df.mut$Sequence
		meme = src.motif
	}
	assertFileExist(meme)
	
	tmp.fa = tempfile(fileext=".fa")
	#fimo = "fimo.mut.txt"
	if(any(is.null(names(fa)))) names(fa) = sprintf("Seq%d", 1:length(fa))
	writeXStringSet(DNAStringSet(fa), tmp.fa)
	cmd=sprintf("runMemeFimo.Peak.sh -g hg38 -p %f %s %s > %s", pv, meme, tmp.fa, des)
	system(cmd)
	unlink(tmp.fa)

	result = read.delim(des, header=TRUE, stringsAsFactors=FALSE, comment.char = "#", blank.lines.skip = TRUE)


	df.maxScore = NULL
	for( motifName in motifNameL ){
		tmp = result[result$Alt_Name == motifName,]
		# maximum score per sequence
		tmp.max = aggregate(Score ~ SeqName, data=tmp, max)
		rownames(tmp.max) = tmp.max$SeqName

		# initialization; default value in case of no instance < pv
		maxScore = rep(0, length(fa))
		# filling max motif score
		names(maxScore) = names(fa)
		maxScore[rownames(tmp.max)] = tmp.max$Score

		if(is.null(df.maxScore)){
			df.maxScore = maxScore
		}else{
			df.maxScore = cbind(df.maxScore, maxScore)
		}
	}
	df.maxScore  = data.frame(df.maxScore)
	colnames(df.maxScore) = motifNameL
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

		# generate scrambled sequence
		seq.scram = shuffleDNA(s)
		s.mut = rbind(s.mut, data.frame(Sequence = seq.scram, MotifMutated = "Scramble", Mask=sprintf("%d-%d", 1, nchar(s))))

		# merge to the result data.frame
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
			# merge overlapping masks -> not used
			#df.mask = mergeMasks(df.mask)
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
# e.g. 
#	for a row, [ TRUE, FALSE, TRUE ], and column names [ C1, C2, C3 ]
#	the resulting combination is "C1,C3"
#	This conversion is done for all rows
#	For rows with FALSE only, "None" is assigned as a name
boolTableToVector=function(tf){
	v = apply(tf, 1, function(x) paste(colnames(tf)[unlist(x)], collapse=","))
	v[v==""] = "None"
	return(v)
}

## Count frequency comma-separted column names whose elements are TRUE
# Input:
#	Data.frame of TRUE/FALSE
# Output:
#	A vector of frequencies whose names are all possible combinations of column names
#	** Combinations of zero occurances are also reported i.e. 2^n rows output for n-column input
getComboFrequency=function(tf){
	# Generate a vector of comma-separted column names and count frequency
	comboNames = boolTableToVector(tf)
	comboFreq = table(comboNames)

	# Prepare all possible combinations of column names
	tmp.tf = makeComboBoolTable(ncol(tf))
	tmp = apply(tmp.tf, 1, function(x) paste(colnames(tf)[unlist(x)], collapse=","))
	tmp[tmp==""] = "None"

	# Fill the frequency values
	result = rep(0, length(tmp))
	names(result) = tmp
	result[names(comboFreq)] = comboFreq

	return(result)
}

## Convert motifScanList object into a data.frame of the number of anchors having each category of motif instances
# Input: a list of motif scan results
# Output: a data.frame of N x 5 where N is the number of factors
#	Columns are:
#		TF: name of transcription factor (or motif)
#		Stringent / Marginal: the number of anchors having Stringent / Marginal motif instances
#		None: the number of anchors without even marginal instance
#		StrOrMarginal: the number of anchors having either Stringent or Marginal instance
motifScanListToCntTable=function(mscan){
	df.cnt = data.frame(TF=names(mscan), Stringent=0, Marginal=0, None=0)
	rownames(df.cnt) = names(mscan)
	for( motifName in names(mscan) ){
		tmp = table(mscan[[motifName]]$status)
		for( status in names(tmp) ) df.cnt[motifName, status] = tmp[status]
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

	# identify missing motif and drop anchors
	status.filtered = motifScanListToStatusTable(motifScanL.filtered)
	stopifnot(all(rownames(df.cobind) == rownames(status.filtered)))
	stopifnot(all(colnames(df.cobind) == colnames(status.filtered)))
	index.missingMotif = apply(! df.cobind == (status.filtered != "None"), 1, any)
	anchorNames = rownames(df.cobind)[!index.missingMotif]

	df.cobind = df.cobind[anchorNames,]
	for(motifName in names(motifScanL.filtered)){
		# motifName = names(motifScanL.filtered)[1]
		tmp = motifScanL.filtered[[motifName]]
		tmp[["instance"]] = tmp[["instance"]][tmp[["instance"]]$Name %in% anchorNames,]
		tmp[["status"]] = tmp[["status"]][anchorNames]
		tmp[["maxScore"]] = tmp[["maxScore"]][anchorNames]
		motifScanL.filtered[[motifName]] = tmp
	}

	resultL=list( motifScan = motifScanL.filtered, cobind = df.cobind )
	return(resultL)
} 


## Intersect cobind and motif scan results
# - Drop motif instances without TF cobinding for each TF
# - Drop anchors without a desired motif instance but with TF cobinding only
# Input:
#	- cobind: N x M data.frame where N is the number of anchors and M is the number of TFs (or motifs)
#	- motifScan: a list of data.frames for motif instances and status, each of which are the same sied data.frames with cobind
#		motifScan[["instance"]]: N x M data.frame containing all motif instance in comma-separated range format
#		motifScan[["status"]]: N x M data.frame containing motif instance status, {Stringent, Marginal, None}
# Ouptut:
#	A list containing
#	result[["cobind"]]: filtered cobind data.frame
#	result[["motifScan"]]: filtered motifScan list
intersectCobindMotif=function(cobind, motifScan){
	if(FALSE){
		cobind = df.cobind
		motifScan = df.fimo
	}

	# mask out motif instance without corresponding TF cobinding
	for( motifName in colnames(cobind) ){
		# anchor names that has TF overlap
		anchorName.noCobind = rownames(cobind)[!cobind[,motifName]]

		# intersection with TF overlap
		motifScan[["instance"]][anchorName.noCobind, motifName] = ""
		motifScan[["status"]][anchorName.noCobind, motifName] = "None"
	}

	# drop anchors desired motif instances with a corresponding TF binding
	index.missingMotif = which(apply( cobind & (motifScan[["status"]] == "None"), 1, any))
	cobind = cobind[-index.missingMotif, ]
	motifScan[["instance"]] = motifScan[["instance"]][-index.missingMotif, ]
	motifScan[["status"]] = motifScan[["status"]][-index.missingMotif, ]
	
	result=list( cobind = cobind, motifScan = motifScan )
	return(result)
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
				if(n.overlap > 0){
					par(mar=c(2,2,2,0.5), las=1)
					hist(tmp[,7], plot=TRUE, main=NULL)
				}else{
					# handling non overlap
					par(mar = c(0,0,0,0))
					plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
					text(x = 0.5, y = 0.5, n.overlap, cex=2)
				}
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

## Conver a vector of comma-separated ranges to numerical data.frame
# Input:
#	- A charactor vector r where
#		r[1] : comma-separted ranges i.e. "1-5,10-5" (1-base coordinate)
#		names(r): unique names
# Output:
# 	- a N x 3 data.frame of instances (1-base coordinate)
#		where row is each range
#		column is name / start / end
rangeListToBed=function(r){
	# r = df[,1]; names(r) = rownames(df)
	if(FALSE){
		r = df[,1]
		names(r) = rownames(df)
	}
	tmp = lapply(r, rangeToMask)
	nameL = unlist(sapply(names(tmp), function(x) rep(x, nrow(tmp[[x]]))))
	instances = do.call(rbind, tmp)
	instances = data.frame(Name = nameL, instances)
	return(instances)
}

## Convert a data.frame of motif instances in range format to mask
# instance / mask: all 1-base coordinate
# Input: N x M data.frame where
#	row is anchors
#	column is motif
#	elements are comma-separated ranges i.e. "1-5,10-15"
# Output: a list of mask data.frames of length M
#	
instanceFimoToMask=function(df){
	# df = df.fimo.filtered[[1]]
	maskL = list()
	for( motifName in colnames(df) ){
		# motifName = colnames(df)[3]
		instances = df[,motifName]
		names(instances) = rownames(df)
		instances = instances[ instances != "" ]
		stopifnot(length(instances) > 0)

		bed = rangeListToBed(instances)
		bed[,2] = bed[,2] - 1
		bed = bed[ order(bed$Name, bed$Start, bed$End), ]
		# merge to make mask
		mask = bedTools.1in("mergeBed", bed, direct=FALSE, debug=FALSE)
		mask[,2] = mask[,2] + 1
		colnames(mask) = c("Name","Start","End")

		maskL[[motifName]] = mask
	}

	return(maskL)
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


# make fimo results data.frames
makeFimoDataFrame=function(fimoL, anchorNameL, motifNameL){
	df.instance=data.frame(matrix("", length(anchorNameL), length(motifNameL)))
	rownames(df.instance) = anchorNameL
	colnames(df.instance) = motifNameL
	df.status=df.instance

	for(motifName in names(fimoL) ){
		# motifName = motifNameL[2]
		write(sprintf("Formatting %s", motifName), stderr())
		instance = fimoL[[motifName]][["instance"]]
		status = fimoL[[motifName]][["status"]]

		motifLen = nchar(instance[1,"Seq"])
		instance = instance[order(instance[,1], instance[,2]), c("Name","Offset")]
		instance = data.frame(Name = instance$Name,
					Start=instance$Offset,
					End=instance$Offset + motifLen -1)
		instance.groupByAnchor = lapply(instance %>% split(f = as.factor(.$Name), drop=TRUE), function(x) x[,2:3])
		instance.range = sapply(instance.groupByAnchor, maskToRange)

		df.instance[names(instance.range), motifName] = instance.range
		df.status[names(status), motifName] = status
	}

	return(list(instance = df.instance, status = df.status))
}


## Scan MEME motif using two p-values, stringent and marginal
# Will return 1-base coordinate
scanMemeMPRA=function(fa, meme, pv1, pv2, motifNameL, des.fimo="fimo.txt", verbose=FALSE){
	stopifnot(file.exists(meme))

	anchorNameL=NULL
	if(length(fa)==1 && is.character(fa) && file.exists(fa)){
		src = fa
		anchorNameL = names(readDNAStringSet(fa))
	}else if(class(fa)=="DNAStringSet"){
		src = tempfile()
		writeXStringSet(fa, src)
		anchorNameL = names(fa)
	}else{
		stopifnot("Currently, fasta file or DNAStringSet are allowed as input")
	}
	stopifnot(all(is.character(anchorNameL)))
	N.anchor = length(anchorNameL)

	# destination folder if necessary
	system(sprintf("mkdir -p %s", dirname(des.fimo)))

	cmd=sprintf("runMemeFimo.Peak.sh -g hg38 -p %e %s %s > %s", pv2, meme, src, des.fimo)
	system(cmd)
	if(class(fa)=="DNAStringSet") unlink(src)

	result = read.delim(des.fimo, header=TRUE, stringsAsFactors=FALSE, comment.char = "#", blank.lines.skip = TRUE)

	## Convert meme scan to legacy format of Homer
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
	motifScanL=list()
	for(motifName in motifNameL ){
		# motifName = motifNameL[2]
		write(sprintf("Formatting %s", motifName), stderr())
		result.select = result[result$Alt_Name == motifName,]
		stopifnot(nrow(result.select) > 0)
		# disjoint set of stringent vs marginal anchors and instances
		# i.e. marginal results doesn't have any anchors containing stringent instances
		
		# stringent instances
		result.stringent = result.select[result.select$PValue < pv1,]
		# anchor names having stringent instances
		anchor.stringent = unique(result.stringent$SeqName)
		# marginal instancees from anchors without stringent instances
		result.marginal = result.select[result.select$PValue >= pv1 & ! result.select$SeqName %in% anchor.stringent ,]

		# Select lowest p-value instance for each anchor as marginal instance
		# Select all if ties
		if(nrow(result.marginal)>0){
			result.marginal = as.data.frame(result.marginal %>% group_by(SeqName) %>% filter(PValue==min(PValue)))
		}else{
			write("Empty", stderr())
		}
		if(FALSE){
			# Test code 
			head(result.marginal[order(result.marginal$SeqName),],10)
			head(tmp[order(tmp$SeqName),],10)
		}
		instance.stringent = result.stringent[,c("SeqName","Start","Score","Direc","Seq")]
		instance.marginal = result.marginal[,c("SeqName","Start","Score","Direc","Seq")]
		# status vector
		status = rep("None", N.anchor)
		names(status) = anchorNameL
		status[unique(result.stringent$SeqName)] = "Stringent"
		status[unique(result.marginal$SeqName)] = "Marginal"

		tmp = list()
		tmp[["instance"]] = rbind(instance.stringent, instance.marginal)
		colnames(tmp[["instance"]]) = c("Name","Offset","Score","Direc","Seq")
		tmp[["status"]] = status
		motifScanL[[motifName]] = tmp
	}

	return(motifScanL)
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
# Output (example)
#	Sequence1    1 GACCTAGGACATTGCATACAGCTGATAGGACGAGCTTTGCAGCTGTAGCAGGCTGTGACAGGGGAAGACACCTTCCTGGCTAAGCAGTGTCCGAAACTCC  100
#	Mismatch                                                                                                           
#	Sequence2    1 GACCTAGGACATTGCATACAGCTGATAGGACGAGCTTTGCAGCTGTAGCAGGCTGTGACAGGGGAAGACACCTTCCTGGCTAAGCAGTGTCCGAAACTCC  100
#	Mask1                                                                                                              
#
#	Sequence1  101 TGTAATAGGACCGGGTCTGCTGGACCTGAGAGTGAAAGTGAGAGTGAGAGTGAGAGTGTGCCGCAAGGGAGGAAATGGGAGGGAAAGCGTCAAAACCAAC  200
#	Mismatch                                   * * *******     **  *  * *    *                    **** ** * *          
#	Sequence2  101 TGTAATAGGACCGGGTCTGCTGGACCTGTGGGGAGGTAAGAGAGATAGGGTAATAGTGAGCCGCAAGGGAGGAAATGGGGAACAGTGAGGCAAAACCAAC  200
#	Mask1                                     |----|----------|--------------|                 |------------|          
#
#	Sequence1  201 TCCTTTGGAGTGCATGATAAAAATTTTAAAGAAAGGATTTAGAGGTGATTGTGGGATGAAACTGGATGTT  270
#	Mismatch                                                                             
#	Sequence2  201 TCCTTTGGAGTGCATGATAAAAATTTTAAAGAAAGGATTTAGAGGTGATTGTGGGATGAAACTGGATGTT  270
#	Mask1
#
# Note:
#	If there is any overalapping masks, later ones will simply overwrites previous output string
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

compareCRSvsWT=function(crs, index, width=100, adapterLen=15){
	if(FALSE){
		crs=read.delim("mpra.design.txt", header=TRUE, stringsAsFactors=FALSE)

		compareCRSvsWT(crs, 100)
	}

	#index=6
	index.wt = which(crs$AnchorName == crs[index, "AnchorName"] & crs[, "MotifMutated"]=="WT")
	seq1 = crs[index.wt,"Sequence"]
	seq2 = crs[index,"Sequence"]
	mask = crs[index,"Mask"]
	write(sprintf("Anchor name: %s
Motif mutatedd: %s", crs[index, "AnchorName"], crs[index, "MotifMutated"]), stderr())

	if(crs[index, "MotifMutated"]=="WT" ){
		write(sprintf("This is WT seqeuence; skip"), stderr())
		return()
	}
	compareTwoStrings(substr(seq1, adapterLen+1, nchar(seq1)-adapterLen), substr(seq2, adapterLen+1, nchar(seq2)-adapterLen), mask, width=100)

}

## Draw a barplot of TF combo frequency comparing original TF coobinding vs motif scan intersection
# Input:
#	df: N x 3 data.frame where
#		Rows are comma-separted TF combinations
#		Columns are TF cobinding frequency before and after intersecting with motif scan
#			Cobind / Intersect / Downsample
#	remove Zero:
#		If TRUE, zero frequency combo will be ignored
drawComboFreqPlot=function(df, removeZero=TRUE){
	if(FALSE){
		df = df.freq
		removeZero=TRUE
	}
	library(ggplot2)
	library(reshape2)

	if(removeZero) df = df[!apply(df, 1, function(x) max(x)==0),]
	df.melt = reshape2::melt(as.matrix(df))
	colnames(df.melt) = c("Combo", "Step", "Frequency")
	df.melt$Step = factor(df.melt$Step, levels=c("Downsample", "Intersect", "Cobind"))
	df.melt$Combo = factor(df.melt$Combo, levels=rownames(df)[order(df$Cobind)])
	g = ggplot(df.melt, aes(x = Combo, y = Frequency, fill = Step)) +
		geom_col(position = "dodge") +
		scale_fill_discrete(breaks=rev(levels(df.melt$Step))) +
		geom_text(
			aes(label = Frequency),
			colour = "black", size = 3, hjust='left',
			position = position_dodge(1)
		) +
		ylim(0, max(df.melt$Frequency) * 1.2) +
		coord_flip()

	return(g)
}


## Calculate the total number of CRS to design except for PRDM1-only sites
# PRDM1-only sites are not counted here considering potential downsampling
# Input: a named vector of combo frequency
#	Each name is a comma-separated TF names
# Output: Total CRS
calcCrsCount=function(freq){
	# x = df.freq
	comboL = setdiff(names(freq)[freq > 0],c("None"))
	freq = freq[comboL]
	n = sum(sapply( comboL, function(x) 2^length(strsplit(x,",")[[1]])+1 ) * freq)
	return(n)
}




## Function to perform downsampling for each combinations of TF
# Input:
#	df: N x M boolean data.frame of cobinding 
#		row is anchor / column is TF (motif)
#	downsampleL: a list of downsampling setting
#		downsampleL[[combo name]] = <number of sites to sample>
# Output:
#	Downsampled df with a subset of rows
downsampleCombo = function(df, downsampleL){
	if(FALSE){
		df = df.cobind.filtered
	}

	# split by combo
	stopifnot(class(df) == "data.frame")
	comboNameL = boolTableToVector(df)
	dfL = split(df, comboNameL)
	
	# split into two groups
	# - to keep
	# - to downsample
	dfL.keep = dfL[! names(dfL) %in% names(downsampleL) ]
	dfL.sample = dfL[ names(dfL) %in% names(downsampleL) ]
	for( name in names(downsampleL) ){
		# name = names(downsampleL)[1]
		depth = downsampleL[[name]]
		tmp = dfL.sample[[name]][ sample(nrow(dfL.sample[[name]]), depth), ]
		dfL.keep[[name]] = tmp
	}
	# list -> data.frame
	df = do.call(rbind, setNames(dfL.keep, NULL))

	return(df)
}

## Not being used: too slow
# Merge overlapping mask ranges
# Input: a N x 2 data.frame of Start and End coordinate in 1-base
# Output: a M x 2 data.frame of merged Start and End coordinate in 1-base
mergeMasks=function(df){
	require(valr)
	tb = as_tibble(data.frame(chrom = "chr", start = df[,1], end = df[,2]))
	tb = bed_merge(tb)
	df.merged = as.data.frame(tb)[,-1]
	return(df.merged)
}




# Export MPRA design product
# Input:
#	- desPrefix : output prefix
#	- anchor : N x 6 data.frame of anchors in BED format used for MPRA design (after intersecting cobind x motif scan)
#	- crs : N_crs x 4 data.frame of MPRA design product. Output from generateMutatedSeqSet
#		where N_crs is the number of designed CRS
#	- motifStatus.raw: N x M data.frame of motif status (Stringent / Marginal / None) before intersecting cobind x motif scan
#		where M is the nurber of motifs
#	- motifStatus.filtered: N x M data.frame of motif status after intersecting cobind x motif scan
#	- maxMotifScore: N_crs x M data.frame of maximum motif score
# Output:
#	All sequences are flanked by adapters specified in "config"
#	- <desPrefix>.anchor.txt
#	- <desPrefix>.design.txt
#	- <desPrefix>.fa
#	- <desPrefix>.xlsx
exportDesign=function(desPrefix, config, anchor, crs, cobind, motifStatus.raw, motifStatus.filtered, maxMotifScore){
	if(FALSE){
		desPrefix = "mpra"
		anchor = scan.bed[anchorNameL.filtered,]
		crs = df.mut
		cobind = df.cobind[anchorNameL.filtered,]
		motifStatus.raw = df.fimo[["status"]][anchorNameL.filtered,]
		motifStatus.filtered = df.fimo.filtered[["status"]][anchorNameL.filtered,]
		maxMotifScore = maxMotifScore.crs
	}

	# validation
	stopifnot(all(rownames(anchor) == rownames(motifStatus.raw)))
	stopifnot(all(rownames(anchor) == rownames(motifStatus.filtered)))
	stopifnot(all(rownames(crs) == rownames(maxMotifScore)))
	stopifnot(all(c("adapter5","adapter3") %in% names(config)))

	# Final CRS data sheet
	crs$Sequence = sprintf("%s%s%s", config[["adapter5"]], crs$Sequence, config[["adapter3"]])
	df.crs = data.frame(Id = rownames(crs), crs, maxMotifScore)

	# Final anchor data sheet
	colnames(anchor) = c("Chr","Start","End","Name","Null","Direc")
	colnames(cobind) = sprintf("Cobind.%s", colnames(cobind))
	colnames(motifStatus.raw) = sprintf("Motif.%s", colnames(motifStatus.raw))
	colnames(motifStatus.filtered) = sprintf("Motif_Filtered.%s", colnames(motifStatus.filtered))
	df.anchor = data.frame(
		anchor,
		cobind,
		motifStatus.raw,
		motifStatus.filtered
	)

	## Export Output

	# 1. anchor sheet in txt file
	write.table(df.anchor, sprintf("%s.anchor.txt", desPrefix), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
	# 2. MPRA CRS sheet in txt file
	write.table(df.crs, sprintf("%s.design.txt", desPrefix), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

	# 3. excel file
	# - design parameters
	# - Two comprehensive data sheet from #1 & #2
	df.config = data.frame(Item=character(), Value=character(), Value2=character())
	df.config["Title",] = c("Title", config[["main_title"]], "")
	df.config["Anchor",] = c("Anchor", config[["src_anchor"]], "")
	df.config["TF Cobind",] = c("TF Cobind", config[["src_cobind"]], "")
	df.config["Motif file",] = c("Motif file", config[["src_motif"]], "")
	df.config["Scan Window",] = c("Scan Window", config[["scan_window"]], "")
	df.config["Genome",] = c("Genome", config[["genome"]], "")
	df.config["P-value stringent",] = c("P-value stringent", config[["pv1"]], "")
	df.config["P-value marginal",] = c("P-value marginal", config[["pv2"]], "")
	
	df.config["Motif Names",] = c("Motif Names", "", "")
	motifNameL = config[["motiflist"]]
	for( i in 1:length(motifNameL) ){
		motifName = motifNameL[i]
		df.config[sprintf("Motif%d", i),] = c("", motifName, "")
	}

	downsampleL = config[["downsample"]]
	df.config["Downsampling",] = c("Downsampling", "", "")
	for( comboName in names(downsampleL) ) df.config[comboName,] = c("", comboName, sprintf("%d", downsampleL[[comboName]]))
	
	df.config["Mutation Mode",] = c("Mutation Mode", config[["mutation_mode"]], "")
	df.config["adapter5",] = c("adapter5", adapter5, "")
	df.config["adapter3",] = c("adapter3", adapter3, "")


	wb = createWorkbook("MPRA-ChIP")
	addWorksheet(wb, "Parameters")
	writeData(wb, sheet = 1, df.config)
	addWorksheet(wb, "Anchor")
	writeData(wb, sheet = 2, df.anchor)
	addWorksheet(wb, "Design")
	writeData(wb, sheet = 3, df.crs)
	saveWorkbook(wb, sprintf("%s.xlsx", desPrefix), overwrite = TRUE)


	# 3. fasta file of all CRS
	fa.final = DNAStringSet(df.crs$Sequence)
	names(fa.final) = df.crs$Id
	writeXStringSet(fa.final, sprintf('%s.crs.fa', desPrefix), width=100)

}
