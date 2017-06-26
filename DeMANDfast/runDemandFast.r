####################################
#
#	Alexander Lachmann 7/6/2015
#	DeMAND reimagined
#	
#	original algorithm: Jung Hoon
#	
#	
#
############################################
##############################################
runDeMANDFast <- function (dobj, fgIndex=NULL, bgIndex=NULL){
	write("RunDeMANDfast: Preparing Input", stderr())
	fgIndex = fgIndex
	bgIndex = bgIndex

	x = dobj

	expData <- x@exp
	expAnno <- x@anno
	inputNetwork <- x@network
	platformGene <- unique(expAnno[, 2])

	ww = which(inputNetwork[,1] %in% platformGene & inputNetwork[,2] %in% platformGene)
	interactome = inputNetwork[ww,]

	analGene <- intersect(as.vector(interactome), platformGene)

	rsd = apply(expData,1,sd)/apply(expData,1,mean)
	oo = rev(order(rsd))

	expData = expData[oo,]
	expAnno = expAnno[oo,]
	mm = match(unique(expAnno[,2]), expAnno[,2])

	expData = expData[mm,]
	rownames(expData) = expAnno[mm,2]

	tempExp = expData

	expData = t(apply(expData[,c(fgIndex, bgIndex)],1,rank, ties.method='first'))

	fgIndex = 1:length(fgIndex)
	bgIndex = (length(bgIndex)+1):(length(bgIndex)+length(fgIndex))

	# To make random network
	message("RunDeMANDfast: Prepare Network and compute Kernel-bandwidth....")
	p1 <- sample(interactome[,1])
	p2 <- sample(interactome[,2])

	#remove self-loops
	pOut <- which(p1 == p2)
	#remove self loops, they are however allowed in the later p-value analysis
	permutedInteractome <- cbind(p1, p2)[-pOut,]

	#compute the kernel-bandwidth
	sigmaFG = round(kernelDensityBandwidth(expData[,fgIndex]))^2
	names(sigmaFG) = rownames(expData)
	sigmaBG = round(kernelDensityBandwidth(expData[,bgIndex]))^2
	names(sigmaBG) = rownames(expData)
	
	message("RunDeMANDfast: Measure dysregulation of the interactions.....")
	message("Going Through the Interactome")
	#compute KLD for all edges in network
	edgeKLD = getSigKLDedge(interactome, expData, sigmaFG, sigmaBG, fgIndex, bgIndex)
        message("RunDeMANDfast: Make a null distribution for KL divergence.....")
	#compute null distribution with randomized network
	nullmodel = getSigKLDedge(permutedInteractome, expData, sigmaFG, sigmaBG, fgIndex, bgIndex)

	pval = getKLDpvalue(edgeKLD[,3], nullmodel[,3])
	edgeKLD = cbind(edgeKLD, pval)

	#get integrated p-value
	#intPval = getCombinedPvalue(edgeKLD, expData)
	#in the demand code the actual expression is used and not the rank transformed. Also the full expression matrix is used rather than the
	#control and treatment subsamples.
	message("RunDeMANDfast: Estimate dysregulation of the genes.....")
	intPval = getCombinedPvalue(edgeKLD, tempExp)

	#correct for multiple hypothesis testing
	intPvalAdjustedp <- p.adjust(intPval, "fdr")

	#make a nice data.frame with the genes and their significance
	finalPrediction <- data.frame(moaGene = names(intPval), Pvalue = intPval, adjustedPvalue = intPvalAdjustedp)
	finalPrediction <- finalPrediction[sort(intPval, decreasing = F, index.return = T)$ix, ]

	# To update moa slot in the demand object
	dobj@moa <- finalPrediction

	#that's what we get
	return(dobj)
}
