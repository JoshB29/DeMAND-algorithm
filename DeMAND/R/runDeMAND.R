#' DEMAND Function
#'
#'  This function is based on the realization that drugs affect the protein activity of their targets, but not necessarily their mRNA expression levels. In contrast, the change in protein activity directly affects the mRNA expression levels of downstream genes. Based on this hypothesis, DeMAND identifies drug MoA by comparing gene expression profiles following drug perturbation with control samples, and computing the change in the individual interactions within a pre-determined integrated transcriptional and post-translational regulatory model (interactome). For each edge in the interactome we determine the two-dimensional probability distribution of the gene expression levels both in the control state, and following drug treatment. Any changes in the probability distribution are estimated using the Kullback-Leibler (KL) divergence, from which we determine the statistical significance of the dysregulation of each edge. In the second step of DeMAND, we interrogate each gene independently to determine whether its interactions are enriched in dysregulated ones, suggesting that it is a candidate mechanism of action.
#'
#' @param dobj Instance of class demand
#' @param fgIndex Sample indices of Drug treated samples
#' @param bgIndex Sample indices of DMSO treated samples
#' @return Objet of class demand with updated moa slot
#' @docType methods
#' @examples
#' data(inputExample, package="demand")
#' dobj <- demandClass(exp=bcellExp, anno=bcellAnno, network=bcellNetwork)
#' dobj <- runDeMAND(dobj, fgIndex=caseIndex, bgIndex=controlIndex)
#' show(dobj)
#' @export


runDeMAND <- function (x, fgIndex=NULL, bgIndex=NULL){

	# Parameter check
	if (is.null(fgIndex) | is.null(bgIndex)) stop("Please provide sample (column of expression data) indices of case/control samples")
	if(length(fgIndex) < 3 | length(bgIndex) < 3) stop("The number of samples in each class should be at least three")
	if((length(fgIndex) >= 3 & length(fgIndex) < 6) | (length(bgIndex) >= 3 & length(bgIndex) < 6)) warning("it will run but please note that DeMAND needs minimum six samples (in each class) to show the best performance")

	message("Preparing input")
	expData <- x@exp
	expAnno <- x@anno
	inputNetwork <- x@network
	platformGene <- unique(expAnno[, 2])

    	tempIdx <- apply(inputNetwork, 1, tempF1 <- function(yy) {
        	if (length(intersect(yy[1:2], platformGene)) == 2) {return(1)}
        	else {return(0)}
    		})

	# To reduce the networks using genes appear in the expression data
    	interactome <- inputNetwork[which(tempIdx ==1), 1:2]
    	analGene <- intersect(unique(as.vector(interactome)), platformGene)

	# To preprocess the expression data. if there are multiple probes for one gene, select one with the maximum coefficient of variation
	tempIdx <- c()
	flush <- apply(as.matrix(analGene), 1, tempF2 <- function(g) {
		probeIdx <- which(expAnno[, 2] == g)
        	if (length(probeIdx) == 1) {tempIdx <<- c(tempIdx, probeIdx)}
        	else {
			tempExpMat <- expData[probeIdx, ]
			tempCV <- apply(tempExpMat, 1, tempF3 <- function(yy) {
				rs <- sd(yy, na.rm=T)/mean(yy, na.rm=T)
                		return(rs)
				})
			tempIdx <<- c(tempIdx, probeIdx[which(tempCV == max(tempCV, na.rm=T))[1]])
        		}
		})
    	expData <- expData[tempIdx, ]
    	rownames(expData) <- analGene

	# internal function definition
	getCombinedPvalue <- function(g) {

		# get neighbor information of a given gene
		nbsIdx <- nbsIdxAll[[g]]
		if (length(nbsIdx) < 2) {return(1)}
		nbs <- setdiff(as.character(edgeKLD[nbsIdx, 1:2]), g)
        	if (length(nbs) < 2) {return(1)}

		# get p-values of edges surrounded using Fisher's method
        	edgeKLDnb <- edgeKLD[nbsIdx, ]
        	nPvec <- apply(as.matrix(nbs), 1, tempF4 <- function(nb) {
			edgeKLDnb[which(edgeKLDnb[, 1] == nb | edgeKLDnb[,2] == nb)[1], 4]
        		})
		
		# p-value integration by Fisher's method
        	nPvec <- as.numeric(nPvec)
        	nPvec[which(nPvec == 0)] <- 10^-20
        	intpChisq <- -2 * sum(log(nPvec))
        
		# Correcting a possible depedency structure among the interactions using Brown's method
		gExpVec <- expData[g, ]
		residualExp <- list()
		flush <- apply(as.matrix(nbs), 1, tempF5 <- function(n) {
			y <- as.numeric(expData[n, ])
			x <- as.numeric(gExpVec)
            		xp <- cbind(x, rep(1, length(x)))
            		res <- y - (xp %*% ginv(xp) %*% y)
            		res <- as.vector(res)
            		residualExp[[n]] <<- res
        		})
        	combNbs <- combn(nbs, 2)
        	combNbsCorr <- apply(combNbs, 2, tempF6 <- function(nn) {
	            	ppiIdxN1 <- ppiIdx[[g]][which(ppiIdx[[g]][, 1] == nn[1]), 2]
            		ppiIdxN2 <- ppiIdx[[g]][which(ppiIdx[[g]][, 1] == nn[2]), 2]
            		if (ppiIdxN1 == 1 | ppiIdxN2 == 1) {return(0)}
            		else {
                		resExpN1 <- residualExp[[nn[1]]]
	          	      	resExpN2 <- residualExp[[nn[2]]]
		                resCorr <- abs(as.vector(cor(resExpN1, resExpN2)))
		            }
	        	    return(resCorr)
		 	})
        	meanChisq <- 2 * length(nbs)
        	covChisq <- combNbsCorr * (3.25 + 0.75 * combNbsCorr)
        	varChisq <- 4 * length(nbs) + 2 * sum(covChisq)
	        f = (2 * (meanChisq)^2)/varChisq
	        c = varChisq/(2 * (meanChisq))
	        correctedChisq <- intpChisq/c
	        rs <- pchisq(correctedChisq, df = f, lower.tail = FALSE)
	        return(rs)
	    	}

	# To get bandwith for Gaussian kernel
	kernelDensityBandwidth <- function(data, varargin) {
        	prop <- 1.06
	        dim <- 1
	        sig <- sd(data)
	        iqrSig <- 0.7413 * IQR(data)
	        if (max(iqrSig) == 0) {iqrSig <- sig}
        	sig = min(sig, iqrSig)
	        n = length(data)
	        h = prop * sig * n^(-1/(4 + dim))
	        return(h)
    		}

	# To get probability from multivariate normal distribution
    	mvn <- function(arg, mu, covariancematrix) {
        	x = arg[1]
	        y = arg[2]
        	y = exp(-0.5 * (((x - mu[1])/sqrt(covariancematrix[1,1]))^2 + ((y - mu[2])/sqrt(covariancematrix[2, 2]))^2))
        	y = y/(2 * pi * sqrt(covariancematrix[1, 1]) * sqrt(covariancematrix[2, 2]))
        	return(y)
    		}

	# To get probability density (2-dimensional) from expression data using Gaussian kernel
    	getPDM <- function(xdata, ydata, xrange, yrange, sigmaX, sigmaY, covarianceMatrix) {
        	sumMvnPmatrix <- matrix(0, ncol = max(yrange), nrow = max(xrange))
	        colnames(sumMvnPmatrix) <- as.character(xrange)
	        rownames(sumMvnPmatrix) <- as.character(rev(yrange))
        	for (n in c(1:length(xdata))) {
			muX <- xdata[n]
            		muY <- ydata[n]
            		mvnPmatrix <- globalProbMatrix[[paste(muX, "-", muY, "-", sigmaX, "-", sigmaY, sep = "")]]
            		if (is.null(mvnPmatrix)) {
                		mvnPmatrix <- matrix(0, ncol = max(yrange), nrow = max(xrange))
                		colnames(mvnPmatrix) <- as.character(xrange)
                		rownames(mvnPmatrix) <- as.character(rev(yrange))
               			for (i in xrange) {
                  			for (j in yrange) {
                    				tempMvnP <- mvn(c(i, j), c(muX, muY), covarianceMatrix)
                    				mvnPmatrix[as.character(j), as.character(i)] <- tempMvnP
                  				}
                			}
                		globalProbMatrix[[paste(muX, "-", muY, "-", sigmaX, "-", sigmaY, sep = "")]] <<- mvnPmatrix
                		sumMvnPmatrix <- sumMvnPmatrix + mvnPmatrix
            			}else {
                		sumMvnPmatrix <- sumMvnPmatrix + mvnPmatrix
            			}
        		}
        	return(sumMvnPmatrix)
    	}
	
    	getSigKLDedge <- function(yy, bgIndex, fgIndex) {

		# To get two genes interacted
	        edge1 <- yy[1]
	        edge2 <- yy[2]

		# To get expression profiles of two genes
        	tempX <- expData[edge1, c(bgIndex, fgIndex)]
        	tempY <- expData[edge2, c(bgIndex, fgIndex)]

		# To rank transformation
        	x <- rank(tempX, ties.method = "first")
        	y <- rank(tempY, ties.method = "first")

		# To make 2-dimensional matrix
	        bgIndexTemp <- c(1:length(bgIndex))
        	fgIndexTemp <- c((length(bgIndex) + 1):(length(bgIndex) + length(fgIndex)))
        	xdataBG <- x[bgIndexTemp]
        	xdataFG <- x[fgIndexTemp]
        	ydataBG <- y[bgIndexTemp]
        	ydataFG <- y[fgIndexTemp]
        	xrange <- c(1:max(x))
        	yrange <- c(1:max(y))

		# To get bandwith with optimization
        	sigmaXbg <- round(kernelDensityBandwidth(xdataBG))
		sigmaYbg <- round(kernelDensityBandwidth(ydataBG))
        	if (sigmaXbg == 0) {sigmaXbg <- 0.1}
		if (sigmaXbg <= 0.5) {sigmaXbg <- 0.1}
        	if (sigmaYbg == 0) {sigmaYbg <- 0.1}
        	if (sigmaYbg <= 0.5) {sigmaYbg <- 0.1}
        	covarianceMatrixBG <- matrix(c(sigmaXbg^2, 0, 0, sigmaYbg^2), nrow = 2, byrow = T)
        	
		sigmaXfg <- round(kernelDensityBandwidth(xdataFG))
        	sigmaYfg <- round(kernelDensityBandwidth(ydataFG))
		if (sigmaXfg == 0) {sigmaXfg <- 0.1}
        	if (sigmaXfg <= 0.5) {sigmaXfg <- 0.1}
		if (sigmaYfg == 0) {sigmaYfg <- 0.1}
        	if (sigmaYfg <= 0.5) {sigmaYfg <- 0.1}
        	covarianceMatrixFG <- matrix(c(sigmaXfg^2, 0, 0, sigmaYfg^2), nrow = 2, byrow = T)

		# To get probability distribution of case/control using Gaussian kernel method
        	sumMvnPmatrixBG <- getPDM(xdataBG, ydataBG, xrange, yrange, sigmaXbg, sigmaYbg, covarianceMatrixBG)
        	sumMvnPmatrixFG <- getPDM(xdataFG, ydataFG, xrange, yrange, sigmaXfg, sigmaYfg, covarianceMatrixFG)
        	sumMvnPmatrixBG[which(sumMvnPmatrixBG == 0)] <- 10^-100 # smoothing
        	sumMvnPmatrixFG[which(sumMvnPmatrixFG == 0)] <- 10^-100 # smoothing

		# Normalization
        	sumMvnPmatrixBG <- sumMvnPmatrixBG/sum(sumMvnPmatrixBG)
        	sumMvnPmatrixFG <- sumMvnPmatrixFG/sum(sumMvnPmatrixFG)

		# To measure KL-divergence between the two distributions, and make the value symmetric
        	kld1 <- sum(sumMvnPmatrixBG * log(sumMvnPmatrixBG/sumMvnPmatrixFG))
        	kld2 <- sum(sumMvnPmatrixFG * log(sumMvnPmatrixFG/sumMvnPmatrixBG))
        	kld <- kld1 + kld2
        	rs <- c(edge1, edge2, kld)
        	return(rs)
    		}

	# To get statistical significance of the KLD value using null distribution from the randomized network
	getKLDpvalue <- function(kld) {
        	kldValue <- as.numeric(kld[3])
	        rs <- length(which(null > kldValue))/length(null)
	        return(rs)
		}

	# To reduce redundancy in probability calculation 
    	globalProbMatrix <- list()
    
	message("Make a null distribution for KL divergence.....")

	# To make random network
	p1 <- sample(analGene, length(analGene), replace = T)
    	p2 <- sample(analGene, length(analGene), replace = T)
    	pOut <- which(p1 == p2)
    	if (length(pOut) > 0) {p1 <- p1[-pOut]; p2 <- p2[-pOut]}
    	permuteInteractome <- cbind(p1, p2)

    	# To get null distribution of KLD values
	tempNull <- apply(permuteInteractome, 1, getSigKLDedge, bgIndex, fgIndex)
    	null <- as.numeric(tempNull[3, ])
    
	message("Measure dysregulation of the interactionss.....")

	# To get KLD for all the edges
    	KLDmat <- apply(interactome, 1, getSigKLDedge, bgIndex, fgIndex)
    	KLDmat <- t(KLDmat)

	# To get pvalue (of KLD) for all the edges
    	KLDpvec <- apply(KLDmat, 1, getKLDpvalue)
    	edgeKLD <- cbind(KLDmat, KLDpvec)
    	colnames(edgeKLD) <- c("gene1", "gene2", "KLD", "KLD.p")
    
	# Indexing for optimization
	nbsIdxAll <- list()
    	flush <- apply(as.matrix(analGene), 1, tempF7 <- function(g) {
        	nbsIdxAll[[g]] <<- which(edgeKLD[, 1] == g | edgeKLD[, 2] == g)
    		})

    	ppiIdx <- list()
    	flush <- apply(as.matrix(analGene), 1, tempF8 <- function(g) {
        	nbsIdx <- nbsIdxAll[[g]]
        	nbs <- setdiff(as.character(edgeKLD[nbsIdx, 1:2]), g)
        	nbsPpiIdx <- apply(as.matrix(nbs), 1, tempF9 <- function(n, g) {
            	rs <- inputNetwork[, "ppi"][which((interactome[, 1] == n & interactome[, 2] == g) | (interactome[, 1] == g & interactome[, 2] == n))][1]
            	return(rs)
        	}, g)
        	ppiIdx[[g]] <<- cbind(nbs, nbsPpiIdx)
    		})

	# To combine p-values
    	message("Estimate dysregulation of the genes.....")
    	intPval <- apply(as.matrix(analGene), 1, getCombinedPvalue)
    	intPvalAdjustedp <- p.adjust(intPval, "fdr")
    	finalPrediction <- data.frame(moaGene = analGene, Pvalue = intPval, adjustedPvalue = intPvalAdjustedp)
    	finalPrediction <- finalPrediction[sort(intPval, decreasing = F, index.return = T)$ix, ]

	# To update moa slot in the demand object
	x@moa <- finalPrediction
	return(x)
	}

