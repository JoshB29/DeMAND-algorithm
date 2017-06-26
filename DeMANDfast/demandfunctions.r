####################################
#
#   Alexander Lachmann 7/6/2015
#   DeMAND reimagined
#   
#   original algorithm: Jung Hoon
#   
#   Collection of rewritten DeMAND functions.
#   Functionality unchanged from original.
#
############################################
##############################################


# To get bandwith for Gaussian kernel, as matrix operation
kernelDensityBandwidth <- function(data) {
	prop <- 1.06
    sig <- apply(data,1,sd)

    iqrSig <- 0.7413 * apply(data,1,IQR)
    iqrSig[which(iqrSig == 0)] = sig[which(iqrSig == 0)]

    sig = apply(cbind(sig,iqrSig),1,min)

    n = ncol(data)
    return(prop * sig * n^(-1/5))
}


# To get probability from multivariate normal distribution
mvn <- function(x,y, mu, sigma1, sigma2) {
    rs = exp(-0.5 * (((x - mu[1])/sqrt(sigma1))^2 + ((y - mu[2])/sqrt(sigma2))^2))
    rs = rs/(2 * pi * sqrt(sigma1) * sqrt(sigma2))
    return(rs)
}

getPDM <- function(expData, r1, r2, xdata, ydata, sig1, sig2) {
    mvnMatrix = rep(0, ncol(expData)^2)
    for(i in 1:length(xdata)){
        mvnMatrix = mvnMatrix + mvn(r1,r2, c(xdata[i], ydata[i]), sig1, sig2)
    }
    mvnMatrix[which(mvnMatrix == 0)] <- 10^-100
    return(mvnMatrix/sum(mvnMatrix))
}

#compute KLD for all edges of interactome
getSigKLDedge <- function(interactome, expData, sigmaFG, sigmaBG, fgIndex, bgIndex) {

    r1 = rep(1:ncol(expData), ncol(expData))
    r2 = rep(1:ncol(expData), each = ncol(expData), times = 1)

    edgeKLD = list()

    for(i in 1:nrow(interactome)){

        fgPD = getPDM(expData,r1, r2, expData[interactome[i,1],fgIndex], expData[interactome[i,2],fgIndex], sigmaFG[interactome[i,1]], sigmaFG[interactome[i,2]])
        bgPD = getPDM(expData,r1, r2, expData[interactome[i,1],bgIndex], expData[interactome[i,2],bgIndex], sigmaBG[interactome[i,1]], sigmaBG[interactome[i,2]])

        #fgPD = getPDM(r1, r2, expData[interactome[i,1],fgIndex], expData[interactome[i,2],fgIndex], 1, 1)
        #bgPD = getPDM(r1, r2, expData[interactome[i,1],bgIndex], expData[interactome[i,2],bgIndex], 1, 1)

        kld1 <- sum(fgPD * log(fgPD/bgPD))
        kld2 <- sum(bgPD * log(bgPD/fgPD))
        kld <- kld1 + kld2
        edgeKLD[[length(edgeKLD)+1]] = c(interactome[i,1], interactome[i,2], kld)
    	message(i) #To tell you how far in the interactome you are
	}
    return(do.call(rbind, edgeKLD))
}

# To get statistical significance of the KLD value using null distribution from the randomized network
getKLDpvalue <- function(kld, nullmodel) {
    kldValue <- as.numeric(kld)
    null <- as.numeric(nullmodel)

    rs = c()
    for(i in 1:length(kld)){
         rs <- c(rs, length(which(null > kldValue[i]))/length(null))
    }
    return(rs)
}

# this method combined p-values with the brown's method. It takes dependencies of variables into consideration.
getCombinedPvalue <- function(edgeKLD, expData) {

    genes = table(c(edgeKLD[,1:2]))

    result = rep(0, length(genes))
    names(result) = names(genes)

    #genes that have only one edge in the network are by default not significant
    result[which(genes < 2)] = 1

    genes = names(genes)[which(genes >= 2)]

    for(g in genes){
        w1 = which(edgeKLD[,1] == g)
        w2 = which(edgeKLD[,2] == g)
        
        neighborList = unique(c(edgeKLD[w1,2], edgeKLD[w2,1]))

        pvals = as.numeric(edgeKLD[c(w1,w2),4])
        pvals[which(pvals == 0)] <- 10^-20
        intpChisq <- -2 * sum(log(pvals))

        x = expData[g,]
        xp <- cbind(x, rep(1, length(x)))
        xp <- xp %*% ginv(xp)

        residualList = list()
        for(neighbor in neighborList){

            y <- expData[neighbor,]

            #magical way to compute residuals, does the same as lm(y~x) but is indeed faster. What a time to be alive!
            residuals <- y - (xp %*% y)
            residualList[[length(residualList)+1]] = residuals
        }
        residuals <- do.call(cbind,residualList)

        residualCorrelation = cor(residuals)
        residualCorrelation = c(residualCorrelation[upper.tri(residualCorrelation)])

        #math magic 
        meanChisq <- 2 * length(neighborList)
        covChisq <- residualCorrelation * (3.25 + 0.75 * residualCorrelation)
        varChisq <- 4 * length(neighborList) + 2 * sum(covChisq)
        f = (2 * (meanChisq)^2)/varChisq
        c = varChisq/(2 * (meanChisq))
        correctedChisq <- intpChisq/c
        result[g] <- pchisq(correctedChisq, df = f, lower.tail = FALSE)
    }

    return(result)
}












