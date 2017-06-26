library("DeMAND")
source("runDemandFast.r")
source("demandfunctions.r")

data(inputExample)


dobj <- demandClass(exp=bcellExp, anno=bcellAnno, network=bcellNetwork)
printDeMAND(dobj)

res <- runDeMANDFast(dobj, fgIndex=caseIndex, bgIndex=controlIndex)
printDeMAND(res)

