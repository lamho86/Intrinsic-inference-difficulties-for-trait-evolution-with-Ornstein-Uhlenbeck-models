### July 23 2014
### Lam Ho (lamho86@gmail.com)
# aic: Akaike information criterion
# bic: Bayesian information criterion
# saic: SURFACE like AIC
# sbic: SURFACE like BIC
# mbic: modified BIC


OUshifts <- function(y, phy, method=c("mbic","aic","bic","saic","sbic"), nmax, check.pruningwise = TRUE) {
	if (!inherits(phy, "phylo")) 
      	stop("object \"phy\" is not of class \"phylo\".")
	if (check.pruningwise)
		phy = reorder(phy, "pruningwise")
	method = match.arg(method)
	options(warn = -1)
	n <- length(phy$tip.label)
	N <- dim(phy$edge)[1]
	ROOT <- n + 1L
	anc <- phy$edge[, 1]
	des <- phy$edge[, 2]
	el <- phy$edge.length
	v <- matrix(0,n + phy$Nnode,n)
	for (i in 1:N) {
		if (des[i] <= n) v[des[i],des[i]] = 1
		v[anc[i],] = v[anc[i],]+v[des[i],]		
	}

	check <- function(model) {
		### return identifiability check and modified BIC penalty
		if (length(model)==0) return(c(1,log(n)))
		checkpar = rep(ROOT,n)
		for (i in N:1) 
			if (i %in% model) checkpar[which(v[des[i],]==1)] = i
		if (length(levels(as.factor(checkpar))) == length(model)+1) {
			numchange = length(model)
			pen = 3*log(n) + (2*numchange - 1)*log(n)
			for (j in 1:numchange)
				pen = pen + log(length(which(as.factor(checkpar)==levels(as.factor(checkpar))[j])))
			return(c(1,pen)) 
			} else return(c(0,0))
	}	

	score <- function(model,method) {
		if (length(model)==0) {			
			logLik = phylolm(y~1,phy=phy,model="OUfixedRoot")$logLik
			if (method=="aic") score = -2*logLik + 2*(length(model)+3)
			if (method=="bic") score = -2*logLik + log(n)*(length(model)+3)
			if (method=="mbic") score = -2*logLik + 3*log(n)
			if (method=="saic") score = -2*logLik + 2*(length(model)+3)
			if (method=="sbic") score = -2*logLik + log(n)*(length(model)+3)
			return(score)
		}
		X = v[ROOT,]
		for (i in 1:length(model)) X = cbind(X,v[des[model[[i]]],])
		ch = check(model)
		if (ch[1]==0) return(Inf)
		logLik = phylolm(y~X-1,phy=phy,model="OUfixedRoot")$logLik
		if (method=="aic") score = -2*logLik + 2*(length(model)+3)
		if (method=="bic") score = -2*logLik + log(n)*(length(model)+3)
		if (method=="mbic") score = -2*logLik + ch[2]
		if (method=="saic") score = -2*logLik + 2*(2*length(model)+3)
		if (method=="sbic") score = -2*logLik + log(n)*(2*length(model)+3)
		return(score)
	}

	curmodel = list()	
	curscore = score(curmodel,method)
	flag = 0
	promodel = curmodel
	proscore = curscore
	while (flag==0) {		
		for (i in 1:N) {
			if (sum(which(curmodel==i))>0) {	
				pos = which(curmodel==i)
				tempmodel = curmodel[-pos]				
				tempscore = score(tempmodel,method)								
				if (tempscore < proscore) {
					promodel = tempmodel
					proscore = tempscore
					flag = flag + 1
				}
			}
			if (sum(which(curmodel==i))==0) {
				if (length(curmodel) < nmax) {
					tempmodel = curmodel
					tempmodel[[length(tempmodel)+1]] = i
					tempscore = score(tempmodel,method)							
					if (tempscore < proscore) {
						promodel = tempmodel
						proscore = tempscore
						flag = flag + 1
					}
				}
				if (length(curmodel)>=1)
					for (j in 1:length(curmodel)) 
						if ((anc[i]==des[curmodel[[j]]])||(des[i]==anc[curmodel[[j]]])||(anc[i]==anc[curmodel[[j]]])) {
							tempmodel = curmodel
							tempmodel[j] = i
							tempscore = score(tempmodel,method)									
							if (tempscore < proscore) {
							promodel = tempmodel 
							proscore = tempscore
							flag = flag + 1
							}
						}	
			}		
		}
		if (flag > 0) flag = 0 else flag = 1
		curmodel = promodel
		curscore = proscore	
	}
	results = list(y = y, phy = phy, score=curscore, nmax = nmax, nshift = length(curmodel))
	if (results$nshift==0) results$pshift = NULL else results$pshift = unlist(curmodel)
	class(results) = "OUshifts"
	return(results)
	}

################################################
plot.OUshifts <-function(object, show.data = TRUE) {
	if (show.data)
		layout(matrix(1:2,1,2))
	plot(object$phy, show.tip.label=F)	
	edgelabels(edge=object$pshift,pch=8,frame="n",cex=1.1,lwd=2,col="red")
	if (show.data) {		
		object$y = as.vector((object$y-mean(object$y))/sd(object$y))
		barplot(object$y, horiz=TRUE, xlab="Normalized trait", names.arg="")
		}
	}






