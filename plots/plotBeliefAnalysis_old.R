library(scales)

almost.equal <- function (x, y, tolerance=.Machine$double.eps^0.5, 
		na.value=TRUE)  { 
	answer <- rep(na.value, length(x)) 
	test <- !is.na(x) 
	answer[test] <- abs(x[test] - y) < tolerance 
	answer 
}

plotAllofAll_belief3 <- function(prefix, nodes=31){
	links <- c(0.3,0.5,0.7)
	edges <- c(4,8)
	methods <- c("belief", "conf-perfect")
	acc <- 0.6
	stddev <- 0.12
	maxSpeed <- 20
	for (edge in edges){
		epsilons <- seq(0.07,0.1,0.001)
		pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/belief-acc-BA-e",edge,".pdf",sep=""))
		plotSuccessOnEpsilon(prefix=prefix, nodes=nodes, netType="barabasi-albert", accuracy=acc, acstdv=stddev,
				methods=methods, edge=edge, link=links[1], colours=rainbow(3), yrange=c(0,1), epsilons=epsilons)
		dev.off()
		pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/belief-speed-BA-e",edge,".pdf",sep=""))
		plotTimeOnEpsilon(prefix=prefix, nodes=nodes, netType="barabasi-albert", accuracy=acc, acstdv=stddev,
				methods=methods, edge=edge, link=links[1], colours=rainbow(3), yrange=c(1,maxSpeed), epsilons=epsilons)
		dev.off()
		pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/belief-trade-BA-e",edge,".pdf",sep=""))
		plotSpeedAccuracyTradeoff(prefix=prefix, nodes=nodes, netType="barabasi-albert", accuracy=acc, acstdv=stddev, addOptim=TRUE,
				methods=methods, edge=edge, link=links[1], colours=rainbow(3), xrange=c(1,maxSpeed), yrange=c(0,1), epsilons=epsilons)
		dev.off()
	}
	for (link in links){
		if (link == 0.3){
			epsilons <- seq(0.09,0.11,0.001)
		} else if (link == 0.5){
			epsilons <- seq(0.07,0.09,0.001)
		} else if (link == 0.7){
			epsilons <- seq(0.05,0.07,0.001)
		}
		
		linkstr <- gsub("\\.", "", toString(link))
		pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/belief-acc-ER-l",linkstr,".pdf",sep=""))
		plotSuccessOnEpsilon(prefix=prefix, nodes=nodes, netType="erdos-renyi", accuracy=acc, acstdv=stddev,
				methods=methods, edge=edges[1], link=link, colours=rainbow(3), yrange=c(0,1), epsilons=epsilons)
		dev.off()
		pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/belief-speed-ER-l",linkstr,".pdf",sep=""))
		plotTimeOnEpsilon(prefix=prefix, nodes=nodes, netType="erdos-renyi", accuracy=acc, acstdv=stddev,
				methods=methods, edge=edges[1], link=link, colours=rainbow(3), yrange=c(1,maxSpeed), epsilons=epsilons)
		dev.off()
		pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/belief-trade-ER-l",linkstr,".pdf",sep=""))
		plotSpeedAccuracyTradeoff(prefix=prefix, nodes=nodes, netType="erdos-renyi", accuracy=acc, acstdv=stddev, addOptim=TRUE,
				methods=methods, edge=edges[1], link=link, colours=rainbow(3), xrange=c(1,maxSpeed), yrange=c(0,1), epsilons=epsilons)
		dev.off()
	}
}

plotAllofAll <- function(prefix, nodes=31, conservativeEps=0.03, plts=c(1,2,3,4)){
	links <- c(0.3,0.5,0.7,0.2)
	edges <- c(4,8,12,16)
	#methods <- c("belief", "conf-perfect")
	methods <- c("belief-acc", "belief-log05", "belief-05log", "conf-perfect")
	rain <- rainbow(5)
	epsilons <- seq(0.01,0.11,0.005)
	acc <- 0.6
	stddev <- 0.12
	maxSpeed <- 25
	
	optimEps <- c(0.075, 0.065, 0.06, 0.055)
	l <- 0
	for (edge in edges){
		l <- l+1
		if (1 %in% plts){
			pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/belief-acc-BA-e",edge,".pdf",sep=""))
			plotSuccessOnEpsilon(prefix=prefix, nodes=nodes, netType="barabasi-albert", accuracy=acc, acstdv=stddev,
					methods=methods, edge=edge, link=links[1], colours=rain, yrange=c(0,1), epsilons=epsilons)
			dev.off()
		}
		if (2 %in% plts){
			pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/belief-speed-BA-e",edge,".pdf",sep=""))
			plotTimeOnEpsilon(prefix=prefix, nodes=nodes, netType="barabasi-albert", accuracy=acc, acstdv=stddev, bxplt=TRUE,
					methods=methods, edge=edge, link=links[1], colours=rain, yrange=c(1,maxSpeed), epsilons=epsilons)
			dev.off()
		}
		if (3 %in% plts){
			pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/belief-trade-BA-e",edge,".pdf",sep=""))
			plotSpeedAccuracyTradeoff(prefix=prefix, nodes=nodes, netType="barabasi-albert", accuracy=acc, acstdv=stddev, addOptim=TRUE,
					methods=methods, edge=edge, link=links[1], colours=rain, xrange=c(1,maxSpeed), yrange=c(0.5,1), epsilons=epsilons,
					optimalEps=optimEps[l], conservativeEps=conservativeEps)
			dev.off()
		}
		if (4 %in% plts){
			pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/belief-both-BA-e",edge,".pdf",sep=""))
			plotSuccTimeOnEpsilon(prefix=prefix, nodes=nodes, netType="barabasi-albert", accuracy=acc, acstdv=stddev, 
					methods=methods, edge=edge, link=links[1], colours=rain, yrangel=c(0.5,1), yranger=c(0,maxSpeed), epsilons=epsilons,
					optimalEps=optimEps[l], conservativeEps=conservativeEps)
			dev.off()
		}
	}
	optimEps <- c(0.095, 0.06, 0.045, 0.11)
	l <- 0
	for (link in links){
		l<-l+1
		linkstr <- gsub("\\.", "", toString(link))
		if (1 %in% plts){
			pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/belief-acc-ER-l",linkstr,".pdf",sep=""))
			plotSuccessOnEpsilon(prefix=prefix, nodes=nodes, netType="erdos-renyi", accuracy=acc, acstdv=stddev,
					methods=methods, edge=edges[1], link=link, colours=rain, yrange=c(0,1), epsilons=epsilons)
			dev.off()
		}
		if (2 %in% plts){
			pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/belief-speed-ER-l",linkstr,".pdf",sep=""))
			plotTimeOnEpsilon(prefix=prefix, nodes=nodes, netType="erdos-renyi", accuracy=acc, acstdv=stddev, bxplt=TRUE,
					methods=methods, edge=edges[1], link=link, colours=rain, yrange=c(1,maxSpeed), epsilons=epsilons)
			dev.off()
		}
		if (3 %in% plts){
			pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/belief-trade-ER-l",linkstr,".pdf",sep=""))
			plotSpeedAccuracyTradeoff(prefix=prefix, nodes=nodes, netType="erdos-renyi", accuracy=acc, acstdv=stddev, addOptim=TRUE,
					methods=methods, edge=edges[1], link=link, colours=rain, xrange=c(1,maxSpeed), yrange=c(0.5,1), epsilons=epsilons,
					optimalEps=optimEps[l], conservativeEps=conservativeEps)
			dev.off()
		}
		if (4 %in% plts){
			pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/belief-both-ER-l",linkstr,".pdf",sep=""))
			plotSuccTimeOnEpsilon(prefix=prefix, nodes=nodes, netType="erdos-renyi", accuracy=acc, acstdv=stddev, 
					methods=methods, edge=edges[1], link=link, colours=rain, yrangel=c(0.5,1), yranger=c(0,maxSpeed), epsilons=epsilons,
					optimalEps=optimEps[l], conservativeEps=conservativeEps)
			dev.off()
		}
	}
}

parseAllFiles <- function(prefix, nodes=31, links=c(0.3,0.5,0.7,0.2), edges=c(4,8,12,16), 
		methods=c("belief-acc", "belief-log05", "belief-05log", "conf-perfect"), methods_id=c(1,2,3,4),
		epsilons=seq(0.01,0.11,0.005) ){
	#methods <- c("belief", "conf-perfect")
	acc <- 0.6
	stddev <- 0.12
	update = "belief-up"
	limitToAvoidNumericalIssues <- 30

	finalTable <- data.frame()
	
	netTypeID = 1 #E-R net
	netType="erdos-renyi"
	edge=edges[1]
	for (link in links){
		m <- 0
		for (method in methods){
			m <- m+1
			for (epsilon in epsilons){
				filename <- paste(prefix,"out_net-", netType, "_nodes-",nodes,"_edges-",edge,
						"_link-",link,"_model-",method,"_up-",update,"_acc-",acc,"_acstdv-",stddev,"_eps-",format(epsilon, nsmall=3),".txt",sep="")
				print(filename)
				data <- read.table(filename, header=T)
				eff <- nrow(data[ data$pos == nodes | data$neg == nodes, ]) / nrow(data) 
				succ <- nrow(data[ data$pos == nodes, ]) / nrow(data) 
				time <-  mean(data[ (data$pos == nodes | data$neg == nodes) & data$iter < limitToAvoidNumericalIssues, 'iter']) 
				row <- c(netTypeID,methods_id[m],link,epsilon,eff,succ,time)
				finalTable <- rbind(finalTable, row)
			}
		}
	}
	
	netTypeID = 2 #B-A net
	netType="barabasi-albert"
	link=links[1]
	for (edge in edges){
		m <- 0
		for (method in methods){
			m <- m+1
			for (epsilon in epsilons){
				filename <- paste(prefix,"out_net-", netType, "_nodes-",nodes,"_edges-",edge,
						"_link-",link,"_model-",method,"_up-",update,"_acc-",acc,"_acstdv-",stddev,"_eps-",format(epsilon, nsmall=3),".txt",sep="")
				print(filename)
				data <- read.table(filename, header=T)
				eff <- nrow(data[ data$pos == nodes | data$neg == nodes, ]) / nrow(data) 
				succ <- nrow(data[ data$pos == nodes, ]) / nrow(data) 
				time <-  mean(data[ (data$pos == nodes | data$neg == nodes) & data$iter < limitToAvoidNumericalIssues, 'iter']) 
				row <- c(netTypeID,methods_id[m],edge,epsilon,eff,succ,time)
				finalTable <- rbind(finalTable, row)
			}
		}
	}
	
	write.table(finalTable, paste(prefix,"/finalTable.txt",sep=""), sep="\t", col.names=FALSE, row.names=FALSE)
	return(finalTable)
}

plotCompareInitConfFunctions <- function(){
	pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/compareInitConfFuncs.pdf",sep=""))
	
	linewidth=4
	par(mar = c(4, 5.5, 0, 0) + 0.3)  # Leave space for z axis
	par(pty="s")
	plot(c(0,1),c(0,2), type='n', xlab=expression(alpha[i]), ylab=expression(c[i]^{t[0]} ), cex.lab=2, cex.axis=1.5 )
	x <- seq(0,1,0.01)
	lines(x, log(x/(1-x)), col="red", lwd=linewidth, lty=1)
	lines(x, log(2*x)/log(2), col="blue", lwd=linewidth, lty=2)
	legend('topleft',c(expression(log(frac(alpha[i],1-alpha[1]))), expression(frac(log(2*alpha[i]), log(2)))), col=c('red','blue'), lty=c(1,2), lwd=linewidth, cex=1.5, bty = "n" )
	
	dev.off() 
}

plotSuccessOnEpsilon <- function(prefix, nodes=31, netType="barabasi-albert", accuracy=0.6, acstdv=0.12,
		methods=c("belief", "conf-perfect"), edge=4, link=0.3, colours=rainbow(2), yrange=c(0,1),
		epsilons=seq(0.01,0.12,0.02) ) {
	pntTypes = c(1,4,5,3,8,7,0,2,6,9,10,11,12,13)
	
	xmin <- min(epsilons)
	xmax <- max(epsilons)
	xrange <- (xmax - xmin)
	plot(c(xmin,xmax), yrange, type='n', xlab="Belief Epsilon", ylab="Success rate", xaxt = "n")
	axis(1, at=epsilons)
	grid()
	if (length(methods)>length(colours)){
		print(paste("WARNING! - There are ", length(methods), " methods to display and you're using only ", length(colours), " colours.", sep=""))
	}
	update = "belief-up"
	
	legTxt <- c()
	m <- 0
	for (method in methods){
		m <- m+1
		dataFilt <- c()
		for (epsilon in epsilons){
			filename <- paste(prefix,"out_net-", netType, "_nodes-",nodes,"_edges-",edge,
					"_link-",link,"_model-",method,"_up-",update,"_acc-",accuracy,"_acstdv-",acstdv,"_eps-",format(epsilon, nsmall=3),".txt",sep="")
			print(filename)
			data <- read.table(filename, header=T)
			dataFilt <- c(dataFilt, nrow(data[ data$pos == nodes, ]) / nrow(data) )
			if (epsilon == epsilons[1]) {
				ltxt <- paste(method,update,sep=" ")
				ltxt <- gsub("perfect", "weigh", ltxt)
				ltxt <- gsub("rand", "rule", ltxt)
				ltxt <- gsub("optim-up", "with-up", ltxt)
				legTxt <- c(legTxt, ltxt)
			}
		}
		
		points( epsilons, dataFilt, pch=pntTypes[m], col=colours[m], cex=1, lwd=2, type='b')
		print(dataFilt)
	}
	legend('bottomright', legTxt, pch=pntTypes[1:length(methods)], col=colours[1:length(methods)], lty=0, cex=1, lwd=2, bg='white' )
}

 

plotSuccTimeOnEpsilon <- function(prefix, nodes=31, netType="barabasi-albert", accuracy=0.6, acstdv=0.12,
		methods=c("belief", "conf-perfect"), edge=4, link=0.3, colours=rainbow(2), yrangel=c(0,1), yranger=c(0,30),
		epsilons=seq(0.01,0.12,0.02), optimalEps=0.065, conservativeEps=0.03 ) {
	pntTypes = c(1,4,5,3,8,7,0,2,6,9,10,11,12,13)
	limitToAvoidNumericalIssues <- 30
	par(mar = c(5, 4, 0, 4) + 0.3)  # Leave space for z axis
	
	xmin <- min(epsilons)
	xmax <- max(epsilons)
	xrange <- (xmax - xmin)
	plot(c(xmin,xmax), yrangel, type='n', xlab="Belief Epsilon", ylab="Success rate", xaxt = "n")
	axis(1, at=epsilons)
	grid()
	lines(c(optimalEps,optimalEps), c(min(yrangel)-0.1,max(yrangel)+0.1), lty=2, lwd=2)
	lines(c(conservativeEps,conservativeEps), c(min(yrangel)-0.1,max(yrangel)+0.1), lty=2, lwd=2)
	text(	optimalEps+0.001* if(optimalEps>conservativeEps){1}else{-1}, min(yrangel)+(max(yrangel)-min(yrangel))*0.02,
			substitute( epsilon == eps, list(eps=optimalEps) ), cex=1, adj = c(if(optimalEps>conservativeEps){0}else{1}, NA)
	)
	text(	conservativeEps+0.001* if(optimalEps>conservativeEps){-1}else{1}, min(yrangel)+(max(yrangel)-min(yrangel))*0.02,
			substitute( epsilon == eps, list(eps=conservativeEps) ), cex=1, adj = c(if(optimalEps>conservativeEps){1}else{0}, NA)
	)
#	text(	optimalEps+0.001,
#			min(yrangel)+(max(yrangel)-min(yrangel))*0.02,
#			substitute(
##					paste( epsilon == oeps %->% paste(integral(P(k,G[list(n==nu,p[e]==pe)]),0,1/epsilon),italic("d"),k%~~%0.72)),
#					paste( epsilon == oeps %->% paste(integral(P(k,G),0,1/epsilon))[paste(italic("d"),k)] %~~% 0.72),
#					list(nu=nodes, pe = link, mnum=round(1/optimalEps, digits=1), oeps=optimalEps)
#			),
#			cex=1,
#			adj = c( if(optimalEps < (xmin + (xrange/2))){0}else{1} , NA)
#	)
	if (length(methods)>length(colours)){
		print(paste("WARNING! - There are ", length(methods), " methods to display and you're using only ", length(colours), " colours.", sep=""))
	}
	update = "belief-up"
	
	legTxt <- c()
	m <- 0
	for (method in methods){
		m <- m+1
		effVect <- c()
		succVect <- c()
		timeVect <- c()
		for (epsilon in epsilons){
			filename <- paste(prefix,"out_net-", netType, "_nodes-",nodes,"_edges-",edge,
					"_link-",link,"_model-",method,"_up-",update,"_acc-",accuracy,"_acstdv-",acstdv,"_eps-",format(epsilon, nsmall=3),".txt",sep="")
			print(filename)
			data <- read.table(filename, header=T)
			effVect <- c(effVect, nrow(data[ data$pos == nodes | data$neg == nodes, ]) / nrow(data) )
			succVect <- c(succVect, nrow(data[ data$pos == nodes, ]) / nrow(data) )
			timeVect <- c(timeVect, mean(data[ (data$pos == nodes | data$neg == nodes) & data$iter < limitToAvoidNumericalIssues, 'iter']) )
			if (epsilon == epsilons[1]) {
				ltxt <- paste(method,update,sep=" ")
				ltxt <- gsub("perfect", "weigh", ltxt)
				ltxt <- gsub("rand", "rule", ltxt)
				ltxt <- gsub("optim-up", "with-up", ltxt)
				legTxt <- c(legTxt, ltxt)
			}
		}
		
#		points( epsilons, succVect, pch=pntTypes[m], col = alpha(colours[m], 0.5), cex=1, lwd=2, type='b')
		points( epsilons, succVect, pch=pntTypes[m], col=alpha(colours[m], effVect), cex=1, lwd=2, type='p')
		for (p in seq(1,length(succVect)-1)){
			points( epsilons[p:(p+1)], succVect[p:(p+1)], col=alpha(colours[m], effVect[(p+1)]) , type='b', pch="", lwd=2)
		}
		points( epsilons[which( almost.equal(epsilons,optimalEps) )], succVect[which( almost.equal(epsilons,optimalEps) )], pch=1, lwd=2, col=alpha(colours[m], 1), cex=2, type='p')
		points( epsilons[which( almost.equal(epsilons,conservativeEps) )], succVect[which( almost.equal(epsilons,conservativeEps) )], pch=0, lwd=2, col=alpha(colours[m], 1), cex=2, type='p')
		print(succVect)
	}
	
	## Right axis plot
	par(new=T)
	plot(c(xmin,xmax), yranger, type='n', xaxt = "n", axes = FALSE, bty = "n", xlab = "", ylab = "")
	axis(side=4, at = pretty(yranger) )
	mtext("Mean RT", side=4, line=3)
	m <- 0
	for (method in methods){
		m <- m+1
		effVect <- c()
		timeVect <- c()
		for (epsilon in epsilons){
			filename <- paste(prefix,"out_net-", netType, "_nodes-",nodes,"_edges-",edge,
					"_link-",link,"_model-",method,"_up-",update,"_acc-",accuracy,"_acstdv-",acstdv,"_eps-",format(epsilon, nsmall=3),".txt",sep="")
			print(filename)
			data <- read.table(filename, header=T)
			effVect <- c(effVect, nrow(data[ data$pos == nodes | data$neg == nodes, ]) / nrow(data) )
			timeVect <- c(timeVect, mean(data[ (data$pos == nodes | data$neg == nodes) & data$iter < limitToAvoidNumericalIssues, 'iter']) )
		}
		
		points( epsilons, timeVect, pch=pntTypes[m], col=alpha(colours[m], effVect), cex=1, lwd=2, type='p')
		for (p in seq(1,length(timeVect)-1)){
			points( epsilons[p:(p+1)], timeVect[p:(p+1)], col=alpha(colours[m], effVect[(p+1)]) , type='b', pch="", lwd=2, lty=3)
		}
		points( epsilons[which( almost.equal(epsilons,optimalEps) )], timeVect[which( almost.equal(epsilons,optimalEps) )], pch=1, lwd=2, col=alpha(colours[m], 1), cex=2, type='p')
		points( epsilons[which( almost.equal(epsilons,conservativeEps) )], timeVect[which( almost.equal(epsilons,conservativeEps) )], pch=0, lwd=2, col=alpha(colours[m], 1), cex=2, type='p')
		print(timeVect)
	}
	
	legend('right', legTxt, pch=pntTypes[1:length(methods)], col=colours[1:length(methods)], lty=0, cex=1, lwd=2, bg='white' )
}


plotTimeOnEpsilon <- function(prefix, nodes=31, netType="barabasi-albert", accuracy=0.6, acstdv=0.12,
		methods=c("belief", "conf-perfect"), edge=4, link=0.3, colours=rainbow(2), yrange=c(1,20), bxplt=FALSE,
		epsilons=seq(0.01,0.12,0.02) ) {
	pntTypes = c(1,4,5,3,8,7,0,2,6,9,10,11,12,13)
	limitToAvoidNumericalIssues <- 30
	
	xmin <- min(epsilons)
	xmax <- max(epsilons)
	xrange <- (xmax - xmin)
	plot(c(xmin,xmax), yrange, type='n', xlab="Belief Epsilon", ylab="Mean RT", xaxt = "n")
	axis(1, at=epsilons)
	grid()
	if (length(methods)>length(colours)){
		print(paste("WARNING! - There are ", length(methods), " methods to display and you're using only ", length(colours), " colours.", sep=""))
	}
	update = "belief-up"
	
	legTxt <- c()
	m <- 0
	for (method in methods){
		m <- m+1
		effVect <- c()
		dataFilt <- c()
		for (epsilon in epsilons){
			filename <- paste(prefix,"out_net-", netType, "_nodes-",nodes,"_edges-",edge,
					"_link-",link,"_model-",method,"_up-",update,"_acc-",accuracy,"_acstdv-",acstdv,"_eps-",format(epsilon, nsmall=3),".txt",sep="")
			print(filename)
			data <- read.table(filename, header=T)
			if (bxplt){
				dataFilt <- append( dataFilt, list(data[ (data$pos == nodes | data$neg == nodes) & data$iter < limitToAvoidNumericalIssues, 'iter']) )
			} else {
				effVect <- c(effVect, nrow(data[ data$pos == nodes | data$neg == nodes, ]) / nrow(data) )
				dataFilt <- c(dataFilt, mean(data[ (data$pos == nodes | data$neg == nodes) & data$iter < limitToAvoidNumericalIssues, 'iter']) )
			}
			if (epsilon == epsilons[1]) {
				ltxt <- paste(method,update,sep=" ")
				ltxt <- gsub("perfect", "weigh", ltxt)
				ltxt <- gsub("rand", "rule", ltxt)
				ltxt <- gsub("optim-up", "with-up", ltxt)
				legTxt <- c(legTxt, ltxt)
			}
		}
		if (bxplt){
			#print(length(dataFilt))
			#print(positions)
			size <- xrange/70
			positions = epsilons - (m-1)*size*1.5
			boxplot( dataFilt, at=positions, boxwex=size, add=T, col=colours[m], axes=F, outcex=0.5)
		} else {
			text( epsilons, dataFilt+(max(yrange)/15), labels=effVect, col=colours[m])
			points( epsilons, dataFilt, pch=pntTypes[m], col=colours[m], cex=1, lwd=2, type='b' )
			print(dataFilt)
		}
		
	}
	if (bxplt){
		legend('topright', legTxt, fill=colours[1:length(methods)], cex=1, bg='white')
	} else {
		legend('topright', legTxt, pch=pntTypes[1:length(methods)], col=colours[1:length(methods)], lty=0, cex=1, lwd=2, bg='white' )
	}
}


plotSpeedAccuracyTradeoff <- function(prefix, nodes=31, netType="barabasi-albert", accuracy=0.6, acstdv=0.12,
		methods=c("belief", "conf-perfect"), edge=4, link=0.3, colours=rainbow(3), xrange=c(1,20), yrange=c(0,1), bxplt=FALSE,
		epsilons=seq(0.01,0.12,0.02), addOptim=TRUE, optimalEps=0.065, conservativeEps=0.03 ) {
	
	pntTypes = c(1,4,5,3,8,7,0,2,6,9,10,11,12,13)
	limitToAvoidNumericalIssues <- 30
	
	plot(xrange, yrange, type='n', ylab="Accuracy", xlab="Mean RT")
	grid()
	if (length(methods)+addOptim>length(colours)){
		print(paste("WARNING! - There are ", length(methods), " methods to display and you're using only ", length(colours), " colours.", sep=""))
	}
	update = "belief-up"
	
	legTxt <- c()
	m <- 0
	for (method in methods){
		m <- m+1
		effVect <- c()
		speedVect <- c()
		accVect <- c()
		for (epsilon in epsilons){
			filename <- paste(prefix,"out_net-", netType, "_nodes-",nodes,"_edges-",edge,
					"_link-",link,"_model-",method,"_up-",update,"_acc-",accuracy,"_acstdv-",acstdv,"_eps-",format(epsilon, nsmall=3),".txt",sep="")
			print(filename)
			data <- read.table(filename, header=T)
			if (bxplt){
				speedVect <- append( speedVect, list(data[ (data$pos == nodes | data$neg == nodes) & data$iter < limitToAvoidNumericalIssues, 'iter']) )
				accVect <- c(accVect, nrow(data[ data$pos == nodes, ]) / nrow(data) )
			} else {
				effVect <- c(effVect, nrow(data[ data$pos == nodes | data$neg == nodes, ]) / nrow(data) )
#				if (effVect[length(effVect)] > 0.3 ){
					speedVect <- c(speedVect, mean(data[ (data$pos == nodes | data$neg == nodes) & data$iter < limitToAvoidNumericalIssues, 'iter']) )
					accVect <- c(accVect, nrow(data[ data$pos == nodes, ]) / nrow(data) )
#				} else {
#					print("Never converging, point is skipped")
#				}
			}
			if (epsilon == epsilons[1]) {
				ltxt <- paste(method,update,sep=" ")
				ltxt <- gsub("perfect", "weigh", ltxt)
				ltxt <- gsub("rand", "rule", ltxt)
				ltxt <- gsub("optim-up", "with-up", ltxt)
				legTxt <- c(legTxt, ltxt)
			}
		}
		if (bxplt){
			#print(length(dataFilt))
			#print(positions)
			size <- xrange/70
			positions = epsilons - (m-1)*size*1.5
			boxplot( speedVect, at=positions, boxwex=size, add=T, col=colours[m], axes=F, outcex=0.5)
		} else {
			points( speedVect, accVect, pch=pntTypes[m], col=alpha(colours[m], effVect), cex=1, lwd=2, type='p')
			for (p in seq(1,length(speedVect)-1)){
				points( speedVect[p:(p+1)], accVect[p:(p+1)], col=alpha(colours[m], effVect[(p+1)]) , type='b', pch="", lwd=2, lty=3)
			}
			points( speedVect[which( almost.equal(epsilons,optimalEps) )], accVect[which( almost.equal(epsilons,optimalEps) )], pch=1, lwd=2, col=alpha(colours[m], 1), cex=2, type='p')
			points( speedVect[which( almost.equal(epsilons,conservativeEps) )], accVect[which( almost.equal(epsilons,conservativeEps) )], pch=0, lwd=2, col=alpha(colours[m], 1), cex=2, type='p')
						
#			points( speedVect, accVect, pch=pntTypes[m], col=colours[m], cex=1, lwd=2, type='b' )
			print(effVect)
			print(accVect)
			print(speedVect)
		}
		
	}
	
	if (addOptim){
		filename <- paste(prefix,"out_net-", netType, "_nodes-",nodes,"_edges-",edge,
				"_link-",link,"_model-conf-perfect_up-optim-up_acc-",accuracy,"_acstdv-",acstdv,"_eps-",format(epsilons[1], nsmall=3),".txt",sep="")
		print(filename)
		data <- read.table(filename, header=T)
		speed <- mean(data[ data$pos == nodes | data$neg == nodes, 'iter'])
		acc <- nrow(data[ data$pos == nodes, ]) / nrow(data)
		points( speed, acc, pch=pntTypes[length(methods)+1], col=colours[length(methods)+1], cex=1, lwd=2, type='b' )
		abline(v=speed, col=colours[length(methods)+1], lty=2, lwd=1)
		abline(h=acc, col=colours[length(methods)+1], lty=2, lwd=1)
		print(acc)
		print(speed)
	}
	
	if (bxplt){
		legend('bottomright', legTxt, fill=colours[1:length(methods)], cex=1, bg='white')
	} else {
		legend('bottomright', legTxt, pch=pntTypes[1:length(methods)], col=colours[1:length(methods)], lty=0, cex=1, lwd=2, bg='white' )
	}
	
	

}

