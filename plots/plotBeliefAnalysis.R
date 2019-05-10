library(scales)

almost.equal <- function (x, y, tolerance=.Machine$double.eps^0.5, 
		na.value=TRUE)  { 
	answer <- rep(na.value, length(x)) 
	test <- !is.na(x) 
	answer[test] <- abs(x[test] - y) < tolerance 
	answer 
}

plotAllofAll <- function(prefix, parse=TRUE, conservativeEps=0.02, plts=c(1,2), T_MAX=30){
	links <- c(0.2,0.4,0.6,0.8)
	edges <- c(3,7,11,15)
	nodesList <- seq(11,47,12)
	
	#methods <- c("belief", "conf-perfect")
	methods <- c("belief-init", "conf-perfect")
	rain <- rainbow(3)
	epsilons <- seq(0.01,0.15,0.005)
	acc <- 0.6
	stddev <- 0.12
	maxSpeedOnPlot <- 25
	
	if (parse){
		fullTable <- parseAllFiles(prefix, nodeList=nodesList, links=links, edges=edges, 
				methods=methods, methods_id=seq(1,length(methods)), epsilons=epsilons, T_MAX=T_MAX )
	} else {
		fullTable <- read.table(paste(prefix,"/finalTable.txt",sep=""), header=TRUE)
	}
	
	for (nodes in nodesList){
		optimEps <- c(0.075, 0.065, 0.06, 0.055)
		l <- 0
		for (edge in edges){
			l <- l+1
			if (edge >= nodes){ next }
			if (1 %in% plts){
				pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/belief-trade-n",nodes,"-BA-e",edge,".pdf",sep=""))
				plotSpeedAccuracyTradeoff(fullTable=fullTable, prefix=prefix, nodes=nodes, netType="barabasi-albert", accuracy=acc, acstdv=stddev, addOptim=TRUE,
						methods=methods, edge=edge, link=links[1], colours=rain, xrange=c(1,maxSpeedOnPlot), yrange=c(0.5,1), epsilons=epsilons,
						optimalEps=optimEps[l], conservativeEps=conservativeEps, T_MAX=T_MAX)
				dev.off()
			}
			if (2 %in% plts){
				pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/belief-both-n",nodes,"-BA-e",edge,".pdf",sep=""))
				plotSuccTimeOnEpsilon(fullTable=fullTable, prefix='none', nodes=nodes, netType="barabasi-albert", accuracy=acc, acstdv=stddev, 
						methods=methods, edge=edge, link=links[1], colours=rain, yrangel=c(0.5,1), yranger=c(0,maxSpeedOnPlot), epsilons=epsilons,
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
				pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/belief-trade-n",nodes,"-ER-l",linkstr,".pdf",sep=""))
				plotSpeedAccuracyTradeoff(fullTable=fullTable, prefix=prefix, nodes=nodes, netType="erdos-renyi", accuracy=acc, acstdv=stddev, addOptim=TRUE,
						methods=methods, edge=edges[1], link=link, colours=rain, xrange=c(1,maxSpeedOnPlot), yrange=c(0.5,1), epsilons=epsilons,
						optimalEps=optimEps[l], conservativeEps=conservativeEps, T_MAX=T_MAX)
				dev.off()
			}
			if (2 %in% plts){
				pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/belief-both-n",nodes,"-ER-l",linkstr,".pdf",sep=""))
				plotSuccTimeOnEpsilon(fullTable=fullTable, prefix=prefix, nodes=nodes, netType="erdos-renyi", accuracy=acc, acstdv=stddev, 
						methods=methods, edge=edges[1], link=link, colours=rain, yrangel=c(0.5,1), yranger=c(0,maxSpeedOnPlot), epsilons=epsilons,
						optimalEps=optimEps[l], conservativeEps=conservativeEps)
				dev.off()
			}
		}
	}
}

parseAllFiles <- function(prefix, nodeList=c(11,47), links=c(0.3,0.5,0.7,0.2), edges=c(4,8,12,16), 
		methods=c("belief-acc", "belief-log05", "belief-05log", "conf-perfect"), methods_id=c(1,2,3,4),
		epsilons=seq(0.01,0.11,0.005), T_MAX=30 ){
	#methods <- c("belief", "conf-perfect")
	acc <- 0.6
	stddev <- 0.12
	update = "belief-up"
	updateID = 2

	finalTable <- data.frame()
	
	for (nodes in nodeList){
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
					eff <- nrow(data[ data$iter <= T_MAX & (data$pos == nodes | data$neg == nodes), ]) / nrow(data) 
					succ <- nrow(data[ data$iter <= T_MAX & data$pos == nodes, ]) / nrow(data) 
					time <-  mean(data[ (data$pos == nodes | data$neg == nodes) & data$iter < T_MAX, 'iter']) 
					row <- c(nodes, netTypeID, methods_id[m], updateID, link, epsilon, eff, succ, time)
					finalTable <- rbind(finalTable, row)
				}
			}
		}
		
		netTypeID = 2 #B-A net
		netType="barabasi-albert"
		link=links[1]
		for (edge in edges){
			if (edge >= nodes){ next }
			m <- 0
			for (method in methods){
				m <- m+1
				for (epsilon in epsilons){
					filename <- paste(prefix,"out_net-", netType, "_nodes-",nodes,"_edges-",edge,
							"_link-",link,"_model-",method,"_up-",update,"_acc-",acc,"_acstdv-",stddev,"_eps-",format(epsilon, nsmall=3),".txt",sep="")
					print(filename)
					data <- read.table(filename, header=T)
					eff <- nrow(data[ data$iter <= T_MAX & (data$pos == nodes | data$neg == nodes), ]) / nrow(data) 
					succ <- nrow(data[ data$iter <= T_MAX & data$pos == nodes, ]) / nrow(data) 
					time <-  mean(data[ (data$pos == nodes | data$neg == nodes) & data$iter <= T_MAX, 'iter']) 
					row <- c(nodes, netTypeID, methods_id[m], updateID, edge, epsilon, eff, succ, time)
					finalTable <- rbind(finalTable, row)
				}
			}
		}
	}
	
	colnames(finalTable) <- c("nodes", "netTypeID", "methodID", "updateID", "netParam", "epsilon", "eff", "succ", "time")
	write.table(finalTable, paste(prefix,"/finalTable.txt",sep=""), sep="\t", col.names=TRUE, row.names=FALSE)
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

plotSuccTimeOnEpsilon <- function(fullTable=NULL, prefix='none', nodes=31, netType="barabasi-albert", accuracy=0.6, acstdv=0.12,
		methods=c("belief", "conf-perfect"), edge=4, link=0.3, colours=rainbow(2), yrangel=c(0,1), yranger=c(0,30),
		epsilons=seq(0.01,0.12,0.02), optimalEps=0.065, conservativeEps=0.03) {
	
	pntTypes = c(1,4,5,3,8,7,0,2,6,9,10,11,12,13)
	par(mar = c(5, 4, 0, 4) + 0.3)  # Leave space for z axis
	xmin <- min(epsilons)
	xmax <- max(epsilons)
	xrange <- (xmax - xmin)
	plot(c(xmin,xmax), yrangel, type='n', xlab="Belief Epsilon", ylab="Success rate", xaxt = "n")
	axis(1, at=epsilons)
	grid()
	
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
	
	if (is.null(fullTable)){
		fullTable <- parseAllFiles(prefix, nodeList=c(nodes), links=links, edges=edges, 
				methods=methods, methods_id=seq(1,length(methods)), epsilons=epsilons )
	}
	
	if (netType == "erdos-renyi"){ netID <- 1; netParam <- link }
	if (netType == "barabasi-albert"){ netID <- 2; netParam <- edge }
	
	subVect <- fullTable[ fullTable['nodes']==nodes & fullTable['netTypeID']==netID & fullTable['methodID']==1 & 
					fullTable['updateID']==2 & fullTable['netParam']==netParam, c("epsilon", "time", "succ", "eff")]
	conservativeEps <- max(subVect[1/(nodes-1)-subVect[,'epsilon']>0, 'epsilon'])
	
	subVect <- subVect[ subVect['succ'] == max(subVect[,'succ']), ]
	subVect <- subVect[ subVect['time'] == min(subVect[,'time']), ]
	optimalEps <- subVect[,'epsilon']
	
	lines(c(optimalEps,optimalEps), c(min(yrangel)-0.1,max(yrangel)+0.1), lty=2, lwd=2)
	text(	optimalEps+0.001* if(optimalEps>conservativeEps){1}else{-1}, min(yrangel)+(max(yrangel)-min(yrangel))*0.02,
			substitute( epsilon == eps, list(eps=optimalEps) ), cex=1, adj = c(if(optimalEps>conservativeEps){0}else{1}, NA)
	)
	lines(c(conservativeEps,conservativeEps), c(min(yrangel)-0.1,max(yrangel)+0.1), lty=2, lwd=2)
	text(	conservativeEps+0.001* if(optimalEps>conservativeEps){-1}else{1}, min(yrangel)+(max(yrangel)-min(yrangel))*0.02,
			substitute( epsilon == eps, list(eps=conservativeEps) ), cex=1, adj = c(if(optimalEps>conservativeEps){1}else{0}, NA)
	)
	
	
	legTxt <- c()
	m <- 0
	for (method in methods){
		m <- m+1
		ltxt <- paste(method,update,sep=" ")
		ltxt <- gsub("perfect", "weigh", ltxt)
		ltxt <- gsub("rand", "rule", ltxt)
		ltxt <- gsub("optim-up", "with-up", ltxt)
		legTxt <- c(legTxt, ltxt)
		
		succVect <- fullTable[ fullTable['nodes']==nodes & fullTable['netTypeID']==netID & fullTable['methodID']==m & 
						fullTable['updateID']==2 & fullTable['netParam']==netParam, c("epsilon", "succ", "eff")]
		points( succVect[,1], succVect[,2], pch=pntTypes[m], col=alpha(colours[m], succVect[,3]), cex=1, lwd=2, type='p')
		for (p in seq(1,length(succVect[,1])-1)){
			points( succVect[p:(p+1),1], succVect[p:(p+1),2], col=alpha(colours[m], succVect[(p+1),3]) , type='b', pch="", lwd=2)
		}
		points( succVect[which( almost.equal(succVect[,1],optimalEps) ),1], succVect[which( almost.equal(succVect[,1],optimalEps) ),2], pch=1, lwd=2, col=alpha(colours[m], 1), cex=2, type='p')
		#points( subVect['epsilon'], subVect['succ'], pch=1, lwd=2, col=alpha(colours[m], 1), cex=2, type='p')
		points( succVect[which( almost.equal(succVect[,1],conservativeEps) ),1], succVect[which( almost.equal(succVect[,1],conservativeEps) ),2], pch=0, lwd=2, col=alpha(colours[m], 1), cex=2, type='p')
		#print(succVect)
	}
	
	## Right axis plot
	par(new=T)
	plot(c(xmin,xmax), yranger, type='n', xaxt = "n", axes = FALSE, bty = "n", xlab = "", ylab = "")
	axis(side=4, at = pretty(yranger) )
	mtext("Mean RT", side=4, line=3)
	m <- 0
	for (method in methods){
		m <- m+1
		
		timeVect <- fullTable[ fullTable['nodes']==nodes & fullTable['netTypeID']==netID & fullTable['methodID']==m & 
						fullTable['updateID']==2 & fullTable['netParam']==netParam, c("epsilon", "time", "eff")]
		points( timeVect[,1], timeVect[,2], pch=pntTypes[m], col=alpha(colours[m], timeVect[,3]), cex=1, lwd=2, type='p')
		for (p in seq(1,length(timeVect[,1])-1)){
			points( timeVect[p:(p+1),1], timeVect[p:(p+1),2], col=alpha(colours[m], timeVect[(p+1),3]) , type='b', pch="", lwd=2, lty=3)
		}
		points( timeVect[which( almost.equal(timeVect[,1],optimalEps) ),1], timeVect[which( almost.equal(timeVect[,1],optimalEps) ),2], pch=1, lwd=2, col=alpha(colours[m], 1), cex=2, type='p')
		#points( subVect['epsilon'], subVect['time'], pch=1, lwd=2, col=alpha(colours[m], 1), cex=2, type='p')
		points( timeVect[which( almost.equal(timeVect[,1],conservativeEps) ),1], timeVect[which( almost.equal(timeVect[,1],conservativeEps) ),2], pch=0, lwd=2, col=alpha(colours[m], 1), cex=2, type='p')
		#print(timeVect)
	}
	
	legend('right', legTxt, pch=pntTypes[1:length(methods)], col=colours[1:length(methods)], lty=0, cex=1, lwd=2, bg='white' )
}

plotSpeedAccuracyTradeoff <- function(fullTable=NULL, prefix='none', nodes=31, netType="barabasi-albert", accuracy=0.6, acstdv=0.12,
		methods=c("belief", "conf-perfect"), edge=4, link=0.3, colours=rainbow(3), xrange=c(1,20), yrange=c(0,1), bxplt=FALSE,
		epsilons=seq(0.01,0.12,0.02), addOptim=TRUE, optimalEps=0.065, conservativeEps=0.03, T_MAX=30) {
	
	pntTypes = c(1,4,5,3,8,7,0,2,6,9,10,11,12,13)
	
	plot(xrange, yrange, type='n', ylab="Accuracy", xlab="Mean RT")
	grid()
	if (length(methods)+addOptim>length(colours)){
		print(paste("WARNING! - There are ", length(methods), " methods to display and you're using only ", length(colours), " colours.", sep=""))
	}
	update = "belief-up"

	if (is.null(fullTable)){
		fullTable <- parseAllFiles(prefix, nodes=nodes, links=links, edges=edges, 
				methods=methods, methods_id=seq(1,length(methods)), epsilons=epsilons )
	}
	if (netType == "erdos-renyi"){ netID <- 1; netParam <- link }
	if (netType == "barabasi-albert"){ netID <- 2; netParam <- edge }
	
	subVect <- fullTable[ fullTable['nodes']==nodes & fullTable['netTypeID']==netID & fullTable['methodID']==1 & 
					fullTable['updateID']==2 & fullTable['netParam']==netParam, c("epsilon", "time", "succ", "eff")]
	conservativeEps <- max(subVect[1/(nodes-1)-subVect[,'epsilon']>0, 'epsilon'])
	
	subVect <- subVect[ subVect['succ'] == max(subVect[,'succ']), ]
	subVect <- subVect[ subVect['time'] == min(subVect[,'time']), ]
	optimalEps <- subVect[,'epsilon']
	
	legTxt <- c()
	m <- 0
	for (method in methods){
		m <- m+1
		ltxt <- paste(method,update,sep=" ")
		ltxt <- gsub("perfect", "weigh", ltxt)
		ltxt <- gsub("rand", "rule", ltxt)
		ltxt <- gsub("optim-up", "with-up", ltxt)
		legTxt <- c(legTxt, ltxt)
		
		subVect <- fullTable[ fullTable['nodes']==nodes & fullTable['netTypeID']==netID & fullTable['methodID']==m & 
						fullTable['updateID']==2 & fullTable['netParam']==netParam, c("epsilon", "succ", "time", "eff")]
		points( subVect[,"time"], subVect[,"succ"], pch=pntTypes[m], col=alpha(colours[m], subVect[,"eff"]), cex=1, lwd=2, type='p')
		for (p in seq(1,length(subVect[,1])-1)){
			points( subVect[p:(p+1),"time"], subVect[p:(p+1),"succ"], col=alpha(colours[m], subVect[(p+1),"eff"]) , type='b', pch="", lwd=2, lty=3)
		}
		points( subVect[which( almost.equal(subVect[,"epsilon"],optimalEps) ), "time"], subVect[which( almost.equal(subVect[,"epsilon"],optimalEps) ),"succ"], pch=1, lwd=2, col=alpha(colours[m], 1), cex=2, type='p')
		points( subVect[which( almost.equal(subVect[,"epsilon"],conservativeEps) ), "time"], subVect[which( almost.equal(subVect[,"epsilon"],conservativeEps) ),"succ"], pch=0, lwd=2, col=alpha(colours[m], 1), cex=2, type='p')
				
	}
	
	if (addOptim){
		filename <- paste(prefix,"out_net-", netType, "_nodes-",nodes,"_edges-",edge,
				"_link-",link,"_model-conf-perfect_up-optim-up_acc-",accuracy,"_acstdv-",acstdv,"_eps-",format(epsilons[1], nsmall=3),".txt",sep="")
		print(filename)
		data <- read.table(filename, header=T)
		speed <- mean(data[ data$iter <= T_MAX & (data$pos == nodes | data$neg == nodes), 'iter'])
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

