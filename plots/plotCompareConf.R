library(data.table)
#library(MASS)
library(lattice)
require(grid)

#############################################################################
################################ AGGREGATES #################################
#############################################################################

simpleVsDDMtoPdf <- function(prefix){
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/succ-simple-vs-ddm.pdf")
	plotSuccessOnNodes(prefix, nodes_list=seq(11, 31, 4), numRuns=100,
			netType="full", accuracy=0.6, acstdv=0.3, agentTypes=c('simple', 'DDM'),
			methods=c("conf-perfect"), colours=rainbow(10), yrange=c(0.8,1))
	dev.off()
}

allDecModelsToPdf <- function(prefix){
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/succ-dec-models.pdf")
	plotSuccessOnNodes(prefix, nodes_list=seq(11, 31, 4), numRuns=100,
			netType="full", accuracy=0.6, acstdv=0.3, agentTypes=c('DDM'),
			methods=c("best-acc", "majority-rand", "conf-perfect", "log-odds"),
			colours=rainbow(10), yrange=c(0,1))
	dev.off()
}

allSuccessMethodPlots <- function(prefix){
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/succ-onbest-ER-maj-rand.pdf")
	diffMatrixSuccessMethod(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.1, 0.6, 0.1), 
		numRuns=100, netType="erdos-renyi", accuracy=0.7, acstdv=0.3, brks=seq(-0.20,0.20,0.005),
		methods=c("majority-rand", "best-acc"))
	dev.off()
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/succ-onbest-ER-conf.pdf")
	diffMatrixSuccessMethod(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.1, 0.6, 0.1), 
		numRuns=100, netType="erdos-renyi", accuracy=0.7, acstdv=0.3, brks=seq(-0.20,0.20,0.005),
		methods=c("conf-perfect", "best-acc"))
	dev.off()

	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/succ-onbest-BA-maj-rand.pdf")
	diffMatrixSuccessMethod(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.1, 0.6, 0.1), 
		numRuns=100, netType="barabasi-albert", accuracy=0.7, acstdv=0.3, brks=seq(-0.20,0.20,0.005),
		methods=c("majority-rand", "best-acc"))
	dev.off()
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/succ-onbest-BA-conf.pdf")
	diffMatrixSuccessMethod(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.1, 0.6, 0.1), 
		numRuns=100, netType="barabasi-albert", accuracy=0.7, acstdv=0.3, brks=seq(-0.20,0.20,0.005),
		methods=c("conf-perfect", "best-acc"))
	dev.off()
}

allTimePlots <- function(prefix){
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/time-nets-maj-rand.pdf", width=17, height=7)
	plotTimeOnNodes(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.2, 0.6, 0.1), numRuns=100,
			netTypes=c("erdos-renyi", "barabasi-albert"), accuracy=0.7, acstdv=0.3,
			methods=c("majority-rand"), colours=rainbow(12))
	dev.off()
	
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/time-nets-conf.pdf", width=17, height=7)
	plotTimeOnNodes(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.2, 0.6, 0.1), numRuns=100,
			netTypes=c("erdos-renyi", "barabasi-albert"), accuracy=0.7, acstdv=0.3,
			methods=c("conf-perfect"), colours=rainbow(12))
	dev.off()
}

#############################################################################
################################# SUCCESS ###################################
#############################################################################

plotSuccessOnNodes <- function(prefix, nodes_list=seq(11, 31, 4), numRuns=100,
		netType="full", accuracy=0.6, acstdv=0.3, agentTypes=c('simple', 'DDM'),
		methods=c("best-acc", "majority-rand", "conf-perfect"), colours=rainbow(10), yrange=c(0,1)){
	combolength <- length(methods) + (length(agentTypes) > 1)
	offset <- 0.02*combolength
	xmin <- min(nodes_list)
	xmax <- max(nodes_list)
	xrange <- (xmax - xmin)
	plot(c(xmin-(xrange*offset),xmax+(xrange*offset)),yrange, type='n', xlab="Number of nodes", ylab="Success rate", xaxt = "n")
	axis(1, at=nodes_list)
	grid()
	if (combolength>length(colours)){
		print(paste("WARNING! - There are ", combolength, " methods to display and you're using only ", length(colours), " colours.", sep=""))
	}
	
	legTxt <- c()
	for (nodes in nodes_list){
		dataFilt <- c()
		for (agentType in agentTypes){
			for (method in methods){
				if (agentType == "simple" && method != 'conf-perfect') { ## getting only conf-perfect for simple agents 
					next
				}
				filename <- paste(prefix,"out_net-", netType,"_",agentType, "_nodes-",nodes,"_model-",method,
						"_acc-",accuracy,"_acstdv-",acstdv,".txt",sep="")
				print(filename)
				data <- read.table(filename, header=T)
				dt <- data.table(data)
				tmp <- dt[, list(pos=sum(pos[pos == nodes])/nodes, neg=sum(neg[neg == nodes])/nodes, nRuns=length(exp[pos == nodes | neg == nodes])), by=exp]
				tmp$metric <- tmp$pos / tmp$nRuns
				dataFilt <- c(dataFilt, tmp[,'metric'])
				if (nodes == nodes_list[1]) {
					ltxt <- paste(agentType,method,sep=" ")
					legTxt <- c(legTxt, ltxt)
				}
			}
		}
		size <- xrange/40
		positions = nodes + (seq(0,combolength-1) - floor((combolength-1)/2))*size*1.5
		print(length(dataFilt))
		print(positions)
		boxplot( dataFilt, at=positions, boxwex=size, add=T, col=colours[1:combolength], axes=F)
		
	}
	legend('bottomright', legTxt, fill=colours[1:combolength], cex=1.5, bg='white')
}

diffMatrixSuccessMethod <- function(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.1, 0.6, 0.1), 
		numRuns=100, netType="barabasi-albert", accuracy=0.7, acstdv=0.3, brks=seq(-0.15,0.15,0.015),
		methods=c( "conf-perfect", "best-acc")){
	if (netType == "barabasi-albert"){
		netParam <- edges
		netParamText <- "Num. edges"
	}
	if (netType == "erdos-renyi"){
		netParam <- links
		netParamText <- "Linking probability"
	}
	
	xrange <- range(netParam)
	#plot(xrange,range(nodes_list), type='n', xlab="Number of nodes", ylab="Success rate")
#	grid()
	net <- netType
	sumSign <- +1
	dataMat <- matrix(0, nrow=length(nodes_list), ncol=length(netParam))
	for (method in methods){
		dataFrame <- data.frame()
		n <- 0
		for (nodes in nodes_list){
			n <- n+1
			j <- 0
			for (parVal in netParam){
#					if (method == "best-acc" && net != netTypes[1]) { ## getting only best-acc for the first net-type (because it's the same for others)
#						next
#					}
#					if (net != "erdos-renyi" && link != links[1]) { ## looping on linking-probablity on if net is erdos-renyi
#						next
#					}
#					if (net != "barabasi-albert" && edge != edges[1]) { ## looping on edges on if net is barabasi-albert
#						next
#					}
				j <- j +1
				link <- links[1]
				edge <- edges[1]
				if (net == "erdos-renyi"){
					link <- parVal
				}
				if (net == "barabasi-albert"){
					edge <- parVal
				}
				filename <- paste(prefix,"out_net-", net, "_nodes-",nodes,"_edges-",edge,
						"_link-",link,"_model-",method,"_acc-",accuracy,"_acstdv-",acstdv,".txt",sep="")
#					print(filename)
				data <- read.table(filename, header=T)
				dt <- data.table(data)
				tmp <- dt[, list(pos=sum(pos[pos == nodes])/nodes, neg=sum(neg[neg == nodes])/nodes, nRuns=length(exp[pos == nodes | neg == nodes])), by=exp]
				tmp$metric <- tmp$pos / tmp$nRuns
				
				line <- c(nodes, parVal, mean(tmp$metric) )
				if (length(dataFrame) == 0){
					dataFrame <- line
				} else {
					dataFrame <- rbind(dataFrame, line )
				}
				
				dataMat[n,j] = dataMat[n,j] + (sumSign * mean(tmp$metric))
			}
		}
		sumSign <- -1
	}
	
	positives <- c(-1,-1)
	for (n in seq(1,nrow(dataMat))){
		for (j in seq(1,ncol(dataMat))){
			if (!is.nan(dataMat[n,j]) && dataMat[n,j] > 0) {positives <- rbind(positives, c(n,j))}
		}
	}
	
	#kde2d( dataMat[,1], dataMat[,2], dataMat[,3] )
	cuts <- cut(dataMat, breaks = brks, include.lowest = TRUE)
	#print(cuts)
	#print(levels(cuts))
	#colors <- colorRampPalette(c("white", "red"))(length(levels(cuts)))
	#colors <- colorRampPalette(c("red", "white"))(length(levels(cuts)))
	colors <- colorRampPalette(c('red', 'orange', 'yellow', 'green'))(length(levels(cuts)))
#	colors <- colorRampPalette(c('red', 'orange', 'yellow', 'green', 'cyan', 'blue'))(length(levels(cuts)))
	#print(colors)
	cls <- rep(colors, times = table(cuts))
	
	print(levelplot(dataMat, xlab="nodes", ylab=netParamText, at = brks, 
			colorkey = list(col = colors, at = brks), cuts = (length(colors)), 
			aspect=1, col.regions=colors, #row.values=nodes_list, column.values=edges,
			scales=list(x=list(at=seq(1,length(nodes_list)), labels=nodes_list),
						y=list(at=seq(1,length(netParam)), labels=netParam)),
			panel = function(...){
				panel.fill(col = "black")
				panel.levelplot(...)
				if (!is.null(nrow(positives))) grid.points(positives[,1], positives[,2], pch = '*', gp=gpar(col='blue', lwd=2, cex=2))
			},
	))
	print( paste("model:", method, " net:", net, " acc:", accuracy, " acstdv:", acstdv, sep="") )	
	print( dataMat )
}

diffMatrixSuccessFullnet <- function(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.1, 0.6, 0.1), 
		numRuns=100, netType="barabasi-albert", accuracy=0.7, acstdv=0.3, brks=seq(-0.15,0.15,0.005),
		method="conf-perfect"){
	if (netType == "barabasi-albert"){
		netParam <- edges
		netParamText <- "Num. edges"
	}
	if (netType == "erdos-renyi"){
		netParam <- links
		netParamText <- "Linking probability"
	}
	
	xrange <- range(netParam)
	#plot(xrange,range(nodes_list), type='n', xlab="Number of nodes", ylab="Success rate")
#	grid()
	net <- netType
	dataMat <- matrix(0, nrow=length(nodes_list), ncol=length(netParam))
	dataFrame <- data.frame()
	positives <- c(-1,-1)
	n <- 0
	for (nodes in nodes_list){
		n <- n+1
		
		### Assigning as base value for N-nodes the mean-success from the "full" network 
		filename <- paste(prefix,"out_net-", "full", "_nodes-",nodes,"_edges-",edges[1],
				"_link-",links[1],"_model-",method,"_acc-",accuracy,"_acstdv-",acstdv,".txt",sep="")
#			print(filename)
		data <- read.table(filename, header=T)
		dt <- data.table(data)
		tmp <- dt[, list(pos=sum(pos[pos == nodes])/nodes, neg=sum(neg[neg == nodes])/nodes, nRuns=length(exp[pos == nodes | neg == nodes])), by=exp]
		tmp$metric <- tmp$pos / tmp$nRuns
		dataMat[n,] = -mean(tmp$metric)
		
		j <- 0
		for (parVal in netParam){
			j <- j +1
			link <- links[1]
			edge <- edges[1]
			if (net == "erdos-renyi"){
				link <- parVal
			}
			if (net == "barabasi-albert"){
				edge <- parVal
			}
			filename <- paste(prefix,"out_net-", net, "_nodes-",nodes,"_edges-",edge,
					"_link-",link,"_model-",method,"_acc-",accuracy,"_acstdv-",acstdv,".txt",sep="")
#					print(filename)
			data <- read.table(filename, header=T)
			dt <- data.table(data)
			tmp <- dt[, list(pos=sum(pos[pos == nodes])/nodes, neg=sum(neg[neg == nodes])/nodes, nRuns=length(exp[pos == nodes | neg == nodes])), by=exp]
			tmp$metric <- tmp$pos / tmp$nRuns
			
			dataMat[n,j] = dataMat[n,j] + mean(tmp$metric)
			if (!is.nan(dataMat[n,j]) && dataMat[n,j] > 0) {positives <- rbind(positives, c(n,j))}
		}
	}
	
	cuts <- cut(dataMat, breaks = brks, include.lowest = TRUE)
	#colors <- colorRampPalette(c("white", "red"))(length(levels(cuts)))
	colors <- colorRampPalette(c('red', 'orange', 'yellow', 'green'))(length(levels(cuts)))
	cls <- rep(colors, times = table(cuts))
	
	print(levelplot(dataMat, xlab="nodes", ylab=netParamText, at = brks, 
					colorkey = list(col = colors, at = brks), cuts = (length(colors)), 
					aspect=1, col.regions=colors, #row.values=nodes_list, column.values=edges,
					scales=list(x=list(at=seq(1,length(nodes_list)), labels=nodes_list),
							y=list(at=seq(1,length(netParam)), labels=netParam)),
					panel = function(...){
						panel.fill(col = "black")
						panel.levelplot(...)
						if (!is.null(nrow(positives))) grid.points(positives[,1], positives[,2], pch = '*', gp=gpar(col='blue', lwd=2, cex=2))
#				if (length(cl)>0){ for (i in 1:length(cl)){ panel.lines( cbind(cl[[i]]$y, cl[[i]]$x), col='blue') } }
#				#panel.text(data$rhoA, data$rhoB, substring(as.character(round(mat,2)), 2), cex=cexsize)
#				#panel.text(data$rhoA, data$rhoB, round(mat,2), cex=cexsize)
					},
			))
	print( paste("model:", method, " net:", net, " acc:", accuracy, " acstdv:", acstdv, sep="") )
	print( dataMat )
#	return( dataMat )
}

#############################################################################
################################# TIME ######################################
#############################################################################

plotTimeOnNodes <- function(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.2, 0.6, 0.1), numRuns=100,
		netTypes=c("full", "erdos-renyi", "barabasi-albert"),accuracy=0.7, acstdv=0.3,
		methods=c("best-acc", "conf-perfect"), colours=rainbow(12)){
	combolength <- length(methods)*length(netTypes)
	methodscount <- length(methods)
	if ('best-acc' %in% methods){
		combolength <- combolength - (length(netTypes)-1)
		methodscount <- methodscount-1  
	}
	if ('erdos-renyi' %in% netTypes){ combolength <- combolength + (length(links)-1)*methodscount }
	if ('barabasi-albert' %in% netTypes){ combolength <- combolength + (length(edges)-1)*methodscount }
	offset <- 0.008*combolength
	xmin <- min(nodes_list)
	xmax <- max(nodes_list)
	xrange <- (xmax - xmin)
	plot(c(xmin-(xrange*offset),xmax+(xrange*offset)),c(1,7), type='n', xlab="Number of nodes", ylab="Number of iterations", xaxt = "n")
	axis(1, at=nodes_list)
	grid()
	if (combolength>length(colours)){
		print(paste("WARNING! - There are ", combolength, " methods to display and you're using only ", length(colours), " colours.", sep=""))
	}
	
	legTxt <- c()
	for (nodes in nodes_list){
		dataFilt <- c()
		for (method in methods){
			for (net in netTypes){
				for (edge in edges){
					for (link in links){
						if (method == "best-acc" && (net != netTypes[1] || link != links[1] || edge != edges[1]) ) { ## getting only best-acc for the first net-type (because it's the same for others)
							next
						}
						if (net != "erdos-renyi" && link != links[1]) { ## looping on linking-probablity on if net is erdos-renyi
							next
						}
						if (net != "erdos-renyi") {link <- 0.1}
						if (net != "barabasi-albert" && edge != edges[1]) { ## looping on edges on if net is barabasi-albert
							next
						}
						filename <- paste(prefix,"out_net-", net, "_nodes-",nodes,"_edges-",edge,
								"_link-",link,"_model-",method,"_acc-",accuracy,"_acstdv-",acstdv,".txt",sep="")
						print(filename)
						data <- read.table(filename, header=T)
						dt <- data.table(data)
						dataFilt <- c(dataFilt, dt[, list(metric=mean(iter[pos == nodes | neg == nodes]) ), by=exp][,'metric'])
						if (nodes == nodes_list[1]) {
							ltxt <- paste(net,method,sep=" ")
							if (net == "erdos-renyi"){
								ltxt <- paste(ltxt, " p:",link,sep="")
							}
							if (net == "barabasi-albert"){
								ltxt <- paste(ltxt, " m:",edge,sep="")
							}
							legTxt <- c(legTxt, ltxt)
						}
					}
				}
			}
		}
		size <- xrange/100
		positions = nodes + (seq(0,combolength-1) - floor((combolength-1)/2))*size*1.5
		print(length(dataFilt))
		print(positions)
		boxplot( dataFilt, at=positions, boxwex=size, add=T, col=colours[1:combolength], axes=F)
		
		#aggr <- dt[, list(pos=sum(pos)/(nodes*numRuns), neg=sum(neg)/(nodes*numRuns), best=length(exp[bestDDM == 1])/numRuns), by=exp]
		#boxshift <- 10
		#boxplot( aggr[,'pos'],  at=nodes-boxshift, boxwex=25, add=T)
		#boxplot( aggr[,'best'], at=nodes+boxshift, boxwex=25, col='red', add=T)
	}
#	joinlegend <- function(v1, v2){
#		leg <- c()
#		for (i in v1){
#			for (j in v2){
#				leg <- c(leg, paste(i,j,sep=" "))
#			}
#		}
#		return(leg)
#	}
#	legTxt <- joinlegend(methods,netTypes)
	legend('topleft', legTxt, fill=colours[1:combolength], cex=0.8, bg='white')
}


