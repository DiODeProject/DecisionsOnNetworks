library(data.table)
#library(MASS)
library(lattice)
require(grid)

#############################################################################
################################ AGGREGATES #################################
#############################################################################

plotAllOfAll <- function(prefix, nodes=c(100, 500, 1000), links=seq(0.2,0.8,0.2), drift=0.1, methods=c("conf","M-rand","M-bias","M-inhib","CDCI"), const_methods=c("bestDDM", "confDDM", "FNconf", "FNmaj")){
	colours=rainbow(length(c(methods, const_methods)))
	plotAllEffectivity(prefix=prefix, nodes=nodes, links=links, drift=drift, pdfout=paste(prefix,"im-effectivity.pdf",sep=""), methods=methods, colours=colours)
	plotAllSuccess(prefix=prefix, nodes=nodes, links=links, drift=drift, pdfout=paste(prefix,"im-success.pdf",sep=""), methods=methods, const_methods=const_methods, colours=colours)
	plotAllTime(prefix=prefix, nodes=nodes, links=links, drift=drift, pdfout=paste(prefix,"im-time.pdf",sep=""), methods=methods, colours=colours)
}

plotAllEffectivity <- function(prefix, nodes=c(100, 500, 1000), links=seq(0.2,0.8,0.2), drift=0.1, pdfout, methods=c("conf","M-rand","M-bias","M-inhib"),colours=rainbow(10)){
	pdf(pdfout)
	lines.per.page <- 2
	graph.per.line <- 2
	par(mfrow=c(lines.per.page,graph.per.line), xpd=FALSE, cex.axis=0.5)
	for (node in nodes){
		plotEffectivity(prefix=prefix, nodes=node, links=links, drift=drift, methods=methods, colours=colours)
		title(paste("Effectivity for ", node, " nodes", sep=""))
	}
	for (link in links){
		plotEffectivityOnNodes(prefix=prefix, nodes_list=nodes, link=link, drift=drift, methods=methods, colours=colours)
		title(paste("Effectivity for ", link, " edges", sep=""))
	}
	dev.off()
}

allEffectivityPlots <- function(prefix){
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/eff-nets-maj-rand.pdf")
	plotEffectivityOnNodes(prefix, nodes_list=c(11, 19, 31), edges=seq(2, 10, 2), links=seq(0.1, 0.6, 0.1), numRuns=100,
			netTypes=c("full", "erdos-renyi", "barabasi-albert"), accuracy=0.7, acstdv=0.3,
			methods=c("majority-rand"), colours=rainbow(12))
	dev.off()
	
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/eff-nets-maj-bias.pdf")
	plotEffectivityOnNodes(prefix, nodes_list=c(11, 19, 31), edges=seq(2, 10, 2), links=seq(0.1, 0.6, 0.1), numRuns=100,
			netTypes=c("full", "erdos-renyi", "barabasi-albert"), accuracy=0.7, acstdv=0.3,
			methods=c("majority-bias"), colours=rainbow(12))
	dev.off()
	
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/eff-nets-maj-inhibit.pdf")
	plotEffectivityOnNodes(prefix, nodes_list=c(11, 19, 31), edges=seq(2, 10, 2), links=seq(0.1, 0.6, 0.1), numRuns=100,
			netTypes=c("full", "erdos-renyi", "barabasi-albert"), accuracy=0.7, acstdv=0.3,
			methods=c("majority-inhibit"), colours=rainbow(12))
	dev.off()
	
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/eff-nets-conf-perfect.pdf")
	plotEffectivityOnNodes(prefix, nodes_list=c(11, 19, 31), edges=seq(2, 10, 2), links=seq(0.1, 0.6, 0.1), numRuns=100,
			netTypes=c("full", "erdos-renyi", "barabasi-albert"), accuracy=0.7, acstdv=0.3,
			methods=c("conf-perfect"), colours=rainbow(12))
	dev.off()
}

allSuccessPlots <- function(prefix){
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/succ-onfull-ER-maj-rand.pdf")
	diffMatrixSuccessFullnet(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.1, 0.6, 0.1), 
			numRuns=100, netType="erdos-renyi", accuracy=0.7, acstdv=0.3, brks=seq(-0.10,0.10,0.005),
			method="majority-rand" )
	dev.off()
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/succ-onfull-BA-maj-rand.pdf")
	diffMatrixSuccessFullnet(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.1, 0.6, 0.1), 
			numRuns=100, netType="barabasi-albert", accuracy=0.7, acstdv=0.3, brks=seq(-0.10,0.10,0.005),
			method="majority-rand" )
	dev.off()
	
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/succ-onfull-ER-conf.pdf")
	diffMatrixSuccessFullnet(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.1, 0.6, 0.1), 
			numRuns=100, netType="erdos-renyi", accuracy=0.7, acstdv=0.3, brks=seq(-0.10,0.10,0.005),
			method="conf-perfect" )
	dev.off()
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/succ-onfull-BA-conf.pdf")
	diffMatrixSuccessFullnet(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.1, 0.6, 0.1), 
			numRuns=100, netType="barabasi-albert", accuracy=0.7, acstdv=0.3, brks=seq(-0.10,0.10,0.005),
			method="conf-perfect" )
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

plotAllSuccess <- function(prefix, nodes=c(100, 500, 1000), links=seq(0.2,0.8,0.2), drift=0.1, pdfout, 
		methods=c("conf","M-rand","M-bias","M-inhib"),
		#methods=c("conf","M-rand","bestDDM", "confDDM", "FNconf", "FNmaj"),
		const_methods=c("bestDDM", "confDDM", "FNconf", "FNmaj"), colours=rainbow(10)){
	pdf(pdfout)
	lines.per.page <- 2
	graph.per.line <- 2
	par(mfrow=c(lines.per.page,graph.per.line), xpd=FALSE, cex.axis=0.5)
	for (node in nodes){
		plotSuccess(prefix=prefix, nodes=node, links=links, drift=drift, methods=methods, const_methods=const_methods, colours=colours)
		title(paste("Success rate for ", node, " nodes", sep=""))
	}
	for (link in links){
		plotSuccessOnNodes(prefix=prefix, nodes_list=nodes, link=link, drift=drift, methods=c(methods,const_methods), colours=colours)
		title(paste("Success rate for ", link, " edges", sep=""))
	}
	dev.off()
}

plotAllTime <- function(prefix, nodes=c(100, 500, 1000), links=seq(0.2,0.8,0.2), drift=0.1, pdfout, methods=c("conf","M-rand","M-bias","M-inhib"), colours=rainbow(10)){
	pdf(pdfout)
	lines.per.page <- 2
	graph.per.line <- 2
	par(mfrow=c(lines.per.page,graph.per.line), xpd=FALSE, cex.axis=0.5)
	for (node in nodes){
		plotTime(prefix=prefix, nodes=node, links=links, drift=drift, methods=methods, colours=colours)
		title(paste("Iterations for ", node, " nodes", sep=""))
	}
	for (link in links){
		plotTimeOnNodes(prefix=prefix, nodes_list=nodes, link=link, drift=drift, methods=methods, colours=colours)
		title(paste("Iterations for ", link, " edges", sep=""))
	}
	dev.off()
}

#############################################################################
############################## EFFECTIVITY ##################################
#############################################################################

plotEffectivity <- function(prefix, nodes=20, links=seq(0.2,0.8,0.2), drift=0.1, numRuns=100, methods=c("conf","M-rand","M-bias","M-inhib"), colours=rainbow(12)){
	offset <- 0.02*length(methods)
	xmin <- min(links)
	xmax <- max(links)
	xrange <- (xmax - xmin)
	plot(c(xmin-(xrange*offset),xmax+(xrange*offset)),c(0,1), type='n', xlab="Edges per new node", ylab="Effectivity")
	grid()
	for (link in links){
		filename <- paste(prefix,"out_nodes-",nodes,"_link-",link,"_drift-",0.1,".txt",sep="")
		data <- read.table(filename, header=T)
		
		#print( nrow(data[ data[,'pos']==nodes | data[,'neg']==nodes ,]) )
		dt <- data.table(data)
		
		dataFilt <- c()
		for (method in methods){
			if (method == "conf"){
				dataFilt <- c(dataFilt, dt[, list(eff=length(exp[pos == nodes | neg == nodes])/numRuns ), by=exp][,'eff'])
			} else if (method == "M-rand"){
				dataFilt <- c(dataFilt, dt[, list(eff=length(exp[posMR == nodes | negMR == nodes])/numRuns ), by=exp][,'eff'])
			} else if (method == "M-bias"){
				dataFilt <- c(dataFilt, dt[, list(eff=length(exp[posMB == nodes | negMB == nodes])/numRuns ), by=exp][,'eff'])
			} else if (method == "M-inhib"){
				dataFilt <- c(dataFilt, dt[, list(eff=length(exp[posMI == nodes | negMI == nodes])/numRuns ), by=exp][,'eff'])
			} else if (method == "CDCI"){
				dataFilt <- c(dataFilt, dt[, list(eff=length(exp[posCDCI == nodes | negCDCI == nodes])/numRuns ), by=exp][,'eff'])
			}
			#dt[, list(pos=sum(pos)/(nodes*100), neg=sum(neg)/(nodes*100)), by=exp]
		}
		size <- xrange/40
		positions = link + (seq(0,length(methods)-1) - floor((length(methods)-1)/2))*size*1.5
		#print(positions)
		boxplot( dataFilt, at=positions, boxwex=size, add=T, col=colours[1:length(methods)], axes=F)
	}
	legend('bottomleft', methods, fill=colours[1:length(methods)], cex=0.5, bg='white')
	#return (data)
}

plotEffectivityOnNodes <- function(prefix, nodes_list=c(11, 19, 31), edges=seq(2, 10, 2), links=seq(0.1, 0.6, 0.1), numRuns=100,
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
	offset <- 0.01*combolength
	xmin <- min(nodes_list)
	xmax <- max(nodes_list)
	xrange <- (xmax - xmin)
	plot(c(xmin-(xrange*offset),xmax+(xrange*offset)),c(0,1), type='n', xlab="Number of nodes", ylab="Effectivity")
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
						if (net != "barabasi-albert" && edge != edges[1]) { ## looping on edges on if net is barabasi-albert
							next
						}
						filename <- paste(prefix,"out_net-", net, "_nodes-",nodes,"_edges-",edge,
								"_link-",link,"_model-",method,"_acc-",accuracy,"_acstdv-",acstdv,".txt",sep="")
						#print(filename)
						data <- read.table(filename, header=T)
						dt <- data.table(data)
						dataFilt <- c(dataFilt, dt[, list(eff=length(exp[pos == nodes | neg == nodes])/numRuns ), by=exp][,'eff'])
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
		size <- xrange/80
		positions = nodes + (seq(0,combolength-1) - floor((combolength-1)/2))*size*1.5
		#print(length(dataFilt))
		#print(positions)
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
	legend('bottomright', legTxt, fill=colours[1:combolength], cex=0.5, bg='white')
}

#############################################################################
################################# SUCCESS ###################################
#############################################################################

plotSuccess <- function(prefix, nodes=20, links=seq(0.2,0.8,0.2), drift=0.1, numRuns=100, 
		methods=c("best-acc", "conf-perfect"), const_methods=c(),
		colours=rainbow(10)){
	offset <- 0.02*length(methods)
	xmin <- min(links)
	xmax <- max(links)
	tmp_size <- (xmax-xmin)/40
	xmin <- xmin - floor((length(methods)-1)/2)*tmp_size*1.5 - tmp_size*8
	xrange <- (xmax - xmin)
	plot(c(xmin-(xrange*offset),xmax+(xrange*offset)),c(0,1), type='n', xlab="Edges per new node", ylab="Success rate")
	if (length(const_methods)>0){
		rect(par("usr")[1], par("usr")[3], par("usr")[1] + xrange*0.9/(length(links)+1), par("usr")[4], col=rgb(0.1, 1, 0.1, alpha=0.5), border='NA')
	}
	grid()
	all_data <- data.frame()
	
	if (length(c(methods,const_methods))>length(colours)){
		print(paste("WARNING! - There are more than ", length(colours), " methods and you're using only ", length(colours), " colours.", sep=""))
	}
	
	for (link in links){
		filename <- paste(prefix,"out_nodes-",nodes,"_link-",link,"_drift-",0.1,".txt",sep="")
		data <- read.table(filename, header=T)
		## creating the full-list
		if (is.element("bestDDM", const_methods) | is.element("confDDM", const_methods) | is.element("FNconf", const_methods) | is.element("FNmaj", const_methods)){
			all_data <- rbind(all_data, data)
		}
		
		dt <- data.table(data)
		dataFilt <- c()
		for (method in methods){
			
			if (method == "conf"){
				tmp <- dt[, list(pos=sum(pos[pos == nodes])/nodes, neg=sum(neg[neg == nodes])/nodes, nRuns=length(exp[pos == nodes | neg == nodes])), by=exp]
				tmp$metric <- tmp$pos / tmp$nRuns
				#print(tmp)
				dataFilt <- c(dataFilt, tmp[,'metric'])
			} else if (method == "M-rand"){
				tmp <- dt[, list(pos=sum(posMR[posMR == nodes])/nodes, neg=sum(negMR[negMR == nodes])/nodes, nRuns=length(exp[posMR == nodes | negMR == nodes])), by=exp]
				tmp$metric <- tmp$pos / tmp$nRuns
				dataFilt <- c(dataFilt, tmp[,'metric'])
			} else if (method == "M-bias"){
				tmp <- dt[, list(pos=sum(posMB[posMB == nodes])/nodes, neg=sum(negMB[negMB == nodes])/nodes, nRuns=length(exp[posMB == nodes | negMB == nodes])), by=exp]
				tmp$metric <- tmp$pos / tmp$nRuns
				dataFilt <- c(dataFilt, tmp[,'metric'])
			} else if (method == "M-inhib"){
				tmp <- dt[, list(pos=sum(posMI[posMI == nodes])/nodes, neg=sum(negMI[negMI == nodes])/nodes, nRuns=length(exp[posMI == nodes | negMI == nodes])), by=exp]
				tmp$metric <- tmp$pos / tmp$nRuns
				dataFilt <- c(dataFilt, tmp[,'metric'])
			} else if (method == "CDCI"){
				tmp <- dt[, list(pos=sum(posCDCI[posCDCI == nodes])/nodes, neg=sum(negCDCI[negCDCI == nodes])/nodes, nRuns=length(exp[posCDCI == nodes | negCDCI == nodes])), by=exp]
				tmp$metric <- tmp$pos / tmp$nRuns
				dataFilt <- c(dataFilt, tmp[,'metric'])
			} else if (method == "bestDDM"){
				# this converges always (therefore, we divide by all runs, i.e. numRuns) 
				dataFilt <- c(dataFilt, dt[, list(metric=length(exp[bestDDM == 1])/numRuns), by=exp][,'metric'])
			} else if (method == "confDDM"){
				# this converges always (therefore, we divide by all runs, i.e. numRuns) 
				dataFilt <- c(dataFilt, dt[, list(metric=length(exp[confDDM == 1])/numRuns), by=exp][,'metric'])
			} else if (method == "FNconf"){
				dataFilt <- c(dataFilt, dt[, list(metric=length(exp[FNconf == 1])/numRuns), by=exp][,'metric'])
			} else if (method == "FNmaj"){
				dataFilt <- c(dataFilt, dt[, list(metric=(length(exp[FNmaj == 1])+0.5*length(exp[FNmaj == 0]))/numRuns), by=exp][,'metric'])
			}
		}
		size <- xrange/40
		positions = link + (seq(0,length(methods)-1) - floor((length(methods)-1)/2))*size*1.5
		boxplot( dataFilt, at=positions, boxwex=size, add=T, col=colours[1:length(methods)], axes=F)
		
#		aggr <- dt[, list(pos=sum(pos)/(nodes*numRuns), neg=sum(neg)/(nodes*numRuns), best=length(exp[bestDDM == 1])/numRuns), by=exp]
#		#print(aggr)
#		boxshift <- 0.03
#		boxplot( aggr[,'pos'],  at=link-boxshift, boxwex=0.08, add=T)
#		boxplot( aggr[,'best'], at=link+boxshift, boxwex=0.08, col='red', add=T)
	}
	
	if (length(const_methods) > 0){
		all_data_dt <- data.table(all_data)
		dataFilt <- c()
		## Plotting the constant methods
		for (method in const_methods){
			if (method == "bestDDM"){
				dataFilt <- c(dataFilt, all_data_dt[, list(metric=length(exp[bestDDM == 1])/(numRuns*length(links))), by=list(exp, seed)][,'metric'])
			} else if (method == "confDDM"){
				dataFilt <- c(dataFilt, all_data_dt[, list(metric=length(exp[confDDM == 1])/(numRuns*length(links))), by=list(exp, seed)][,'metric'])
			} else if (method == "FNconf"){
				dataFilt <- c(dataFilt, all_data_dt[, list(metric=length(exp[FNconf == 1])/(numRuns*length(links))), by=list(exp, seed)][,'metric'])
			} else if (method == "FNmaj"){
				dataFilt <- c(dataFilt, all_data_dt[, list(metric=(length(exp[FNmaj == 1]) + 0.5*length(exp[FNmaj == 0]))/(numRuns*length(links))), by=list(exp, seed)][,'metric'])
			}
		}
		positions = xmin + ( (seq(0,length(const_methods)-1) - floor((length(const_methods)-1)/2))*size*1.5 ) - size*3
		boxplot( dataFilt, at=positions, boxwex=size, add=T, col=colours[(length(methods)+1):(length(methods)+length(const_methods))], axes=F)
	}
	
#	if (is.element("bestDDM", const_methods) | is.element("confDDM", const_methods) | is.element("FNconf", const_methods) | is.element("FNmaj", const_methods)){
#		for (link in links){
#			filename <- paste(prefix,"out_nodes-",nodes,"_link-",link,"_drift-",0.1,".txt",sep="")
#			data <- read.table(filename, header=T)
#			all_data <- rbind(all_data, data)
#		}
#		all_data_dt <- data.table(all_data)
#		size <- 0.02
#		values <- all_data_dt[, list(metric=length(exp[bestDDM == 1])/(numRuns*length(links))), by=list(exp, seed)][,'metric']
#		print(paste("values:" , values))
#		boxplot( values, at=0, boxwex=size, add=T, col='red', axes=F)
#	}
	
	legend('bottomleft', c(methods, const_methods), fill=colours[1:(length(methods)+length(const_methods))], cex=0.5, bg='white')
}

plotSuccessOnNodes <- function(prefix, nodes_list=seq(11, 31, 4), edges=c(2, 6, 8), links=c(0.1, 0.3, 0.5), drift=0.1, numRuns=100,
		netTypes=c("full", "erdos-renyi", "barabasi-albert"),accuracy=0.7, acstdv=0.3,
		methods=c("best-acc", "conf-perfect"), colours=rainbow(10)){
	combolength <- length(methods)*length(netTypes)
	methodscount <- length(methods)
	if ('best-acc' %in% methods){
		combolength <- combolength - (length(netTypes)-1)
		methodscount <- methodscount-1  
	}
	if ('erdos-renyi' %in% netTypes){ combolength <- combolength + (length(links)-1)*methodscount }
	if ('barabasi-albert' %in% netTypes){ combolength <- combolength + (length(edges)-1)*methodscount }
	offset <- 0.02*combolength
	xmin <- min(nodes_list)
	xmax <- max(nodes_list)
	xrange <- (xmax - xmin)
	plot(c(xmin-(xrange*offset),xmax+(xrange*offset)),c(0,1), type='n', xlab="Number of nodes", ylab="Success rate", xaxt = "n")
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
						if (net != "barabasi-albert" && edge != edges[1]) { ## looping on edges on if net is barabasi-albert
							next
						}
						filename <- paste(prefix,"out_net-", net, "_nodes-",nodes,"_edges-",edge,
								"_link-",link,"_model-",method,"_acc-",accuracy,"_acstdv-",acstdv,".txt",sep="")
						print(filename)
						data <- read.table(filename, header=T)
						dt <- data.table(data)
						tmp <- dt[, list(pos=sum(pos[pos == nodes])/nodes, neg=sum(neg[neg == nodes])/nodes, nRuns=length(exp[pos == nodes | neg == nodes])), by=exp]
						tmp$metric <- tmp$pos / tmp$nRuns
						dataFilt <- c(dataFilt, tmp[,'metric'])
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
		size <- xrange/40
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


plotTime <- function(prefix, nodes=20, links=seq(0.2,0.8,0.2), drift=0.1, methods=c("conf","M-rand","M-bias","M-inhib"), colours=rainbow(10)){
	offset <- 0.02*length(methods)
	xmin <- min(links)
	xmax <- max(links)
	xrange <- (xmax - xmin)
	plot(c(xmin-(xrange*offset),xmax+(xrange*offset)),c(0,5), type='n', xlab="Edges per new node", ylab="Number of iterations")
	grid()
	for (link in links){
		# Read the data
		filename <- paste(prefix,"out_nodes-",nodes,"_link-",link,"_drift-",0.1,".txt",sep="")
		data <- read.table(filename, header=T)
		
		# Filter: keeping only converged runs
		#data[ data[,'pos']==nodes | data[,'neg']==nodes ,]
		
		# Computing mean
		dt <- data.table(data)
		dataFilt <- c()
		for (method in methods){
			if (method == "conf"){
				dataFilt <- c(dataFilt, dt[, list(metric=mean(iter[pos == nodes | neg == nodes]) ), by=exp][,'metric'])
			} else if (method == "M-rand"){
				dataFilt <- c(dataFilt, dt[, list(metric=mean(iterMR[posMR == nodes | negMR == nodes]) ), by=exp][,'metric'])
			} else if (method == "M-bias"){
				dataFilt <- c(dataFilt, dt[, list(metric=mean(iterMB[posMB == nodes | negMB == nodes]) ), by=exp][,'metric'])
			} else if (method == "M-inhib"){
				dataFilt <- c(dataFilt, dt[, list(metric=mean(iterMI[posMI == nodes | negMI == nodes]) ), by=exp][,'metric'])
			} else if (method == "CDCI"){
				dataFilt <- c(dataFilt, dt[, list(metric=mean(iterCDCI[posCDCI == nodes | negCDCI == nodes]) ), by=exp][,'metric']/10)
			}
		}
		# Plotting
		size <- xrange/40
		positions = link + (seq(0,length(methods)-1) - floor((length(methods)-1)/2))*size*1.5
		boxplot( dataFilt, at=positions, boxwex=size, add=T, col=colours[1:length(methods)], axes=F)
		
		#aggr <- dt[, list(time=mean(iter) ), by=exp]		
		#boxplot( aggr[,'time'], at=link, boxwex=0.1, add=T)
	}
	legend('topleft', methods, fill=colours[1:length(methods)], cex=0.5, bg='white')
}

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


