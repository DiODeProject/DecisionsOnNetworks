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

allEffectivityPlots_old <- function(prefix){
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

allEffectivityPlots_explore <- function(prefix, nodes=c(11, 23, 31), acc=0.6){
	for (n in nodes){
		pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/eff-upds-erdos-n",n,".pdf",sep=""))
		plotEffectivityOnNetParam(prefix, nodes_list=c(n), param=seq(0.3, 0.8, 0.2), xlabel="Parameter p_e (linking probability)", numRuns=1,
				netType="erdos-renyi", accuracy=acc, acstdv=0.12, methods=c("majority-rand", "majority-bias", "majority-inhibit", "conf-perfect"),
				edge=2, link=0.3, updates=c("no-up", "theta-up", "theta-norm"), colours=rainbow(12), bxplt=F, truncated='true' )
		title(paste("Effectivity for ", n, " nodes on Erdos-Renyi network", sep=""))
		dev.off()
		
		pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/eff-upds-barabasi-n",n,".pdf",sep=""))
		plotEffectivityOnNetParam(prefix, nodes_list=c(n), param=seq(2, 10, 4), xlabel="Parameter m (number of edges)", numRuns=1,
				netType="barabasi-albert", accuracy=acc, acstdv=0.12, methods=c("majority-rand", "majority-bias", "majority-inhibit", "conf-perfect"),
				edge=2, link=0.3, updates=c("no-up", "theta-up", "theta-norm"), colours=rainbow(12), bxplt=F, truncated='true' )
		title(paste("Effectivity for ", n, " nodes on Barabasi-Albert network", sep=""))
		dev.off()
	}
}

allEffectivityPlots <- function(prefix, nodes=c(11, 23, 31), acc=0.6){
	colours <- rainbow(6)
	for (n in nodes){
		pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/eff-upds-erdos-n",n,".pdf",sep=""))
		plotEffectivityOnNetParam(prefix, nodes_list=c(n), param=seq(0.3, 0.8, 0.1), xlabel="Parameter p_e (linking probability)", numRuns=1,
				netType="erdos-renyi", accuracy=acc, acstdv=0.12, 
				#				methods=c("majority-rand", "conf-perfect"), updates=c("no-up", "optim-up", "belief-up"), combineMethods=TRUE,
				methods=c("majority-rand", "majority-rand", "conf-perfect", "conf-perfect", "conf-perfect", "conf-perfect"), updates=c("no-up", "optim-up", "no-up", "optim-up", "belief-up", "theta-up"), combineMethods=FALSE,
				edge=2, link=0.3, colours=colours, bxplt=F, truncated='true' )
		title(paste("Effectivity for ", n, " nodes on Erdos-Renyi network", sep=""))
		dev.off()
		
		pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/eff-upds-barabasi-n",n,".pdf",sep=""))
		plotEffectivityOnNetParam(prefix, nodes_list=c(n), param=seq(2, 10, 2), xlabel="Parameter m (number of edges)", numRuns=1,
				netType="barabasi-albert", accuracy=acc, acstdv=0.12,
				#				methods=c("majority-rand", "conf-perfect"), updates=c("no-up", "optim-up", "belief-up"), combineMethods=TRUE,
				methods=c("majority-rand", "majority-rand", "conf-perfect", "conf-perfect", "conf-perfect", "conf-perfect"), updates=c("no-up", "optim-up", "no-up", "optim-up", "belief-up", "theta-up"), combineMethods=FALSE,
				edge=2, link=0.3, colours=colours, bxplt=F, truncated='true' )
		title(paste("Effectivity for ", n, " nodes on Barabasi-Albert network", sep=""))
		dev.off()
	}
}

allSuccessPlots <- function(prefix, nodes=c(11, 23, 31), acc=0.6, acstdv=0.12){
	colours <- rainbow(6)
	for (n in nodes){
		pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/succ-methods-ER-n",n,".pdf",sep=""))
		plotSuccessOnNetParam(prefix, nodes_list=c(n), param=seq(0.3, 0.8, 0.1), xlabel="Parameter p_e (linking probability)",
				netType="erdos-renyi", accuracy=acc, acstdv=acstdv, 
#				methods=c("majority-rand", "conf-perfect"), updates=c("no-up", "optim-up", "belief-up"), combineMethods=TRUE,
				methods=c("majority-rand", "majority-rand", "conf-perfect", "conf-perfect", "conf-perfect", "conf-perfect"), updates=c("no-up", "optim-up", "no-up", "optim-up", "belief-up", "theta-up"), combineMethods=FALSE,
				edge=2, link=0.3, colours=colours, yrange=c(0.5,1), truncated="true", bxplt=FALSE )
		title(paste("Success rate for ", n, " nodes on Erdos-Renyi network", sep=""))
		dev.off()
		
		pdf(paste("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/succ-methods-BA-n",n,".pdf",sep=""))
		plotSuccessOnNetParam(prefix, nodes_list=c(n), param=seq(2, 10, 2), xlabel="Parameter m (number of edges)",
				netType="barabasi-albert", accuracy=acc, acstdv=acstdv, 
#				methods=c("majority-rand", "conf-perfect"), updates=c("no-up", "optim-up", "belief-up"), combineMethods=TRUE,
				methods=c("majority-rand", "majority-rand", "conf-perfect", "conf-perfect", "conf-perfect", "conf-perfect"), updates=c("no-up", "optim-up", "no-up", "optim-up", "belief-up", "theta-up"), combineMethods=FALSE,
				edge=2, link=0.3, colours=colours, yrange=c(0.5,1), truncated="true", bxplt=FALSE )
		title(paste("Success rate for ", n, " nodes on Barabasi-Albert network", sep=""))
		dev.off()
	}
}

fullNetSuccRate <- function(prefix){
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/succ-full.pdf")
	plotSuccessOnNodes(prefix, nodes_list=seq(11, 31, 4), edges=c(2, 6, 8), links=c(0.1, 0.3, 0.5), 
			netTypes=c("full"),accuracy=0.6, acstdv=0.12, updates=c("no-up"), truncated="true", bxplt = FALSE,
			methods=c("best-acc", "majority-rand", "conf-perfect"), colours=rainbow(3), yrange=c(0.7, 1))
	dev.off()
}

allSuccessPlotsMtxVsFull <- function(prefix, acc=0.6, stdev=0.12){
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/succ-onfull-ER-maj-rand.pdf")
	diffMatrixSuccessFullnet(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.3, 0.8, 0.1), 
			netType="erdos-renyi", accuracy=acc, acstdv=stdev, brks=seq(-0.06,0.06,0.0025), 
			update="optim-up", method="majority-rand", truncated='true', singleRun=TRUE )
	dev.off()
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/succ-onfull-BA-maj-rand.pdf")
	diffMatrixSuccessFullnet(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.3, 0.8, 0.1), 
			netType="barabasi-albert", accuracy=acc, acstdv=stdev, brks=seq(-0.06,0.06,0.0025),
			update="optim-up", method="majority-rand", truncated='true', singleRun=TRUE )
	dev.off()
	
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/succ-onfull-ER-conf.pdf")
	diffMatrixSuccessFullnet(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.3, 0.8, 0.1), 
			netType="erdos-renyi", accuracy=acc, acstdv=stdev, brks=seq(-0.06,0.06,0.0025),
			update="optim-up", method="conf-perfect", truncated='true', singleRun=TRUE )
	dev.off()
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/succ-onfull-BA-conf.pdf")
	diffMatrixSuccessFullnet(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.3, 0.8, 0.1), 
			netType="barabasi-albert", accuracy=acc, acstdv=stdev, brks=seq(-0.06,0.06,0.0025),
			update="optim-up", method="conf-perfect", truncated='true', singleRun=TRUE )
	dev.off()
}

allSuccessMethodPlots <- function(prefix, acc=0.6, stdev=0.12){
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/succ-onbest-ER-maj-rand.pdf")
	diffMatrixSuccessMethod(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.3, 0.8, 0.1), 
		netType="erdos-renyi", accuracy=acc, acstdv=stdev, brks=seq(-0.15,0.15,0.005),
		update="optim-up", methods=c("majority-rand", "best-acc"), truncated='true', singleRun=TRUE)
	dev.off()
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/succ-onbest-ER-conf.pdf")
	diffMatrixSuccessMethod(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.3, 0.8, 0.1), 
		netType="erdos-renyi", accuracy=acc, acstdv=stdev, brks=seq(-0.15,0.15,0.005),
		update="optim-up", methods=c("conf-perfect", "best-acc"), truncated='true', singleRun=TRUE)
	dev.off()

	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/succ-onbest-BA-maj-rand.pdf")
	diffMatrixSuccessMethod(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.3, 0.8, 0.1), 
		netType="barabasi-albert", accuracy=acc, acstdv=stdev, brks=seq(-0.15,0.15,0.005),
		update="optim-up", methods=c("majority-rand", "best-acc"), truncated='true', singleRun=TRUE)
	dev.off()
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/succ-onbest-BA-conf.pdf")
	diffMatrixSuccessMethod(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.3, 0.8, 0.1), 
		netType="barabasi-albert", accuracy=acc, acstdv=stdev, brks=seq(-0.15,0.15,0.005),
		update="optim-up", methods=c("conf-perfect", "best-acc"), truncated='true', singleRun=TRUE)
	dev.off()
}

allTimePlots <- function(prefix, acc=0.6, stdev=0.12){
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/time-nets-erdos.pdf", width=17, height=7)
#	plotTimeOnNodes(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.3, 0.8, 0.1), 
#			netTypes=c("erdos-renyi", "barabasi-albert"), accuracy=acc, acstdv=stdev, truncated='true',
#			updates=c("optim-up"), methods=c("majority-rand"), colours=rainbow(12), yrange=c(1,20), singleRun=TRUE)
	plotTimeTwoMethodsOnNodes(prefix=prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), 
			links=seq(0.3, 0.8, 0.1), netTypes=c("erdos-renyi"), accuracy=acc, acstdv=stdev, truncated='true', 
			updates=c("theta-up", "optim-up"), methods=c("majority-rand","conf-perfect"), colours=rainbow(6), 
			yrange=c(0,8), singleRun=TRUE)
	dev.off()
	
#	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/time-nets-conf-theta.pdf", width=17, height=7)
#	plotTimeOnNodes(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.3, 0.8, 0.1), 
#			netTypes=c("erdos-renyi", "barabasi-albert"), accuracy=acc, acstdv=stdev, truncated='true',
#			updates=c("theta-up"), methods=c("conf-perfect"), colours=rainbow(12), yrange=c(1,20), singleRun=TRUE)
#	dev.off()
	
	pdf("/Users/joefresna/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/time-nets-barabasi.pdf", width=17, height=7)
	plotTimeTwoMethodsOnNodes(prefix=prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), 
			links=seq(0.3, 0.8, 0.1), netTypes=c("barabasi-albert"), accuracy=acc, acstdv=stdev, truncated='true', 
			updates=c("theta-up", "optim-up"), methods=c("majority-rand","conf-perfect"), colours=rainbow(6), 
			yrange=c(0,8), singleRun=TRUE)
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

plotEffectivityOnNetParam <- function(prefix, nodes_list=c(23), param=seq(2, 10, 4), xlabel="Parameter m (number of edges)",numRuns=100,
		netType="barabasi-albert", accuracy=0.6, acstdv=0.12, methods=c("majority-rand", "conf-perfect"), edge=2, link=0.3,
		updates=c("no-up", "theta-up", "theta-norm"), colours=rainbow(10), yrange=c(0,1), bxplt=TRUE, truncated='true', combineMethods=FALSE ) {
	
	pntTypes = c(1,4,5,3,8,7,0,2,6,9,10,11,12,13)
	
	if (combineMethods){
		combolength <- length(methods)
		methodscount <- length(methods)
		if ('best-acc' %in% methods){
			combolength <- combolength - (length(netTypes)-1)
			methodscount <- methodscount-1  
		}
		combolength <- combolength*length(updates)
		updates_loop <- updates
	} else {
		combolength <- length(methods)
		updates_loop <- c(1)
	}
	
	
	offset <- 0.01*combolength
	xmin <- min(param)
	xmax <- max(param)
	xrange <- (xmax - xmin)
	plot(c(xmin-(xrange*offset),xmax+(xrange*offset)), yrange, type='n', xlab=xlabel, ylab="Effectivity", xaxt = "n")
	axis(1, at=param)
	grid()
	if (combolength>length(colours)){
		print(paste("WARNING! - There are ", combolength, " methods to display and you're using only ", length(colours), " colours.", sep=""))
	}
	
	legTxt <- c()
	for (netPar in param){
		dataFilt <- c()
		for (nodes in nodes_list){
			m <- 0
			for (method in methods){
				m <- m+1
				for (update in updates_loop){
					if (method == "best-acc" && (netType != netTypes[1] || link != links[1] || edge != edges[1]) ) { ## getting only best-acc for the first net-type (because it's the same for others)
						next
					}
					if (netType == "erdos-renyi") {link <- netPar}
					if (netType == "barabasi-albert") {edge <- netPar}
					if (!combineMethods) {
						update <- updates[m]
					}
					filename <- paste(prefix,"out_net-", netType, "_nodes-",nodes,"_edges-",edge,
							"_link-",link,"_model-",method,"_up-",update,"_acc-",accuracy,"_acstdv-",acstdv,"_trunc-",truncated,".txt",sep="")
#							print(filename)
					data <- read.table(filename, header=T)
					dt <- data.table(data)
					dataFilt <- c(dataFilt, dt[, list(eff=length(exp[pos == nodes | neg == nodes])/numRuns ), by=exp][,'eff'])
					if (netPar == param[1]) {
						ltxt <- paste(method,update,sep=" ")
						ltxt <- gsub("perfect", "weigh", ltxt)
						ltxt <- gsub("rand", "rule", ltxt)
						ltxt <- gsub("optim-up", "with-up", ltxt)
						legTxt <- c(legTxt, ltxt)
					}
				}
			}
		}
		size <- xrange/70
		positions = netPar + (seq(0,combolength-1) - floor((combolength-1)/2))*size*1.5
		#print(length(dataFilt))
		#print(positions)
		if (bxplt){
			boxplot( dataFilt, at=positions, boxwex=size, add=T, col=colours[1:combolength], axes=F, outcex=0.5)
		} else {
			print( paste(netType, netPar) )
			print( data.frame(lapply(dataFilt, mean, na.rm=TRUE)) )
			points( positions, lapply(dataFilt, mean, na.rm=TRUE), pch=pntTypes[1:combolength], col=colours[1:combolength], cex=1, lwd=2 )
		}
	}
	if (bxplt){
		legend('bottomright', legTxt, fill=colours[1:combolength], cex=1.5, bg='white')
	} else {
		legend('bottomright', legTxt, pch=pntTypes[1:combolength], col=colours[1:combolength], lty=0, cex=1, lwd=2, bg='white' )
	}
}

plotEffectivityOnNodes <- function(prefix, nodes_list=c(11, 19, 31), edges=seq(2, 10, 2), links=seq(0.3, 0.8, 0.1), numRuns=100,
		netTypes=c("full", "erdos-renyi", "barabasi-albert"),accuracy=0.7, acstdv=0.3, methods=c("best-acc", "conf-perfect"), 
		updates=c("no-up","theta-up","theta-norm"), colours=rainbow(12)){
	combolength <- length(methods)*length(netTypes)*length(updates)
	methodscount <- length(methods)*length(updates)
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
						for (update in updates){
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
									"_link-",link,"_model-",method,"_up-",update,"_acc-",accuracy,"_acstdv-",acstdv,".txt",sep="")
							#print(filename)
							data <- read.table(filename, header=T)
							dt <- data.table(data)
							dataFilt <- c(dataFilt, dt[, list(eff=length(exp[pos == nodes | neg == nodes])/numRuns ), by=exp][,'eff'])
							if (nodes == nodes_list[1]) {
								ltxt <- paste(net,method,update,sep=" ")
								if (net == "erdos-renyi"){
									ltxt <- paste(ltxt, " p:",link,sep="")
								}
								if (net == "barabasi-albert"){
									ltxt <- paste(ltxt, " m:",edge,sep="")
								}
#								if (update){
#									ltxt <- paste(ltxt, " w/ update",sep="")
#								} else {
#									ltxt <- paste(ltxt, " w/o update",sep="")
#								}
								ltxt <- gsub("perfect", "weigh", ltxt)
								ltxt <- gsub("rand", "rule", ltxt)
								legTxt <- c(legTxt, ltxt)
							}
						}
					}
				}
			}
		}
		size <- xrange/70
		positions = nodes + (seq(0,combolength-1) - floor((combolength-1)/2))*size*1.5
		#print(length(dataFilt))
		#print(positions)
		if (bxplt){
			boxplot( dataFilt, at=positions, boxwex=size, add=T, col=colours[1:combolength], axes=F, outcex=0.5)
		} else {
			#return(dataFilt)
			points( positions, lapply(dataFilt, mean, na.rm=TRUE), pch=pntTypes[1:combolength], col=colours[1:combolength], cex=1, lwd=2 )
		}
		
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
	if (bxplt){
		legend('bottomright', legTxt, fill=colours[1:combolength], cex=1.5, bg='white')
	} else {
		legend('bottomright', legTxt, pch=pntTypes[1:combolength], col=colours[1:combolength], lty=0, cex=1, lwd=2, bg='white' )
	}
}

#############################################################################
################################# SUCCESS ###################################
#############################################################################

plotSuccessOnNetParam <- function(prefix, nodes_list=c(23), param=seq(2, 10, 4), xlabel="Parameter m (number of edges)",
		netType="barabasi-albert", accuracy=0.7, acstdv=0.3, methods=c("majority-rand", "conf-perfect"), edge=2, link=0.3,
		updates=c("no-up", "theta-up", "theta-norm"), colours=rainbow(10), yrange=c(0,1), truncated="true", bxplt=TRUE,
		combineMethods=FALSE) {
	pntTypes = c(1,4,5,3,8,7,0,2,6,9,10,11,12,13)
	
	if (combineMethods){
		combolength <- length(methods)
		methodscount <- length(methods)
		if ('best-acc' %in% methods){
			combolength <- combolength - (length(netTypes)-1)
			methodscount <- methodscount-1  
		}
		if ('belief-up' %in% updates){
			combolength <- combolength*(length(updates)-1) +1
		} else {		
			combolength <- combolength*length(updates)
		}
		updates_loop <- updates
	} else {
		combolength <- length(methods)
		updates_loop <- c(1)
	}
	
	offset <- 0.01*combolength
	xmin <- min(param)
	xmax <- max(param)
	xrange <- (xmax - xmin)
	plot(c(xmin-(xrange*offset),xmax+(xrange*offset)), yrange, type='n', xlab=xlabel, ylab="Success rate", xaxt = "n")
	axis(1, at=param)
	grid()
	if (combolength>length(colours)){
		print(paste("WARNING! - There are ", combolength, " methods to display and you're using only ", length(colours), " colours.", sep=""))
	}
	
	legTxt <- c()
	for (netPar in param){
		dataFilt <- c()
		for (nodes in nodes_list){
			m <- 0
			for (method in methods){
				m <- m+1
				for (update in updates_loop){
					if (method == "best-acc" && (netType != netTypes[1] || link != links[1] || edge != edges[1]) ) { ## getting only best-acc for the first net-type (because it's the same for others)
						next
					}
					if (combineMethods && update == "belief-up" && method != "conf-perfect" ) { ## getting only best-acc for the first net-type (because it's the same for others)
						next
					}
					if (!combineMethods){
						update <- updates[m] 
					}
					if (netType == "erdos-renyi") {link <- netPar}
					if (netType == "barabasi-albert") {edge <- netPar}
					filename <- paste(prefix,"out_net-", netType, "_nodes-",nodes,"_edges-",edge,
							"_link-",link,"_model-",method,"_up-",update,"_acc-",accuracy,"_acstdv-",acstdv,"_trunc-",truncated,".txt",sep="")
#							print(filename)
					data <- read.table(filename, header=T)
					if (bxplt){
						dt <- data.table(data)
						tmp <- dt[, list(pos=sum(pos[pos == nodes])/nodes, neg=sum(neg[neg == nodes])/nodes, nRuns=length(exp[pos == nodes | neg == nodes])), by=exp]
						tmp$metric <- tmp$pos / tmp$nRuns
						dataFilt <- c(dataFilt, tmp[,'metric'])
					} else {
						#dataFilt <- c(dataFilt, nrow(data[ data$pos == nodes, ]) / nrow( data[ data$pos == nodes | data$neg == nodes, ] ) )
						dataFilt <- c(dataFilt, nrow(data[ data$pos == nodes, ]) / nrow(data) )
					}
					if (netPar == param[1]) {
						ltxt <- paste(method,update,sep=" ")
						ltxt <- gsub("perfect", "weigh", ltxt)
						ltxt <- gsub("rand", "rule", ltxt)
						ltxt <- gsub("optim-up", "with-up", ltxt)
						legTxt <- c(legTxt, ltxt)
					}
				}
			}
		}
		size <- xrange/70
		positions = netPar + (seq(0,combolength-1) - floor((combolength-1)/2))*size*1.5
		if (bxplt){
			#print(length(dataFilt))
			#print(positions)
			boxplot( dataFilt, at=positions, boxwex=size, add=T, col=colours[1:combolength], axes=F, outcex=0.5)
		} else {
			points( positions, dataFilt, pch=pntTypes[1:combolength], col=colours[1:combolength], cex=1, lwd=2 )
			print(dataFilt)
		}
		
	}
	if (bxplt){
		legend('bottomright', legTxt, fill=colours[1:combolength], cex=1, bg='white')
	} else {
		legend('bottomright', legTxt, pch=pntTypes[1:combolength], col=colours[1:combolength], lty=0, cex=1, lwd=2, bg='white' )
	}
}

plotSuccessOnNodes <- function(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 4), links=seq(0.3, 0.8, 0.1), 
		netTypes=c("full", "erdos-renyi", "barabasi-albert"),accuracy=0.7, acstdv=0.3, methods=c("best-acc", "conf-perfect"),
		updates=c("no-up","theta-up","theta-norm"), colours=rainbow(10), yrange=c(0,1), truncated="true", bxplt=TRUE ){
	
	pntTypes = c(1,4,5,3,8,7,0,2,6,9,10,11,12,13)
	combolength <- length(methods)*length(netTypes)
	methodscount <- length(methods)
	if ('best-acc' %in% methods){
		combolength <- combolength - (length(netTypes)-1)
		methodscount <- methodscount-1  
	}
	if ('erdos-renyi' %in% netTypes){ combolength <- combolength + (length(links)-1)*methodscount }
	if ('barabasi-albert' %in% netTypes){ combolength <- combolength + (length(edges)-1)*methodscount }
	combolength <- combolength*length(updates)
	offset <- 0.01*combolength
	xmin <- min(nodes_list)
	xmax <- max(nodes_list)
	xrange <- (xmax - xmin)
	plot(c(xmin-(xrange*offset),xmax+(xrange*offset)), yrange, type='n', xlab="Number of nodes", ylab="Success rate", xaxt = "n")
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
						for (update in updates){
							if (method == "best-acc" && (net != netTypes[1] || link != links[1] || edge != edges[1]) ) { ## getting only best-acc for the first net-type (because it's the same for others)
								next
							}
							if (net != "erdos-renyi" && link != links[1]) { ## looping on linking-probablity on if net is erdos-renyi
								next
							}
							if (net != "barabasi-albert" && edge != edges[1]) { ## looping on edges on if net is barabasi-albert
								next
							}
							if (net != "erdos-renyi") {link <- 0.3}
							filename <- paste(prefix,"out_net-", net, "_nodes-",nodes,"_edges-",edge,
									"_link-",link,"_model-",method,"_up-",update,"_acc-",accuracy,"_acstdv-",acstdv,"_trunc-",truncated,".txt",sep="")
#							print(filename)
							data <- read.table(filename, header=T)
							if (bxplt){
								dt <- data.table(data)
								tmp <- dt[, list(pos=sum(pos[pos == nodes])/nodes, neg=sum(neg[neg == nodes])/nodes, nRuns=length(exp[pos == nodes | neg == nodes])), by=exp]
								tmp$metric <- tmp$pos / tmp$nRuns
								dataFilt <- c(dataFilt, tmp[,'metric'])
							} else {
								dataFilt <- c(dataFilt, nrow(data[ data$pos == nodes, ]) / nrow( data[ data$pos == nodes | data$neg == nodes, ] ) )
							}
							if (nodes == nodes_list[1]) {
								if (length(methods)>1 && length(netTypes)>1) {
									ltxt <- paste(net,method,sep=" ")
								} else {
									if (length(methods)==1) {
										ltxt <- net 
									} else {
										ltxt <- method
									}
								}
								if (length(updates)>1){
									ltxt <- paste(ltxt,update,sep=" ")
								}
								if (net == "erdos-renyi"){
									ltxt <- paste(ltxt, " p:",link,sep="")
								}
								if (net == "barabasi-albert"){
									ltxt <- paste(ltxt, " m:",edge,sep="")
								}
								ltxt <- gsub("perfect", "weigh", ltxt)
								ltxt <- gsub("rand", "rule", ltxt)
								legTxt <- c(legTxt, ltxt)
							}
						}
					}
				}
			}
		}
		size <- xrange/70
		positions = nodes + (seq(0,combolength-1) - floor((combolength-1)/2))*size*1.5
		if (bxplt){
			#print(length(dataFilt))
			#print(positions)
			boxplot( dataFilt, at=positions, boxwex=size, add=T, col=colours[1:combolength], axes=F, outcex=0.5)
		} else {
			points( positions, dataFilt, pch=pntTypes[1:combolength], col=colours[1:combolength], cex=1, lwd=2 )
			print(dataFilt)
		}
		
	}
	if (bxplt){
		legend('bottomright', legTxt, fill=colours[1:combolength], cex=1.5, bg='white')
	} else {
		legend('bottomright', legTxt, pch=pntTypes[1:combolength], col=colours[1:combolength], lty=0, cex=1.5, lwd=2, bg='white' )
	}
}

diffMatrixSuccessMethod <- function(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.3, 0.8, 0.1), 
		netType="barabasi-albert", accuracy=0.7, acstdv=0.3, brks=seq(-0.15,0.15,0.015), update="theta-up",
		methods=c( "conf-perfect", "best-acc"), truncated='true', singleRun=TRUE){
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
						"_link-",link,"_model-",method,"_up-",update,"_acc-",accuracy,"_acstdv-",acstdv,"_trunc-",truncated,".txt",sep="")
#					print(filename)
				data <- read.table(filename, header=T)
				if (singleRun) {
					#dataMat[n,j] = dataMat[n,j] + (sumSign * nrow(data[ data$pos == nodes, ]) / nrow( data[ data$pos == nodes | data$neg == nodes, ] ) )
					dataMat[n,j] = dataMat[n,j] + (sumSign * nrow(data[ data$pos == nodes, ]) / nrow( data ) )
				} else { 
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
		netType="barabasi-albert", accuracy=0.7, acstdv=0.3, brks=seq(-0.15,0.15,0.005), update="theta-up",
		method="conf-perfect", truncated='true', singleRun=TRUE){
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
				"_link-",links[1],"_model-",method,"_up-",update,"_acc-",accuracy,"_acstdv-",acstdv,"_trunc-",truncated,".txt",sep="")
#			print(filename)
		data <- read.table(filename, header=T)
		if (singleRun){
			#dataMat[n,] = -nrow(data[ data$pos == nodes, ]) / nrow( data[ data$pos == nodes | data$neg == nodes, ] )
			dataMat[n,] = -nrow(data[ data$pos == nodes, ]) / nrow( data )
		} else {
			dt <- data.table(data)
			tmp <- dt[, list(pos=sum(pos[pos == nodes])/nodes, neg=sum(neg[neg == nodes])/nodes, nRuns=length(exp[pos == nodes | neg == nodes])), by=exp]
			tmp$metric <- tmp$pos / tmp$nRuns
			dataMat[n,] = -mean(tmp$metric)
		}
		
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
					"_link-",link,"_model-",method,"_up-",update,"_acc-",accuracy,"_acstdv-",acstdv,"_trunc-",truncated,".txt",sep="")
#					print(filename)
			data <- read.table(filename, header=T)
			if (singleRun){
				#dataMat[n,j] = dataMat[n,j] + (nrow(data[ data$pos == nodes, ]) / nrow( data[ data$pos == nodes | data$neg == nodes, ] ))
				dataMat[n,j] = dataMat[n,j] + (nrow(data[ data$pos == nodes, ]) / nrow( data ))
			} else {
				dt <- data.table(data)
				tmp <- dt[, list(pos=sum(pos[pos == nodes])/nodes, neg=sum(neg[neg == nodes])/nodes, nRuns=length(exp[pos == nodes | neg == nodes])), by=exp]
				tmp$metric <- tmp$pos / tmp$nRuns
				
				dataMat[n,j] = dataMat[n,j] + mean(tmp$metric)
			}
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

plotTimeOnNodes <- function(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.2, 0.6, 0.1), 
		netTypes=c("full", "erdos-renyi", "barabasi-albert"),accuracy=0.7, acstdv=0.3, updates=c("theta-up", "theta-norm"),
		methods=c("best-acc", "conf-perfect"), colours=rainbow(12), yrange=c(1,7), truncated='true', singleRun=TRUE){
	combolength <- length(methods)*length(netTypes)
	methodscount <- length(methods)
	if ('best-acc' %in% methods){
		combolength <- combolength - (length(netTypes)-1)
		methodscount <- methodscount-1  
	}
	if ('erdos-renyi' %in% netTypes){ combolength <- combolength + (length(links)-1)*methodscount }
	if ('barabasi-albert' %in% netTypes){ combolength <- combolength + (length(edges)-1)*methodscount }
	combolength <- combolength*length(updates)
	offset <- 0.008*combolength
	xmin <- min(nodes_list)
	xmax <- max(nodes_list)
	xrange <- (xmax - xmin)
	plot(c(xmin-(xrange*offset),xmax+(xrange*offset)), yrange, type='n', xlab="Number of nodes", ylab="Number of iterations", xaxt = "n")
	axis(1, at=nodes_list)
	grid()
	if (combolength>length(colours)){
		print(paste("WARNING! - There are ", combolength, " methods to display and you're using only ", length(colours), " colours.", sep=""))
	}
	
	legTxt <- c()
	for (nodes in nodes_list){
		dataFilt <- c()
		dataList <- list()
		for (method in methods){
			for (net in netTypes){
				for (edge in edges){
					for (link in links){
						for (update in updates){
							if (method == "best-acc" && (net != netTypes[1] || link != links[1] || edge != edges[1]) ) { ## getting only best-acc for the first net-type (because it's the same for others)
								next
							}
							if (net != "erdos-renyi" && link != links[1]) { ## looping on linking-probablity on if net is erdos-renyi
								next
							}
							if (net != "erdos-renyi") {link <- 0.3}
							if (net != "barabasi-albert" && edge != edges[1]) { ## looping on edges on if net is barabasi-albert
								next
							}
							filename <- paste(prefix,"out_net-", net, "_nodes-",nodes,"_edges-",edge,
									"_link-",link,"_model-",method,"_up-",update,"_acc-",accuracy,"_acstdv-",acstdv,"_trunc-",truncated,".txt",sep="")
							#print(filename)
							data <- read.table(filename, header=T)
							if (singleRun){
								dataList <- append( dataList, list(data[ data$pos == nodes | data$neg == nodes, 'iter']) )
							} else {
								dt <- data.table(data)
								dataFilt <- c(dataFilt, dt[, list(metric=mean(iter[pos == nodes | neg == nodes]) ), by=exp][,'metric'])
							}
							if (nodes == nodes_list[1]) {
								ltxt <- paste(net,method,sep=" ")
								if (net == "erdos-renyi"){
									ltxt <- paste(ltxt, " p:",link,sep="")
								}
								if (net == "barabasi-albert"){
									ltxt <- paste(ltxt, " m:",edge,sep="")
								}
								if (length(updates)>1){
									ltxt <- paste(ltxt,update,sep=" ")
								}
								ltxt <- gsub("perfect", "weigh", ltxt)
								ltxt <- gsub("rand", "rule", ltxt)
								legTxt <- c(legTxt, ltxt)
							}
						}
					}
				}
			}
		}
		size <- xrange/100
		positions = nodes + (seq(0,combolength-1) - floor((combolength-1)/2))*size*1.5
		if (singleRun){ dataFilt <- dataList }
#		print(length(dataFilt))
#		print(positions)
		boxplot( dataFilt, at=positions, boxwex=size, add=T, col=colours[1:combolength], axes=F, outcex=0.5)
		
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

riffle <- function(a, b) { 
	mlab <- min(length(a), length(b)) 
	seqmlab <- seq(length=mlab) 
	c(rbind(a[seqmlab], b[seqmlab]), a[-seqmlab], b[-seqmlab]) 
} 

plotTimeTwoMethodsOnNodes <- function(prefix, nodes_list=seq(11, 31, 4), edges=seq(2, 10, 2), links=seq(0.3, 0.8, 0.1), 
		netTypes=c("erdos-renyi", "barabasi-albert"), accuracy=0.6, acstdv=0.12, updates=c("theta-up", "optim-up"),
		methods=c("best-acc", "conf-perfect"), colours=rainbow(12), yrange=c(1,7), truncated='true', singleRun=TRUE,
		ttest=TRUE){
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
	plot(c(xmin-(xrange*offset),xmax+(xrange*offset)), yrange, type='n', xlab="Number of nodes", ylab="Number of iterations", xaxt = "n")
	axis(1, at=nodes_list)
	grid()
	if (combolength/2>length(colours)){
		print(paste("WARNING! - There are ", combolength, " methods to display and you're using only ", length(colours), " colours.", sep=""))
	}
	
	legTxt <- c()
	for (nodes in nodes_list){
		dataFilt <- c()
		dataList <- list()
		pval <- c()
		for (net in netTypes){
			for (edge in edges){
				for (link in links){
					meth_i <- 0
					for (method in methods){
						meth_i <- meth_i +1
						update <- updates[meth_i]
						if (method == "best-acc" && (net != netTypes[1] || link != links[1] || edge != edges[1]) ) { ## getting only best-acc for the first net-type (because it's the same for others)
							next
						}
						if (net != "erdos-renyi" && link != links[1]) { ## looping on linking-probablity on if net is erdos-renyi
							next
						}
						if (net != "erdos-renyi") {link <- 0.3}
						if (net != "barabasi-albert" && edge != edges[1]) { ## looping on edges on if net is barabasi-albert
							next
						}
						filename <- paste(prefix,"out_net-", net, "_nodes-",nodes,"_edges-",edge,
								"_link-",link,"_model-",method,"_up-",update,"_acc-",accuracy,"_acstdv-",acstdv,"_trunc-",truncated,".txt",sep="")
						#print(filename)
						data <- read.table(filename, header=T)
						if (singleRun){
							#return( list(data[ data$pos == nodes | data$neg == nodes, 'iter'])  )
							#data[ data$pos != nodes & data$neg != nodes, 'iter'] = 21
							#dataList <- append( dataList, list(data[ , 'iter']) )
							dataList <- append( dataList, list(data[ data$pos == nodes | data$neg == nodes, 'iter']) )
						} else {
							dt <- data.table(data)
							dataFilt <- c(dataFilt, dt[, list(metric=mean(iter[pos == nodes | neg == nodes]) ), by=exp][,'metric'])
						}
						if (nodes == nodes_list[1]) {
							ltxt <- paste(net,method,sep=" ")
							if (net == "erdos-renyi"){
								ltxt <- paste(ltxt, " p:",link,sep="")
							}
							if (net == "barabasi-albert"){
								ltxt <- paste(ltxt, " m:",edge,sep="")
							}
							if (length(updates)>1){
								ltxt <- paste(ltxt,update,sep=" ")
							}
							ltxt <- gsub("perfect", "weigh", ltxt)
							ltxt <- gsub("rand", "rule", ltxt)
							legTxt <- c(legTxt, ltxt)
						}
						if (meth_i > 1 & ttest){
							indx <- length(dataList)
							#return(dataList)
							ttst <- t.test( x=data.frame(dataList[indx-1]), y = data.frame(dataList[indx]),
									alternative = "greater", mu = 0, paired = FALSE, var.equal = FALSE,
									conf.level = 0.95)
							#return(ttst)
							pval <- c(pval, ttst['p.value'])
						}
					}
				}
			}
		}
		size <- xrange/100
		positions = nodes + (seq(0,combolength-1) - floor((combolength-1)/2))*size*1.5
		if (singleRun){ dataFilt <- dataList }
#		print(length(dataFilt))
#		print(positions)
		cols <- riffle( colours[1:(combolength/2)], rep("grey", (combolength/2)) )
		cols_line <- riffle( rep("black", (combolength/2)), colours[1:(combolength/2)] )
#		cols_line <- "black"
		boxplot( dataFilt, at=positions, boxwex=size, add=T, col=cols, axes=F, outcex=0.5, border=cols_line)
		
		odd <- function(x) x%%2 != 0 
		significant <- pval < 0.001
		#points( (positions[odd(1:length(positions))]+size)*significant, y=rep(-0.1, length(positions)/2), pch='*', col=colours[1:(combolength/2)], cex=2 )
		#return(pval)
		
	}
	legend('topright', legTxt, cex=0.7, pt.bg=cols, col=cols_line, pch=22, pt.cex=1.1, lty=-1, bg='white')
}


