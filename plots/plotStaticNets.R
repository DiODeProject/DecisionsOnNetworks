# Author: Andreagiovanni Reina - a.reina@sheffield.ac.uk - University of Sheffield
###################################################################################
library(data.table)
library(MASS)
library(RColorBrewer)

setyrange <- function(yaxis){
	if(yaxis == "degree" ){
		rng <- c(0,50)
	}else{ if(yaxis == "degree-stddev" ){
			rng <- c(0,5)
		}else{ if(yaxis == "time" ){
				rng <- c(0,5)
			}else{ if(yaxis == "confidence" ){
					rng <- c(0,23)
			}else{ if(yaxis == "confidence-stddev" ){
					rng <- c(0,10)
				}else{ if(yaxis == "effectivity" ){
						rng <- c(0.,1)
#						rng <- c(0.,1)
					}else{if(yaxis == "success" ){
#							rng <- c(0.,1)
							rng <- c(0.8,1)
						}else{
					rng <- c(0,1)
				}}}}}}}
	return(rng)
}

plotAll <-function(resdir){
	pltdir <- "/Users/joefresna/DecisionsOnNetworks/results/synch_1-plt/"
	system(paste("mkdir -p ", pltdir))
	
	#yaxes <- c("effectivity", "degree", "degree-scaled", "degree-stddev", "success", "clustering", "clustering-stddev", "time")
	yaxes <- c("effectivity", "success", "time")
	darkcols <- brewer.pal(8, "Dark2")
	#yaxes <- c("time")
#	for (yaxis in yaxes){
#		rng <- setyrange(yaxis)
#		#for (speed in speeds){
#			pdf(paste("/Users/joefresna/DecisionsOnNetworks/results/static_lin2-plt/",yaxis,"_BA-on-edge.pdf",sep=""))
#			plotOnNetParam(prefix=resdir, nodes_list=seq(20, 100, 20), link_list=0.2, edge_list=seq(3, 12, 3), range_list=c(0), xlabel="Link edges",
#					netType="barabasi-albert", accuracy=0.6, acstdv=0.12, methods=c("conf-perfect"), pbound='false', T_MAX=1000, xpar='edge',
##					updates=c("finite-time"), colours=rainbow(5), yaxis=yaxis, yrange=rng, epsilon=0.05, legNodes=FALSE, env=1, quorum=1)
##					updates=c("optim-up", "finite-time"), colours=rainbow(5), yaxis=yaxis, yrange=rng, epsilon=0, legNodes=TRUE, env=1, quorum=1)
#					updates=c("optim-up", "belief-up"), colours=rainbow(5), yaxis=yaxis, yrange=rng, epsilon=0, legNodes=TRUE, env=1, quorum=1)
#			dev.off()
#			pdf(paste("/Users/joefresna/DecisionsOnNetworks/results/static_lin2-plt/",yaxis,"_ER-on-link.pdf",sep=""))
#			plotOnNetParam(prefix=resdir, nodes_list=seq(20, 100, 20), link_list=seq(0.2, 0.8, 0.2), edge_list=3, range_list=c(0), xlabel="Link probability",
#					netType="erdos-renyi", accuracy=0.6, acstdv=0.12, methods=c("conf-perfect"), pbound='false', T_MAX=1000, xpar='link',
##					updates=c("finite-time"), colours=rainbow(5), yaxis=yaxis, yrange=rng, epsilon=0.05, legNodes=FALSE, env=1, quorum=1)
##					updates=c("optim-up", "finite-time"), colours=rainbow(5), yaxis=yaxis, yrange=rng, epsilon=0, legNodes=TRUE, env=1, quorum=1)
#					updates=c("optim-up", "belief-up"), colours=rainbow(5), yaxis=yaxis, yrange=rng, epsilon=0, legNodes=TRUE, env=1, quorum=1)
#			dev.off()
#			pdf(paste("/Users/joefresna/DecisionsOnNetworks/results/static_lin2-plt/",yaxis,"_SP-on-range.pdf",sep=""))
#			plotOnNetParam(prefix=resdir, nodes_list=seq(20, 100, 20), link_list=0.2, edge_list=3, range_list=seq(0.20, 0.35, 0.05), xlabel="Communication range",
#					netType="space", accuracy=0.6, acstdv=0.12, methods=c("conf-perfect"), pbound='false', T_MAX=1000, xpar='range',
##					updates=c("finite-time"), colours=rainbow(5), yaxis=yaxis, yrange=rng, epsilon=0.05, legNodes=FALSE, env=1, quorum=1)
##					updates=c("optim-up", "finite-time"), colours=rainbow(5), yaxis=yaxis, yrange=rng, epsilon=0, legNodes=TRUE, env=1, quorum=1)
#					updates=c("optim-up", "belief-up"), colours=rainbow(5), yaxis=yaxis, yrange=rng, epsilon=0, legNodes=TRUE, env=1, quorum=1)
#			dev.off()
#		#}
#	}
	for (yaxis in yaxes){
		rng <- setyrange(yaxis)
		#for (speed in speeds){
			nodes_list <- seq(20, 100, 10)
			#if (yaxis == 'time') {nodes_list <- seq(20, 100, 20)}
			pdf(paste(pltdir,yaxis,"_BA-on-nodes.pdf",sep=""))
			res <- plotOnNodes(prefix=resdir, nodes_list=nodes_list, link_list=0.2, range_list=c(0), edge_list=c(3,9), #edge_list=if(yaxis == 'time'){c(3,12)}else{seq(3, 12, 3)}, 
					netType="barabasi-albert", accuracy=0.1, acstdv=0.20, methods=c("log-odds-distr"), pbound='false', T_MAX=1000, xpar='edge', boxptime=FALSE, xlabel="Nodes", 
#					updates=c("finite-time"), colours=rainbow(5), yaxis=yaxis, yrange=rng, epsilon=0.05, legNodes=FALSE, env=1, quorum=1)
#					updates=c("optim-up", "finite-time"), colours=rainbow(5), yaxis=yaxis, yrange=rng, epsilon=0, legNodes=TRUE, env=1, quorum=1)
					updates=c("optim-up", "no-up"), colours=darkcols, yaxis=yaxis, yrange=rng, epsilon=0, legNodes=TRUE, env=1, quorum=1)
			if (yaxis == yaxes[1]){ write.table(res, file = paste(resdir,"/finalTable_BA.txt",sep=""), sep="\t", col.names=TRUE, row.names=FALSE) }
			dev.off()
			pdf(paste(pltdir,yaxis,"_ER-on-nodes.pdf",sep=""))
			res <- plotOnNodes(prefix=resdir, nodes_list=nodes_list, edge_list=3, range_list=c(0), link_list=c(0.2,0.6),#link_list=if(yaxis == 'time'){c(0.2,0.8)}else{seq(0.2, 0.8, 0.2)}, 
					netType="erdos-renyi", accuracy=0.1, acstdv=0.20, methods=c("log-odds-distr"), pbound='false', T_MAX=1000, xpar='link', boxptime=FALSE, xlabel="Nodes",
#					updates=c("finite-time"), colours=rainbow(5), yaxis=yaxis, yrange=rng, epsilon=0.05, legNodes=FALSE, env=1, quorum=1)
#					updates=c("optim-up", "finite-time"), colours=rainbow(5), yaxis=yaxis, yrange=rng, epsilon=0, legNodes=TRUE, env=1, quorum=1)
					updates=c("optim-up", "no-up"), colours=darkcols, yaxis=yaxis, yrange=rng, epsilon=0, legNodes=TRUE, env=1, quorum=1)
			if (yaxis == yaxes[1]){ write.table(res, file = paste(resdir,"/finalTable_ER.txt",sep=""), sep="\t", col.names=TRUE, row.names=FALSE) }
			dev.off()
			pdf(paste(pltdir,yaxis,"_SP-on-nodes.pdf",sep=""))
			if (yaxis == 'time') rng <- c(0,30)
			res <- plotOnNodes(prefix=resdir, nodes_list=nodes_list, link_list=0.2, edge_list=3, range_list=0, #range_list=if(yaxis == 'time'){c(0.2,0.35)}else{seq(0.20, 0.35, 0.05)},
					netType="rgg-fixed-degree", accuracy=0.1, acstdv=0.20, methods=c("log-odds-distr"), pbound='false', T_MAX=1000, xpar='range', boxptime=FALSE, xlabel="Nodes", 
#					updates=c("finite-time"), colours=rainbow(5), yaxis=yaxis, yrange=rng, epsilon=0.05, legNodes=FALSE, env=1, quorum=1)
#					updates=c("optim-up", "finite-time"), colours=rainbow(5), yaxis=yaxis, yrange=rng, epsilon=0, legNodes=TRUE, env=1, quorum=1)
					updates=c("optim-up", "no-up"), colours=darkcols, yaxis=yaxis, yrange=rng, epsilon=0, legNodes=TRUE, env=1, quorum=1)
			if (yaxis == yaxes[1]){ write.table(res, file = paste(resdir,"/finalTable_SP.txt",sep=""), sep="\t", col.names=TRUE, row.names=FALSE) }
			dev.off()
		#}
	}
	return(0)
	
	rng <- c(0,20)
	yaxis <- "time"
	nodes_list <- seq(20, 100, 20)
	pdf(paste(pltdir,yaxis,"_BA-on-nodes_b.pdf",sep=""))
	plotOnNodes(prefix=resdir, nodes_list=nodes_list, link_list=0.2, edge_list=seq(3, 12, 3), range_list=c(0), xlabel="Nodes",
			netType="barabasi-albert", accuracy=0.6, acstdv=0.12, methods=c("conf-perfect"), pbound='false', T_MAX=1000, xpar='edge', boxptime=TRUE,
			updates=c("optim-up", "belief-up"), colours=darkcols, yaxis=yaxis, yrange=rng, epsilon=0, legNodes=TRUE, env=1, quorum=1)
	dev.off()
	#return(res)
	pdf(paste(pltdir,yaxis,"_ER-on-nodes_b.pdf",sep=""))
	plotOnNodes(prefix=resdir, nodes_list=nodes_list, link_list=seq(0.2, 0.8, 0.2), edge_list=3, range_list=c(0), xlabel="Nodes",
			netType="erdos-renyi", accuracy=0.6, acstdv=0.12, methods=c("conf-perfect"), pbound='false', T_MAX=1000, xpar='link', boxptime=TRUE,
			updates=c("optim-up", "belief-up"), colours=darkcols, yaxis=yaxis, yrange=rng, epsilon=0, legNodes=TRUE, env=1, quorum=1)
	dev.off()
	pdf(paste(pltdir,yaxis,"_SP-on-nodes_b.pdf",sep=""))
	plotOnNodes(prefix=resdir, nodes_list=nodes_list, link_list=0.2, edge_list=3, range_list=seq(0.20, 0.35, 0.05), xlabel="Nodes",
			netType="space", accuracy=0.6, acstdv=0.12, methods=c("conf-perfect"), pbound='false', T_MAX=1000, xpar='range', boxptime=TRUE,
			updates=c("optim-up", "belief-up"), colours=darkcols, yaxis=yaxis, yrange=rng, epsilon=0, legNodes=TRUE, env=1, quorum=1)
	dev.off()
	for (yaxis in yaxes){
		rng <- setyrange(yaxis)
		#for (speed in speeds){
			pdf(paste("/Users/joefresna/DecisionsOnNetworks/results/static_lin1-plt/",yaxis,"_BA-on-edge.pdf",sep=""))
			plotOnNetParam(prefix=resdir, nodes_list=seq(20, 100, 20), link_list=0.2, edge_list=seq(3, 12, 3), xlabel="Link edges",numRuns=1000,
					netType="barabasi-albert", accuracy=0.6, acstdv=0.12, methods=c("conf-perfect"), pbound='false', T_MAX=1000, xpar='edge',
#					updates=c("belief-up"), colours=rainbow(14), yaxis=yaxis, yrange=rng, epsilon=0.01, legNodes=FALSE, env=1, quorum=1)
					updates=c("optim-up", "belief-up", "finite-time"), colours=rainbow(15), yaxis=yaxis, yrange=rng, epsilon=0.01, legNodes=FALSE, env=1, quorum=1)
			dev.off()
			pdf(paste("/Users/joefresna/DecisionsOnNetworks/results/static_lin1-plt/",yaxis,"_ER-on-link.pdf",sep=""))
			plotOnNetParam(prefix=resdir, nodes_list=seq(20, 100, 20), link_list=seq(0.2, 0.8, 0.2), edge_list=3, xlabel="Link probability",numRuns=1000,
					netType="erdos-renyi", accuracy=0.6, acstdv=0.12, methods=c("conf-perfect"), pbound='false', T_MAX=1000, xpar='link',
#					updates=c("belief-up"), colours=rainbow(14), yaxis=yaxis, yrange=rng, epsilon=0.01, legNodes=FALSE, env=1, quorum=1)
					updates=c("optim-up", "belief-up", "finite-time"), colours=rainbow(15), yaxis=yaxis, yrange=rng, epsilon=0.01, legNodes=FALSE, env=1, quorum=1)
			dev.off()
		#}
	}
	return(0)
#	for (yaxis in yaxes){
#		rng <- setyrange(yaxis)
#		for (range in ranges){
#			pdf(paste("/Users/joefresna/DecisionsOnNetworks/results/dynamic-env033/",yaxis,"_r",format(range, nsmall=2),".pdf",sep=""))
#			plotEffectivityOnNetParam(prefix=resdir, nodes_list=seq(11,47,12), ranges=c(range), speeds=speeds, xlabel="Agents' speed", numRuns=1000,
#					netType="space", accuracy=0.6, acstdv=0.12, methods=c("conf-perfect"), pbound='true', T_MAX=200, xpar='speed',
#					updates=c("optim-up", "belief-up"), colours=rainbow(10), bxplt=FALSE, truncated='true', combineMethods=FALSE, epsilon=0.33,
#					legNodes=F, yaxis=yaxis, yrange=rng)
#			dev.off()
#		}
#	}
}

plotFiniteTime <-function(resdir){
	#yaxes <- c("effectivity", "degree", "degree-scaled", "degree-stddev", "success", "clustering", "clustering-stddev", "time")
	yaxes <- c("effectivity", "success", "time", "confidence")
	#yaxes <- c("time")
	for (yaxis in yaxes){
		rng <- setyrange(yaxis)
		#for (speed in speeds){
#		pdf(paste("/Users/joefresna/DecisionsOnNetworks/results/static_lin2-plt/",yaxis,"_BA-on-edge.pdf",sep=""))
#		plotOnNetParam(prefix=resdir, nodes_list=seq(20, 100, 20), link_list=0.2, edge_list=seq(3, 12, 3), range_list=c(0), xlabel="Link edges",
#				netType="barabasi-albert", accuracy=0.6, acstdv=0.12, methods=c("conf-perfect"), pbound='false', T_MAX=1000, xpar='edge',
#				updates=c("finite-time"), colours=rainbow(5), yaxis=yaxis, yrange=rng, epsilon=0, legNodes=TRUE, env=1, quorum=1)
#		dev.off()
		pdf(paste("/Users/joefresna/DecisionsOnNetworks/results/static_lin5-plt/",yaxis,"_ER-on-link.pdf",sep=""))
		plotOnNetParam(prefix=resdir, nodes_list=seq(20, 100, 20), link_list=seq(0.2, 0.8, 0.2), edge_list=3, range_list=c(0), xlabel="Link probability",
				netType="erdos-renyi", accuracy=0.6, acstdv=0.12, methods=c("conf-perfect"), pbound='false', T_MAX=1000, xpar='link',
				updates=c("finite-time"), colours=rainbow(5), yaxis=yaxis, yrange=rng, epsilon=0.9, legNodes=TRUE, env=1, quorum=1)
		dev.off()
#		pdf(paste("/Users/joefresna/DecisionsOnNetworks/results/static_lin2-plt/",yaxis,"_SP-on-range.pdf",sep=""))
#		plotOnNetParam(prefix=resdir, nodes_list=seq(20, 100, 20), link_list=0.2, edge_list=3, range_list=seq(0.20, 0.35, 0.05), xlabel="Communication range",
#				netType="space", accuracy=0.6, acstdv=0.12, methods=c("conf-perfect"), pbound='false', T_MAX=1000, xpar='range',
#				updates=c("optim-up", "belief-up"), colours=rainbow(5), yaxis=yaxis, yrange=rng, epsilon=0, legNodes=TRUE, env=1, quorum=1)
#		dev.off()
		#}
	}
	for (yaxis in yaxes){
		rng <- setyrange(yaxis)
		#for (speed in speeds){
		nodes_list <- seq(20, 100, 10)
		if (yaxis == 'time') {nodes_list <- seq(20, 100, 20)}
#		pdf(paste("/Users/joefresna/DecisionsOnNetworks/results/static_lin2-plt/",yaxis,"_BA-on-nodes.pdf",sep=""))
#		plotOnNodes(prefix=resdir, nodes_list=nodes_list, link_list=0.2, edge_list=seq(3, 12, 3), range_list=c(0), xlabel="Nodes",
#				netType="barabasi-albert", accuracy=0.6, acstdv=0.12, methods=c("conf-perfect"), pbound='false', T_MAX=1000, xpar='edge',
#				updates=c("finite-time"), colours=rainbow(5), yaxis=yaxis, yrange=rng, epsilon=0, legNodes=TRUE, env=1, quorum=1)
#		dev.off()
		pdf(paste("/Users/joefresna/DecisionsOnNetworks/results/static_lin5-plt/",yaxis,"_ER-on-nodes.pdf",sep=""))
		plotOnNodes(prefix=resdir, nodes_list=nodes_list, link_list=seq(0.2, 0.8, 0.2), edge_list=3, range_list=c(0), xlabel="Nodes",
				netType="erdos-renyi", accuracy=0.6, acstdv=0.12, methods=c("conf-perfect"), pbound='false', T_MAX=1000, xpar='link',
				updates=c("finite-time"), colours=rainbow(5), yaxis=yaxis, yrange=rng, epsilon=0.9, legNodes=TRUE, env=1, quorum=1)
		dev.off()
#		pdf(paste("/Users/joefresna/DecisionsOnNetworks/results/static_lin2-plt/",yaxis,"_SP-on-nodes.pdf",sep=""))
#		plotOnNodes(prefix=resdir, nodes_list=nodes_list, link_list=0.2, edge_list=3, range_list=seq(0.20, 0.35, 0.05), xlabel="Nodes",
#				netType="space", accuracy=0.6, acstdv=0.12, methods=c("conf-perfect"), pbound='false', T_MAX=1000, xpar='range',
#				updates=c("optim-up", "belief-up"), colours=rainbow(5), yaxis=yaxis, yrange=rng, epsilon=0, legNodes=TRUE, env=1, quorum=1)
#		dev.off()
		#}
	}
}

plotEpsilons <-function(resdir){
	#yaxes <- c("effectivity", "degree", "degree-scaled", "degree-stddev", "success", "clustering", "clustering-stddev", "time")
	yaxes <- c("effectivity", "success", "time")
	#yaxes <- c("time")
#	epsilons <- c(0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00)
	#epsilons <- seq(0.01, 0.1, 0.01)
	epsilons <- seq(0.001, 0.016, 0.003)
#	edges_list <- c(5, 10, 15, 20)
	edges_list <- seq(3, 12, 3)
	link_list <- c(0.2, 0.4, 0.6, 0.8)
	nodes_list <- seq(20, 100, 20) 
	#nodes_list=c(20, 30, 40, 50, 60, 70, 80, 100) #c(50, 100, 150, 200)
	for (yaxis in yaxes){
		rng <- setyrange(yaxis)
		res <- data.frame()
		for (edges in edges_list){
			pdf(paste("/Users/joefresna/DecisionsOnNetworks/results/static3-plt/",yaxis,"_BA-e",edges,".pdf",sep=""))
#			dp <- plotOnEpsilon(prefix=resdir, nodes_list=c(100), range=0.05, speed=0.05, xlabel="Epsilon", numRuns=100,
#					netType="space", accuracy=0.6, acstdv=0.12, methods=c("conf-perfect"), pbound='false', T_MAX=1000, xpar='range',
#					updates=c("optim-up", "belief-up", "finite-time"), colours=rainbow(10), bxplt=FALSE, truncated='true', combineMethods=FALSE, 
#					epsilons=epsilons, legNodes=F, yaxis=yaxis, yrange=rng, env=env, stddev=45)
			dp <- plotOnEpsilon(prefix=resdir, nodes_list=nodes_list, edges=edges, link=0.2, xlabel="Epsilon", numRuns=1000,
					netType="barabasi-albert", accuracy=0.6, acstdv=0.12, methods=c("conf-perfect"), pbound='false', T_MAX=1000, xpar='epsilon',
#					updates=c("optim-up", "belief-up", "finite-time"), colours=rainbow(15), bxplt=FALSE, truncated='true', combineMethods=FALSE, 
					updates=c("finite-time"), colours=rainbow(5), bxplt=FALSE, truncated='true', combineMethods=FALSE, 
					epsilons=epsilons, legNodes=F, yaxis=yaxis, yrange=rng, env=1, quorum=1)
			dev.off()
			#dp <- dp[dp[,'update']=='belief-up',]
			dp[,2] <- edges
			colnames(dp)[2] <- 'env'
			res <- rbind(res, dp)
		}
		for (link in link_list){
			pdf(paste("/Users/joefresna/DecisionsOnNetworks/results/static3-plt/",yaxis,"_ER-l",link,".pdf",sep=""))
			dp <- plotOnEpsilon(prefix=resdir, nodes_list=nodes_list, edges=3, link=link, xlabel="Epsilon", numRuns=1000,
					netType="erdos-renyi", accuracy=0.6, acstdv=0.12, methods=c("conf-perfect"), pbound='false', T_MAX=1000, xpar='epsilon',
#					updates=c("optim-up", "belief-up", "finite-time"), colours=rainbow(15), bxplt=FALSE, truncated='true', combineMethods=FALSE, 
					updates=c("finite-time"), colours=rainbow(5), bxplt=FALSE, truncated='true', combineMethods=FALSE, 
					epsilons=epsilons, legNodes=F, yaxis=yaxis, yrange=rng, env=1, quorum=1)
			dev.off()
			#dp <- dp[dp[,'update']=='belief-up',]
			dp[,2] <- link
			colnames(dp)[2] <- 'env'
			res <- rbind(res, dp)
		}
	}
	write.table(res, file = paste(resdir,"/finalTable.txt",sep=""), sep="\t", col.names=TRUE, row.names=FALSE)
	return(res)
#	for (yaxis in yaxes){
#		pdf(paste("/Users/joefresna/DecisionsOnNetworks/results/deadlock/ER-",yaxis,".pdf",sep=""))
#		plotEffectivityOnNetParam(prefix=resdir, nodes_list=c(11,23,35,47), param=seq(0.2,0.8,0.2), xlabel="Link probability", netType="erdos-renyi", methods=c("conf-perfect"), edge=0.20,
#				updates=c("optim-up"), bxplt=F, combineMethods=T, legNodes=T, yaxis=yaxis, yrange=if(yaxis == "degree" ){c(0,20)}else{if(yaxis == "degree-stddev" ){c(0,5)}else{c(0,1)}})
#		dev.off()
#	}
#	for (yaxis in yaxes){
#		pdf(paste("/Users/joefresna/DecisionsOnNetworks/results/deadlock/combo-",yaxis,".pdf",sep=""))
#		plotEffCombo(prefix=resdir, nodes_list=c(11,23,35,47), netPar1=seq(0.2, 0.8, 0.2), netPar2=seq(0.2, 0.5, 0.05), netTypes=c("erdos-renyi","space"), accuracy=0.6, acstdv=0.12, 
#				method="conf-perfect", T_MAX=50, update="optim-up", colours=rainbow(10), yaxis=yaxis, yrange=if(yaxis == "degree" ){c(0,20)}else{if(yaxis == "degree-stddev" ){c(0,5)}else{c(0,1)}})
#		dev.off()
#	}
}

plotOnNetParam <- function(prefix, nodes_list=seq(11,47,12), link_list=seq(0.2, 0.8, 0.2), edge_list=seq(5, 20, 5), range_list=c(1), xlabel="Link probability",
		netType="erdos-renyi", accuracy=0.6, acstdv=0.12, methods=c("conf-perfect"), pbound='true', T_MAX=50, xpar='range',
		updates=c("optim-up"), colours=rainbow(10), yrange=c(0,1), epsilon=0.1, legNodes=FALSE, yaxis="effectivity", env='none', quorum=1) {
	
	pntTypes = c(1,4,5,3,8,7,0,2,6,9,10,11,12,13)
	
	combolength <- length(updates)*length(nodes_list)
	if (xpar == 'link'){
		xpars <- link_list
	}
	if (xpar == 'edge'){
		xpars <- edge_list
	}
	if (xpar == 'range'){
		xpars <- range_list
	}
	
	xmin <- min(xpars)
	xmax <- max(xpars)
	xrange <- (xmax - xmin)
	if (yaxis == 'time' ){ offset <- 0.008*combolength} else {offset <- 0}
	ylabname <- simpleCap( gsub("effectivity", "Convergence", gsub("success", "Group accuracy", yaxis)) )
	plot(c(xmin-(xrange*offset), xmax+(xrange*offset)), yrange, type='n', xlab=xlabel, ylab=ylabname, xaxt = "n")
	axis(1, at=xpars)
	grid()
	legTxt <- c()
	collist <- c()
	collistB <- c()
	dataPlot <- data.frame()
#	for (edge in edge_list){
#	for (link in link_list){
#	for (range in range_list){
	for (netPar in xpars){
		dataFilt <- c()
		m <- 0
		for (nodes in nodes_list){
			for (method in methods){
				for (update in updates){
					m <- m+1
					#						filename <- paste(prefix,"out_net-", netType, "_nodes-",nodes,"_edges-",edge,"_link-",link,"_bound-",pbound,
					#								if(env!='none'){paste("_env-",env,sep="")}else{""},"_model-",method,"_up-",update,"_acc-",accuracy,"_acstdv-",acstdv,"_eps-",
					#								format(epsilon, nsmall=2),".txt",sep="")
#					if  (netType == 'erdos-renyi'){
#						netPar <- link
#					}
#					if  (netType == 'barabasi-albert'){
#						netPar <- edge
#					}
#					if  (netType == 'space'){
#						netPar <- range
#					}
					filename <- paste(prefix,"out_net-", netType, "_nodes-",nodes,"_link-",if(netType == 'space'){format(netPar,nsmall=2)}else{netPar},"_bound-",pbound,
							if(env!='none'){paste("_env-",env,sep="")}else{""},"_model-",method,"_up-",update,"_acc-",accuracy,"_acstdv-",acstdv,"_eps-",
							epsilon,".txt",sep="")
					
					#print(filename)
					data <- read.table(filename, header=T)
					#						if (xpar == 'link'){
					#							xparv <- link
					#							#legv <- speed
					#						} else {
					#							xparv <- edge
					#							#legv <- range
					#						}
					effectivity <- nrow(data[ data$iter <= T_MAX & (data$pos >= (nodes*quorum) | data$neg >= (nodes*quorum) ) , ]) / nrow(data)
					success <- nrow(data[ data$pos >= (nodes*quorum), ]) / nrow(data)
					times <- data[ data$pos == nodes | data$neg == nodes, 'iter']
					dataPlot <- rbind(dataPlot, c(netPar, m, which(update == updates)[1], nodes, effectivity,
									success, mean(times), mean(data$deg), mean(data$deg)/nodes, mean(data$dstd), mean(data$clust), 
									mean(data$cstd), mean(data$conf), mean(data$confsd), 
									1.96*sqrt( effectivity * (1-effectivity) / nrow(data) ),
									1.96*sqrt( success * (1-success) / nrow(data) ),
									1.96*sd(times)/sqrt(length(times)) ))
					
					if (yaxis == 'time' ){
						dataFilt <- append( dataFilt, list(data[ data$pos == nodes | data$neg == nodes, 'iter']) )
					}
					
					if (netPar == xpars[1]) {
						#ltxt <- paste(update," ",substr(xpar,1,1),":",legv," n:",nodes,sep="")
						ltxt <- paste(update," n:",nodes,sep="")
						ltxt <- gsub("perfect", "weigh", ltxt)
						#ltxt <- gsub("rand", "rule", ltxt)
						#ltxt <- gsub("optim-up", "with-up", ltxt)
						legTxt <- c(legTxt, ltxt)
						if (legNodes){
							if (m %% 2 == 0){
								collist <- append(collist, 'white' )
								collistB <- append(collistB, colours[m/length(updates)] )
							} else {
								collist <- append(collist, colours[(m+1)/length(updates)])
								collistB <- append(collistB, 'black')
							}
						}
					}
				}
			}
		}
		if (yaxis == 'time' ){
			size <- xrange/80
			positions = netPar + (seq(0,combolength-1) - floor((combolength-1)/2))*size*1.5
			#print(positions)
			if (legNodes){
				boxplot( dataFilt, at=positions, boxwex=size, add=T, col=collist, axes=F, outcex=0.5, border=collistB, lwd=c(1,1.5), notch=TRUE)
			} else {
				boxplot( dataFilt, at=positions, boxwex=size, add=T, col=colours[1:combolength], axes=F, outcex=0.5)
			}
			
			
		}
	}
	colnames(dataPlot) <- c("netPar", "loop", "update", "nodes", "effectivity", "success", "time", "degree", "degree-scaled", "degree-stddev", "clustering", 
			"clustering-stddev", "confidence", "confidence-stddev", "effectivity-stddev", "success-stddev", "time-stddev")
	#print(colnames(dataPlot))
	if (yaxis != 'time' ){
		for (m in seq(1,length(legTxt))){
			if (legNodes){
				points( dataPlot[dataPlot[,2]==m,1], dataPlot[dataPlot[,2]==m, yaxis], pch=pntTypes[(m+(m%%2))/length(updates)], col=colours[(m+(m%%2))/length(updates)], cex=1, lwd=2, type='b', lty=(((m-1) %% length(updates))+1) )
			} else {
				points( dataPlot[dataPlot[,2]==m,1], dataPlot[dataPlot[,2]==m, yaxis], pch=pntTypes[m], col=colours[m], cex=1, lwd=2, type='b', lty=m)
			}
			if (yaxis=="confidence" || yaxis=="success"){
				arrows(dataPlot[dataPlot[,2]==m,1], dataPlot[dataPlot[,2]==m, yaxis]-dataPlot[dataPlot[,2]==m, paste(yaxis,"-stddev",sep="")], 
						dataPlot[dataPlot[,2]==m,1], dataPlot[dataPlot[,2]==m, yaxis]+dataPlot[dataPlot[,2]==m, paste(yaxis,"-stddev",sep="")], 
						length=0.05, angle=90, code=3, col=colours[(m+(m%%2))/length(updates)], lty=(((m-1) %% length(updates))+1) )
			}
		}
	}
	
	legpos <- 'bottomright'
	legpos2 <- 'bottom'
	if (yaxis == 'time' ){ legpos <- 'topright'; legpos2 <- 'top' }
	if (yaxis == 'clustering' || yaxis == 'deegree' || yaxis == 'degree-scaled'){ legpos <- 'topleft'; legpos2 <- 'top'  }
	if (legNodes){
		if (yaxis == 'time' ){
			legend(legpos, gsub("belief-up", "Belief consensus", gsub("optim-up", "Bayes update", updates)), pch=22, col=c('black','gray'), pt.bg=c('gray','white'), lty=-1, cex=1, lwd=2, bg='white' )
		} else {			
			legend(legpos, gsub("belief-up", "Belief consensus", gsub("optim-up", "Bayes update", updates)), pch=-1, col='black', lty=seq(1,length(updates)), cex=1, lwd=2, bg='white' )
		}
		txtL <- paste("Nodes: ", nodes_list, sep="")
		legend(legpos2, txtL, pch=pntTypes[1:(combolength/length(updates))], col=colours[1:(combolength/length(updates))], lty=-1, cex=1, lwd=2, bg='white' )
	} else {
		legend(legpos, legTxt, pch=pntTypes[1:combolength], col=colours[1:combolength], lty=0, cex=1, lwd=2, bg='white' )
	}
	#return(dataPlot)
}

plotOnNodes <- function(prefix, nodes_list=seq(11,47,12), link_list=seq(0.2, 0.8, 0.2), edge_list=seq(5, 20, 5), range_list=c(1), xlabel="Link probability",
		netType="erdos-renyi", accuracy=0.6, acstdv=0.12, methods=c("conf-perfect"), pbound='true', T_MAX=50, xpar='range', boxptime=FALSE,
		updates=c("optim-up"), colours=rainbow(10), yrange=c(0,1), epsilon=0.1, legNodes=FALSE, yaxis="effectivity", env='none', quorum=1) {
	
	pntTypes = c(1,4,5,3,8,7,0,2,6,9,10,11,12,13)
	
	combolength <- length(updates)*length(link_list)*length(edge_list)*length(range_list)
	
	xmin <- min(nodes_list)
	xmax <- max(nodes_list)
	xrange <- (xmax - xmin)
	if (yaxis == 'time' ){ offset <- 0.008*combolength} else {offset <- 0}
	ylabname <- simpleCap( gsub("effectivity", "Convergence", gsub("success", "Group accuracy", gsub("time", "Timesteps", yaxis))) )
	plot(c(xmin-(xrange*offset), xmax+(xrange*offset)), yrange, type='n', xlab=xlabel, ylab=ylabname, xaxt = "n")
	axis(1, at=nodes_list)
	grid()
	legTxt <- c()
	collist <- c()
	collistB <- c()
	dataPlot <- data.frame()
	for (nodes in nodes_list){
		dataFilt <- c()
		m <- 0
		for (edge in edge_list){
			for (link in link_list){
				for (range in range_list){
					for (method in methods){
						for (update in updates){
							m <- m+1
							if  (netType == 'erdos-renyi'){
								netPar <- link
							}
							if  (netType == 'barabasi-albert'){
								netPar <- edge
							}
							if  (netType == 'space' || netType == 'rgg-fixed-degree'){
								netPar <- range
							}
#							filename <- paste(prefix,"out_net-", netType, "_nodes-",nodes,"_link-",if(netType == 'space'){format(netPar,nsmall=2)}else{netPar},"_bound-",pbound,
#									if(env!='none'){paste("_env-",env,sep="")}else{""},"_model-",method,"_up-",update,"_acc-",accuracy,"_acstdv-",acstdv,"_eps-",
#									epsilon,".txt",sep="")
							filename <- paste(prefix,"out_net-", netType, "_nodes-",nodes,"_link-",if(netType == 'space'){format(netPar,nsmall=2)}else{netPar},
									"_model-",method,"_up-",update,"_driftbase-",accuracy,"_driftrange-",format(acstdv,nsmall=2),".txt",sep="")
							print(filename)
							data <- read.table(filename, header=T)
							xparv <- nodes
							effectivity <- nrow(data[ data$iter <= T_MAX & (data$pos >= (nodes*quorum) | data$neg >= (nodes*quorum) ) , ]) / nrow(data)
							success <- nrow(data[ data$pos >= (nodes*quorum), ]) / nrow(data)
							times <- data[ data$pos == nodes | data$neg == nodes, 'iter']
							dataPlot <- rbind(dataPlot, c(netPar, m, which(update == updates)[1], nodes, effectivity,
											success, mean(times), mean(data$deg), mean(data$deg)/nodes, mean(data$dstd), mean(data$clust), 
											mean(data$cstd), mean(data$conf), mean(data$confsd), 
											1.96*sqrt( effectivity * (1-effectivity) / nrow(data) ),
											1.96*sqrt( success * (1-success) / nrow(data) ),
											1.96*sd(times)/sqrt(length(times)),
											effectivity - 1.96*sqrt( effectivity * (1-effectivity) / nrow(data) ), effectivity + 1.96*sqrt( effectivity * (1-effectivity) / nrow(data) ),
											success - 1.96*sqrt( success * (1-success) / nrow(data) ), success + 1.96*sqrt( success * (1-success) / nrow(data) ),
											mean(times) - 1.96*sd(times)/sqrt(length(times)), mean(times) + 1.96*sd(times)/sqrt(length(times))
											))
							
							if (yaxis == 'time' && boxptime){
								dataFilt <- append( dataFilt, list(data[ data$pos == nodes | data$neg == nodes, 'iter']) )
							}
							
							if (xparv == nodes_list[1]) {
								#ltxt <- paste(update," ",substr(xpar,1,1),":",legv," n:",nodes,sep="")
								ltxt <- paste(update," p:",netPar,sep="")
								ltxt <- gsub("perfect", "weigh", ltxt)
								#ltxt <- gsub("rand", "rule", ltxt)
								#ltxt <- gsub("optim-up", "with-up", ltxt)
								legTxt <- c(legTxt, ltxt)
								if (legNodes){
									if (m %% 2 == 0){
										collist <- append(collist, 'white' )
										collistB <- append(collistB, colours[m/length(updates)] )
									} else {
										collist <- append(collist, colours[(m+1)/length(updates)])
										collistB <- append(collistB, 'black')
									}
								}
							}
						}
					}
				}
			}
		}
		if (yaxis == 'time' && boxptime){
			size <- xrange/80
			positions = xparv + (seq(0,combolength-1) - floor((combolength-1)/2))*size*1.5
			#print(positions)
			if (legNodes){
				boxplot( dataFilt, at=positions, boxwex=size, add=T, col=collist, axes=F, outcex=0.5, border=collistB, lwd=c(1,1.5), notch=TRUE)
			} else {
				boxplot( dataFilt, at=positions, boxwex=size, add=T, col=colours[1:combolength], axes=F, outcex=0.5)
			}
		}
	}
	#print(dataPlot)
	colnames(dataPlot) <- c("netPar", "loop", "update", "nodes", "effectivity", "success", "time", "degree", "degree-scaled", "degree-stddev", "clustering", 
			"clustering-stddev", "confidence", "confidence-stddev", "effectivity-stddev", "success-stddev", "time-stddev",
			"emin", "epl", "smin", "spl", "tmin", "tpl" )
	#print(colnames(dataPlot))
	if (yaxis != 'time' || !boxptime){
		for (m in seq(1,length(legTxt))){
			if (legNodes){
				points( dataPlot[dataPlot[,2]==m,'nodes'], dataPlot[dataPlot[,2]==m, yaxis], pch=pntTypes[(m+(m%%2))/length(updates)], col=colours[(m+(m%%2))/length(updates)], cex=1, lwd=2, type='b', lty=(((m-1) %% length(updates))+1) )
			} else {
				points( dataPlot[dataPlot[,2]==m,'nodes'], dataPlot[dataPlot[,2]==m, yaxis], pch=pntTypes[m], col=colours[m], cex=1, lwd=2, type='b', lty=m )
			}
			if (yaxis=="confidence" || yaxis=="success" || yaxis=="time" || yaxis=="effectivity"){
				arrows(dataPlot[dataPlot[,2]==m,'nodes'], dataPlot[dataPlot[,2]==m, yaxis]-dataPlot[dataPlot[,2]==m, paste(yaxis,"-stddev",sep="")], 
						dataPlot[dataPlot[,2]==m,'nodes'], dataPlot[dataPlot[,2]==m, yaxis]+dataPlot[dataPlot[,2]==m, paste(yaxis,"-stddev",sep="")], 
						length=0.05, angle=90, code=3, col=colours[(m+(m%%2))/length(updates)] ) #lty=(((m-1) %% length(updates))+1))
			}
		}
	}
	
	legpos <- 'bottomright'
	legpos2 <- 'bottom'
	if (yaxis == 'time' ){ legpos <- 'topright'; legpos2 <- 'top' }
	if (yaxis == 'clustering' || yaxis == 'deegree' || yaxis == 'degree-scaled'){ legpos <- 'topleft'; legpos2 <- 'top'  }
	if (legNodes){
		#legend(legpos, paste("Nodes",nodes_list), pch=pntTypes[1:length(nodes_list)], col=colours[1:length(nodes_list)], lty=0, cex=1, lwd=2, bg='white' )
		if (yaxis == 'time' && boxptime){
			legend(legpos, gsub("belief-up", "Belief consensus", gsub("optim-up", "Bayes update", updates)), pch=22, col=c('black','gray'), pt.bg=c('gray','white'), lty=-1, cex=1, lwd=2, bg='white' )
		} else {
			legend(legpos, gsub("belief-up", "Belief consensus", gsub("optim-up", "Bayes update", updates)), pch=-1, col='black', lty=seq(1,length(updates)), cex=1, lwd=2, bg='white' )			
		}
		if  (netType == 'erdos-renyi'){
			txtL <- paste("Link P: ", link_list, sep="")
		}
		if  (netType == 'barabasi-albert'){
			txtL <- paste("Edges: ", edge_list, sep="")
		}
		if  (netType == 'space' || netType == 'rgg-fixed-degree'){
			txtL <- paste("Range: ", range_list, sep="")
		}
		legend(legpos2, txtL, pch=pntTypes[1:(combolength/length(updates))], col=colours[1:(combolength/length(updates))], lty=-1, cex=1, lwd=2, bg='white' )
	} else {
		legend(legpos, legTxt, pch=pntTypes[1:combolength], col=colours[1:combolength], lty=0, cex=1, lwd=2, bg='white' )
	}
	return(dataPlot)
}

plotOnEpsilon <- function(prefix, nodes_list=seq(11,47,12), range=0.05, edges=5, link=0.2, xlabel="Epsilon", numRuns=100,
		netType="space", accuracy=0.6, acstdv=0.12, methods=c("conf-perfect"), pbound='false', T_MAX=50, xpar='epsilon',
		updates=c("optim-up"), colours=rainbow(10), yrange=c(0,1), bxplt=FALSE, truncated='true', combineMethods=FALSE, epsilons=c(0.1),
		legNodes=FALSE, yaxis="effectivity", env='none', quorum=1) {
	
	pntTypes = c(1,4,5,3,8,7,0,2,6,9,10,11,12,13)
	
	combolength <- length(updates)*length(nodes_list)
	if (xpar == 'epsilon'){
		xpars <- epsilons
	} else {
		xpars <- epsilons
	}
	
	xmin <- min(xpars)
	xmax <- max(xpars)
	xrange <- (xmax - xmin)
	if (yaxis == 'time' ){ offset <- 0.008*combolength} else {offset <- 0}
	plot(c(xmin-(xrange*offset), xmax+(xrange*offset)), yrange, type='n', xlab=xlabel, ylab=simpleCap(yaxis), xaxt = "n")
	axis(1, at=xpars)
	grid()
	legTxt <- c()
	dataPlot <- data.frame()
	for (epsilon in epsilons){
		dataFilt <- c()
		m <- 0
		for (nodes in nodes_list){
			for (method in methods){
				for (update in updates){
					m <- m+1
					if  (netType == 'erdos-renyi'){
						netPar <- link
					}
					if  (netType == 'barabasi-albert'){
						netPar <- edges
					}
					if  (netType == 'space'){
						netPar <- range
					}
					filename <- paste(prefix,"out_net-", netType, "_nodes-",nodes,"_link-",if(netType == 'space'){format(netPar,nsmall=2)}else{netPar},"_bound-",pbound,
							if(env!='none'){paste("_env-",env,sep="")}else{""},"_model-",method,"_up-",update,"_acc-",accuracy,"_acstdv-",acstdv,"_eps-",
							format(epsilon, nsmall=3),".txt",sep="")
#					filename <- paste(prefix,"out_net-", netType, "_nodes-",nodes,"_edges-",edges,"_link-",link,"_bound-",pbound,
#							if(env!='none'){paste("_env-",env,sep="")}else{""},"_model-",method,"_up-",update,"_acc-",accuracy,"_acstdv-",acstdv,"_eps-",
#							format(if(update=='optim-up'){epsilons[1]}else{epsilon}, nsmall=2),".txt",sep="")
					#print(filename)
					data <- read.table(filename, header=T)
					xparv <- epsilon
					dataPlot <- rbind(dataPlot, c(xparv, m, which(update == updates)[1], nodes, nrow(data[ data$iter <= T_MAX & (data$pos >= (nodes*quorum) | data$neg >= (nodes*quorum) ) , ]) / nrow(data) ,
									nrow(data[ data$pos >= (nodes*quorum), ]) / nrow(data), mean(data$deg), mean(data$deg)/nodes, mean(data$dstd), mean(data$clust), 
									mean(data$cstd), mean(data$conf), mean(data$confsd) ))
					
					if (yaxis == 'time' ){
						dataFilt <- append( dataFilt, list(data[ data$pos >= (nodes*quorum) | data$neg >= (nodes*quorum), 'iter']) )
					}
					
					if (xparv == xpars[1]) {
						#ltxt <- paste(update," ",substr(xpar,1,1),":",legv," n:",nodes,sep="")
						ltxt <- paste(update," n:",nodes,sep="")
						ltxt <- gsub("perfect", "weigh", ltxt)
						#ltxt <- gsub("rand", "rule", ltxt)
						#ltxt <- gsub("optim-up", "with-up", ltxt)
						legTxt <- c(legTxt, ltxt)
					}
				}
			}
		}
		size <- xrange/80
		positions = xparv + (seq(0,combolength-1) - floor((combolength-1)/2))*size*1.5
		#print(positions)
		if (yaxis == 'time' ){
			boxplot( dataFilt, at=positions, boxwex=size, add=T, col=colours[1:combolength], axes=F, outcex=0.5)
		}
	}
	colnames(dataPlot) <- c("netPar", "loop", "update", "nodes", "effectivity", "success", "degree", "degree-scaled", "degree-stddev", "clustering", 
			"clustering-stddev", "confidence", "confidence-stddev")
	#print(colnames(dataPlot))
	#return(dataPlot)
	if (yaxis != 'time' ){
		for (m in seq(1,length(legTxt))){
			points( dataPlot[dataPlot[,2]==m,1], dataPlot[dataPlot[,2]==m, yaxis], pch=pntTypes[m], col=colours[m], cex=1, lwd=2, type='b', lty=((m %% length(updates))+1) )
			if (yaxis=="confidence"){
				arrows(dataPlot[dataPlot[,2]==m,1], dataPlot[dataPlot[,2]==m, yaxis]-dataPlot[dataPlot[,2]==m, paste(yaxis,"-stddev",sep="")], 
						dataPlot[dataPlot[,2]==m,1], dataPlot[dataPlot[,2]==m, yaxis]+dataPlot[dataPlot[,2]==m, paste(yaxis,"-stddev",sep="")], 
						length=0.05, angle=90, code=3, col=colours[m])
			}
		}
	}
	
	if (legNodes){
		legend('bottomright', paste("Nodes",nodes_list), pch=pntTypes[1:length(nodes_list)], col=colours[1:length(nodes_list)], lty=0, cex=1, lwd=2, bg='white' )
	} else {
		legpos <- 'bottomright'
		if (yaxis == 'time' ){ legpos <- 'topright' }
		if (yaxis == 'clustering' || yaxis == 'deegree'){ legpos <- 'topleft' }
		legend(legpos, legTxt, pch=pntTypes[1:combolength], col=colours[1:combolength], lty=0, cex=1, lwd=2, bg='white' )
	}
	return(dataPlot)
}


plotEffCombo <- function(prefix, nodes_list=c(23), netPar1=seq(0.2, 0.8, 0.2), netPar2=seq(0.25, 0.5, 0.05), netTypes=c("erdos-renyi","space"), accuracy=0.6, acstdv=0.12, 
		method="conf-perfect", T_MAX=50, update="optim-up", colours=rainbow(10), yrange=c(0,1), yaxis="effectivity") {
	
	pntTypes = c(1,4,5,3,8,7,0,2,6,9,10,11,12,13)
	
	xmin <- min(netPar1)
	xmax <- max(netPar1)
	xrange <- (xmax - xmin)
	if(netTypes[1]=="space"){
		xlabel1 <- "Communication Range"
		xlabel2 <- "Link Probability"
	}else{
		xlabel1 <- "Link Probability"
		xlabel2<- "Communication Range"
	}
	plot(c(xmin,xmax), yrange, type='n', xlab=xlabel1, ylab=simpleCap(yaxis), xaxt = "n")
	axis(1, at=netPar1)
#	axis(3, at=netPar1[1] + (netPar1[length(netPar1)] - netPar1[1]) * (netPar2 - netPar2[1])/(netPar2[length(netPar2)] - netPar2[1]), labels=netPar2 )
	axis(3, at=netPar2)
	mtext(side = 3, line = 3, xlabel2)
	grid()
	
	for (netType in netTypes){
		dataPlot <- data.frame()
		dataFilt <- c()
		for (np1 in netPar1){
			for (np2 in netPar2){
				if ( (netType == netTypes[1] && np2 != netPar2[1]) || (netType == netTypes[2] && np1 != netPar1[1]) ) {
					next
				}
				m <- 0
				for (nodes in nodes_list){
					m <- m+1
					if (netTypes[1] == "space"){
						rangepar <- np1
						linkpar <- np2
					} else {
						rangepar <- np2
						linkpar <- np1
					}
					if (netType == netTypes[1]){
						netPar <- np1
					} else {
						netPar <- netPar1[1] + (netPar1[length(netPar1)] - netPar1[1]) * (np2 - netPar2[1])/(netPar2[length(netPar2)] - netPar2[1])
						netPar <- np2
					}
					
					filename <- paste(prefix,"out_net-", netType, "_nodes-",nodes,"_range-",format(rangepar, nsmall=2),
							"_link-",linkpar,"_model-",method,"_up-",update,"_acc-",accuracy,"_acstdv-",acstdv,"_eps-0.1.txt",sep="")
					#print(filename)
					data <- read.table(filename, header=T)
					dataPlot <- rbind(dataPlot, c(netPar, m, nrow(data[ data$iter <= T_MAX & (data$pos == nodes | data$neg == nodes) , ]) / nrow(data) , mean(data$deg), mean(data$deg)/nodes, mean(data$dstd), mean(data$clust), mean(data$cstd) ))
				}
			}
		}
		colnames(dataPlot) <- c("netPar", "nodes", "effectivity", "degree", "degree-scaled", "degree-stddev", "clustering", "clustering-stddev")
		#print(colnames(dataPlot))
		#return(dataPlot)
		for (m in seq(1,length(nodes_list))){
			points( dataPlot[dataPlot[,2]==m,1], dataPlot[dataPlot[,2]==m, yaxis], pch=pntTypes[m], col=colours[m], cex=1, lwd=2, type='b', lty=((netType == netTypes[2])+1) )
		}
	}
	
	legend('bottomright', paste("Nodes",nodes_list), pch=pntTypes[1:length(nodes_list)], col=colours[1:length(nodes_list)], lty=0, cex=1, lwd=2, bg='white' )
	legend('bottomleft', netTypes, pch=-1, col='black', lty=c(1,2), cex=1, lwd=2, bg='white' )
	
}

repeatConfPlots <- function(n=1){
	for (i in seq(1,n)){
		confidenceOverTime(paste("/Users/joefresna/DecisionsOnNetworks/data/conf-log",i,".txt",sep=""), step=5, pdfout=paste("/Users/joefresna/DecisionsOnNetworks/results/dynamic-env/conf-",i,".pdf",sep=""))
	}
}

confidenceOverTime <- function(filename, step=1, pdfout='none'){
	if (pdfout != 'none'){
		pdf(pdfout)
	}
	data <- read.table(filename, header=F)
	plot( c(1,max(data[,1])), c(0,max(data[,-c(1,2)])), type='n', xlab="timestep", ylab="confidence" )
	axis(4, at=seq(0,1,0.2) * max(data[,-c(1,2)]), labels=seq(0,1,0.2), col='blue' )
	lines( data[,1], data[,2]*max(data[,-c(1,2)]), lwd=2, col='blue' )
	for (t in seq(1,nrow(data), step )){
		bplt <- boxplot( as.vector(data[data[,1]==t,-c(1,2)], mode="numeric"), at=t, boxwex=step, axes=F, add=T, outline=FALSE )
		points( rep(t,length(bplt$out)), bplt$out, type="p", pch=1, cex=0.1)
	}
	if (pdfout != 'none'){		dev.off()	}
}


simpleCap <- function(x) {
	s <- strsplit(x, " ")[[1]]
	paste(toupper(substring(s, 1,1)), substring(s, 2),
			sep="", collapse=" ")
}

