# Author: Andreagiovanni Reina - a.reina@sheffield.ac.uk - University of Sheffield
###################################################################################
library(data.table)

plotAll <-function(resdir){
	yaxes <- c("effectivity", "degree", "degree-scaled", "degree-stddev", "clustering", "clustering-stddev")
#	for (yaxis in yaxes){
#		pdf(paste("/Users/joefresna/DecisionsOnNetworks/results/deadlock/space-",yaxis,".pdf",sep=""))
#		plotEffectivityOnNetParam(prefix=resdir, nodes_list=seq(11,47,12), param=seq(0.25,0.5,0.05), xlabel="Communication range", netType="space", methods=c("conf-perfect"), link=0.2, 
#				updates=c("optim-up"), bxplt=F, combineMethods=T, legNodes=T, yaxis=yaxis, yrange=if(yaxis == "degree" ){c(0,20)}else{if(yaxis == "degree-stddev" ){c(0,5)}else{c(0,1)}})
#		dev.off()
#	}
#	for (yaxis in yaxes){
#		pdf(paste("/Users/joefresna/DecisionsOnNetworks/results/deadlock/ER-",yaxis,".pdf",sep=""))
#		plotEffectivityOnNetParam(prefix=resdir, nodes_list=c(11,23,35,47), param=seq(0.2,0.8,0.2), xlabel="Link probability", netType="erdos-renyi", methods=c("conf-perfect"), edge=0.20,
#				updates=c("optim-up"), bxplt=F, combineMethods=T, legNodes=T, yaxis=yaxis, yrange=if(yaxis == "degree" ){c(0,20)}else{if(yaxis == "degree-stddev" ){c(0,5)}else{c(0,1)}})
#		dev.off()
#	}
	for (yaxis in yaxes){
		pdf(paste("/Users/joefresna/DecisionsOnNetworks/results/deadlock/combo-",yaxis,".pdf",sep=""))
		plotEffCombo(prefix=resdir, nodes_list=c(11,23,35,47), netPar1=seq(0.2, 0.8, 0.2), netPar2=seq(0.2, 0.5, 0.05), netTypes=c("erdos-renyi","space"), accuracy=0.6, acstdv=0.12, 
				method="conf-perfect", T_MAX=50, update="optim-up", colours=rainbow(10), yaxis=yaxis, yrange=if(yaxis == "degree" ){c(0,20)}else{if(yaxis == "degree-stddev" ){c(0,5)}else{c(0,1)}})
		dev.off()
	}
}

plotEffectivityOnNetParam <- function(prefix, nodes_list=c(23), param=seq(2, 10, 4), xlabel="Communication range",numRuns=100,
		netType="barabasi-albert", accuracy=0.6, acstdv=0.12, methods=c("majority-rand", "conf-perfect"), edge=2, link=0.3, T_MAX=50,
		updates=c("no-up", "theta-up", "theta-norm"), colours=rainbow(10), yrange=c(0,1), bxplt=FALSE, truncated='true', combineMethods=FALSE, 
		legNodes=FALSE, yaxis="effectivity") {
	
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
	plot(c(xmin-(xrange*offset),xmax+(xrange*offset)), yrange, type='n', xlab=xlabel, ylab=simpleCap(yaxis), xaxt = "n")
	axis(1, at=param)
	grid()
	if (combolength>length(colours)){
		print(paste("WARNING! - There are ", combolength, " methods to display and you're using only ", length(colours), " colours.", sep=""))
	}
	
	legTxt <- c()
	dataPlot <- data.frame()
	for (netPar in param){
		dataFilt <- c()
		m <- 0
		for (nodes in nodes_list){
			if (netPar >= nodes) {next}
			for (method in methods){
				for (update in updates_loop){
					m <- m+1
					if (method == "best-acc" && (netType != netTypes[1] || link != links[1] || edge != edges[1]) ) { ## getting only best-acc for the first net-type (because it's the same for others)
						next
					}
					if (netType == "space") {
						filename <- paste(prefix,"out_net-", netType, "_nodes-",nodes,"_range-",format(netPar, nsmall=2),
								"_link-",link,"_model-",method,"_up-",update,"_acc-",accuracy,"_acstdv-",acstdv,"_eps-0.1.txt",sep="")
					} else {
						if (netType == "erdos-renyi") {link <- netPar}
						if (netType == "barabasi-albert") {edge <- netPar}
						if (!combineMethods) {
							update <- updates[m]
						}
						filename <- paste(prefix,"out_net-", netType, "_nodes-",nodes,"_range-",edge,
								"0_link-",link,"_model-",method,"_up-",update,"_acc-",accuracy,"_acstdv-",acstdv,"_eps-0.1.txt",sep="")
					}
					#print(filename)
					data <- read.table(filename, header=T)
					#print(ncol(dataPlot))
					if (bxplt){
						dt <- data.table(data)
						dataFilt <- c(dataFilt, dt[, list(eff=length(exp[ iter <= T_MAX & (pos == nodes | neg == nodes) ])/numRuns ), by=exp][,'eff'])
					} else {
						#print(netPar)
						#print(m)
						#print(nrow(data[ data$iter <= T_MAX & (data$pos == nodes | data$neg == nodes) , ]))
						#print(nrow(data))
						dataPlot <- rbind(dataPlot, c(netPar, m, nrow(data[ data$iter <= T_MAX & (data$pos == nodes | data$neg == nodes) , ]) / nrow(data) , mean(data$deg), mean(data$deg)/nodes, mean(data$dstd), mean(data$clust), mean(data$cstd) ))
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
		#print(length(dataFilt))
		#print(positions)
		if (bxplt){
			boxplot( dataFilt, at=positions, boxwex=size, add=T, col=colours[1:combolength], axes=F, outcex=0.5)
		}
#		else {
#			print( paste(netType, netPar) )
#			print( data.frame(lapply(dataFilt, mean, na.rm=TRUE)) )
#			points( positions, lapply(dataFilt, mean, na.rm=TRUE), pch=pntTypes[1:combolength], col=colours[1:combolength], cex=1, lwd=2 )
#		}
	}
	
	if (bxplt){
		legend('bottomright', legTxt, fill=colours[1:combolength], cex=1.5, bg='white')
	} else {
		
		colnames(dataPlot) <- c("netPar", "nodes", "effectivity", "degree", "degree-scaled", "degree-stddev", "clustering", "clustering-stddev")
		#print(colnames(dataPlot))
		#return(dataPlot)
		for (m in seq(1,length(legTxt))){
			points( dataPlot[dataPlot[,2]==m,1], dataPlot[dataPlot[,2]==m, yaxis], pch=pntTypes[m], col=colours[m], cex=1, lwd=2, type='b' )
		}
		
		if (legNodes){
			legend('bottomright', paste("Nodes",nodes_list), pch=pntTypes[1:length(nodes_list)], col=colours[1:length(nodes_list)], lty=0, cex=1, lwd=2, bg='white' )
		} else {
			legend('bottomright', legTxt, pch=pntTypes[1:combolength], col=colours[1:combolength], lty=0, cex=1, lwd=2, bg='white' )
		}
	}
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

simpleCap <- function(x) {
	s <- strsplit(x, " ")[[1]]
	paste(toupper(substring(s, 1,1)), substring(s, 2),
			sep="", collapse=" ")
}

