library(data.table)

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

plotEffectivity <- function(prefix, nodes=20, links=seq(0.2,0.8,0.2), drift=0.1, numRuns=100, methods=c("conf","M-rand","M-bias","M-inhib"), colours=rainbow(10)){
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

plotEffectivityOnNodes <- function(prefix, nodes_list=c(100), link=0.2, drift=0.1, numRuns=100, methods=c("conf","M-rand","M-bias","M-inhib"), colours=rainbow(10)){
	offset <- 0.02*length(methods)
	xmin <- min(nodes_list)
	xmax <- max(nodes_list)
	xrange <- (xmax - xmin)
	offset <- 0.08
	plot(c(xmin-(xrange*offset),xmax+(xrange*offset)),c(0,1), type='n', xlab="Number of nodes", ylab="Effectivity")
	grid()
	for (nodes in nodes_list){
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
		positions = nodes + (seq(0,length(methods)-1) - floor((length(methods)-1)/2))*size*1.5
		boxplot( dataFilt, at=positions, boxwex=size, add=T, col=colours[1:length(methods)], axes=F)
		
		#aggr <- dt[, list(eff=length(exp[pos == nodes | neg == nodes])/100 ), by=exp]
		#dt[, list(pos=sum(pos)/(nodes*100), neg=sum(neg)/(nodes*100)), by=exp]
	}
	legend('bottomright', methods, fill=colours[1:length(methods)], cex=0.5, bg='white')
	#return (data)
}

#############################################################################
################################# SUCCESS ###################################
#############################################################################

plotSuccess <- function(prefix, nodes=20, links=seq(0.2,0.8,0.2), drift=0.1, numRuns=100, 
		methods=c("conf","M-rand","M-bias","M-inhib"), const_methods=c("bestDDM", "confDDM", "FNconf", "FNmaj"),
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

plotSuccessOnNodes <- function(prefix, nodes_list=c(100, 500, 1000), link=0.2, drift=0.1, numRuns=100, methods=c("conf","M-rand","M-bias","M-inhib","bestDDM", "confDDM"), colours=rainbow(10)){
	offset <- 0.02*length(methods)
	xmin <- min(nodes_list)
	xmax <- max(nodes_list)
	xrange <- (xmax - xmin)
	plot(c(xmin-(xrange*offset),xmax+(xrange*offset)),c(0,1), type='n', xlab="Number of nodes", ylab="Success rate")
	grid()
	if (length(c(methods))>length(colours)){
		print(paste("WARNING! - There are more than ", length(colours), " methods and you're using only ", length(colours), " colours.", sep=""))
	}
	
	for (nodes in nodes_list){
		filename <- paste(prefix,"out_nodes-",nodes,"_link-",link,"_drift-",0.1,".txt",sep="")
		data <- read.table(filename, header=T)
		
		dt <- data.table(data)
		dataFilt <- c()
		for (method in methods){
			if (method == "conf"){
				tmp <- dt[, list(pos=sum(pos[pos == nodes])/nodes, neg=sum(neg[neg == nodes])/nodes, nRuns=length(exp[pos == nodes | neg == nodes])), by=exp]
				tmp$metric <- tmp$pos / tmp$nRuns
				dataFilt <- c(dataFilt, tmp[,'metric'])
			} else if (method == "M-rand"){
				tmp <- dt[, list(pos=sum(posMR[posMR == nodes])/nodes, neg=sum(neg[negMR == nodes])/nodes, nRuns=length(exp[posMR == nodes | negMR == nodes])), by=exp]
				tmp$metric <- tmp$pos / tmp$nRuns
				dataFilt <- c(dataFilt, tmp[,'metric'])
			} else if (method == "M-bias"){
				tmp <- dt[, list(pos=sum(posMB[posMB == nodes])/nodes, neg=sum(neg[negMB == nodes])/nodes, nRuns=length(exp[posMB == nodes | negMB == nodes])), by=exp]
				tmp$metric <- tmp$pos / tmp$nRuns
				dataFilt <- c(dataFilt, tmp[,'metric'])
			} else if (method == "M-inhib"){
				tmp <- dt[, list(pos=sum(posMI[posMI == nodes])/nodes, neg=sum(neg[negMI == nodes])/nodes, nRuns=length(exp[posMI == nodes | negMI == nodes])), by=exp]
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
		positions = nodes + (seq(0,length(methods)-1) - floor((length(methods)-1)/2))*size*1.5
		boxplot( dataFilt, at=positions, boxwex=size, add=T, col=colours[1:length(methods)], axes=F)
		
		#aggr <- dt[, list(pos=sum(pos)/(nodes*numRuns), neg=sum(neg)/(nodes*numRuns), best=length(exp[bestDDM == 1])/numRuns), by=exp]
		#boxshift <- 10
		#boxplot( aggr[,'pos'],  at=nodes-boxshift, boxwex=25, add=T)
		#boxplot( aggr[,'best'], at=nodes+boxshift, boxwex=25, col='red', add=T)
	}
	legend('bottomright', methods, fill=rainbow(length(methods)), cex=0.5, bg='white')
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

plotTimeOnNodes <- function(prefix, nodes_list=c(100, 500, 1000), link=0.2, drift=0.1, numRuns=100, methods=c("conf","M-rand","M-bias","M-inhib"), colours=rainbow(10)){
	#if (length(methods)>3) {offset <- 0.08} else {offset<-0}
	offset <- 0.02*length(methods)
	xmin <- min(nodes_list)
	xmax <- max(nodes_list)
	xrange <- (xmax - xmin)
	plot(c(xmin-(xrange*offset),xmax+(xrange*offset)),c(0,5), type='n', xlab="Number of nodes", ylab="Number of iterations")
	grid()
	for (nodes in nodes_list){
		# Read the data
		filename <- paste(prefix,"out_nodes-",nodes,"_link-",link,"_drift-",0.1,".txt",sep="")
		data <- read.table(filename, header=T)
		
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
		positions = nodes + (seq(0,length(methods)-1) - floor((length(methods)-1)/2))*size*1.5
		boxplot( dataFilt, at=positions, boxwex=size, add=T, col=colours[1:length(methods)], axes=F)
		
		#aggr <- dt[, list(time=mean(iter) ), by=exp]
		#boxplot( aggr[,'time'],  at=nodes, boxwex=25, add=T)
	}
	legend('bottomright', methods, fill=colours[1:length(methods)], cex=0.5, bg='white')
}


