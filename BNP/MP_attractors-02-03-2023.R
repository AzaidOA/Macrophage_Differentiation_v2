library(BoolNet)
library(data.table)
library(gtools)
#library(BoolNetPerturb)
library(stringr)
library(dplyr)
for (f in list.files("BNP")) { source( paste('BNP/', f, sep='') ) }

######################### Loading Paths ##################################
net.path <- "model/MP_full_boolean_network.V2_01-20-2023.csv"
net.path <- "model/MP_full_boolean_network.V2_01-30-2023.csv"
#net.path <- "model/MP_full_boolean_network.V2_11-11-22.csv"
net.path <- "model/MP_reduced_network.csv"
#net.path <- "model/MP_full_boolean_network.V2_10-31-22.csv"
#attr.table.path <- "data/MP_attr_v2_ginsim.csv"
lab <- read.csv("model/MP_label_rules-09-04-22.csv")
lab <- read.csv("model/MP_label_rules.csv")

#tag.date <- "16-01-23"
#replace.labels <- list('il6/M0'='il6', 'il6/M0/M2'='M2', 'il6/M2b'='M2b', 'il6/M0/M2b'='M2b', 'M0/M2'='M2', 'M1/M1M2'='M1M2', 'M1/M1M2/NoLabel'='NoLabel', 'M1/M1M2/NoLabel'='NoLabel', 'M1M2/NoLabel'='NoLabel', 'il6/M2'='M2', 'M2/M2b'='M2b', 'M2/M2b'='M2b', 'il6/M2/M2b'='M2b', 'M0/M2/M2b'='M2b', 'M2/M2d'='M2d', 'M0/M2a'='M2a', 'M0/M2b'='M2b', 'M0/M2/M2c'='M2c')
######################### Getting Attractors ##################################
net <- loadNetwork(file = net.path)
attr <- getAttractors(network = net, method = "sat.restricted", maxAttractorLength = 1)
#attr.table <- attractorToDataframe(attr, Boolean=TRUE)
#attr.table <- fread(input = attr.table.path)
#attr.table <- apply(attr.table, 1, function(attractor){
#  attractor <- as.data.table(t(attractor))
#  asterik_deconvolution(attractor = attractor)
#})
#attr.table <- rbindlist(attr.table)
#attr.table[,attractor:=1:.N]
####################
#AO: Only works for attractors of fixed point
labelAttractors <- function(attr, label.rules, node.names=NULL, sep="/") {
  # takes an attractors object created by BoolNet
  # returns a list of the labels for each attractor in order.
  # If an attractor has multiple states it will return a label for each state.
  if (!is(attr, "AttractorInfo")) { stop("Error: non-valid attractor") }
  node.names <- attr$stateInfo$genes
  res <- list()

  #for (i in 1:length(attr$attractors)) {
  #    label <- sapply(attr$attractors[[i]]$involvedStates, function(state) {
  #    state <- int2binState(state, node.names) #state to binary
  #    l <- labelState(state, node.names, label.rules, sep='') #label
  #  })
  #  if (!is.null(sep)) { label <- paste(label, collapse=sep) }
  #  res <- append(res, list(label))
  #}

  for(i in 1:length(attr$attractors)){
    state <- attr$attractors[[i]]$involvedStates
    state <- int2binState(state, node.names) #state to binary
    label <- labelState(state, node.names, label.rules, sep='') #label
    if (!is.null(sep)) { label <- paste(label, collapse=sep) }
    res <- append(res, list(label))
  }
  unlist(res)
}

####################### Getting labels #########################################
labels.attr <- labelAttractors(attr, lab, net$genes)
attr.table <- as.data.table(attractorToDataframe(attr=attr,Boolean=T))
attr.table[,label:=labels.attr[attractor]]

attr.table[IC_e == 1 & IL1B_e == 1 & LPS_e == 1 & EGFR_e == 1 & IFNG_e == 0 & GMCSF_e == 0 & NECA_e == 0 & IL4_e == 0 & IL6_e == 0 & IL10_e == 0]

# Analysis on first network about the estimulus of m2b macrophages.
est.cols <- c('attractor',net$genes[net$genes %like% '_e'])
m2b.estimulus <- unique(attr.table[label == 'M2b',..est.cols])
m2b.estimulus[,total.genes := rowSums(m2b.estimulus[,!'attractor'])]
m2b.estimulus[total.genes == 4]

# As we can see here EGFR_e is always on!.

attr.table[attractor == 1871,]
# IL4 stimulus
tmp.state <- rep(0,length(net$genes)); names(tmp.state) <- net$genes
tmp.state['IL4_e'] <- 1

getPathToAttractor(network=net,state = tmp.state)

last.state <- getPathToAttractor(network=net,state = tmp.state)
last.state <- last.state[nrow(last.state),]
last.state <- getPathToAttractor(network=net,state = last.state)
last.state <- last.state[nrow(last.state),]
getPathToAttractor(network=net,state = last.state)

getPathToAttractor(network=net,state = tmp.state)[,c('IL4_e','IL10_out','STAT6','SOCS1','IL12_out','IL6_out','NFKB','VEGF_out','FCGR')]

# LPS and IgG stimulus
tmp.state <- rep(0,length(net$genes)); names(tmp.state) <- net$genes
tmp.state[c('LPS_e','IgG_e','IL1B_e','EGFR_e')] <- 1
m2b.state <- getPathToAttractor(network=net,state = tmp.state)
getPathToAttractor(network=net,state = m2b.state[nrow(m2b.state),])

m2b.list <- list(tmp.state)
m2b.attr <- getAttractors(network = net, method = 'chosen', startStates=m2b.list)
m2b.attr.table <- as.data.table(attractorToDataframe(attr=m2b.attr,Boolean=T))

getPathToAttractor(network=net,state = tmp.state)[,c('IL12_out','TLR4','TRAF2','STAT1','STAT5','NFKB','IL10_out','IL6_out','IL10R','STAT3_IL10','FCGR','PPARG','STAT6','VEGF_out')]

# I think it is NFKB what is activating IL12 in a M2b enviroment.
m2b.test.net <- fixGenes(network=net,fixIndices='NFKB',values=0)
getPathToAttractor(network=m2b.test.net,state = tmp.state)[,c('IL12_out','STAT1','STAT5','NFKB','IL10_out','IL10R','STAT3_IL10','FCGR','PPARG')]

# I think that the facor upstream is FCGR who is activating NFKB
m2b.test.net <- fixGenes(network=net,fixIndices='FCGR',values=0)
getPathToAttractor(network=m2b.test.net,state = tmp.state)[,c('FCGR','IL10_out','IL6_out','STAT3_IL6','IL6R','IL12_out','STAT6','mir155','NFKB','VEGF_out')]
# NOTE: It wasn't, IL12 is still expressing.

# I have changed the rule of the node il12_out as the one used in the first network. This works but this rule it's not true at all.

# I have read that mir146 can increased PPARG and inhibit NFKB. The NFKB inihibition works,
# so I will set an OE in PPARG to see if this could work too.
tmp.state <- rep(0,length(net$genes)); names(tmp.state) <- net$genes
tmp.state[c('LPS_e','IgG_e','IL1B_e','EGFR_e','PPARG')] <- 1
m2b.test.net <- fixGenes(network=net,fixIndices='PPARG',values=1)
getPathToAttractor(network=m2b.test.net,state = tmp.state)[,c('IL12_out','STAT1','STAT5','NFKB','IL10_out','IL10R','STAT3_IL10','FCGR','PPARG')]
## It does!!

# Now by TLR4 inhibition.
tmp.state <- rep(0,length(net$genes)); names(tmp.state) <- net$genes
tmp.state[c('LPS_e','IgG_e','IL1B_e','EGFR_e')] <- 1
m2b.test.net <- fixGenes(network=net,fixIndices=c('TRAF2','FCGR'),values=c(0,0))
getPathToAttractor(network=m2b.test.net,state = tmp.state)[,c('IL12_out','TLR4','STAT1','STAT5','NFKB','IL10_out','IL10R','STAT3_IL10','FCGR','PPARG')]

# It's the rule, but based in the searched carried out we have concluded before this

estimulus <- colnames(attr.table)[colnames(attr.table) %like% '_e']
estimulus.m2b <- unique(attr.table[label=='M2b',..estimulus])
#################### PLotting attractors ####################################
images.dir <- paste0("images/",tag.date)
if(!dir.exists(images.dir)) dir.create(images.dir)
file.attr.pdf <- paste0(images.dir,"/MP_attr_all")


# Select key nodes for plotting
node.subset <- c("STAT1","mir155","SOCS1","IL4_e","STAT5", "NFKB", "STAT6", "STAT3_IL10", "STAT3_IL6", "STAT4", "IL12_out", "IL10_out", "IL6_out", "VEGF_out","IFNa_out","IFNG_out","TNFa_out","IFNb_out")
attr.table.short <- as.data.frame(subset(attr.table, select = c("attractor","label","state", node.subset) ))

# select unique attractors of subset
attrs <- list()
for (n in 1:max(attr.table$attractor)) {
  at <- attr.table.short[attr.table.short$attractor==n,]
  at <- subset(at, select = -c(attractor,state) )
  at[nrow(at) + 1,] = NA
  row.names(at) <- NULL
  attrs[[n]] <- at
}
attrs <- unique(attrs)

for (n in unique(sapply(attrs, nrow))) {
  n <-  n-1
  print(n)
  attrs_f <- Filter(function(x) (nrow(x))-1==n, attrs)
  attrs_f <- bind_rows(attrs_f)
  label_f <- attrs_f$label
  print(table(label_f))
  title = paste(c('Attractors of size ',n), collapse='')
  attrs_f <- t(as.matrix(sapply(  subset(attrs_f, select=-c(label))  , as.numeric)))
  colnames(attrs_f) <- label_f
  if (n==1) {
    attrs_f <- attrs_f[,! is.na(colnames(attrs_f))]
    attrs_f <- attrs_f[ , order(colnames(attrs_f))]
  }
  pdf(paste0(file.attr.pdf,"_",n,".pdf"))
  heatmap(attrs_f, main=title,
          col=c('#fb8072','#b3de69'), cexCol=0.75, cexRow=1,
          Colv = NA, Rowv = NA, scale="none" )
  dev.off()
}
