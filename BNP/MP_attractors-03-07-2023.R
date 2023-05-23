### -------------------------- Loading Libraries -------------------------- ###
library(BoolNet)
library(data.table)
library(gtools)
#library(BoolNetPerturb)
library(stringr)
library(dplyr)
# for (f in list.files("BNP")) { source( paste('BNP/', f, sep='') ) }
source('BNP/BNP_Dataframe-10-08-22.R')
source('BNP/BNP_Helper-10-08-22.R')
source('BNP/BNP_Label-10-08-22.R')
### -------------------------- Arguments -------------------------- ###
# ---> File of the network rules to use.
net.path <- "model/MP_full_boolean_network.V2_02-22-2023.csv"
net.path <- "model/MP_full_boolean_network.V2_02-26-2023.csv"
net.path <- "model/MP_full_boolean_network.V2_02-28-2023.csv"
# ---> File of the labeling rules to use.
lab.path <- "model/MP_label_rules-09-04-22.csv"
# ---> Maximun size of the attractors.
attr.len <- 1

### -------------------------- Getting Attractors -------------------------- ###
# ---> Loading the network.
net <- loadNetwork(file = net.path)

attr <- getAttractors(network = net, method = "sat.restricted", maxAttractorLength = attr.len)
attr.table <- as.data.table(attractorToDataframe(attr=attr,Boolean=T))

### -------------------------- Labeling The Attractors -------------------------- ###
# ---> Loading the labels
lab <- read.csv(lab.path)
# Applying labeling rules.
labels.attr <- labelAttractors(attr, lab, net$genes)
attr.table[,label:=labels.attr[attractor]]
# ---> Getting a summary of the attractor object.
attr.summmary <- attr.table[,.(states=.N),by=.(label,attractor)][
  ,.(attractors=.N),by=.(label,states)]

# ---> Getting a summary of the obteined attractors.
# "Label" is the tagged dads by the labeling rules.
# "States" is the classification of the length of the attractors.
# "Attractors" are the number of attractors of that length.
cat('Summary of the attractors obtained:\n')
print(attr.summmary)
tmp.attr.summ <- attr.summmary
### -------------------------- Environment analysis -------------------------- ###
# ---> Adding an extra column to identify attractors faster.
collapsed.attr <- apply(X=attr.table[,!c('state','attractor','label')],
  MARGIN=1,FUN=paste,collapse='')
attr.table[,collapsed:=as.character(collapsed.attr)]

# ---> Making analysis for M1 macrophages: IFNG, GMCSF and LPS stimulus.
cat('Making analysis for M1 macrophages: IFNG, GMCSF and LPS stimulus.\n')
# Creating initital state where only the stimulus is activated.
initial.state <- rep(0,length(net$genes)); names(initial.state) <- net$genes
initial.state[c('IFNG_e','GMCSF_e','LPS_e')] <- 1
# Nodes of interest for this enviroment.
nodes.of.int <- c('IL12_out','STAT3_IL10','IL10_out','IL6_out','SOCS1','STAT1','STAT5','NFKB','CXCL10_out','VEGF_out')

path.to.attr <- getPathToAttractor(network=net,state = initial.state)
last.state <- tail(path.to.attr,1)
cat('First iteration.\n')
path.to.attr <- rbind(path.to.attr[-nrow(path.to.attr),],
  getPathToAttractor(network=net,state = last.state))
last.state <- tail(path.to.attr,1)
cat('Second iteration:\n')
# NOTE: Getting an error in this iteration.
path.to.attr <- rbind(path.to.attr[-nrow(path.to.attr),],
  getPathToAttractor(network=net,state = last.state))

last.state <- as.character(apply(X=last.state,MARGIN=1,FUN=paste,collapse=''))
tmp.check <- attr.table[collapsed == last.state,.N == 0]
if(tmp.check)
  cat('There is no coincidence in the attractors table. This is a cyclic attractor.\n')

cat('Path followed by the attractor using the nodes of interest:\n')
print(path.to.attr[,nodes.of.int])
cat('Observation:
As we can see the IL12 expression has beeen totally
.\n')

# ---> Making analysis for M2a macrophages: IL4 stimulus.
cat('Makling analysis for M2a macrophages: IL4 stimulus.\n')
# Creating initital state where only the stimulus is activated.
initial.state <- rep(0,length(net$genes)); names(initial.state) <- net$genes
initial.state['IL4_e'] <- 1
# Nodes of interest for this enviroment.
nodes.of.int <- c('IL4_e','IL10_out','SOCS1','STAT6','STAT1','IL12_out','IL6_out',
  'NFKB','VEGF_out','FCGR')


path.to.attr <- getPathToAttractor(network=net,state = initial.state)
last.state <- tail(path.to.attr,1)
cat('First iteration.\n')
path.to.attr <- rbind(path.to.attr[-nrow(path.to.attr),],
  getPathToAttractor(network=net,state = last.state))
last.state <- tail(path.to.attr,1)
cat('Second iteration:\n')
# NOTE: Getting an error in this iteration.
path.to.attr <- rbind(path.to.attr[-nrow(path.to.attr),],
  getPathToAttractor(network=net,state = last.state))

last.state <- as.character(apply(X=last.state,MARGIN=1,FUN=paste,collapse=''))
tmp.check <- attr.table[collapsed == last.state,.N == 0]
if(tmp.check)
  cat('There is no coincidence in the attractors table. This is a cyclic attractor.\n')

cat('Path followed by the attractor using the nodes of interest:\n')
print(path.to.attr[,nodes.of.int])
cat('Observation:
As we can see in the pathway followed by the attractor, SOCS1 and STAT6 are interplaying, causing these to be inhibited and activated cyclically. But as we can see the rule for m2b macrophages is being fulfilled, just one of its markers, STAT6, is being activated and inhibited.\n')

# ---> Making analysis for M2b macrophages: LPS, IgG, IL1B and EGFR estimulus.
cat('Making analysis for M2b macrophages: LPS, IgG, IL1B and EGFR estimulus.\n')
# Creating initital state where only the stimulus is activated.
initial.state <- rep(0,length(net$genes)); names(initial.state) <- net$genes
initial.state[c('LPS_e','IgG_e','IL1B_e','EGFR_e')] <- 1
# Nodes of interest for this enviroment.
nodes.of.int <- c('IL10_out','IL6_out','STAT6','VEGF_out','FCGR','IL12_out',
  'STAT4','NFKB','IFNGR','IFNb_out','IFNG_out','TNFa_out')

path.to.attr <- getPathToAttractor(network=net,state = initial.state)
last.state <- tail(path.to.attr,1)
cat('First iteration.\n')
path.to.attr <- rbind(path.to.attr[-nrow(path.to.attr),],
  getPathToAttractor(network=net,state = last.state))
last.state <- tail(path.to.attr,1)
cat('Second iteration:\n')
# NOTE: Getting an error in this iteration.
path.to.attr <- rbind(path.to.attr[-nrow(path.to.attr),],
  getPathToAttractor(network=net,state = last.state))

last.state <- as.character(apply(X=last.state,MARGIN=1,FUN=paste,collapse=''))
tmp.check <- attr.table[collapsed == last.state,.N == 0]
if(tmp.check)
  cat('There is no coincidence in the attractors table. This is a cyclic attractor.\n')

cat('Path followed by the attractor using the nodes of interest:\n')
print(path.to.attr[,nodes.of.int])
cat('Observation:
As we can see in the attractor pathway, IL12 is being secreted due to the activation of NFKB. This is causing STAT4 to activate and therefore the TNF pathway, the IFN pathway and the IL12 autocrine loop. This completely inhibits M2 polarization.\n')
