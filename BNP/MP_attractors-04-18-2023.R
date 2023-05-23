### -------------------------- Loading Libraries -------------------------- ###
library(BoolNet)
library(data.table)
library(gtools)
#library(BoolNetPerturb)
library(stringr)
library(dplyr)
library(R.utils)
library(alluvial)
# for (f in list.files("BNP")) { source( paste('BNP/', f, sep='') ) }
source('BNP/BNP_Dataframe-10-08-22.R')
source('BNP/BNP_Helper-03-12-2023.R')
source('BNP/BNP_Label-03-12-2023.R')
source('BNP/BNP_PerDynamic.R')
### -------------------------- Arguments -------------------------- ###
# ---> File of the network rules to use.
net.path <- "model/MP_full_boolean_network.V2_02-28-2023.csv"
# ---> File of the labeling rules to use.
lab.path <- "model/MP_label_rules-09-04-22.csv"
# ---> Maximun size of the attractors.
attr.len <- 3

### -------------------------- Getting Attractors -------------------------- ###
# ---> Loading the network.
net <- loadNetwork(file = net.path)
single.attr <- attr
attr <- getAttractors(network = net, method = "sat.restricted", maxAttractorLength = attr.len)
# attr <- getAttractors(network = net, method = "sat.exhaustive")
attr.table <- as.data.table(attractorToDataframe(attr=attr,Boolean=T))
# tmp.attr.table <- attr.table[1:.N]
### -------------------------- Labeling The Attractors -------------------------- ###
# ---> Loading the labels
lab <- read.csv(lab.path)
# Applying labeling rules.
labels.attr <- labelAttractors(attr, lab, net$genes, sep="/")
attr.table[,label:=labels.attr[attractor]]
# ---> Getting a summary of the attractor object.
attr.summary <- attr.table[,.(states=.N),by=.(label,attractor)][
  ,.(attractors=.N),by=.(label,states)]

# ---> Getting a summary of the obteined attractors.
# "Label" is the tagged dads by the labeling rules.
# "States" is the classification of the length of the attractors.
# "Attractors" are the number of attractors of that length.
cat('Summary of the attractors obtained:\n')
print(attr.summary)
# tmp.attr.summ <- attr.summary


feats.of.int <- c('IL12_out','IL10_out','IL6_out','VEGF_out','STAT1','STAT5','NFKB','STAT6','FCGR','STAT3_IL10','STAT3_IL6','SOCS1','SOCS3')
# ---> Assesing the cyclic attractors.

tmp.feats.of.int <- c(feats.of.int,'IFNG_e','IFNG_out','STAT4')
attr.of.int <- 'M2d/M2d/il6'
tmp.attr.table <- attr.table[label == attr.of.int, ..tmp.feats.of.int]
unique(tmp.attr.table)
# Why is STAT 1 activating?

tmp.feats.of.int <- c(feats.of.int)
attr.of.int <- 'M2c/M2/M0'
tmp.attr.table <- attr.table[label == attr.of.int, ..tmp.feats.of.int]
unique(tmp.attr.table)
# SOCS3 is triggering that IL10 being cycled, but even when IL10 is off STAT3_10 another maker of m2c is on, and any other conflicted marker is on.

# NOTE: M2c

tmp.feats.of.int <- c(feats.of.int,'IL10R')
attr.of.int <- 'M2/M0/M0'
tmp.attr.table <- attr.table[label == attr.of.int, ..tmp.feats.of.int]
unique(tmp.attr.table)

# Why is activated 1IL10? There is no stimulus
# Posible noising attractor.

tmp.feats.of.int <- c(feats.of.int)
attr.of.int <- 'M0/M2/M0'
tmp.attr.table <- attr.table[label == attr.of.int, ..tmp.feats.of.int]
unique(tmp.attr.table)
attr.table[label == attr.of.int][1:3]
# Why is activated 1IL10? There is no stimulus
# Posible noising attractor.

tmp.feats.of.int <- c(feats.of.int)
attr.of.int <- 'M0/M2/M1'
tmp.attr.table <- attr.table[label == attr.of.int, ..tmp.feats.of.int]
unique(tmp.attr.table)
attr.table[label == attr.of.int][1:3]
# There are no stimulus for IL12 and IL10

tmp.feats.of.int <- c(feats.of.int)
attr.of.int <- 'M2d/M2c/il6'
tmp.attr.table <- attr.table[label == attr.of.int, ..tmp.feats.of.int]
unique(tmp.attr.table)
attr.table[label == attr.of.int][1:3]
# This should be M2d or TAMs, due what is happening SOCS1 is regulating the VEGF production yn a clycler manner, but always the attractor is secreting or IL10 or IL6 what are cytoines produced by M2d.
# NOTE: M2d

tmp.feats.of.int <- c(feats.of.int)
attr.of.int <- 'M0/M0/M1'
tmp.attr.table <- attr.table[label == attr.of.int, ..tmp.feats.of.int]
unique(tmp.attr.table)
attr.table[label == attr.of.int][1:3]
# This should be M1, SOCS1 is inhibiting ina cycled manner the STAT1 or STAT5 that TFs for IL12_out

# NOTE: M1

tmp.feats.of.int <- c(feats.of.int)
attr.of.int <- 'M2d/M2d/M1'
tmp.attr.table <- attr.table[label == attr.of.int, ..tmp.feats.of.int]
unique(tmp.attr.table)
attr.table[label == attr.of.int][1:3]
attr.of.int <- c('M2/M0/M0','M0/M2/M0','M0/M2/M1','M2d/M2/M1','NoLabel/M2/M1','M2/M2/M1','M2d/M2d/M1','M2d/M0/M1')
estimulus <- c(net$genes[net$genes %like% '_e'],'label')
estimulus.attr <- attr.table[label %chin% attr.of.int, ..estimulus]
estimulus.attr <- unique(estimulus.attr)
non.estimulus <- setdiff(net$genes,estimulus)
for (x in non.estimulus) {
  estimulus.attr <-cbind(estimulus.attr,tmp.tag=0)
  colnames(estimulus.attr)[ncol(estimulus.attr)] <- x
}
label.tags <- estimulus.attr$label
estimulus.attr$label <- NULL

estimulus.attr <- apply(X=estimulus.attr,MARGIN=1,FUN=function(x){
  as.data.table(t(x))
})
names(estimulus.attr) <- label.tags

initial.attr <- apply(
  X=attr.table[,!c('attractor','state','label')]
  ,MARGIN=1,
  FUN=function(x){
    as.data.frame(t(x))
})

# ---> Asynchronous
attr.asyn <- getAttractors(
  network = net, method = "chosen",
  type='asynchronous',startStates=initial.attr)
attr.table.asyn <- as.data.table(attractorToDataframe(attr=attr.asyn,Boolean=T))
# Applying labeling rules.
labels.attr.asyn <- labelAttractors(attr.asyn, lab, net$genes, sep="/")
attr.table.asyn[,label:=labels.attr.asyn[attractor]]

attr.summary.asyn <- attr.table.asyn[,.(states=.N),by=.(label,attractor)][
  ,.(attractors=.N),by=.(label,states)]
attr.summary.asyn

# ---> Asynchronous
set.seed(534234)
attr.asyn <- lapply(X=1:length(initial.attr), FUN=function(tmp.x){
  tmp.attr <- initial.attr[tmp.x]
  print(tmp.x)
  tmp.attr.asyn <- tryCatch(expr= {
    withTimeout({getAttractors(
      network = net, method = "chosen",
      type='asynchronous',startStates=tmp.attr)}, timeout=3)
  },
    error=function(i) NA,
    TimeoutException = function(ex) NA)
  if(!is.na(tmp.attr.asyn)){
    tmp.attr.table.asyn <- as.data.table(attractorToDataframe(attr=tmp.attr.asyn,Boolean=T))
    tmp.labels.attr.asyn <- labelAttractors(tmp.attr.asyn, lab, net$genes, sep="/")
    if(is.null(dim(tmp.labels.attr.asyn))){
        tmp.attr.table.asyn[,label:=tmp.labels.attr.asyn]
    } else {
        tmp.attr.table.asyn[,label:=tmp.labels.attr.asyn[attractor]]
    }

  } else{
    tmp.attr.table.asyn <- as.data.table(matrix(ncol=length(net$genes) + 3, nrow=1))2
    colnames(tmp.attr.table.asyn) <- c('attractor','state',net$genes,'label')
  }
  tmp.attr.table.asyn[,attractor:=tmp.x]
  return(tmp.attr.table.asyn)
})
attr.table.asyn <- rbindlist(attr.asyn)
attr.summary.asyn <- attr.table.asyn[,.(states=.N),by=.(label,attractor)][
  ,.(attractors=.N),by=.(label,states)]
attr.summary.asyn

tmp.feats.of.int <- c('state',feats.of.int,'HIF1a','STAT3_IL6')
attr.of.int <- 'M2b/M2b/M2b/M2b/M2b/M2b/M2b/M2b/M2b/M2b/M2b/M2b/M2/M2/M2/M2/M2/M2/M2/M2/M2/M2/M2/M2'
tmp.attr.table.asyn <- attr.table.asyn[label == attr.of.int, ..tmp.feats.of.int]
unique(tmp.attr.table.asyn)

cat('VEGF is being activated due to the IL6 out, trigerring the phenotybe being classified as M2 although it has all the M2b markers, this being cycled due SOCS3 that can inhibit IL6R->STAT3_IL6->VEGF')

# Transition between synchronus and asynchronous
tmp.file.name <- 'images/MP2_transition_syn_vs_asyn.pdf'
transition.table <- merge(x=attr.table[,.(attractor=1:.N,syn.lab=label)],
  y=unique(attr.table.asyn[,.(attractor,asyn.lab=label)]),
  by='attractor')

transition.table <- transition.table[,.(freq=.N),
  by=.(syn.lab,asyn.lab)]

transition.table[is.na(transition.table)] <- 'NoAttractor'

pdf(tmp.file.name, height=20)
alluvial(transition.table[,.(syn.lab,asyn.lab)],
freq=transition.table$freq)
dev.off()
