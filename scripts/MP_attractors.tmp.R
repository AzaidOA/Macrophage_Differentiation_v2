### -------------------------- Loading Libraries -------------------------- ###
library(BoolNet)
library(data.table)
library(gtools)
#library(BoolNetPerturb)
library(ComplexHeatmap)
library(stringr)
library(dplyr)
library(R.utils)
library(alluvial)
library(parallel)
# for (f in list.files("BNP")) { source( paste('BNP/', f, sep='') ) }
source('BNP/BNP_Dataframe.R')
source('BNP/BNP_Helper.R')
source('BNP/BNP_Label-03-12-2023.R')
source('BNP/BNP_PerDynamic.R')
source('BNP/BNP_PerFunctions.R')

### ------------------------------ Functions ------------------------------ ###
simplifyLabel <- function(old, sep='/', mark='*', replace.labels=NULL) {
  # separate abd simplify
  new <- str_split(old,sep)
  new <- sort(unique(new[[1]]))
  new <- paste(new, collapse=sep)
  # replace special cases
  if (! is.null(replace.labels)) {
    # order by descending size to avoid problems
    rl <- names(replace.labels)
    rl <- rl[order(nchar(rl), rl, decreasing=T)]
    replace.labels <- replace.labels[rl]
    # replace all
    for (key in names(replace.labels)) {
      new <- str_replace( new, pattern=key, replace.labels[[key]] )
    }
  }
  # mark the ones that have changed
  if (! is.null(mark)) {
    if (new!=old) {
      new <- paste(c(new,mark), collapse='')
    }
  }
  new
}

### -------------------------- Arguments -------------------------- ###
# ---> File of the network rules to use.
net.path <- "model/MP_full_boolean_network.V2_02-28-2023.csv"
# ---> File of the labeling rules to use.
lab.path <- "model/MP_label_rules-09-04-22.csv"
# ---> Maximun size of the attractors.
attr.len <- 3
# ---> Environments
env <- read.csv("model/MP_environment.csv")

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
fwrite(attr.table,file=paste0('data/MP_attr.csv'))
# ---> Getting a summary of the obteined attractors.
# "Label" is the tagged dads by the labeling rules.
# "States" is the classification of the length of the attractors.
# "Attractors" are the number of attractors of that length.
cat('Summary of the attractors obtained:\n')
print(attr.summary)
tmp.output.file <- paste0('data/MP_attr_summary.csv')
fwrite(x=attr.summary, file=tmp.output.file)
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


### ---> Simplifiying attractors
replace.labels <- c(
  'M2c/M2c'='M2c','M2d/M2d'='M2d','M2/M2b'='M2b','M0/M0'='M0')
attr.replace <- attr.table[label %chin% names(replace.labels)]
attr.replace[,label:=paste0(replace.labels[label],'*')]


### ---> Plotting Attractors
file.attr.pdf <- "images/MP_attr_all"
# Select key nodes for plotting
node.subset <- c('IL12_out','IL10_out','IL6_out','VEGF_out','STAT1','STAT5','NFKB','STAT6','FCGR','STAT3_IL10','STAT3_IL6','SOCS1','SOCS3')
attr.table.short <- subset(attr.table, select = c("attractor","label","state", node.subset) )

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
  # n <-  n-1
  print(n)
  attrs_f <- Filter(function(x) (nrow(x))==n, attrs)
  # lapply(attrs_f,function(x) {
  #   x[ , (order(colnames(x[,-'label']))+1)]
  # }
  
  attrs_f <- bind_rows(attrs_f, .id='attractor')
  label_f <- attrs_f$label
  print(table(label_f))
  title <-  paste(c('Attractors of size ',n), collapse='')
  attrs_f <- t(as.matrix(sapply(  subset(attrs_f, select=-c(label))  , as.numeric)))
  colnames(attrs_f) <- label_f
  # if (n==1) {
  attrs_f <- attrs_f[,! is.na(colnames(attrs_f))]
  attrs_f <- attrs_f[ , order(colnames(attrs_f))]
  # }
  # attrs_f <- t(attrs_f)
  tmp.split <- attrs_f[1,]
  tmp.split <- match(tmp.split, unique(tmp.split))
  attrs_f <- attrs_f[-1,]
  # attrs_f <- t(attrs_f)
  tmp.output.file <- paste0(file.attr.pdf,'_',n,'.pdf')
  if(n == 3){
    n.pages <- ceiling(ncol(attrs_f)/24)
    pdf(tmp.output.file)
    start.ind <- 1; end.ind <- 24
    for (tmp.p in 1:n.pages){
      h1 <- Heatmap(attrs_f[,start.ind:end.ind], name=title, cluster_rows = FALSE,
                    col=c('#fb8072','#b3de69'), cluster_columns = FALSE,
                    column_split=tmp.split[start.ind:end.ind], row_dend_reorder = FALSE,
                    column_dend_reorder = FALSE,
                    show_heatmap_legend = FALSE, show_column_dend=FALSE,
                    show_parent_dend_line=FALSE,
                    column_title = NULL
      )
      print(h1)
      start.ind <- start.ind + 24; end.ind <- end.ind + 24;
      if(end.ind > ncol(attrs_f)) end.ind <- ncol(attrs_f)
    }
    dev.off()
    next
  }
  pdf(tmp.output.file)
  h1 <- Heatmap(attrs_f, name=title, cluster_rows = FALSE,
                col=c('#fb8072','#b3de69'), cluster_columns = FALSE,
                column_split=tmp.split, row_dend_reorder = FALSE,
                column_dend_reorder = FALSE,
                show_heatmap_legend = FALSE, show_column_dend=FALSE,
                show_parent_dend_line=FALSE,
                column_title = NULL
  )
  print(h1)
  dev.off()
}
invert = TRUE




# ---> Asynchronous
initial.attr <- apply(
  X=attr.table[,!c('attractor','state','label')]
  ,MARGIN=1,
  FUN=function(x){
    as.data.frame(t(x))
  })
# attr.asyn <- getAttractors(
#   network = net, method = "chosen",
#   type='asynchronous',startStates=initial.attr)
# attr.table.asyn <- as.data.table(attractorToDataframe(attr=attr.asyn,Boolean=T))
# # Applying labeling rules.
# labels.attr.asyn <- labelAttractors(attr.asyn, lab, net$genes, sep="/")
# attr.table.asyn[,label:=labels.attr.asyn[attractor]]
#
# attr.summary.asyn <- attr.table.asyn[,.(states=.N),by=.(label,attractor)][
#   ,.(attractors=.N),by=.(label,states)]
# attr.summary.asyn

### ------------------------ Asynchronous ------------------------ ###
# ---> Asynchronous
set.seed(534234)
attr.asyn <- lapply(X=1:length(initial.attr), FUN=function(tmp.x){
  tmp.attr <- initial.attr[tmp.x]
  print(tmp.x)
  tmp.attr.asyn <- tryCatch(expr= {
    withTimeout({getAttractors(
      network = net, method = "chosen",
      type='asynchronous',startStates=tmp.attr)}, timeout=5)
  },
  error=function(i) NA,
  TimeoutException = function(ex) NA)
  if(all(!is.na(tmp.attr.asyn))){
    tmp.attr.table.asyn <- as.data.table(attractorToDataframe(attr=tmp.attr.asyn,Boolean=T))
    tmp.labels.attr.asyn <- labelAttractors(tmp.attr.asyn, lab, net$genes, sep="/")
    if(is.null(dim(tmp.labels.attr.asyn))){
      tmp.attr.table.asyn[,label:=tmp.labels.attr.asyn]
    } else {
      tmp.attr.table.asyn[,label:=tmp.labels.attr.asyn[attractor]]
    }
    
  } else {
    tmp.attr.table.asyn <- as.data.table(matrix(ncol=length(net$genes) + 3, nrow=1))
    colnames(tmp.attr.table.asyn) <- c('attractor','state',net$genes,'label')
  }
  tmp.attr.table.asyn[,attractor:=tmp.x]
  return(tmp.attr.table.asyn)
})
attr.table.asyn <- rbindlist(attr.asyn)

# Saving asyn attractators.
tmp.output.file <- paste0('data/MP_asyn_attr.csv')
fwrite(attr.table.asyn,file=tmp.output.file)

tmp.output.file <- paste0('data/MP_asyn_attr_summary.csv')
attr.summary.asyn <- attr.table.asyn[,.(states=.N),by=.(label,attractor)][
  ,.(attractors=.N),by=.(label,states)]
fwrite(x=attr.summary.asyn, file=tmp.output.file)
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


# attr <- getAttractors(net, method="sat.exhaustive")
async <- verifySyncronousVsAsyncronous(net, attr, lab)
async

### --------------------------- Mutants analysis --------------------------- ###
file.mut <- 'data/MP_mut_label.csv'
mutants <- perturbNetworkFixedNodes(net, label.rules=lab, 
                                    returnDataFrame='occurrence', method="sat.restricted",
                                    maxAttractorLength = attr.len)

# Making mutant analysis over an iteration of each gene.
# NOTE: Error generated by canonical why is this happening? Ask if deactivating is ok
mutants <- mclapply(X=net$genes, mc.cores=100,
  FUN=function(tmp.gene){
    tmp.mutants <- perturbNetworkFixedNodes(net, label.rules=lab, verbose=T,
                                        genes=rep(tmp.gene,2), values=c(0,1),
                                        returnDataFrame='occurrence', method="sat.restricted",
                                        maxAttractorLength = attr.len, canonical=F)
    return(tmp.mutants)
    })

mutants <- Reduce(function(x, y){
  tmp.merged <- merge(x, y, by = 0, all = TRUE)
  rownames(tmp.merged) <- tmp.merged$Row.names
  tmp.merged$Row.names <- NULL
  tmp.merged
  }, 
  mutants)
mutants[is.na(mutants)] <- 0
net.socs1 <- fixGenes(network = net, fixIndices = 'SOCS1', values = 0)
attr.socs1 <- getAttractors(network = net.socs1, method = "sat.restricted", maxAttractorLength = attr.len)
mutants.socs1 <- attractorListToDataframe(mutants, sep=sep, returnDataFrame=returnDataFrame)

perturbNetworkFixedNodes(net, label.rules=lab, verbose=T,
                         genes='IFNa_out', values=0,
                         returnDataFrame='occurrence', method="sat.restricted",
                         maxAttractorLength = 3)
tmp.attr <- getAttractors(net,method="sat.restricted",maxAttractorLength = 3)
Error in if (stateMatrix[j, i] < smallestVal[j]) { : missing value where TRUE/FALSE needed

### --------------------- Environment analysis --------------------- ###
env <- read.csv("model/MP_environment.csv")

file.env.attr <- "data/MP_env_attr_basin.csv"

if (! file.exists(file.env.attr)) {
  env.attr <- perturbNetworkFixedNodes(net, label.rules=lab,
                                       returnDataFrame='basinSize',
                                       genes  = rep( list(colnames(env)), times=nrow(env) ),
                                       values = lapply( split(env,seq_along(env[,1])), as.list),
                                       names  = rownames(env))
  setDT(env.attr, keep.rownames = TRUE)
  colnames(env.attr)[[1]] <- 'label'
  env.attr
  label, simplifyLabel, replace=replace.labels)
env.attr <- env.attr %>% group_by(label) %>% summarize_all(sum)
write.csv(env.attr, file.env.attr, row.names=F)
} else {
  env.attr <- read.csv(file = file.env.attr, row.names=1, stringsAsFactors=F)
  env.attr[env.attr == 0] <- NA
}

file.env.pdf <- "images/MP_env_attr.pdf"

normalize = F
if (normalize) {
  env.attr <- env.attr/env.attr
  color = '#a65628'
} else {
  colfunc <- colorRampPalette(c('#fee391', '#a65628'))
  color <- colfunc(10)
}

#pdf(file.env.pdf)
heatmap(t(as.matrix( subset(env.attr ))),
        main="Macrophage by environment",
        xlab="label", ylab="Environment",
        col=color, cexCol=0.75, cexRow=0.75,
        Colv = NA, Rowv = NA, scale="none",
)
#dev.off()
