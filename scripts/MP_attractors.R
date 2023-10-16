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
lab.path <- "model/MP_label_rules.tmp.csv"
# ---> Maximun size of the attractors.
attr.len <- 3
# ---> Environments
env.file <- "model/MP_environment.csv"
# ---> Simplied labels.
replace.labels <- c('M2/M2b'='M2b','M0/M2/M2c'='M2c','il6/M2d'='M2d-like','M1/M2d'='M2d-like','M1/M2/M2d'='M2d-like','M0/M1/M2d'='M2d-like','M0/M1/M2'='Discarded','M0/M2'='M2c','M0/M1'='M1','il6/M2c/M2d'='M2d-like','M1/M2'='Discarded','M1/M2/NoLabel'='Discarded')
### -------------------------- Loading data -------------------------- ###
# ---> Loading the network.
net <- loadNetwork(file = net.path)
# ---> Loading the labels.
lab <- read.csv(lab.path)
# ---> Loading environments.
env <- read.csv(env.file)
### -------------------------- Getting Attractors -------------------------- ###
# ---> Getting attractors.
attr <- getAttractors(network = net, method = "sat.restricted", maxAttractorLength = attr.len)
# ---> Transforing attractors to data table format.
attr.table <- as.data.table(attractorToDataframe(attr=attr,Boolean=T))
### -------------------------- Labeling The Attractors -------------------------- ###
# Applying labeling rules.
labels.attr <- labelAttractors(attr, lab, net$genes, sep="/")
attr.table[,label:=labels.attr[attractor]]
# ---> Getting a summary of the attractor object.
attr.summary <- attr.table[,.(states=.N),by=.(label,attractor)][
  ,.(attractors=.N),by=.(label,states)]
fwrite(attr.table,file=paste0('data/MP_attr.csv'))
# ---> Getting a summary of the obteined attractors.
# "Label" is the tagged names by the labeling rules.
# "States" is the classification of the length of the attractors.
# "Attractors" are the number of attractors of that length.
cat('Summary of the attractors obtained:\n')
print(attr.summary)
tmp.output.file <- paste0('data/MP_attr_summary.csv')
# fwrite(x=attr.summary, file=tmp.output.file)

### ---> Simplifiying attractors
attr.table$label <- sapply(attr.table$label,simplifyLabel,replace.labels=replace.labels)

# NOTE: Removing discarded labeled attractors. These attractors will be reported in the supplementary files, but they wont be included in firther analysis
attr.table <- attr.table[!label %like% 'Discarded']

attr.summary.r <- attr.table[,.(states=.N),by=.(label,attractor)][
  ,.(attractors=.N),by=.(label,states)]
cat('Summary of the attractors obtained after replacing:\n')
print(attr.summary.r)
tmp.output.file <- paste0('data/MP_attr_summary_replacing.csv')
fwrite(x=attr.summary.r, file=tmp.output.file)


### ---> Plotting Attractors
file.attr.pdf <- "images/MP_attr_all"
# Select key nodes for plotting
node.subset <- c('IL12_out','IL10_out','IL6_out','VEGF_out','STAT1','STAT5','NFKB','STAT6','FCGR','STAT3_IL10','STAT3_IL6','SOCS1','SOCS3')
attr.table.short <- subset(attr.table, select = c("attractor","label","state", node.subset) )

# select unique attractors of subset
attrs <- list()
for (n in sort(unique(attr.table$attractor))) {
    at <- attr.table.short[attr.table.short$attractor==n,]
    at <- subset(at, select = -c(attractor,state) )
    at[nrow(at) + 1,] = NA
    row.names(at) <- NULL
    attrs[[length(attrs) + 1]] <- at
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
        start.ind <- 1 + 24; end.ind <- 24 + 24;
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


### ------------------------ Asynchronous ------------------------ ###
set.seed(534234)

initial.attr <- apply(
  X=attr.table[,!c('attractor','state','label')], MARGIN=1,
  FUN=function(x) {
    as.data.frame(t(x))
  })

pb <- txtProgressBar(min = 0,
                     max = length(initial.attr),
                     style = 3, 
                     width = 50,
                     char = "=") 

attr.asyn <- lapply(X=1:length(initial.attr), FUN=function(tmp.x){
  tmp.attr <- initial.attr[tmp.x]
#   print(tmp.x)
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
  setTxtProgressBar(pb, tmp.x)
  return(tmp.attr.table.asyn)
})
close(pb)
attr.table.asyn <- rbindlist(attr.asyn)
attr.table.asyn$label <- sapply(attr.table.asyn$label,simplifyLabel,replace.labels=replace.labels)

# Saving asyn attractators.
tmp.output.file <- paste0('data/MP_asyn_attr.csv')
fwrite(attr.table.asyn,file=tmp.output.file)

tmp.output.file <- paste0('data/MP_asyn_attr_summary.csv')
attr.summary.asyn <- attr.table.asyn[,.(states=.N),by=.(label,attractor)][
  ,.(attractors=.N),by=.(label,states)]
fwrite(x=attr.summary.asyn, file=tmp.output.file)
attr.summary.asyn

feats.of.int <- c('IL12_out','IL10_out','IL6_out','VEGF_out','STAT1','STAT5','NFKB','STAT6','FCGR','STAT3_IL10','STAT3_IL6','SOCS1','SOCS3')
tmp.feats.of.int <- c('state',feats.of.int,'HIF1a','STAT3_IL6','IFNG_out')
attr.of.int <- 'M2b*'
tmp.attr.table.asyn <- attr.table.asyn[label == attr.of.int, ..tmp.feats.of.int]
unique(tmp.attr.table.asyn)
cat('VEGF is being activated due to the IL6 out, trigerring the phenotybe being classified as M2 although it has all the M2b markers, this being cycled due SOCS3 that can inhibit IL6R->STAT3_IL6->VEGF\n\n')

# Transition between synchronus and asynchronous
tmp.file.name <- 'images/MP2_transition_syn_vs_asyn.pdf'
transition.table <- merge(x=attr.table[,.(attractor=1:.N,syn.lab=label)],
  y=unique(attr.table.asyn[,.(attractor,asyn.lab=label)]),
  by='attractor')

transition.table <- transition.table[,.(freq=.N),
  by=.(syn.lab,asyn.lab)]

# transition.table[is.na(transition.table)] <- 'NoAttractor'

pdf(tmp.file.name, height=20)
alluvial(transition.table[,.(syn.lab,asyn.lab)],
freq=transition.table$freq)
dev.off()

tmp.file.name <- paste0('data/MP_transition_table_syn_vs_asyn.csv')
fwrite(x=transition.table,file=tmp.file.name)


# NOTE: Check this function, because it should do the same as this analysis but it is integrated in BNP, it must be corrected.
# async <- verifySyncronousVsAsyncronous(net, attr, lab)
# async
### --------------------------- Mutants Analysis --------------------------- ###
file.mut <- 'data/MP_mut_label.csv'


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
mutants <- as.data.table(mutants,keep.rownames = 'label')
mutants[mutants == 0] <- NA
tmp.mutants <- mutants[1:nrow(mutants),1:ncol(mutants)]
# Simplifying labels in order to callapse the table
mutants$label <- sapply(mutants$label,simplifyLabel,replace.labels=replace.labels)
mutants <- mutants %>% group_by(label) %>% summarize_all(sum)
write.csv(x=mutants, file=file.mut, row.names=F)

mutants <- as.data.frame(mutants)
rownames(mutants) <- mutants$label; mutants$label <- NULL

file.mut.pdf = 'images/MP_mut_label.pdf'
mutants <- mutants/mutants
color <- c('#bebada')

pdf(file.mut.pdf)
heatmap(t(as.matrix( mutants )),
        main="Macrophage mutants", 
        xlab="label", ylab="Mutant",
        col=color, cexCol=0.75, cexRow=0.5,
        Colv = NA, Rowv = NA, scale="none",
        )
dev.off()

### --------------------- Environment analysis --------------------- ###

file.env.attr <- "data/MP_env_attr.csv"

env.attr <- perturbNetworkFixedNodes(net, label.rules=lab, 
                                        # returnDataFrame='basinSize', 
                                        genes  = rep( list(colnames(env)), times=nrow(env) ),
                                        values = lapply( split(env,seq_along(env[,1])), as.list),
                                        names  = rownames(env), 
                                        returnDataFrame='occurrence', 
                                        method="sat.restricted",
                                        maxAttractorLength = attr.len)

setDT(env.attr, keep.rownames = 'label')
env.attr
env.attr$label <- sapply(env.attr$label, simplifyLabel, replace.labels=replace.labels)
env.attr <- env.attr %>% group_by(label) %>% summarize_all(sum)
write.csv(env.attr, file.env.attr, row.names=F)

env.attr <- as.data.frame(env.attr)
rownames(env.attr) <- env.attr$label; env.attr$label <- NULL
env.attr[env.attr == 0] <- NA
colfunc <- colorRampPalette(c('#fee391', '#a65628'))
color <- colfunc(10)

file.env.pdf <- "images/MP_env_attr.pdf"
pdf(file.env.pdf)
heatmap(t(as.matrix( subset(env.attr ))),
        main="Macrophage by environment", 
        xlab="label", ylab="Environment",
        col=color, cexCol=0.75, cexRow=0.75,
        Colv = NA, Rowv = NA, scale="none",
       )
dev.off()

# Getting list of attractors
file.env.attr <- "data/MP_env_attr_list.csv"

env.attr <- perturbNetworkFixedNodes(net, label.rules=lab, 
                                        genes  = rep( list(colnames(env)), times=nrow(env) ),
                                        values = lapply( split(env,seq_along(env[,1])), as.list),
                                        names  = rownames(env), 
                                        returnDataFrame='attrList', 
                                        method="sat.restricted",
                                        maxAttractorLength = attr.len)
env.lab <- unlist(lapply(env.attr,labelAttractors, label.rules=lab, node.names=net$genes, sep="/"))
env.attr <- mclapply(env.attr, mc.cores=4, attractorToDataframe, Boolean=T)

# Initialize a variable to keep track of the last value in the sequence
for (i in 1:length(env.attr)){
  if (i == 1) {
    tmp.last <- tail(env.attr[[i]]$attractor,1)
    next
  }
  tmp.seq <- (tmp.last+1):(tmp.last+tail(env.attr[[i]]$attractor,1))
  env.attr[[i]]$ID <- tmp.seq[env.attr[[i]]$attractor]
  tmp.last <- tail(env.attr[[i]]$attractor,1)
}

env.attr <- rbindlist(env.attr, idcol='envriroment')
env.attr[,label:=env.lab[attractor]]
env.attr$simply_label <- sapply(env.attr$label, simplifyLabel, replace.labels=replace.labels)


# Create a list of data frames with a sequence column
list_of_data_frames <- list(
  data.frame(ID = 1:5, Value = 11:15),
  data.frame(ID = c(1:3,3,4), Value = 21:25),
  data.frame(ID = 1:6, Value = 31:36)
)

