# Main: MP_attractor.tmp.R
# Date: 11-06-2023
# Update:
# The output file created by sink() is conditional on interactive jobs only.
if (!interactive()){
  out.file <- 'scripts/MP_attractors.out'
  sink(out.file)
}
### -------------------------- Loading Libraries -------------------------- ###
cat('### -------------------------- Loading Libraries -------------------------- ###\n')
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
source('BNP/BNP_Label.R')
source('BNP/BNP_PerDynamic.R')
source('BNP/BNP_PerFunctions.R')

cat('The libraries have been loaded successfully!\n')
### ------------------------------ Loading Functions ------------------------------ ###
cat('### ------------------------------ Loading Functions ------------------------------ ###\n')
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
# ---> Function to remove all not used variables.
cleanUp <- function(vars.to.keep){
  vars.to.remove <- ls(envir = .GlobalEnv)[!ls(envir = .GlobalEnv) %chin% vars.to.keep]
  rm(list=vars.to.remove, envir = .GlobalEnv)
  cat('Environment has been cleaned up...\n')
}

cat('The functions have been loaded successfully!\n')
### -------------------------- Arguments -------------------------- ###
cat('### -------------------------- Arguments -------------------------- ###\n')
# ---> File of the network rules to use.
net.file <- "model/MP_full_boolean_network.csv"
# ---> File of the labeling rules to use.
lab.file <- "model/MP_label_rules.tmp.csv"
# ---> Maximun size of the attractors.
attr.len <- 3
# ---> Environments
env.file <- "model/MP_M1_environments.csv"
# ---> Simplied labels.
replace.labels <- c('M2/M2b'='M2b','M0/M2/M2c'='M2c','il6/M2d'='M2d-like','M1/M2d'='M2d-like','M1/M2/M2d'='M2d-like','M0/M1/M2d'='M2d-like','M0/M1/M2'='Discarded','M0/M2'='M2c','M0/M1'='M1','il6/M2c/M2d'='M2d-like','M1/M2'='Discarded','M1/M2/NoLabel'='Discarded')

if(!file.exists(net.file)) stop('The path to the file with the network does not exist!!!')
if(!file.exists(lab.file)) stop('The path to the file with the label rules does not exist!!!')
if(!file.exists(env.file)) stop('The path to the file with the environments does not exist!!!')

cat('Boolean network to use: ',net.file,'\n')
cat('Label rules to use: ',lab.file,'\n')
cat('Maximum length of the attractors to calculate: ',attr.len,'\n')
cat('Enviroments to use: ',env.file,'\n')

### -------------------------- Loading data -------------------------- ###
cat('### -------------------------- Loading data -------------------------- ###\n')
# ---> Loading the network.
net <- loadNetwork(file = net.file)
# ---> Loading the labels.
lab <- read.csv(lab.file)
# ---> Loading environments.
env <- read.csv(env.file)
cat('Data loaded successfully!\n')
base.vars <- c(ls(),'base.vars')
### --------------------------  Attractors analysis -------------------------- ###
cat('### -------------------------- Attractors analysis -------------------------- ###\n')
attr.file <- 'data/MP_attr.csv'
if (!file.exists(attr.file)){
  # ---> Getting attractors.
  attr <- getAttractors(network = net, method = "sat.restricted", maxAttractorLength = attr.len)
  # ---> Transforing attractors to data table format.
  attr.table <- as.data.table(attractorToDataframe(attr=attr,Boolean=T))
  ### ---> Labeling The Attractors
  # Applying labeling rules.
  labels.attr <- labelAttractors(attr, lab, net$genes, sep="/")
  attr.table[,label:=labels.attr[attractor]]
  ### ---> Simplifiying attractors
  attr.table$simply_label <- sapply(attr.table$label,simplifyLabel,replace.labels=replace.labels)
  # Saving table of labeled attractors.
  fwrite(attr.table,file=attr.file)
  # NOTE: Removing discarded labeled attractors. These attractors will be reported in the supplementary files, but they wont be included in firther analysis
  attr.table <- attr.table[!simply_label %like% 'Discarded']
  cat('Attractors have been calculated...\n')
} else {
  attr.table <- fread(file=attr.file)
  # NOTE: Removing discarded labeled attractors. These attractors will be reported in the supplementary files, but they wont be included in firther analysis
  attr.table <- attr.table[!simply_label %like% 'Discarded']
  cat('Attractors have been loaded...\n')
}
print(attr.table)
cat('\n')
# ---> Getting a summary of the attractor object.
attr.summary.file <- 'data/MP_attr_summary.csv'
if(!file.exists(attr.summary.file)){
  attr.summary <- attr.table[,.(states=.N),by=.(label,attractor)][
  ,.(attractors=.N),by=.(label,states)]
  fwrite(x=attr.summary, file=attr.summary.file)
  cat('Attractors summary have been generated..\n')
} else {
  attr.summary <- fread(attr.summary.file)
  cat('Attractors summary have been loaded..\n')
}
# "Label" is the tagged names by the labeling rules.
# "States" is the classification of the length of the attractors.
# "Attractors" are the number of attractors of that length.
cat('Summary of the attractors obtained:\n')
print(attr.summary)
cat('\n')

# ---> Getting a summary of the attractor object: Simplified labels.
attr.summary.r.file <- 'data/MP_attr_summary_simplified.csv'
if(!file.exists(attr.summary.r.file)){
  attr.summary.r <- attr.table[,.(states=.N),by=.(simply_label,attractor)][
  ,.(attractors=.N),by=.(simply_label,states)]
  attr.summary.r <- attr.summary.r[order(states)]
  fwrite(x=attr.summary.r, file=attr.summary.r.file)
  cat('Attractors summary with simplified labels have been generated..\n')
} else {
  attr.summary.r <- fread(attr.summary.r.file)
  cat('Attractors summary with simplified labels have been loaded..\n')
}
cat('Summary of the attractors obtained after replacing:\n')
print(attr.summary.r)
cat('\n')


### ---> Plotting Attractors
file.attr.pdf <- "images/MP_attr_all"
if(!all(file.exists(paste0(file.attr.pdf,'_',1:attr.len,'.pdf')))){
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
} 
cat('Attractors have been plotted...\n')

cat('Attractor analysis has been successfully completed!\n\n')
base.vars <- c(base.vars,'attr.table')
cleanUp(base.vars)
### ------------- Asynchronous network analysis ------------- ###
cat('### ------------- Asynchronous network analysis ------------- ###\n')
seed.to.use <- 534234
cat('Setting seed for analysis replication: ',seed.to.use,'\n')
set.seed(seed.to.use)

# ---> Calculating asynchronous attractors from synchronous atroctors.
asyn.attr.file <- 'data/MP_asyn_attr.csv'
if(!file.exists(asyn.attr.file)){
  # Getting list of synchronous attractors.
  initial.attr <- apply(
    X=attr.table[,!c('attractor','state','label','simply_label')], MARGIN=1,
    FUN=function(x) {
      as.data.frame(t(x))
    })
  # Start progress bar.
  pb <- txtProgressBar(min = 0,max = length(initial.attr),style = 3,width = 50,char = "=") 

  # Applying asynchronous network to each of the synchronous attractors.
  attr.asyn <- lapply(X=1:length(initial.attr), FUN=function(tmp.x){
    tmp.attr <- initial.attr[tmp.x]
    # If after 5 seconds it is not possible to obtain an asynchronous attractor, move on to the next attractor.
    tmp.attr.asyn <- tryCatch(expr= {
      withTimeout({getAttractors(
        network = net, method = "chosen",
        type='asynchronous',startStates=tmp.attr)}, timeout=5)
    },
      error=function(i) NA,
      TimeoutException = function(ex) NA)
    # Returning asynchronous attractor in data table format.
    if(all(!is.na(tmp.attr.asyn))){
      tmp.attr.table.asyn <- as.data.table(attractorToDataframe(attr=tmp.attr.asyn,Boolean=T))
      # Labeling asynchronous attractor.
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
    # Progress bar tracking.
    setTxtProgressBar(pb, tmp.x)
    return(tmp.attr.table.asyn)
  })
  # Closing the progress bar.
  close(pb)
  # Joining list of attractors in a single data table.
  attr.table.asyn <- rbindlist(attr.asyn)
  # Simplifying the labels of asynchronous attractors.
  attr.table.asyn$simply_label <- sapply(attr.table.asyn$label,simplifyLabel,replace.labels=replace.labels)
  # Saving asynchronous attractators.
  fwrite(attr.table.asyn,file=asyn.attr.file)
  cat('Asynchronous attractors have been calculated...\n')
} else {
  attr.table.asyn <- fread(asyn.attr.file)
  cat('Asynchronous attractors have been loaded...\n')
}
print(attr.table.asyn)
cat('\n')
# ---> Calculating summary of asynchronous attractors.
asyn.attr.summ.file <- 'data/MP_asyn_attr_summary.csv'
if(!file.exists(asyn.attr.summ.file)){
  attr.summary.asyn <- attr.table.asyn[,.(states=.N),by=.(label,attractor)][
    ,.(attractors=.N),by=.(label,states)]
  attr.summary.asyn <- attr.summary.asyn[order(states,label)]
  fwrite(x=attr.summary.asyn, file=asyn.attr.summ.file)
  cat('Summary of asynchronous attractors has been generated...\n')
} else {
  attr.summary.asyn <- fread(asyn.attr.summ.file)
  cat('Summary of asynchronous attractors has been loaded...\n')
}
cat('Summary of asynchronous attractors:\n')
print(attr.summary.asyn)
cat('\n')

# ---> Calculating summary of asynchronous attractors with simplified labels.
asyn.attr.summ.r.file <- 'data/MP_asyn_attr_summary_simplified.csv'
if(!file.exists(asyn.attr.summ.r.file)){
  attr.summary.asyn.r <- attr.table.asyn[,.(states=.N),by=.(simply_label,attractor)][
    ,.(attractors=.N),by=.(simply_label,states)]
  attr.summary.asyn.r <- attr.summary.asyn.r[order(states,simply_label)]
  fwrite(x=attr.summary.asyn.r, file=asyn.attr.summ.r.file)
  cat('Summary of asynchronous attractors with simplified labels has been generated...\n')
} else {
  attr.summary.asyn.r <- fread(asyn.attr.summ.r.file)
  cat('Summary of asynchronous attractors with simplified labels has been loaded...\n')
}
cat('Summary of asynchronous attractors:\n')
print(attr.summary.asyn)
cat('\n')

feats.of.int <- c('IL12_out','IL10_out','IL6_out','VEGF_out','STAT1','STAT5','NFKB','STAT6','FCGR','STAT3_IL10','STAT3_IL6','SOCS1','SOCS3')
tmp.feats.of.int <- c('state',feats.of.int,'HIF1a','STAT3_IL6','IFNG_out')
attr.of.int <- 'M2b*'
tmp.attr.table.asyn <- attr.table.asyn[simply_label == attr.of.int, ..tmp.feats.of.int]
unique(tmp.attr.table.asyn)
cat('VEGF is being activated due to the IL6 out, trigerring the phenotybe being classified as M2 although it has all the M2b markers, this being cycled due SOCS3 that can inhibit IL6R -> STAT3_IL6 -> VEGF\n\n')

# ---> Creating alluvial transition between synchronous and asynchronous attractors.
syn.asyn.table.file <- 'data/MP_transition_table_syn_vs_asyn.csv'
if(!file.exists(syn.asyn.table.file)){
  # Creating transition table between synchronous and asynchronous attractors.
  transition.table <- merge(x=attr.table[,.(attractor=1:.N,syn.lab=simply_label)],
    y=unique(attr.table.asyn[,.(attractor,asyn.lab=simply_label)]),
    by='attractor')
  transition.table <- transition.table[,.(freq=.N),
    by=.(syn.lab,asyn.lab)]
  transition.table <- transition.table[order(syn.lab,asyn.lab)]
  # Saving transition table.
  fwrite(x=transition.table,file=syn.asyn.table.file)
  cat('Transition table between synchronous and asynchronous attractors has been generated...\n')
} else {
  transition.table <- fread(syn.asyn.table.file)
  cat('Transition table between synchronous and asynchronous attractors has been loaded...\n')
}
cat('Transition table:\n')
print(transition.table)
cat('\n')

# Plotting the transition table in an alluvial
syn.asyn.plot.file <- 'images/MP_transition_syn_vs_asyn.pdf'
if(!file.exists(syn.asyn.plot.file)){
  # Identifier colors of asynchronous labels.
  lab.col <- c('M0'='#a4a5a4','M1'='#fbcb98','M2'='#afc2ff','M2a'='#cffecd','M2b'='#cdfefe','M2b*'='#cdfefe','M2c'='#ffc0cb',
  'M2d'='#cc98fe','il6'='#939292')
  lab.col <- lab.col[transition.table$asyn.lab]
  # Plotting alluvial.
  pdf(syn.asyn.plot.file, height=20)
  alluvial(transition.table[,.(syn.lab,asyn.lab)],
    freq=transition.table$freq,
    col=lab.col)
  dev.off()
  cat('Alluvial transition between synchronous and asynchronous attractors has been generated...\n')
}

cat('Asynchronous network analysis has been successfully completed!\n\n')
cleanUp(base.vars)

# NOTE: Check this function, because it should do the same as this analysis but it is integrated in BNP, it must be corrected.
# async <- verifySyncronousVsAsyncronous(net, attr, lab)
# async
### --------------------- Environment analysis --------------------- ###
cat('### --------------------- Environment analysis --------------------- ###')
cat('Environments to use:\n')
print(env)
cat('\n')

# ---> Calculating the number of attractors for each of the environments.
env.attr.file <- "data/MP_env_attr.csv"
if(!file.exists(env.attr.file)){
  env.attr <- perturbNetworkFixedNodes(net, label.rules=lab, 
                                        genes  = rep( list(colnames(env)), times=nrow(env) ),
                                        values = lapply( split(env,seq_along(env[,1])), as.list),
                                        names  = rownames(env), 
                                        returnDataFrame='occurrence', 
                                        method="sat.restricted",
                                        maxAttractorLength = attr.len)

  setDT(env.attr, keep.rownames = 'label')
  env.attr
  env.attr <- env.attr %>% group_by(label) %>% summarize_all(sum)
  fwrite(x=env.attr, file=env.attr.file)
  cat('Number of attractors generated in the environments has been calculated...\n')
} else {
  env.attr <- fread(env.attr.file)
  cat('Number of attractors generated in the environments has been loaded...\n')
}
cat('Number of attractors in the environments:\n')
env.attr
cat('\n')

# ---> Calculating the number of attractors for each of the environments with simplified labels.
env.attr.r.file <- "data/MP_env_attr_simply.csv"
if(!file.exists(env.attr.r.file)){
  env.attr$simply_label <- sapply(env.attr$label, simplifyLabel, replace.labels=replace.labels)
  env.attr$label <- NULL
  env.attr <- env.attr %>% group_by(simply_label) %>% summarize_all(sum)
  fwrite(x=env.attr, file=env.attr.r.file)
  cat('Number of attractors generated in the environments with simplified labels has been calculated...\n')
} else {
  env.attr <- fread(env.attr.r.file)
  cat('Number of attractors generated in the environments with simplified labels has been loaded...\n')
}
cat('Number of attractors in the environments with simplified labels:\n')
print(env.attr)
cat('\n')

# ---> Plotting a heatmap with the number of attractors obtained in the environments with the simplified labels.
env.heatmap.file <- "images/MP_env_attr.pdf"
if(!file.exists(env.heatmap.file)){
  env.attr <- as.data.frame(env.attr)
  rownames(env.attr) <- env.attr$simply_label; env.attr$simply_label <- NULL
  env.attr[env.attr == 0] <- NA
  colfunc <- colorRampPalette(c('#fee391', '#a65628'))
  color <- colfunc(10)
  # Plotting heatmap.
  pdf(env.heatmap.file)
  heatmap(t(as.matrix( subset(env.attr ))),
        main="Macrophage by environment", 
        xlab="label", ylab="Environment",
        col=color, cexCol=0.75, cexRow=0.75,
        Colv = NA, Rowv = NA, scale="none",
        )
  dev.off()
  cat('Heatmap with the number of attractors obtained in the environments with simplified labels has been generated...\n\n')
}

# ---> Calculating the list of attractors obtained in the environments.
env.attr.list.file <- "data/MP_env_attr_list.csv"
if(!file.exists(env.attr.list.file)){
  env.attr.list <- perturbNetworkFixedNodes(net, label.rules=lab, 
                                          genes  = rep( list(colnames(env)), times=nrow(env) ),
                                          values = lapply( split(env,seq_along(env[,1])), as.list),
                                          names  = rownames(env), 
                                          returnDataFrame='attrList', 
                                          method="sat.restricted",
                                          maxAttractorLength = attr.len)
  env.lab <- unlist(lapply(env.attr.list,labelAttractors, label.rules=lab, node.names=net$genes, sep="/"))
  env.attr.list <- mclapply(env.attr.list, mc.cores=4, attractorToDataframe, Boolean=T)

  # Initialize a variable to keep track of the last value in the attractor sequence.
  for (i in 1:length(env.attr.list)){
    if (i == 1) {
      tmp.last <- tail(env.attr.list[[i]]$attractor,1)
      next
    }
    tmp.seq <- (tmp.last+1):(tmp.last+tail(env.attr.list[[i]]$attractor,1))
    env.attr.list[[i]]$attractor <- tmp.seq[env.attr.list[[i]]$attractor]
    tmp.last <- tail(env.attr.list[[i]]$attractor,1)
  }
  # Formatting table with the list of attractors, and applying the simplified labels.
  env.attr.list <- rbindlist(env.attr.list, idcol='environment')
  env.attr.list[,label:=env.lab[attractor]]
  env.attr.list$simply_label <- sapply(env.attr.list$label, simplifyLabel, replace.labels=replace.labels)
  # Saving the list of attractors.
  fwrite(x=env.attr.list, file=env.attr.list.file)
  cat('The list of attractors obtained in the environments has been calculated...\n')
} else {
  env.attr.list <- fread(env.attr.list.file)
  cat('The list of attractors obtained in the environments has been loaded...\n')
}
cat('List of attractors obtained in the environments:\n')
print(env.attr.list)
cat('\n')

cat('Environment analysis has been successfully completed!\n\n')
cleanUp(base.vars)
### --------------------------- Mutants Analysis --------------------------- ###
cat('### --------------------------- Mutants Analysis --------------------------- ###\n')
# ---> Obtaining the number of attractors of each mutated gene.
mut.file <- 'data/MP_mut_label.csv'
if(!file.exists(mut.file)){
  # Making mutant analysis over an iteration of each gene.
  # NOTE: Error generated by canonical why is this happening? Ask if deactivating is ok
  mutants <- mclapply(X=net$genes, mc.cores=50, FUN=function(tmp.gene){
      tmp.mutants <- perturbNetworkFixedNodes(net, label.rules=lab,
        genes=rep(tmp.gene,2), values=c(0,1),
        returnDataFrame='occurrence', method="sat.restricted",
        maxAttractorLength = attr.len, canonical=F)
      return(tmp.mutants)
      })
  # Merging list of mutants into a data table.
  mutants <- Reduce(function(x, y){
    tmp.merged <- merge(x, y, by = 0, all = TRUE)
    rownames(tmp.merged) <- tmp.merged$Row.names
    tmp.merged$Row.names <- NULL
    tmp.merged
    }, 
    mutants)
  mutants <- as.data.table(mutants,keep.rownames = 'label')
  # Calculating number of attractors per label.
  mutants <- as.data.table(mutants %>% group_by(label) %>% summarize_all(sum))
  # Saving table of mutants.
  fwrite(x=mutants, file=mut.file)
  cat('Number of attractors of each mutated gene has been calculated...\n')
} else {
  mutants <- fread(mut.file)
  cat('Number of attractors of each mutated gene has been loaded...\n')
}
cat('Number of attractors of each mutated gene:\n')
print(mutants)
cat('\n')

# ---> Obtaining the number of attractors of each mutated gene with simplified labels.
mut.r.file <- 'data/MP_mut_label_simply.csv'
if(!file.exists(mut.r.file)){
  # Applying simplified labels.
  mutants$simply_label <- sapply(mutants$label,simplifyLabel,replace.labels=replace.labels)
  mutants$label <- NULL
  # Summarizing attractors by simplified label.
  mutants <- as.data.table(mutants %>% group_by(simply_label) %>% summarize_all(sum))
  # Saving table of mutants.
  fwrite(x=mutants, file=mut.r.file)
  cat('Number of attractors of each mutated gene with simplified labels has been calculated...\n')
} else {
  mutants <- fread(mut.r.file)
  cat('Number of attractors of each mutated gene with simplified labels has been loaded...\n')
}
cat('Number of attractors of each mutated gene with simplified labels:\n')
print(mutants)
cat('\n')

# ---> Plotting a heatmap with the number of attractors obtained of each mutated gene with the simplified labels.
mut.heatmap.file  <-  'images/MP_mut_label.pdf'
if(!file.exists(mut.heatmap.file)){
  # Format tablet to be able to plot as heatmap. All numbers are normalized to 1.
  mutants[mutants == 0] <- NA
  mutants <- as.data.frame(mutants)
  rownames(mutants) <- mutants$simply_label; mutants$simply_label <- NULL
  mutants <- mutants/mutants
  color <- c('#bebada')
  # Plotting heatmap.
  pdf(mut.heatmap.file)
  heatmap(t(as.matrix( mutants )),
          main="Macrophage mutants", 
          xlab="label", ylab="Mutant",
          col=color, cexCol=0.75, cexRow=0.5,
          Colv = NA, Rowv = NA, scale="none",
          )
  dev.off()
  cat('Heatmap with the number of attractors obtained of each mutated gene with simplified labels has been generated...\n\n')
}

# ---> Calculating the list of attractors obtained from each mutated gene.
mutants.list.file <- 'data/MP_mut_label_list.csv'
if(!file.exists(mutants.list.file)){
  # Making mutant analysis over an iteration of each gene.
  # NOTE: Error generated by canonical why is this happening? Ask if deactivating is ok
  mut.vec <- c(paste(net$genes,'1'),paste(net$genes,'0')); names(mut.vec) <- mut.vec
  mutants.list <- mclapply(X=mut.vec, mc.cores=50,
    FUN=function(tmp.string){
      tmp.gene <- unlist(str_split(tmp.string,' '))[1]
      tmp.mut <- unlist(str_split(tmp.string,' '))[2]
      tmp.mutants <- perturbNetworkFixedNodes(net, label.rules=lab,
        genes=tmp.gene, values=tmp.mut,
        returnDataFrame='attrList', method="sat.restricted",
        maxAttractorLength = attr.len, canonical=F)
      tmp.mutants <- tmp.mutants[[1]]
      return(tmp.mutants)
      })
  mutants.list.backup <- mutants.list[1:length(mutants.list)]

  mutants.lab <- unlist(lapply(mutants.list,labelAttractors, label.rules=lab, node.names=net$genes, sep="/"))
  mutants.list <- mclapply(mutants.list, mc.cores=50, attractorToDataframe, Boolean=T)
  # Initialize a variable to keep track of the last value in the attractor sequence.
  for (i in 1:length(mutants.list)){
    if (i == 1) {
      tmp.last <- tail(mutants.list[[i]]$attractor,1)
      next
    }
    tmp.seq <- (tmp.last+1):(tmp.last+tail(mutants.list[[i]]$attractor,1))
    mutants.list[[i]]$attractor <- tmp.seq[mutants.list[[i]]$attractor]
    tmp.last <- tail(mutants.list[[i]]$attractor,1)
  }
  # Formatting table with the list of attractors, and applying the simplified labels.
  mutants.list <- rbindlist(mutants.list, idcol='mutant')
  mutants.list[,label:=mutants.lab[attractor]]
  mutants.list$simply_label <- sapply(mutants.list$label, simplifyLabel, replace.labels=replace.labels)
  # Saving the list of attractors.
  fwrite(x=mutants.list, file=mutants.list.file)
  cat('The list of attractors obtained from each mutated gene has been calculated...\n')
} else {
  mutants.list <- fread(mutants.list.file)
  cat('The list of attractors obtained from each mutated gene has been loaded...\n')
}
cat('List of attractors obtained from each mutated gene:\n')
print(mutants.list)
cat('\n')
sessionInfo()

if (!interactive()) sink()