
################################################################
#### Prosopis invasions on South American dry lands
#### Brunno Oliveira, 2015
#### Universidade Federal do Rio Grande do Norte - Brasil
################################################################


#set working directory
#setwd("/media/brunno/FAT/Prosopis")
setwd("F:/Prosopis")

gc()

rm(list=ls())

# problem loading rJava solved
# write in terminal:
# $sudo updatedb
# $locate libjvm.so
# You will see something like it:
# /usr/lib/debug/usr/lib/jvm/java-7-openjdk-amd64/jre/lib/amd64/server/libjvm.so
# /usr/lib/debug/usr/lib/jvm/java-7-openjdk-amd64/jre/lib/amd64/zero/libjvm.so
# /usr/lib/jvm/java-7-openjdk-amd64/jre/lib/amd64/jamvm/libjvm.so
# /usr/lib/jvm/java-7-openjdk-amd64/jre/lib/amd64/server/libjvm.so
# /usr/lib/jvm/java-7-openjdk-amd64/jre/lib/amd64/zero/libjvm.so
# Move the lobjvm.so to /usr/lib/
# $sudo ln -s /usr/lib/jvm/java-7-openjdk-amd64/jre/lib/amd64/server/libjvm.so /usr/lib/  

# In window environment:
# Sys.setenv(JAVA_HOME="C:\\Arquivos de Programas\\Java\\jre1.8.0_60\\")
# options(java.home="C:\\Arquivos de Programas\\Java\\jre1.8.0_60\\")

# load libraries
library(raster)
library(maptools)
library(dismo)
library(rJava)
library(ENMeval)
library(ecospat)

# If from the beginning
load('Prosopis_data.RData')
# If not
# load('Prosopis_ENMs.RData')

getFCs <- function(html) {
  htmlRead <- readLines(html)
  featureTypes <- htmlRead[grep("Feature types", htmlRead)]
  substr(featureTypes, start=21, stop=nchar(featureTypes)-4)
}

##############################################################################################
#################################### ENMs ####################################################
##############################################################################################

### Maxent modeling
# workstation
options(java.parameters = "-Xmx30g") 
# my laptop
options(java.parameters = "-Xmx7g")

Sys.setenv(NOAWT=TRUE)

### define runs and background environments
runs <- list(jul_nat=jul_nat,jul_inv=jul_inv, # juliflora
             pal_nat=pal_nat,pal_inv=pal_inv) # pallida

envs <- list(BgNatJul,BgInv,
             BgNatPal,BgInv)

# For projection
predictor <- BgAll
  
### DEFAULT maxent models

## create dir to store results of ENM Default
dir.create('Results/ENMDefault')

for(i in 1:length(runs)){
  dir.create(paste('Results/ENMDefault/',names(runs[i]),sep=''))
}

### run maxent models with default settings
for(i in 1:length(runs)) {
  
  cat('\r','run',i,'of',length(runs),'-',names(runs[i]),'\n')
  
  me <- maxent(envs[[i]],runs[[i]][,1:2], 
               path=paste('Results/ENMDefault/', names(runs[i]), sep=''), overight=T)
  r <- predict(me, predictor, overwrite=TRUE, progress='text')
  writeRaster(r, filename=paste('Results/ENMDefault/', names(runs[i]),sep=''), format="ascii", overwrite=T)
}

save.image('Prosopis_ENMs.RData')

### get FT from default models

def.results <- list(NA)

for(i in 1:length(runs)) {
  
  cat('\r','run',i,'of',length(runs),'-',names(runs[i]),'\n')
  
  def.results[[i]] <- getFCs(paste('Results/ENMDefault/',names(runs[i]),"/maxent.html",sep=''))
  def.results[[i]] <- strsplit(def.results[[i]], " ")[[1]]
}

def.results <- lapply(def.results, function(x) gsub("hinge", "H", x))
def.results <- lapply(def.results, function(x) gsub("linear", "L", x))
def.results <- lapply(def.results, function(x) gsub("product", "P", x))
def.results <- lapply(def.results, function(x) gsub("threshold", "T", x))
def.results <- lapply(def.results, function(x) gsub("quadratic", "Q", x))

def.results <- lapply(def.results, function(x) paste(x, collapse = ""))
names(def.results) <- names(runs)

write.csv(def.results, "default_results.csv")


### Evaluating ecological niche models using ENMeval 

## create dir to store results of ENMeval
dir.create('Results/ENMeval')

for(i in 1:length(runs)){
  dir.create(paste('Results/ENMeval/',names(runs[i]),sep=''))
}


### run ENMeval
coc.results <- list()

for(i in 1:length(runs)) {
  
cat('\r','run',i,'of',length(runs),'-',names(runs[i]),'\n')

coc.results[[i]] <- ENMevaluate(occ=runs[[i]][,1:2], env=envs[[i]], RMvalues=seq(0.5, 5, 0.5), 
                                fc=c('L', 'H', 'LQ', 'LQH', 'LQHP', 'LQHPT'), method='block', parallel=TRUE, numCores = 20)
}


dir.create('Results/ENMeval')


# write results ENMeval
for (i in 1:length(runs)) {
  write.csv(coc.results[[i]]@results,paste('Results/ENMeval/',names(runs[i]),'.csv',sep=''))
}

# project ENMs from ENMeval
for(i in 1:length(runs)) {
  
  cat('\r','run',i,'of',length(runs),'-',names(runs[i]),'\n')
  
  aicmods <- which(coc.results[[i]]@results$AICc == min(na.omit(coc.results[[i]]@results$AICc)))[1] # best model
  aicmods <- coc.results[[i]]@results[aicmods,]
  FC<- as.character(aicmods$features[1]) # get feature class used in the best model
  
  rm<-aicmods$rm # get RM used in the best model
  
  maxent.args <- make.args(RMvalues = rm, fc = FC)
  
  me <- maxent(envs[[i]], runs[[i]][,1:2], args=maxent.args[[1]],
               path = paste('Results/ENMeval/', names(runs[i]), sep=''), overight=T)
  r <- predict(me, predictor, overwrite=TRUE, progress = 'text')
  writeRaster(r, filename=paste('Results/ENMeval/', names(runs[i]), sep=''), format="ascii", overwrite=T)
}


save.image('Prosopis_ENMs.RData')

load('Prosopis_ENMs.RData')

##################################################
### Supp Figure
### Create plot of maps comparing default and AICc models ####

head.names1 <- c('','P. juliflora', 'P. pallida')
head.names2 <- c('', 'Default', 'AIC', 'Default', 'AIC')

run.plot <- names(runs[c(1,1,3,3,2,2,4,4)])

pdf("maps_default_vs_enmeval.pdf", width = 8, height = 4)
layout(matrix(c(1,2,2,3,3,seq(4,18)), nrow=4, byrow = T),widths = c(1,3,3,3,3),heights = c(.5,.5,3,3))
par(mar=c(0,0,0,0))

for(i in 1:3) { plot.new(); text(0.5,0.5, paste(head.names1[i]), font=3, cex=1.5) }

for(i in 1:5) { plot.new(); text(0.5,0.5, paste(head.names2[i]), cex=1.5) }

plot.new(); text(0.5,0.5, 'Native', cex=1.5)

image(raster('Results/ENMDefault/jul_nat.asc'), legend=F, axes=F, main='', box=F, col=colorRampPalette(c("#3d46a9","#49bff0","#58c53f","#f2f327", "#ec2236"))(40))
image(raster('Results/ENMeval/jul_nat.asc'), legend=F, axes=F, main='', box=F, col=colorRampPalette(c("#3d46a9","#49bff0","#58c53f","#f2f327", "#ec2236"))(40))

image(raster('Results/ENMDefault/pal_nat.asc'), legend=F, axes=F, main='', box=F, col=colorRampPalette(c("#3d46a9","#49bff0","#58c53f","#f2f327", "#ec2236"))(40))
image(raster('Results/ENMeval/pal_nat.asc'), legend=F, axes=F, main='', box=F, col=colorRampPalette(c("#3d46a9","#49bff0","#58c53f","#f2f327", "#ec2236"))(40))

plot.new(); text(0.5,0.5, 'Invasive', cex=1.5)

image(raster('Results/ENMDefault/jul_inv.asc'), legend=F, axes=F, main='', box=F, col=colorRampPalette(c("#3d46a9","#49bff0","#58c53f","#f2f327", "#ec2236"))(40))
image(raster('Results/ENMeval/jul_inv.asc'), legend=F, axes=F, main='', box=F, col=colorRampPalette(c("#3d46a9","#49bff0","#58c53f","#f2f327", "#ec2236"))(40))

image(raster('Results/ENMDefault/pal_inv.asc'), legend=F, axes=F, main='', box=F, col=colorRampPalette(c("#3d46a9","#49bff0","#58c53f","#f2f327", "#ec2236"))(40))
image(raster('Results/ENMeval/pal_inv.asc'), legend=F, axes=F, main='', box=F, col=colorRampPalette(c("#3d46a9","#49bff0","#58c53f","#f2f327", "#ec2236"))(40))

dev.off()


##################################################
### Figure 3
### Create plot of maps ONLY with AICc models ####

# Compare maps fitted with default and AIC features

head.names1 <- c('','P. juliflora', 'P. pallida')

run.plot <- names(runs[c(1,1,3,3,2,2,4,4)])

pdf("Figure3.pdf", width = 5, height = 4)
layout(matrix(c(seq(1,9)), nrow=3, byrow = T),widths = c(1,3,3),heights = c(.5,3,3))

par(mar=c(0,0,0,0))

for(i in 1:3) { plot.new(); text(0.5,0.5, paste(head.names1[i]), font=3, cex=1.5) }

plot.new(); text(0.5,0.5, 'Native', cex=1.5)

image(raster('Results/ENMeval/jul_nat.asc'), legend=F, axes=F, main='', box=F, col=colorRampPalette(c("#3d46a9","#49bff0","#58c53f","#f2f327", "#ec2236"))(40))

image(raster('Results/ENMeval/pal_nat.asc'), legend=F, axes=F, main='', box=F, col=colorRampPalette(c("#3d46a9","#49bff0","#58c53f","#f2f327", "#ec2236"))(40))

plot.new(); text(0.5,0.5, 'Invasive', cex=1.5)

image(raster('Results/ENMeval/jul_inv.asc'), legend=F, axes=F, main='', box=F, col=colorRampPalette(c("#3d46a9","#49bff0","#58c53f","#f2f327", "#ec2236"))(40))

image(raster('Results/ENMeval/pal_inv.asc'), legend=F, axes=F, main='', box=F, col=colorRampPalette(c("#3d46a9","#49bff0","#58c53f","#f2f327", "#ec2236"))(40))

dev.off()



##################################################
#### Do some basic exploration of the results ####
##################################################


# View ENMevaluation object:
tabres<- read.csv(paste('Results/ENMeval/',names(runs[i]),'.csv',sep=''))
                       
# Which settings minimum AICc?
aicmods <- which(tabres$AICc == min(na.omit(tabres$AICc)))
View(tabres[aicmods,])

# Visualize how data were partitioned
# Background points:
plot(envs[[i]][[1]])
points(coc.results[[i]]@bg.pts, col= coc.results[[i]]@bg.grp, cex=.75)

# Occurrence localities:
plot(coc.results[[i]]@predictions[[aicmods]])
points(coc.results[[i]]@occ.pts, pch=16, col= coc.results[[i]]@occ.grp, cex=.75)

# View predictions in geographic space for best and default models
par(mfrow=c(1,2))
#fxn to rescale cell values between 0 and 1
rasterRescale<-function(r){
  ((r-cellStats(r,"min"))/(cellStats(r,"max")-cellStats(r,"min")))
}

newmap<-rasterRescale(coc.results[[i]]@predictions[[aicmods]])

plot(newmap)

##################################################
### Create table with Best models

for(i in 1:length(runs)){
# View ENMevaluation object:
if(i==1){
  tabres<- read.csv(paste('Results/ENMeval/',names(runs[i]),'.csv',sep=''))
  
  # Which settings minimum AICc?
  aicmods <- which(tabres$AICc == min(na.omit(tabres$AICc)))
  aicmods <- tabres[aicmods,]
  }
  if(i>1){
    tabres<- read.csv(paste('Results/ENMeval/',names(runs[i]),'.csv',sep=''))
    
    # Which settings minimum AICc?
    x <- which(tabres$AICc == min(na.omit(tabres$AICc)))
    x <- tabres[x,]
    aicmods <- rbind(aicmods, x)
  }
}

aicmods <- cbind(run = names(runs), aicmods)
View(aicmods)
write.csv(aicmods, "AICc_results.csv")

##################################################
### Plot ENMeval models comparizons for each run

for(i in 1:length(runs)){
pdf(paste('Results/ENMeval/',names(runs[i]),'.pdf',sep=''), width = 8, height = 5.5)
par(mfrow=c(2,3))
eval.plot(coc.results[[i]]@results)
eval.plot(coc.results[[i]]@results, 'Mean.AUC', legend = F)
eval.plot(coc.results[[i]]@results, 'Mean.AUC.DIFF', legend = F)
eval.plot(coc.results[[i]]@results, 'Mean.OR10', legend = F)
eval.plot(coc.results[[i]]@results, 'Mean.ORmin', legend = F)

# Note that using the results table for plotting enables creativity via the full range of R plotting options. For example:
# Let's see how AUCtrain values compare with delta.AICc, depending on feature classes and regularization multiplier settings...
plot(coc.results[[i]]@results$Mean.AUC, coc.results[[i]]@results$delta.AICc, bg=coc.results[[i]]@results$features, pch=21, 
     cex= coc.results[[i]]@results$rm/2, ylim=c(0,30), xlab = "Mean.AUC", ylab = 'delta.AICc', cex.lab = 1.5)
#legend("topright", legend=unique(coc.results[[i]]@results$features), pt.bg=coc.results[[i]]@results$features, pch=21)
#mtext("Circle size proportional to regularization multiplier value")
dev.off()
}



##################################################
### Plot delta.AICc for different settings
### Supp Figure
### multiple graphs

pdf(paste('ENMeval_res','.pdf',sep=''), width = 7, height = 6.5)
layout(matrix(c(1,2,3,4,5,6,7,8,9), nrow=3, byrow = T),widths = c(1.5,3,3),heights = c(.5,3,3))
par(mar=c(0,0,0,0),mai=c(0,0,0,0))
plot.new(); text(0.5,0.5, '', cex=1.5)
plot.new(); text(0.5,0.5, "Prosopis juliflora", font=3, cex=1.5)
plot.new(); text(0.5,0.5, "Prosopis pallida", font=3, cex=1.5)
plot.new(); text(0.5,0.5, 'Native', cex=1.5)
eval.plot(coc.results[[1]]@results, legend = F)
eval.plot(coc.results[[3]]@results)
plot.new(); text(0.5,0.5, 'Invasive', cex=1.5)
eval.plot(coc.results[[2]]@results, legend = F)
eval.plot(coc.results[[4]]@results, legend = F)
dev.off()


###################################################