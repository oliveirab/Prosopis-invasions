################################################################
#### Prosopis invasions on South American dry lands
#### Brunno Oliveira, 2015
#### Universidade Federal do Rio Grande do Norte - Brasil
################################################################

#### 1) LOAD FILES (PREDICTORS AND SPPS OCCURRENCES)
#### 2) COMPARE NATIVE AND INVADED NICHES USING PCA
#### 3) 

#set working directory
setwd("/media/brunno/FAT/Prosopis")
#setwd("F:/Prosopis")

gc()

rm(list=ls())

#load First time
#load('Prosopis_data.RData')

#load Not first time
load('Prosopis_analysis.RData')


# problem loading rJava solved
# I just went to Software center and intalled rjava
# In other situation this worked:
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


### Shape files
# WE - polygon Wesern Hemisphere
# WE_c -  polygon Wesern Hemisphere contour
# pal.nat.bg.shp - P. pallida Bg sensu Pasiecznik 2001
# jul.nat.bg.shp - P. juliflora Bg following Pasiecznik 2001
# NatFrame - Native Bg (union pal.nat.bg.shp and jul.nat.bg.shp)
# InvFrame - Invasive Bg (Brazil)
# AllFrame - All Bg (crop WE from Soulth of USA to Soulth Brazil)
# World - World shp

### climatic predictors

# BgAll
# BgNat
# BgNatJul
# BgNatPal
# BgInv

# resolution ~ 0.04166Dg ~ 5km

### Species occurrence data

#sp.prosopis
#jul_nat
#pal_nat
#jul_inv_raw
#pal_inv_raw
#jul_inv
#pal_inv

#################################################################################################
################################# ANALYSIS - selection of parameters  ###########################

### ANALYSIS - selection of parameters 

# selection of the type of analysis.
# If PROJ =F, the models are calibrated on both ranges.
# If PROJ =T, the models are calibrated on species 1 range only and projected to range 2. 
PROJ = F

#number of interation for the tests of equivalency and similarity
iterations<-100

#resolution of the gridding of the climate space
R=100


#################################################################################################
################################# ANALYSIS #################################

# Comparizons:
# First we compare Prosopis complex between ranges to test if the species show niche conservatism or niche shifts
# Secondly, we test the same prediction separatly for P. juliflora and P. pallida
# Thirdly, we compare differences between species within their native and invaded ranges.
# Finally, we compare whether the misidentification in P. juliflora generated different models.

# 1) P.juliflora-P.pallida complex between native and invaded ranges
# 2) P.juliflora between native and invaded ranges
# 3) P.juliflora between native and invaded ranges (raw)
# 4) P.pallida between native and invaded ranges


runs <- list(sp.prosopis=sp.prosopis, #1
             sp.juliflora=sp.juliflora, #2
             sp.juliflora.raw=sp.juliflora.raw, #3
             sp.pallida=sp.pallida) #4
             

clims1 <- list(sampBgNat, #1
               sampBgNatJul, #2
               sampBgNatJul, #3
               sampBgNatPal) #4
               
               

clims2 <- list(sampBgInv, #1
               sampBgInv, #2
               sampBgInv, #3
               sampBgInv) #4
               
               
               
### Create table for store results

res.col <- c('unfilling','stability','expansion','I','D',
             'identity (I)','similarity A->B (I)','similarity A<-B (I)',
             'identity (D)','similarity A->B (D)','similarity A<-B (D)',
             'PCA')
res.row <- c('Prosopis spp (native vs invasive)', 
             'P. juliflora (native vs invasive)', 
             'P. juliflora (native vs invasive raw)', 
             'P. pallida (native vs invasive)', 
             'Native range  (P. juliflora vs P. pallida)', 
             'Invasive range  (P. juliflora vs P. pallida)',
             'Invasive range - raw vs reanalyzed (P. juliflora)')

results <- as.data.frame(matrix(ncol=length(res.col),nrow=length(res.row)))
colnames(results) <- res.col
rownames(results) <- res.row

### Create dir to stores results and graphs

dir.create('Results/identity_test')
dir.create('Results/overlap_graphs')
dir.create('Results/PCA_graphs')
dir.create('Results/similarity_test_A-B')
dir.create('Results/similarity_test_B-A')
dir.create('Results/summary_overlap')

################################# 1) Comparison between ranges  #################################

for (i in 1:length(runs)) {
  
cat('\r','round',i,'of',length(runs),'-',"piura",'\n')

spp <- runs[[i]]
clim1<-clims1[[i]]
clim2<-clims2[[i]]
clim12<-rbind(clim1,clim2)

# sample environmental values for species occurences
occ.sp1<-na.exclude(ecospat.sample.envar(dfsp=spp,colspxy=1:2,colspkept=1:3,dfvar=clim1,
                                        colvarxy=1:2,colvar="all",resolution=resolution))
occ.sp2<-na.exclude(ecospat.sample.envar(dfsp=spp,colspxy=1:2,colspkept=1:3,dfvar=clim2,
                                         colvarxy=1:2,colvar="all",resolution=resolution))

# selection of variables to include in the analyses
Xvar<-names(BgAll)
nvar<-length(Xvar)

#################################################################################################
#################################### PCA-ENV ####################################################
#################################################################################################

row.w.1.env<-1-(nrow(clim1)/nrow(clim12))  # prevalence of clim1
row.w.2.env<-1-(nrow(clim2)/nrow(clim12))  # prevalence of clim2
row.w.env<-c(rep(row.w.1.env, nrow(clim1)),rep(row.w.2.env, nrow(clim2)),rep(0, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))

# global dataset for the analysis and rows for each sub dataset
data.env.occ<-rbind(clim1[Xvar],clim2[Xvar],occ.sp1[Xvar],occ.sp2[Xvar])
row.clim1<-1:nrow(clim1)
row.clim2<-(nrow(clim1)+1):(nrow(clim1)+nrow(clim2))
row.clim12<-1:(nrow(clim1)+nrow(clim2))
row.sp1<-(nrow(clim1)+nrow(clim2)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1))
row.sp2<-(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+nrow(occ.sp2))

# measures niche overlap along the two first axes of a PCA calibrated on all the pixels of the study areas
pca.cal <-dudi.pca(data.env.occ,row.w = row.w.env, center = T, scale = T, scannf = F, nf = 2)

# measure % explained by the first two PCs
pca.scores <- 100 * pca.cal$eig/sum(pca.cal$eig)
pca.scores <- pca.scores[1] + pca.scores[2]

results$PCA[i] <- round(pca.scores, 2)

# predict the scores on the axes
scores.clim12<- pca.cal$li[row.clim12,]
scores.clim1<- pca.cal$li[row.clim1,]
scores.clim2<- pca.cal$li[row.clim2,]
scores.sp1<- pca.cal$li[row.sp1,]
scores.sp2<- pca.cal$li[row.sp2,]

# calculation of occurence density and test of niche equivalency and similarity 
z1<- ecospat.grid.clim.dyn(scores.clim12,scores.clim1,scores.sp1, th.sp=0, R)
z2<- ecospat.grid.clim.dyn(scores.clim12,scores.clim2,scores.sp2, th.sp=0, R)

# Plot niche categories and species density created by ecospat.grid.clim.dyn
pdf(paste('Results/PCA_graphs/',res.row[i],'_PCA','.pdf',sep=''), width = 6, height = 5.5)
plot.PCA(pca.cal$co,pca.cal$eig)
dev.off()

# Calculate niche expansion, st <- stability and unfilling
dyn.indexes <- ecospat.niche.dyn.index (z1, z2, intersection=NA)$dynamic.index.w
results[i,1] <- dyn.indexes[[3]]
results[i,2] <- dyn.indexes[[2]]
results[i,3] <- dyn.indexes[[1]]

# Plot niche categories and species density created by ecospat.grid.clim.dyn
pdf(paste('Results/overlap_graphs/',res.row[i],'.pdf',sep=''), width = 7*.7, height = 5)
ecospat.plot.niche.dyn (z1, z2,quant=.25,interest=1,title=res.row[i])
legend('topleft',legend=c(paste("Unfilling:",round(results[i,1],3)),
                        paste("Stability:",round(results[i,2],3)),
                        paste("Expansion:",round(results[i,3],3))),
                        bty='n')
dev.off()


# calculation of niche overlap - I and D indexes
niche_over <- ecospat.niche.overlap (z1, z2, cor=F)

results$I[i] <- niche_over$I
results$D[i] <- niche_over$D

# test of niche equivalency and similarity according to Warren et al. 2008
cat('\n','running niche equivalency test')
a<-ecospat.niche.equivalency.test(z1, z2, rep=100) 
cat('\n','running niche similarity test A->B','\n')
b<-ecospat.niche.similarity.test(z1, z2, rep=100)
cat('\n','running niche similarity test A<-B','\n')
b2<-ecospat.niche.similarity.test(z2, z1, rep=100)


results[i,6] <- a$p.I
results[i,9] <- a$p.D

results[i,7] <- b$p.I
results[i,8] <- b2$p.I

results[i,10] <- b$p.D
results[i,11] <- b2$p.D


write.csv(a$sim$I,paste(sep='','Results/identity_test/',res.row[i],'_null_ident_I.csv')) # null distribuion Identity test - I
write.csv(a$sim$D,paste(sep='','Results/identity_test/',res.row[i],'_null_ident_D.csv')) # null distribuion Identity test - D

write.csv(b$sim$I,paste(sep='','Results/similarity_test_A-B/',res.row[i],'_null_sim_A-B_I.csv')) # null distribuion Similarity test A -> B - I
write.csv(b$sim$D,paste(sep='','Results/similarity_test_A-B/',res.row[i],'_null_sim_A-B_D.csv')) # null distribuion Similarity test A -> B - D

write.csv(b2$sim$I,paste(sep='','Results/similarity_test_B-A/',res.row[i],'_null_sim_B-A_I.csv')) # null distribuion Similarity test A <- B - I
write.csv(b2$sim$D,paste(sep='','Results/similarity_test_B-A/',res.row[i],'_null_sim_B-A_D.csv')) # null distribuion Similarity test A <- B - D

save.image('Prosopis_analysis.RData')

#plot Summary

pdf(paste('Results/summary_overlap/',res.row[i],'summary_D','.pdf',sep=''), width = 8, height = 4)
layout(matrix(c(1,1,1,1,2,3,4,5), 2, 4, byrow = F))
ecospat.plot.niche.dyn (z1, z2,quant=0,title=paste("piura"))
plot.new(); text(0.5,0.5,paste("Niche overlap:","\n","D=",round(results[i,5],3),"\n","\n",
                               "Niche dynamics indices:","\n",
                               "Unfilling:",round(results[i,1],3),"\n",
                               "Stability:",round(results[i,2],3),"\n",
                               "Expansion:",round(results[i,3],3),"\n"))
plot.overlap.test.brunno(obs=results[i,5], sim=read.csv(paste('Results/similarity_test_A-B/',sep='',"piura",'_null_sim_A-B_D.csv'))[,2],
                         p=results[i,10],title='Similarity A->B',type='D')
plot.overlap.test.brunno(obs=results[i,5], sim=read.csv(paste('Results/identity_test/',sep='',"piura",'_null_ident_D.csv'))[,2],
                         p=results[i,9],title="Equivalency",type='D')
plot.overlap.test.brunno(obs=results[i,5], sim=read.csv(paste('Results/similarity_test_B-A/',sep='',"piura",'_null_sim_B-A_D.csv'))[,2],
                         p=results[i,11],title='Similarity B->A',type='D')
dev.off()

pdf(paste('Results/summary_overlap/',res.row[i],'summary_I','.pdf',sep=''), width = 8, height = 4)
layout(matrix(c(1,1,1,1,2,3,4,5), 2, 4, byrow = F))
ecospat.plot.niche.dyn (z1, z2,quant=0,title=paste("piura"))
plot.new(); text(0.5,0.5,paste("Niche overlap:","\n","I=",round(results[i,4],3),"\n","\n",
                               "Niche dynamics indices:","\n",
                               "Unfilling:",round(results[i,1],3),"\n",
                               "Stability:",round(results[i,2],3),"\n",
                               "Expansion:",round(results[i,3],3),"\n"))
plot.overlap.test.brunno(obs=results[i,4], sim=read.csv(paste('Results/similarity_test_A-B/',sep='',"piura",'_null_sim_A-B_D.csv'))[,2],
                         p=results[i,7],title='Similarity A->B',type='I')
plot.overlap.test.brunno(obs=results[i,4], sim=read.csv(paste('Results/identity_test/',sep='',"piura",'_null_ident_D.csv'))[,2],
                         p=results[i,6],title="Equivalency",type='I')
plot.overlap.test.brunno(obs=results[i,4], sim=read.csv(paste('Results/similarity_test_B-A/',sep='',"piura",'_null_sim_B-A_D.csv'))[,2],
                         p=results[i,8],title='Similarity B->A',type='I')
dev.off()

}

save.image('Prosopis_analysis.RData')

#############  2) Comparison between species in native and invaded ranges ################
##########################################################################################
# 5) P.juliflora vs P. pallida within their native range
# 6) P.juliflora vs P. pallida within their invaded range
# 7) P.juliflora raw vs P. juliflora within their invaded range

spp1 <- list(jul_nat=jul_nat, #5
             jul_inv=jul_inv, #6
             jul_inv=jul_inv) #7
             
spp2 <- list(pal_nat=pal_nat, #5
             pal_inv=pal_inv, #6
             jul_inv_raw=jul_inv_raw) #7

clims1 <- list(sampBgNatJul, #5
               sampBgInv, #6
               sampBgInv) #7
               
clims2 <- list(sampBgNatPal, #5
               sampBgInv, #6
               sampBgInv) #7
               
### Loop ####

for (i in 1:length(spp1)) { # start from n=5
  
  j = i+4
  
  cat('\r','round',i,'of',length(spp1),'-',res.row[j],'\n')

  sp1 <- spp1[[i]]
  sp2 <- spp2[[i]]
  clim1<-clims1[[i]]
  clim2<-clims2[[i]]
  clim12<-rbind(clim1,clim2)
  
  # sample environmental values for species occurences
  occ.sp1<-na.exclude(ecospat.sample.envar(dfsp=sp1,colspxy=1:2,colspkept=1:3,dfvar=clim1,
                                           colvarxy=1:2,colvar="all",resolution=resolution))
  occ.sp2<-na.exclude(ecospat.sample.envar(dfsp=sp2,colspxy=1:2,colspkept=1:3,dfvar=clim2,
                                           colvarxy=1:2,colvar="all",resolution=resolution))
  
  # selection of variables to include in the analyses
  Xvar<-names(BgAll)
  nvar<-length(Xvar)
  
  #################################################################################################
  #################################### PCA-ENV ####################################################
  #################################################################################################
  
  row.w.1.env<-1-(nrow(clim1)/nrow(clim12))  # prevalence of clim1
  row.w.2.env<-1-(nrow(clim2)/nrow(clim12))  # prevalence of clim2
  row.w.env<-c(rep(row.w.1.env, nrow(clim1)),rep(row.w.2.env, nrow(clim2)),rep(0, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))
  
  # global dataset for the analysis and rows for each sub dataset
  data.env.occ<-rbind(clim1[Xvar],clim2[Xvar],occ.sp1[Xvar],occ.sp2[Xvar])
  row.clim1<-1:nrow(clim1)
  row.clim2<-(nrow(clim1)+1):(nrow(clim1)+nrow(clim2))
  row.clim12<-1:(nrow(clim1)+nrow(clim2))
  row.sp1<-(nrow(clim1)+nrow(clim2)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1))
  row.sp2<-(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+nrow(occ.sp2))
  
  # measures niche overlap along the two first axes of a PCA calibrated on all the pixels of the study areas
  pca.cal <-dudi.pca(data.env.occ,row.w = row.w.env, center = T, scale = T, scannf = F, nf = 2)
  
  # measure % explained by the first two PCs
  pca.scores <- 100 * pca.cal$eig/sum(pca.cal$eig)
  pca.scores <- pca.scores[1] + pca.scores[2]
  
  results$PCA[j] <- round(pca.scores, 2)
  
  # predict the scores on the axes
  scores.clim12<- pca.cal$li[row.clim12,]
  scores.clim1<- pca.cal$li[row.clim1,]
  scores.clim2<- pca.cal$li[row.clim2,]
  scores.sp1<- pca.cal$li[row.sp1,]
  scores.sp2<- pca.cal$li[row.sp2,]
  
  # calculation of occurence density and test of niche equivalency and similarity 
  z1<- ecospat.grid.clim.dyn(scores.clim12,scores.clim1,scores.sp1, th.sp=0, R)
  z2<- ecospat.grid.clim.dyn(scores.clim12,scores.clim2,scores.sp2, th.sp=0, R)
  
# Plot niche categories and species density created by ecospat.grid.clim.dyn
pdf(paste('Results/PCA_graphs/',res.row[j],'_PCA','.pdf',sep=''), width = 6, height = 5.5)
plot.PCA(pca.cal$co,pca.cal$eig)
dev.off()

# Calculate niche expansion, st <- stability and unfilling
dyn.indexes <- ecospat.niche.dyn.index (z1, z2, intersection=NA)$dynamic.index.w
results[j,1] <- dyn.indexes[[3]]
results[j,2] <- dyn.indexes[[2]]
results[j,3] <- dyn.indexes[[1]]

# Plot niche categories and species density created by ecospat.grid.clim.dyn
pdf(paste('Results/overlap_graphs/',res.row[j],'.pdf',sep=''), width = 7*.7, height = 5)
ecospat.plot.niche.dyn (z1, z2,interest=1,quant=.25,title=paste(res.row[j]))
legend('topleft',legend=c(paste("Unfilling:",round(results[j,1],3)),
                           paste("Stability:",round(results[j,2],3)),
                           paste("Expansion:",round(results[j,3],3))),
       bty='n')
dev.off()


# calculation of niche overlap - I and D indexes
niche_over <- ecospat.niche.overlap (z1, z2, cor=F)

results$I[j] <- niche_over$I
results$D[j] <- niche_over$D

# test of niche equivalency and similarity according to Warren et al. 2008
cat('\r','running niche equivalency test','\n')
a<-ecospat.niche.equivalency.test(z1, z2, rep=100) 
cat('\r','running niche similarity test A->B','\n')
b<-ecospat.niche.similarity.test(z1, z2, rep=100)
cat('\r','running niche similarity test A<-B','\n')
b2<-ecospat.niche.similarity.test(z2, z1, rep=100)


results[j,6] <- a$p.I
results[j,9] <- a$p.D

results[j,7] <- b$p.I
results[j,8] <- b2$p.I

results[j,10] <- b$p.D
results[j,11] <- b2$p.D


write.csv(a$sim$I,paste(sep='','Results/identity_test/',res.row[j],'_null_ident_I.csv')) # null distribuion Identity test - I
write.csv(a$sim$D,paste(sep='','Results/identity_test/',res.row[j],'_null_ident_D.csv')) # null distribuion Identity test - D

write.csv(b$sim$I,paste(sep='','Results/similarity_test_A-B/',res.row[j],'_null_sim_A-B_I.csv')) # null distribuion Similarity test A -> B - I
write.csv(b$sim$D,paste(sep='','Results/similarity_test_A-B/',res.row[j],'_null_sim_A-B_D.csv')) # null distribuion Similarity test A -> B - D

write.csv(b2$sim$I,paste(sep='','Results/similarity_test_B-A/',res.row[j],'_null_sim_B-A_I.csv')) # null distribuion Similarity test A <- B - I
write.csv(b2$sim$D,paste(sep='','Results/similarity_test_B-A/',res.row[j],'_null_sim_B-A_D.csv')) # null distribuion Similarity test A <- B - D

save.image('Prosopis_analysis.RData')

#plot Summary

pdf(paste('Results/summary_overlap/',res.row[j],'summary_D','.pdf',sep=''), width = 8, height = 4)
layout(matrix(c(1,1,1,1,2,3,4,5), 2, 4, byrow = F))
ecospat.plot.niche.dyn (z1, z2,quant=0,title=paste(res.row[j]))
plot.new(); text(0.5,0.5,paste("Niche overlap:","\n","D=",round(results[j,5],3),"\n","\n",
                               "Niche dynamics indices:","\n",
                               "Unfilling:",round(results[j,1],3),"\n",
                               "Stability:",round(results[j,2],3),"\n",
                               "Expansion:",round(results[j,3],3),"\n"))
plot.overlap.test.brunno(obs=results[j,5], sim=read.csv(paste('Results/similarity_test_A-B/',sep='',res.row[j],'_null_sim_A-B_D.csv'))[,2],
                         p=results[j,10],title='Similarity A->B',type='D')
plot.overlap.test.brunno(obs=results[j,5], sim=read.csv(paste('Results/identity_test/',sep='',res.row[j],'_null_ident_D.csv'))[,2],
                         p=results[j,9],title="Equivalency",type='D')
plot.overlap.test.brunno(obs=results[j,5], sim=read.csv(paste('Results/similarity_test_B-A/',sep='',res.row[j],'_null_sim_B-A_D.csv'))[,2],
                         p=results[j,11],title='Similarity B->A',type='D')
dev.off()

pdf(paste('Results/summary_overlap/',res.row[j],'summary_I','.pdf',sep=''), width = 8, height = 4)
layout(matrix(c(1,1,1,1,2,3,4,5), 2, 4, byrow = F))
ecospat.plot.niche.dyn (z1, z2,quant=0,title=paste(res.row[j]))
plot.new(); text(0.5,0.5,paste("Niche overlap:","\n","I=",round(results[j,4],3),"\n","\n",
                               "Niche dynamics indices:","\n",
                               "Unfilling:",round(results[j,1],3),"\n",
                               "Stability:",round(results[j,2],3),"\n",
                               "Expansion:",round(results[j,3],3),"\n"))
plot.overlap.test.brunno(obs=results[j,4], sim=read.csv(paste('Results/similarity_test_A-B/',sep='',res.row[j],'_null_sim_A-B_D.csv'))[,2],
                         p=results[j,7],title='Similarity A->B',type='I')
plot.overlap.test.brunno(obs=results[j,4], sim=read.csv(paste('Results/identity_test/',sep='',res.row[j],'_null_ident_D.csv'))[,2],
                         p=results[j,6],title="Equivalency",type='I')
plot.overlap.test.brunno(obs=results[j,4], sim=read.csv(paste('Results/similarity_test_B-A/',sep='',res.row[j],'_null_sim_B-A_D.csv'))[,2],
                         p=results[j,8],title='Similarity B->A',type='I')
dev.off()

}

save.image('Prosopis_analysis.RData')

write.csv(results,'results_COUE.csv')



######################### 3) Supp Analysis: Comparison Piura vs Invasive range  #################################

dir.create('Results/Piura_comparizons')

results_piura <- as.data.frame(matrix(ncol=length(res.col),nrow=1))
colnames(results_piura) <- res.col
rownames(results_piura) <- 'Piura'

spp1 <- list(piura_square = cbind(read.csv("Data/prosopis_occr_piura_square.csv")[,1:2],sp="Piura_occr"),
             piura_desert=cbind(read.csv("Data/prosopis_occr_piura_desert.csv")[,1:2],sp="Piura_occr"),
             piura_MCP=sampPiura[,1:2]) 


# see points
AllFrame<-readShapePoly("Data/All_Rectangle_world") # All Bg
proj4string(AllFrame)


Piura.sq<- readShapePoly("Data/Piura_square") # Piura square

Piura.shp <- readShapePoly("Data/Piura") # Piura State
crs(Piura.shp) <-"+init=epsg:32718"
Piura.shp <- spTransform(Piura.shp, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

plot(Piura.sq, axes=T)
plot(Piura.shp, add=T)
points(spp1[[1]][,-3],cex=.5, col='red')
points(spp1[[2]][,-3],cex=.5, col='black')


### Vamos marcar as ocorrencias que estão em piura.
piura_occr=read.csv("Data/prosopis_occr_piura_desert.csv")[,1:2]
#piura_MCP=cbind(sampPiura[,1:2],sp='Piura_occr')

spp <- sp.prosopis
spp$sp <- as.character(spp$sp)

for(i in 1:nrow(spp)){
  if(any(spp$x[i]==piura_occr$X)){
    if(any(spp$y[i]==piura_occr$Y)){
      spp$sp[i] <- "piura"
    }
  }
}

unique(spp$sp)
table(spp$sp)

row.piura <- which(spp$sp=='piura') # estas são as linhas de spp que estão em piura desert

# Agora vamos fazer o mesmo para a região de piura

# 1) piura desert
bg.nat <- cbind(sampBgNat[,1:2],sp='bg')
bg.nat$sp <- as.character(bg.nat$sp)

for(i in 1:nrow(sampBgNat)){
  cat(paste(round(i/nrow(sampBgNat),2)*100,'%'),'\r')
  if(any(round(bg.nat$x[i],2)==round(sampBgPiura.d$x,2))){
    if(any(round(bg.nat$y[i],2)==round(sampBgPiura.d$y,2))){
      bg.nat$sp[i] <- "piura"
    }
  }
}

unique(bg.nat$sp)
table(bg.nat$sp)

row.bg.pd <- which(bg.nat$sp=='piura') # estas são as linhas de piura desert que estão em bgnat

# 2) piura desert + codilera
bg.nat <- cbind(sampBgNat[,1:2],sp='bg')
bg.nat$sp <- as.character(bg.nat$sp)

for(i in 1:nrow(sampBgNat)){
  cat(paste(round(i/nrow(sampBgNat),2)*100,'%'),'\r')
  if(any(round(bg.nat$x[i],2)==round(sampBgPiura.dc$x,2))){
    if(any(round(bg.nat$y[i],2)==round(sampBgPiura.dc$y,2))){
      bg.nat$sp[i] <- "piura"
    }
  }
}

unique(bg.nat$sp)
table(bg.nat$sp)

row.bg.pdc <- which(bg.nat$sp=='piura') # estas são as linhas de piura desert + codilera que estão em bgnat

### fazer a analise

clim1<-sampBgNat
clim2<-sampBgInv
clim12<-rbind(clim1,clim2)

# sample environmental values for species occurences
occ.sp1<-na.exclude(ecospat.sample.envar(dfsp=spp,colspxy=1:2,colspkept=1:3,dfvar=clim1,
                                         colvarxy=1:2,colvar="all",resolution=resolution))
occ.sp2<-na.exclude(ecospat.sample.envar(dfsp=spp,colspxy=1:2,colspkept=1:3,dfvar=clim2,
                                         colvarxy=1:2,colvar="all",resolution=resolution))

# selection of variables to include in the analyses
Xvar<-names(BgAll)
nvar<-length(Xvar)

#################################################################################################
#################################### PCA-ENV ####################################################
#################################################################################################

row.w.1.env<-1-(nrow(clim1)/nrow(clim12))  # prevalence of clim1
row.w.2.env<-1-(nrow(clim2)/nrow(clim12))  # prevalence of clim2
row.w.env<-c(rep(row.w.1.env, nrow(clim1)),rep(row.w.2.env, nrow(clim2)),
             rep(0, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))

# global dataset for the analysis and rows for each sub dataset
data.env.occ<-rbind(clim1[Xvar],clim2[Xvar],occ.sp1[Xvar],occ.sp2[Xvar])
row.clim1<-1:nrow(clim1)
row.clim2<-(nrow(clim1)+1):(nrow(clim1)+nrow(clim2))
row.clim12<-1:(nrow(clim1)+nrow(clim2))
row.sp1<-(nrow(clim1)+nrow(clim2)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1))
row.sp2<-(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+nrow(occ.sp2))

row.piura2<-row.sp1[row.piura]
row.bg2<-row.clim1[row.bg.pd]

# measures niche overlap along the two first axes of a PCA calibrated on all the pixels of the study areas
pca.cal <-dudi.pca(data.env.occ,row.w = row.w.env, center = T, scale = T, scannf = F, nf = 2)

# predict the scores on the axes
scores.clim12<- pca.cal$li[row.clim12,]
scores.clim1<- pca.cal$li[row.clim1,]
scores.clim2<- pca.cal$li[row.clim2,]
scores.sp1<- pca.cal$li[row.sp1,]
scores.sp2<- pca.cal$li[row.sp2,]
scores.piura <- pca.cal$li[row.piura2,]
scores.bg <- pca.cal$li[row.bg2,]

# calculation of occurence density and test of niche equivalency and similarity 
z1<- ecospat.grid.clim.dyn(scores.clim12,scores.clim1,scores.sp1, th.sp=0, R)
z2<- ecospat.grid.clim.dyn(scores.clim12,scores.clim2,scores.sp2, th.sp=0, R)


# Plot niche categories and species density created by ecospat.grid.clim.dyn
pdf(paste('Results/overlap_graphs/','piura_overlap.pdf',sep=''), width = 7*.7, height = 5)
ecospat.plot.niche.dyn (z1, z2,interest=2,quant=.25,title='')
points(scores.bg$Axis1,scores.bg$Axis2,col = 'lightgray')
#points(scores.piura$Axis1,scores.piura$Axis2,col = 'black')
dev.off()











######################### 4) Supp Analysis: comparing native range conditions: P. jul vs P. pal  #################################

occ.sp1<-na.exclude(ecospat.sample.envar(dfsp=jul_nat,colspxy=1:2,colspkept=1:3,dfvar=sampBgNatJul,
                                         colvarxy=1:2,colvar="all",resolution=resolution))[,4:9]
occ.sp2<-na.exclude(ecospat.sample.envar(dfsp=pal_nat,colspxy=1:2,colspkept=1:3,dfvar=sampBgNatPal,
                                         colvarxy=1:2,colvar="all",resolution=resolution))[,4:9]

for(i in 3:length(sampBgNatPal)){ # a partir de 3 por causa do XY
  #tmp <- rbind(cbind(env=sampBgNatJul[i],sp='jul'),cbind(env=sampBgNatPal[i],sp='pal'))
  #colnames(tmp) <- c('env','sp')
  #tmp$env <- scale(tmp$env)
  
  #supp_env[i-2] <- summary(aov(env~sp, data=tmp))
  #res.t <- pairwise.t.test(tmp$env,tmp$sp,p.adj = "bonf")
  #res.t <- t.test(sampBgNatJul[i], sampBgNatPal[i],var.equal=T)
  
  res.t <- t.test(occ.sp1[,i-2], occ.sp2[,i-2])
  p.t <- res.t$p.value
  #p.t <- p.adjust(res.t$p.value,"bonferroni",6)
  
  supp_env[i-2,'P-value'] <- signif(p.t,3)
  supp_env[i-2,'t'] <- round(res.t$statistic,2)
  supp_env[i-2,'P.juliflora'] <- paste(round(mean(sampBgNatJul[,i]),2),' (',round(sd(sampBgNatJul[,i]),2),')',sep='')
  supp_env[i-2,'P.pallida'] <- paste(round(mean(sampBgNatPal[,i]),),' (',round(sd(sampBgNatPal[,i]),2),')',sep='')
  supp_env[i-2,'95.CI'] <- paste(round(res.t$conf.int[1],2),round(res.t$conf.int[2],2),sep='-')
}
supp_env <- data.frame(supp_env)
supp_env

write.csv(supp_env,'table S1.csv')

pdf(paste('Supp_image_diff_env.pdf',sep=''), width = 10, height = 6,useDingbats=F)
par(mfrow=c(2,3),mar=c(5,3,5,1))
for(i in 3:length(sampBgNatPal)){ # a partir de 3 por causa do XY
  tmp <- rbind(cbind(env=sampBgNatJul[i],sp='jul'),cbind(env=sampBgNatPal[i],sp='pal'))
  colnames(tmp) <- c('env','sp')
  tmp$env <- scale(tmp$env)
  
  boxplot(env~sp, data=tmp, main = 
            paste(names(sampBgNat)[i],'\n','t = ',supp_env$t[i-2], '\n','p = ', supp_env$p[i-2]), 
          names=c("",""))
  mtext(c('P. juliflora', 'P. pallida'),at=1:2,line=c(1.5,1.5),side=1, font=3)
  }
dev.off()

##################################################################################################
## Supp analysis time lag - human effect constraints
library(rgdal)

# load species distributions
julnatpred <- raster('Results/ENMeval/jul_nat.asc')
palnatpred <- raster('Results/ENMeval/pal_nat.asc')
#plot(julnatpred)

# load cattle density
x <- new("GDALReadOnlyDataset", '/media/brunno/FAT/GIS/cattle/cattle/totcor/glbctd1t0503m/hdr.adf')
getDriver(x)
getDriverLongName(getDriver(x))
xx<-asSGDF_GROD(x,output.dim=c(200, 200))
r <- raster(xx)
plot(r)

# load goat density
x <- new("GDALReadOnlyDataset", '/media/brunno/FAT/GIS/goats/goats/totcor/glbgtd1t0503m/hdr.adf')
getDriver(x)
getDriverLongName(getDriver(x))
xx<-asSGDF_GROD(x,output.dim=c(200, 200))
goats <- raster(xx)
goats <- crop(goats, extent(BgInv))
plot(goats)

# load global cattle density from Robinson et al. 2014 Plos
gcat <- raster('/media/brunno/FAT/GIS/CATTLE (2)/Glb_Cattle_CC2006_AD.tif')
gcat <- crop(gcat, extent(BgInv))
plot(gcat,col=colorRampPalette(c("#3d46a9","#49bff0","#58c53f","#f2f327", "#ec2236"))(40))

# load global ruminant data from  Robinson et al. 2014 Plos
grum <- raster('/media/brunno/FAT/GIS/Ruminat_prod_systemsv5/Ruminant_production_systems.tif')
grum <- crop(grum, extent(BgInv))
plot(grum)

# extract suitability values in points
jul_r2 <- extract(julnatpred, jul_inv[,1:2])

# extract cattle/goat values in points
cat_r <- extract(r, jul_inv[,1:2])
gcat_r <- extract(gcat, jul_inv[,1:2])
grum_r <- extract(grum, jul_inv[,1:2])
goat_r <- extract(goats, jul_inv[,1:2])

plot(gcat_r,jul_r2)
plot(grum_r,jul_r2)
plot(log(goat_r+1),jul_r2); summary(lm(jul_r2~log(goat_r+1)));abline(lm(jul_r2~log(goat_r+1)))


# extract suitability values in points
pal_r2 <- extract(palnatpred, pal_inv[,1:2])

# extract cattle values in points
cat_r <- extract(r, pal_inv[,1:2])
gcat_r <- extract(gcat, pal_inv[,1:2])
goat_r <- extract(goats, pal_inv[,1:2])


plot(gcat_r,pal_r2)
plot(log(goat_r+1),pal_r2); summary(lm(pal_r2~log(goat_r+1)));abline(lm(pal_r2~log(goat_r+1)))


