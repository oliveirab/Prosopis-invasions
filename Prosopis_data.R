
################################################################
#### Prosopis invasions on South American dry lands
#### Brunno Oliveira, 2015
#### Universidade Federal do Rio Grande do Norte - Brasil
################################################################

#### 1) LOAD PREDICTOR VARIABLES
#### 2) CROP PREDICTOR VARIABLES
#### 3) FIX OCCURRENCES (Split native and invasive occurrences for Prosopis complex, P julifoa and P pallida)
#### 4) ANALYSE OUTPUT of raw vs reanalyzed vouchers
#### 5) CREATE FIGURE OF SPECIES DISTRIBUTIONS
#### 6) SAVE RESULTS

#set working directory
setwd("F:/Prosopis")

gc()
rm(list=ls())

##########Loading packages################
library(raster)
library(maptools)
library(ecospat)
library(rgeos)
library(rgdal)

source("vif_func.R")
source("functions.R")

### READ ###
#BIO1 = Annual Mean Temperature
#BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
#BIO3 = Isothermality (#BIO2/#BIO7) (* 100)
#BIO4 = Temperature Seasonality (standard deviation *100)
#BIO5 = Max Temperature of Warmest Month
#BIO6 = Min Temperature of Coldest Month
#BIO7 = Temperature Annual Range (#BIO5-#BIO6)
#BIO8 = Mean Temperature of Wettest Quarter
#BIO9 = Mean Temperature of Driest Quarter
#BIO10 = Mean Temperature of Warmest Quarter
#BIO11 = Mean Temperature of Coldest Quarter
#BIO12 = Annual Precipitation
#BIO13 = Precipitation of Wettest Month
#BIO14 = Precipitation of Driest Month
#BIO15 = Precipitation Seasonality (Coefficient of Variation)
#BIO16 = Precipitation of Wettest Quarter
#BIO17 = Precipitation of Driest Quarter
#BIO18 = Precipitation of Warmest Quarter
#BIO19 = Precipitation of Coldest Quarter


# Load shape files ###### 
WE<-readShapePoly("Data/WE_countries_crop") # WE
WE_c<-readShapePoly("Data/WE_contour_crop") # WE contour
pal.nat.bg.shp <- readShapePoly("Data/pallida_native_bg_eco") # P. pallida Bg
jul.nat.bg.shp <- readShapePoly("Data/juliflora_native_bg_eco") # P. juliflora Bg
NatFrame<-gUnion(jul.nat.bg.shp,pal.nat.bg.shp) # Native Bg
InvFrame<-readShapePoly("Data/Brazil") # Invasive Bg
AllFrame<-readShapePoly("Data/All_Rectangle_world") # All Bg
World <- readShapePoly("Data/Wrld_shp") # World shp

piura.bg.shp <- readShapePoly("Data/Piura") # Piura Bg
crs(piura.bg.shp) <-"+init=epsg:32718"
piura.bg.shp <- spTransform(piura.bg.shp, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

pdesert.bg.shp <- readShapePoly("Data/piura_desert") # Piura desert Bg

pdesertcord.bg.shp <- readShapePoly("Data/piura_tumbes+cordileira") # Piura desert + cordilera Bg

###### Get predictors ######
climatelayers <- getData('worldclim', var='bio', res=2.5) # 2.5'min' ~ 5km ~ 0.041666 Dg

# selection of variables to include in the analyses
# We used temperature and precipitation variables from the WorldClim database (Hijmans et al. 2005) 
# representing average and extreme climatic conditions, as well as measures of temporal climatic variability, 
# to characterize the native and invaded realized niches.

# We used the variance inflation factor (VIF) to detect collinearity (Marquardt 1970). 
# A VIF greater than 10 is a signal that the model has a collinearity problem (Chatterjee and Hadi 2006). 
# We excluded the variables with large VIF values (greater than 10) one by one using a stepwise
# procedure. We repeated this procedure until all strongly correlated variables (i.e. with VIF ? 10) were excluded

#Extract random 10000 random values from variables
samp<-data.frame(sampleRandom(climatelayers, 10000, na.rm=T))

Xvar <- vif_func(samp, thresh=3, trace=F)

climatelayers <- climatelayers[[Xvar]]

BgAll <- crop(climatelayers, extent(AllFrame)); BgAll <- mask(BgAll, AllFrame)

BgNat <- crop(climatelayers, extent(NatFrame)); BgNat <- mask(BgNat, NatFrame)
BgNatJul <- crop(climatelayers, extent(jul.nat.bg.shp)); BgNatJul <- mask(BgNatJul, jul.nat.bg.shp); BgNatJul <- crop(BgNatJul, extent(BgAll))
BgNatPal <- crop(climatelayers, extent(pal.nat.bg.shp)); BgNatPal <- mask(BgNatPal, pal.nat.bg.shp); BgNatPal <- crop(BgNatPal, extent(BgAll))
BgInv <- crop(climatelayers, extent(InvFrame)); BgInv <- mask(BgInv, InvFrame)
BgPiura  <- crop(climatelayers, extent(piura.bg.shp)); BgPiura <- mask(BgPiura, piura.bg.shp); BgPiura <- crop(BgPiura, extent(BgAll))
BgPiura.d  <- crop(climatelayers, extent(pdesert.bg.shp)); BgPiura.d <- mask(BgPiura.d, pdesert.bg.shp); BgPiura.d <- crop(BgPiura.d, extent(BgAll))
BgPiura.dc  <- crop(climatelayers, extent(pdesertcord.bg.shp)); BgPiura.dc <- mask(BgPiura.dc, pdesertcord.bg.shp); BgPiura.dc <- crop(BgPiura.dc, extent(BgAll))


rm(climatelayers)

resolution <- res(BgAll[[1]])[1] # ~ 0.04166Dg ~ 5km

# Check if is everything ok...
#plot(BgNat[[1]])
#plot(BgInv[[1]])
#plot(BgAll[[1]])


# Load occurrence files

#install.packages("stringr", dependencies=TRUE)
require(stringr)

juliflora <- read.csv('Data/occurrence/P_juliflora_occr.csv',encoding='UTF-8') # all occurrences
juliflora$y <- as.numeric(str_trim(juliflora$y)) # remove white spaces from vector
pallida <- read.csv('Data/occurrence/P_pallida_occr.csv') # all occurrences
sp.prosopis <- rbind(juliflora,pallida) # all occurrences

sp.prosopis <- SpatialPoints(sp.prosopis[,c('x','y')])

teste <- extract(BgAll[[1]],sp.prosopis) # clean data : remove points that fall in the ocean
sp.prosopis <- sp.prosopis[!is.na(teste)]
sp.prosopis <- as.data.frame(sp.prosopis)
sp.prosopis.raw <- data.frame(sp.prosopis,sp='Complex') # all native and invaded records

# Prosopis native occurrences
pointx <- SpatialPoints(sp.prosopis[,c('x','y')])

teste <- extract(BgNat[[1]],pointx) # clean data : remove points that fall outside the invasive range
sp.nat <- as.data.frame(pointx[!is.na(teste)])
sp.nat <- data.frame(sp.nat,sp='Complex_native') # all native records

# Prosopis invasive occurrences
pointx <- SpatialPoints(sp.prosopis[,c('x','y')])

teste <- extract(BgInv[[1]],pointx) # clean data : remove points that fall outside the invasive range
sp.inv <- as.data.frame(pointx[!is.na(teste)])
sp.inv.raw <- data.frame(sp.inv,sp='Complex_invasive') # sp invasive raw

# P. juliflora native occurrences
juliflora <- SpatialPoints(juliflora[,c('x','y')])

teste <- extract(BgNatJul[[1]],juliflora) # clean data 
jul_nat <- as.data.frame(juliflora[!is.na(teste)])
jul_nat <- data.frame(jul_nat,sp='juliflora_native')

# P. pallida native occurrences
pallida <- SpatialPoints(pallida[,c('x','y')])

teste <- extract(BgNatPal[[1]],pallida) # clean data 
pal_nat <- as.data.frame(pallida[!is.na(teste)])
pal_nat <- data.frame(pal_nat,sp='pallida_native')

### Invasive occurrences
## reaxamined
invasive.occ <- read.csv('Data/occurrence/Prosopis_reexamined_voushers.csv')

jul_inv <- subset(invasive.occ,Species=='Prosopis juliflora')
jul_inv <- SpatialPoints(na.omit(jul_inv)[,c('Long','Lat')])
teste <- jul_inv[complete.cases(over(jul_inv,InvFrame)),] # clean data 
jul_inv <- as.data.frame(teste)
jul_inv <- data.frame(jul_inv,sp='juliflora_invasive')
names(jul_inv) <- names(jul_nat)

pal_inv <- subset(invasive.occ,Species=='Prosopis pallida')
pal_inv <- SpatialPoints(na.omit(pal_inv)[,c('Long','Lat')])
teste <- pal_inv[complete.cases(over(pal_inv,InvFrame)),] # clean data 
pal_inv <- as.data.frame(teste)
pal_inv <- data.frame(pal_inv,sp='pallida_invasive')
names(pal_inv) <- names(pal_nat)

# Complex Invasive (reanalyzed voucher)
sp.inv <- rbind(jul_inv,pal_inv)

# Complex prosopis (with reanalyzed vouchers)
sp.prosopis <- rbind(sp.nat,sp.inv)

## Raw data invasive
# P. juliflora invasive raw
teste <- subset(invasive.occ,Original.determination=='Prosopis juliflora')
teste <- SpatialPoints(na.omit(teste)[,c('Long','Lat')])
teste <- teste[complete.cases(over(teste,InvFrame)),] # clean data 
jul_inv_raw <- as.data.frame(teste)
jul_inv_raw <- data.frame(jul_inv_raw,sp='juliflora_invasive')
names(jul_inv_raw) <- names(jul_inv)

# P. pallida invasive raw (DONT WORK BECAUSE THERE ARE FEW OCCURRENCES)
teste <- subset(invasive.occ,Original.determination=='Prosopis pallida')
teste <- SpatialPoints(na.omit(teste)[,c('Long','Lat')])
teste <- teste[complete.cases(over(teste,InvFrame)),] # clean data 
pal_inv_raw <- as.data.frame(teste)
pal_inv_raw <- data.frame(pal_inv_raw,sp='pallida_invasive')
names(pal_inv_raw) <- names(pal_inv)

# Complex Invasive (raw)
sp.inv.raw <- rbind(jul_inv_raw,pal_inv_raw)

# Complex prosopis (raw)
sp.prosopis.raw <- rbind(sp.nat,sp.inv.raw)

### Summary occr before desaggregation
spps <- c('P. juliflora', 'P. pallida')
regions <- c('Native', 'Invasive Brazil raw','Invasive Brazil')

summary.occr <- as.data.frame(matrix(ncol=length(spps),nrow=length(regions)))
colnames(summary.occr) <- spps
rownames(summary.occr) <- regions
summary.occr[1,1] <- nrow(jul_nat)
summary.occr[2,1] <- nrow(jul_inv_raw)
summary.occr[3,1] <- nrow(jul_inv)
summary.occr[1,2] <- nrow(pal_nat)
summary.occr[2,2] <- nrow(pal_inv_raw)
summary.occr[3,2] <- nrow(pal_inv)
print(summary.occr)

write.csv(summary.occr,'summary_occr_before.csv')

# clean data : # Remove occurrences closer than a minimum distance to each other (remove aggregation). 
sp.prosopis <- ecospat.occ.desaggregation(df=sp.prosopis,colxy=1:2, 
                                          colvar=3, min.dist=resolution,plot=F)  
jul_nat <- ecospat.occ.desaggregation(df=jul_nat,colxy=1:2, 
                                      colvar=3, min.dist=resolution,plot=F)  
pal_nat <- ecospat.occ.desaggregation(df=pal_nat,colxy=1:2, 
                                      colvar=3, min.dist=resolution,plot=F)  
jul_inv_raw <- ecospat.occ.desaggregation(df=jul_inv_raw,colxy=1:2, 
                                          colvar=3, min.dist=resolution,plot=F)  
pal_inv_raw <- ecospat.occ.desaggregation(df=pal_inv_raw,colxy=1:2, 
                                          colvar=3, min.dist=resolution,plot=F)  
jul_inv <- ecospat.occ.desaggregation(df=jul_inv,colxy=1:2, 
                                      colvar=3, min.dist=resolution,plot=F)  
pal_inv <- ecospat.occ.desaggregation(df=pal_inv,colxy=1:2, 
                                      colvar=3, min.dist=resolution,plot=F)  

### Summary occr after desaggregation
spps <- c('P. juliflora', 'P. pallida')
regions <- c('Native', 'Invasive raw', 'Invasive Brazil')

summary.occr <- as.data.frame(matrix(ncol=length(spps),nrow=length(regions)))
colnames(summary.occr) <- spps
rownames(summary.occr) <- regions
summary.occr[1,1] <- nrow(jul_nat)
summary.occr[2,1] <- nrow(jul_inv_raw)
summary.occr[3,1] <- nrow(jul_inv)
summary.occr[1,2] <- nrow(pal_nat)
summary.occr[2,2] <- nrow(pal_inv_raw)
summary.occr[3,2] <- nrow(pal_inv)
print(summary.occr)

write.csv(summary.occr,'summary_occr_agrregated.csv')

sp.juliflora.raw <- rbind(jul_nat,jul_inv_raw)
sp.pallida.raw <- rbind(pal_nat,pal_inv_raw)
sp.juliflora <- rbind(jul_nat,jul_inv)
sp.pallida <- rbind(pal_nat,pal_inv)
sp.prosopis <- rbind(sp.nat,sp.inv)

# Introdution sites XY

intro_sites <- rbind(data.frame(readShapePoints("Data/Serra_Talhada_point.shp"))[,1:3],
                     data.frame(readShapePoints("Data/Angicos_point.shp"))[,1:3])

colnames(intro_sites) <- c('place','x','y')
intro_sites$place <- c("Serra_Talhada","Angicos")

# Piura XY

piura <- data.frame(readShapePoints("Data/Piura_XY.shp"))[,1:3]
colnames(piura) <- c('place','x','y')
piura$place <- "Piura"


# get values from BgAll
sampBgAll<-na.exclude(as.data.frame(BgAll,xy=T))
# get values from BgNat
sampBgNat<-na.exclude(as.data.frame(BgNat,xy=T))
# get values from BgJulNat
sampBgNatJul<-na.exclude(as.data.frame(BgNatJul,xy=T))
# get values from BgJulNat
sampBgNatPal<-na.exclude(as.data.frame(BgNatPal,xy=T))
# get values from BgInv
sampBgInv<-na.exclude(as.data.frame(BgInv,xy=T))
# get values from BgPiura.desert
sampBgPiura<-na.exclude(as.data.frame(BgPiura,xy=T))
# get values from BgPiura.desert
sampBgPiura.d<-na.exclude(as.data.frame(BgPiura.d,xy=T))
# get values from BgPiura.desert+cordilera
sampBgPiura.dc<-na.exclude(as.data.frame(BgPiura.dc,xy=T))


# get values from occurrences a MCP in Piura region
occ_MCP_Piura <- readShapePoly("Data/occurrences_in_piura_MCP.shp")

sampPiura  <- crop(BgAll, extent(occ_MCP_Piura)) ; sampPiura <- mask(sampPiura, occ_MCP_Piura)

sampPiura <- na.exclude(as.data.frame(sampPiura,xy=T))



######################
### Save image
save.image('Prosopis_data.RData')

#################################################################################################
########################################## Figure 1 #############################################

#pch = symbols
#lwd = line width
#cex = symbol size
#col = border color
#bg = fill color

library(plyr)


pdf('Figure1.pdf',width=10.5, height=9.5)
# Plot the First map
par(las=2, ps=12)
plot(AllFrame, lwd=2, axes=F);plot(WE, lwd=1, add=T);plot(WE_c, lwd=2, add=T)
points(jul_nat[,-3], pch=21, lwd=1, cex=1.3, col='white', bg='black')
points(pal_nat[,-3], pch=21, lwd=1, cex=1, col='black', bg='white')
points(jul_inv[,-3], pch=24, lwd=1, cex=1.1, col='white', bg='black')
points(pal_inv[,-3], pch=24, lwd=1, cex=.8, col='black', bg='white')

#Add Scale bar
scalebar_b(xc=-42, yc=13, len=1000/111.32, units="Km", ndivs=1, subdiv=1000, t.cex = .8)

#Add North
north.arrow(xb=-42, yb=22, len=.7, lab="N",cex.lab=.8)

# Plot the second map as inset 
#   x_left, x_right, y_bottom, y_top are portions of the main chart X and Y spans.
par(fig=c(0.05, 0.45, 0.05,0.525), new=T, las=1, ps=9)
plot(WE_c);#plot(NatFrame,add=T,lty=2);plot(InvFrame,add=T,lt=2)
points(piura[,-1], pch=23, lwd=1.5, cex=1.6, col='white', bg='black')
points(intro_sites[,-1], pch=23, lwd=1.5, cex=1.3, col='black', bg='white')


# Return par(fig=) to full size for next plot
par(fig=c(0,1,0,1))
# Save figure
dev.off()
