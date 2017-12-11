#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°#
# Harvey and MacDougall  
# Non-interacting impacts of fertilization and habitat area on plant diversity via contrasting assembly mechanisms.#
#                                           Script                                                #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°#


#########################################################################
################# DATA STRUCTURE AND VISUALIZATION
#########################################################################

rm(list=ls())

#°°°°°°°°°°°°°°°°°°°°#
# Load libraries     #
#°°°°°°°°°°°°°°°°°°°°#

library(vegan)
library(scales)
library(reshape)
library(nlme)
library(Hmisc)
library(indicspecies)
library(adespatial)
library(sciplot)

#°°°°°°°°°°°°°°°°°°°°#
# Load data          #
#°°°°°°°°°°°°°°°°°°°°#

Dat <- "~/Documents/Research/1.Projects/PHD/Data/R_scripts_data/" #Working directory for scripts and data
fig.path <- "~/Documents/Research/1.Projects/PHD/redaction/Global.change.Plant/5.Diversity_Distribution/Figures/" #Working directory for scripts and data

Data = read.delim(paste0(Dat,"Plant(12-13-14)_20161107.txt"))
treatment = read.delim(paste0(Dat,"treatments.txt"))
space = read.delim(paste0(Dat,"space.txt"))
biomass = read.delim(paste0(Dat,"PlantBiomass(12-13-14).txt")) 
light = read.delim(paste0(Dat,"env.txt")) 


#...Remove Mainland data (not used for this study)
Data2 = Data[Data$size !="mainland",]

#...Fix the variable structure
{ 
#.Island as a discrete factor
island.order = c(1:36)
island.order = as.character(island.order)
Data2$Island = factor(as.character(Data2$Island),levels=island.order,labels=island.order)
#.Sampling year
years.order = c("2012","2013","2014")
Data2$Years = factor(as.character(Data2$Years),levels=years.order,labels=c("2012","2013","2014"))
#.Island size
size.order = c("25","100","400")
Data2$size = factor(as.character(Data2$size),levels=size.order,labels=c("25","100","400"))
#.Sampling quadrats, within each island (see Methods)
quadrat.order = c("A","B","C","D","E")
Data2$quadrat = factor(as.character(Data2$quadrat),levels=quadrat.order,labels=c("A","B","C","D","E"))
#.Sampling plots, withine ach island (only for 2013)
plot.order = c("1","2","3","4")
Data2$plot = factor(as.character(Data2$plot),levels=plot.order,labels=c("1","2","3","4"))
#.Nitrogen addition (change labelling)
levels(Data2$nitrogen)= c("F-","F+")
#.Mowing treatment (change labelling)
levels(Data2$mowing) = c("D-","D+")
#.Distance from mainland (change to character - technicality for the next step - will be reverted to numeric after)
Data2$dist.mainland = as.character(Data2$dist.mainland)
#.Free the R environment of all the transitory objects
rm(island.order,years.order,size.order,quadrat.order,plot.order)

}

#°°°°°°°°°°°°°°°°°°°°#
# 2 spatial scales   #
#°°°°°°°°°°°°°°°°°°°°#

{ 
#Separate 2012, 2013 and 2014 for convenience
Data.2012 = Data2[which(Data2$Years == 2012),]
Data.2013 = Data2[which(Data2$Years == 2013),]
Data.2014 = Data2[which(Data2$Years == 2014),]

#Average each 1m2 Plot within a Quadrat for 2012 and 2013
Data.2012.2 = cast(melt(Data.2012),Years + Island + quadrat + size + nitrogen + mowing + dist.mainland ~ variable, mean)
Data.2012.2 = data.frame(Data.2012.2)
Data.2013.2 = cast(melt(Data.2013),Years + Island + quadrat + size + nitrogen + mowing + dist.mainland ~ variable, mean)
Data.2013.2 = data.frame(Data.2013.2)

#Sum each 1 m2 Plot within a Quadrat for 2014
Data.2014.2 = cast(melt(Data.2014),Years + Island + quadrat + size + nitrogen + mowing + dist.mainland ~ variable, sum)
Data.2014.2 = data.frame(Data.2014.2)

####Bind them back together - Data3 will be the dataset used for Quadrat-level analyses
Data3 = rbind(Data.2012.2,Data.2013.2,Data.2014.2)

#####Pool patch-level information - Data4 will be the dataset used for Patch-level analyses
Data4 = cast(melt(Data3),Years + Island + size + nitrogen + mowing + dist.mainland ~ variable, mean)
Data4 = data.frame(Data4)

#Turn back variables as numerical 
Data4$dist.mainland = as.numeric(Data4$dist.mainland)
Data4$Years = as.numeric(Data4$Years)
Data3$dist.mainland = as.numeric(Data3$dist.mainland)
Data3$Years = as.numeric(Data3$Years)

}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°#
# SAC curves                 #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°#

{ 

#Because the sampling design was different across all years we here test how it changed across years and spatial scale

##Across years - Patch-level

pdf(paste0(fig.path,"SAC_PATCH(20170327.pdf"),width=10,height=5)

par(mfrow=c(1,2))
# #2012
# large.sample=specaccum(Data4[which(Data4$Years==1 & Data4$size==400),7:109],conditioned=FALSE,method="random")
# plot(large.sample,add=FALSE,ci=2,ci.type="bar",xvar=c("sites"),col="darkgreen",xlab="",ylab="Number species",ylim=c(0,60))
# medium.sample=specaccum(Data4[which(Data4$Years==1 & Data4$size==100),7:109],conditioned=FALSE,method="random")
# plot(medium.sample,add=TRUE,ci=2,ci.type="bar",xvar=c("sites"),col="darkred",ylim=c(0,60))
# small.sample=specaccum(Data4[which(Data4$Years==1 & Data4$size==25),7:109],conditioned=FALSE,method="random")
# plot(small.sample,add=TRUE,ci=2,ci.type="bar",xvar=c("sites"),col="darkblue",ylim=c(0,60))
#2013
large.sample=specaccum(Data4[which(Data4$Years==2 & Data4$size==400),7:109],conditioned=FALSE,method="random")
plot(large.sample,add=FALSE,ci=2,ci.type="bar",xvar=c("sites"),col="darkgreen",xlab="",ylab="Number species",ylim=c(0,60))
medium.sample=specaccum(Data4[which(Data4$Years==2 & Data4$size==100),7:109],conditioned=FALSE,method="random")
plot(medium.sample,add=TRUE,ci=2,ci.type="bar",xvar=c("sites"),col="darkred",ylim=c(0,60))
small.sample=specaccum(Data4[which(Data4$Years==2 & Data4$size==25),7:109],conditioned=FALSE,method="random")
plot(small.sample,add=TRUE,ci=2,ci.type="bar",xvar=c("sites"),col="darkblue",ylim=c(0,60))
#2014
large.sample=specaccum(Data4[which(Data4$Years==3 & Data4$size==400),7:109],conditioned=FALSE,method="random")
plot(large.sample,add=FALSE,ci=2,ci.type="bar",xvar=c("sites"),col="darkgreen",xlab="",ylab="Number species",ylim=c(0,60))
medium.sample=specaccum(Data4[which(Data4$Years==3 & Data4$size==100),7:109],conditioned=FALSE,method="random")
plot(medium.sample,add=TRUE,ci=2,ci.type="bar",xvar=c("sites"),col="darkred",ylim=c(0,60))
small.sample=specaccum(Data4[which(Data4$Years==3 & Data4$size==25),7:109],conditioned=FALSE,method="random")
plot(small.sample,add=TRUE,ci=2,ci.type="bar",xvar=c("sites"),col="darkblue",ylim=c(0,60))
dev.off()


##Across years - at Quadrat-level

pdf(paste0(fig.path,"SAC_QUADRAT(20170327.pdf"),width=10,height=5)
par(mfrow=c(1,2))
#2012
# large.sample=specaccum(Data3[which(Data3$Years==1 & Data3$size==400),7:109],conditioned=FALSE,method="random")
# plot(large.sample,add=FALSE,ci=2,ci.type="bar",xvar=c("sites"),col="darkgreen",xlab="",ylab="Number species")
# medium.sample=specaccum(Data3[which(Data3$Years==1 & Data3$size==100),7:109],conditioned=FALSE,method="random")
# plot(medium.sample,add=TRUE,ci=2,ci.type="bar",xvar=c("sites"),col="darkred")
# small.sample=specaccum(Data3[which(Data3$Years==1 & Data3$size==25),7:109],conditioned=FALSE,method="random")
# plot(small.sample,add=TRUE,ci=2,ci.type="bar",xvar=c("sites"),col="darkblue")
#2013
large.sample=specaccum(Data3[which(Data3$Years==2 & Data3$size==400),7:109],conditioned=FALSE,method="random")
plot(large.sample,add=FALSE,ci=2,ci.type="bar",xvar=c("sites"),col="darkgreen",xlab="",ylab="Number species",ylim=c(0,60))
medium.sample=specaccum(Data3[which(Data3$Years==2 & Data3$size==100),7:109],conditioned=FALSE,method="random")
plot(medium.sample,add=TRUE,ci=2,ci.type="bar",xvar=c("sites"),col="darkred",ylim=c(0,60))
small.sample=specaccum(Data3[which(Data3$Years==2 & Data3$size==25),7:109],conditioned=FALSE,method="random")
plot(small.sample,add=TRUE,ci=2,ci.type="bar",xvar=c("sites"),col="darkblue",ylim=c(0,60))
#2014
large.sample=specaccum(Data3[which(Data3$Years==3 & Data3$size==400),7:109],conditioned=FALSE,method="random")
plot(large.sample,add=FALSE,ci=2,ci.type="bar",xvar=c("sites"),col="darkgreen",xlab="",ylab="Number species",ylim=c(0,60))
medium.sample=specaccum(Data3[which(Data3$Years==3 & Data3$size==100),7:109],conditioned=FALSE,method="random")
plot(medium.sample,add=TRUE,ci=2,ci.type="bar",xvar=c("sites"),col="darkred",ylim=c(0,60))
small.sample=specaccum(Data3[which(Data3$Years==3 & Data3$size==25),7:109],conditioned=FALSE,method="random")
plot(small.sample,add=TRUE,ci=2,ci.type="bar",xvar=c("sites"),col="darkblue",ylim=c(0,60))


dev.off()


}



#°°°°°°°°°°°°°°°°°°°°#
# Remove 2012        #
#°°°°°°°°°°°°°°°°°°°°#

#The design in 2012 was sensibly different (see Appendix) and the mowing treatment had not yet been applied
#Thus, we will treat 2012 independantly for analyses and focus mainly on 2013-2014 for the main analyses.

#remove year 2012 for Quadrat-level and Patch-level data
Data.2012.quadrat = Data3[Data3$Years==1,] #2012 at quadrat-level
Data.2012.patch = Data4[Data4$Years==1,] #2012 at patch-level
Data3 = Data3[Data3$Years!=1,] #Include 2013 and 2014 at quadrat-level
Data4 = Data4[Data4$Years!=1,] #Include 2013 and 2014 at patch-level

#...Add diversity metrics to all datasets
Data.2012.quadrat$plant.richness = specnumber(Data.2012.quadrat[,8:110])
Data.2012.patch$plant.richness = specnumber(Data.2012.patch[,7:109])
Data3$plant.richness = specnumber(Data3[,8:110])
Data4$plant.richness = specnumber(Data4[,7:109])

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°#
# NMDS PREM                  #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°#

#NMDS to visualize changes in time

x = rbind(Data.2012.patch[,7:109],Data4[,7:109])
Years2 = as.factor(c(Data.2012.patch$Years,as.character(Data4$Years)))
colvec = c("red","blue","green")

MDS.mod = metaMDS(decostand(Comp.mat0,"hell"), k=3,autotransform=FALSE,distance="jaccard")
pdf(paste0(fig.path,"NMDS_TEMP.pdf"),width=5,height = 5)
ordiplot(MDS.mod,type="n",choices=c(1,2),main=NULL,xlab="", ylab="")
ordihull(MDS.mod,groups=Years2,show.groups="1",col="red",label=F,lwd=3,lty=1)
ordihull(MDS.mod,groups=Years2,show.groups="2",col="blue",label=F,lwd=3,lty=2)
ordihull(MDS.mod,groups=Years2,show.groups="3",col="green",label=F,lwd=3,lty=2)
points(MDS.mod, dis="sites", pch = 16  , col = colvec[Years2])
dev.off()


#°°°°°°°°°°°°°°°°°°°°°°°°°°°°#
# Extract spatial variables  #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°#

#...Distance between quadrats and islands


#..Load function to converte geographical coordinates into meters
source(paste0(fig.path,"Geographical_distance.R"))

#..Calculate pairwise distance between quadrats within each patch

{ 
#..Organize space dataset per quadrat (currently at plot level)
space$Island = as.factor(space$Island)
space$plot = as.factor(space$plot)
space2 = cast(melt(space[,-2]), Island + quadrat ~ variable, mean)

#..Generate distance matrices
dist.quad = data.frame()
quad.mat = list()
nearest.quad = data.frame()
for(i in 1:36){
  df.patch = data.frame(name = paste(space2[which(space2$Island==i),1],space2[which(space2$Island==i),2],sep="."),
                        lat = space2[which(space2$Island==i),6],
                        lon = space2[which(space2$Island==i),5])
  #Calculate mean distance between quadrats for each patch
  dist.mat = GeoDistanceInMetresMatrix(df.patch)
  mean.dist.quad = mean(dist.mat)
  sd.dist.quad = sd(dist.mat)
  output = data.frame(mean.dist.quad,sd.dist.quad)
  dist.quad = rbind.data.frame(dist.quad,output)
  #Save the distance matrix for each patch into a list
  quad.mat[[paste0("Island", i)]] <- dist.mat
  #Calculate distance to nearest quadrat within each patch
  dist.mat[which(dist.mat == 0)] = NA
  output2 = data.frame(nearest.quadrat = apply(dist.mat,1,function(x) min(x,na.rm=T)))
  nearest.quad = rbind.data.frame(nearest.quad,output2) }

}
  
#..Calculate currencies of interest for further analysis
nearest.quad = rep(nearest.quad$nearest.quadrat,2) #rep 2 times for 2013 and 2014
mean.dist.quad = rep(dist.quad$mean.dist.quad,2)
sd.dist.quad = rep(dist.quad$sd.dist.quad,2)
lala = sapply(quad.mat,mean)
tapply(lala,Data4$size[1:36],mean)
tapply(lala,Data4$size[1:36],sd)
max(unlist(quad.mat))
min(unlist(quad.mat)[which(unlist(quad.mat)>0)])

#..Calculate pairwise distance between patches 

{ 
#..Organize space dataset per island (currently at quadrat level)
space3 = cast(melt(space2), Island ~ variable, mean)

df.island = data.frame(name = space3$Island,
                      lat = space3[,3],
                      lon = space3[,2])

dist.mat.island = GeoDistanceInMetresMatrix(df.island)
str(dist.mat.island)
dist.mat.island = matrix(data=dist.mat.island,nrow=36,ncol=36)
str(dist.mat.island)
dist.mat.island[which(dist.mat.island == 0)] = NA
#test[lower.tri(test,diag=T)] = NA 
nearest.island = apply(dist.mat.island,1,function(x) min(x,na.rm=T))
}

#..Calculate currencies of interest for further analysis
nearest.island = rep(nearest.island,2)


#...PCNM analysis to extract intrinsic spatial structure

{ 
x = rbind(Data.2012.patch[,7:109],Data4[,7:109])
pcnm1=pcnm(dist(treatment[,8:9]))

#..select most important pcnm axes for...
library(adespatial)
#...species composition for each year
forward.sel(x[Years2==1,],scores(pcnm1)) #1-3-5
forward.sel(x[Years2==2,],scores(pcnm1)) #none
forward.sel(x[Years2==3,],scores(pcnm1)) #2-3-18
ordisurf(treatment[,8:9], scores(pcnm1, choi=1), bubble = 4, main=NULL,xlab="Longitude",ylab="Latitude")
ordisurf(treatment[,8:9], scores(pcnm1, choi=2), bubble = 4, main=NULL,xlab="Longitude",ylab="Latitude")
ordisurf(treatment[,8:9], scores(pcnm1, choi=3), bubble = 4, main=NULL,xlab="Longitude",ylab="Latitude")
#...plant species richness for each year
forward.sel(Data.2012.patch$plant.richness,scores(pcnm1))#2-6
forward.sel(Data4$plant.richness[Data4$Years==2],scores(pcnm1))#3
forward.sel(Data4$plant.richness[Data4$Years==3],scores(pcnm1)) #21-18-6
pdf(paste0(fig.path,"PCNM2.pdf"))
ordisurf(treatment[,8:9], scores(pcnm1, choi=2), bubble = 4, main=NULL,xlab="Longitude",ylab="Latitude")
dev.off()
ordisurf(treatment[,8:9], scores(pcnm1, choi=3), bubble = 4, main=NULL,xlab="Longitude",ylab="Latitude")
ordisurf(treatment[,8:9], scores(pcnm1, choi=21), bubble = 4, main=NULL,xlab="Longitude",ylab="Latitude")

}

#..Extract selected axes

pcnm.var1 = rep(scores(pcnm1)[,1],2)
pcnm.var2 = rep(scores(pcnm1)[,2],2)
pcnm.var3 = rep(scores(pcnm1)[,3],2)
pcnm.var21 = rep(scores(pcnm1)[,21],2)



#########################################################################
################# DATA ANALYSES
#########################################################################


#°°°°°°°°°°°°°°°°°°°°°°°°°°°°#
# Statitics - plant richness #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°#

###################################
######.....QUADRAT-LEVEL - 2013-2014

with(Data3,plot(density((plant.richness))))
summary(Data3$plant.richness)

#Initial model
mod1 = with(Data3,lme(plant.richness ~ size*nitrogen*mowing*Years + dist.mainland + nearest.quad,random = ~ as.factor(Years)|Island/quadrat,method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim")))
summary(mod1)$tTable
anova(mod1)
summary(mod1)$AIC

{ 
# Incremental model simplification procedure (see Methods) 

#Effect of sampling year (Starting with 4-way interactions)
mod1.2 = with(Data3,update(mod1, ~. - size:nitrogen:mowing:Years))
anova(mod1,mod1.2)


mod1.3 = with(Data3,update(mod1.2, ~. - nitrogen:mowing:Years))
anova(mod1.2,mod1.3)


mod1.4 = with(Data3,update(mod1.3, ~. - size:mowing:Years))
anova(mod1.3,mod1.4)


mod1.5 = with(Data3,update(mod1.4, ~. - size:nitrogen:Years))
anova(mod1.4,mod1.5)


mod1.6 = with(Data3,update(mod1.5, ~. - size:Years)) 
anova(mod1.5,mod1.6) 


mod1.7 = with(Data3,update(mod1.6, ~. - mowing:Years))
anova(mod1.6,mod1.7)


mod1.8 = with(Data3,update(mod1.7, ~. - nitrogen:Years))
anova(mod1.7,mod1.8)

mod1.8.1 = with(Data3,update(mod1.8, ~. - Years))
anova(mod1.8.1,mod1.8)

#Interactions among perturbations (starting with three-way interaction)
mod1.9 = with(Data3,update(mod1.8, ~. - size:nitrogen:mowing))
anova(mod1.8.1,mod1.9)

mod1.10 = with(Data3,update(mod1.9, ~. - nitrogen:mowing))
anova(mod1.9,mod1.10)

mod1.11 = with(Data3,update(mod1.10, ~. - nitrogen:size))
anova(mod1.10,mod1.11)

mod1.12 = with(Data3,update(mod1.11, ~. - mowing:size))
anova(mod1.12,mod1.11)

mod1.13 = with(Data3,update(mod1.12, ~. - mowing))
anova(mod1.13,mod1.12)

mod1.14 = with(Data3,update(mod1.13, ~. - size))
anova(mod1.14,mod1.13)
anova(mod1.14)

mod1.15 = with(Data3,update(mod1.14, ~. - dist.mainland))
anova(mod1.14,mod1.15)
anova(mod1.15)

mod1.16 = with(Data3,update(mod1.15, ~. - nearest.quad))
anova(mod1.15,mod1.16)
anova(mod1.16)

}

####Final model with REML estimator

mod1.17 = with(Data3,lme(plant.richness ~ nitrogen + Years,random = ~ as.factor(Years)|Island/quadrat,method="REML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim")))
summary(mod1.17)$tTable
anova(mod1.17)

#               numDF denDF  F-value p-value
# (Intercept)     1   119 939.2878  <.0001
# nitrogen        1    34  27.3873  <.0001
# Years           1   119   4.2660  0.0411


# Value Std.Error  DF   t-value      p-value
# (Intercept) 11.3292766 1.1184588 119 10.129364 9.111027e-18
# nitrogenF+  -3.8991857 0.7450734  34 -5.233291 8.553736e-06
# Years        0.8053733 0.3899317 119  2.065421 4.105463e-02

##Extract model prediction + figues
#Extract and plot model prediction + CI (http://stackoverflow.com/questions/14358811/extract-prediction-band-from-lme-fit)
{ 
X = rep(seq(2,3,len=2),times=1,each=120)
new.dat = data.frame(Years=X,nitrogen= Data3$nitrogen)
new.dat$pred = predict(mod1.15,newdata=new.dat,level=0)
Designmat = model.matrix(eval(eval(mod1.15$call$fixed)[-2]),new.dat[-ncol(new.dat)])
predvar = diag(Designmat %*% mod1.15$varFix %*% t(Designmat))
new.dat$SE = sqrt(predvar)

colo.s = c("#99999930", "#E69F0030", "#56B4E930", "#009E7330","#CC79A730")
colo.l = c("#999999", "#E69F00", "#56B4E9", "#009E73","#CC79A7")
with(Data3,plot(Years,plant.richness,type="n",xaxt="n"))#xaxt="n",yaxt="n"
axis(side = 1, at =c(2,3),labels=c(2013,2014))

x1 = as.numeric(new.dat$Years[which(new.dat$nitrogen=="F+")])
y = new.dat$pred[which(new.dat$nitrogen=="F+")]
y1.1 = new.dat$pred[which(new.dat$nitrogen=="F+")]+new.dat$SE[which(new.dat$nitrogen=="F+")]
y2.1 = new.dat$pred[which(new.dat$nitrogen=="F+")]-new.dat$SE[which(new.dat$nitrogen=="F+")]
errbar(x1,y,y1.1,y2.1,add=T,errbar.col="red",col="red",lty=1,pch=16)
points(new.dat$Years[which(new.dat$nitrogen=="F+")],Data3$plant.richness[which(Data3$nitrogen=="F+")],lwd=2,col=alpha("red",0.15),pch=16)

x1 = as.numeric(new.dat$Years[which(new.dat$nitrogen=="F-")])
y = new.dat$pred[which(new.dat$nitrogen=="F-")]
y1.1 = new.dat$pred[which(new.dat$nitrogen=="F-")]+new.dat$SE[which(new.dat$nitrogen=="F-")]
y2.1 = new.dat$pred[which(new.dat$nitrogen=="F-")]-new.dat$SE[which(new.dat$nitrogen=="F-")]
errbar(x1,y,y1.1,y2.1,add=T,errbar.col="blue",col="blue",lty=1,pch=16)
points(x1,Data3$plant.richness[which(Data3$nitrogen=="F-")],lwd=2,col=alpha("blue",0.15),pch=16)

}

pdf(paste0(fig.path,"quadrat_level(13-14).pdf"),width=5,height = 5)
bargraph.CI(nitrogen,plant.richness,legend=T,data=Data3)
bargraph.CI(Years,plant.richness,legend=T,data=Data3)
dev.off()

###################################
######.....Patch-level - 2013-2014

with(Data4,plot(density((plant.richness))))
summary(Data4$plant.richness)

#Initial model
mod3 = with(Data4,lme(plant.richness ~ size*nitrogen*mowing*Years + dist.mainland + nearest.island + pcnm.var1 + pcnm.var2 + pcnm.var21 + pcnm.var3,random = ~ as.factor(Years)|Island,method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim")))
summary(mod3)$tTable
anova(mod3)
summary(mod3)$AIC

# Incremental model simplification procedure 
#Effect of sampling year (Starting with 4-way interactions)

{ 
mod3.2 = with(Data4,update(mod3, ~. - size:nitrogen:mowing:Years))
anova(mod3,mod3.2)

mod3.3 = with(Data4,update(mod3.2, ~. - nitrogen:mowing:Years))
anova(mod3.2,mod3.3)

mod3.4 = with(Data4,update(mod3.3, ~. - size:mowing:Years))
anova(mod3.3,mod3.4)

mod3.5 = with(Data4,update(mod3.4, ~. - size:nitrogen:Years))
anova(mod3.4,mod3.5)

mod3.6 = with(Data4,update(mod3.5, ~. - size:Years)) 
anova(mod3.5,mod3.6) 

mod3.7 = with(Data4,update(mod3.6, ~. - mowing:Years))
anova(mod3.6,mod3.7)

mod3.8 = with(Data4,update(mod3.7, ~. - nitrogen:Years))
anova(mod3.7,mod3.8)

mod3.8.1 = with(Data4,update(mod3.8, ~. - Years))
anova(mod3.8.1,mod3.8)

#Interactions among perturbations (starting with three-way interaction)
mod3.9 = with(Data4,update(mod3.8.1, ~. - size:nitrogen:mowing))
anova(mod3.8.1,mod3.9)


mod3.10 = with(Data4,update(mod3.9, ~. - nitrogen:mowing))
anova(mod3.9,mod3.10)

mod3.11 = with(Data4,update(mod3.10, ~. - nitrogen:size))
anova(mod3.10,mod3.11)

mod3.12 = with(Data4,update(mod3.11, ~. - mowing:size))
anova(mod3.12,mod3.11)

#Spatial variables

mod3.13 = with(Data4,update(mod3.11, ~. - nearest.island ))
anova(mod3.11,mod3.13)
anova(mod3.13)

mod3.14 = with(Data4,update(mod3.13, ~. - pcnm.var1 ))
anova(mod3.13,mod3.14)
anova(mod3.14)

mod3.15 = with(Data4,update(mod3.14, ~. - dist.mainland))
anova(mod3.14,mod3.15)
anova(mod3.15)

mod3.16 = with(Data4,update(mod3.15, ~. -pcnm.var21 ))
anova(mod3.15,mod3.16)

mod3.17 = with(Data4,update(mod3.15, ~. -pcnm.var3 ))
anova(mod3.15,mod3.17)
anova(mod3.17)

 }

#Rerun Final models with REML estimators
#SR
mod3.18 = with(Data4,lme(plant.richness ~ size*mowing + nitrogen + pcnm.var2 + pcnm.var21,random = ~ as.factor(Years)|Island,method="REML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim")))
summary(mod3.18)$tTable
anova(mod3.18)


# Value Std.Error DF    t-value      p-value
# (Intercept)      14.856699  1.483849 36 10.0122681 6.012159e-12
# size100           4.152456  1.814920 27  2.2879560 3.019504e-02
# size400           6.004959  1.841169 27  3.2614929 2.998496e-03
# mowingD+          2.485606  1.951651 27  1.2735913 2.136637e-01
# nitrogenF+       -3.737489  1.078396 27 -3.4657840 1.783508e-03
# pcnm.var2        -8.773189  3.567325 27 -2.4593186 2.061391e-02
# pcnm.var21       -7.144865  3.545988 27 -2.0149153 5.397079e-02
# size100:mowingD+ -4.048231  2.702713 27 -1.4978396 1.457759e-01
# size400:mowingD+  2.407967  2.633784 27  0.9142615 3.686724e-01


# numDF denDF   F-value p-value
# (Intercept)     1    36 1126.7401  <.0001
# size            2    27   17.3914  <.0001
# mowing          1    27    0.8101  0.3760
# nitrogen        1    27   20.7724  0.0001
# pcnm.var2       1    27    8.9716  0.0058
# pcnm.var21      1    27    2.6084  0.1179
# size:mowing     2    27    3.1070  0.0610

pdf(paste0(fig.path,"patch_level(13-14).pdf"),width=5,height = 5)
bargraph.CI(size,plant.richness,legend=T,data=Data4)
bargraph.CI(nitrogen,plant.richness,legend=T,data=Data4)
dev.off()


#Extract and plot model prediction + CI (http://stackoverflow.com/questions/14358811/extract-prediction-band-from-lme-fit)
{ 
X= rep(seq(2,3,len=2),times=1,each=72)
new.dat3 = data.frame(nitrogen= Data4$nitrogen,mowing=Data4$mowing,size=Data4$size)
new.dat3$pred = predict(mod3.15,newdata=new.dat3,level=0)
Designmat3 = model.matrix(eval(eval(mod3.15$call$fixed)[-2]),new.dat3[-ncol(new.dat3)])
predvar3 = diag(Designmat3 %*% mod3.15$varFix %*% t(Designmat3))
new.dat3$SE = sqrt(predvar3)

#Fertilized
colo.s = c("#99999930", "#E69F0030", "#56B4E930", "#009E7330","#CC79A730")
colo.l = c("#999999", "#E69F00", "#56B4E9", "#009E73","#CC79A7")
with(Data4,plot(as.numeric(size),plant.richness,type="n",xaxt="n"))#xaxt="n",yaxt="n"
axis(side = 1, at =c(1,2,3),labels=c(25,100,400))

#Fertilized - UnMowed
x1 = as.numeric(new.dat3$size[which(new.dat3$mowing=="D-" & new.dat3$nitrogen=="F+")])
y=new.dat3$pred[which(new.dat3$mowing=="D-" & new.dat3$nitrogen=="F+")]
y1.1 = new.dat3$pred[which(new.dat3$mowing=="D-" & new.dat3$nitrogen=="F+")]+new.dat3$SE[which(new.dat3$mowing=="D-" & new.dat3$nitrogen=="F+")]
y2.1 = new.dat3$pred[which(new.dat3$mowing=="D-" & new.dat3$nitrogen=="F+")]-new.dat3$SE[which(new.dat3$mowing=="D-" & new.dat3$nitrogen=="F+")]
errbar(x1,y,y1.1,y2.1,add=T,errbar.col="red",col="red",lty=1,pch=16)
points(new.dat3$size[which(new.dat3$mowing=="D-" & new.dat3$nitrogen=="F+")],Data4$plant.richness[which(Data4$mowing=="D-" & Data4$nitrogen=="F+")],lwd=2,col=alpha("red",0.15),pch=16)

#Fertilized - mowed
x1 = as.numeric(new.dat3$size[which(new.dat3$mowing=="D+" & new.dat3$nitrogen=="F+")])
y=new.dat3$pred[which(new.dat3$mowing=="D+" & new.dat3$nitrogen=="F+")]
y1.1 = new.dat3$pred[which(new.dat3$mowing=="D+" & new.dat3$nitrogen=="F+")]+new.dat3$SE[which(new.dat3$mowing=="D+" & new.dat3$nitrogen=="F+")]
y2.1 = new.dat3$pred[which(new.dat3$mowing=="D+" & new.dat3$nitrogen=="F+")]-new.dat3$SE[which(new.dat3$mowing=="D+" & new.dat3$nitrogen=="F+")]
errbar(x1,y,y1.1,y2.1,add=T,errbar.col="green",col="green",lty=1,pch=17)
points(new.dat3$size[which(new.dat3$mowing=="D+" & new.dat3$nitrogen=="F+")],Data4$plant.richness[which(Data4$mowing=="D+" & Data4$nitrogen=="F+")],lwd=2,col=alpha("green",0.15),pch=17)

#Unfertilized
with(Data4,plot(as.numeric(size),plant.richness,type="n",xaxt="n"))#xaxt="n",yaxt="n"
axis(side = 1, at =c(1,2,3),labels=c(25,100,400))

#Unfertilized - Unmowed
x1 = as.numeric(new.dat3$size[which(new.dat3$mowing=="D-" & new.dat3$nitrogen=="F-")])
y=new.dat3$pred[which(new.dat3$mowing=="D-" & new.dat3$nitrogen=="F-")]
y1.1 = new.dat3$pred[which(new.dat3$mowing=="D-" & new.dat3$nitrogen=="F-")]+new.dat3$SE[which(new.dat3$mowing=="D-" & new.dat3$nitrogen=="F-")]
y2.1 = new.dat3$pred[which(new.dat3$mowing=="D-" & new.dat3$nitrogen=="F-")]-new.dat3$SE[which(new.dat3$mowing=="D-" & new.dat3$nitrogen=="F-")]
errbar(x1,y,y1.1,y2.1,add=T,errbar.col="red",col="red",lty=1,pch=16)
points(new.dat3$size[which(new.dat3$mowing=="D-" & new.dat3$nitrogen=="F-")],Data4$plant.richness[which(Data4$mowing=="D-" & Data4$nitrogen=="F-")],lwd=2,col=alpha("red",0.15),pch=16)

#Unfertilized - mowed
x1 = as.numeric(new.dat3$size[which(new.dat3$mowing=="D+" & new.dat3$nitrogen=="F-")])
y=new.dat3$pred[which(new.dat3$mowing=="D+" & new.dat3$nitrogen=="F-")]
y1.1 = new.dat3$pred[which(new.dat3$mowing=="D+" & new.dat3$nitrogen=="F-")]+new.dat3$SE[which(new.dat3$mowing=="D+" & new.dat3$nitrogen=="F-")]
y2.1 = new.dat3$pred[which(new.dat3$mowing=="D+" & new.dat3$nitrogen=="F-")]-new.dat3$SE[which(new.dat3$mowing=="D+" & new.dat3$nitrogen=="F-")]
errbar(x1,y,y1.1,y2.1,add=T,errbar.col="green",col="green",lty=1,pch=17)
points(new.dat3$size[which(new.dat3$mowing=="D+" & new.dat3$nitrogen=="F-")],Data4$plant.richness[which(Data4$mowing=="D+" & Data4$nitrogen=="F-")],lwd=2,col=alpha("green",0.15),pch=17)

}


###################################
######.....QUADRAT-LEVEL - 2012

with(Data.2012.quadrat,plot(density((plant.richness))))
summary(Data.2012.quadrat$plant.richness)

#Initial model
mod4 = with(Data.2012.quadrat,lme(plant.richness ~ size*nitrogen + dist.mainland, random = ~ 1|Island/quadrat,method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim")))
summary(mod4)$tTable
summary(mod4)$AIC

{ 
mod4.1 = with(Data.2012.quadrat,lme(plant.richness ~ size+nitrogen + dist.mainland, random = ~ 1|Island/quadrat,method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim")))
summary(mod4.1)$tTable
summary(mod4.1)$AIC

mod4.2 = with(Data.2012.quadrat,lme(plant.richness ~ nitrogen + dist.mainland, random = ~ 1|Island/quadrat,method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim")))
summary(mod4.2)$tTable
summary(mod4.2)$AIC

}

#Rerun Final models with REML estimators
#SR
mod4.3 = with(Data.2012.quadrat,lme(plant.richness ~  dist.mainland, random = ~ 1|Island/quadrat,method="REML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim")))
summary(mod4.3)$tTable
anova(mod4.3)

#                  numDF denDF   F-value p-value
# (Intercept)       1    96 153.66183  <.0001
# dist.mainland     1    34   4.62629  0.0387

# Value  Std.Error DF   t-value      p-value
# (Intercept)    3.567572285 0.42950549 96  8.306232 6.362815e-13
# dist.mainland -0.004840173 0.00225032 34 -2.150881 3.868068e-02

pdf(paste0(fig.path,"quadrat_level(12).pdf"),width=5,height = 5)
plot(plant.richness~dist.mainland,data=Data.2012.quadrat,pch=16)
#lo <- loess(plant.richness~dist.mainland,data=Data.2012.quadrat)
smoothingSpline = with(Data.2012.quadrat,smooth.spline(dist.mainland,plant.richness,spar=1))
#lines(predict(lo), col='red', lwd=2)
lines(smoothingSpline,col='red',lwd=2)
dev.off()

# Smoothing Parameter  spar= 1  lambda= 0.1299019
# Equivalent Degrees of Freedom (Df): 3.041563
# Penalized Criterion (RSS): 210.2833
# GCV: 2.769731


###################################
######.....Patch-LEVEL - 2012
with(Data.2012.patch,plot(density((log(plant.richness+1)))))
summary(Data.2012.patch$plant.richness)

#Initial model
nearest = nearest.island[1:36]
pcnm.var1.2 = pcnm.var1[1:36]
pcnm.var2.2 =pcnm.var2[1:36]
pcnm.var21.2 = pcnm.var21[1:36]
pcnm.var3.2 = pcnm.var3[1:36]
mod5 = with(Data.2012.patch,lme(plant.richness ~ size*nitrogen + dist.mainland + nearest + pcnm.var1.2 + pcnm.var2.2 + pcnm.var21.2 + pcnm.var3.2, random = ~ 1|Island,method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim")))
summary(mod5)$tTable
summary(mod5)$AIC

{ 
#Remove two-way interaction
mod5.1 = with(Data.2012.patch,update(mod5, ~. - size:nitrogen))
anova(mod5,mod5.1)
anova(mod5.1)

#Remove smallest F-values first
mod5.2 = with(Data.2012.patch,update(mod5.1, ~. - pcnm.var21.2))
anova(mod5.1,mod5.2)
anova(mod5.2)

mod5.3 = with(Data.2012.patch,update(mod5.2, ~. - pcnm.var3.2))
anova(mod5.2,mod5.3)
anova(mod5.3)

mod5.4 = with(Data.2012.patch,update(mod5.3, ~. -size))
anova(mod5.3,mod5.4)
anova(mod5.4)

mod5.5 = with(Data.2012.patch,update(mod5.4, ~. -nearest))
anova(mod5.4,mod5.5)
anova(mod5.5)

mod5.6 = with(Data.2012.patch,update(mod5.5, ~. -nitrogen))
anova(mod5.5,mod5.6)

mod5.7 = with(Data.2012.patch,update(mod5.5, ~. -pcnm.var1.2))
anova(mod5.5,mod5.7)
anova(mod5.7)

}

#Rerun Final models with REML estimators
#SR
mod5.8 = with(Data.2012.patch,lme(plant.richness ~ nitrogen + dist.mainland +  pcnm.var2.2, random = ~ 1|Island,method="REML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim")))
summary(mod5.8)$tTable
anova(mod5.8)

#                   Value   Std.Error DF     t-value    p-value
# (Intercept)    3.6772150476 1.134193140 32  3.24214185 0.00277213
# nitrogenF+     1.7573189778 0.922834938 32  1.90426143 0.06590017
# dist.mainland  0.0005153477 0.005591779 32  0.09216166 0.92714416
# pcnm.var2.2   -8.8937435374 3.390996215 32 -2.62275242 0.01324685

#                 numDF denDF   F-value p-value
# (Intercept)       1    32 104.30893  <.0001
# nitrogen          1    32   1.97854  0.1692
# dist.mainland     1    32   3.07813  0.0889
# pcnm.var2.2       1    32   6.87883  0.0132


#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°#
# Statitics - Beta-diversity       #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°#


#############################
#  PERMANOVA-PERMDISP
#################

###################
#Extract presence-absence matrix

#Quadrat-level
Comp.mat1 = Data3[,8:110]
Comp.mat1[Comp.mat1>0] = 1

#Patch-level
Comp.mat2 = Data4[,7:109]
Comp.mat2[Comp.mat2>0] = 1

###################
#Quadrat-level

PERM.mod1 = with(Data3,adonis(Comp.mat1 ~ size*nitrogen*mowing*Years + dist.mainland + nearest.quad,strata= quadrat %in% Island,method="jaccard",permutation=999))

PERM.mod1
# Number of permutations: 999
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# size                         2     1.347  0.6736  4.2171 0.02557  0.001 ***
#   nitrogen                     1     3.949  3.9491 24.7245 0.07496  0.001 ***
#   mowing                       1     0.979  0.9789  6.1284 0.01858  0.001 ***
#   Years                        1     4.264  4.2641 26.6966 0.08094  0.001 ***
#   dist.mainland                1     0.908  0.9075  5.6817 0.01723  0.001 ***
#   nearest.quad                 1     0.315  0.3154  1.9745 0.00599  0.015 *  
#   size:nitrogen                2     0.914  0.4571  2.8621 0.01735  0.001 ***
#   size:mowing                  2     1.096  0.5480  3.4307 0.02080  0.001 ***
#   nitrogen:mowing              1     0.395  0.3947  2.4713 0.00749  0.006 ** 
#   size:Years                   2     0.483  0.2413  1.5107 0.00916  0.034 *  
#   nitrogen:Years               1     0.961  0.9608  6.0156 0.01824  0.001 ***
#   mowing:Years                 1     0.595  0.5946  3.7227 0.01129  0.001 ***
#   size:nitrogen:mowing         2     0.908  0.4542  2.8438 0.01724  0.001 ***
#   size:nitrogen:Years          2     0.385  0.1925  1.2052 0.00731  0.194    
# size:mowing:Years            2     0.384  0.1922  1.2030 0.00729  0.202    
# nitrogen:mowing:Years        1     0.216  0.2164  1.3548 0.00411  0.139    
# size:nitrogen:mowing:Years   2     0.401  0.2006  1.2562 0.00762  0.130    
# Residuals                  214    34.181  0.1597         0.64882           
# Total                      239    52.682                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#PERMDISP
dist1 = vegdist(Comp.mat1,"jaccard")

DISP.mod1 = betadisper(dist1,Data3$size)
anova(DISP.mod1)

DISP.mod2 = betadisper(dist1,Data3$nitrogen)
mean(DISP.mod2$distances[which(Data3$nitrogen=="F-")])

###################
#Patch-level

PERM.mod2 = with(Data4,adonis(Comp.mat2 ~ size*nitrogen*mowing*Years + dist.mainland + nearest.island + pcnm.var1 + pcnm.var2 + pcnm.var21 + pcnm.var3,strata = Island,method="jaccard",permutations=999))
PERM.mod2

#                            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# size                        2    0.9549 0.47745  4.1007 0.07284  0.001 ***
#   nitrogen                    1    1.1022 1.10225  9.4669 0.08408  0.001 ***
#   mowing                      1    0.4101 0.41006  3.5219 0.03128  0.001 ***
#   Years                       1    1.5026 1.50261 12.9055 0.11463  0.001 ***
#   dist.mainland               1    0.3018 0.30180  2.5921 0.02302  0.001 ***
#   nearest.island              1    0.1579 0.15788  1.3560 0.01204  0.001 ***
#   pcnm.var1                   1    0.4011 0.40114  3.4453 0.03060  0.001 ***
#   pcnm.var2                   1    0.4117 0.41167  3.5357 0.03140  0.001 ***
#   pcnm.var21                  1    0.1867 0.18673  1.6038 0.01424  0.001 ***
#   pcnm.var3                   1    0.1106 0.11060  0.9499 0.00844  0.001 ***
#   size:nitrogen               2    0.3867 0.19333  1.6605 0.02950  0.001 ***
#   size:mowing                 2    0.4520 0.22599  1.9410 0.03448  0.001 ***
#   nitrogen:mowing             1    0.2014 0.20140  1.7298 0.01536  0.001 ***
#   size:Years                  2    0.2228 0.11138  0.9566 0.01699  0.380    
# nitrogen:Years              1    0.2600 0.25996  2.2327 0.01983  0.031 *  
#   mowing:Years                1    0.1998 0.19983  1.7163 0.01524  0.077 .  
# size:nitrogen:mowing        2    0.3761 0.18806  1.6152 0.02869  0.001 ***
#   size:nitrogen:Years         2    0.2135 0.10677  0.9170 0.01629  0.432    
# size:mowing:Years           2    0.2124 0.10621  0.9122 0.01620  0.440    
# nitrogen:mowing:Years       1    0.0722 0.07221  0.6202 0.00551  0.733    
# size:nitrogen:mowing:Years  2    0.0822 0.04108  0.3528 0.00627  0.998    
# Residuals                  42    4.8901 0.11643         0.37304           
# Total                      71   13.1088                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#PERMDISP
dist2 = vegdist(Comp.mat2,"jaccard")

DISP.mod1 = betadisper(dist2,Data4$size)
anova(DISP.mod1)

DISP.mod2 = betadisper(dist2,interaction(Data4$size,Data4$mowing))
anova(DISP.mod2)

DISP.mod3 = betadisper(dist2,Data4$nitrogen)
anova(DISP.mod3)
mean(DISP.mod3$distances[which(Data4$nitrogen=="F+")])


#NMDS to visualize effect of Nitrogen addition
MDS.mod1 = metaMDS(Comp.mat2, k=3,autotransform=FALSE,distance="jaccard")
ordiplot(MDS.mod1,type="n",choices=c(1,2),main=NULL,xlab="", ylab="")
ordihull(MDS.mod1,groups=Data4$nitrogen,show.groups="F+",col="gray88",label=F,lwd=3,lty=1)
ordihull(MDS.mod1,groups=Data4$nitrogen,show.groups="F-",col="black",label=F,lwd=3,lty=2)
text(MDS.mod1, display="species", col="black",cex=0.5)
points(MDS.mod1, display="sites", col=c("red"),cex=0.5)

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°#
# Alpha-Beta relationship          #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°#

#Based on Catano, Dickson, and Myers, “Dispersal and Neutral Sampling Mediate Contingent Effects of Disturbance on Plant Beta-Diversity.”

###############
#Quadrat-level
Comp.mat1 = Data3[,8:110]
Comp.mat1[Comp.mat1>0] = 1

#PERMDISP
dist2 = vegdist(Comp.mat1,"jaccard")
dist2.null = raupcrick(Comp.mat1,chase=F,nsimul=999)

fert.mod = betadisper(dist2,Data3$nitrogen,type="centroid")
size.mod = betadisper(dist2,Data3$size,type="centroid")
fert.mod.null = betadisper(dist2.null,Data3$nitrogen,type="centroid")
size.mod.null = betadisper(dist2.null,Data3$size,type='centroid')
permutest(size.mod.null)

#Extract distances from PERMDISP (now raup crick and jaccard have the same range value)
beta.fert.dist.null = fert.mod.null$distances
beta.fert.dist = fert.mod$distances
beta.size.dist.null = size.mod.null$distances
beta.size.dist = size.mod$distances


#Calculate Log Response Ratio for Fertilization effect
beta.fert.effect = log(beta.fert.dist[which(Data3$nitrogen=="F+")]/beta.fert.dist[which(Data3$nitrogen=="F-")])
beta.fert.effect.null = log(beta.fert.dist.null[which(Data3$nitrogen=="F+")]/beta.fert.dist.null[which(Data3$nitrogen=="F-")])

#***
#At quadrat-level there is an uneven number of quadrats in largest (120 total) versus smallest (48 total) patches
#therefore we resampled 10 000 times to generate the figures. 

mean.size.effect.null = 0
mean.size.effect = 0
sd.size.effect.null = 0
sd.size.effect = 0
for(i in 1:10000){
  
  beta.size.effect = log(sample(beta.size.dist[which(Data3$size=="400")],48)/beta.size.dist[which(Data3$size=="25")])
  beta.size.effect.null = log(sample(beta.size.dist.null[which(Data3$size=="400")],48)/beta.size.dist.null[which(Data3$size=="25")])
  
  mean.size.effect[i] = mean(beta.size.effect)
  sd.size.effect[i] = sd(beta.size.effect)
  mean.size.effect.null[i] = mean(beta.size.effect.null)
  sd.size.effect.null[i] = sd(beta.size.effect.null)
}

#FIGURE

#Calculate 95% IC 
error.fert <- qnorm(0.975)*sd(beta.fert.effect)/sqrt(36)
error.fert.null <- qnorm(0.975)*sd(beta.fert.effect.null)/sqrt(36)

error.size <- qnorm(0.975)*mean(sd.size.effect)/sqrt(48)
error.size.null <- qnorm(0.975)*mean(sd.size.effect.null)/sqrt(48)

#Generate figure
pdf(paste0(fig.path,"Beta_Alpha_Quadrat.pdf"),width=5,height=5)

plot(beta.fert.effect.null,type="n",ylab="Effect size",xaxt="n",main="Nitrogen (neutral sampling)_Quadrat-level",ylim=c(-1,1.2))
abline(h=0,lwd=1,col="gray")
#points(x=50,mean(beta.size.effect))
errbar(x=20,mean(beta.fert.effect),mean(beta.fert.effect)+error.fert,mean(beta.fert.effect)-error.fert,add=T,errbar.col="blue",col="blue",lty=1,pch=16)
errbar(x=100,mean(beta.fert.effect.null),mean(beta.fert.effect.null)+error.fert.null,mean(beta.fert.effect.null)-error.fert.null,add=T,errbar.col="orange",col="orange",lty=1,pch=16)

plot(beta.size.effect.null,type="n",ylab="Effect size",xaxt="n",main="Patch size (selection)_Quadrat-level",ylim=c(-1,1.2))
abline(h=0,lwd=1,col="gray")
#points(x=50,mean(beta.size.effect))
errbar(x=10,mean(mean.size.effect),mean(mean.size.effect)+error.size,mean(mean.size.effect)-error.size,add=T,errbar.col="blue",col="blue",lty=1,pch=16)
errbar(x=40,mean(mean.size.effect.null),mean(mean.size.effect.null)+error.size.null,mean(mean.size.effect.null)-error.size.null,add=T,errbar.col="orange",col="orange",lty=1,pch=16)

dev.off()

#############
#Patch-level

Comp.mat2 = Data4[,7:109]
Comp.mat2[Comp.mat2>0] = 1

#PERMDISP
dist2 = vegdist(Comp.mat2,"jaccard")
dist2.null = raupcrick(Comp.mat2,chase=T,nsimul=999)

fert.mod = betadisper(dist2,Data4$nitrogen,type="centroid")
size.mod = betadisper(dist2,Data4$size,type="centroid")
fert.mod.null = betadisper(dist2.null,Data4$nitrogen,type="centroid")
size.mod.null = betadisper(dist2.null,Data4$size,type='centroid')


#Extract distance from centroid (beta-div) from PERMDISP models 
beta.fert.dist.null = fert.mod.null$distances 
beta.fert.dist = fert.mod$distances
beta.size.dist.null = size.mod.null$distances
beta.size.dist = size.mod$distances

#Log Response Ratio
beta.fert.effect = log(beta.fert.dist[which(Data4$nitrogen=="F+")]/beta.fert.dist[which(Data4$nitrogen=="F-")])
beta.fert.effect.null = log(beta.fert.dist.null[which(Data4$nitrogen=="F+")]/beta.fert.dist.null[which(Data4$nitrogen=="F-")])

beta.size.effect = log(beta.size.dist[which(Data4$size=="400")]/beta.size.dist[which(Data4$size=="25")])
beta.size.effect.null = log(beta.size.dist.null[which(Data4$size=="400")]/beta.size.dist.null[which(Data4$size=="25")])

#Calculate 95% IC
error.fert <- qnorm(0.975)*sd(beta.fert.effect)/sqrt(36)
error.fert.null <- qnorm(0.975)*sd(beta.fert.effect.null)/sqrt(36)

error.size <- qnorm(0.975)*sd(beta.size.effect)/sqrt(24)
error.size.null <- qnorm(0.975)*sd(beta.size.effect.null)/sqrt(24)

#Figure
pdf(paste0(fig.path,"Beta_Alpha_Patch.pdf"),width=5,height = 5)

plot(beta.fert.effect.null,type="n",ylab="Effect size",xaxt="n",main="Nitrogen (neutral sampling)_Patch-level",ylim=c(-1,1.2))
abline(h=0,lwd=1,col="gray")
#points(x=50,mean(beta.size.effect))
errbar(x=5,mean(beta.fert.effect),mean(beta.fert.effect)+error.fert,mean(beta.fert.effect)-error.fert,add=T,errbar.col="blue",col="blue",lty=1,pch=16)
errbar(x=30,mean(beta.fert.effect.null),mean(beta.fert.effect.null)+error.fert.null,mean(beta.fert.effect.null)-error.fert.null,add=T,errbar.col="orange",col="orange",lty=1,pch=16)

plot(beta.size.effect.null,type="n",ylab="Effect size",xaxt="n",main="Patch size (selection)_Patch-level",ylim=c(-1,1.2))
abline(h=0,lwd=1,col="gray")
#points(x=50,mean(beta.size.effect))
errbar(x=5,mean(beta.size.effect),mean(beta.size.effect)+error.size,mean(beta.size.effect)-error.size,add=T,errbar.col="blue",col="blue",lty=1,pch=16)
errbar(x=20,mean(beta.size.effect.null),mean(beta.size.effect.null)+error.size.null,mean(beta.size.effect.null)-error.size.null,add=T,errbar.col="orange",col="orange",lty=1,pch=16)

dev.off()


#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°#
# Changes in species composition   #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°#

#############################
#INDICATOR SPECIES ANALYSES
#################

############
#Quadrat-level
size.nitrogen = interaction(Data3$nitrogen,Data3$size)
size.mowing = interaction(Data3$mowing,Data3$size)
mowing.nitrogen = interaction(Data3$nitrogen,Data3$mowing)
mowing.nitrogen.size = interaction(Data3$nitrogen,Data3$mowing,Data3$size)

indval = multipatt(Comp.mat1,Data3$nitrogen,control=how(nperm=999),func="IndVal.g")
summary(indval,indvalcomp=TRUE)
round(head(indval$str),2)

indval = multipatt(Comp.mat1,Data3$mowing,control=how(nperm=999))
summary(indval,indvalcomp=TRUE)

indval = multipatt(Comp.mat1,Data3$size,control=how(nperm=999),duleg=T)
summary(indval,indvalcomp=TRUE)

indval = multipatt(Comp.mat1,mowing.nitrogen,control=how(nperm=999),duleg=T)
summary(indval,indvalcomp=TRUE)

indval = multipatt(Comp.mat1,size.nitrogen,control=how(nperm=999),duleg=T)
summary(indval,indvalcomp=TRUE)

indval = multipatt(Comp.mat1,size.mowing,control=how(nperm=999),duleg=T)
summary(indval,indvalcomp=TRUE)

indval = multipatt(Comp.mat1,mowing.nitrogen.size,control=how(nperm=999),duleg=T)
summary(indval,indvalcomp=TRUE)

############
#Patch-level
size.nitrogen = interaction(Data4$nitrogen,Data4$size)
size.mowing = interaction(Data4$mowing,Data4$size)
mowing.nitrogen = interaction(Data4$nitrogen,Data4$mowing)
mowing.nitrogen.size = interaction(Data4$nitrogen,Data4$mowing,Data4$size)

indval = multipatt(Comp.mat2,Data4$nitrogen,control=how(nperm=999))
summary(indval,indvalcomp=TRUE)

indval = multipatt(Comp.mat2,Data4$mowing,control=how(nperm=999))
summary(indval,indvalcomp=TRUE)

indval = multipatt(Comp.mat2,Data4$size,control=how(nperm=999),duleg=T)
summary(indval,indvalcomp=TRUE)

indval = multipatt(Comp.mat2,mowing.nitrogen,control=how(nperm=999),duleg=T)
summary(indval,indvalcomp=TRUE)

indval = multipatt(Comp.mat2,size.nitrogen,control=how(nperm=999),duleg=T)
summary(indval,indvalcomp=TRUE)

indval = multipatt(Comp.mat2,size.mowing,control=how(nperm=999),duleg=T)
summary(indval,indvalcomp=TRUE)

indval = multipatt(Comp.mat2,mowing.nitrogen.size,control=how(nperm=999),duleg=T)
summary(indval,indvalcomp=TRUE)


#RDA to visualize plant composition constrained on by the effects of the treatments

##############
#Quadrat-level
env1 = cbind(Data3[,c(4:6)])
rda.mod1 = with(Data3,(rda(Comp.mat1 ~ size + nitrogen + mowing)))
(fit = envfit(rda.mod1, env1, perm = 999,display="sites"))
plot(rda.mod1,type="n",scaling=3)
orditorp(rda.mod1, display = "species",scaling=3,pch=16,air=0.5,cex=1,pcol=alpha("black",0.15))
# plot(fit, add = T)
text(rda.mod1, display="bp",labels=c("Small","Large","Fertilized","Mowed"), col=alpha("blue",1),cex=0.8)
#Small here is for 100m2 and not for 25 m2


##############
#Patch-level
env2 = cbind(Data4[,c(3:5)])
rda.mod2 = with(Data4,(rda(Comp.mat2 ~ size + nitrogen + mowing)))
(fit = envfit(rda.mod2, env2, perm = 999,display="sites"))
plot(rda.mod2,type="n",scaling=3)
orditorp(rda.mod2, display = "species",scaling=3,pch=16,air=0.5,cex=1,pcol=alpha("black",0.15))
# plot(fit, add = T)
text(rda.mod2, display="bp",labels=c("Small","Large","Fertilized","Mowed"), col=alpha("blue",1),cex=0.8)
#Small here is for 100m2 and not for 25 m2

#########################################################################
################# COMPLEMENATARY INFORMATION (BIOMASS AND LIGHT)
#########################################################################

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°#
# Biomass analysis                 #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°#
library(car)
biomass$nearest.island = rep(nearest.island[1:36],3)

biomass$total = biomass$F+biomass$G
biomass$ratio = logit(biomass$F/biomass$total)

forward.sel(biomass$ratio[biomass$Years=="2014"],scores(pcnm1))#2-17-5
forward.sel(biomass$total[biomass$Years=="2012"],scores(pcnm1))#20-5
forward.sel(biomass$total[biomass$Years=="2013"],scores(pcnm1))
forward.sel(biomass$total[biomass$Years=="2014"],scores(pcnm1))#8

biomass$pcnm2 = rep(pcnm.var2[1:36],3)
biomass$pcnm20 = rep(scores(pcnm1)[,20],3)
biomass$pcnm8 = rep(scores(pcnm1)[,8],3)

#fix structure
biomass$Size = factor(biomass$Size)
biomass$Years = factor(biomass$Years)
str(biomass)

#Total biomass
Temp.tot = lm(log(total+1) ~ Size*Nitrogen*Mowing*Years 
               + Dist + nearest.island + pcnm2 + pcnm8 + pcnm20, data=biomass)
shapiro.test(Temp.tot$res)
summary(Temp.tot)
anova(Temp.tot)

# Response: log(total + 1)
#                             Df  Sum Sq Mean Sq F value    Pr(>F)    
# Size                        2  1.0062  0.5031  2.2215   0.11649    
# Nitrogen                    1  5.7271  5.7271 25.2900 4.026e-06 ***
# Mowing                      1  0.2367  0.2367  1.0451   0.31038    
# Years                       2  0.2505  0.1252  0.5531   0.57783    
# Dist                        1  0.3035  0.3035  1.3402   0.25116    
# nearest.island              1  0.7477  0.7477  3.3017   0.07375 .  
# pcnm2                       1  0.0403  0.0403  0.1781   0.67438    
# pcnm8                       1  0.1315  0.1315  0.5808   0.44872    
# pcnm20                      1  0.8209  0.8209  3.6250   0.06128 .  
# Size:Nitrogen               2  1.8579  0.9289  4.1021   0.02093 *  
# Size:Mowing                 2  0.4240  0.2120  0.9361   0.39729    
# Nitrogen:Mowing             1  0.2207  0.2207  0.9746   0.32713    
# Size:Years                  4  1.7374  0.4344  1.9181   0.11775    
# Nitrogen:Years              2  0.9636  0.4818  2.1276   0.12723    
# Mowing:Years                2  0.2703  0.1351  0.5967   0.55357    
# Size:Nitrogen:Mowing        2  1.4082  0.7041  3.1092   0.05124 .  
# Size:Nitrogen:Years         4  2.7940  0.6985  3.0845   0.02168 *  
#   Size:Mowing:Years           4  0.6968  0.1742  0.7692   0.54909    
# Nitrogen:Mowing:Years       2  0.1147  0.0574  0.2532   0.77702    
# Size:Nitrogen:Mowing:Years  4  1.3877  0.3469  1.5320   0.20312    
# Residuals                  66 14.9462  0.2265   


bargraph.CI(Nitrogen,total,xlab="Nitrogen addition",ylab="Total plant biomass (g)",data=biomass)



#Biomass ratio
Temp.ratio = lm(ratio ~ Size*Nitrogen*Mowing*Years 
               + Dist + nearest.island + pcnm2 + pcnm8 + pcnm20, data=biomass)
shapiro.test(Temp.ratio$res)
summary(Temp.ratio)
anova(Temp.ratio)

# Response: ratio
#                             Df  Sum Sq Mean Sq F value    Pr(>F)    
# Size                        2   3.601   1.801  0.9592  0.388494    
# Nitrogen                    1   0.804   0.804  0.4283  0.515077    
# Mowing                      1  19.536  19.536 10.4067  0.001955 ** 
# Years                       2 271.435 135.718 72.2949 < 2.2e-16 ***
# Dist                        1   0.057   0.057  0.0305  0.861840    
# nearest.island              1  13.628  13.628  7.2593  0.008936 ** 
# pcnm2                       1   1.716   1.716  0.9141  0.342526    
# pcnm8                       1   4.827   4.827  2.5712  0.113597    
# pcnm20                      1   6.417   6.417  3.4183  0.068957 .  
# Size:Nitrogen               2   2.809   1.404  0.7481  0.477225    
# Size:Mowing                 2   0.595   0.297  0.1585  0.853784    
# Nitrogen:Mowing             1   3.215   3.215  1.7128  0.195161    
# Size:Years                  4   6.103   1.526  0.8127  0.521560    
# Nitrogen:Years              2  21.399  10.699  5.6994  0.005209 ** 
# Mowing:Years                2   4.109   2.054  1.0943  0.340754    
# Size:Nitrogen:Mowing        2  20.311  10.156  5.4097  0.006674 ** 
# Size:Nitrogen:Years         4  12.334   3.083  1.6425  0.174050    
# Size:Mowing:Years           4   7.776   1.944  1.0355  0.395618    
# Nitrogen:Mowing:Years       2  10.619   5.309  2.8282  0.066302 .  
# Size:Nitrogen:Mowing:Years  4   3.656   0.914  0.4869  0.745303    
# Residuals                  66 123.900   1.877    


bargraph.CI(Mowing,ratio,xlab="Defoliation",ylab="logit(Forb:Total)",data=biomass)
bargraph.CI(Years,ratio,group=NULL,legend=F,xlab="Sampling year",ylab="logit(Forb:Total)",data=biomass)
bargraph.CI(Nitrogen,ratio,group=Years,legend=FALSE,xlab="Nitrogen addition",ylab="Forb:Total biomass",data=biomass)
bargraph.CI(Mowing,ratio,group=Years,legend=FALSE,xlab="Mowing",ylab="Forb:Total biomass",data=biomass)
plot(ratio~nearest.island,data=biomass)
lines(predict(Temp.ratio)) #artifact! 

#Change in litter biomass

Temp.lit = lm(L ~ Size*Nitrogen*Mowing*Years + 
                 + Dist + nearest.island + pcnm2 + pcnm8 + pcnm20,data=biomass)
shapiro.test(Temp.lit$res)
summary(Temp.lit)
anova(Temp.lit)
# 
# 
# Response: L
#                             Df Sum Sq Mean Sq F value    Pr(>F)    
# Size                        2   5785  2892.3  1.6072  0.208196    
# Nitrogen                    1  12900 12900.0  7.1682  0.009353 ** 
# Mowing                      1  10978 10978.5  6.1005  0.016108 *  
# Years                       2  42614 21306.8 11.8397 4.037e-05 ***
# Dist                        1   2199  2198.6  1.2217  0.273039    
# nearest.island              1     37    36.5  0.0203  0.887138    
# pcnm2                       1    292   291.7  0.1621  0.688518    
# pcnm8                       1   2257  2256.9  1.2541  0.266826    
# pcnm20                      1   5785  5785.1  3.2146  0.077565 .  
# Size:Nitrogen               2   8204  4102.2  2.2795  0.110339    
# Size:Mowing                 2   3150  1575.0  0.8752  0.421560    
# Nitrogen:Mowing             1   1173  1173.1  0.6519  0.422347    
# Size:Years                  4   3526   881.5  0.4898  0.743157    
# Nitrogen:Years              2    851   425.6  0.2365  0.790067    
# Mowing:Years                2   8727  4363.3  2.4246  0.096362 .  
# Size:Nitrogen:Mowing        2   6358  3179.1  1.7666  0.178904    
# Size:Nitrogen:Years         4   6300  1574.9  0.8751  0.483685    
# Size:Mowing:Years           4   4171  1042.8  0.5795  0.678543    
# Nitrogen:Mowing:Years       2   2097  1048.3  0.5825  0.561349    
# Size:Nitrogen:Mowing:Years  4    976   244.1  0.1356  0.968622    
# Residuals                  66 118774  1799.6  



bargraph.CI(Mowing,L,xlab="Mowing year",ylab="Litter biomass (g)",data=biomass)
bargraph.CI(Years,L,xlab="Sampling year",ylab="Litter biomass (g)",data=biomass)
bargraph.CI(Mowing,L,group=Years,legend=TRUE,xlab="Defoliation",ylab="Litter biomass (g)",x.leg=7,y.leg=115,data=biomass)



#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°#
# Light analysis                 #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°#

#...At the quadrat-level
## Average all nested 1 m2 plots
light1 = light[1:480,] #remove mainland plots
light1$plot = as.factor(light1$plot)
str(light1)
light2 = cast(melt(light1),Island + quadrat ~ variable, mean)
desired_order <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36")
light2$Island <- factor(as.character(light2$Island), levels=desired_order )
light3 <- as.data.frame(light2[order(light2$Island),])
light= data.frame(light3,size=Data3$size,nitrogen=Data3$nitrogen,mowing=Data3$mowing)

## ANOVA
plot(density(light$light.av))
plot(density(log(light$light.av)))
model = lm(log(light.av) ~ size*nitrogen*mowing,data=light)
shapiro.test(model$res)
summary(model)
anova(model)

#                        Df Sum Sq Mean Sq F value  Pr(>F)    
# size                   2   2.65    1.33   3.530 0.03093 *  
# nitrogen               1  56.03   56.03 149.052 < 2e-16 ***
# mowing                 1   0.41    0.41   1.100 0.29538    
# size:nitrogen          2   1.21    0.61   1.612 0.20171    
# size:mowing            2   3.43    1.72   4.566 0.01137 *  
# nitrogen:mowing        1   1.23    1.23   3.270 0.07189 .  
# size:nitrogen:mowing   2   4.45    2.22   5.917 0.00312 ** 
# Residuals            228  85.71    0.38                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

bargraph.CI(nitrogen,light.av,xlab="Nitrogen addition",ylab="Light availability",data=light)
bargraph.CI(size,group=mowing,light.av,xlab="Nitrogen addition",ylab="Light availability",data=light,legend=T)



#...At the patch-level

## Average all nested 4 m2 quadrats
light4 = cast(melt(light1),Island ~ variable, mean)
desired_order <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36")
light4$Island <- factor(as.character(light4$Island), levels=desired_order )
light4 <- as.data.frame(light4[order(light4$Island),])
light.patch= data.frame(light4,size=Data4$size[1:36],nitrogen=Data4$nitrogen[1:36],mowing=Data4$mowing[1:36],
                        cv = with(light,tapply(light.av,Island,sd)/tapply(light.av,Island,mean)))

##ANOVA
plot(density(light.patch$cv))
plot(density(log(light.patch$cv)))
model = lm(cv ~ size*nitrogen*mowing,data=light.patch)
shapiro.test(model$res)
summary(model)
anova(model)

#                       Df Sum Sq Mean Sq F value Pr(>F)  
# size                  2 0.0183 0.00915   0.331 0.7212  
# nitrogen              1 0.1217 0.12165   4.404 0.0466 *
# mowing                1 0.0305 0.03049   1.104 0.3039  
# size:nitrogen         2 0.2015 0.10076   3.648 0.0414 *
# size:mowing           2 0.0142 0.00708   0.256 0.7759  
# nitrogen:mowing       1 0.0238 0.02384   0.863 0.3621  
# size:nitrogen:mowing  2 0.0931 0.04657   1.686 0.2065  
# Residuals            24 0.6629 0.02762                        


bargraph.CI(nitrogen,cv,xlab="Nitrogen addition",ylab="Light CV",data=light.patch)

capture.output(summary(model),file="light_CV.doc")

######THE END########





