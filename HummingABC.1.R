library(igraph)
library(reshape2)
library(moments)
library(ggplot2)
library(pROC)
library(PresenceAbsence)
library(tidyverse)
library(knitr)
library(msm)

#######################################
#Define necessary functions
#######################################


#Connectance function
conn1<-function(net){
  n<-nrow(net)
  conn<-length(which(net!=0))/(n*(n-1))
  return(conn)
}

#Functions to build a network based on overlaps in utility curves

myoverlap<-function(a,b){
  a = data.frame( value = a, Source = "a" )
  b = data.frame( value = b, Source = "b" )
  d = rbind(a, b)
  library(scales) 
  d$value <- rescale( d$value, to = c(0,2*pi))
  
  a <- d[d$Source == "a", 1]
  b <- d[d$Source == "b", 1]     
  
  library(overlap)
  
  # generate kernel densities
  da <- density(a, from=min(c(a,b)), to=max(c(a,b)), adjust = 1)
  db <- density(b, from=min(c(a,b)), to=max(c(a,b)), adjust = 1)
  
  # Compute overlap coefficient
  return(overlapTrue(da$y,db$y))
}



birdnet<-function(u){#u is an array of smoothed points of utilized flower resource traits for each bird species (each column)
  mat<-matrix(rep(0,ncol(u)^2),ncol=ncol(u))
  for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
      mat[i,j]<-myoverlap(u[,i],u[,j])
    }
  }
  diag(mat)<-0
  g<-graph.adjacency(mat,weighted=TRUE,mode="undirected")
  cl<-fastgreedy.community((g))
  mod<-max(cl$modularity)
  return(list(mat=mat,cl=cl,g=g,mod=mod))
}



#Load input data
#1. Resource abundance (flowers) mapped over one trait (Corolla length): Over different months and elevation points
#2. Inteaction data: Each data of bird species visiting a flower species:Over different months and elevation points
#3. Flower morphology data
#4. Bird morphology data
interact<-read.csv("HummingbirdInteractions.csv")
interact<-interact%>%select(month,year,ele,Hummingbird,Iplant_Double)
colnames(interact)[c(3,4,5)]<-c("elevation","bird","plant")
bins<-c(1300,1500,1700,1900,2100,2300,2500)
interact$elevation<-.bincode(interact$elevation,bins,FALSE)
interact<-na.omit(interact)
interact$Y<-1

bird.obs<-interact%>%group_by(bird)%>%summarise(Y=sum(Y))%>%
  ungroup()%>%filter(Y>200)

birdlist<-as.vector(unlist((interact%>%group_by(bird)%>%summarise(Y=sum(Y))%>%
                              ungroup()%>%filter(Y>200)%>%select(bird))))
birdlist<-birdlist[order(birdlist)]
nsp<-length(birdlist)

plantlist<-unique(interact$plant)

interact<-interact[which(interact$bird%in%birdlist),]


fl.abun<-read.csv("FlowerTransectClean.csv")
abund<-fl.abun%>%select(month,year,Elevation.Begin,Iplant_Double,Total_Flowers)
colnames(abund)[c(3,4,5)]<-c("elevation","plant","abundance")
abund$elevation<-as.factor(abund$elevation)
levels(abund$elevation)<-c(1,2,3,4,5,6)
abund<-abund%>%group_by(month,year,elevation,plant)%>%summarise(abundance=sum(abundance))
abund<-as.data.frame(na.omit(abund))
abund<-abund[which(abund$plant%in%plantlist),]
abund$abundance<-log(exp(1)+abund$abundance)

flowermorph<-read.csv("FlowerMorphology.csv")
flower.trait<-flowermorph[,c(1,3)]
colnames(flower.trait)<-c("plant","corolla")
flower.trait<-flower.trait[which(flower.trait$plant%in%plantlist),]
add.plant<-setdiff(plantlist,flower.trait$plant)
add.plant<-cbind(add.plant,runif(length(add.plant),min(flower.trait$corolla),max(flower.trait$corolla)))
flower.trait[(nrow(flower.trait)+1):(nrow(flower.trait)+nrow(add.plant)),]<-add.plant     
flower.trait$corolla<-as.numeric(flower.trait$corolla)
flower.trait[,2]<-round(flower.trait[,2],2)


birdmorph<-read.csv("HummingbirdMorphology.csv")
birdbill<-birdmorph[,c(28,10)]
birdbill[23,1]<-birdlist[6]
birdbill[,2]<-round(birdbill[,2],2)
birdbill<-birdbill[match(birdlist,birdbill$English),]
colnames(birdbill)<-c("bird","culmen")

birdwing<-birdmorph[,c(28,16)] #Wing area is shown to affect flight ability correlate negatively with #elevation. (Altshuler et.al. 2004, Stiles et.al. 2005)
birdwing[23,1]<-birdlist[6]
birdwing[,2]<-round(birdwing[,2],2)
colnames(birdwing)<-c("bird","wing")
birdwing<-birdwing[which(birdwing$bird%in%birdlist),]

birdforage<-read.csv("HummingbirdSpecies.csv")
birdforage<-birdforage[,c(1,3)]
birdforage[23,1]<-birdlist[6]
colnames(birdforage)<-c("bird","role")
birdforage<-birdforage[match(birdlist,birdforage$bird),]


interact$corolla<-flower.trait$corolla[match(interact$plant,flower.trait$plant)]
interact$culmen<-birdbill$culmen[match(interact$bird,birdbill$bird)]
totaldat<-interact


#Get the utility curves of all hummingbird species as a function of flower trait values
#These utility curves are based on pooled data over years and locations

int.sum<-totaldat%>%group_by(plant, bird,month,year,elevation)%>%
  summarize(Y=sum(Y), corolla=mean(corolla))%>%ungroup()


int.sum<-merge(int.sum,abund,all.x=TRUE)


miss<-which(is.na(int.sum$abundance))
for(i in 1:length(miss)){
  plant<-int.sum[which(int.sum$plant==int.sum[i,1]),]
  if(!all(is.na(plant$abundance))){
    int.sum$abundance[i]<-mean(plant$abundance,na.rm=TRUE)
  }
}
miss<-which(is.na(int.sum$abundance))

for(j in 1:length(miss)){
  i<-miss[j]
  samp<-int.sum[i,2:4]
  temp<-int.sum[which(int.sum$month==as.numeric(samp[1]) & int.sum$year==as.numeric(samp[2]) & int.sum$elevation==as.numeric(samp[3])),8]
  if(!all(is.na(temp))){
    int.sum$abundance[i]<-mean(temp,na.rm=TRUE)
  }else{
    int.sum$abundance[i]<-int.sum$Y[i]
  }
}

int.sum$Y<-int.sum$Y/int.sum$abundance



resources<-int.sum%>%select(corolla,month,year,elevation,abundance)
#unique timepoints
tps<-unique(resources[,2:4])
tps<-tps[order(tps$year,tps$month,tps$elevation),]


bird.obs<-int.sum%>%select(bird,month,year,elevation,Y)%>%
              group_by(bird,month,year,elevation)%>%
              summarize(Y=sum(Y))%>%ungroup()%>%
              complete(bird,nesting(month,year,elevation),fill=list(Y=0))%>%
              mutate(Y=as.numeric(Y>0))

obs<-bird.obs%>%pivot_wider(names_from=bird,values_from=Y)%>%
                arrange(year,month,elevation)%>%select(-month,-year,-elevation)

obs<-as.data.frame(obs)
##########################################################################
##########################################################################

getsens<-function(obs,pred){
  dat<-data.frame(ID=1:length(obs),obs=obs,pred=pred)
  thresh<-unlist(optimal.thresholds(dat,opt.methods = 3,na.rm=TRUE)[2])
  cmdat<-cmx(DATA=dat,threshold = thresh)
  return(sensitivity(cmdat,st.dev=FALSE))
}
#simulate a competition dynamic from resources and utility curves
#Utility curve has vectors of 500 for each bird species

sim.comp<-function(resources,utility,a.vec){
  #Create an alpha matrix
  alphas<-birdnet(utility)$mat
  alphas.sp<-list()
  for(i1 in 1:ncol(utility)){
    mat<-alphas
    mat[i1,]<-mat[i1,]*a.vec[i1]
    mat[i1,i1]<-1
    alphas.sp<-c(alphas.sp,list(mat))
  }


  kdat<-matrix(nrow=nrow(tps),ncol=ncol(utility))
  occ<-matrix(nrow=nrow(tps),ncol=ncol(utility))
  AUCs<-vector(length=ncol(utility))

  #sensdat<-vector(length=ncol(utility))
  #Estimate carrying capacities
  for(i2 in 1:nrow(tps)){
    dat<-resources%>%filter(month==tps[i2,1],year==tps[i2,2],elevation==tps[i2,3])
    ks<-c()
    for(i3 in 1:ncol(utility)){
      freqs<-table(utility[,i3])
      dat1<-merge(dat,freqs,by.x="corolla",by.y="Var1",all.x=TRUE)
      ks<-c(ks,sum(dat1$abundance*dat1$Freq,na.rm=TRUE))
    }
    for(i4 in 1:ncol(utility)){
      pops.k<-solve(alphas.sp[[i4]],ks)
      pops.k<-pops.k+abs(min(pops.k))
      occ[i2,i4]<-pops.k[i4]
    }
    
  }
  occ<-apply(occ,2,function(i5) i5/max(i5))
  for(i6 in 1:ncol(utility)){
    AUCs[i6]<-roc(obs[,i6],occ[,i6],quiet=TRUE)$auc
  #  sensdat[k]<-getsens(obs[,k],occ[,k])
  }
return(AUCs)
}


######################################################################################
#Parameterize utility curves with two parameters for shape
#Assume that utility curves follow a bounded normal distribution.
#The bounds are fixed by the range of flower traits. 
#The mean and variance of the distribution are assumed to be the functions of bird trait (culmen length)
#Mean preferred corolla= mean corolla + p1*(culmen-mean culmen)
#Variance has a quadriatic relationship with (culmen-min(culmen)) and thus has 3 parameters
#So, we have in total 4+ (no. of bird species) parameters to estimate 
culs<-birdbill$culmen
cors<-flower.trait$corolla
mid.cor<-min(cors)+(0.5*diff(range(cors)))
mid.cul<-min(culs)+(0.5*diff(range(culs)))

#Build a utility curve based on parameters of mean and variances functions
build.util<-function(p1,Vb,nsp){
  mat<-matrix(nrow=500,ncol=length(birdlist))
  m1<-mid.cor+p1*(culs-mid.cul)   
  v1<-5+Vb*(culs-min(culs))
  #v1<-Va*((culs-min(culs))^2)+Vb*(culs-min(culs))+Vc
  for(j1 in 1:nsp){
    probs<-dtnorm(cors,m1[j1],v1[j1],min(cors),max(cors))
    probs[is.na(probs)]<-0
    mat[,j1]<-c(cors,sample(cors,(500-length(cors)),replace=TRUE,prob=probs))
  }
  return(mat)
}

p1.max<-diff(range(cors))/diff(range(culs))
v1.max<-20/diff(range(culs))




#mat[,i]<-c(cors,sample(cors,(500-length(cors)),replace=TRUE,prob=probs))


########################################################################
#Write an ABC algorithm 

#priors:
#a.vec: all elemements have U[0,1]
#p1: U[0,(range of corollas/range of culmens)]
#The quadriatic coefficients for variance function: V=Va*(culmen^2)+Vb*(Culmen)+Vc
#Va~U[-1,0.5]
#Vb~U[-sqrt(20),sqrt(40)]
#Vc~U[5,25]

#Start the cluster and split the loop on #cores
library(parallel)
library(foreach)
library(iterators)
library(doParallel)


numcores<-detectCores()
clust<-makeCluster(min(3,numcores))
clusterExport(clust,c("myoverlap","birdnet","build.util","resources","obs","sim.comp",
                      "cors","mid.cor","mid.cul","culs","p1.max","v1.max","nsp","birdlist"))
registerDoParallel(clust)




paramdat<-foreach(k1=1:50,.combine=rbind,.packages=c("msm","pROC","igraph","reshape2","tidyverse"))%dopar%
{

  
  p1<-runif(1,0,p1.max)
  Vb<-runif(1,0,v1.max)
  #Va<-runif(1,-0.1,0.1)
  #Vc<-runif(1,5,25)
  #Vb<-runif(1,((5-961*Va-Vc)/31),((25-961*Va-Vc)/31))
  a.vec<-runif(nsp,0,1)
  
  res<-sim.comp(resources,build.util(p1,Vb,nsp),a.vec)
  return(c(p1,Vb,a.vec,res))
}

stopCluster(clust)


write.csv(paramdat,"firstparam.csv")
######################


#Second filter
param.main<-paramdat
param.main<-param.main[,-1]

# Choose the rows where all AUCs are above 0.6
good<-apply(param.main[,8:12],1,min)
good.ind<-which(good>0.6)
paramset<-param.main[good.ind,1:7]

#Add new parameter values in the neighborhood of the filtered values.
newset<-paramset
for(i1 in 1:30){
  newset<-rbind(newset,apply(paramset,c(1,2),function(i2) rnorm(1,i2,i2*0.1)))
}



#Set the computation
numcores<-detectCores()
clust<-makeCluster(numcores)
clusterExport(clust,c("myoverlap","birdnet","build.util","resources","obs","sim.comp","newset",
                      "cors","mid.cor","mid.cul","culs","p1.max","v1.max","nsp","birdlist"))
registerDoParallel(clust)




paramdat2<-foreach(k1=1:nrow(newset),.combine=rbind,
                   .packages=c("msm","pROC","igraph","reshape2","tidyverse"))%dopar%
                   {
                     
                     p1<-newset[k1,1]
                     Vb<-newset[k1,2]
                     a.vec<-unlist(newset[k1,3:7])
                     
                     res<-sim.comp(resources,build.util(p1,Vb,nsp),a.vec)
                     return(c(p1,Vb,a.vec,res))
                   }

stopCluster(clust)
write.csv(paramdat2,"secondparam.csv")


#Third filter
param.main<-paramdat2
param.main<-param.main[,-1]

# Choose the rows where all AUCs are above 0.7
good<-apply(param.main[,8:12],1,min)
good.ind<-which(good>0.7)
paramset<-param.main[good.ind,1:7]

#Add new parameter values in the neighborhood of the filtered values.
newset<-paramset
for(i1 in 1:30){
  newset<-rbind(newset,apply(paramset,c(1,2),function(i2) rnorm(1,i2,i2*0.1)))
}



#Set the computation
numcores<-detectCores()
clust<-makeCluster(numcores)
clusterExport(clust,c("myoverlap","birdnet","build.util","resources","obs","sim.comp","newset",
                      "cors","mid.cor","mid.cul","culs","p1.max","v1.max","nsp","birdlist"))
registerDoParallel(clust)




paramdat3<-foreach(k1=1:nrow(newset),.combine=rbind,
                   .packages=c("msm","pROC","igraph","reshape2","tidyverse"))%dopar%
                   {
                     
                     p1<-newset[k1,1]
                     Vb<-newset[k1,2]
                     a.vec<-unlist(newset[k1,3:7])
                     
                     res<-sim.comp(resources,build.util(p1,Vb,nsp),a.vec)
                     return(c(p1,Vb,a.vec,res))
                   }

stopCluster(clust)
write.csv(paramdat3,"thirdparam.csv")



