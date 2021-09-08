library(igraph)
library(reshape)
library(moments)
library(ggplot2)
library(pROC)
library(PresenceAbsence)
library(tidyverse)
library(knitr)
library(kableExtra)


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


####################################################################################

#Set up data

interact<-read.csv("HummingbirdInteractions.csv")
interact<-interact%>%select(month,year,ele,Hummingbird,Iplant_Double)
colnames(interact)[c(3,4,5)]<-c("elevation","bird","plant")
bins<-c(1300,1500,1700,1900,2100,2300,2500)
interact$elevation<-.bincode(interact$elevation,bins,FALSE)
interact<-na.omit(interact)
interact$Y<-1

bird.obs<-interact%>%group_by(bird)%>%summarise(Y=sum(Y))%>%ungroup()%>%filter(Y>100)

birdlist<-as.vector(unlist((interact%>%group_by(bird)%>%summarise(Y=sum(Y))%>%
                              ungroup()%>%filter(Y>100)%>%select(bird))))

birdlist<-birdlist[order(birdlist)]

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
colnames(birdforage)<-c("bird","role")
birdforage<-birdforage[match(birdlist,birdforage$bird),]


interact$corolla<-as.numeric(flower.trait$corolla[match(interact$plant,flower.trait$plant)])
totaldat<-merge(interact,abund,all.x=TRUE)



miss.p<-unique(totaldat[which(is.na(totaldat$abundance)),4])
for(i in miss.p){
  ind<-which(totaldat$plant==i)
  plant<-totaldat[ind,]
  if(!all(is.na(plant$abundance))){
      miss<-which(is.na(totaldat[ind,8]))
      totaldat[ind[miss],8]<-mean(plant$abundance,na.rm=TRUE)
  }
}

totaldat<-na.omit(totaldat)


#######################################################################################
#Set up a test data from years 2015 and 2016

test.dat<-totaldat%>%filter(year>2014)%>%
  select(month,year,elevation,bird,Y)%>%
  group_by(month,year,elevation,bird)%>%
  summarize(Y=mean(Y))%>%ungroup()


testdat<-as.data.frame(as_tibble(test.dat)%>%complete(bird,nesting(month,year,elevation),fill=list(Y=0)))
test.bird<-unique(testdat$bird)
samples.test<-unique(testdat[,2:4])

int.sum<-as.data.frame(totaldat%>%filter(year>2014)%>%
                         group_by(plant,bird,month,year,elevation)%>%
                         summarize(Y=sum(Y), corolla=mean(corolla),abundance=mean(abundance))%>%ungroup())

int.sum$Y<-int.sum$Y/int.sum$abundance

##################################################################


######################################################################################
#Set up a training data to estimate elevation frequencies and  utility curves
totaldat.tr<-totaldat[totaldat$year<2015,]

#Bootstrap the interaction data to get distributions of alphas and Ks

finset<-vector(mode="list",length=length(test.bird))

for(i1 in 1:200){
  
  boots<-sample(1:nrow(totaldat.tr),nrow(totaldat.tr),replace=TRUE)
  totaldat.b<-totaldat.tr[boots,]

  
#1. Probability of occurnace is predicted by elevations only. Create the predictions based
#on the observations from 2013 and 2014. 
train.el<-totaldat.tr
train.el$Y<-train.el$Y/train.el$abundance

train.el<-as.data.frame(train.el%>%filter(bird%in%test.bird)%>%
                          group_by(bird,elevation)%>%
                          summarize(Y=sum(Y))%>%ungroup())

train.el<-as.data.frame(as_tibble(train.el)%>%complete(bird,elevation,fill=list(Y=0)))
train.el<-as.data.frame(cast(train.el,bird~elevation))
train.el[,-1]<-train.el[,-1]/rowSums(train.el[,-1])
train.el<-melt(train.el)
train.el<-train.el[order(train.el$bird,train.el$variable),]
colnames(train.el)<-c("bird","elevation","prob")

#############################################
#Plot elevation-wise observations for species
train.el2<-train.el
train.el2$elevation<-as.numeric(train.el2$elevation)
ggplot(train.el2,aes(x=elevation,y=prob))+geom_point()+geom_line()+facet_wrap(vars(bird),ncol=3)+
  ylab("Frequency of observations")+
  theme(axis.title= element_text(face = "bold", size = 12),
        strip.text = element_text(face = "bold", size = 11))
###############################################
  results<-merge(testdat,train.el,all.x=TRUE)
  colnames(results)[ncol(results)]<- "el"
  results<-results[order(results$bird,results$month,results$year,results$elevation),]
  
##################################################################################
#Estimate alpha and carrying capacities
  
 # use training data (years 2013 and 2014)
  
  int.sum1<-as.data.frame(totaldat.b%>%group_by(plant,bird,month,year,elevation)%>%
          summarize(Y=sum(Y), corolla=mean(corolla),abundance=mean(abundance))%>%ungroup())
  
  int.sum1$Y<-int.sum1$Y/int.sum1$abundance
  
 # corol<-unique(totaldat$corolla)
  #cors<-length(corol)
  
  utility.s<-int.sum1%>%
    select(corolla,bird,month,year,elevation,Y)%>%
    group_by(corolla,bird)%>%
    summarize(Y=sum(Y))%>%ungroup()
  utility.s<-as.data.frame(utility.s)
  
  utility<-matrix(0,500,length(birdlist))
  for(i in 1:length(birdlist)){
    dat<-utility.s[utility.s$bird==birdlist[i],c(1,3)]
    #utility[1:cors,i]<-corol
    if(nrow(dat)>1){
      utility[1:nrow(dat),i]<-dat[,1]
      utility[(nrow(dat)+1):500,i]<-sample(dat[,1],(500-nrow(dat)),replace=TRUE,prob=dat[,2])
    }else{utility[,i]<-rep(dat[1,1],500)}
  }
  
##########################################  
#plot utility curves
  
  #Plot histgrams of utility curves for all species
  df.u<-as.data.frame(cbind(rep(birdlist,each=500),as.vector(utility)))
  colnames(df.u)<-c("species","corolla")
  df.u$species<-as.factor(df.u$species)
  df.u$corolla<-as.numeric(df.u$corolla)
  ggplot(df.u,aes(x=corolla))+geom_histogram(aes(y=..density..), position="identity",binwidth = 2)+
    geom_density(color="red",alpha=0.6)+
    facet_wrap(vars(species),ncol=3,scales="free_y",labeller = label_wrap_gen())+
    ylab("Probability")+xlab("Corolla length")+
    theme(axis.title= element_text(face = "bold", size = 12),
          strip.text = element_text(face = "bold", size = 11))
  
  #plot the correlation
  
  alphamat<-birdnet(utility)$mat
  alphamat<-alphamat+0.05
  colnames(alphamat)<-birdlist
  diag(alphamat)<-1
  
  
  #Estimate carrying capacities for all sampling points in 2015 and 2016
  
  samples<-unique(int.sum[,3:5])
  kdat<-expand.grid(birdlist,1:nrow(samples))
  kdat<-cbind(samples[kdat[,2],],kdat[,1])
  kdat$K<-NA

  for(i in 1:nrow(kdat)){
    
    dat1<-as.data.frame(int.sum%>%filter(month==kdat[i,1],year==kdat[i,2],elevation==kdat[i,3])%>%
                          select(corolla,abundance)%>%group_by(corolla)%>%
                          summarize(abundance=mean(abundance))%>%ungroup())
      
    if(nrow(dat1)>=1){
      dat2<-as.data.frame(table(utility[,which(birdlist==kdat[i,4])]))
      ans<-merge(dat1,dat2,by.x="corolla",by.y="Var1",all.x=TRUE)
      #if(!all(is.na(ans$Freq))){
        kdat[i,5]<-sum(ans$abundance*ans$Freq,na.rm=TRUE)#/sum(ans$Freq,na.rm=TRUE)
      #}else{kdat[i,5]<-0}
    }else{
      kdat[i,5]<-0
    }
  }
  
  
  colnames(kdat)<-c("month","year","elevation","bird","K")
  


# model 2. Probability of occurance depends on elevation preference and the amount of resources available
  kdat.train<-kdat%>%filter(bird%in%test.bird)
  results<-merge(results,kdat.train,all.x=TRUE)
  results<-results[order(results$bird,results$month,results$year,results$elevation),]
  
  colnames(results)[ncol(results)]<-"K"

  results$K1<-results$K
  
  results$K1[which(results$el==0)]<-0
  
  colnames(results)[ncol(results)]<-"K+el"
  
  
#3. LV model where the prob. of occurance depends on equilibrium population size of competing species.
  
  eldat<-as.data.frame(cast(train.el,bird~elevation))
  abd2<-c()
  abd<-c()
  
  for(i in 1:nrow(samples.test)){
    el<-samples.test[i,3]
    begin.list<-eldat$bird[which(eldat[,(el+1)]!=0)]
    begin.list<-begin.list[order(begin.list)]
    ind<-match(begin.list,birdlist)
    alphamat1<-alphamat[ind,ind]
    Klist<-kdat.train%>%
      filter(month==samples.test[i,1],year==samples.test[i,2],elevation==el,bird%in%begin.list)%>%
      select(K)
    Klist2<-kdat.train%>%
      filter(month==samples.test[i,1],year==samples.test[i,2],elevation==el)%>%select(K)
    
    pops<-solve(alphamat1,unlist(Klist))
    pops[pops<1]<-0
    pop1<-rep(0,length(eldat$bird))
    pop1[ind]<-pops
    #pops<-pops+abs(min(pops))
    #pops<-(pops/max(pops))+0.1
    
    pops2<-solve(alphamat,unlist(Klist2))
    pops2[pops2<1]<-0
    #pop2<-rep(0,length(eldat$bird))
    #pop2[ind]<-pops2
    #pops2<-pops2+abs(min(pops2))
    #pops2<-(pops2/max(pops2))+0.1
    
    abd<-rbind(abd,cbind(samples.test[rep(i,nrow(eldat)),],eldat$bird,pops2))
    #abd2<-rbind(abd2,cbind(samples.test[rep(i,length(begin.list)),],begin.list,pops))
    abd2<-rbind(abd2,cbind(samples.test[rep(i,nrow(eldat)),],eldat$bird,pop1))
    
  }
  

  colnames(abd)<-c("month","year","elevation","bird","pop")
  colnames(abd2)<-c("month","year","elevation","bird","pop")
  #abd2<-as.data.frame(abd2%>%complete(bird,nesting(month,year,elevation),fill=list(pop=0)))
  
  
  results<-merge(results,abd,all.x=TRUE)
  colnames(results)[ncol(results)]<-"LV"
  results<-merge(results,abd2,all.x=TRUE)
  colnames(results)[ncol(results)]<-"LV+el"
  results<-results[order(results$bird,results$month,results$year,results$elevation),]
  

#Evaluate models 

no_models<-ncol(results)-5
aucdat<-matrix(nrow=length(unique(results$bird)),ncol=no_models)

for(j in 1:length(test.bird)){
  for(k in 1:no_models){
    
    res<-results[which(results$bird==test.bird[j]),c(5,k+5)]
    aucdat[j,k]<-roc(res[,1],res[,2])$auc
    
  }
}

for(i2 in 1:length(test.bird)){
  finset[[i2]]<-rbind(finset[[i2]],aucdat[i2,])
}

}


findat<-c()


for(i in 1:length(test.bird)){
  a<-finset[[i]]
  colnames(a)<-colnames(results)[6:ncol(results)]
  a<-melt(a)[,-1]
  findat<-rbind(findat,cbind(rep(test.bird[i],nrow(a)),a))
}

colnames(findat)<-c("bird","model","auc")
findat$bird<-as.factor(findat$bird)
findat$model<-as.factor(findat$model)

#OR load the data:
findat<-read.table("modelp_15sp.csv")

best.m<-c()

for(i in 1:length(test.bird)){
  dat<-findat[findat$bird==test.bird[i],2:3]
  means<-aggregate(auc~model,dat,mean)
  means<-means[order(means$auc,decreasing=T),]
  ao<-aov(auc~model,dat)
  tuk<-as.data.frame(TukeyHSD(ao,conf.level = 0.99)$model)
  tuk1<-cbind(matrix(unlist(strsplit(rownames(tuk),split="-")),ncol=2,byrow=T),tuk)
  
  tukset<-tuk1[which(tuk1[,1]==means[1,1] | tuk1[,2]==means[1,1]),]
  mods<-(means[1,])
  j<-2
  while(j<5){
    pval<-tukset[which(tukset[,1]==means[j,1] | tukset[,2]==means[j,1]),6]
    if(pval<0.005){
      j<-5
    }else{mods<-rbind(mods,means[j,])
        j<-j+1  
    }
  }
  
  best.m<-rbind(best.m,cbind(rep(test.bird[i],nrow(mods)),mods))
}

best.m$auc<-round(best.m$auc,3)
best.m$model<-as.character(best.m$model)
best.m$model[best.m$model=="LV"]<-"Comp"
best.m$model[best.m$model=="LV+el"]<-"Comp + El"




best.m$bill<-birdbill[match(best.m$`rep(test.bird[i], nrow(mods))`,birdbill[,1]),2]
best.m$fora<-birdforage[match(best.m$`rep(test.bird[i], nrow(mods))`,birdforage[,1]),2]

best.m1 <-aggregate(best.m[,2:3],by=list(best.m[,1],best.m$bill,best.m$fora), paste)
new1<-c()
for(i in 1:nrow(best.m1)){
  b1<-unlist(best.m1[i,4])
  b2<-unlist(best.m1[i,5])
  d1<-c()
  d2<-c()
  for(j in 1:length(b1)){
    d1<-paste(d1,b1[[j]])
    d2<-paste(d2,b2[[j]])
  }
  new1<-rbind(new1,c(d1,d2))
}
best.m1<-cbind(best.m1[1:3],new1)
colnames(best.m1)<-c("Species","Culmen Length","Foraging Behavior","Model","Mean AUC")

write.csv(best.m,"bestm_15sp.csv")

best.m1%>%kbl()%>%kable_classic(full_width=F)%>%column_spec(4:5,width="2cm")%>%save_kable(file="full_model.html")

#Find the performance of competition models:
bes.comp<-findat[findat$model%in%c("LV","LV+el"),]
comp.means<- aggregate(bes.comp$auc,by=list(bes.comp$bird,bes.comp$model),mean)
comp05<-aggregate(bes.comp$auc,by=list(bes.comp$bird,bes.comp$model),function(x) quantile(x,0.05))
comp95<-aggregate(bes.comp$auc,by=list(bes.comp$bird,bes.comp$model),function(x) quantile(x,0.95))
com.per<-cbind(comp.means,comp05$x,comp95$x)
colnames(com.per)<-c("Species","Model","Mean AUC", "5th Percentile", "95th Percentile")

comp.per1<-com.per[com.per$Model=="LV+el",-2]
comp.per1[,-1]<-round(comp.per1[,-1],3)
comp.per1$new<-paste("(",comp.per1$`5th Percentile`,",",comp.per1$`95th Percentile`,")",sep="")
colnames(comp.per1)[5]<-"95% Confidence Interval"
comp.per2<-com.per[com.per$Model=="LV",-2]

comp.per1[,c(1,2,5)]%>%kbl(caption="Performance of the competition model")%>%kable_classic(full_width=F)%>%
  row_spec(0,bold=TRUE)%>%save_kable(file="png")

ggplot(bes.comp,aes(x=model,y=auc))+ geom_boxplot() + 
  facet_wrap(vars(bird),ncol=3)

ggplot(findat,aes(x=model,y=auc))+ geom_boxplot() + 
facet_wrap(vars(bird),ncol=3)

write.csv(best.m,"bestm_9sp.csv")
write.table(findat,"modelp_9sp.csv")
