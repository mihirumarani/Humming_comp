library(igraph)
library(reshape)
library(moments)
library(ggplot2)
library(pROC)
library(PresenceAbsence)
library(tidyverse)
library(knitr)
library(kableExtra)
library(nnet)



#############
#Function to get equilibrium sp. abundances
get.abun=function(alphamat,k){
  
  pops=rep(0,length(birdlist))
  
  sp1=which(k>0)
  
  if(length(sp1)==1){
    
    pops[sp1]=k[sp1]
    
  }else{
    
    set=sp1
    
    while(length(set)>1){

      ans=solve(alphamat[set,set],k[set])
      
      if(sum(ans<=0)==0){
        
        pops[set]=ans
        
        break
        
      }else{
        set=set[which(ans>0)]
      }
    }
  }
  
  return(pops)
}
  


TSS = function(cmx){
  require(PresenceAbsence)
  sensitivity(cmx, st.dev=F)+specificity(cmx, st.dev=F)-1
}


####################################################################################


#Set up data

interact<-read.csv("HummingbirdInteractions.csv")%>%as_tibble()
interact<-interact%>%select(month,year,ele,Hummingbird,Iplant_Double)
colnames(interact)[c(3,4,5)]<-c("elevation","bird","plant")
bins<-c(1300,1500,1700,1900,2100,2300,2500)
interact$elevation<-.bincode(interact$elevation,bins,FALSE)
interact<-na.omit(interact)
interact$Y<-1

birdlist=interact%>%
  group_by(bird)%>%
  summarize(Y=sum(Y))%>%
  ungroup()%>%
  filter(Y>100)%>%
  pull(bird)

interact=interact%>%
          filter(bird %in% birdlist)


plantlist<-unique(interact$plant)



fl.abun<-read.csv("FlowerTransectClean.csv")%>%as_tibble()
abund<-fl.abun%>%select(month,year,Elevation.Begin,Iplant_Double,Total_Flowers)
colnames(abund)[c(3,4,5)]<-c("elevation","plant","abundance")
abund$elevation<-as.factor(abund$elevation)
levels(abund$elevation)<-c(1,2,3,4,5,6)
abund<-abund%>%
        group_by(month,year,elevation,plant)%>%
        summarise(abundance=sum(abundance))%>%
        ungroup()%>%
        filter(plant %in% plantlist)%>%
        drop_na()%>%
        mutate(abundance=log(abundance+1))%>%
        mutate(elevation=as.factor(elevation))


totaldat=interact%>%
  mutate(elevation=as.factor(elevation))%>%
  left_join(abund)%>%
  mutate(time=(100*year)+month)%>%
  select(-c(month,year))
 # add.plant<-setdiff(plantlist,flower.trait$plant)
 # add.plant<-cbind(add.plant,runif(length(add.plant),min(flower.trait$corolla),max(flower.trait$corolla)))
 # flower.trait[(nrow(flower.trait)+1):(nrow(flower.trait)+nrow(add.plant)),]<-add.plant     
 # flower.trait$corolla<-as.numeric(flower.trait$corolla)
 # flower.trait[,2]<-round(flower.trait[,2],2)


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
birdwing=birdwing[-10,]

birdforage<-read.csv("HummingbirdSpecies.csv")
birdforage<-birdforage[,c(1,3)]
colnames(birdforage)<-c("bird","role")
birdforage<-birdforage[match(birdlist,birdforage$bird),]






miss.p<-unlist(unique(totaldat[which(is.na(totaldat$abundance)),3]))
for(i in miss.p){
  ind<-which(totaldat$plant==i)
  plant<-totaldat[ind,]
  if(!all(is.na(plant$abundance))){
    miss<-which(is.na(totaldat[ind,5]))
    totaldat[ind[miss],5]<-mean(plant$abundance,na.rm=TRUE)
  }
}

totaldat<-totaldat%>%drop_na()

birdlist=totaldat%>%select(bird)%>%unique()%>%pull(bird)
plantlist=totaldat%>%select(plant)%>%unique()%>%arrange(desc(plant))%>%pull(plant)
#######################################################################################
#Set up a test data from years 2015 and 2016

test.dat<-totaldat%>%slice(c(grep("2015",time),grep("2016",time)))%>%
  select(time,elevation,bird,Y)%>%
  group_by(time,elevation,bird)%>%
  summarize(Y=mean(Y))%>%ungroup()


testdat<-as_tibble(test.dat)%>%
        complete(bird,nesting(time,elevation),fill=list(Y=0))
test.bird<-unique(testdat$bird)
samples.test<-unique(testdat[,2:3])

int.sum=totaldat%>%slice(c(grep("2015",time),grep("2016",time)))%>%
                         group_by(plant,bird,time,elevation)%>%
                         summarize(Y=sum(Y),abundance=mean(abundance))%>%ungroup()%>%
                        mutate(Y=Y/abundance)  


##################################################################


######################################################################################
#Set up a training data to estimate elevation frequencies and  utility curves
totaldat1<-totaldat%>%
              slice(c(grep("2013",time),grep("2014",time)))
  
  
#Bootstrap the samples to estimate the parameters

reps=1:100

results=NULL



for(r in reps){

  totaldat.tr=totaldat1%>%
              group_by(bird)%>%
              slice(sample((1:n()),round(n()*0.9)))%>%
              ungroup()
  
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
  
  
  train.el$prob[train.el$prob!=0]<-1
  ##################################################################################
  #Estimate alpha and carrying capacities
  
  # use training data (years 2013 and 2014)
  
  int.sum1<-totaldat.tr%>%
            group_by(plant,bird,time,elevation)%>%
            summarize(Y=sum(Y),abundance=mean(abundance))%>%
            ungroup()%>%
            mutate(Y=Y/abundance)
  

  # corol<-unique(totaldat$corolla)
  #cors<-length(corol)
  
  utility.s<-int.sum1%>%
    select(plant,bird,time,elevation,Y)%>%
    group_by(plant,bird)%>%
    summarize(Y=sum(Y))%>%ungroup()%>%
    right_join(
      tibble(expand.grid(bird=birdlist,
                         plant=plantlist))
    )%>%
    replace_na(list(Y=0))

  
  util.raw=matrix(nrow=length(plantlist),ncol=length(birdlist))
  

  for(i in 1:length(birdlist)){
    util.raw[,i]=utility.s%>%filter(bird==birdlist[i])%>%arrange(desc(plant))%>%pull(Y)
    
  }
  
  util.raw=as.data.frame(apply(util.raw,2,function(x) return(x/sum(x))))
  
  names(util.raw)=birdlist
  
  util.raw=util.raw%>%
          mutate(plant=plantlist)
  
  
  birdnet=matrix(NA,length(birdlist),length(birdlist))
  
  for(i in 1:length(birdlist)){
    for(j in 1:length(birdlist)){
      a=util.raw[,i]
      b=util.raw[,j]
      birdnet[i,j]=sum(sapply(1:nrow(util.raw),function(x) return(
        ifelse((a[x]<=b[x]),a[x],b[x]))))
    }
  }
  
  
  alphamat=birdnet
  colnames(alphamat)<-birdlist
  

##############################################################################################################  
 #Estimate carrying capacities for all sampling points in 2015 and 2016
  samples<-int.sum%>%select(time,elevation)%>%unique()
  
  kdat=NULL
  
 for(j in 1:nrow(samples)){

   
   for(i in 1:length(birdlist)){
     
     dat=int.sum%>%
       filter(bird==birdlist[i],
              time==unlist(samples[j,1]),
              elevation==unlist(samples[j,2]))
     
     if(nrow(dat)>0){
       
       dat2=util.raw%>%
         filter(plant %in% dat$plant)%>%
         select(birdlist[i],plant)
       names(dat2)[1]="weights"
       
       
       kdat=bind_rows(kdat,
                      tibble(
                        time=unlist(samples[j,1]),
                        elevation=unlist(samples[j,2]),
                        bird=birdlist[i],
                        k=dat%>%
                          inner_join(dat2)%>%
                          mutate(x=sum(abundance*weights))%>%
                          pull(x)))
       

      }else{
        kdat=bind_rows(kdat,
                       tibble(
                         time=unlist(samples[j,1]),
                         elevation=unlist(samples[j,2]),
                         bird=birdlist[i],
                         k=0))
     }
     
     }

 }
        
  kdat=kdat%>%unique()
 
  
  colnames(kdat)<-c("time","elevation","bird","K")
  
 ##########################################
  #Elevation-wise frequency of species
  eldat<-as.data.frame(cast(train.el,bird~elevation))
  

  
####################################################################  
#Create prediction probabilities under different assembly processes
  predictions=kdat%>%select(time,elevation,bird)%>%
              inner_join(birdwing)
  
  
  #Model 1a: Species presence/absence is explained only by elevation
  predictions=predictions%>%
                inner_join(train.el)%>%
              rename(el1=prob)
  
  #Model 1b: Species presence/absence is explained only regression between
  #wing size and elevation
  # wing.el1=train.el%>%
  #           inner_join(birdwing)%>%
  #           filter(prob==1)%>%
  #           mutate(elevation=as.factor(elevation))
  # el.glm=glm(wing.el1$elevation~wing.el$wing, 
  #            family=poisson) #Model performance is terrible!
  # 
  # #Try regression between max elevation and wing
  # wing.el2=train.el%>%
  #       inner_join(birdwing)%>%
  #       filter(prob==1)%>%
  #       group_by(bird,wing)%>%
  #       slice_max(elevation)%>%
  #       ungroup()
  # 
  # el.glm2=glm(wing.el2$elevation~wing.el2$wing,family='poisson')
  
#Model 2: Species presence/absence is explained by resources
  predictions=predictions%>%
                inner_join(kdat)%>%
                mutate(K=K/max(K))%>%
                mutate(el_K=K*el1)
  
#Model 3:  LV model where the prob. of occurrence depends on equilibrium 
  #population size of competing species.
  lv=NULL
  lv_dash=NULL
  
  for(i in 1:nrow(samples.test)){
    
    ks=kdat%>%filter(time==samples.test$time[i],
                     elevation==samples.test$elevation[i])%>%
      pull(K)
    
    ksdash=ks*(predictions%>%filter(time==predictions$time[i],
                                   elevation==predictions$elevation[i])%>%
                 pull(el1))
    
    pop=get.abun(alphamat,ks)
    popdash=get.abun(alphamat,ksdash)
    
    lv=bind_rows(lv,
                 tibble(time=samples.test$time[i],
                 elevation=samples.test$elevation[i],
                 bird=birdlist,
                 lv=pop))
    
    lv_dash=bind_rows(lv_dash,
                      tibble(
                        time=samples.test$time[i],
                 elevation=samples.test$elevation[i],
                 bird=birdlist,
                 el_lv=popdash))
  }
  
  
predictions=predictions%>%
              inner_join(lv)%>%
              inner_join(lv_dash)

#Compare observations
predictions=predictions%>%
  inner_join(testdat)



#####################################################################
#Model performances

predictions=as.data.frame(predictions)


for(i in 1:length(birdlist)){
  for(j in 1:5){
    
    dat=predictions%>%
        filter(bird==birdlist[i])
    
    dat=dat[,c(4+j,10)]
    
    AUC=as.numeric(roc(dat[,2],dat[,1])$auc)
    
    tssdat=NULL
    
    for(k in seq(0.01,1,0.01)){
      
      cmd=cmx(data.frame(ID=1:nrow(dat),obs=dat[,2], mod=dat[,1]),threshold=k)
      
      tssdat=bind_rows(tssdat,
                       tibble(threshold=k,
                              TSS=TSS(cmd)))
      
    }
    
    results=bind_rows(results,
                      tibble(
                        reps=r,
                        bird=birdlist[i],
                        model=names(predictions[4+j]),
                        auc=AUC,
                        threshold=tssdat$threshold,
                        tss=tssdat$TSS)
                      )
    
  }
}
}

  saveRDS(results,"humming_res.rds")
  
  #Plot the results
  
  theme_set(theme_bw())
  theme_update(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = 'orange')
  )
  
  best.aucs=results%>%
            select(reps,bird,model,auc)%>%
            unique()
  
  relabel.model = as_labeller(c(el1 = "Elevation", 
                                  K = "Resources", 
                                  el_K = "Resources+Elevation",
                                  lv="Competition",
                                  el_lv="Competition+Elevation"))
  
  best.aucs%>%
    rowwise()%>%
    mutate(hypothesis=factor(relabel.model(model)))%>%
    ggplot(aes(hypothesis,auc))+
    geom_boxplot()+
    theme(axis.text.x = element_text(face='bold',size=10,angle = 90))+
    facet_wrap(vars(bird))
  
  
  #####################################################################################
  #ABC model
  #DO NOT RUN
  
  #Start the cluster and split the loop on #cores
  library(parallel)
  library(foreach)
  library(iterators)
  library(doParallel)
  
  numcores<-detectCores()
  clust<-makeCluster(3)
  clusterExport(clust,c("eldat","alphamat","kdat","samples.test","testdat","test.bird"))
  registerDoParallel(clust)
  
  
  dat<-foreach(j1=1:40000,.combine=rbind,.packages=c("tidyverse","PresenceAbsence","pROC"))%dopar%
  {
    
  abd<-c()
  betas<-runif(9,0,1)
  
  for(i in 1:nrow(samples.test)){
    begin<-which(eldat[,(1+as.numeric(samples.test[i,2]))]==1)
    alphamat1<-alphamat[begin,begin]
    Klist<-kdat%>%
      filter(time==samples.test[i,1],elevation==samples.test[i,2])%>%
      pull(K)
    
    beta1<-betas[begin]
    
    pops<-solve(alphamat1,((Klist)/beta1))
    pops<-pops+abs(min(pops))
    pops<-(pops/max(pops))+0.1
    
    abd1<-rep(0,9)
    abd1[begin]<-pops
    
    abd<-rbind(abd,cbind(samples.test[rep(i,9),],eldat$bird,abd1))
  }
  abd<-abd[order(abd$`eldat$bird`,abd$month,abd$year,abd$elevation),]
  aucdat<-c()
  for(j in 1:9){
    obs<-testdat$Y[testdat$bird==test.bird[j]]
    preds<-abd$abd1[abd$`eldat$bird`==test.bird[j]]
    aucdat<-c(aucdat,roc(obs,preds)$auc)
  }
  return(c(betas,aucdat))
  
  }
  stopCluster(clust)
  
#write.csv(dat,"hum_beta.csv")  
master1<-read.csv("hum_beta.csv")
master1<-master1[,-1]
master1<-dat  
m1<-master1[,c(1,10)]
m1<-m1[m1[,2]>0.7,]
par(mfrow=c(3,3))
hist(m1[,1],main=test.bird[1],xlab=expression(beta),xlim=c(0,1))
max.ind<-apply(master1[,10:18],2,which.max)
for(i in 2:9){
  mi<-master1[,c(i,9+i)]
  mi<-mi[mi[,2]>0.7,]
  if(length(mi)>3){
  hist(mi[,1],main=test.bird[i],xlab=expression(beta),xlim=c(0,1))
  }else{hist(runif(1000,0,1),main=paste(test.bird[i],"NS!"),xlim=c(0,1))}
}


max.beta<-vector(length=9)
for(i in 1:9){
  master.i<-master1[,c(i,i+9)]
  master.i<-master.i[master.i[,2]>0.6,]
  max.beta[i]<-master1[max.ind[i],i]
  
}
  
  
  
  
  