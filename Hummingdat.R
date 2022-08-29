library(igraph)
library(reshape)
library(moments)
library(ggplot2)
library(pROC)
library(PresenceAbsence)
library(tidyverse)


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
totaldat<-interact

#Get the utility curves of all hummingbird species as a function of flower trait values
#These utility curves are based on pooled data over years and locations

int.sum<-totaldat%>%group_by(plant, bird,month,year,elevation)%>%
  summarize(Y=sum(Y), corolla=mean(corolla))%>%ungroup()


int.sum<-merge(int.sum,abund,all.x=TRUE)


#Fill in the missing values
miss<-which(is.na(int.sum$abundance))
for(i in 1:length(miss)){
  plant<-int.sum[which(int.sum$plant==int.sum[i,1]),]
  if(!all(is.na(plant$abundance))){
    int.sum$abundance[i]<-mean(plant$abundance,na.rm=TRUE)
  }
}

int.sum<-na.omit(int.sum)
int.sum$Y<-int.sum$Y/int.sum$abundance

  



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



#EStimate utility curves
utilset<-list()
dat<-int.sum%>%select(elevation,bird,Y,corolla)

for(i in 1:length(birdlist)){
  ut.bird<-dat%>%filter(bird==birdlist[i])

  mat1<-matrix(0,500,length(unique(int.sum$elevation)))
  
  for(i1 in 1:length(unique(int.sum$elevation))){
    ut<-ut.bird%>%filter(elevation==i1)%>%
                  group_by(bird,corolla)%>%
                  summarize(Y=sum(Y))%>%ungroup()%>%
                  select(corolla,Y)
    ut<-as.data.frame(ut)
    
      if(nrow(ut)>1){
        mat1[1:nrow(ut),i1]<-ut[,1]
        mat1[(nrow(ut)+1):500,i1]<-sample(ut[,1],(500-nrow(ut)),replace=TRUE,prob=ut[,2])
        }else{mat1[,i1]<-rep(flower.trait$corolla,length.out=500)}
  }
  utilset<-c(utilset,list(mat1))
}



utility.s<-int.sum1%>%
  select(corolla,bird,month,year,elevation,Y)%>%
  group_by(corolla,bird)%>%
  summarize(Y=sum(Y))%>%ungroup()
utility.s<-as.data.frame(utility.s)

utility<-matrix(0,500,length(birdlist))
for(i in 1:length(birdlist)){
  dat<-utility.s[utility.s$bird==birdlist[i],c(1,3)]
  if(nrow(dat)>1){
    utility[1:nrow(dat),i]<-dat[,1]
    utility[(nrow(dat)+1):500,i]<-sample(dat[,1],(500-nrow(dat)),replace=TRUE,prob=dat[,2])
  }else{utility[,i]<-rep(dat[1,1],500)}
}

#Assume a Gaussian utility function of which mean and variances are estimated form empirical data
utility.g<-matrix(0,500,length(birdlist))
res.trait<-unique(int.sum$corolla)

for(i in 1:length(birdlist)){
  dat<-utility.s[utility.s$bird==birdlist[i],c(1,3)]
  
  if(nrow(dat)>1){
    mean<-sum(dat$corolla*dat$Y)/sum(dat$Y)
    var<-sqrt((sum(dat$Y)*sum((dat$corolla^2)*dat$Y)-(mean^2))/(sum(dat$Y)*(sum(dat$Y)-1)))
    utility.g[,i]<-c(res.trait,sample(res.trait,(500-length(res.trait)),replace=TRUE,
                                      prob=dnorm(res.trait,mean,var)))
  }else{
    utility.g[,i]<-c(res.trait,rep(dat[1,1],(500-length(res.trait))))
  }
}

#Build a competition matrix

alphamat<-birdnet(utility)$mat
alphamat<-alphamat+0.05
colnames(alphamat)<-birdlist
diag(alphamat)<-1

alphamat.g<-birdnet(utility.g)$mat
alphamat.g<-alphamat+0.05
colnames(alphamat.g)<-birdlist
diag(alphamat.g)<-1

#Plot histgrams of utility curves for all species
df.u<-as.data.frame(cbind(rep(birdlist,each=500),as.vector(utility)))
colnames(df.u)<-c("species","corolla")
df.u$species<-as.factor(df.u$species)
df.u$corolla<-as.numeric(df.u$corolla)
ggplot(df.u,aes(x=corolla))+geom_density()+facet_wrap(vars(species),ncol=3)+
  ylab("Probability")+xlab("Corolla length")+
  theme(axis.title= element_text(face = "bold", size = 12),
        strip.text = element_text(face = "bold", size = 11))

#plot the correlation between utility curves and bird traits.
ut.m<-apply(utility,2,mean)
ut.s<-apply(utility,2,sd)
df1<-as.data.frame(cbind(birdbill[,2],ut.m,ut.s))

ggplot(df1,aes(x=V1,y=ut.m))+geom_point()+geom_smooth(method="lm")+xlab("Culmen length")+
  ylab("Mean corolla length")
ggplot(df1,aes(x=V1,y=ut.s))+geom_point()+geom_smooth(method="lm")+xlab("Culmen length")+
  ylab("SD in corolla length")

bill.d<-dist(birdbill[,2])

df2<-c()
for(i in 1:nrow(birdbill)){
  d1<-abs(birdbill[,2]-birdbill[i,2])
  comp<-alphamat[i,]
  df2<-rbind(df2,cbind(rep(birdbill[i,1],nrow(birdbill)),d1,comp))
}
df2<-as.data.frame(df2)
colnames(df2)<-c("bird","trd","comp")
df2$bird<-as.factor(df2$bird)
df2$trd<-as.numeric(df2$trd)
df2$comp<-as.numeric(df2$comp)
df2$comp<--log(df2$comp)

ggplot(df2,aes(x=trd,y=comp))+geom_point()+geom_smooth(method="lm")+
  facet_wrap(vars(bird),ncol=3)+xlab(expression(paste("|",Delta,"z","|")))+
  ylab("-log(resource overlap)")+
  theme(axis.title= element_text(face = "bold", size = 12))


#Get the resource amounts

#For hypothesis 2a: Just get the total amount of resources at every sample point:This assumes that species have equal preferences
#for all species.
abund1<-int.sum%>%group_by(month,year,elevation)%>%summarize(abundance=sum(abundance))%>%ungroup()

#For hypothesis 2b: standardize the resource abundance with resource utility for each bird species, i.e. get the carrying capacity
samples<-unique(int.sum[,3:5])
kdat<-expand.grid(birdlist,1:nrow(samples))
kdat<-cbind(samples[kdat[,2],],kdat[,1])
kdat.g<-kdat
ks<-c()
kg<-c()

for(i in 1:nrow(samples)){
  for(j in 1:length(birdlist)){
    dat1<-int.sum[which(int.sum$month==samples[i,1] & int.sum$year==samples[i,2] & 
                          int.sum$elevation==samples[i,3] & int.sum$bird==birdlist[j]),7:8]
    if(nrow(dat1)>=1){
      dat2<-as.data.frame(table(utility[,j]))
      #dat2.g<-as.data.frame(table(utility.g[,j]))
      ans<-merge(dat1,dat2,by.x="corolla",by.y="Var1",all.x=TRUE)
     # ans.g<-merge(dat1,dat2.g,by.x="corolla",by.y="Var1",all.x=TRUE)
      ks<-c(ks,sum(ans$abundance*ans$Freq,na.rm=TRUE))#/sum(ans$Freq))
      #kg<-c(kg,sum(ans.g$abundance*ans.g$Freq,na.rm=TRUE)/sum(ans$Freq))
    }else{
      ks<-c(ks,0)
      #kg<-c(kg,0)
    }
  }
}

kdat$K<-ks
colnames(kdat)<-c("month","year","elevation","bird","K")

kdat$K[is.na(kdat$K)]<-0

kdat.g$K<-kg
colnames(kdat.g)<-c("month","year","elevation","bird","K")
kdat.g$K[is.na(kdat.g$K)]<-0

for(i in 1:length(birdlist)){
  ind1<-which(kdat$bird==birdlist[i])
  kdat$K[ind1]<-kdat$K[ind1]/max(kdat$K[ind1])
}


#######################################################################################
#Set up a test data from years 2015 and 2016
test.dat<-totaldat%>%filter(year>2014)%>%
  select(month,year,elevation,bird,Y)%>%
  group_by(month,year,elevation,bird)%>%
  summarize(Y=mean(Y))%>%ungroup()


testdat<-as.data.frame(as_tibble(test.dat)%>%complete(bird,nesting(month,year,elevation),fill=list(Y=0)))
test.bird<-unique(testdat$bird)
samples.test<-unique(testdat[,2:4])

results<-testdat
results<-results[order(results$bird,results$month,results$year,results$elevation),]


##################################################################
#Testing different models of occurances
#Model 1a:

train.el<-as.data.frame(totaldat%>%filter(year<2015,bird%in%test.bird)%>%
                          select(bird,month,year,elevation,Y)%>%group_by(bird,elevation)%>%
                          summarize(Y=sum(Y))%>%ungroup())
train.el<-as.data.frame(as_tibble(train.el)%>%complete(bird,elevation,fill=list(Y=0)))
train.el<-as.data.frame(cast(train.el,bird~elevation))
train.el[,-1]<-train.el[,-1]/rowSums(train.el[,-1])
train.el<-melt(train.el)
train.el<-train.el[order(train.el$bird,train.el$variable),]
colnames(train.el)<-c("bird","elevation","prob")

results<-merge(testdat,train.el,all.x=TRUE)
colnames(results)[6]<- "Model1a"
results<-results[order(results$bird,results$month,results$year,results$elevation),]

#model1b: Only wing_area/elevation explains the observed data:

#Create a linear regression model to estimate the relationship between elevation and wing-size.
#Use the model with residual error to create 100 sets of predictions of elevation points on the test data and create a probability distribution of getting each elevation point for a given value of wing-size.
#Two sub-methods: 1) Use continous elevation data 2) Use categorical elevation data

testdat1b<-merge(testdat,birdwing,all.x=TRUE)
testdat1b<-testdat1b[order(testdat1b$bird,testdat1b$month,testdat1b$year,testdat1b$elevation),]
test.wing<-data.frame(wing=testdat1b$wing)

interact1<-read.csv("HummingbirdInteractions.csv")
interact1<-interact1%>%select(month,year,ele,Hummingbird,Iplant_Double)
colnames(interact1)[c(3,4,5)]<-c("elevation","bird","plant")
interact1<-merge(interact1,birdwing,all.x=TRUE)

lm1<-lm(elevation~wing,data=interact1)
mean.res1<-summary(lm1)$sigma
sd.res1<-sd(abs(resid(lm1)))

dat.el<-merge(train.el,birdwing,all.x=TRUE)

lm2<-lm(as.numeric(elevation)~wing,data=dat.el,weights=prob)
mean.res2<-summary(lm2)$sigma
sd.res2<-sd(abs(resid(lm2)))


pred1<-matrix(nrow=nrow(test.wing),ncol=100)
pred2<-matrix(nrow=nrow(test.wing),ncol=100)

for(i in 1:100){
  preds1<-predict(lm1,newdata=test.wing)+rnorm(nrow(test.wing),mean.res1,sd.res1)

  pred1[,i]<-preds1
}

for(i in 1:100){
  preds1<-predict(lm1,newdata=test.wing)+rnorm(nrow(test.wing),mean.res1,sd.res1)
  bin1<-bins
  bin1[c(1,7)]<-c(min(bins[1],min(pred1)-1),max((max(pred1)+1),bins[7]))
  pred1[,i]<-.bincode(preds1,bin1,right=TRUE)
  
  preds2<-(predict(lm2,newdata=test.wing)+rnorm(nrow(test.wing),mean.res2,sd.res2))
  ra<-range(preds2)
  preds2<-preds2+1-ra[1]
  ratio<-5.5/diff(range(preds2))
  preds2<-round(1+(ratio*(preds2-1)))
  pred2[,i]<-preds2
}

prob1<-c()
prob2<-c()
for(j in 1:nrow(test.wing)){
  el<-testdat1b$elevation[j]
  prob1<-c(prob1,sum(pred1[j,]==el)/100)
  prob2<-c(prob2,sum(pred2[j,]==el)/100)
}

results$Model1b.1<-prob1
results$Model1b.2<-prob2


#Model2: Resource avaiability predicts the species presences

kdat.train<-kdat%>%filter(year>2014, bird%in%test.bird)
results<-merge(results,kdat.train,all.x=TRUE)
results<-results[order(results$bird,results$month,results$year,results$elevation),]

colnames(results)[ncol(results)]<-"Model2"

#Combine model 1 with model 2

results$Model1a.2<-results$Model2

results$Model1a.2[which(results$Model1a==0)]<-0
results$Model1b1.2<-results$Model2
results$Model1b1.2[which(results$Model1b.1==0)]<-0
results$Model1b2.2<-results$Model2
results$Model1b2.2[which(results$Model1b.2==0)]<-0



#Model3a: Species occupancy is driven by both the resource availability and competition (apply L-V equations)
#Test data from 2015 and 2016

abd<-matrix(ncol=5)


for(i in 1:nrow(samples.test)){
  Klist<-kdat.train%>%filter(month==samples.test[i,1],year==samples.test[i,2],elevation==samples.test[i,3])%>%
    select(K)
  pops<-solve(alphamat,unlist(Klist))
  pops[pops<0]<-0
  pops<-pops+abs(min(pops))
  pops<-pops/max(pops)
  abd<-rbind(abd,cbind(as.matrix(samples.test[rep(i,length(pops)),]),test.bird,pops))
}

abd<-as.data.frame(abd)
abd<-abd[-1,]
abd[,1]<-as.numeric(abd[,1])
abd[,2]<-as.numeric(abd[,2])
abd[,3]<-as.numeric(abd[,3])
abd[,5]<-as.numeric(abd[,5])

colnames(abd)<-c("month","year","elevation","bird","pop")

results<-merge(results,abd,all.x=TRUE)
colnames(results)[ncol(results)]<-"Model3"
results<-results[order(results$bird,results$month,results$year,results$elevation),]


#Model3a: Species occupancy is driven by elevation and LV equations.
#Test data from 2015 and 2015

eldat<-as.data.frame(cast(train.el,bird~elevation))
abd2<-matrix(ncol=5)

for(i in 1:nrow(samples.test)){
  el<-samples.test[i,3]
  begin.list<-eldat$bird[which(eldat[,(el+1)]!=0)]
  begin.list<-begin.list[order(begin.list)]
  ind<-match(begin.list,birdlist)
  alphamat1<-alphamat[ind,ind]
  Klist<-kdat.train%>%
    filter(month==samples.test[i,1],year==samples.test[i,2],elevation==el,bird%in%begin.list)%>%
    select(K)
  pops<-solve(alphamat1,unlist(Klist))
  pops<-pops+abs(min(pops))
  pops<-pops/max(pops)
  abd2<-rbind(abd2,cbind(as.matrix(samples.test[rep(i,length(begin.list)),]),begin.list,pops))
}

abd2<-as.data.frame(abd2)
abd2<-abd2[-1,]
abd2[,1]<-as.numeric(abd2[,1])
abd2[,2]<-as.numeric(abd2[,2])
abd2[,3]<-as.numeric(abd2[,3])
abd2[,5]<-as.numeric(abd2[,5])

colnames(abd2)<-c("month","year","elevation","bird","pop")
abd2<-as.data.frame(abd2%>%complete(bird,nesting(month,year,elevation),fill=list(pop=0)))


results<-merge(results,abd2,all.x=TRUE)
colnames(results)[ncol(results)]<-"Model13"
results<-results[order(results$bird,results$month,results$year,results$elevation),]

txt<-c("Elevation freq.","Wing size cont.","wing size disc.","Carrying capacity",
       "Elevation freq + Carrying Capacity","wing size cont.+ Carrying Capacity",
       "wing size disc. + Carrying Capacity","Lotka-Volterra model","Elevation freq.+LV model")
txt2<-c("Elevation","Elevation","Elevation","Resources","Elevation + Resources","Elevation + Resources","Elevation + Resources","Competition","Competition")
Model_info<-data.frame(colnames(results[6:14]),Description=txt)
kable(Model_info)

AUC_table<-data.frame(bird=test.bird,culmen=birdbill$culmen,role=birdforage$role)
no_models<-ncol(results)-5
aucdat<-matrix(nrow=nrow(AUC_table),ncol=no_models)

for(i in 1:length(test.bird)){
for(j in 1:no_models){
  
    res<-results[which(results$bird==test.bird[i]),c(5,j+5)]
    aucdat[i,j]<-roc(res[,1],res[,2])$auc
  }
}

colnames(aucdat)<-colnames(AUC_table[6,ncol(AUC_table)])
AUC_table<-cbind(AUC_table,aucdat)
kable(data.frame(AUC_table[],Total_data=bird.obs$Y))

best.auc<-apply(AUC_table[,-c(1:3)],1,max)
best.model<-apply(AUC_table[,-c(1:3)],1,which.max)


Winners<-data.frame(AUC_table[,1:3],Total_data=bird.obs$Y,Best_model=Model_info[best.model,2],AUC=best.auc)

Winners<-Winners[order(Winners$culmen,Winners$role),]


kable(Winners)

results1<-data.frame(ID=(1:nrow(results)),results[,-(1:4)])

thresh<-c()
for(i in 1:nrow(Model_info)){
  res1<-results1[,c(1,2,i+2)]
  thresh<-c(thresh,optimal.thresholds(DATA=res1,opt.methods = 3,na.rm=TRUE)[,2])
}



TSS = function(cmx){
  require(PresenceAbsence)
  sensitivity(cmx, st.dev=F)+specificity(cmx, st.dev=F)-1
}


pCC.dat<-matrix(nrow=length(test.bird),ncol=no_models)
rownames(pCC.dat)<-test.bird
colnames(pCC.dat)<-colnames(results1)[3:11]
sens.dat<-pCC.dat
spec.dat<-pCC.dat
kap.dat<-pCC.dat
tss.dat<-pCC.dat

t.dat<-rep(0.5,nrow(Model_info))

for(i in 1:length(test.bird)){
  for(j in 1:no_models){
    thr<-t.dat[j]
    res<-results1[which(results$bird==test.bird[i]),c(1,2,j+2)]
    res<-na.omit(res)
    cmdat<-cmx(DATA=res,threshold=thr)
    pCC.dat[i,j]<-pcc(cmdat,st.dev=F)
    sens.dat[i,j]<-sensitivity(cmdat,st.dev=F)
    spec.dat[i,j]<-specificity(cmdat,st.dev=F)
    kap.dat[i,j]<-kappa(cmdat,st.dev=F)
    tss.dat[i,j]<-TSS(cmdat)
  }
}

kable(pCC.dat, caption="Correct classification rates")
kable(sens.dat, caption="Sensitivity")
kable(spec.dat, caption="Specificity")
kable(kap.dat, caption="kappa statistic")
kable(tss.dat, caption= "True Skill statistic")

best.pcc<-apply(pCC.dat,1,which.max)
best.sens<-apply(sens.dat,1,which.max)
best.spec<-apply(spec.dat,1,which.max)
best.kap<-apply(kap.dat,1,which.max)
best.tss<-apply(tss.dat,1,which.max)


fin.win<-data.frame(AUC_table[,1:3],Total_data=bird.obs$Y,AUC_winners=Model_info[best.model,2],
                    pCC_winners=Model_info[best.pcc,2],
                    sens_winners=Model_info[best.sens,2],
                    Spec_winners=Model_info[best.spec,2],
                    Kappa_winners=Model_info[best.kap,2],
                    TSS_winners=Model_info[best.tss,2]
)
fin.win<-fin.win[order(fin.win$culmen,fin.win$role),]
kable(fin.win)

#Second try: I use separate optimal thresholds for each models. Optimal threshold is chosen to maximize the sum of sensitivity and specificity.

pCC.dat<-matrix(nrow=length(test.bird),ncol=no_models)
rownames(pCC.dat)<-test.bird
colnames(pCC.dat)<-colnames(results1)[3:11]
sens.dat<-pCC.dat
spec.dat<-pCC.dat
kap.dat<-pCC.dat
tss.dat<-pCC.dat

t.dat<-thresh

for(i in 1:length(test.bird)){
  for(j in 1:no_models){
    thr<-t.dat[j]
    res<-results1[which(results$bird==test.bird[i]),c(1,2,j+2)]
    res<-na.omit(res)
    cmdat<-cmx(DATA=res,threshold=thr)
    pCC.dat[i,j]<-pcc(cmdat,st.dev=F)
    sens.dat[i,j]<-sensitivity(cmdat,st.dev=F)
    spec.dat[i,j]<-specificity(cmdat,st.dev=F)
    kap.dat[i,j]<-kappa(cmdat,st.dev=F)
    tss.dat[i,j]<-TSS(cmdat)
  }
}

best.pcc<-apply(pCC.dat,1,which.max)
best.sens<-apply(sens.dat,1,which.max)
best.spec<-apply(spec.dat,1,which.max)
best.kap<-apply(kap.dat,1,which.max)
best.tss<-apply(tss.dat,1,which.max)


fin.win<-data.frame(AUC_table[,1:3],Total_data=bird.obs$Y,AUC_winners=Model_info[best.model,2],
                    pCC_winners=Model_info[best.pcc,2],
                    sens_winners=Model_info[best.sens,2],
                    Spec_winners=Model_info[best.spec,2],
                    Kappa_winners=Model_info[best.kap,2],
                    TSS_winners=Model_info[best.tss,2]
)
fin.win<-fin.win[order(fin.win$culmen,fin.win$role),]
kable(fin.win)

