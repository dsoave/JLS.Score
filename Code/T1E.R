#############  R Code for Score tests for scale effects, with application to genomic analysis
#############  David Soave, dsoave@wlu.ca, Date: June 2020

source('../Functions/funcs1.R')
library(parallel)
library(xtable)

### Data generating and testing function

sim1<-function(seedi=i,covZ="Bin"){
  
  print(seedi)
  set.seed(111*seedi)
  
  XG<-rbinom(n,1,pX)
  ZG<-NULL
  if(covZ=="Bin"){
    ZG<-rbinom(n,1,pZ)
  }else if(covZ=="Norm"){
    ZG<-rnorm(n,0,1)
  }
  
  if(sum(XG)<2){XG[sample(1:n,2)]=1}
  if(covZ=="Bin"&sum(ZG)<2){ZG[sample(1:n,2)]=1}

  if(dist=="logX2"){
    y<-bZ*ZG+log(rchisq(n,2))
  } else if(dist=="t4"){
    y<-bZ*ZG+rt(n,4)
  } else if(dist=="t10"){
    y<-bZ*ZG+rt(n,10)
  } else if(dist=="logNorm"){
    y<-bZ*ZG+rlnorm(n,meanlog = 2, sdlog = .75)
  } else if(dist=="Normal"){
    y<-bZ*ZG+rnorm(n)
  }
  

  W=XG;Z=cbind(XG,ZG)

  Scale.Score.perm<-Scale_Test(y=y,X=NULL,W=XG,Z=cbind(XG,ZG))
  
  VN<-Scale.Score.perm$VN_p
  VPN<-Scale.Score.perm$VPN_p
  VPL<-Scale.Score.perm$VPL_p
  
  c(VN,VPN,VPL)  
}

sim1gz<-function(seedi=i,covZ="Bin"){
  
  print(seedi)
  set.seed(111*seedi)
  
  XG<-rbinom(n,1,pX)
  ZG<-NULL
  if(covZ=="Bin"){
    ZG<-rbinom(n,1,pZ)
  }else if(covZ=="Norm"){
    ZG<-rnorm(n,0,1)
  }
  
  if(sum(XG)<2){XG[sample(1:n,2)]=1}
  if(covZ=="Bin"&sum(ZG)<2){ZG[sample(1:n,2)]=1}
  
  if(dist=="logX2"){
    y<-bZ*ZG+log(rchisq(n,2))*exp(gZ*ZG)
  } else if(dist=="t4"){
    y<-bZ*ZG+rt(n,4)*exp(gZ*ZG)
  } else if(dist=="t10"){
    y<-bZ*ZG+rt(n,10)*exp(gZ*ZG)
  } else if(dist=="logNorm"){
    y<-bZ*ZG+rlnorm(n,meanlog = 2, sdlog = .75)*exp(gZ*ZG)
  } else if(dist=="Normal"){
    y<-bZ*ZG+rnorm(n)*exp(gZ*ZG)
  }
  
  
  W=XG;Z=cbind(XG,ZG)
  
  Scale.Score.perm<-Scale_Test(y=y,X=NULL,W=XG,Z=cbind(XG,ZG))
  
  VN<-Scale.Score.perm$VN_p
  VPN<-Scale.Score.perm$VPN_p
  VPL<-Scale.Score.perm$VPL_p
  
  c(VN,VPN,VPL)  
}


#############################################################################
### Generate Results for Table 1 Web Tables 1-9 (S1-S9)

##### Common simulation parameters
reps=100000

parms<-matrix(rbind(c(100,0.1,0.2,0.3,0.5),
                    c(500,0.05,0.1,0.2,0.3),
                    c(1000,0.01,0.05,0.1,0.3)),nrow=3)

##### Table 1
Tablei="1"
dist<-"Normal"
covZd="Bin";pZ<-0.3;bZ<-1

for (j in 1:3){
  
  n<- parms[j,1]
  pXv<-parms[j,2:5]
  
  for (pX in pXv){
    
    print(pX)
    pvals<-data.frame(do.call(cbind,mclapply(1:reps,function(i,covZ){sim1(i,covZ=covZd)},mc.cores=20)))
    
    power.05<-apply(pvals,1,function(x){sum(x<0.05)/length(x)});power.05
    power.01<-apply(pvals,1,function(x){sum(x<0.01)/length(x)});power.01
    power.001<-apply(pvals,1,function(x){sum(x<0.001)/length(x)});power.001
    
    objsave<-rbind(power.05,power.01,power.001);objsave
    
    saveRDS(objsave, paste("../Results/Table",Tablei,"n",n,"pX",pX,".rds",sep="_"))
    
  }
}




##### Web Table 1
Tablei="S1"
dist<-"Normal"
covZd="Norm";bZ<-1


for (j in 1:3){
  
  n<- parms[j,1]
  pXv<-parms[j,2:5]
  
  for (pX in pXv){
    
    print(pX)
    pvals<-data.frame(do.call(cbind,mclapply(1:reps,function(i,covZ){sim1(i,covZ=covZd)},mc.cores=20)))
    
    power.05<-apply(pvals,1,function(x){sum(x<0.05)/length(x)});power.05
    power.01<-apply(pvals,1,function(x){sum(x<0.01)/length(x)});power.01
    power.001<-apply(pvals,1,function(x){sum(x<0.001)/length(x)});power.001
    
    objsave<-rbind(power.05,power.01,power.001);objsave
    
    saveRDS(objsave, paste("../Results/Table",Tablei,"n",n,"pX",pX,".rds",sep="_"))
    
  }
}


for (j in 1:3){
  
  n<- parms[j,1]
  pXv<-parms[j,2:5]
  
  mat1<-NULL
  for (pX in pXv){
    print(pX)
    mat1<-cbind(mat1,readRDS(paste("../Results/Table",Tablei,"n",n,"pX",pX,".rds",sep="_")));mat1               
    alphas<-c(0.05,0.01,0.001)
    
  }
  print(xtable(cbind(alphas,mat1),digits=c(0,3,4,4,4,4,4,4,4,4,4,4,4,4)),include.rownames=FALSE,include.colnames=FALSE)
}



##### Web Table 2
Tablei="S2"
dist<-"t10"
covZd="Bin";pZ<-0.3;bZ<-1


for (j in 1:3){
  
  n<- parms[j,1]
  pXv<-parms[j,2:5]
  
  for (pX in pXv){
    
    print(pX)
    pvals<-data.frame(do.call(cbind,mclapply(1:reps,function(i,covZ){sim1(i,covZ=covZd)},mc.cores=20)))
    
    power.05<-apply(pvals,1,function(x){sum(x<0.05)/length(x)});power.05
    power.01<-apply(pvals,1,function(x){sum(x<0.01)/length(x)});power.01
    power.001<-apply(pvals,1,function(x){sum(x<0.001)/length(x)});power.001
    
    objsave<-rbind(power.05,power.01,power.001);objsave
    
    saveRDS(objsave, paste("../Results/Table",Tablei,"n",n,"pX",pX,".rds",sep="_"))
    
  }
}



##### Web Table 3
gc()
Tablei="S3"
dist<-"t10"
covZd="Norm";bZ<-1


for (j in 1:3){
  
  n<- parms[j,1]
  pXv<-parms[j,2:5]
  
  for (pX in pXv){
    
    print(pX)
    pvals<-data.frame(do.call(cbind,mclapply(1:reps,function(i,covZ){sim1(i,covZ=covZd)},mc.cores=20)))
    
    power.05<-apply(pvals,1,function(x){sum(x<0.05)/length(x)});power.05
    power.01<-apply(pvals,1,function(x){sum(x<0.01)/length(x)});power.01
    power.001<-apply(pvals,1,function(x){sum(x<0.001)/length(x)});power.001
    
    objsave<-rbind(power.05,power.01,power.001);objsave
    
    saveRDS(objsave, paste("../Results/Table",Tablei,"n",n,"pX",pX,".rds",sep="_"))
    
  }
}



##### Web Table 4
gc()
Tablei="S4"
dist<-"t4"
covZd="Bin";pZ<-0.3;bZ<-1


for (j in 1:3){
  
  n<- parms[j,1]
  pXv<-parms[j,2:5]
  
  for (pX in pXv){
    
    print(pX)
    pvals<-data.frame(do.call(cbind,mclapply(1:reps,function(i,covZ){sim1(i,covZ=covZd)},mc.cores=20)))
    
    power.05<-apply(pvals,1,function(x){sum(x<0.05)/length(x)});power.05
    power.01<-apply(pvals,1,function(x){sum(x<0.01)/length(x)});power.01
    power.001<-apply(pvals,1,function(x){sum(x<0.001)/length(x)});power.001
    
    objsave<-rbind(power.05,power.01,power.001);objsave
    
    saveRDS(objsave, paste("../Results/Table",Tablei,"n",n,"pX",pX,".rds",sep="_"))
    
  }
}


##### Web Table 5
gc()
Tablei="S5"
dist<-"t4"
covZd="Norm";bZ<-1


for (j in 1:3){
  
  n<- parms[j,1]
  pXv<-parms[j,2:5]
  
  for (pX in pXv){
    
    print(pX)
    pvals<-data.frame(do.call(cbind,mclapply(1:reps,function(i,covZ){sim1(i,covZ=covZd)},mc.cores=20)))
    
    power.05<-apply(pvals,1,function(x){sum(x<0.05)/length(x)});power.05
    power.01<-apply(pvals,1,function(x){sum(x<0.01)/length(x)});power.01
    power.001<-apply(pvals,1,function(x){sum(x<0.001)/length(x)});power.001
    
    objsave<-rbind(power.05,power.01,power.001);objsave
    
    saveRDS(objsave, paste("../Results/Table",Tablei,"n",n,"pX",pX,".rds",sep="_"))
    
  }
}




##### Web Table 6
gc()
Tablei="S6"
dist<-"logX2"
covZd="Bin";pZ<-0.3;bZ<-1


for (j in 1:3){
  
  n<- parms[j,1]
  pXv<-parms[j,2:5]
  
  for (pX in pXv){
    
    print(pX)
    pvals<-data.frame(do.call(cbind,mclapply(1:reps,function(i,covZ){sim1(i,covZ=covZd)},mc.cores=20)))
    
    power.05<-apply(pvals,1,function(x){sum(x<0.05)/length(x)});power.05
    power.01<-apply(pvals,1,function(x){sum(x<0.01)/length(x)});power.01
    power.001<-apply(pvals,1,function(x){sum(x<0.001)/length(x)});power.001
    
    objsave<-rbind(power.05,power.01,power.001);objsave
    
    saveRDS(objsave, paste("../Results/Table",Tablei,"n",n,"pX",pX,".rds",sep="_"))
    
  }
}


##### Web Table 7
gc()
Tablei="S7"
dist<-"logX2"
covZd="Norm";bZ<-1


for (j in 1:3){
  
  n<- parms[j,1]
  pXv<-parms[j,2:5]
  
  for (pX in pXv){
    
    print(pX)
    pvals<-data.frame(do.call(cbind,mclapply(1:reps,function(i,covZ){sim1(i,covZ=covZd)},mc.cores=20)))
    
    power.05<-apply(pvals,1,function(x){sum(x<0.05)/length(x)});power.05
    power.01<-apply(pvals,1,function(x){sum(x<0.01)/length(x)});power.01
    power.001<-apply(pvals,1,function(x){sum(x<0.001)/length(x)});power.001
    
    objsave<-rbind(power.05,power.01,power.001);objsave
    
    saveRDS(objsave, paste("../Results/Table",Tablei,"n",n,"pX",pX,".rds",sep="_"))
    
  }
}



##### Web Table 8
gc()
Tablei="S8"
dist<-"logNorm"
covZd="Bin";pZ<-0.3;bZ<-1

j=1
pX=0.1
for (j in 1:3){
  
  n<- parms[j,1]
  pXv<-parms[j,2:5]
  
  for (pX in pXv){
    
    print(pX)
    pvals<-data.frame(do.call(cbind,mclapply(1:reps,function(i,covZ){sim1(i,covZ=covZd)},mc.cores=20)))
    
    power.05<-apply(pvals,1,function(x){sum(x<0.05)/length(x)});power.05
    power.01<-apply(pvals,1,function(x){sum(x<0.01)/length(x)});power.01
    power.001<-apply(pvals,1,function(x){sum(x<0.001)/length(x)});power.001
    
    objsave<-rbind(power.05,power.01,power.001);objsave
    
    saveRDS(objsave, paste("../Results/Table",Tablei,"n",n,"pX",pX,".rds",sep="_"))
    
  }
}


##### Web Table 9
gc()
Tablei="S9"
dist<-"logNorm"
covZd="Norm";bZ<-1


for (j in 1:3){
  
  n<- parms[j,1]
  pXv<-parms[j,2:5]
  
  for (pX in pXv){
    
    print(pX)
    pvals<-data.frame(do.call(cbind,mclapply(1:reps,function(i,covZ){sim1(i,covZ=covZd)},mc.cores=20)))
    
    power.05<-apply(pvals,1,function(x){sum(x<0.05)/length(x)});power.05
    power.01<-apply(pvals,1,function(x){sum(x<0.01)/length(x)});power.01
    power.001<-apply(pvals,1,function(x){sum(x<0.001)/length(x)});power.001
    
    objsave<-rbind(power.05,power.01,power.001);objsave
    
    saveRDS(objsave, paste("../Results/Table",Tablei,"n",n,"pX",pX,".rds",sep="_"))
    
  }
}






##### Web Table 10
gc()
Tablei="1gz"
dist<-"Normal"
covZd="Bin";pZ<-0.3;bZ<-1;gZ<-0.5

for (j in 1:3){
  
  n<- parms[j,1]
  pXv<-parms[j,2:5]
  
  for (pX in pXv){
    
    print(pX)
    pvals<-data.frame(do.call(cbind,mclapply(1:reps,function(i,covZ){sim1gz(i,covZ=covZd)},mc.cores=20)))
    
    power.05<-apply(pvals,1,function(x){sum(x<0.05)/length(x)});power.05
    power.01<-apply(pvals,1,function(x){sum(x<0.01)/length(x)});power.01
    power.001<-apply(pvals,1,function(x){sum(x<0.001)/length(x)});power.001
    
    objsave<-rbind(power.05,power.01,power.001);objsave
    
    saveRDS(objsave, paste("../Results/Table",Tablei,"n",n,"pX",pX,".rds",sep="_"))
    
  }
}



##### Web Table 11
gc()
Tablei="S4gz"
dist<-"t4"
covZd="Bin";pZ<-0.3;bZ<-1;gZ<-0.5


for (j in 1:3){
  
  n<- parms[j,1]
  pXv<-parms[j,2:5]
  
  for (pX in pXv){
    
    print(pX)
    pvals<-data.frame(do.call(cbind,mclapply(1:reps,function(i,covZ){sim1gz(i,covZ=covZd)},mc.cores=20)))
    
    power.05<-apply(pvals,1,function(x){sum(x<0.05)/length(x)});power.05
    power.01<-apply(pvals,1,function(x){sum(x<0.01)/length(x)});power.01
    power.001<-apply(pvals,1,function(x){sum(x<0.001)/length(x)});power.001
    
    objsave<-rbind(power.05,power.01,power.001);objsave
    
    saveRDS(objsave, paste("../Results/Table",Tablei,"n",n,"pX",pX,".rds",sep="_"))
    
  }
}



##### Web Table 12
gc()
Tablei="S8gz"
dist<-"logNorm"
covZd="Bin";pZ<-0.3;bZ<-1;gZ<-0.5

j=1
pX=0.1
for (j in 1:3){
  
  n<- parms[j,1]
  pXv<-parms[j,2:5]
  
  for (pX in pXv){
    
    print(pX)
    pvals<-data.frame(do.call(cbind,mclapply(1:reps,function(i,covZ){sim1gz(i,covZ=covZd)},mc.cores=20)))
    
    power.05<-apply(pvals,1,function(x){sum(x<0.05)/length(x)});power.05
    power.01<-apply(pvals,1,function(x){sum(x<0.01)/length(x)});power.01
    power.001<-apply(pvals,1,function(x){sum(x<0.001)/length(x)});power.001
    
    objsave<-rbind(power.05,power.01,power.001);objsave
    
    saveRDS(objsave, paste("../Results/Table",Tablei,"n",n,"pX",pX,".rds",sep="_"))
    
  }
}


#############################################################################
### Produce *.tex files to be imported using in manuscript latex


######### Run each of the following #commented lines (one at a time) and then the subsequent for loop

#dist<-"Normal";covZd="Bin";pZ<-0.3;bZ<-1;Tablei="1"
#dist<-"Normal";covZd="Norm";bZ<-1;Tablei="S1"

#dist<-"t10";covZd="Bin";pZ<-0.3;bZ<-1;Tablei="S2"
#dist<-"t10";covZd="Norm";bZ<-1;Tablei="S3"

#dist<-"t4";covZd="Bin";pZ<-0.3;bZ<-1;Tablei="S4"
#dist<-"t4";covZd="Norm";bZ<-1;Tablei="S5"

#dist<-"logX2";covZd="Bin";pZ<-0.3;bZ<-1;Tablei="S6"
#dist<-"logX2";covZd="Norm";bZ<-1;Tablei="S7"

#dist<-"logNorm";covZd="Bin";pZ<-0.3;bZ<-1;Tablei="S8"
#dist<-"logNorm";covZd="Norm";bZ<-1;Tablei="S9"

#dist<-"Normal";covZd="Bin";pZ<-0.3;bZ<-1;Tablei="1gz"
#dist<-"t4";covZd="Bin";pZ<-0.3;bZ<-1;Tablei="S4gz"
#dist<-"logNorm";covZd="Bin";pZ<-0.3;bZ<-1;Tablei="S8gz"


for (j in 1:3){
  
  n<- parms[j,1]
  pXv<-parms[j,2:5]
  
  mat1<-NULL
  for (pX in pXv){
    print(pX)
    mat1<-cbind(mat1,readRDS(paste("../Results/Table",Tablei,"n",n,"pX",pX,".rds",sep="_")));mat1               
    alphas<-c(0.05,0.01,0.001)
    
  }
  idx<-paste(Tablei,"n",n,sep="_")
  print(xtable(cbind(alphas,mat1),digits=c(0,3,4,4,4,4,4,4,4,4,4,4,4,4)),type='latex',file=paste("../Results/table_",idx,".tex",sep=""),only.contents=TRUE, hline.after=NULL,include.colnames=FALSE,include.rownames=FALSE,sanitize.text.function = function(x) x)
  
}
