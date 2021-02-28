#############  R Code for Score tests for scale effects, with application to genomic analysis
#############  David Soave, dsoave@wlu.ca, Date: June 2020

source('../Functions/funcs1.R')

### Generate Results for Web Figure 2 (S2a-S2d)

Figures<-c("S2a","S2b","S2c","S2d")

n<-2000
nn=500
siglev<-5e-8
pA=0.3
pE1=0.3

BE1=0.3
bXE<-seq(0,1,.1)

for(Figurei in Figures){
  
  print(Figurei)
  if(Figurei %in% c("S2a")){
    BG=0.01;bs1<-1
  }else if(Figurei %in% c("S2b")){
    BG=0.1;bs1<-1
  }else if(Figurei %in% c("S2c")){
    BG=0.01;bs1<-(-1)
  }else if(Figurei %in% c("S2d")){
    BG=0.1;bs1<-(-1)
  }
  
  BE1<-bs1*0.3
  
  set.seed(538785)
  powerC=NULL
  
  for(bxe1 in bXE){  
    
    bxe1<-bs1*bxe1
    print(c(BG,bxe1))
    
    sim1 <-replicate(
      nn,
      expr = {
        
        genocount<-rmultinom(1,size=n,prob=c(pA*pA, 2*pA*(1-pA), (1-pA)*(1-pA)))
        XG<-c(rep(0, genocount[1]), rep(1, genocount[2]), rep(2,genocount[3]))
        XG<-sample(XG,size=length(XG),replace=FALSE)
        E1=rbinom(n,1,prob=pE1)
        y=BG*XG+BE1*E1+bxe1*XG*E1+rnorm(n,0,1)
        

          Scale.Score.perm<-Scale_Test(y=y,X=NULL,W=XG,Z=cbind(XG))
          
          VN<-Scale.Score.perm$VN_p
          VPN<-Scale.Score.perm$VPN_p
          VPL<-Scale.Score.perm$VPL_p
          
          c(VN,VPN,VPL)  
          
      }
    )
    
    powerC.t<-apply(sim1,1,function(x){sum(x<siglev)/length(x)})
    powerC<-rbind(powerC,c(bxe1,powerC.t));print(powerC)
    
  }
  
  powerC<-data.frame(powerC)
  names(powerC)<-c("biX","VN","VPN","VPL")
  head(powerC)
  saveRDS(powerC, paste("../Results/Figure",Figurei,".rds",sep="_"))
}



### Create Images for Web Figure 2 (S2a-S2d)

X11(width=8, height=7,type="cairo")
par(mar=c(5, 4+1, 4-1, 2) + 0.1) 

for(Figurei in Figures){
  
  print(Figurei)
  
  powerC<-readRDS(paste("../Results/Figure",Figurei,".rds",sep="_"));powerC
  biX<-powerC[,1]

  plot(biX,powerC[,3],type='l',ylab='Power',xlab=expression(paste("Interaction Effect ", beta[xz] )),ylim=c(0,1),col=gray(.1),lty=1,lwd=3,cex.lab=1.5,cex.axis=1.5)
  lines(biX,powerC[,3],col=gray(.4),lty=3,lwd=3)
  lines(biX,powerC[,4],col=4,lty=2,lwd=3)
  
  if(Figurei %in% c("S2a")){
    title(expression(paste(beta[x],"=0.01, ",beta[xz],">0")),cex.main=2)
  }else if(Figurei %in% c("S2b")){
    title(expression(paste(beta[x],"=0.1, ",beta[xz],">0")),cex.main=2)
  }else if(Figurei %in% c("S2c")){
    title(expression(paste(beta[x],"=0.01, ",beta[xz],"<0")),cex.main=2)
  }else if(Figurei %in% c("S2d")){
    title(expression(paste(beta[x],"=0.1, ",beta[xz],"<0")),cex.main=2)
  }
  dev.copy2pdf(device = x11,file=paste("../Results/F",Figurei,".pdf",sep=""), out.type = "pdf")
  
} 
