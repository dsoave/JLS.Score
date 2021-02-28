#############  R Code for Score tests for scale effects, with application to genomic analysis
#############  David Soave, dsoave@wlu.ca, Date: June 2020

source('../Functions/funcs1.R')

### Generate Results for Figure 1 (1a-1e) and Web Figure 1 (S1a-S1e)

n<-1000
nn=500
a<-0
covZ<-"Norm";pZ<-NA;bZ<-1

Figures<-c("1a","1b","1c","1d","1e","S1a","S1b","S1c","S1d","S1e")

for(Figurei in Figures){
  
  print(Figurei)
  if(Figurei %in% c("1a","S1a")){
    dist<-'Norm';gX<-seq(0,0.26,.26/10)
  }else if(Figurei %in% c("1b","S1b")){
    dist<-'t10';gX<-seq(0,0.32,.32/10)
  }else if(Figurei %in% c("1c","S1c")){
    dist<-'t4';gX<-seq(0,0.5,.5/10)
  }else if(Figurei %in% c("1d","S1d")){
    dist<-'logX2';gX<-seq(0,0.32,.32/10)
  }else if(Figurei %in% c("1e","S1e")){
    dist<-'logNorm';gX<-seq(0,0.55,.55/10)
  }
  
  if(Figurei %in% c("1a","1b","1c","1d","1e")){
    pX=0.3
  }else if(Figurei %in% c("S1a","S1b","S1c","S1d","S1e")){
    pX=0
    gX<-gX*.6
  }
  
  set.seed(538785)
  powerC=NULL
  
  for(g1 in gX){
    
    print(g1)
    sim1 <-replicate(
      nn,
      expr = {
        
        if(pX==0){
          XG<-rnorm(n,0,1)
        }else if(pX>0){
          XG<-rbinom(n,1,pX)
        }
        
        ZG<-NULL
        if(covZ=="Bin"){
          ZG<-a*XG+rbinom(n,1,pZ)
        }else if(covZ=="Norm"){
          ZG<-a*XG+rnorm(n,0,1)
        }
        
        if(covZ=="Bin"&sum(XG)<2){XG[sample(1:n,2)]=1}
        if(dist=="logX2"){
          y<-bZ*ZG+log(rchisq(n,2))*exp(g1*XG)
        } else if(dist=="t4"){
          y<-bZ*ZG+rt(n,4)*exp(g1*XG)
        } else if(dist=="t2"){
          y<-bZ*ZG+rt(n,2)*exp(g1*XG)
        } else if(dist=="t10"){
          y<-bZ*ZG+rt(n,10)*exp(g1*XG)
        } else if(dist=="logNorm"){
          y<-bZ*ZG+rlnorm(n,meanlog = 2, sdlog = .75)*exp(g1*XG)
        } else if(dist=="Norm"){
          y<-bZ*ZG+rnorm(n)*exp(g1*XG)
        }
      
        
        Scale.Score.perm<-Scale_Test(y=y,X=NULL,W=XG,Z=cbind(XG,ZG))
        
        VN<-Scale.Score.perm$VN_p
        VPN<-Scale.Score.perm$VPN_p
        VPL<-Scale.Score.perm$VPL_p
        
        c(VN,VPN,VPL)  
      }
    )
    
    powerC.t<-apply(sim1,1,function(x){sum(x<0.05)/length(x)})
    powerC<-rbind(powerC,c(g1,powerC.t));print(powerC)
    
  }
  
  powerC<-data.frame(powerC)
  names(powerC)<-c("g1","VN","VPN","VPL")
  head(powerC)
  saveRDS(powerC, paste("../Results/Figure",Figurei,".rds",sep="_"))
}





### Create Images for Figure 1 (1a-1e) and Web Figure 1 (S1a-S1e)

X11(width=8, height=7,type="cairo")
par(mar=c(5, 4+1, 4-1, 2) + 0.1) 


for(Figurei in Figures){
  
  print(Figurei)
  
  powerC<-readRDS(paste("../Results/Figure",Figurei,".rds",sep="_"));powerC
  gammaX<-powerC[,1]
  plot(gammaX,powerC[,3],type='l',ylab='Power',xlab=expression(paste("Scale Effect ", gamma[X] )),ylim=c(0,1),col=gray(.1),lty=1,lwd=3,cex.lab=1.5,cex.axis=1.5)
  lines(gammaX,powerC[,2],col=gray(.4),lty=3,lwd=3)
  lines(gammaX,powerC[,4],col=4,lty=2,lwd=3)
  
  
  if(Figurei %in% c("1a","S1a")){
    title(expression(paste(N(0,1) )),cex.main=2)
  }else if(Figurei %in% c("1b","S1b")){
    title(expression(paste("Student's ", t(10) )),cex.main=2)
  }else if(Figurei %in% c("1c","S1c")){
    title(expression(paste("Student's ", t(4) )),cex.main=2)
  }else if(Figurei %in% c("1d","S1d")){
    title(expression(paste("Extreme Value")),cex.main=2)
  }else if(Figurei %in% c("1e","S1e")){
    title(expression(paste("log-Normal")),cex.main=2)
  }
  dev.copy2pdf(device = x11,file=paste("../Results/F",Figurei,".pdf",sep=""), out.type = "pdf")
  
} 




X11(width=11.3, height=2.5,type="cairo")
par(xpd=T, mar=par()$mar+c(0,0,0,0))
plot(1:5,1:5,axes=F,xlab='',ylab='',col='white')
legend(0.8,.5,c(expression(paste("V"[N] )),expression(paste("V"[P(N)] )),expression(paste("V"[P(L)] ))), lty = c(3,1,2),lwd=c(8),
       col=c(gray(.4),gray(.1),4)
       ,horiz=T,box.lwd=0,cex=2.5)
dev.copy2pdf(device = x11,file=paste("../Results/legend_FS2.pdf",sep=""), out.type = "pdf")
dev.off()


### Create Figure Legend
X11(width=5.5, height=4.5,type="cairo")
par(xpd=T, mar=par()$mar+c(0,0,0,0))
plot(1:5,1:5,axes=F,xlab='',ylab='',col='white')
legend(0.8,5,c(expression(paste("V"[N] )),expression(paste("V"[P(N)] )),expression(paste("V"[P(L)] ))), lty = c(3,1,2),lwd=c(8),
       col=c(gray(.4),gray(.1),4)
       ,horiz=F,box.lwd=NULL,cex=2.5)
dev.copy(png,file=paste("../Results/legend_FV.png",sep=""))
dev.off()
