#############  R Code for Score tests for scale effects, with application to genomic analysis
#############  David Soave, dsoave@wlu.ca, Date: June 2020

library(quantreg)

#### Scale Score test:  
Scale_Test<-function(y,X=NULL,W=NULL,Z=NULL){
  ## MLE's (under H0):
  y=y
  X=X
  W=W
  Z=Z
  I1=matrix(rep(1,length(y)))
  X0<-cbind(I1,Z)
  
  r1=0;r2=0
  n=length(y)
  
  pred.m0q<-predict(rq(lm(y~X0-1)))
  pred.m0<-predict((lm(y~X0)))
  
  G0hat=1/2*log(1/n*sum((y-pred.m0)^2));G0hat
  exp(G0hat)
  
  bi=exp(G0hat)
  ei0q=(y-pred.m0q)/bi
  ei0=(y-pred.m0)/bi
  ei<-ei0
  
  

  ### Permutation Approximation Scale Tests VPN, VPL
  ###################################################
    sx<-t(W-mean(W))%*%(W-mean(W))
  
    alpha_i<-(ei0^2-1)
    R<-cor(alpha_i,W)
    testi<-sqrt((n-2)*R^2/(1-R^2))
    VPN_p<-2*(1-pt(testi,df=(n-2)))
    VPN_p
    
    alpha_i<-(abs(ei0q))
    R<-cor(alpha_i,W);R
    testi<-sqrt((n-2)*R^2/(1-R^2))
    VPL_p<-2*(1-pt(testi,df=(n-2)))
    VPL_p
 

    
  ### SCale Score Test with Model Based (Gausian) Covriance Matrix
  ### (Obvserved or Expected (Fisher) Information Matrix)
 
    U10=NULL;U20=NULL
    DBB=NULL;DBG=NULL;DGB=NULL;DGG=NULL
    DBB0=NULL;DBG0=NULL;DGB0=NULL;DGG0=NULL
    IBB=NULL;IBG=NULL;IGB=NULL;IGG=NULL
    IBB0=NULL;IBG0=NULL;IGB0=NULL;IGG0=NULL
    
    # Observed Information Matrix (not currently used in the results)
    {     
    if(!is.null(X)){
      r1=dim(data.matrix(X))[2]
      U10=t(ei*X/bi)%*%I1;U10
      DBB=t(ei^2*X)%*%X/bi^2;DBB
      DBB0=t(ei^2*X)%*%X0/bi^2;DBB0
      DBG0=t(-ei*(1-ei^2)*X)%*%I1/bi;DBG0
      IBB=(t(X)%*%(X)/bi^2);IBB
      IBB0=t(X/bi^2)%*%X0;IBB0
      IBG0=t(2*t(ei)%*%X/bi);IBG0
    }
    if(!is.null(X)&!is.null(W)){
      DBG=t(-ei*(1-ei^2)*X)%*%W/bi;DBG
      DGB=t(DBG);DGB
      IBG=(2*t(ei*X)%*%(W)/bi);IBG
      IGB=t(IBG);IGB
    }
    if(!is.null(W)){
      r2=dim(data.matrix(W))[2]
      U20=t((1-ei^2)*-W)%*%I1;U20
      DGG=t((1-ei^2)^2*W)%*%W;DGG
      DGB0=t(-ei*(1-ei^2)*W)%*%X0/bi;DGB0
      DGG0=t((1-ei^2)^2*W)%*%I1;DGG0
      IGG=(2*t((ei^2)*W)%*%(W));IGG
      IGB0=(2*t((ei)*W)%*%X0/bi);IGB0
      IGG0=t(2*t(ei^2)%*%W);IGG0
    }
    
    ## Score vector U (under H0):
    U0=matrix(c(U10,U20));U0
    
    ## Variance of Score vector Var(U) = I11 - (I10) (I00^-1) (I01) (Under H0):
    I11=(rbind(cbind(IBB,IBG),cbind(IGB,IGG)));I11
    I10=(rbind(cbind(IBB0, IBG0),cbind(IGB0, IGG0)));I10
    I01=t(I10)
    
    IB0B0=t(X0/bi^2)%*%X0; IB0B0
    IB0G0=2*t(ei*X0/bi)%*%I1; IB0G0   #IB0G0=0; IB0G0
    IG0B0=t(IB0G0)
    IG0G0=2*t(ei)%*%ei;IG0G0
    I00=(rbind(cbind(IB0B0,IB0G0),cbind(IG0B0,IG0G0)));I00
    I=(rbind(cbind(I11, I10),cbind(I01,I00)));I
    
    
    ## Observed Information Matrix
    VarU0=I11-I10%*%solve(I00)%*%I01;VarU0
    IV0=solve(VarU0);IV0
    S0=t(U0)%*%IV0%*%U0;S0
    df=r1+r2
    p=1-pchisq(S0,df);p
    # You may add VNobs_p=p to the results list at end of function if you would like 
    # to include the score test result using model-based observed information matrix
    }
    
    
    ## Fisher Expected Information -> Variance of Score vector Var(U) = I11 - (I10) (I00^-1) (I01) (Under H0):
    {
    if(!is.null(X)){
      r1=dim(data.matrix(X))[2]
      U10=t(ei*X/bi)%*%I1;U10
      DBB=t(ei^2*X)%*%X/bi^2;DBB
      DBB0=t(ei^2*X)%*%X0/bi^2;DBB0
      DBG0=0*t(-ei*(1-ei^2)*X)%*%I1/bi;DBG0
      IBB=(t(X)%*%(X)/bi^2);IBB
      IBG=(0*t(ei*X)%*%(W)/bi);IBG
      IBB0=t(X/bi^2)%*%X0;t(t(IBB0))
      IBG0=t(0*t(ei)%*%X/bi);IBG0
    }
    if(!is.null(X)&!is.null(W)){
      DBG=0*t(-ei*(1-ei^2)*X)%*%W/bi;DBG
      DGB=t(DBG);DGB
      IBG=(0*t(ei*X)%*%(W)/bi);IBG
      IGB=t(IBG);IGB
    }
    if(!is.null(W)){
      r2=dim(data.matrix(W))[2]
      U20=t((1-ei^2)*-W)%*%I1;U20
      DGG=t((1-ei^2)^2*W)%*%W;DGG
      DGB0=0*t(-ei*(1-ei^2)*W)%*%X0/bi;DGB0
      DGG0=t((1-ei^2)^2*W)%*%I1;DGG0
      IGG=(2*t((ei^2)*W)%*%(W));IGG
      IGG=(2*t(W)%*%(W));IGG
      IGB0=(0*t((ei)*W)%*%X0/bi);IGB0
      IGG0=t(2*t(I1^2)%*%W);IGG0
    }
    
    ## Score vector U (under H0):
    U0=matrix(c(U10,U20));U0
    
    ## Variance of Score vector Var(U) = I11 - (I10) (I00^-1) (I01) (Under H0):
    I11=(rbind(cbind(IBB,IBG),cbind(IGB,IGG)));I11
    I10=(rbind(cbind(IBB0, IBG0),cbind(IGB0,IGG0)));I10
    I01=t(I10)
    
    IB0B0=t(X0/bi^2)%*%X0; IB0B0
    IB0G0=0*t(ei*X0/bi)%*%I1; IB0G0   #IB0G0=0; IB0G0
    IG0B0=t(IB0G0)
    IG0G0=2*t(I1)%*%I1;IG0G0
    I00=(rbind(cbind(IB0B0,IB0G0),cbind(IG0B0,IG0G0)));I00
    I=(rbind(cbind(I11, I10),cbind(I01,I00)));I
    
    ## Variance of Score vector Var(U) = I11 - (I10) (I00^-1) (I01) (Under H0):
    VarU0F=I11-I10%*%solve(I00)%*%I01;VarU0F
    IV0=solve(VarU0F);IV0
    SF=t(U0)%*%IV0%*%U0;SF
    df=r1+r2
    pF=1-pchisq(SF,df);pF
  }
  
  

  list(VN_p=pF,VPN_p=VPN_p,VPL_p=VPL_p)
}



