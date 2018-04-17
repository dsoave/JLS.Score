#' Score Tests for Joint Location and Scale Models
#'
#' This function performs the joint location and scale (JLS) or scale-only score tests (Soave, Lawless an Awadalla, submitted) to simultaneously test for mean and variance (or variance-only) associations with covariates.
#' @param y response vector
#' @param X vector (or matrix) of covariate(s) to be tested for "location" effect(s) on the response. If a matrix, each column represents a variable and each row represents an observation. This may be the same as W.
#' @param W vector (or matrix) of covariate(s) to be tested for "location" effect(s) on the response. If a matrix, each column represents a variable and each row represents an observation. This may be the same as X.
#' @param Z vector (or matrix) of covariate(s) to be included as adjustment variables in the "location" portion of the model. If a matrix, each column represents a variable and each row represents an observation.
#' @keywords JLS Score test
#' @export
#' @author David Soave
#' @details No missing data are allowed - function will return an "error". Outcome must be quantitative and covariates may be discrete (categorical) or continuous.
#' @return a table consisting of test statistics, degrees of freedom and p-vaules for the score test using the observed (W_obs) and expected (W_exp) information covariance matrix, and the robust covariance estimators V_A(D,I), V_B(D,E(I)), and V_C(E(D),E(I)).
#' @references Soave, D., Lawless, J.F., and Awadalla, P. Score tests for association in location and scale models. Submitted.
#' @examples
#' #################################################################################
#' ## Example simulating data from Table 1, row 3 (Soave, Lawless, and Awadalla, submitted)
#' #################################################################################
#'
#' #### Simulation parameters
#' n<-1000 # sample size
#' pX<-0.1 # covariate frequency
#'
#' reps<-100000  # number of simulation replicates
#'
#' #################################################################################
#' ## Simulation replicates
#' #################################################################################
#'
#' sims <-replicate(
#'  reps,
#'  expr = {
#'    XG<-rbinom(n,1,pX)
#'    y<-rnorm(length(XG),0,1)
#'    result_J<-JLS.Score(y=y,X=XG,W=XG,Z=NULL)
#'    result_S<-JLS.Score(y=y,X=NULL,W=XG,Z=XG)
#'    c(result_J[,3],result_S[,3])
#'  }
#' )
#'
#' #################################################################################
#' ## Power estimates
#' #################################################################################
#'
#' 
#' alpha<-0.01
#' power<-rbind(apply(sims,1,function(x){sum(x<alpha)/length(x)}))
#' power.J<-power[,1:5]
#' power.S<-power[,6:10]
#' names(power.J)=names(power.S)=c("W_obs","W_exp","V_A","V_B","V_C")
#' power.J
#' power.S

JLS.Score <-function(y,X=NULL,W=NULL,Z=NULL){

  ## check if there is missing data
  if(sum(is.na(cbind(y,X,W,Z))) > 0)  stop("missing value(s) not allowed")

  I1=matrix(rep(1,length(y)))
  X0<-cbind(I1,Z)
  
  r1=0;r2=0
  n=length(y)
  
  pred.m0<-predict((lm(y~X0)))
  #pred.m0<-predict(rq(lm(y~X0-1)))
  G0hat=1/2*log(1/n*sum((y-pred.m0)^2));G0hat
  exp(G0hat)
  
  bi=exp(G0hat)
  ei=(y-pred.m0)/bi
  ei<-qnorm((rank(ei)-0.5)/length(ei))
  U10=NULL;U20=NULL
  DBB=NULL;DBG=NULL;DGB=NULL;DGG=NULL
  DBB0=NULL;DBG0=NULL;DGB0=NULL;DGG0=NULL
  IBB=NULL;IBG=NULL;IGB=NULL;IGG=NULL
  IBB0=NULL;IBG0=NULL;IGB0=NULL;IGG0=NULL
  
  
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
  
  pU2X<-1-pchisq(U20+sum(W),sum(W));pU2X
  
  ## Score vector U (under H0):
  U0=matrix(c(U10,U20));U0
  
  D11=(rbind(cbind(DBB,DBG),cbind(DGB,DGG)));D11
  D10=(rbind(cbind(DBB0,DBG0),cbind(DGB0,DGG0)));D10
  D01=t(D10);D01
  
  DB0B0=t(ei*X0)%*%(ei*X0)/bi^2;DB0B0
  DB0G0=t(-ei*X0)%*%(1-ei^2)/bi;DB0G0
  DG0B0=t(DB0G0);DG0B0
  DG0G0=t(1-ei^2)%*%(1-ei^2);DG0G0
  D00=(rbind(cbind(DB0B0,DB0G0),cbind(DG0B0,DG0G0)));D00
  
  D=(rbind(cbind(D11, D10),cbind(D01,D00)));D
  
  
  ## Variance of Score vector Var(U) = I11 - (I10) (I00^-1) (I01) (Under H0):
  I11=(rbind(cbind(IBB,IBG),cbind(IGB,IGG)));I11
  
  I10=(rbind(cbind(IBB0, IBG0),cbind(IGB0, IGG0)));I10
  I01=t(I10)
  I01
  
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
  U0
  p=1-pchisq(S0,df);p
  
  US<-NA
  if(min(diag(IV0))>0){US=U0*sqrt(diag(IV0));US}
  
  ## Robust Observed Information Matrix
  VarR1a=D11-I10%*%solve(I00)%*%t(D10)-D10%*%solve(I00)%*%t(I10)+I10%*%solve(I00)%*%D00%*%solve(I00)%*%t(I10)  ;   VarR1a
  IVR1=solve(VarR1a)
  SR1=t(U0)%*%IVR1%*%U0;SR1
  pR1=1-pchisq(SR1,df);pR1
  
  USR=U0*sqrt(diag(IVR1));USR
  
  #VarR2=(VarU0)%*%((solve(I)%*%D%*%solve(I))[1:df,1:df])%*%(VarU0);VarR2  ### Exact same as VarR1
  
  ## Fisher Expected Information -> Variance of Score vector Var(U) = I11 - (I10) (I00^-1) (I01) (Under H0):
  
  if(!is.null(X)){
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
    DGG=t((1-ei^2)^2*W)%*%W;DGG
    DGB0=0*t(-ei*(1-ei^2)*W)%*%X0/bi;DGB0
    DGG0=t((1-ei^2)^2*W)%*%I1;DGG0
    IGG=(2*t((ei^2)*W)%*%(W));IGG
    IGG=(2*t(W)%*%(W));IGG
    IGB0=(0*t((ei)*W)%*%X0/bi);IGB0
    IGG0=t(2*t(I1^2)%*%W);IGG0
  }
  
  I11=(rbind(cbind(IBB,IBG),cbind(IGB,IGG)));I11
  
  
  I10=(rbind(cbind(IBB0, IBG0),cbind(IGB0,IGG0)));I10
  I01=t(I10)
  I01
  
  IB0B0=t(X0/bi^2)%*%X0; IB0B0
  IB0G0=0*t(ei*X0/bi)%*%I1; IB0G0   #IB0G0=0; IB0G0
  IG0B0=t(IB0G0)
  IG0G0=2*t(I1)%*%I1;IG0G0
  I00=(rbind(cbind(IB0B0,IB0G0),cbind(IG0B0,IG0G0)));I00
  
  I=(rbind(cbind(I11, I10),cbind(I01,I00)));I
  
  
  ## Variance of Score vector Var(U) = I11 - (I10) (I00^-1) (I01) (Under H0):
  
  VarU0F=I11-I10%*%solve(I00)%*%I01;VarU0F
  IV0=solve(VarU0F);IV0
  
  SF=t(U0)%*%IV0%*%U0;S0
  df=r1+r2
  pF=1-pchisq(SF,df);pF
  
  USF<-U0*sqrt(diag(IV0));USF
  
  ## robust (b)
  ## Robust Fisher Information Matrix
  VarR1b=D11-I10%*%solve(I00)%*%t(D10)-D10%*%solve(I00)%*%t(I10)+I10%*%solve(I00)%*%D00%*%solve(I00)%*%t(I10)  ;   VarR1b
  IVR1=solve(VarR1b)
  SR1b=t(U0)%*%IVR1%*%U0;SR1b
  pR1b=1-pchisq(SR1b,df);pR1b
  USRb=U0*sqrt(diag(IVR1));USRb
  
  
  
  ## Robust (c)
  D11=(rbind(cbind(DBB,DBG),cbind(DGB,DGG)));D11
  D10=(rbind(cbind(DBB0,DBG0),cbind(DGB0,DGG0)));D10
  D01=t(D10);D01
  
  DB0B0=t(ei*X0)%*%(ei*X0)/bi^2;DB0B0
  DB0G0=0*t(-ei*X0)%*%(1-ei^2)/bi;DB0G0
  DG0B0=t(DB0G0);DG0B0
  DG0G0=t(1-ei^2)%*%(1-ei^2);DG0G0
  D00=(rbind(cbind(DB0B0,DB0G0),cbind(DG0B0,DG0G0)));D00
  
  D=(rbind(cbind(D11, D10),cbind(D01,D00)));D
  
  
  ## Robust Fisher (c) Information Matrix
  VarR1c=D11-I10%*%solve(I00)%*%t(D10)-D10%*%solve(I00)%*%t(I10)+I10%*%solve(I00)%*%D00%*%solve(I00)%*%t(I10)  ;   VarR1c
  
  IVR1=solve(VarR1c)
  SR1c=t(U0)%*%IVR1%*%U0;SR1c
  pR1c=1-pchisq(SR1c,df);pR1c
  
  USRc=U0*sqrt(diag(IVR1));USRc
  
  dataf<-data.frame(rbind(c(S0,df,p),c(SF,df,pF),c(SR1,df,pR1),c(SR1b,df,pR1b),c(SR1c,df,pR1c)))
  names(dataf)=c("Statistic","df","p-value")
  rownames(dataf)=c("W_obs","W_exp","V_A","V_B","V_C")
  return(dataf)
  
  #list(ScoreJ=S0,ScoreF=SF,df=df,pval=p,pvalF=pF,ScoreR=SR1,pvalR=pR1,pvalRb=pR1b,pvalRc=pR1c,US=US,USF=USF,USR=USR,USRb=USRb,USRc=USRc,VarU0=VarU0,VarU0F=VarU0F,VarR1a=VarR1a,VarR1b=VarR1b,VarR1c=VarR1c,U0=U0,ScoreRb=SR1b,ScoreRc=SR1c)
  
}
