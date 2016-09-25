###############################################
###function: Simulate under the inference model
##############################################
simulate.baypass <- function(omega.mat,nsnp=1000,beta.coef=NA,beta.pi=c(1,1),pop.trait=0,sample.size=100,pi.maf=0.05,suffix="sim",remove.fixed.loci=FALSE,coverage=NA){
#coverage = matrix with npop colums or vector of length npop with coverage =>activate poolseq data
#sample.size=matrix with npop colums or vector of length npop with count (if a matrix and poolseq data: poolsize is set to the colMax of sample.size)
 require(mvtnorm)
 npop=nrow(omega.mat)
 
 if(sum(is.na(beta.coef))==0){
  if(length(pop.trait)!=npop){stop("Trait dimension must have the same size as the rank of omega.mat")}
  simu.cov=TRUE
  nsnp.asso=length(beta.coef)
  if(nsnp>0){beta.coef=c(rep(0,nsnp),beta.coef)}
  nsnp=nsnp + nsnp.asso
  alpha.cov = beta.coef %*% t(pop.trait)
 }else{simu.cov=FALSE}
 
 if(length(sample.size)==1){
  NN=matrix(sample.size,nsnp,npop)
  poolsize=rep(sample.size,npop)
  }else{
  sample.size=as.matrix(sample.size)
  if(ncol(sample.size)==1){#c'est un vecteur
    if(nrow(sample.size)!=npop){stop("Sample size dimension must be of length 1 or have the same size as the rank of omega.mat or must be a matrix with npop columns")}
    NN=matrix(rep(sample.size,nsnp),nsnp,npop,byrow=TRUE)
    poolsize=as.numeric(sample.size)
  }else{
   tmp.snp=sample(1:nrow(sample.size),nsnp,replace=TRUE)
   NN=sample.size[tmp.snp,]
   poolsize=apply(sample.size,2,max)
  }}

 if(length(coverage)==1){
  if(is.na(coverage)){
  poolseq=FALSE
  }else{poolseq=TRUE ; NN.coverage=matrix(coverage,nsnp,npop)}
  }else{
   poolseq=TRUE
   coverage=as.matrix(coverage)
   if(ncol(coverage)==1){#c'est un vecteur
    if(nrow(coverage)!=npop){
     stop("Coverage dimension must be of length 1, or have the same size as the rank of omega.mat or must be a matrix with npop columns")
    }
   NN.coverage=matrix(rep(coverage,nsnp),nsnp,npop,byrow=TRUE)
  }else{
   tmp.snp=sample(1:nrow(coverage),nsnp,replace=TRUE)
   NN.coverage=coverage[tmp.snp,]
  }
 }
 
 if(poolseq){
  YY.coverage=NN.coverage*0
  NN=matrix(rep(poolsize,nsnp),nsnp,npop,byrow=TRUE)
  }

 if(length(beta.pi)!=2){stop("beta.pi must be of length 2")}else{
  if(sum(beta.pi==1)==2){Pi=runif(nsnp,pi.maf,1-pi.maf)}else{
   Pi=rbeta(nsnp,beta.pi[1],beta.pi[2])
   Pi[Pi<pi.maf]=pi.maf ; Pi[Pi>1-pi.maf]=1-pi.maf
  }
 }
 
 ALPHA=YY=matrix(0,nsnp,npop)
 for(i in 1:nsnp){
   mat=Pi[i]*(1-Pi[i])*omega.mat
   ALPHA[i,]=rmvnorm(1,rep(Pi[i],npop),mat)
 }
 if(simu.cov){ALPHA=ALPHA + alpha.cov}
 ALPHA_tr=ALPHA ; ALPHA_tr[ALPHA>1]=1 ; ALPHA_tr[ALPHA<0]=0 
 for(i in 1:nsnp){
  for(j in 1:npop){
    YY[i,j]=rbinom(1,size=NN[i,j],prob=ALPHA_tr[i,j])
    if(poolseq){
     YY.coverage[i,j]=rbinom(1,size=NN.coverage[i,j],prob=YY[i,j]/NN[i,j])
    }
  }
  if(i%%(nsnp/10)==0){cat(i,"SNP simulated out of",nsnp,"\n")}
  }
  
 if(remove.fixed.loci){
  ##Attention Bayevn n'accepte pas les locus fixe!
 if(poolseq){tmp.freq=rowSums(YY.coverage)/rowSums(NN.coverage)}else{tmp.freq=rowSums(YY)/rowSums(NN)}
 snp.sel=tmp.freq>0 & tmp.freq<1
 nsnp=sum(snp.sel)
 YY=YY[snp.sel,] ; NN=NN[snp.sel,] ; Pi=Pi[snp.sel] ; ALPHA=ALPHA[snp.sel,] 
 if(poolseq){YY.coverage=YY.coverage[snp.sel,];NN.coverage=NN.coverage[snp.sel,]}
 if(simu.cov){beta.coef=beta.coef[snp.sel]}
 cat("Number of SNPs removed: ",sum(!snp.sel),"\n")
 }

 write.table(file=paste("pi.",suffix,sep=""),Pi,quote=F,col.names=F,row.names=F)
 write.table(file=paste("alpha.",suffix,sep=""),ALPHA,quote=F,col.names=F,row.names=F)

 if(simu.cov){
  write.table(file=paste("betacoef.",suffix,sep=""),beta.coef,quote=F,col.names=F,row.names=F)
  write.table(file=paste("pheno.",suffix,sep=""),t(pop.trait),quote=F,col.names=F,row.names=F)
 }

 mat_nichmnv=cbind(YY,NN)
 all2.pos=2*(1:npop)
 mat_nichmnv[,all2.pos-1]=YY ; mat_nichmnv[,all2.pos]=NN-YY
 write.table(file=paste("G.",suffix,sep="") ,mat_nichmnv,quote=F,col.names=F,row.names=F)

 mat_bayenv=rbind(YY,NN)
 all2.pos=2*(1:nsnp)
 mat_bayenv[all2.pos-1,]=YY ; mat_bayenv[all2.pos,]=NN-YY
 mat_bayenv=cbind(mat_bayenv,rep("",nsnp))
 write.table(file=paste("bayenv_freq.",suffix,sep=""),mat_bayenv,sep="\t",quote=F,col.names=F,row.names=F)

 if(poolseq){
   mat_nichmnv=cbind(YY.coverage,NN.coverage)
   all2.pos=2*(1:npop)
   mat_nichmnv[,all2.pos-1]=YY.coverage ; mat_nichmnv[,all2.pos]=NN.coverage-YY.coverage
   write.table(file=paste("Gpool.",suffix,sep="") ,mat_nichmnv,quote=F,col.names=F,row.names=F)

  mat_bayenv=rbind(YY.coverage,NN.coverage)
  all2.pos=2*(1:nsnp)
  mat_bayenv[all2.pos-1,]=YY.coverage ; mat_bayenv[all2.pos,]=NN.coverage-YY.coverage
  mat_bayenv=cbind(mat_bayenv,rep("",nsnp))
  write.table(file=paste("bayenv_freq_pool.",suffix,sep=""),mat_bayenv,sep="\t",quote=F,col.names=F,row.names=F)
  
  write.table(file=paste("poolsize.",suffix,sep=""),t(poolsize),quote=F,col.names=F,row.names=F)
 }


 if(simu.cov){
  if(poolseq){
  list(Y.pool=YY.coverage,N.pool=NN.coverage,Y.sim=YY,N.sim=NN,Pi.sim=Pi,alpha.sim=ALPHA,omega.sim=omega.mat,betacoef.sim=beta.coef)
  }else{
  list(Y.sim=YY,N.sim=NN,Pi.sim=Pi,alpha.sim=ALPHA,omega.sim=omega.mat,betacoef.sim=beta.coef)
  }
 }else{
  if(poolseq){
   list(Y.pool=YY.coverage,N.pool=NN.coverage,Y.sim=YY,N.sim=NN,Pi.sim=Pi,alpha.sim=ALPHA,omega.sim=omega.mat)
  }else{
   list(Y.sim=YY,N.sim=NN,Pi.sim=Pi,alpha.sim=ALPHA,omega.sim=omega.mat)
  }
}
}

###############################################
###function: transform geno input files into YY and NN
##############################################

geno2YN<-function(genofile){
 data=read.table(genofile)
 npop=ncol(data)/2 ; all1=seq(1,2*npop,2) ; all2=all1+1
 YY=data[,all1]
 NN=YY+data[,all2]
 list(YY=YY,NN=NN)
}

###############################################
###function: Compute FMD distance (Forstner and Moonen, 2003) between two covariance matrices
##############################################
fmd.dist<-function(mat1,mat2){
 require(geigen)
 return(sqrt(sum(log(geigen(mat1,mat2)$values)**2)))
}


