


####ZIGDM_EM usage

####library(DMBC)
####library(Rlab)

#' @export
#'


####First step: generate pi
Bernoulli <- function(c.est, row) {

  #pi = exp(c.est %*% row)/(1+exp(c.est %*% row))
  pi = 0
  Delta = rbern(dim(c.est)[1], pi)
  return (Delta)
}


####Second step: generate Z
Zvalue_Generator <- function(Delta,Test_result, row2, row3) {

  Beta_mean = exp(Test_result$alpha.est %*% row2)/(1+exp(Test_result$alpha.est %*% row2))
  Beta_dispersion = exp(Test_result$beta.est %*% row3)/(1+exp(Test_result$beta.est %*% row3))
  Data_Matrix = cbind(as.matrix(Delta),Beta_mean,Beta_dispersion)

  test = apply(Data_Matrix, 1, RandomGenerate_beta)
  return(test)
}

RandomGenerate_beta <- function(input) {
    delta = input[1]
    beta_mean =input[2]
    beta_dispersion = input[3]

   if (delta ==1) {
      Z = 0
    }else{
      shape1= beta_mean*(1/beta_dispersion)-beta_mean
      shape2 = (1-beta_mean)*(1/beta_dispersion -1)
      Z = rbeta(1, shape1, shape2, ncp = 0)
    }
  return(Z)
  }

####Third step: generate p for multinormail distribution
multinormial_P <- function(Z) {
  P= list()
  P[1] = Z[1]
  if(length(Z)>1) {
    for (j in 2:length(Z)) {
      multiply =1
      for(k in 1:(j-1)) {
        multiply = multiply*(1-Z[k])
      }
      P[j] = Z[j]*multiply
    }
  }
  P= unlist(P)
  return(P)
}






################ Run below functions before calling
## The function ZIGDM.EM.PAR2 estimates three sets of regression coefficients "c", "alpha", and "beta" in the ZIGDM regression model.

## Input
# Y: taxa count matrix.
# Xc, Xa, Xb: design matrices for probability of absence model, mean model and dispersion model, respectively. The first column should contain "1" for the intercept. For intercept only model, the three matrices are simply contain vectors of 1.
# c0: initial values for regression coefficients "c" organized in a matrix (dim: K x d, K+1 is the #taxa (i.e. K+1 = ncol(Y)), d is the number of columns in the corresponding design matrix)
# alpha0: initial values for regression coefficients "alpha" organized in a matrix (dim: K x d)
# beta0: initial values for regression coefficients "beta" organized in a matrix (dim: K x d)
# tol: convergence tolerance
# max.iter: maximum number of iterations

## output
# c.est: regression coefficient for the probability of absence organized in a matrix. The probabilities of absence for all the taxa in sample i can then be calculated as exp(c.est %*% Xc[i,])/(1+exp(c.est %*% Xc[i,]))
# alpha.est: regression coefficient for the Beta mean parameter organized in a matrix. The Beta means for all the taxa in sample i  can then be calculated as exp(alpha.est %*% Xa[i,])/(1+exp(alpha.est %*% Xa[i,]))
# beta.est: regression coefficient for the Beta dispersion parameter organized in a matrix. The Beta dispersions for all the taxa in sample i can then be calculated as exp(beta.est %*% Xb[i,])/(1+exp(beta.est %*% Xb[i,]))


## Two important remarks:
## 1. Sort the taxa as Y = Y[,order( colSums(Y), decreasing = TRUE )] before run the analysis
## 2. If you don't want to have zero-inflation on a certain subset of taxa, just set the elements of the corresponding columns in c0 to -Inf.
##    So to fit the GDM model, all elements in c0 matrix should be set to -Inf

ZIGDM.EM.PAR2 <- function(Y, Xc, Xa, Xb, c0, alpha0, beta0, tol=0.0001, max.iter=1000){
  Y = Y[,order( colSums(Y), decreasing = TRUE )]
  CONV = 0
  CONV.iter = max.iter
  n = nrow(Y)
  K = ncol(Y)-1
  dc = ncol(Xc)
  da = ncol(Xa)
  db = ncol(Xb)

  c.last = c0; alpha.last = alpha0; beta.last = beta0;
  c.now = c0; alpha.now = alpha0; beta.now = beta0;

  Del.R = matrix(NA, n, K)
  A.R = matrix(NA, n, K)
  B.R = matrix(NA, n, K)

  for(l in 1:max.iter){

    # E-step

    # print(paste("====== ", l, "th ======", sep=""))
    for(i in 1:n){

      tmp = exp(c.last %*% Xc[i,])
      tmp[is.na(tmp)] = 0
      pv = as.numeric( tmp/(1+tmp) )
      pv[is.infinite(tmp) & tmp>0] =1  # positive inf

      tmp = exp(alpha.last %*% Xa[i,])
      mv = as.numeric( tmp/(1+tmp) )
      tmp = exp(beta.last %*% Xb[i,])
      sv = as.numeric( tmp/(1+tmp) )

      av = (1/sv - 1) * mv
      bv = (1/sv - 1) * (1-mv)

      par.post = rP.Y.par(pv, av, bv, Y[i,])

      Del.R[i, ] = par.post$pv.post
      tmp = par.post$av.post + par.post$bv.post
      #A.R[i, ] = digamma(par.post$av.post) - digamma(tmp)
      #B.R[i, ] = digamma(par.post$bv.post) - digamma(tmp)

      A.R[i, ] = as.matrix(digamma(par.post$av.post) - digamma(tmp))
      B.R[i, ] = as.matrix(digamma(par.post$bv.post) - digamma(tmp))

    }



    # M-step
    for(j in 1:K){

      if( !is.infinite(c.last[j,1]) ){
        c.now[j, ] = Logit.optim(Del.R[,j], Xc, c.last[j, ])
      }

      tmp = Beta.optim.PAR2(Del.R[,j], A.R[,j], B.R[,j], Xa, Xb, alpha.last[j, ], beta.last[j, ])
      alpha.now[j, ] = tmp[1:da]; beta.now[j, ] = tmp[-(1:da)]

    }


    diff = 0
    if( sum(!is.infinite(c.now))>0 ){

      diff = diff + sum(abs(c.now[!is.infinite(c.now)]-c.last[!is.infinite(c.now)]))

    }

    diff = diff + sum(abs(alpha.now-alpha.last))
    diff = diff + sum(abs(beta.now-beta.last))
    if( diff < tol ){

      CONV = 1; CONV.iter = l;
      break;

    }else{

      c.last = c.now; alpha.last = alpha.now; beta.last = beta.now;
    }

  }

  return(list(c.est = c.now, alpha.est = alpha.now, beta.est = beta.now, CONV=CONV, CONV.iter = CONV.iter))

}



rP.Y.par <- function(pv, av, bv, Y){

  K = length(Y)-1
  N = sum(Y)

  pv.post = rep(0, K)

  av.prim = av + Y[1:K]
  bv.prim = bv
  for(j in 1:K){
    bv.prim[j] = bv.prim[j] + (N - sum(Y[1:j]))

    if(Y[j]==0){

      if(beta(av[j], bv[j])==0){

        pv.post[j] = pv[j]/(pv[j] + (1-pv[j]))

      }else{
        tmp = pv[j] + (1-pv[j])*beta(as.numeric(av.prim[j]), bv.prim[j])/beta(av[j], bv[j])
        #tmp = pv[j] + (1-pv[j])*beta(av.prim[j], bv.prim[j])/beta(av[j], bv[j])
        if(tmp==0){
          pv.post[j] = 1
        }else{
          pv.post[j] = pv[j]/tmp
        }

      }


    }

  }

  return(list(pv.post = pv.post, av.post = av.prim, bv.post = bv.prim))
}

lbeta2 <- function(a,b){
  a[a<0] = 0
  b[b<0] = 0
  return(lbeta(a,b))
}

Logit.neg.loglik <- function(par, data){

  tmp = as.numeric( exp(  data$X %*% par ) )
  p = tmp/(1+tmp)
  tmp = data$Del * log(p) + (1-data$Del) * log(1-p)
  index = which(p==0 | p==1)
  if(length(index)>0){
    tmp[index] = 0
  }

  return( -sum( tmp) )

}

Logit.neg.score <- function(par, data){

  tmp = as.numeric( exp(  data$X %*% par ) )
  p = tmp/(1+tmp)
  return( -colSums( (data$Del - p) * data$X ) )

}

Logit.optim <- function(Del, X, c.ini){

  Logit.par.ini = c.ini
  Logit.data = list(Del=Del, X=X)

  return( optim(par=Logit.par.ini, fn=Logit.neg.loglik, gr=Logit.neg.score, data = Logit.data, method="BFGS")$par)


}



Beta.neg.loglik.PAR2 <- function(par, data){

  da = ncol(data$Xa)
  alpha = par[1:da]
  beta = par[-(1:da)]


  tmp = as.numeric( exp( data$Xa %*% alpha) )
  mu.tmp = as.numeric( tmp/(1+tmp) )

  tmp = as.numeric( exp( data$Xb %*% beta) )
  sigma.tmp = as.numeric( tmp/(1+tmp) )

  a = (1/sigma.tmp - 1) * mu.tmp
  b = (1/sigma.tmp - 1) * (1-mu.tmp)
  a[a<0] = 0
  b[b<0] = 0

  return( -sum( (1-data$Del) * ( -lbeta2(a,b) + data$A * (a-1) + data$B *(b-1) )  ) )

}

Beta.neg.score.PAR2 <- function(par, data){

  da = ncol(data$Xa)
  alpha = par[1:da]
  beta = par[-(1:da)]

  tmp.a = as.numeric( exp( data$Xa %*% alpha) )
  mu.tmp = as.numeric( tmp.a/(1+tmp.a) )

  tmp.b = as.numeric( exp( data$Xb %*% beta) )
  sigma.tmp = as.numeric( tmp.b/(1+tmp.b) )

  a = (1/sigma.tmp - 1) * mu.tmp
  b = (1/sigma.tmp - 1) * (1-mu.tmp)
  a[a<0] = 0
  b[b<0] = 0

  a.a = (1/tmp.b) * mu.tmp * (1/(1+tmp.a))
  a.b = -(1/tmp.b) * mu.tmp
  b.a = - a.a
  b.b = -(1/tmp.b) *  (1/(1+tmp.a))

  one = (1-data$Del)*(digamma(a+b)-digamma(a) + data$A)
  two = (1-data$Del)*(digamma(a+b)-digamma(b) + data$B)

  return( -c( colSums( (one*a.a + two*b.a) * data$Xa ), colSums( (one*a.b + two*b.b) * data$Xb ) ) )
}

Beta.optim.PAR2 <- function(Del, A, B, Xa, Xb, alpha.ini, beta.ini){

  Beta.par.ini = c(alpha.ini, beta.ini)
  Beta.data = list(Del=Del, A=A, B=B, Xa=Xa, Xb=Xb)

  return( optim(par=Beta.par.ini, fn=Beta.neg.loglik.PAR2, gr=Beta.neg.score.PAR2, data = Beta.data, method="BFGS")$par)

}

