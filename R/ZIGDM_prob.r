#' @title Calculate probability based on zero inflated generalized Dirichlet-multinomial distribution
#'
#' @description
#' This function estimates parameters from zero inflated generalized Dirichlet-multinomial distribution, and test the pro.
#' @param Y Count table of microbiome
#' @param Xc design matrices for probability of absence model. The first column should contain "1" for the intercept. For intercept only model, the three matrices are simply contain vectors of 1.
#' @param Xa design matrices for mean model.  The first column should contain "1" for the intercept. For intercept only model, the three matrices are simply contain vectors of 1.
#' @param Xb design matrices for dispersion model. The first column should contain "1" for the intercept. For intercept only model, the three matrices are simply contain vectors of 1.
#' @param c0 c0: initial values for regression coefficients "c" organized in a matrix (dim: K x d, K+1 is the #taxa (i.e. K+1 = ncol(Y)), d is the number of columns in the corresponding design matrix)
#' @param alpha0: initial values for regression coefficients "alpha" organized in a matrix (dim: K x d)
#' @param beta0: initial values for regression coefficients "beta" organized in a matrix (dim: K x d)
#' @param tol: convergence tolerance
#' @param max.iter: maximum number of iterations

#' @return A vector containning the probablity of multinomial distribution
#' @export
#' @examples
####estimate probability of of absent matrix, mean and dispersion matrix from Y based on ZIGDM model
###Test_result <- ZIGDM.EM.PAR2(Y=Y, Xc=Xc, Xa=Xa, Xb=Xb, c0=c0, alpha0=alpha0, beta0=beta0, tol=tol, max.iter=max.iter)
#### using probability of of absent matrix to random generate present and absent value


ZIGDM_prob <- function(Test_result, row) {
 ####row is Xc[1,]
  Delta <- Bernoulli(c.est = Test_result$c.est, row = row)
  #### based on Delta, generate Z, row2 is Xa[1,], row3 is Xb[1,]. since Xc[1,]=Xa[1,]=Xb[1,],row is used for all
  Z_value <- Zvalue_Generator(Delta = Delta, Test_result = Test_result, row2 = row, row3 = row)
  #### based on Z, generate p vector for multinomial

  MultinomialP <- multinormial_P(Z=Z_value)
  MultinomialP[length(MultinomialP)+1] <- 1- sum(MultinomialP)  #### add the K+1 taxa back

  return (MultinomialP)

  #probability = apply(Y,1,dmultinom, prob = MultinomialP)

  #return (probability)
}

