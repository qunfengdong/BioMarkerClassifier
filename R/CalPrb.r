#' @title Calculate Likelihood based on Dirichlet-multinomial distribution estimated parameters
#'
#' @description
#' This function estimates parameters from Dirichlet-multinomial distribution.
#' @param FS_out An object from the FeatureSelection()
#' @param testSet A test set in data frame or matrix form. The colnames should have the same bacteria (features) as in the training set.
#' @param col_start An index indicating at which column is the beginning of bacteria (features) data. The default is the 3rd column.
#' @param type_col An index indicating at which column is group/type variable. The default is the 2nd column.
#' @param HighestRank The top number of features inclueded in model. The default is all the features left after filtering.
#' @return A data frame with 17 columns, each row represents a model estimation output.
#' @export
#' @examples
#' data(training)
#'
#' #### Take one row as testSet ####
#' idx <- sample(1:nrow(training),1)
#' test <- training[idx,]
#' train <- training[-idx,]
#'
#' CalPrb(FS(train),test) # This may take up to one minute

CalPrb <- function(FS_out=FS_out,testSet=testSet,col_start=3,type_col=2,HighestRank=nrow(FS_out$Feature)){

  Genus <- FS_out$CountData
  SortP <- FS_out$Feature

  Disease <- levels(Genus[,type_col])

  TotalType1 <- nrow(Genus[Genus[,type_col] == Disease[1], ])
  TotalType2 <- nrow(Genus[Genus[,type_col] == Disease[2], ])
  Prior1 <- TotalType1/(TotalType1+TotalType2)
  Prior2 <- TotalType2/(TotalType1+TotalType2)

  lh <- list()
  for(rank in 2:HighestRank) {
  #for(rank in 2:5) {
print(rank)
    ######## choose the signature taxa and merge the training data

    #NameList = as.vector(P_table[as.numeric(as.character(P_table$P_Wilcoxon))<0.2,]$Genera)
    #select the NameList based on the rank
    NameList <- as.vector(rownames(SortP)[1:rank])
    NewDF <- Genus[,colnames(Genus) %in% NameList]
    NewDF2 <- Genus[,!colnames(Genus)%in% NameList]

    col_end2 <- dim(NewDF2)[2]
    NewDF2$Others <- rowSums(NewDF2[, col_start:col_end2])
    #NewDFTotal = data.frame(NewDF2[, 1:(col_start-1)],NewDF, NewDF2$Others)
    NewDFTotal <- data.frame(NewDF2[, 1:(col_start-1)],NewDF)
    ###### merge in test row

    NewTestSignature <- testSet[, colnames(Genus) %in% NameList]
    NewTestNonSignature <- testSet[, !colnames(Genus)%in% NameList]
    NewTestNonSignature$Others <- rowSums(NewTestNonSignature[, col_start:col_end2])
    NewTestTotal <- data.frame(NewTestNonSignature[, 1:(col_start-1)],NewTestSignature)
    #NewTestTotal = data.frame(NewTestNonSignature[, 1:(col_start-1)],NewTestSignature, NewTestNonSignature$Others)


    rep2_Type1 <- NewDFTotal[NewDFTotal[,type_col] ==Disease[1],]
    rep2_Type2 <- NewDFTotal[NewDFTotal[,type_col] ==Disease[2],]


    ############Estimate Zero Inflated Generalized Dirchlet-Multinomial parameters
    tol <-0.0001
    max.iter<-1000

     ###### Estimate the parameters from the control data
    Xc <- as.matrix(rep(1, dim(rep2_Type1)[1]))
    Xa <- Xc
    Xb <- Xc
    c0 <- as.matrix(rep(1, dim(rep2_Type1)[2]-col_start))
    alpha0 <- c0
    beta0 <- c0


    ZIGDM_est1 <-ZIGDM.EM.PAR2(rep2_Type1[,-(1:(col_start-1))], Xc, Xa, Xb, c0, alpha0, beta0, tol=0.0001, max.iter=1000)

 #   ZIGDM_est1 <-ZIGDM.EM.PAR2(rep2_Type1[,-c((1:(col_start-1)), ncol(rep2_Type1))], Xc, Xa, Xb, c0, alpha0, beta0, tol=0.0001, max.iter=1000)
    ###p_Type1 <- ZIGDM_prob(rep2_Type1[,-(1:(col_start-1))], Xc, Xa, Xb, c0, alpha0, beta0, tol=0.0001, max.iter=1000)
    #fit3 <- dirmult(rep2_Type1[,-(1:(col_start-1))],epsilon=10^(-4),trace=FALSE)

    ###### Estimate the paramenters from the baseline data
    Xc <- as.matrix(rep(1, nrow(rep2_Type2)))
    Xa <- Xc
    Xb <- Xc
    c0 <- as.matrix(rep(1, ncol(rep2_Type1)-col_start))
    alpha0 <- c0
    beta0 <- c0

    ZIGDM_est2<-ZIGDM.EM.PAR2(rep2_Type2[,-(1:(col_start-1))], Xc, Xa, Xb, c0, alpha0, beta0, tol=0.0001, max.iter=1000)

    ##p_Type2 <- ZIGDM_prob(rep2_Type2[,-(1:(col_start-1))], Xc, Xa, Xb, c0, alpha0, beta0, tol=0.0001, max.iter=1000)
    #fit4 <- dirmult(rep2_Type2[,-(1:(col_start-1))],epsilon=10^(-4),trace=FALSE)


    #########Calculate the likelihood of the test sample being Type1 and Type2
######error p_Type1 has 7 number as 7 taxa? while only 2 taxa in
    #p_Type1 <- ZIGDM_prob(ZIGDM_est1, row = Xc[1,])
    #p_Type2 <- ZIGDM_prob(ZIGDM_est2, row = Xc[1,])

    testlist <- list()

    for (r in 1:nrow(testSet)){
      repeatlist <- list()
      for (rep in 1:10) {
        ##### Calculate zero inflated generalized Dirichlet multinomial probability mass function P(x|Type1)
        p_Type1 <- ZIGDM_prob(ZIGDM_est1, row = Xc[1,])
        p_Type2 <- ZIGDM_prob(ZIGDM_est2, row = Xc[1,])
        normalizeTest = NewTestTotal[r,-(1:(col_start-1))]/sum(NewTestTotal[r,-(1:(col_start-1))])*100
        if(length(normalizeTest) == 2) {
          probability_Type1 = apply(round(normalizeTest[1]),1,dbinom, 100, prob = p_Type1[1])
        }else{
        probability_Type1 = apply(normalizeTest,1,dmultinom, prob = p_Type1)
        }
        lhP_Type1 = probability_Type1*Prior1

        #pdfln_Type1 <- ddirmn(NewTestTotal[r,-(1:(col_start-1))], t(as.matrix(p_Type1)))
        #lh_Type1=exp(pdfln_Type1)
        #lhP_Type1 = lh_Type1*Prior1

        ##### Calculate log of Dirichlet multinomial probability mass function P(x|Type2)
        if(length(normalizeTest) == 2) {
          probability_Type2 = apply(round(normalizeTest[1]),1,dbinom, 100, prob = p_Type2[1])
        }else{
          probability_Type2 = apply(normalizeTest,1,dmultinom, prob = p_Type2)
        }

        lhP_Type2 = probability_Type2*Prior2

        #pdfln_Type2 <- ddirmn(NewTestTotal[r,-(1:(col_start-1))], t(as.matrix(p_Type2)))
        #lh_Type2=exp(pdfln_Type2)
        #lhP_Type2 = lh_Type2*Prior2

        #Pos_Type1 = lh_Type1/(lh_Type1+lh_Type2)
        #Pos_Type2 = lh_Type2/(lh_Type1+lh_Type2)

        PosP_Type1 = lhP_Type1/(lhP_Type1+lhP_Type2)
        PosP_Type2 = lhP_Type2/(lhP_Type1+lhP_Type2)

        #### Create truth labels ####

        if(testSet[r,type_col] == Disease[1]) {
          Type1Label = 1
          Type2Label = 0
        } else if (testSet[r,type_col] == Disease[2]) {
          Type1Label = 0
          Type2Label = 1
        }
      repeatlist[[rep]] <- as.matrix(t(c(Disease[1], Disease[2], rownames(testSet)[r], rank,  PosP_Type1, PosP_Type2, Type1Label, Type2Label,paste(NameList,collapse=";"))))

    }
      RandomRepeat = data.frame(t(sapply(repeatlist,'[')))
      PosP_Type1 = mean(as.numeric(as.character(RandomRepeat$X5)), na.rm = TRUE)
      PosP_Type2 = mean(as.numeric(as.character(RandomRepeat$X6)), na.rm = TRUE)

      testlist[[r]] <- as.matrix(t(c(Disease[1], Disease[2], rownames(testSet)[r], rank,  PosP_Type1, PosP_Type2, Type1Label, Type2Label,paste(NameList,collapse=";"))))




##    lh[[rank]] <- as.matrix(t(c(Disease[1], Disease[2], rownames(testSet), rank, Prior1, Prior2, lh_Type1, lh_Type2, Pos_Type1, Pos_Type2, lhP_Type1, lhP_Type2, PosP_Type1, PosP_Type2, Type1Label, Type2Label,paste(NameList,collapse=";"))))

    lh[[rank]] <- data.frame(t(sapply(testlist,'[')))

  } #end of rank

  lh_table <- data.frame(do.call(rbind,lh))
##  colnames(lh_table) = c("Type1", "Type2", "row", "feature_rank", "Prior1", "Prior2", "lh_Type1", "lh_Type2", "Poster_Type1", "Poster_Type2", "lhP_Type1", "lhP_Type2","Poster_Prio_Type1","Poster_Prio_Type2", "Type1Label", "Type2Label","SelectedFeatures" )

  colnames(lh_table) = c("Type1", "Type2", "row", "feature_rank", "Type1_Posterior_Prb","Type2_Postereior_Prb", "Type1Label", "Type2Label","SelectedFeatures" )
}
  return(lh_table)
}
