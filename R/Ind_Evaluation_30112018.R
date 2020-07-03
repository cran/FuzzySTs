##### Function calculating the evaluation per observation with non-response #####
#################################################################################

# Non response functions #
##########################
#' Calculates the factor delta_jki
#' @param x a particular observation.
#' @return A binary value 0 or 1.
# #' @export

delta_jki <- function(x){
  delta_jki <- !is.na(x)
  return(delta_jki)
} 

#' Calculates the factor Delta_jki
#' @param x a dataset.
#' @param i an observation index.
#' @param K the total number of linguistics in a sub-item.
#' @return The response matrix of binary values (0 or 1) related to the answers of a particular dataset for its corresponding sub-items.
# #' @export

Delta_jki <- function(x, i, K){
  
  # x is the database
  # i is the observation index
  # K is the total number of linguistics in a sub-item
  Delta_jki=matrix(rep(0,K))
  
  for  (k in 1:K)
  {
    Delta_jki[k] = delta_jki(x[i, k])
  }
  return(Delta_jki[,])
} 



#' Calculates the adjusted weight for a given sub-item of a linguistic questionnaire
#' @param x the data set to evaluate.
#' @param i an observation index.
#' @param k a sub-item index.
#' @param b_jk an array referring to the initial weights given to each sub-item of the considered main-item. This array will be afterwards re-calculated.
#' @return A numerical value giving the readjusted weight of the sub-item k of the considered main-item for the observation i.
#' @export
#' @examples data <- matrix(c(3,4,2,3,3,2,4,3,3,4,3,4,4,2,5,3,4,4,3,3,3,4,4,3,
#' 3,3,4,3,3,3,3,4,4,3,5,3,4,3,3,3), ncol = 4)
#' adjusted.weight.SI(data, 7, 1, c(0.5,0.5))
adjusted.weight.SI  <- function(x, i, k, b_jk){
  
  # adjusted.weight.SI is the adjusted weight for a certain sub-item k 
  # x is the database
  # i is the observation index
  # k is the sub-item index
  # b_jk is the sub-items initial weights fixed for a certain main-item j
  # b_jk is a form array b_jk=c(w_j1,...,w_jk,...w_j(m_j)) where w_j1, ... , w_jk, ..., w_j(m_j) are the weights of sub-item 1, ... , k, ..., r etc respectively)
  adjusted.weight.SI <- ( delta_jki(x[i,k]) * b_jk[k] ) / sum ( Delta_jki(x, i, length(b_jk)) * b_jk)
  
  adjusted.weight.SI[is.na(adjusted.weight.SI)] <- 0
  
  return(adjusted.weight.SI)
  
}



#' Calculates the adjusted weight for a given main-item of a linguistic questionnaire
#' @param x the data set to evaluate.
#' @param i an observation index.
#' @param j a main-item index.
#' @param b_j an array referring to the initial weights given to each main-item of the considered main-item. This array will be afterwards re-calculated.
#' @param b_jk a matrix of length(b_j) rows and max(SI) columns expressing the initial weights of each sub-item of a given main-item.
#' @param SI an array representing the total numbers of sub-items per main-item.
#' @return A numerical value giving the readjusted weight of the main-item j for the observation i.
#' @export
#' @examples data <- matrix(c(3,4,2,3,3,2,4,3,3,4,3,4,4,2,5,3,4,4,3,3,3,4,4,3,
#' 3,3,4,3,3,3,3,4,4,3,5,3,4,3,3,3), ncol = 4)
#' data <- as.data.frame(data)
#' MI <- 2
#' SI1 <- 2
#' SI2 <- 2
#' SI <- c(SI1,SI2)
#' b_j <- c(1/2,1/2)
#' b_jk <- matrix(c(0.5,0.5,0.5,0.5),nrow=2) 
#' PA11 <- c(1,2,3,4,5)
#' PA12 <- c(1,2,3,4,5)
#' PA21 <- c(1,2,3,4,5)
#' PA22 <- c(1,2,3,4,5)
#' # ------------------
#' MF111 <- TrapezoidalFuzzyNumber(0,2,2,7)
#' MF112 <- TrapezoidalFuzzyNumber(2,7,7,15)
#' MF113 <- TrapezoidalFuzzyNumber(7,15,15,23)
#' MF114 <- TrapezoidalFuzzyNumber(15,23,23,28)
#' MF115 <- TrapezoidalFuzzyNumber(23,28,28,30)
#' MF11 <- GFUZZ(data, 1, 1, PA11, spec="Identical", breakpoints = 100)
#' # ------------------
#' MF121 <- TrapezoidalFuzzyNumber(0,2,2,7)
#' MF122 <- TrapezoidalFuzzyNumber(2,7,7,15)
#' MF123 <- TrapezoidalFuzzyNumber(7,15,15,23)
#' MF124 <- TrapezoidalFuzzyNumber(15,23,23,28)
#' MF125 <- TrapezoidalFuzzyNumber(23,28,28,30)
#' MF12 <- GFUZZ(data, 1, 2, PA12, spec="Identical", breakpoints = 100)
#' # ------------------
#' MF211 <- TrapezoidalFuzzyNumber(0,2,2,7)
#' MF212 <- TrapezoidalFuzzyNumber(2,7,7,15)
#' MF213 <- TrapezoidalFuzzyNumber(7,15,15,23)
#' MF214 <- TrapezoidalFuzzyNumber(15,23,23,28)
#' MF215 <- TrapezoidalFuzzyNumber(23,28,28,30)
#' MF21 <- GFUZZ(data, 2, 1, PA21, spec="Identical", breakpoints = 100)
#' # ------------------
#' MF221 <- TrapezoidalFuzzyNumber(0,2,2,7)
#' MF222 <- TrapezoidalFuzzyNumber(2,7,7,15)
#' MF223 <- TrapezoidalFuzzyNumber(7,15,15,23)
#' MF224 <- TrapezoidalFuzzyNumber(15,23,23,28)
#' MF225 <- TrapezoidalFuzzyNumber(23,28,28,30)
#' MF22 <- GFUZZ(data, 2, 2, PA22, spec="Identical", breakpoints = 100)
#' # ------------------
#' range <- matrix(c(0,0,0,0,28,28,28,28), ncol=2)
#' adjusted.weight.MI(data, 9, 1, b_j, b_jk, SI)

adjusted.weight.MI  <- function(x, i, j, b_j, b_jk, SI){
  
  # adjusted.weight.MI is the adjusted weight for a certain main-item j 
  # x is the database
  # i is the observation index
  # j is the main-item index
  # b_j is the main-items initial weights fixed (b_j is a form array b_j=c(w_1,...,w_j,...w_r) where w_1, ... , w_j, ..., w_r are the weights of main-item 1, ... , j, ..., r etc respectively)
  # b_jk is the sub-items initial weights fixed 
  # b_jk is a matrix of j rows and k columns b_jk=matrix(c(w_1,...,w_k,...w_(m_j)) where w_1, ... , w_k, ..., w_(m_j) are the weights of sub-item 1, ... , k, ..., m_j etc respectively)
  # b_jk = b_11   b_12   ...    b_1k   ...    b_1(m_j) 
  #        b_21   b_22   ...    b_2k   ...    b_2(m_j) 
  #         .       .    ...      .    ...       . 
  #         .       .    ...      .    ...       . 
  #         .       .    ...      .    ...       . 
  #        b_j1   b_j2   ...    b_jk   ...    b_j(m_j) 
  #         .       .    ...      .    ...       . 
  #         .       .    ...      .    ...       . 
  #         .       .    ...      .    ...       . 
  #        b_r1   b_r2   ...    b_rk   ...    b_r(m_j) 
  # N.B: If, for example, case of 2 main-items don't have the same number of sub-items:
  #             main_item 1 : 5 sub-items
  #             main_item 2 : 3 sub-items
  # b_jk will be
  # b_jk = b_11   b_12   b_13    b_14   b_15 
  #        b_21   b_22   b_23     0      0 
  # SI is the total numbers of sub-items per main-item (used to decompose the database)
  # (SI is an array SI=c(k_1,...,k_j,...k_r) where k_1, ... , k_j, ..., k_r are the total numbers of sub-items per main-item 1, ... , j, ..., r etc respectively)
  
  cum.SI <- c(0,cumsum(SI))
  
  sum.mi <- c(rep(0,length(b_j)))
  for (u in 1:length(b_j))
  {
    x.mi <- cbind(x[,(cum.SI[u]+1):cum.SI[u+1]])
    sum.mi[u] <- sum ( adjusted.weight.SI(x.mi, i, 1:ncol(x.mi), b_jk[u,1:ncol(x.mi)]))  * b_j[u]
  }
  
  xj <- cbind(x[,(cum.SI[j]+1):cum.SI[j+1]])
  adjusted.weight.MI <- ( sum ( adjusted.weight.SI(xj, i, 1:length(xj), b_jk[j,1:ncol(xj)]) ) * b_j[j] ) / sum ( sum.mi )
  
  adjusted.weight.MI[is.na(adjusted.weight.MI)] <- 0
  
  return(adjusted.weight.MI)
  
} 

#' Calculates the indicator of information's rate of the data base for a given unit
#' @param x the data set to evaluate.
#' @param i an observation index.
#' @param b_jk a matrix of length(b_j) rows and max(SI) columns expressing the initial weights of each sub-item of a given main-item. 
#' @param SI an array representing the total numbers of sub-items per main-item.
#' @return A numerical value giving the indicator of information's rate of the complete linguistic questionnaire for a particular observation. Note that the obtained value is interpreted as the more it tends to the value 1, the less the observation i contains missing values.
#' @export
#' @examples data <- matrix(c(3,4,2,3,3,2,4,3,3,4,3,4,4,2,5,3,4,4,3,3,3,4,4,3,
#' 3,3,4,3,3,3,3,4,4,3,5,3,4,3,3,3), ncol = 4)
#' data <- as.data.frame(data)
#' SI1 <- 2
#' SI2 <- 2
#' SI <- c(SI1,SI2)
#' b_jk <- matrix(c(0.5,0.5,0.5,0.5),nrow=2) 
#' Ri(data, 7, b_jk, SI)

Ri <- function(x, i, b_jk, SI){
  
  # Ri is the pertinence indicator of an individual questionnaire regarding the weighted number of missing values in his answers, Ri is between 0 and 1, 1 means no missing values and 0 the opposite.
  # x is the database
  # i is the observation index
  # b_jk is the sub-items initial weights fixed 
  # b_jk is a matrix of j rows and k columns b_jk=matrix(c(w_1,...,w_k,...w_(m_j)) where w_1, ... , w_k, ..., w_(m_j) are the weights of sub-item 1, ... , k, ..., m_j etc respectively)
  # b_jk = b_11   b_12   ...    b_1k   ...    b_1(m_j) 
  #        b_21   b_22   ...    b_2k   ...    b_2(m_j) 
  #         .       .    ...      .    ...       . 
  #         .       .    ...      .    ...       . 
  #         .       .    ...      .    ...       . 
  #        b_j1   b_j2   ...    b_jk   ...    b_j(m_j) 
  #         .       .    ...      .    ...       . 
  #         .       .    ...      .    ...       . 
  #         .       .    ...      .    ...       . 
  #        b_r1   b_r2   ...    b_rk   ...    b_r(m_j) 
  # N.B: If, for example, case of 2 main-items don't have the same number of sub-items:
  #             main_item 1 : 5 sub-items
  #             main_item 2 : 3 sub-items
  # b_jk will be
  # b_jk = b_11   b_12   b_13    b_14   b_15 
  #        b_21   b_22   b_23     0      0 
  # SI is the total numbers of sub-items per main-item (used to decompose the database)
  # (SI is an array SI=c(k_1,...,k_j,...k_r) where k_1, ... , k_j, ..., k_r are the total numbers of sub-items per main-item 1, ... , j, ..., r etc respectively)
  
  
  num.NA <- c(rep(0,length(SI)))
  S1=0
  cum.SI <- c(0,cumsum(SI))
  for (m in 1:length(SI))
  {
    S1=0
    x.mi <- cbind(x[,(cum.SI[m]+1):cum.SI[m+1]])
    for (v in 1:ncol(x.mi))
    {
      S1=S1+is.na(x.mi[i,v])
    }
    num.NA[m] <- S1
  }
  
  
  prod <- c(rep(0,length(SI)))
  for (j in 1:length(SI))
  {
    xj <- cbind(x[,(cum.SI[j]+1):cum.SI[j+1]])
    diff <- c(rep(0,ncol(xj)))
    for(k in 1:SI[j])
    {
      diff[k] = (b_jk[j,k] - adjusted.weight.SI(xj, i, k, b_jk[j,1:ncol(xj)])) * is.na(xj[i,k])
    }
    prod[j] <- (num.NA[j] / SI[j]) * sum(diff)
  }
  
  Ri <- 1 - (sum(prod) / sum(b_jk))
  
  return(Ri)
  
}


#' Calculates the indicator of information's rate of the data base
#' @param x the data set to evaluate.
#' @param p_ind a vector of the relative sampling weights of the units, for which \eqn{length(p_ind) = nrow(data)}. If the weights are not relative, the following expression should be applied on the vector: \deqn{\frac{p_{ind}}{\sum_{i=1}^{n} p_{ind}}.} If no sampling weights are used, the vector of weights is reduced to a vector of values 1, i.e. \eqn{rep(1, nrow(data))}.
#' @param b_jk a matrix of length(b_j) rows and max(SI) columns expressing the initial weights of each sub-item of a given main-item.
#' @param SI an array representing the total numbers of sub-items per main-item.
#' @return A numerical value giving the indicator of information's rate of the complete linguistic questionnaire. Note that the obtained value is interpreted as the more it tends to the value 1, the less the complete questionnaire contains missing values.
#' @export
#' @examples data <- matrix(c(3,4,2,3,3,2,4,3,3,4,3,4,4,2,5,3,4,4,3,3,3,4,4,3,
#' 3,3,4,3,3,3,3,4,4,3,5,3,4,3,3,3), ncol = 4)
#' data <- as.data.frame(data)
#' p_ind <- c(0.1,0.05,0.05,0.2,0.1,0.05,0.1,0.1,0.2,0.05)
#' SI1 <- 2
#' SI2 <- 2
#' SI <- c(SI1,SI2)
#' b_jk <- matrix(c(0.5,0.5,0.5,0.5),nrow=2) 
#' R(data, p_ind, b_jk, SI)

R <- function(x, p_ind, b_jk, SI){
  
  # R is the pertinence indicator of a full questionnaire regarding the weighted number of missing values in all individual answers, R is between 0 and 1, 1 means no missing values and 0 the opposite.
  # x is the database
  # p_ind is the vector of individual weights in the database (p_ind should be relative)
  # b_jk is the sub-items initial weights fixed 
  # b_jk is a matrix of j rows and k columns b_jk=matrix(c(w_1,...,w_k,...w_(m_j)) where w_1, ... , w_k, ..., w_(m_j) are the weights of sub-item 1, ... , k, ..., m_j etc respectively)
  # b_jk = b_11   b_12   ...    b_1k   ...    b_1(m_j) 
  #        b_21   b_22   ...    b_2k   ...    b_2(m_j) 
  #         .       .    ...      .    ...       . 
  #         .       .    ...      .    ...       . 
  #         .       .    ...      .    ...       . 
  #        b_j1   b_j2   ...    b_jk   ...    b_j(m_j) 
  #         .       .    ...      .    ...       . 
  #         .       .    ...      .    ...       . 
  #         .       .    ...      .    ...       . 
  #        b_r1   b_r2   ...    b_rk   ...    b_r(m_j) 
  # N.B: If, for example, case of 2 main-items don't have the same number of sub-items:
  #             main_item 1 : 5 sub-items
  #             main_item 2 : 3 sub-items
  # b_jk will be
  # b_jk = b_11   b_12   b_13    b_14   b_15 
  #        b_21   b_22   b_23     0      0 
  # SI is the total numbers of sub-items per main-item (used to decompose the database)
  # (SI is an array SI=c(k_1,...,k_j,...k_r) where k_1, ... , k_j, ..., k_r are the total numbers of sub-items per main-item 1, ... , j, ..., r etc respectively)
  
  
  num.NA <- matrix(c(rep(0)),nrow(x), length(SI))
  S1=0
  cum.SI <- c(0,cumsum(SI))
  for (i in 1:nrow(x))
  {
    for (m in 1:length(SI))
    {
      S1=0
      x.mi <- cbind(x[,(cum.SI[m]+1):cum.SI[m+1]])
      for (v in 1:ncol(x.mi))
      {
        S1=S1+is.na(x.mi[i,v])
      }
      num.NA[i,m] <- S1
    }
    
  }
  
  
  prod <- matrix(c(rep(0)),nrow(x), length(SI))
  for (i in 1:nrow(x))
  {
    for (j in 1:length(SI))
    {
      xj <- cbind(x[,(cum.SI[j]+1):cum.SI[j+1]])
      diff <- c(rep(0,ncol(xj)))
      for(k in 1:SI[j])
      {
        diff[k] = (b_jk[j,k] - adjusted.weight.SI(xj, i, k, b_jk[j,1:ncol(xj)])) * is.na(xj[i,k])
      }
      prod[i,j] <- p_ind[i] * (num.NA[i,j] / SI[j]) * sum(diff)
    }
  }
  
  R <- 1 - (sum(prod) / sum(b_jk))
  
  return(R)
  
  # This can be written as well as the weighted mean of all the individual indicators Ri as follows:
  
  # Ri.total <- rep(NA, nrow(x))
  # for (i in 1: nrow(x)){
  #   Ri.total[i] <- Ri(x, i, b_jk, SI)
  # }
  # R <- weighted.mean(Ri.total, p_ind)
  # return(R)
}


# Individual evaluations #
##########################

#' Calculates the individual evaluations of a linguistic questionnaire
#' @param Full_Database the data set to evaluate.
#' @param MI a numerical value representing the total number of main-items dividing the linguistic questionnaire.
#' @param bmi an array referring to the initial weights of the main-items.
#' @param SI an array representing the total numbers of sub-items per main-item.
#' @param b_jkt a matrix of MI rows and max(SI) columns expressing the initial weights of each sub-item of a given main-item.
#' @param range a vector of 2 elements giving the range of definition of the produced individual evaluations. The range is usually chosen in the interval between 0 and the maximum of the support set of all the membership functions modelling the data set.
#' @param distance.type type of distance chosen from the family of distances, set by default to the signed distance. The different choices are given by: "Rho1", "Rho2", "Bertoluzza", "Rhop", "Delta.pq", "Mid/Spr", "wabl", "DSGD", "DSGD.G", "GSGD".
#' @param i parameter of the density function of the Beta distribution, fixed by default to i = 1.
#' @param j parameter of the density function of the Beta distribution, fixed by default to j = 1.
#' @param theta a numerical value between 0 and 1, representing a weighting parameter. By default, theta is fixed to 1/3 referring to the Lebesgue space. This measure is used in the calculations of the following distances: d_Bertoluzza, d_mid/spr and d_phi-wabl/ldev/rdev.
#' @param thetas a decimal value between 0 and 1, representing the weight given to the shape of the fuzzy number. By default, thetas is fixed to 1. This parameter is used in the calculations of the d_theta star and the d_GSGD distances.
#' @param p a positive integer such that 1 \eqn{\le} p < infinity, referring to the parameter of the Rho_p and Delta_pq.
#' @param q a decimal value between 0 and 1, referring to the parameter of the metric Delta_pq.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @param spec specification of the fuzzification matrix. The possible values are "Identical" and "Not Identical".
#' @return A data set of individual evaluations, for which the number of observations is exactly the same as the initial data set.
#' @export
#' @examples data <- matrix(c(3,4,2,3,3,2,4,3,3,4,3,4,4,2,5,3,4,4,3,3,3,4,4,3,
#' 3,3,4,3,3,3,3,4,4,3,5,3,4,3,3,3), ncol = 4)
#' data <- as.data.frame(data)
#' MI <- 2
#' SI1 <- 2
#' SI2 <- 2
#' SI <- c(SI1,SI2)
#' b_j <- c(1/2,1/2)
#' b_jk <- matrix(c(0.5,0.5,0.5,0.5),nrow=2) 
#' PA11 <- c(1,2,3,4,5)
#' PA12 <- c(1,2,3,4,5)
#' PA21 <- c(1,2,3,4,5)
#' PA22 <- c(1,2,3,4,5)
#' # ------------------
#' MF111 <- TrapezoidalFuzzyNumber(0,2,2,7)
#' MF112 <- TrapezoidalFuzzyNumber(2,7,7,15)
#' MF113 <- TrapezoidalFuzzyNumber(7,15,15,23)
#' MF114 <- TrapezoidalFuzzyNumber(15,23,23,28)
#' MF115 <- TrapezoidalFuzzyNumber(23,28,28,30)
#' MF11 <- GFUZZ(data, 1, 1, PA11, spec="Identical", breakpoints = 100)
#' # ------------------
#' MF121 <- TrapezoidalFuzzyNumber(0,2,2,7)
#' MF122 <- TrapezoidalFuzzyNumber(2,7,7,15)
#' MF123 <- TrapezoidalFuzzyNumber(7,15,15,23)
#' MF124 <- TrapezoidalFuzzyNumber(15,23,23,28)
#' MF125 <- TrapezoidalFuzzyNumber(23,28,28,30)
#' MF12 <- GFUZZ(data, 1, 2, PA12, spec="Identical", breakpoints = 100)
#' # ------------------
#' MF211 <- TrapezoidalFuzzyNumber(0,2,2,7)
#' MF212 <- TrapezoidalFuzzyNumber(2,7,7,15)
#' MF213 <- TrapezoidalFuzzyNumber(7,15,15,23)
#' MF214 <- TrapezoidalFuzzyNumber(15,23,23,28)
#' MF215 <- TrapezoidalFuzzyNumber(23,28,28,30)
#' MF21 <- GFUZZ(data, 2, 1, PA21, spec="Identical", breakpoints = 100)
#' # ------------------
#' MF221 <- TrapezoidalFuzzyNumber(0,2,2,7)
#' MF222 <- TrapezoidalFuzzyNumber(2,7,7,15)
#' MF223 <- TrapezoidalFuzzyNumber(7,15,15,23)
#' MF224 <- TrapezoidalFuzzyNumber(15,23,23,28)
#' MF225 <- TrapezoidalFuzzyNumber(23,28,28,30)
#' MF22 <- GFUZZ(data, 2, 2, PA22, spec="Identical", breakpoints = 100)
#' # ------------------
#' range <- matrix(c(0,0,0,0,28,28,28,28), ncol=2)
#' ind.eval <- IND.EVAL(data,MI,b_j,SI,b_jk, range = range, distance.type ="DSGD.G")

IND.EVAL <- function(Full_Database,MI,bmi,SI,b_jkt, range, distance.type, i=1, j=1, theta = 1/3, thetas=1, p=2, q=0.5, breakpoints=100, spec="Identical"){

SI=t(t(SI))
SSI=0
crisp_data_MI=matrix(as.numeric(NA), nrow(Full_Database), length(SI))
crisp_data_OBS_NA=NULL
n_obs <- nrow(Full_Database)
bsmi <- matrix(rep(0),nrow=n_obs,ncol=MI)

for (mi in 1:MI)
{
  
  MI_Database = Full_Database[,(SSI+1):(SSI+SI[mi])]
  n_var=ncol(MI_Database)
  mat_res = matrix (rep(0, n_obs*n_var), n_obs, n_var)
  crisp_data = matrix (rep(0, n_obs*n_var), n_obs, n_var)
  adjwgt=matrix(rep(0), nrow=n_obs, ncol=n_var)
  b_jkmi <- b_jkt[mi,][b_jkt[mi,] != 0]
  
  for (varindex in 1:n_var){
    
    #Fuzzifier <- GFUZZ(Full_Database, mi, varindex, get(paste0("PA",mi,varindex)), spec=spec, breakpoints = breakpoints)
    Fuzzifier <- get(paste0("MF",mi,varindex))
    
    for (u in 1:n_obs){
      if (is.na(MI_Database[u,varindex]) == TRUE){ mat_res[u,varindex] <- 0
      } else{
        MF.u <- cbind(Fuzzifier[u,,1], rev(Fuzzifier[u,,2]))
        colnames(MF.u) <- c("L", "U")
        mat_res[u,varindex] <- distance(MF.u, TrapezoidalFuzzyNumber(0,0,0,0), type = distance.type, i=i, j=j, theta = theta, p=p, q=q, breakpoints= breakpoints, thetas = thetas)  
        }
      
      adjwgt[u,varindex]=adjusted.weight.SI(MI_Database,u,varindex,b_jkmi)
      
      crisp_data[u,varindex] <- mat_res[u,varindex] * adjwgt[u,varindex]
      
      bsmi[u,mi] <- adjusted.weight.MI(Full_Database[,1:ncol(Full_Database)],u,mi,bmi,b_jkt,SI)
      
    }
    crisp_data[,varindex] <- crisp_data[,varindex]/range[varindex,2] # step for uniforming the scales of the MFs of the variables
}

  crisp_data_MI[,mi] <- rowSums(crisp_data) * max(range[,2])

  SSI=SSI+SI[mi]

}

crisp_data_OBS_NA <- tcrossprod(crisp_data_MI , bsmi)[,1]

return(crisp_data_OBS_NA)
}

# Global evaluation #
#####################

#' Calculates the weighted mean of the set of individual evaluations
#' @param ind.eval the set of individual evaluations.
#' @param weight a vector of the relative sampling weights of the units, for which \eqn{length(weight) = length(ind.eval)}, set by default to \eqn{rep(1, length(ind.eval))}.
#' @return An integer.
#' @export
#' @examples data <- matrix(c(3,4,2,3,3,2,4,3,3,4,3,4,4,2,5,3,4,4,3,3,3,4,4,3,
#' 3,3,4,3,3,3,3,4,4,3,5,3,4,3,3,3), ncol = 4)
#' data <- as.data.frame(data)
#' MI <- 2
#' SI1 <- 2
#' SI2 <- 2
#' SI <- c(SI1,SI2)
#' b_j <- c(1/2,1/2)
#' b_jk <- matrix(c(0.5,0.5,0.5,0.5),nrow=2) 
#' PA11 <- c(1,2,3,4,5)
#' PA12 <- c(1,2,3,4,5)
#' PA21 <- c(1,2,3,4,5)
#' PA22 <- c(1,2,3,4,5)
#' # ------------------
#' MF111 <- TrapezoidalFuzzyNumber(0,2,2,7)
#' MF112 <- TrapezoidalFuzzyNumber(2,7,7,15)
#' MF113 <- TrapezoidalFuzzyNumber(7,15,15,23)
#' MF114 <- TrapezoidalFuzzyNumber(15,23,23,28)
#' MF115 <- TrapezoidalFuzzyNumber(23,28,28,30)
#' MF11 <- GFUZZ(data, 1, 1, PA11, spec="Identical", breakpoints = 100)
#' # ------------------
#' MF121 <- TrapezoidalFuzzyNumber(0,2,2,7)
#' MF122 <- TrapezoidalFuzzyNumber(2,7,7,15)
#' MF123 <- TrapezoidalFuzzyNumber(7,15,15,23)
#' MF124 <- TrapezoidalFuzzyNumber(15,23,23,28)
#' MF125 <- TrapezoidalFuzzyNumber(23,28,28,30)
#' MF12 <- GFUZZ(data, 1, 2, PA12, spec="Identical", breakpoints = 100)
#' # ------------------
#' MF211 <- TrapezoidalFuzzyNumber(0,2,2,7)
#' MF212 <- TrapezoidalFuzzyNumber(2,7,7,15)
#' MF213 <- TrapezoidalFuzzyNumber(7,15,15,23)
#' MF214 <- TrapezoidalFuzzyNumber(15,23,23,28)
#' MF215 <- TrapezoidalFuzzyNumber(23,28,28,30)
#' MF21 <- GFUZZ(data, 2, 1, PA21, spec="Identical", breakpoints = 100)
#' # ------------------
#' MF221 <- TrapezoidalFuzzyNumber(0,2,2,7)
#' MF222 <- TrapezoidalFuzzyNumber(2,7,7,15)
#' MF223 <- TrapezoidalFuzzyNumber(7,15,15,23)
#' MF224 <- TrapezoidalFuzzyNumber(15,23,23,28)
#' MF225 <- TrapezoidalFuzzyNumber(23,28,28,30)
#' MF22 <- GFUZZ(data, 2, 2, PA22, spec="Identical", breakpoints = 100)
#' # ------------------
#' range <- matrix(c(0,0,0,0,28,28,28,28), ncol=2)
#' ind.eval <- IND.EVAL(data,MI,b_j,SI,b_jk, range = range, distance.type ="DSGD.G")
#' GLOB.mean <- GLOB.EVAL.mean(ind.eval)

GLOB.EVAL.mean <- function(ind.eval,weight=rep(1, length(ind.eval))){

    weighted.mean(ind.eval,w=weight)
  
  }

#' Calculates the number of answers by a specific linguistic of a sub-item
#' @param x the data set to evaluate.
#' @param varindex index of a particular sub-item.
#' @param q index of a particular linguistic term.
#' @param p_ind a vector of the relative sampling weights of the units, for which \eqn{length(p_ind) = nrow(data)}. If the weights are not relative, the following expression should be applied on the vector: \deqn{\frac{p_{ind}}{\sum_{i=1}^{n} p_{ind}}.} If no sampling weights are used, the vector of weights is reduced to a vector of values 1, i.e. \eqn{rep(1, nrow(data))}.
#' @return A positive integer.
# #' @export

n_jkq. <- function(x, varindex, q, p_ind = rep(1, nrow(x))){ #mi, si, q, SI, p_ind = rep(1, nrow(x))){
 
  n_jkq. <- sum(p_ind[which(x[, varindex] == q)])
  
  
  #npar <- p_ind * delta_jki(x[,varindex])
  #n_jkq. <- sum(npar)
  
  #for (i in 1:n_obs){
  # npar <- p_ind[i] * delta_jki(x[i,varindex])
  # n_jki <- n_jki + npar 
  #}
    return(n_jkq.)
} 

#' Calculates the number of answers by a specific sub-item
#' @param x the data set to evaluate.
#' @param varindex index of a particular sub-item.
#' @param PA set of possible linguistic terms.
#' @param p_ind a vector of the relative sampling weights of the units, for which \eqn{length(p_ind) = nrow(data)}. If the weights are not relative, the following expression should be applied on the vector: \deqn{\frac{p_{ind}}{\sum_{i=1}^{n} p_{ind}}.} If no sampling weights are used, the vector of weights is reduced to a vector of values 1, i.e. \eqn{rep(1, nrow(data))}.
#' @return A positive integer.
# #' @export

n_jk.. <- function(x, varindex, PA, p_ind = rep(1, nrow(x))){ #mi, si, SI, p_ind = rep(1, nrow(x))){
  
  mat_tot <- matrix(rep(0), nrow=nrow(x), ncol = max(PA))#get(paste0("PA",mi,si))))
  
  for (q in PA){#get(paste0("PA",mi,si))){
    mat_tot[,q] <-p_ind*match(x[, varindex], q)
  }
  
  n_jk.. <- sum(rowSums(mat_tot, na.rm = TRUE))
  
  return(n_jk..)
} 

#' Calculates the global evaluation of a linguistic questionnaire
#' @param Full_Database the data set to evaluate.
#' @param MI a numerical value representing the total number of main-items dividing the linguistic questionnaire.
#' @param bmi an array referring to the initial weights of the main-items.
#' @param SI an array representing the total numbers of sub-items per main-item.
#' @param b_jkt a matrix of MI rows and max(SI) columns expressing the initial weights of each sub-item of a given main-item. 
#' @param p_ind a vector of the relative sampling weights of the units, for which \eqn{length(p_ind) = nrow(data)}. If the weights are not relative, the following expression should be applied on the vector: \deqn{\frac{p_{ind}}{\sum_{i=1}^{n} p_{ind}}.} If no sampling weights are used, the vector of weights is reduced to a vector of values 1, i.e. \eqn{rep(1, nrow(data))}.
#' @param distance.type type of distance chosen from the family of distances, set by default to the signed distance. The different choices are given by: "Rho1", "Rho2", "Bertoluzza", "Rhop", "Delta.pq", "Mid/Spr", "wabl", "DSGD", "DSGD.G", "GSGD".
#' @param i parameter of the density function of the Beta distribution, fixed by default to i = 1.
#' @param j parameter of the density function of the Beta distribution, fixed by default to j = 1.
#' @param theta a numerical value between 0 and 1, representing a weighting parameter. By default, theta is fixed to 1/3 referring to the Lebesgue space. This measure is used in the calculations of the following distances: d_Bertoluzza, d_mid/spr and d_phi-wabl/ldev/rdev.
#' @param thetas a decimal value between 0 and 1, representing the weight given to the shape of the fuzzy number. By default, thetas is fixed to 1. This parameter is used in the calculations of the d_theta star and the d_GSGD distances.
#' @param p a positive integer such that 1 \eqn{\le} p < infinity, referring to the parameter of the Rho_p and Delta_pq.
#' @param q a decimal value between 0 and 1, referring to the parameter of the metric Delta_pq.
#' @param breakpoints a positive arbitrary integer representing the number of breaks chosen to build the numerical alpha-cuts. It is fixed to 100 by default.
#' @return A data set of individual evaluations, for which the number of observations is exactly the same as the initial data set.
#' @importFrom stats weighted.mean
#' @export
#' @examples data <- matrix(c(3,4,2,3,3,2,4,3,3,4,3,4,4,2,5,3,4,4,3,3,3,4,4,3,
#' 3,3,4,3,3,3,3,4,4,3,5,3,4,3,3,3), ncol = 4)
#' data <- as.data.frame(data)
#' MI <- 2
#' SI1 <- 2
#' SI2 <- 2
#' SI <- c(SI1,SI2)
#' b_j <- c(1/2,1/2)
#' b_jk <- matrix(c(0.5,0.5,0.5,0.5),nrow=2) 
#' PA11 <- c(1,2,3,4,5)
#' PA12 <- c(1,2,3,4,5)
#' PA21 <- c(1,2,3,4,5)
#' PA22 <- c(1,2,3,4,5)
#' # ------------------
#' MF111 <- TrapezoidalFuzzyNumber(0,2,2,7)
#' MF112 <- TrapezoidalFuzzyNumber(2,7,7,15)
#' MF113 <- TrapezoidalFuzzyNumber(7,15,15,23)
#' MF114 <- TrapezoidalFuzzyNumber(15,23,23,28)
#' MF115 <- TrapezoidalFuzzyNumber(23,28,28,30)
#' MF11 <- GFUZZ(data, 1, 1, PA11, spec="Identical", breakpoints = 100)
#' # ------------------
#' MF121 <- TrapezoidalFuzzyNumber(0,2,2,7)
#' MF122 <- TrapezoidalFuzzyNumber(2,7,7,15)
#' MF123 <- TrapezoidalFuzzyNumber(7,15,15,23)
#' MF124 <- TrapezoidalFuzzyNumber(15,23,23,28)
#' MF125 <- TrapezoidalFuzzyNumber(23,28,28,30)
#' MF12 <- GFUZZ(data, 1, 2, PA12, spec="Identical", breakpoints = 100)
#' # ------------------
#' MF211 <- TrapezoidalFuzzyNumber(0,2,2,7)
#' MF212 <- TrapezoidalFuzzyNumber(2,7,7,15)
#' MF213 <- TrapezoidalFuzzyNumber(7,15,15,23)
#' MF214 <- TrapezoidalFuzzyNumber(15,23,23,28)
#' MF215 <- TrapezoidalFuzzyNumber(23,28,28,30)
#' MF21 <- GFUZZ(data, 2, 1, PA21, spec="Identical", breakpoints = 100)
#' # ------------------
#' MF221 <- TrapezoidalFuzzyNumber(0,2,2,7)
#' MF222 <- TrapezoidalFuzzyNumber(2,7,7,15)
#' MF223 <- TrapezoidalFuzzyNumber(7,15,15,23)
#' MF224 <- TrapezoidalFuzzyNumber(15,23,23,28)
#' MF225 <- TrapezoidalFuzzyNumber(23,28,28,30)
#' MF22 <- GFUZZ(data, 2, 2, PA22, spec="Identical", breakpoints = 100)
#' # ------------------
#' range <- matrix(c(0,0,0,0,28,28,28,28), ncol=2)
#' ind.eval <- IND.EVAL(data,MI,b_j,SI,b_jk, range = range, distance.type ="DSGD.G")
#' GLOB <- GLOB.EVAL(data, MI, b_j, SI, b_jk, distance.type ="GSGD")


GLOB.EVAL <- function(Full_Database,MI,bmi,SI,b_jkt, p_ind = rep(1/nrow(Full_Database), nrow(Full_Database)), distance.type, i=1, j=1, theta = 1/3, thetas=1, p=2, q=0.5, breakpoints=100){
  
  SI=t(t(SI))
  SSI=0
  crisp_data_SGD_SI=NULL
  crisp_data_SGD_MI=NULL
  crisp_data_SGD_GLOB=NULL
  crisp_data_SGD_GLOB_NA=NULL
  
  
  for (mi in 1:MI)
  {
    
    MI_Database = Full_Database[,(SSI+1):(SSI+SI[mi])]#Full_Database[,(13+SSI+1):(13+SSI+SI[mi])]
    
    # Weight of each variable in the set #
    b_jk=matrix(rep(1/SI[mi],SI[mi])) # In the case of all the variables in a main item having the same weight, otherwize weights should be entered manually
    
    S1=0
    n_jki=0
    n_obs=nrow(MI_Database)
    n_var=ncol(MI_Database)
    mat_obs = matrix (rep(0, n_obs*max(MI_Database, na.rm = TRUE)), n_obs, max(MI_Database, na.rm = TRUE))
    coeff_glob=NULL
    
    # mat_obs2 = matrix (rep(0, n_obs*max(MI_Database, na.rm = TRUE)), n_obs, max(MI_Database, na.rm = TRUE)) # A voir
    
    for (k in 1:n_var)
    {
      mat_obs = matrix (rep(0, n_obs*max(MI_Database, na.rm = TRUE)), n_obs, max(MI_Database, na.rm = TRUE))
      for(u in 1:n_obs)
      {
        for(lt in 1:max(MI_Database, na.rm = TRUE))
        {
          V = MI_Database[,k]
          if (V[u]==lt && is.na(V[u])==FALSE)
          {
            mat_obs[u,lt]= 1  * p_ind[u] * n_obs * adjusted.weight.MI(Full_Database,u,j=mi,b_j=bmi,b_jk=b_jkt,SI) * adjusted.weight.SI(MI_Database,u,k,b_jk)
          }
          #else if (is.na(V[i])==TRUE)
          # {
          #   mat_obs[i,q]=NA
          #}
        }
      }
      
      n_jki = rowSums(mat_obs)
      mat_n_jkq = colSums(mat_obs)
      mat_glob = mat_n_jkq / length(which(!is.na(MI_Database[,k])))#n_obs
      
      
      S1=0
      for(lt in 1:max(get(paste0("PA", mi, k)), na.rm = TRUE))
      {
        S1 = S1 + ( mat_glob[lt] )*( distance(get(paste0("MF",mi,k,lt)) ,TriangularFuzzyNumber(0,0,0), type = distance.type, i=i, j=j, theta = theta, thetas=thetas, p=p, q=q, breakpoints= breakpoints) )      
      }
      coeff_glob[k]=S1      
      
    }
    
    crisp_data_SGD_SI = (sum(coeff_glob))
    
    crisp_data_SGD_MI[mi] = crisp_data_SGD_SI
    
    SSI=SSI+SI[mi]
    crisp_data_SGD_SI=NULL
    
  }
  
  crisp_data_SGD_GLOB = sum(crisp_data_SGD_MI,na.rm = TRUE)
  
  return(crisp_data_SGD_GLOB)
  
}

