##############################################################
# R implementation of STAGECOACH, Cochran and Ellner (1992)
#
# Notation: transition matrix A = B + F + P 
#   B = births (sexual), age 0 at birth
#   F = fission (vegetative), regarded as a form of survival
#      F is the paper's notation; here it is matrix Fission 
#   P = survival without fission 
#   C = F + P, both forms of survival
#############################################################

rm(list=ls(all=TRUE)); 
library(readxl)
library(Matrix)
library(expm)

USER = "Steve"; ## variable to flag where the data files are. 

if(USER == "Steve"){
    home = "c:/repos/stagecoach"; setwd(home); 
    }else{
    home = "C:/Users/Erin/Desktop"
}    
setwd(home); 
## source("working R.R"); 
   
Caswell_A <- read_excel("Caswell_A.xlsx",sheet = 1)
Caswell_P <- read_excel("Caswell_P.xlsx")
Caswell_B <- read_excel("Caswell_B.xlsx")
Fission <- read_excel("Caswell_F.xlsx")

  
A <- as.matrix(Caswell_P + Caswell_B + Fission)
B <- as.matrix(Caswell_B)
C <- as.matrix(Caswell_P + Fission)
P <- as.matrix(Caswell_P)

##################################################################
## Stage Based Information 
##################################################################

pop_growth <- function(A) {
  # Population growth: Lambda
  # Need to use the transition matrix A  
  # This function gives the right eigenvalue of matrix A
  results <- eigen(A)
  Re(results$values[1])
  
}

stable_stage_dist <- function(A){
  # Stable stage distribution: "w"
  # Need to use transition matrix A
  # COMMENT: distribution is scaled so that it sums to 1
  allvectors <- eigen(A) 
  num <- allvectors$vectors[,1]
  # Must specify column number 1
  den <- sum(num)
  # scale vectors to sum to 1. 
  results <- num/den
  Re(results) 
}

## (Minor) modification by SPE to ensure real-valued vector 
reproductive_value <- function(A){
  # Transpose of matrix A gives left eigenvectors also known as reproductive value
  num <- Re(eigen(t(A))$vectors[,1]) 
  # Must specify colum number 1
  den <- sum(num)
  # scale vectors to sum to 1
  results <- num/den
  Re(results)
}

sensitivy_mat <- function(A){
  # Using both the eigenvectors found above, sensitivity matrix can be calulated
  # transpose of stable stage distribution * reproduction vector
  num <- reproductive_value(A) %*% t(stable_stage_dist(A))
  # stable stage distribution * reproduction vector
  den <- as.numeric(reproductive_value(A) %*% stable_stage_dist(A))
  # divide top and bottom to get sensitivity matrix
  results <- num/den
  Re(results)
}

elasticity_mat <- function(A){
  # Dominant right eigenvalue and sensitivity matrix
  x <- 1/pop_growth(A)
  y <- sensitivy_mat(A)
  
  ## e(i,j) = A(i,j)*s(i,j)/lambda 
  results <- (x * y) * A
  Re(results)
}

################################################## 
## "Sanity check" - compare sens and elas 
## functions with Steve's standard code. All OK. 
#################################################

sens_and_elas = function(L) {
    x <- eigen(L) 
    w <- Re(x$vectors[,1]); 
    lambda= Re(x$values[1]); 
    x <- eigen(t(L))
    v <- Re(x$vectors[,1])
    sens <- outer(v,w)/sum(v*w)
    elas <- (L/lambda) * sens
    return(list(sens=sens, elas=elas))
}     

## Comparison of results. No difference! 
## S1 = sensitivy_mat(A); E1 = elasticity_mat(A); 
## SE2 = sens_and_elas(A); S2 = SE2$sens; E2 = SE2$elas; 
## range(S1-S2); range(E1-E2); 
#############################################

###END OF STAGE BASED INFO###

##################################################################
## AGE BASED INFO ##
##################################################################

n_bj <- function(A,B){
  # Equation 19
  # Stage frequency of newborns at stable stage
  # using the stable stage vector that has been scaled to sum 1
  # first multitply matrix B by the stable stage
  num <- B %*% stable_stage_dist(A)
  # sum matrix B by the stable stage vectors
  den <- sum(B %*% stable_stage_dist(A))
  # divide to get the values of bj. 
  results <-  num/den
  results
}

############## Matches results from fortran code 
age_in_stage <- function(A, B, C){
  # Equation 23
  # Using matrix C which is P + F
  # Identity matrix
  # Need Pop_Growth, N_bj, identity matrix for this function
  Imat <- diag(dim(A)[1])
  lambda <- pop_growth(A)
  bj <- n_bj(A,B)
  num <- rowSums(solve((Imat - (C/lambda)) %^% 2) %*% bj)
  # solve and multiply twice -> need inverse, squared 
  # sum over the rows and multitply by stage frequency birth rate
  den <- rowSums(solve(Imat - (C/lambda)) %*% bj)
  # solve -> need inverse
  results <- num/den
  results
}

############## Matches results from CE92 fortran code 
age_in_stage_SD <- function(A, B, C){
  # Equation 24
  # Standard deviation of mean age in stage
  # need identity matrix
  Imat <- diag(dim(A)[1])
  # need to sum the rows of matrix C * 1/lambda multiply by 2,Solve for this to get the inverse of the matrices. 
  # multiply by bj
  num <- rowSums(solve((Imat - 1/pop_growth(A) * C) %*% (Imat - 1/pop_growth(A) * C) %*% (Imat - 1/pop_growth(A) * C)) %*% (Imat + 1/pop_growth(A) * C) %*% n_bj(A,B))
  # Identity matrix minus the inverse of lambda * matrix C. Sum rows, Solve for this to get inverse. 
  den <- rowSums(solve(Imat - 1/pop_growth(A) * C) %*% n_bj(A,B))
  # take square root to get standard deviation
  # multiply by bj
  results <- sqrt((num/den) - (age_in_stage(A,B,C)) ^2)
  results
}

## Modified by SPE: it was returning one number, not a vector! 
scaled_rep_value <- function(A){
  # Scaling the reproductive value so the first one is 1
  x <- reproductive_value(A)[1]
  results <- reproductive_value(A)/x
  results
}

### Equation 12, stage-specific total fecundity in newborn-equivalents 
gam_i <- function(A,B){
  # Equation 12
  v <- scaled_rep_value(A); 
  results <- t(B)%*%v; 
  results
}

### Equation 11, stage-specific total fecundity in raw offspring numbers 
bet_i <- function(B){
  # Equation 11
  results <- colSums(B); 
  results
}

## Matches results of fortran code for Abar 
pop_gen_time <- function(A, B, C, weighted = TRUE){
  # Equation 26 when weighted == TRUE. Otherwise treats all kids as equals.
  # The fortran uses weighted == TRUE, hence the default here.  
  # Sum the total of mean ages, stable stage vectors, and gamma multiplied 
  # and divide by the sum of stable stage vectors and gamma to produce Abar 
  
  if(weighted) wts = gam_i(A,B); 
  if(!weighted) wts = bet_i(B); 
  
  num <- sum(age_in_stage(A,B,C) * stable_stage_dist(A) * wts) 
  den <- sum(stable_stage_dist(A) * wts)
  results <- num/den
  results
}

########################################################
# Eqn 22, fraction of age t individuals in stage i. 
# Note, this returns a matrix where p_{i,t}
# is the (t,i) entry of the matrix. Original version and
# this one return the same values. 
########################################################
pit <- function(A,B,C, maxAge = 50){
  # Equation 22
  p <- pop_growth(A)
  #lambda
  b <- n_bj(A,B)
  
  Imat <- diag(x = 1, dim(C))
  # Identity matrix
  
  den <- rowSums(solve(Imat - (C / p)) %*% b) 
  # same for all a, so compute it just once 
  
  results <- matrix(NA, nrow = maxAge + 1, ncol = nrow(A))
  Ca = Imat; 
  # This will be C%^%a on each pass through the loop 
  
  for (a in 0 : (maxAge)) {
    # for loop to determine the fraction of newborns in each stage
    num <- (p ^ -a)*(Ca %*% b)
    x <- num/den
    # print(x)
    results[a + 1,] <- x
    Ca <- C %*% Ca 
  }
  results
}
 # range(pit(A,B,C)-pit_original(A,B,C)) # matches 

stable_age_distribution <- function(A,B,C, maxAge = 50){
  # Equation 31
  pp <- pit(A,B,C,maxAge=maxAge)
  results <- apply(pp, 1, weighted.mean, w = stable_stage_dist(A))
  results
}

#### Same results using matrix multiplication 
sad = function(A,B,C, maxAge=50) {
      pp <- pit(A,B,C,maxAge=maxAge)
      w <- stable_stage_dist(A); 
      num = pp%*%w; den=sum(w); 
      return(num/den)
}      
# stable_age_distribution(A,B,C) - sad(A,B,C); # matches 


# Expected remaining lifespan, conditional on current state 
# Matches the fortran output for this case 
life_expectancy <- function(C){
  # Equation 3
  # Must use identity matrix the same size as your C or P matrix.
  # Must use P matrix *if* there is fission in your C matrix
  Imat <- diag(dim(C)[1]) 
  # take the inverse of I-c in order to determine life expectancy
  results <- colSums(solve(Imat - C))
  results
}

# SD of remaining lifespan, conditional on current state 
# Matches the fortran output for this case 
life_expectancy_SD <- function(C){
  # Equation 5
  # Must use identity matrix the same size as your C or P matrix
  # Must use P matrix if there is fission in your C matrix
  Imat <- diag(dim(C)[1])
  #take the inverse and mulitply it through twice and multiply by the sum of (I + C)
  y <- colSums((Imat + C) %*% (solve(Imat - C) %^% 2))
  # take the square-root of (this answer minus life expectancy squared) for the SD
  results <- sqrt(y - (life_expectancy(C))^2) 
  results
}

######### Calculate the D_i matrices, equation 8  
Di_mat <- function(P){
  results <- array(0, dim = c(dim(P),ncol(P)))
  for (i  in 1:ncol(P)) {
    results[,,i] = P; 
    results[,i,i] = 0; 
  }
  return(results) 
}   
# range(Di_mat(P)-Di_mat_original(P));  ## works 

#####################################################################
# Equation 9, subtracting 1 to match the fortran code!!
# Mean time to get from one state to another. 
# So now this is really "how much time does it take to get there?"
# rather than "when do you get there, if you start at time 1?".   
######################################################################
meantime <- function(P){
   D <- Di_mat(P)
  results <- matrix(0, nrow=dim(D[,,1])[1], ncol=dim(D[,,1])[2])
  # Make new matrix 
  Imat <- diag(x=1, nrow(P))
  for (i in 1:dim(results)[1]) {
    den <- solve(Imat - D[,,i]) 
    num <- den %*% den 
    results[i,] <- (num[i,]/den[i,])-1; 
  }
  results[is.nan(results)] <- NA
  # NA used because can't get to that stage
  
  return(results)
}

# Equation 10 ######################################
# The output matches fortran stagecoach where SD is 
# defined, and correctly gives NA otherwise. 
####################################################
meantime_SD <- function(P){
  D <- Di_mat(P)
  m0 <- meantime(P) ### this 'mean time' subtracts 1 to match the fortan 
  m1 = m0 + 1 ### To match the definition of 'mean time' in eqn. 10 of the paper 
  results <- matrix(0, nrow=nrow(P), ncol=ncol(P))
  
  for (i in 1:dim(results)[1]) {
    Imat <- diag(x=1, nrow(P))
    num <- (Imat + D[,,i]) %*% solve((Imat - D[,,i]) %*% (Imat - D[,,i]) %*% (Imat -D[,,i])) 
    den <- solve(Imat - D[,,i]) 
    results[i,] <- (num[i,]/den[i,])
    results[i,] <- sqrt(results[i,] - (m1[i,]) ^ 2) 
  }
  results[is.nan(results)] <- -1
  # -1 used because can't get to that stage
  results[m0 < 0] <- NA  ## because if transition is impossible, SD is undefined.  
  return(results) 
 }

total_lifeSpan <- function(P){
 # Equation 6, matches fortran output 
 LE <- life_expectancy(P)
 # must add each stage number to the correct vector
 MT <- meantime(P)
 results <- matrix(0,nrow = dim(P)[1],ncol = dim(P)[2])
 for (i in 1:dim(P)[1]) {
   for (j in 1:dim(P)[2]) {
     results[i,j] <- MT[i,j] + LE[i] + 1
     if(i < j){
       results[i,j] <- -1
     }
     else{
       results[i,j] <- MT[i,j] + LE[i] + 1
   }
     }
 }
  results[results<0] = NA;
  return(results); 
}


total_lifeSpan_SD <- function(P){
  #Equation 7. Matches the output of the fortran code!  
  MSD <- meantime_SD(P)
  LSD <- life_expectancy_SD(P)
  results <- matrix(0,nrow = nrow(P), ncol = ncol(P))
  for (i in 1:nrow(P)) {
      results[i,] <- LSD[i]^2 + MSD[i,]^2  
  }
  results <- sqrt(results); 
  CTL = total_lifeSpan(P); 
  results[is.na(CTL)] = NA; # getting to target state is impossible
  return(results) 
}

#########################################################
# Survivorship curve conditional on stage at birth
# Default is to compute this for all stages, but optional
# argument newbornTypes can limit it to only the stages 
# that are actually possible for newborns. 
#########################################################

### This is right, and the fortran is wrong, for >1 newborn type 
lx <- function (P, newbornTypes=NULL, max=20) {
  # Equation 2
  if(is.null(newbornTypes)) newbornTypes=c(1:ncol(P));  
  res = matrix(1, max, ncol(P)); 
  for (x in 2:max) {
    res[x,] = colSums(P%^%(x-1)) 
  }
  return(res[,newbornTypes])
}

### This is right, and the fortran is wrong 
lx_pop <- function (A,B,P,max=20) {
  # Table 2
  l <- lx(P,max=max)
  b <- n_bj(A,B)
  results = l %*% b; 
  return(results)
}

### This is also right, and the fortran is wrong 
lx_pop_v2 <- function(A,B,P,max=20) {
    n = n_bj(A,B); 
    lx = numeric(max); lx[1]=1;
    for(k in 2:max) {
            n = P%*%n; 
            lx[k]=sum(n)
    }        
    return(lx)
}    

## Equation 13, fx based on raw counts or newborn equivalents (if weighted==TRUE) 
## This is right, and the fortan is wrong, for >1 newborn type 
fx = function(A,B,C,newbornTypes=NULL,max=20,weighted=FALSE) {
    if(is.null(newbornTypes)) newbornTypes=c(1:ncol(P)); 
    res = matrix(NA, max, ncol(P)); 
    if(weighted) wts = gam_i(A,B); 
    if(!weighted) wts = bet_i(B); 
    gamma_i = gam_i(A,B); 
    for(x in 1:max) {
        Cx1 = C%^%(x-1); 
        num = t(Cx1) %*% wts; 
        den = colSums(Cx1); 
        res[x,]=num/den; res[x,den==0]=0; 
    }    
    return(res[,newbornTypes]); 
} 


### This is right, and the fortran is wrong 
fx_pop = function(A,B,C,max=20,weighted=FALSE) {
    f <- fx(A,B,C,newbornTypes=c(1:ncol(P)),max=max,weighted=weighted)
    b <- n_bj(A,B)
    results <- (f %*% b)
    results
}


##### R0(j) for each initial state j (can specify newborn types).  
### This is right and matches the fortran 
net_rep <- function(B, C, weighted=FALSE,newbornTypes=NULL){
  # Equation 17
  
  if(is.null(newbornTypes)) newbornTypes=c(1:ncol(C)); 
  Imat <- diag(dim(C)[1])
  N <- solve(Imat - C)
  if(weighted) wts = gam_i(A,B); 
  if(!weighted) wts = bet_i(B); 
  results <- t(N)%*%wts; 
  results[newbornTypes]; 
}

### Population R0
### This is right and matches the fortran  
net_rep_pop <- function(A,B,C,weighted=FALSE){
  #Table 2
  R <- net_rep(B, C, weighted=weighted,newbornTypes=c(1:ncol(C)))
  b <- n_bj(A,B)
  results <- sum(R *b)
 results
}

### Mu_1 definition of generation time 
average_age_production <- function(A,B,C,weighted=FALSE, newbornTypes=NULL){
  # Equation 27
  if(is.null(newbornTypes)) newbornTypes=c(1:ncol(C)); 
  if(weighted) wts = gam_i(A,B); 
  if(!weighted) wts = bet_i(B); 
 
  Imat <- diag(dim(C)[1])
  N = solve(Imat-C)
  num = t(N%*%N)%*%wts; 
  den <- net_rep(B,C,weighted=weighted)
  results <- num/den
  results[newbornTypes]
}

### Variance in age of offspring production 
average_age_production_SD <- function(A,B,C,weighted=FALSE, newbornTypes=NULL){
  # Equation 28
  if(is.null(newbornTypes)) newbornTypes=c(1:ncol(C)); 
  if(weighted) wts = gam_i(A,B); 
  if(!weighted) wts = bet_i(B); 
  Imat <- diag(dim(C)[1])
  N = solve(Imat-C); 
  U = (Imat + C)%*%N%*%N%*%N; 
  num = t(U)%*%wts;   
  den <- net_rep(B,C,weighted=weighted)
  results2 = num/den - average_age_production(A,B,C,weighted=weighted)^2 
  results = sqrt(results2); 
  results[newbornTypes]
}

average_age_production_pop <- function(A,B,C,weighted=FALSE){
  if(weighted) wts = gam_i(A,B); 
  if(!weighted) wts = bet_i(B); 
  Imat <- diag(dim(C)[1])
  N = solve(Imat-C); 
  num = t(N%*%N)%*%wts; 
  b <- n_bj(A,B);
  results = sum(num*b)/net_rep_pop(A,B,C,weighted=weighted) 
  results
}    

############### OK down to here, SPE December 6. 

Vx_V1 <- function(A,B,C,newbornType, MAX10 = 10){
  #Equation 32
  #for loop 
  results <- NULL
  v <- age(A,B)
  for (x in 1:MAX10) {
    results <- cbind(results,(colSums(v * C %^%(x-1))[newbornType])/ colSums(C %^% (x-1))[newbornType])
    
  }
  results
}

Vx_V1_pop <- function(A,B,C,newbornType, MAX10 = 10){
  # Table 2, Equation 33
  #RESULTS DO NOT MATCH PAPER
  res <- NULL
  vx <- Vx_V1(A,B,C,newbornType,MAX10)
  n <- n_bj(A,B)[newbornType]
  for (x in 1:MAX10) {
    res <- cbind(res, vx[,x] * n)
  }
  results <- colSums(res)
 results
}

mean_age_residence <- function(C){
 # Equation 29
 # Mean age of residence for each stage 
 #Need to do by newbornType and make NaN = 0
Imat <- diag(dim(C)[1])
# Identity matrix
num <- solve((Imat - C) %^% 2)
den <- solve(Imat - C)
results <- (num/den)
results[is.nan(results)] <- 0
results
}

mean_age_residence_SD <- function(C){
  # Equation 30
  #Need to fix this so it just gives the newbornTypes 1,3,4,5
  
  S <- mean_age_residence(C,newbornType)
  Imat <- diag(dim(C)[1])
  # identity matrix
  num <- (Imat + C) %*% solve((Imat - C) %^% 3)
  den <- solve(Imat - C)
  results <- abs(sqrt((num / den) - S ^ 2))
  results[is.nan(results)] <- 0
 results
  
}

mean_age_residence_pop <- function(A,B,C){
  #Table 2
  bj <- n_bj(A,B)
  Imat <- diag(dim(C)[1])
  num <- solve((Imat - C) %^% 2) %*% bj
  den <- solve(Imat - C) %*% bj
  results <- num/den
  results
  
}


age <- function(A,B){
  # scaling the reproductive value so sum v * bj = 1
  #this becomes first value in the loop for Vx_V1
  num <-  repro_value(A)
  den <- sum(repro_value(A) * n_bj(A,B))
  results <- (num/den) 
  return(results)
}



mean_age_residence_pop_SD <- function(A,B,C){
  #Table 2
  #DOESN'T MATCH OUTPUT
  s <- mean_age_residence_pop(A,B,C)
  b <- n_bj(A,B)
  Imat <- diag(x = 1, dim(C))
  # identity matrix
  num <- ((Imat + C) %*% solve(Imat - C) %^% 3) %*% b
  den <- (solve(Imat - C)) %*% b
  results <- sqrt((num/den) - s ^ 2)
  results
}


Q_mat <- function(P,stage){
 # Equation 14
 # New transition matrix called Q.
 # Use the stage where birth occurs. For Caswell, it is only column 6
  results <- matrix(0, nrow = nrow(P), ncol = ncol(P))
  for (i in 1:dim(P)[1]) {
    for (j in  1: dim(P)[2]){
      if(j == stage){
        results[i,j] <- 0
      }
      else{
        results[i,j] <- P[i,j]
      }
    }
    }
  
 results
}  


maturity_age <- function(P,Q_mat,stage,newbornType){
  # Equation 15
  # Use the transition matrix from above
  q <- Q_mat(P,stage)
  Imat <- diag(x = 1, dim(P))
  # identity matrix
  num <- (solve((Imat - q) %*% (Imat - q)))
  den <- (solve(Imat - q))
  
  results <- num/den 
  results[stage,][newbornType]
}


maturity_age_SD <- function(P,Q_mat,stage,newbornType){
# Equation 16
  #DOESN'T MATCH OUTPUT
 q <- Q_mat(P,stage)
 m <- maturity_age(P,Q_mat,stage,newbornType)
 Imat <- diag(dim(P)[1])
  num <- (Imat + q) %*% (solve((Imat - q) %*% (Imat - q) %*% (Imat - q)))
  den <- solve(Imat - q)
   results <- sqrt(abs(num/den) - m^2)
 results[stage,][newbornType]
}


generation_time <- function(A,B,C,newbornType){
  lam <- pop_growth(A)
  R <- sum(net_rep_pop(A,B,C,newbornType))
  results <- log(R) / log(lam)
  results
}

if(FALSE) {
#########################################
#Running program

Pop_Growth(A)

Stable_Stage(A)

Rep_Value(A)

Sens_Mat(A)

Elas_Mat(A)

N_bj(A,B)

Age_in_stage(A,B,C)

Age_in_stage_sd(A,B,C)

scaled(A)

g(B)

Gammai(A,B)

Age_of_parents(A,B,C)

#

#

Life_expectancy(P)

Life_expectancy_sd(P)

Di_mat(P)

Expected_years(Di)

#

#

lx(P,c(1,3,4,5))

poplx(A,B,P,c(1,3,4,5))

fx(C,c(1,3,4,5))

fx_pop(A,B,C,c(1,3,4,5))

Vx_V1(A,B,C,c(1,3,4,5))

Vx_V1_pop(A,B,C,c(1,3,4,5))

net_rep(A,B,C)

net_rep_pop(A,B,C,c(1,3,4,5))

average_age_production(A,B,C,c(1,3,4,5))

average_age_production_sd(A,B,C,c(1,3,4,5))

Si(C)

Si_sd(C)

pop_Si(A,B,C)

Q_mat(P)

age_first_reproduction(Q)

age_first_reproduction_sd(Q)

generation_time(A,B,C,c(1,3,4,5))

}
