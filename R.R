library(readxl)
Caswell_A <- read_excel("E:/Experiment/Caswell_A.xlsx")
Caswell_P <- read_excel("E:/Experiment/Caswell_P.xlsx")
Caswell_C <- read_excel("E:/Experiment/Caswell_C.xlsx")
Caswell_B <- read_excel("E:/Experiment/Caswell_B.xlsx")
library(readxl)
Caswell_D <- read_excel("E:/Experiment/Caswell_D.xlsx")

library(Matrix)
library(expm)


A <- as.matrix(Caswell_A)
B <- as.matrix(Caswell_B)
C <- as.matrix(Caswell_C)
P <- as.matrix(Caswell_P)
D <- as.matrix(Caswell_D)

##Stage Based Information##



domeigval <- function(A) {
  #Need to use the transition matrix A  
  #funtion accepts data and returns the eigen value of the input data  
  results <- eigen(A)
  return(Re(results$values[1]))
}



stablestagevec <- function(A){
  #Need to use transition matrix A
  #function accepts data and returns scaled stable stage vectors of input data
  allvectors <- eigen(A) 
  num <- allvectors$vectors[,1]
  #Must specify column number
  denom <- sum(num)
  #scale vectors to sum to 1. 
  results <- num/denom
  return(Re(results))
}


leigenval <- function(A){
  #Transpose of matrix A gives you the left eigenvalue
  #function accepts and returns eigen value of transpose input data
  results <- eigen(t(A))
  return(Re(results$values[1]))
}


repval <- function(A){
  #Transpose of matrix A gives left eigenvectors also known as reproductive value
  #function accepts data and returns scaled reproductive value vectors of input data
  num <- eigen(t(A))$vectors[,1]
  #Must specify colum number
  denom <- sum(num)
  #scale vectors to sum to 1
  results <- num/denom
  return(Re(results))
}



sensmat <- function(A){
  #Using both the eigenvectors found above, sensitivity matrix can be calulated
  #function accepts data and returns sensitivity matrix  
  
  #transpose of stable stage distribution * reproduction vector
  num <- repval(A) %*% t(stablestagevec(A))
  #stable stage distribution * reproduction vector
  den <- as.numeric(repval(A) %*% stablestagevec(A))
  #divide top and bottom to get sensitivity matrix
  results <- num/den
  return(Re(results))
}



elasmat <- function(A){
  #function accepts data and returns elasticity matrix
  # dominant eigen value * sensitivity matrix
  num <- domeigval(A)
  den <- sensmat(A)
  # lambda * sensitivity matrix divided by matrix A
  results <- (num/den) * (A)
  return(Re(results))
}

eigC <- function(C){
  #This function returns the eigenvalue for matrix C
  results <- eigen(C)
  return(Re(results$values[1]))
}


###END OF STAGE BASED INFO###

##AGE BASED INFO##


stagefreqnb <- function(A, B){
  #Stage frequency of newborns at stable stage
  #using the stable stage vector that has been scaled to sum 1
  #first multitply matrix B by the stable stage
  num <- B %*% stablestagevec(A) 
  #sum matrix B by the stable stage vectors
  den <- sum(B %*% stablestagevec(A))
  #divide to get the values of bj. 
  results <-  (num/den) 
  return(results)
  }

 

aveage <- function(A, B, C){
  #function accepts data and returns average age at each stage
  #Using matrix C which is P + F
  #Identity matrix
  #Need domeig, stagefreqnb, identity matrix for this function
  I <- diag(dim(C)[1])
  lambda <- domeigval(A)
  bj <- stagefreqnb(A,B)
  num <- rowSums(solve((I - (C/lambda)) %^% 2) %*% bj)
  #solve and multiply twice -> need inverse, squared 
  #sum over the rows and multitply by stage frequency birth rate
  den <- rowSums(solve(I - (C/lambda)) %*% bj)
  #solve -> need inverse
  results <- num/den
  return(Re(results))
}


sdevaveage <- function(A, B, C){
  #function accepts data and returns standard deviation of average ages
  # need identity matrix
  I <- diag(dim(C)[1])
  # need to sum the rows of matrix C * 1/lambda multiply by 2,Solve for this to get the inverse of the matrices. 
  ##multiply by bj
  num <- rowSums(solve((I - 1/domeigval(A) * C) %*% (I - 1/domeigval(A) * C) %*% (I - 1/domeigval(A) * C)) %*% (I + 1/domeigval(A) * C) %*% stagefreqnb(A,B))
  # Identity matrix minus the inverse of lambda * matrix C. Sum rows, Solve for this to get inverse. 
  den <- rowSums(solve(I - 1/domeigval(A) * C) %*% stagefreqnb(A,B))
  #take square root to get standard deviation
  ## multiply by bj
  results <- sqrt((num/den) - (aveage(A,B,C)) ^2)
  return(Re(results))
}

scaled <- function(A){
  #scaling the reproductive value so the first one is 1
  x <- repval(A)[1]
  results <- repval(A) / x
  return(sum(results))
}

Gammai <- function(A,B){
  results <- vector("numeric", length = 6)
  for (i in 1:6) {
    for (j in 1:6) {
      results[i] <- results[i] + (B[j,i] *  scaled(A))
      
    }
    
  }
    
  return(results)
}

ABAR <- function(A, B, C){
  #sum the total of mean ages, stable stage vectors, and gamma multiplied 
  #and divide by the sum of stable stage vectors and gamma in order to produce 
  num <- sum(aveage(A,B,C) * stablestagevec(A) * Gammai(A,B)) 
  den <- sum(stablestagevec(A) * Gammai(A,B))
  results <- num/den
  return(Re(results))
}

#Stable Age Distribution
library(expm)

T <- P

F <- B

w <- stablestagevec(A)
lambda <- domeigval(A)

ageWithinStage <- function (T, F, w, lambda, maxAge = 10, normalise=TRUE) {
  require(expm)
  result <- matrix(NA, nrow=length(w), ncol=maxAge)
  for (i in 0:(maxAge-1)) {
    this.res <- lambda^(-i) * (T %^% i) %*%  F %*% w
    result[, i+1] <- this.res
  }
  if (normalise) {
    I <- diag(length(w))
    denom <- solve(I-1/lambda * T) %*% F %*% w
    result <-  apply(result, 2, function(x) x/denom) 
  }
  result
}


lifexpect <- function(P){
  #Function accepts data and returns vector for L(x)
  #Must use identity matrix the same size as your C or P matrix.
  ##Must use P matrix *if* there is fission in your C matrix
  I <- diag(dim(P)[1]) ## NEED TO REMOVE THIS CONSTANT AND USE THE REAL MATRIX DIMENSIONS
  #take the inverse of I-c in order to determine life expectancy
  results <- colSums(solve(I - P))
  return(Re(results))
}



sdevlife <- function(P){
  #Function accepts data and resturns standard deviation for L(x)
  #Must use identity matrix the same size as your C or P matrix
  ##Must use P matrix if there is fission in your C matrix
  I <- diag(dim(P)[1])
  #take the inverse and mulitply it through twice and multiply by the sum of (I + C)
  y <- colSums((I + P) %*% (solve(I - P) %^% 2))
  #take the squaretoot of this answer minus the life expectancy squared for the SD
  results <- sqrt(y - (lifexpect(P))^2) 
  return(Re(results))
}

Expectedyears <- function(D){
  I <- diag(dim(D)[1])
  num <- solve((I - D) %*% (I - D))
  den <- solve(I - D)
  results <- num/den
  return(results)
}


survivorshipAll <- function (P, newbornType, TLX=0.00000000001,MAX10 = 10 ) {
  #function takes data and provides information for survival(l(x)
  res <- NULL
  for (z in 1:(MAX10*2)) {
    res <- cbind(res, colSums(P %^% (z - 1))[newbornType])
    #if (any(res < TLX)) break
  }
  #for loop in order to calculate the total survivorship (l(x))
  res
}


poplx <- function (A, B, P, newbornType) {
res <- NULL
surv <- survivorshipAll(P, newbornType)
sfnb <- stagefreqnb(A, B)[newbornType]
for (x in 1:20) {
    res <- cbind(res, surv[,x] * sfnb)
  }
res <- colSums(res)
res
}


maternityAll <- function(C, newbornType, TLX=1000,MAX10 = 10 ){
  #Equation 13
    res <- NULL
    for (z in 1:(MAX10*2)) {
      res <- cbind(res, (colSums(C %^% (z - 1)[newbornType] * Gammai(A,B))[newbornType]/colSums(C %^% (z-1)) [newbornType]))
      #if (any(res > TLX)) break
    }
    #for loop in order to calculate the total maternity (f(x))
    res
}

popmx <- function(A,B,C,newbornType){
  #
  res <- NULL
  mx <- maternityAll(C,newbornType)
  sfnb <- stagefreqnb(A,B)[newbornType]
  for (x in 1:20) {
    res <- cbind(res, mx[,x] * sfnb)
  }
  res <- colSums(res)
  res
}


agespecrepAll <- function(A,B,C,newbornType, MAX10 = 10){
  #age specific reproductive value (Vx/V1)
  #Using A,B,C, matrix
  #21/11/2018
agespecrep <- function(A,B,C,pow,newbornType){
  age <- function(A,B){
    #scaling the reproductive value so sum v * bj = 1
    #bj <- stagefreqnb(A,B)
    #this becomes first value in the loop
    num <-  repval(A)
    den <- sum(repval(A) * stagefreqnb(A,B))
    results <- (num/den) 
    return(results)
  }
  v <- age(A,B)
 num <- sapply(newbornType, function (x) colSums(v %*% (C %^% (pow - 1))))
 den <- sapply(newbornType, function (x) colSums(C %^% (pow - 1)))
 results <- num/den
 results[newbornType,1]
}
#for loop for time steps
results <- matrix(age(A,B)[newbornType],nrow=length(newbornType), ncol=1)
for (z in 2:(MAX10*2)) {
  results <- cbind(results, agespecrep(A,B,C, z, newbornType))
  #if (any(res > TLX)) break
}
#for loop in order to calculate the total age specific birth rate (Vx/V1)
results
} 

popvx <- function(A,B,C,newbornType){
  res <- NULL
  vx <- agespecrepAll(A,B,C,newbornType)
  sfnb <- stagefreqnb(A,B)[newbornType]
  for (x in 1:20) {
    res <- cbind(res, vx[,x] * sfnb)
  }
  res <- colSums(res)
  res
}
  


netrep <- function(A,B,C){
  I <- diag(dim(C)[1])
  results <- vector("numeric", length=6)
  for (j in 1:6) {
    for (i in 1:6) {
    results[j] <- results[j] + sum(solve(I - C)[i,j] * Gammai(A,B))
    }
  }
  return(results)
}
  
 


popnetrep <- function(A,B,C,newbornType){
  stnb <- stagefreqnb(A,B)[newbornType]
  np <- netrep(A,B,C,newbornType)
  results <- sum(np * stnb)
 return(results)
}

mustage <- function(A,B,C,newbornType){
  I <- diag(dim(C)[1])
  num <- sapply(newbornType, function(x) colSums(solve((I - C) %^% 2))  %*% colSums(Gammai(A,B)))
  den <- netrep(A,B,C,newbornType)
  results <- num/den
  results[newbornType]
}

stdevmustage <- function(A,B,C,newbornType){
  I <- diag(dim(C)[1])
  num <- sapply(newbornType, function(x) colSums((I + C) %*% (solve(I - C) %^% 3)) %*% colSums(Gammai(A,B)))
  den <- netrep(A,B,C,newbornType)
  results <- sqrt((num/den) - ((mustage(A,B,C,newbornType)) ^ 2))
  results[newbornType]
}





#################START##############################
domeigval(A)

stablestagevec(A)

leigenval(A)

repval(A)

sensmat(A)

elasmat(A)

eigC(C)

stagefreqnb(A,B)

aveage(A,B,C)

sdevaveage(A,B,C)

scaled(A)

Gammai(A,B)

ABAR (A,B,C)



lifexpect(C)

sdevlife(C)

survivorshipAll(P,c(1,3,4,5))

maternityAll(C,c(1,3,4,5))

agespecrepAll(A,B,C,c(1,3,4,5))















