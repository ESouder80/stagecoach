library(readxl)
Caswell_A <- read_excel("C:/Users/Erin/Desktop/Caswell_A.xlsx")
Caswell_P <- read_excel("C:/Users/Erin/Desktop/Caswell_P.xlsx")
Caswell_C <- read_excel("C:/Users/Erin/Desktop/Caswell_C.xlsx")
Caswell_B <- read_excel("C:/Users/Erin/Desktop/Caswell_B.xlsx")
Caswell_D1 <- read_excel("C:/Users/Erin/Desktop/D1.xlsx")
Caswell_D2 <- read_excel("C:/Users/Erin/Desktop/D2.xlsx")
Caswell_D3 <- read_excel("C:/Users/Erin/Desktop/D3.xlsx")
Caswell_D4 <- read_excel("C:/Users/Erin/Desktop/D4.xlsx")
library(readxl)


library(Matrix)
library(expm)


A <- as.matrix(Caswell_A)
B <- as.matrix(Caswell_B)
C <- as.matrix(Caswell_C)
P <- as.matrix(Caswell_P)
D1 <- as.matrix(Caswell_D1)
D2 <- as.matrix(Caswell_D2)
D3 <- as.matrix(Caswell_D3)
D4 <- as.matrix(Caswell_D4)
##Stage Based Information##



Pop_Growth <- function(A) {
  # Need to use the transition matrix A  
  # This function gives the right eigenvalue of matrix A
  results <- eigen(A)
  return(Re(results$values[1]))
}



Stable_Stage <- function(A){
  # Need to use transition matrix A
  # scaled stable stage vectors of input data
  allvectors <- eigen(A) 
  num <- allvectors$vectors[,1]
  # Must specify column number
  denom <- sum(num)
  # scale vectors to sum to 1. 
  results <- num/denom
  return(Re(results))
}


L_eigenval <- function(A){
  # Transpose of matrix A gives you the left eigenvalue
  results <- eigen(t(A))
  return(Re(results$values[1]))
}


Rep_Value <- function(A){
  # Transpose of matrix A gives left eigenvectors also known as reproductive value
  num <- eigen(t(A))$vectors[,1]
  # Must specify colum number
  denom <- sum(num)
  # scale vectors to sum to 1
  results <- num/denom
  return(Re(results))
}



Sens_Mat <- function(A){
  # Using both the eigenvectors found above, sensitivity matrix can be calulated
  # transpose of stable stage distribution * reproduction vector
  num <- Rep_Value(A) %*% t(Stable_Stage(A))
  # stable stage distribution * reproduction vector
  den <- as.numeric(Rep_Value(A) %*% Stable_Stage(A))
  # divide top and bottom to get sensitivity matrix
  results <- num/den
  return(Re(results))
}



Elas_Mat <- function(A){
  # Dominant eigen value * sensitivity matrix
  num <- Pop_Growth(A)
  den <- Sens_Mat(A)
  # lambda * sensitivity matrix divided by matrix A
  results <- (num/den) * (A)
  return(Re(results))
}

Eig_C <- function(C){
  # This function returns the eigenvalue for matrix C
  results <- eigen(C)
  return(Re(results$values[1]))
}


###END OF STAGE BASED INFO###

##AGE BASED INFO##


N_bj <- function(A, B){
  # Equation 19
  # Stage frequency of newborns at stable stage
  # using the stable stage vector that has been scaled to sum 1
  # first multitply matrix B by the stable stage
  num <- B %*% Stable_Stage(A) 
  # sum matrix B by the stable stage vectors
  den <- sum(B %*% Stable_Stage(A))
  # divide to get the values of bj. 
  results <-  (num/den) 
  return(results)
}



age_in_stage <- function(A, B, C){
  # Equation 23
  # Using matrix C which is P + F
  # Identity matrix
  # Need Pop_Growth, N_bj, identity matrix for this function
  I <- diag(dim(A)[1])
  lambda <- domeigval(A)
  bj <- N_bj(A,B)
  num <- rowSums(solve((I - (C/lambda)) %^% 2) %*% bj)
  # solve and multiply twice -> need inverse, squared 
  # sum over the rows and multitply by stage frequency birth rate
  den <- rowSums(solve(I - (C/lambda)) %*% bj)
  # solve -> need inverse
  results <- num/den
  return(Re(results))
}


age_in_stage_sd <- function(A, B, C){
  # Equation 24
  # Standard deviation of mean age in stage
  # need identity matrix
  I <- diag(dim(A)[1])
  # need to sum the rows of matrix C * 1/lambda multiply by 2,Solve for this to get the inverse of the matrices. 
  # multiply by bj
  num <- rowSums(solve((I - 1/Pop_Growth(A) * C) %*% (I - 1/Pop_Growth(A) * C) %*% (I - 1/Pop_Growth(A) * C)) %*% (I + 1/Pop_Growth(A) * C) %*% N_bj(A,B))
  # Identity matrix minus the inverse of lambda * matrix C. Sum rows, Solve for this to get inverse. 
  den <- rowSums(solve(I - 1/Pop_Growth(A) * C) %*% N_bj(A,B))
  # take square root to get standard deviation
  # multiply by bj
  results <- sqrt((num/den) - (age_in_stage(A,B,C)) ^2)
  return(Re(results))
}

scaled <- function(A){
  # Scaling the reproductive value so the first one is 1
  x <- Rep_Value(A)[1]
  results <- Rep_Value(A) / x
  return(sum(results))
}

Gammai <- function(A,B){
  # Equation 12
  # NEED GENERIC
  results <- vector("numeric", length = 6)
  for (i in 1:6) {
    for (j in 1:6) {
      results[i] <- results[i] + (B[j,i] *  scaled(A))
      
    }
    
  }
  
  return(results)
}

age_of_parents <- function(A, B, C){
  # Equation 26
  # sum the total of mean ages, stable stage vectors, and gamma multiplied 
  # and divide by the sum of stable stage vectors and gamma in order to produce 
  num <- sum(age_in_stage(A,B,C) * Stable_Stage(A) * Gammai(A,B)) 
  den <- sum(Stable_Stage(A) * Gammai(A,B))
  results <- num/den
  return(Re(results))
}





age_comp <- function(A,B,C, maxAge = 10){
  #Age composition in stage class i
  #Equation 22
 lam <- Pop_Growth(A)
 bj <- N_bj(A,B)
 I <- diag(dim(C)[1])
 
 for (a in 1: maxAge) {
   num <- lam ^ -a * colSums(C %^% a) * bj
   print(num)
   
 }
 
    
   return(results)
 
}
  
  
stable_age_distribution <- function(A,B,C){
  # Equation 31
  num <- rowSums(com(A,B,C, maxAge = 10) * Stable_Stage(A))
  den <-   rowSums(Stable_Stage(A))
  results <- num/den
  return(results)
  
}


life_expectancy <- function(P){
  # Equation 3
  # Must use identity matrix the same size as your C or P matrix.
  # Must use P matrix *if* there is fission in your C matrix
  I <- diag(dim(P)[1]) 
  # take the inverse of I-c in order to determine life expectancy
  results <- colSums(solve(I - P))
  return(Re(results))
}



life_expectancy_sd <- function(P){
  # Equation 5
  # Must use identity matrix the same size as your C or P matrix
  # Must use P matrix if there is fission in your C matrix
  I <- diag(dim(P)[1])
  #take the inverse and mulitply it through twice and multiply by the sum of (I + C)
  y <- colSums((I + P) %*% (solve(I - P) %^% 2))
  #take the squaretoot of this answer minus the life expectancy squared for the SD
  results <- sqrt(y - (life_expectancy(P))^2) 
  return(Re(results))
}

Di_mat <- function(P){
  # Equation 8
  # NEED GENERIC
  results <- array(0, dim = c(6,6,6))
  for (i  in 1:6) {
    for (l in 1:6) {
      for (k in 1:6) {
        if (i == l) {
          results[k,l,i] <- 0
        } 
        else {
          results[k,l,i] <- P[k,l]
        }
      }
    }
  }
  return(results)
  res <- Di_mat(P)
  print(res)
}


  expected_years <- function(Di){ 
    # Equation 9 
    # Need to use the D transition matrix from above 
    D <- Di
    results <- matrix(0, nrow=dim(D[,,1])[1], ncol=dim(D[,,1])[2])
    # Make new matrix 
    
    for (i in 1:dim(results)[1]) {
      I <- diag(x=1, dim(D[,,i])[1])
      num <- solve((I - D[,,i]) %*% (I - D[,,i])) 
      den <- solve(I - D[,,i]) 
      results[i,] <- (num[i,]/den[i,]) - 1 
    }
    results[is.nan(results)] <- -1
    
    return(results)
  }
  

expected_years_sd <- function(P){
  # Equation 10
  
}

#Life_Span
#Equation 6


survivorshipAll <- function (P, newbornType, TLX=0.00000000001,MAX10 = 10 ) {
  # Equation 2
  res <- NULL
  for (z in 1:(MAX10*2)) {
    res <- cbind(res, colSums(P %^% (z - 1))[newbornType])
    #if (any(res < TLX)) break
  }
  #for loop in order to calculate the total survivorship (l(x))
  res
}


poplx <- function (A, B, P, newbornType) {
  # Table 2
  res <- NULL
  surv <- survivorshipAll(P, newbornType)
  sfnb <- N_bj(A, B)[newbornType]
  for (x in 1:20) {
    res <- cbind(res, surv[,x] * sfnb)
  }
  res <- colSums(res)
  res
}


maternityAll <- function(C, newbornType, TLX=1000,MAX10 = 10 ){
  # Equation 13
  res <- NULL
  for (z in 1:(MAX10*2)) {
    res <- cbind(res, (colSums(C %^% (z - 1)[newbornType] * Gammai(A,B))[newbornType]/colSums(C %^% (z-1)) [newbornType]))
    #if (any(res > TLX)) break
  }
  #for loop in order to calculate the total maternity (f(x))
  res
}

popmx <- function(A,B,C,newbornType){
  # Table 2
  res <- NULL
  mx <- maternityAll(C,newbornType)
  sfnb <- N_bj(A,B)[newbornType]
  for (x in 1:20) {
    res <- cbind(res, mx[,x] * sfnb)
  }
  res <- colSums(res)
  res
}


age_specific_birthAll <- function(A,B,C,newbornType, MAX10 = 10){
  # Equation 32
  # age specific reproductive value (Vx/V1)
  # Using A,B,C, matrix
  # 21/11/2018
  agespecrep <- function(A,B,C,pow,newbornType){
    age <- function(A,B){
      # scaling the reproductive value so sum v * bj = 1
      #this becomes first value in the loop
      num <-  repval(A)
      den <- sum(repval(A) * N_bj(A,B))
      results <- (num/den) 
      return(results)
    }
    v <- age(A,B)
    num <- sapply(newbornType, function (x) colSums(v %*% (C %^% (pow - 1))))
    den <- sapply(newbornType, function (x) colSums(C %^% (pow - 1)))
    results <- num/den
    results[newbornType,1]
  }
  # for loop for time steps
  results <- matrix(age(A,B)[newbornType],nrow=length(newbornType), ncol=1)
  for (z in 2:(MAX10*2)) {
    results <- cbind(results, agespecrep(A,B,C, z, newbornType))
    # if (any(res > TLX)) break
  }
  # for loop in order to calculate the total age specific birth rate (Vx/V1)
  results
} 

popvx <- function(A,B,C,newbornType){
  # Table 2, Equation 33
  res <- NULL
  vx <- agespecrepAll(A,B,C,newbornType)
  sfnb <- N_bj(A,B)[newbornType]
  for (x in 1:20) {
    res <- cbind(res, vx[,x] * sfnb)
  }
  res <- colSums(res)
  res
}



net_rep <- function(A,B,C){
  # Equation 17
  # NEED GENERIC
  I <- diag(dim(C)[1])
  results <- vector("numeric", length=6)
  for (j in 1:6) {
    for (i in 1:6) {
      results[j] <- results[j] + sum(solve(I - C)[i,j] * Gammai(A,B))
    }
  }
  return(results)
}




pop_net_rep <- function(A,B,C,newbornType){
  #Table 2
  stnb <- N_bj(A,B)[newbornType]
  np <- net_rep(A,B,C,newbornType)
  results <- sum(np * stnb)
  return(results)
}

average_age_production <- function(A,B,C,newbornType){
  # Equation 27
  I <- diag(dim(C)[1])
  num <- sapply(newbornType, function(x) colSums(solve((I - C) %^% 2))  %*% colSums(Gammai(A,B)))
  den <- net_rep(A,B,C,newbornType)
  results <- num/den
  results[newbornType]
}

average_age_production_sd <- function(A,B,C,newbornType){
  # Equation 28
  I <- diag(dim(C)[1])
  num <- sapply(newbornType, function(x) colSums((I + C) %*% (solve(I - C) %^% 3)) %*% colSums(Gammai(A,B)))
  den <- net_rep(A,B,C,newbornType)
  results <- sqrt((num/den) - ((mustage(A,B,C,newbornType)) ^ 2))
  results[newbornType]
}

Si <- function(C, newbornType){
 # Equation 29
 # Mean age of residence for each stage 
 #Need to do by newbornType and make NaN = 0
I <- diag(dim(C)[1])
num <- solve((I - C) %^% 2)
den <- solve(I - C)
result <- num/den
return(result)
}

s_d_Si <- function(C){
  # Equation 30
  #Need to fix this
  
  S <- Si(C)
  I <- diag(dim(C)[1])
  num <- (I + C) %*% solve((I - C) %^% 3)
  den <- solve(I - C)
  result <- sqrt((num / den) - S ^ 2)
  return(result)
  
}

pop_Si <- function(A,B,C){
  #Table 2
  bj <- N_bj(A,B)
  I <- diag(dim(C)[1])
  num <- solve((I - C) %^% 2) %*% bj
  den <- solve(I - C) %*% bj
  result <- num/den
  return(result)
  
}

Q <- function(B,P){
 # Equation 14
 # New tranistion matrix called Q. 
 Q <- P
 if ()
  
}

#Age at First Reproduction
# Equation 15

#SD Age at first Reproduction
# Equation 16



