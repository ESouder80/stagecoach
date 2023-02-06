##############################################################
# R implementation of STAGECOACH, Cochran and Ellner (1992)
##############################################################
USER = "Steve"; ## variable to flag where the data files are. 

if(USER == "Steve"){
    home = "c:/repos/stagecoach"; setwd(home); 
    }else{
    home = "C:/Users/Erin/Desktop"
}    
setwd(home); 
source("working R.R"); 

############################################################# 
# Fortran output Caswell.rst lines 1 to 81 
#############################################################

pop_growth(A);   ## matches Dominant Eigenvalue=  .21943D+01

u = stable_stage_dist(A); u/u[1]; # matches fortan output 

u = reproductive_value(A); u/sum(u); # matches fortran output 
scaled_rep_value(A); u/u[1]; # matches 

sensitivy_mat(A); elasticity_mat(A); # matches fortran output 

n_bj(A,B); # gets the right answer; no fortran output.

age_in_stage(A,B,C); # matches lines 76 - 81

age_in_stage_SD(A, B, C); # matches lines 76 - 81 


############################################################# 
# Fortran output Caswell.rst lines 82 to 144 
#############################################################
gam_i(A,B); bet_i(B); ## no analog in the fortran output 

pop_gen_time(A,B,C); ## matches Abar, line 83, with default weighted=TRUE.  

out = pit(A,B,C); apply(out,2,sum);  # should equal vector of all 1's, check. 

u = stable_age_distribution(A,B,C); matrix(u,ncol=5,byrow=TRUE); # matches lines 85 - 89 

life_expectancy(C); life_expectancy_SD(C); # matches lines 91 - 97 

Di_mat(P); # no fortran output, but used in meantime

meantime(P);  # matches lines 99 - 144, "mean time to I" 
              # column j gives mean time to get FROM j, TO other states 
              
meantime_SD(P); # matches lines 99-144, Std. Dev. of time to I 

total_lifeSpan(P);      # What CE92 call 'conditional total lifespan' 
total_lifeSpan_SD(P);   # These two match the fortran output lines 99-144       

############################################################# 
# Fortran output Caswell.rst lines 147 to 244 
#############################################################
       
## The R functions lx, lx_pop, fx, and fx_pop are all correct. 

## When there is only one newborn type, the fortran is right also. 

## When there is more than one, the fortran OK up to age 11 and 
## then column indices get messed up: several are exact duplicates of
## each other, when they should not be. This is now totally obvious
## in caswell.rst. It is unbelievable that nobody noticed it. 

############################################################# 
# Fortran output Caswell.rst lines 246 to 260
#############################################################
net_rep(B,C); # matches fortran net reproductive rates R0(j), with the default weighted=FALSE 

average_age_production(A,B,C,weighted=TRUE,newbornTypes=c(1,3,4,5));  # matches the fortran output for MU1(j) 
age_production_SD(A,B,C,weighted=TRUE,newbornTypes=c(1,3,4,5)); # matches the fortran output for SD(j) 

average_age_production_pop(A,B,C,weighted=TRUE); # matches population output for population MU1, line 259. 
age_production_SD_pop(A,B,C,weighted=TRUE); #  matches STD DEV on line 259 

generation_time(A,B,C,weighted=FALSE); # matches line 260 




