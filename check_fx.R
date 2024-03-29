home = "c:/repos/stagecoach"; setwd(home); 
source("working R.R"); 

 M = c(.10000e+01,  .00000e+00,  .10000e+01,
         .98200e+00,  .19704e+02,  .22755e+01,
         .36171e-01,  .10323e+03,  .10378e+03,
         .10395e-01,  .12940e+03,  .14834e+03,
         .35450e-02,  .16437e+03,  .16232e+03,
         .12619e-02,  .20226e+03,  .16407e+03,
         .40008e-03,  .22198e+03,  .16344e+03,
         .11562e-03,  .23274e+03,  .16264e+03,
         .31391e-04,  .23908e+03,  .16202e+03,
        .81764e-05,  .24303e+03,  .16156e+03,
        .20711e-05,  .24560e+03,  .12517e+03,
       .51472e-06,  .24730e+03,  .12583e+03,
       .12626e-06,  .24846e+03,  .12628e+03,
       .30691e-07,  .24925e+03,  .12658e+03,
       .74134e-08,  .24980e+03,  .12680e+03,
       .17828e-08,  .25018e+03,  .12694e+03,
       .42741e-09,  .25045e+03,  .12705e+03,
       .10224e-09,  .25064e+03,  .12712e+03,
       .24422e-10,  .25077e+03,  .12717e+03,
       .58271e-11,  .25086e+03,  .12720e+03)
M = matrix(M,ncol=3,byrow=TRUE); M; 
        
fx = fx(A,B,C);
fx = data.frame(fx);  

y = M[,2]; 
fit = lm(y~X1+X2+X3+X4+X5+X6-1,data=fx); 
summary(fit); 

Lx = lx(P);
Lx = data.frame(Lx); 
y = M[,1]; 
fit2 = lm(y~X1+X2+X3+X4+X5+X6-1,data=Lx); 
summary(fit2); 

plot(coef(fit),coef(fit2)); abline(0,1,lty=2); 