#include <cstdio>
#include <cmath>
#include <cassert>
#include "asaMap.h"
#include "analysis.h"
#include "readplink.h"
#include <cstring>
#include <cfloat>
#include "kstring.h"
#include <iostream>

#ifdef  EIGEN
#include <Eigen/Dense>
#endif

// if under/overflow causes some value to go below this, use this value
double lower_bound = 1e-20;

// fits a linear model using weighted least squares, 4 observations per individual -
// up to 4 possible states, when geno is known
void getFit(double* start, double* Y, double** covMatrix, const double* weights, int nInd4, int nEnv, int df){
   //  fprintf(stderr,"%s: nInd4:%d nEnv:%d df:%d\n",__FUNCTION__,nInd4,nEnv,df);

   /*
     linear regression. Fits a linear model with weights
     Y is the responce (nInd)
     covMatrix is the design matrix [nInd][nEnv]
     weights is the weights (nInd*nEnv)
     int df is the number of degrees of freedom
     nInd is the number of individuals
     nEnv is the number of predictors (including the intercept)
   */
  
   int nIndW=0;
   char nonZeroWeight[nInd4];
   memset(nonZeroWeight,0,nInd4);
   if(weights==NULL){
     nIndW=nInd4;
     memset(nonZeroWeight,1,nInd4);
   } else{
     for(int i=0;i<nInd4;i++){
       // checks that weights are greater than 0, and keeps those individuals
       if(weights[i]>0){
	 //fprintf(stdout,"%f ",weights[i]);
	 nonZeroWeight[i]=1;
	 nIndW++;
       }
     }
   }
   //  fprintf(stderr,"nnIndW:%d\n",nIndW);
   double yw[nIndW]; //<-stripped phenos scaled by stripped weights
   double xw[nIndW*nEnv]; //<-stripped designs

   int cnt=0;
   for(int i=0;i<nInd4;i++){
     if(nonZeroWeight[i]){
       if(weights!=NULL)
	 yw[cnt] = Y[i]*sqrt(weights[i]);
       else
	 yw[cnt] = Y[i];
       //      fprintf(stdout,"\nyw\t%f\n",yw[cnt]);
       for(int j=0;j<nEnv;j++){
	 if(weights!=NULL)
	   xw[cnt*nEnv+j] = covMatrix[i][j] * sqrt(weights[i]);
	 else
	   xw[cnt*nEnv+j] = covMatrix[i][j];
	   
       }
       cnt++;
     }
   }
    
   double XtX[nEnv*nEnv];
   for(int i=0;i<nEnv*nEnv;i++)
     XtX[i]=0;

   // this is doing the matrix product of (X)^T*W*X 
   for(int x=0;x<nEnv;x++)
     for(int y=0;y<nEnv;y++)
       for(int i=0;i<nIndW;i++)
	 XtX[x*nEnv+y]+=xw[i*nEnv+x]*xw[i*nEnv+y];



#if 0
   //print before inversion
   fprintf(stderr,"BEFORE:\n");
   for(int i=0;i<nEnv;i++){
     for(int j=0;j<nEnv;j++)
       fprintf(stderr,"%f ",XtX[i*nEnv+j]);
     fprintf(stderr,"\n");
   } 
#endif

   // doing inv((X)^T*W*X)
   svd_inverse(XtX, nEnv, nEnv);
   //double workspace[2*nEnv];   
   //matinv(XtX, nEnv, nEnv, workspace);

#if 0
   //print after inversion
   fprintf(stderr,"AFTER:\n");
   for(int i=0;i<nEnv;i++){
     for(int j=0;j<nEnv;j++)
       fprintf(stderr,"%f ",XtX[i*nEnv+j]);
     fprintf(stderr,"\n");
   } 
#endif
       
   double Xt_y[nEnv];
   double invXtX_Xt_y[nEnv]; 
   for(int i=0;i<nEnv;i++)
     Xt_y[i]=0;
   for(int i=0;i<nEnv;i++)
     invXtX_Xt_y[i]=0;

   // doing (X)^T*W*Y
   for(int x=0;x<nEnv;x++)
     for(int i=0;i<nIndW;i++)
       Xt_y[x]+=xw[i*nEnv+x]*yw[i];

   // calculating the coefs: inv((X)^T*W*X)*((X)^T*W*Y)
   for(int x=0;x<nEnv;x++)
     for(int y=0;y<nEnv;y++)
       invXtX_Xt_y[x] += XtX[y*nEnv+x]*Xt_y[y];

   for(int x=0;x<nEnv;x++){
     start[x]=invXtX_Xt_y[x];
   }   
   
   double yTilde[nInd4];
   for(int i=0;i<nInd4;i++)
     yTilde[i] = 0;

   for(int i=0;i<nInd4;i++)
     for(int x=0;x<nEnv;x++)
       yTilde[i] += covMatrix[i][x]*start[x];
       
   double ts=0;
   for(int i=0;i<nInd4;i++){
     // getting the residuals
     double tmp = Y[i]-yTilde[i];
     //    fprintf(stdout,"\nYY[%d]\t%f\n",i,tmp);
     if(weights!=NULL){       
       ts += tmp*tmp*weights[i];
     } else{
       ts += tmp*tmp;
     }

   }


   // why did weights change here...?
   // apparently because of residuals variable, where I wrote out of bounds and ended up writing to weights
        
   //fprintf(stderr,"ts:%f\n",ts);
   if(df==-1){
     start[nEnv] = sqrt(ts/(1.0*(nInd4-nEnv)));
   } else{
     start[nEnv] = sqrt(ts/(1.0*df));
   }
        
 }

#ifdef EIGEN

void getFit2(double* start, double* Y, double** covMatrix, const double* weights, int nInd4, int nEnv, int df){
   //  fprintf(stderr,"%s: nInd4:%d nEnv:%d df:%d\n",__FUNCTION__,nInd4,nEnv,df);

   /*
     linear regression. Fits a linear model with weights
     Y is the responce (nInd)
     covMatrix is the design matrix [nInd][nEnv]
     weights is the weights (nInd*nEnv)
     int df is the number of degrees of freedom
     nInd is the number of individuals
     nEnv is the number of predictors (including the intercept)
   */
  
   int nIndW=0;
   char nonZeroWeight[nInd4];
   memset(nonZeroWeight,0,nInd4);
   if(weights==NULL){
     nIndW=nInd4;
     memset(nonZeroWeight,1,nInd4);
   } else{
     for(int i=0;i<nInd4;i++){
       // checks that weights are greater than 0, and keeps those individuals
       if(weights[i]>0){
	 //fprintf(stdout,"%f ",weights[i]);
	 nonZeroWeight[i]=1;
	 nIndW++;
       }
     }
   }
   
   //  fprintf(stderr,"nnIndW:%d\n",nIndW);
   double yw[nIndW]; //<-stripped phenos scaled by stripped weights
   double xw[nIndW*nEnv]; //<-stripped designs

   int cnt=0;
   for(int i=0;i<nInd4;i++){
     if(nonZeroWeight[i]){
       if(weights!=NULL)
	 yw[cnt] = Y[i]*sqrt(weights[i]);
       else
	 yw[cnt] = Y[i];
       //      fprintf(stdout,"\nyw\t%f\n",yw[cnt]);
       for(int j=0;j<nEnv;j++){
	 if(weights!=NULL)
	   xw[j*nIndW+cnt] = covMatrix[i][j] * sqrt(weights[i]);
	 else
	   xw[j*nIndW+cnt] = covMatrix[i][j];	   
       }
       cnt++;
     }
   }   

   Eigen::Map<Eigen::MatrixXd> WY2(yw,nIndW,1);          
   Eigen::Map<Eigen::MatrixXd> WX2(xw,nIndW,nEnv);

   Eigen::VectorXd x1 = WX2.householderQr().solve(WY2);
   
   for(int x=0;x<nEnv;x++){
     start[x]=x1[x];
   }   
    
   double yTilde[nInd4];
   for(int i=0;i<nInd4;i++)
     yTilde[i] = 0;

   for(int i=0;i<nInd4;i++)
     for(int x=0;x<nEnv;x++)
       yTilde[i] += covMatrix[i][x]*start[x];
       
   double ts=0;
   for(int i=0;i<nInd4;i++){
     // getting the residuals
     double tmp = Y[i]-yTilde[i];
     //    fprintf(stdout,"\nYY[%d]\t%f\n",i,tmp);
     if(weights!=NULL){       
       ts += tmp*tmp*weights[i];
     } else{
       ts += tmp*tmp;
     }

   }

   // why did weights change here...?
   // apparently because of residuals variable, where I wrote out of bounds and ended up writing to weights
        
   //fprintf(stderr,"ts:%f\n",ts);
   if(df==-1){
     start[nEnv] = sqrt(ts/(1.0*(nInd4-nEnv)));
   } else{
     start[nEnv] = sqrt(ts/(1.0*df));
   }
        
 }

#endif


void getFitBin(double* start, double* Y, double** covMatrix, const double* weights, int nInd4, int nEnv, int df){
   //  fprintf(stderr,"%s: nInd4:%d nEnv:%d df:%d\n",__FUNCTION__,nInd4,nEnv,df);

   /*
     Fits logistic model using iterativt weighted least squares (IWLS)
     linear regression. Fits a logistic model with weights
     Y is the binary responce (nInd)
     covMatrix is the design matrix [nInd][nEnv]
     weights is the weights (nInd*nEnv)
     int df is the number of degrees of freedom
     nInd is the number of individuals
     nEnv is the number of predictors (including the intercept)
   */
  
   int nIndW=0;
   // SHOULD I MAKE USE of nonZeroWeight in logistic regression
   char nonZeroWeight[nInd4];
   memset(nonZeroWeight,0,nInd4);
   if(weights==NULL){
     nIndW=nInd4;     
     memset(nonZeroWeight,1,nInd4);
   } else{
     for(int i=0;i<nInd4;i++){
       if(weights[i]>0){
	 //fprintf(stdout,"%f ",weights[i]);
	 // counts non zero weights!
	 nonZeroWeight[i]=1;
	 nIndW++;
       }
     }
   }
   
   //  fprintf(stderr,"nnIndW:%d\n",nIndW);
   double yw[nIndW]; //<-stripped phenos scaled by stripped weights
   double ww[nIndW]; //<-stripped weights
   double xw[nIndW*nEnv]; //<-stripped designs

   int cnt=0;
   for(int i=0;i<nInd4;i++){
     if(nonZeroWeight[i]){
       yw[cnt] = Y[i];
       if(weights==NULL){
	 ww[cnt] = 1;
       } else{
	 ww[cnt] = weights[i];
       }
       for(int j=0;j<nEnv;j++){
	 xw[cnt*nEnv+j] = covMatrix[i][j];
       }
       cnt++;
     }
   }
   
   double mustart[nIndW];
   // if weights do not exists, have all weights be 1     
   for(int i=0;i<nIndW;i++){
     mustart[i] = (ww[i]*yw[i] + 0.5)/(ww[i] + 1);
   }
    
   double eta[nIndW];
   for(int i=0;i<nIndW;i++){
     //link
     eta[i] = log(mustart[i]/(1-mustart[i]));
   }
   
   double mu[nIndW];
   for(int i=0;i<nIndW;i++){
     //linkinv
     mu[i] = 1 / (1 + exp(-eta[i]));
   }

   double mu0[nIndW];
   double muetaval[nIndW];
   double z[nIndW];
   double w[nIndW];  
   double XtX[nEnv*nEnv] = {0};
   double Xt_y[nEnv] = {0};
   double invXtX_Xt_y[nEnv] = {0}; 
   int iter = 0;
   
   // we run this 20 times...
   for(int t=0;t<20;t++){
     
     for(int i=0;i<nIndW;i++){
       // because mu is same as linkinv(eta)
       muetaval[i] = mu[i]*(1-mu[i]); 
     }

     // add nonZeriWeight here!!
     for(int i=0;i<nIndW;i++){
       // can be calculated together as do not depent on eachother
       z[i] = eta[i]+(yw[i]-mu[i])/muetaval[i];       
       w[i] = sqrt( (ww[i]*(muetaval[i]*muetaval[i])) / (mu[i]*(1-mu[i])) );       
     }
      
     // this is doing the matrix product of (X)^T*W*X
     // takes all columns of second matrix and puts on first column of first matrix - stores at first 0 to nEnv-1 values of XtX
     // takes all columns of second matrix and puts on second column of first matrix - stores at nEnv to 2*nEnv-1 values of XtX
     // thereby the same as taking rows of transposed first matrix (cols of org) and putting it on all columns of second matrix
     for(int x=0;x<nEnv;x++){
       for(int y=0;y<nEnv;y++){
	 for(int i=0;i<nIndW;i++){
	   // t(xw) %*% xw
	   //XtX[x*nEnv+y]+=(covMatrix[i][x]*w[i])*(covMatrix[i][y]*w[i]);	   
	   XtX[x*nEnv+y]+=xw[i*nEnv+x]*w[i]*xw[i*nEnv+y]*w[i];
	 }
       }
     }     

#if 0
     //print before inversion
     fprintf(stderr,"BEFORE:\n");
     for(int i=0;i<nEnv;i++){
       for(int j=0;j<nEnv;j++)
	 fprintf(stderr,"%f ",XtX[i*nEnv+j]);
       fprintf(stderr,"\n");
     } 
#endif

     svd_inverse(XtX, nEnv, nEnv);     
     //matinv(XtX, nEnv, nEnv, workspace);
    
#if 0
     //print after inversion
     fprintf(stderr,"AFTER:\n");
     for(int i=0;i<nEnv;i++){
       for(int j=0;j<nEnv;j++)
	 fprintf(stderr,"%f ",XtX[i*nEnv+j]);
       fprintf(stderr,"\n");
     } 
#endif 
    
     // this is doing the matrix product of (X)^T*W*Y
     // takes first column of first matrix and puts on first and only column of second matrix
     // takes second column of first matrix and puts on first and only column of second matrix
     for(int x=0;x<nEnv;x++){
       for(int i=0;i<nIndW;i++){	 
	 Xt_y[x] += xw[i*nEnv+x]*w[i]*z[i]*w[i];
       }       
     }

     // calculating coefficients so inv((X)^T*W*X)*((X)^T*W*Y)
     // the coefficients are stored in start
     for(int x=0;x<nEnv;x++){
       for(int y=0;y<nEnv;y++){
	 invXtX_Xt_y[x] += XtX[y*nEnv+x]*Xt_y[y];
       }
       start[x]=invXtX_Xt_y[x];
     }
      
     // eta
     for(int i=0;i<nIndW;i++){
       // clear values of eta
       eta[i]=0;
       for(int x=0;x<nEnv;x++)
	 eta[i] += xw[i*nEnv+x]*start[x];
     }
     double diff = 0;
     
     // mu
     for(int i=0;i<nIndW;i++){
       //linkinv
       mu[i] = 1 / (1 + exp(-eta[i]));

       // we cannot do this for first iteration
       // this will see diff between previous (mu0) and current iteration (mu)

       // fabs is absolute value for floats, abs is for ints
       diff += fabs(mu[i]-mu0[i]);

       // mu0 has values of previous iteration
       mu0[i] = 1 / (1 + exp(-eta[i]));    
     }

     if(diff<1e-4 & iter>0){
       break;
     }

     //not sure I have to calculate residuals...
     
     // clear those that have +=
     // much faster than setting values 0 than in for loop
     memset(XtX, 0, sizeof(XtX));
     memset(Xt_y, 0, sizeof(Xt_y));
     memset(invXtX_Xt_y, 0, sizeof(invXtX_Xt_y));
     
   }

}

#ifdef  EIGEN

void getFitBin2(double* start, double* Y, double** covMatrix, const double* weights, int nInd4, int nEnv, int df){
   //  fprintf(stderr,"%s: nInd4:%d nEnv:%d df:%d\n",__FUNCTION__,nInd4,nEnv,df);

   /*
     Fits logistic model using iterativt weighted least squares (IWLS)
     linear regression. Fits a logistic model with weights
     Y is the binary responce (nInd)
     covMatrix is the design matrix [nInd][nEnv]
     weights is the weights (nInd*nEnv)
     int df is the number of degrees of freedom
     nInd is the number of individuals
     nEnv is the number of predictors (including the intercept)
   */
  
   int nIndW=0;
   // SHOULD I MAKE USE of nonZeroWeight in logistic regression
   char nonZeroWeight[nInd4];
   memset(nonZeroWeight,0,nInd4);
   if(weights==NULL){
     nIndW=nInd4;     
     memset(nonZeroWeight,1,nInd4);
   } else{
     for(int i=0;i<nInd4;i++){
       if(weights[i]>0){
	 //fprintf(stdout,"%f ",weights[i]);
	 // counts non zero weights!
	 nonZeroWeight[i]=1;
	 nIndW++;
       }
     }
   }
   
   //  fprintf(stderr,"nIndW:%d\n",nIndW);
   double yw[nIndW]; //<-stripped phenos scaled by stripped weights
   double ww[nIndW]; //<-stripped weights
   double xw[nIndW*nEnv]; //<-stripped designs

   int cnt=0;
   for(int i=0;i<nInd4;i++){
     if(nonZeroWeight[i]){
       yw[cnt] = Y[i];
       if(weights==NULL){
	 ww[cnt] = 1;
       } else{
	 ww[cnt] = weights[i];
       }
       for(int j=0;j<nEnv;j++){
	 xw[cnt*nEnv+j] = covMatrix[i][j];
       }
       cnt++;
     }
   }
   
   double mustart[nIndW];
   // if weights do not exists, have all weights be 1     
   for(int i=0;i<nIndW;i++){
     mustart[i] = (ww[i]*yw[i] + 0.5)/(ww[i] + 1);
   }
    
   double eta[nIndW];
   for(int i=0;i<nIndW;i++){
     //link
     eta[i] = log(mustart[i]/(1-mustart[i]));
   }
   
   double mu[nIndW];
   for(int i=0;i<nIndW;i++){
     //linkinv
     mu[i] = 1 / (1 + exp(-eta[i]));
   }

   double mu0[nIndW];
   double muetaval[nIndW];
   double z[nIndW];
   double w[nIndW];
   double zw2[nIndW];
   double xw2[nIndW*nEnv];

   int iter = 0;
   
   // we run this 20 times...
   // PROBABLY SHOULD PUT SOME BREAKING CRITERIA
   for(int t=0;t<20;t++){
     
     for(int i=0;i<nIndW;i++){
       // because mu is same as linkinv(eta)
       muetaval[i] = mu[i]*(1-mu[i]); 
     }

     // add nonZeriWeight here!!
     for(int i=0;i<nIndW;i++){
       // can be calculated together as do not depent on eachother
       z[i] = eta[i]+(yw[i]-mu[i])/muetaval[i];       
       w[i] = sqrt( (ww[i]*(muetaval[i]*muetaval[i])) / (mu[i]*(1-mu[i])) );       
     }    

     for(int i=0;i<nIndW;i++){
       for(int x=0;x<nEnv;x++){
	 // t(xw) %*% xw
	 //XtX[x*nEnv+y]+=(covMatrix[i][x]*w[i])*(covMatrix[i][y]*w[i]);

	 // put first weights on first row and then put product as first col in xw2  
	 xw2[x*nIndW+i]=xw[i*nEnv+x]*w[i];
	 zw2[i]=z[i]*w[i];
       }
     }
          
     Eigen::Map<Eigen::MatrixXd> WY2(zw2,nIndW,1);          
     Eigen::Map<Eigen::MatrixXd> WX2(xw2,nIndW,nEnv);
     
     Eigen::VectorXd x1 = WX2.householderQr().solve(WY2);
              
     for(int x=0;x<nEnv;x++){
       start[x]=x1[x];
     }
          
     // eta
     for(int i=0;i<nIndW;i++){
       // clear values of eta
       eta[i]=0;
       for(int x=0;x<nEnv;x++)
	 eta[i] += xw[i*nEnv+x]*start[x];
     }
     double diff = 0;
     
     // mu
     for(int i=0;i<nIndW;i++){
       //linkinv
       mu[i] = 1 / (1 + exp(-eta[i]));

       // we cannot do this for first iteration
       // this will see diff between previous (mu0) and current iteration (mu)

       // fabs is absolute value for floats, abs is for ints
       diff += fabs(mu[i]-mu0[i]);

       // mu0 has values of previous iteration
       mu0[i] = 1 / (1 + exp(-eta[i]));    
     }

     if(diff<1e-4 & iter>0){
       break;
     }
             
   }

}

#endif

//pat[add/rec][x1/x2][X][g]
char pat[2][5][4][3] = {

  // B1 allele from first pop (B is assumed effect allele)
  {{{0,1,2},
    {0,1,1},
    {0,0,1},
    {0,0,0}
    },{
      // B2 allele from second pop
      {0,0,0},
      {0,0,1},
      {0,1,1},
      {0,1,2}
    },{
      // A1 other allele from first pop
      {2,1,0},
      {1,0,0},
      {1,1,0},
      {0,0,0}}   
  },

  // rec elements
  {
    // why is Rm not here, because it is not in any tests...?
    //R1 rec for allele from first pop
    {{0,0,1},
     {0,0,0},
     {0,0,0},
     {0,0,0}
    },{
      //R2 rec for allele from second pop
      {0,0,0},
      {0,0,0},
      {0,0,0},
      {0,0,1}
    },{
      //Rm rec for allele from each pop
      {0,0,0},
      {0,0,1},
      {0,0,1},
      {0,0,0}
    },{
      //R3 rec for other allele from first pop
      {1,0,0},
      {0,0,0},
      {0,0,0},
      {0,0,0}
    },{
      //R4 rec for other allele from second pop
      {0,0,0},
      {0,0,0},
      {0,0,0},
      {1,0,0}
    }}
};

double dnorm(double x,double mean,double sd){
    
  double fac = 1.0/(sd*sqrt(2.0*M_PI));
  double val = exp(-(((x-mean)*(x-mean))/(2*sd*sd)));
  
  // if val is 0 because exp(-(x-mean)*(x-mean)) is due to underflow, returns low value
  if(val<lower_bound){
    return(lower_bound);
  } else{
    return fac*val;
  }

}

double logdnorm(double x, double mean, double sd){    
  double fac = log(1.0/(sd*sqrt(2.0*M_PI)));
  //log of exp(-(((x-mean)*(x-mean))/(2*sd*sd))), exp causes overflow
  double val = -(((x-mean)*(x-mean))/(2*sd*sd));
  return (fac+val);
}


// from /home/albrecht/github/angsd/fet.c
// getting log of binomial
// log\binom{n}{k}
double lbinom(int n, int k){
  if (k == 0 || n == k) return 0;
  return lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
}


// from /home/albrecht/github/angsd/misc/contamination2.cpp (called ldbinom in there)
double logdbinom(int k, int n,double p){
  return lbinom(n,k)+k*log(p)+(n-k)*log(1-p);
}


double dbinom(int k, int n,double p){
  return exp(lbinom(n,k)+k*log(p)+(n-k)*log(1-p));
}


double logbernoulli(int k, double p){
  // if p is 0 or 1, cannot do log
  // however this because of over/underlow and p i just very close 0 or 1
  if(p>1-lower_bound){
    p = 1-lower_bound;
  } else if(p<lower_bound){
    p = lower_bound;
  } 
    return( log(pow(p,k)*pow(1-p,1-k)) );
}

double bernoulli(int k, double p){
  // if p is 0 or 1, cannot do log
  // however this because of over/underlow and p i just very close 0 or 1
  if(p>1-lower_bound){
    p = 1-lower_bound;
  } else if(p<lower_bound){
    p = lower_bound;
  } 
    return( pow(p,k)*pow(1-p,1-k) );
}


double logLike(double *start,double* pheno,Matrix<double> *design,double *p_sCg,int regression){
  double ret = 0;
#if 0
  for(int i=0;i<p->start_len;i++)
    fprintf(stderr,"%f ",p->start[i]);
  fprintf(stderr,"\n");
#endif
      
  double tmp=0;
  for(int i=0;i<design->dx;i++){
    double m=0;
    for(int j=0;j<design->dy;j++){           
      m +=  design->d[i][j]*start[j];      
    }

    if(regression==0){
      m = dnorm(pheno[i],m,start[design->dy]);
    } else{      
      m = bernoulli(pheno[i],(exp(m)/(exp(m)+1)));
    }
    
    tmp += m*p_sCg[i];
    if((i % 4 )==3){

      // same as taking log of all values and then summing and then flipping sign (plus to minus)
      ret -= log(tmp);    
      tmp = 0;
    }    
    
  }

  return ret;
}

// emil more numerically stable - I think - seems not to be working 28-09-2018
double logLike2(double *start,double* pheno,Matrix<double> *design,double *p_sCg,int regression){
#if 0
  for(int i=0;i<p->start_len;i++)
    fprintf(stderr,"%f ",p->start[i]);
  fprintf(stderr,"\n");
#endif

  double tmp[design->dx];
  double ret = 0;
  for(int i=0;i<design->dx;i++){
    double m = 0;
    for(int j=0;j<design->dy;j++){
      m += design->d[i][j]*start[j];
    }
    
    if(regression==0){
      // density function of normal distribution
      m = logdnorm(pheno[i],m,start[design->dy]);      
    } else{      
      double prob = exp(m)/(exp(m)+1.0);
      // now handling if p is 0 or 1 (makes it very small or large)
      // before I had issues with getting inf or nan values - because were taking log of 0!!
      m = logbernoulli(pheno[i],prob);
    }
    
    tmp[i] = m + log(p_sCg[i]);    
    
  }
    
  //colSumsLog
  double maxval;
  double sum = 0;
  for(int i=0;i<design->dx;i++){
    if(i % 4 == 0){
      maxval = tmp[i];
    } else{
      maxval = std::max(maxval,tmp[i]);     
    }
    // trick to avoid over/underflow - log(exp(log(val1)-log(max)) + ...) + log(max) = (exp(log(val1))/exp(log(max)))*(max) + ...
    // same (exp(log(val1))/exp(log(max)))*(max)
    if(i % 3 == 0){
      sum += log(exp(tmp[i-3]-maxval)+exp(tmp[i-2]-maxval)+exp(tmp[i-1]-maxval)+exp(tmp[i]-maxval)) + maxval;
    }
  }
  
  return(-(sum));
  
}


inline double logLikeP(pars *p){
  return logLike(p->start,p->pheno,p->design,p->p_sCg,p->regression);
}

//double updateEM(pars *p){
double updateEM(double *start,Matrix<double> *design,Matrix<double> *ysCgs,double *pheno,int nInd,double *p_sCg, int regression){
#if 0
  for(int i=0;i<=design->dy;i++)
    fprintf(stderr,"%f ",start[i]);
  fprintf(stderr,"\n");
#endif
    
  double ret;
  for(int i=0;i<design->dx;i++){
    double m = 0;
    for(int j=0;j<design->dy;j++)
      m += design->d[i][j]*start[j];
    // density function of normal distribution
    m = logdnorm(pheno[i],m,start[design->dy]);   
    double tmp = m*p_sCg[i];
    //    fprintf(stderr,"(%lu,%d):%f\n",(size_t)floor(i/4),i %4,tmp);
    ysCgs->d[(size_t)floor(i/4)][i % 4] = tmp;
    
  }
  
  // dx is number of indis, dy is 4
  ysCgs->dx = (size_t) design->dx/4;
  ysCgs->dy = 4;

  // to get rowSums
  double ycGs[ysCgs->dx];
  for(int i=0;i<ysCgs->dx;i++){
    double tmp = 0;
    for(int j=0;j<ysCgs->dy;j++){
      tmp += ysCgs->d[i][j];
    }
    ycGs[i] = tmp;
  }  

  // divide each entry of a row with sum of that row
  int df = ysCgs->dx - design->dy;
  for(int i=0;i<ysCgs->dx;i++){
    double tmp = 0;
    for(int j=0;j<ysCgs->dy;j++){     
      ysCgs->d[i][j] /= ycGs[i];
    }
  }
  //print(ysCgs,stdout);
    
  //need to flatten the weights, which is p(s|y,G,phi,Q,f)
  // so first four values is for first indi, and so on...
  double weigths[ysCgs->dx*ysCgs->dy];
  int a = 0;
  for(int i=0;i<ysCgs->dx;i++)
    for(int j=0;j<ysCgs->dy;j++){
      weigths[a++] = ysCgs->d[i][j];
      // fprintf(stdout,"%f ",weigths[a-1]);
    }

 
   //double resi[nInd];
   getFit(start,pheno,design->d,weigths,nInd*4,design->dy,df);
  
#if 0
  for(int i=0;i<=design->dy;i++)
    fprintf(stderr,"%f ",start[i]);
  fprintf(stderr,"\n");
#endif 
  
  return logLike(start,pheno,design,p_sCg,regression);
}


//double updateEM(pars *p){
double logupdateEM(double *start,Matrix<double> *design,Matrix<double> *ysCgs,double *pheno,int nInd,double *p_sCg, int regression){
#if 0
  for(int i=0;i<=design->dy;i++)
    fprintf(stderr,"%f ",start[i]);
  fprintf(stderr,"\n");
#endif 
  
  double ret;
  for(int i=0;i<design->dx;i++){
    double m = 0;
    for(int j=0;j<design->dy;j++){
      m += design->d[i][j]*start[j];
    }

    if(regression==0){
      // density function of normal distribution
      m = logdnorm(pheno[i],m,start[design->dy]);      
    } else{      
      double prob = exp(m)/(exp(m)+1.0);
      // now handling if p is 0 or 1 (makes it very small or large)
      // before I had issues with getting inf or nan values - because were taking log of 0!!
      m = logbernoulli(pheno[i],prob);
    }
    
    double tmp = m + log(p_sCg[i]);
    //    fprintf(stderr,"(%lu,%d):%f\n",(size_t)floor(i/4),i %4,tmp);
    ysCgs->d[(size_t)floor(i/4)][i % 4] = tmp;
  }

  // dx is number of indis, dy is 4
  ysCgs->dx = (size_t) design->dx/4;
  ysCgs->dy = 4;

  // to get rowSums
  double ycGs[ysCgs->dx];
  for(int i=0;i<ysCgs->dx;i++){
    double tmp = 0;
    double maxval = ysCgs->d[i][0];
    for(int j=0;j<ysCgs->dy;j++){
      // find max - part of trick for doing log(p1+p2+p3+p4)
      maxval = std::max(maxval,ysCgs->d[i][j]);     
    }
    // trick to avoid over/underflow - log(exp(log(val1)-log(max)) + ...) + log(max) = (exp(log(val1))/exp(log(max)))*(max) + ...
    // same (exp(log(val1))/exp(log(max)))*(max)
    ycGs[i] = log(exp(ysCgs->d[i][0]-maxval)+exp(ysCgs->d[i][1]-maxval)+exp(ysCgs->d[i][2]-maxval)+exp(ysCgs->d[i][3]-maxval)) + maxval;  
  }

  // divide each entry of a row with sum of that row
  int df = ysCgs->dx - design->dy;
  for(int i=0;i<ysCgs->dx;i++){
    double tmp = 0;
    for(int j=0;j<ysCgs->dy;j++){          
      ysCgs->d[i][j] -= ycGs[i];
    }    
  }
       
  //need to flatten the weights, which is p(s|y,G,phi,Q,f)
  // so first four values is for first indi, and so on...
  double weigths[ysCgs->dx*ysCgs->dy];
  int a = 0;
  for(int i=0;i<ysCgs->dx;i++)
    for(int j=0;j<ysCgs->dy;j++){
      weigths[a++] = exp(ysCgs->d[i][j]);
      //fprintf(stdout,"%f ",weigths[a-1]);

      //check if issue with weights
      if(exp(ysCgs->d[i][j])!=exp(ysCgs->d[i][j]) or std::isinf(exp(ysCgs->d[i][j]))){
	  fprintf(stderr,"Issue with weights being nan or inf\n");
	  return(-9);
	}
    }   
  
  if(regression==0){
    //double resi[nInd];
#ifdef EIGEN   
    getFit2(start,pheno,design->d,weigths,nInd*4,design->dy,df);
#else
    getFit(start,pheno,design->d,weigths,nInd*4,design->dy,df);
#endif
  } else{
#ifdef EIGEN 
    getFitBin2(start,pheno,design->d,weigths,nInd*4,design->dy,df);
#else
    getFitBin(start,pheno,design->d,weigths,nInd*4,design->dy,df);
#endif    
  }
  
#if 0
  for(int i=0;i<=design->dy;i++)
    fprintf(stderr,"%f ",start[i]);
  fprintf(stderr,"\n");
#endif 
  
  return logLike(start,pheno,design,p_sCg,regression);
}

inline double updateEMP(pars *p){  
  return logupdateEM(p->start,p->design,p->ysCgs,p->pheno,p->len,p->p_sCg,p->regression);
}


#define FOR(i,n) for(i=0; i<n; i++)

// SEEMS like there might be a problem with underflow when values of matrix are very high
// APPARENTLY doing matrix inversion for solving linear equations is not a good idea numerically...
// when using age as covariate, it will work every time if starting value is low (0.4 or below approx.),
// and not if starting value is above 0.4 approx..., what to do about this??

// XtX is array (N*N long so works like matrix) with allele counts wieghted by prob of state
// n and m are number of predictors, space is array where results are put into
int matinv( double x[], int n, int m, double space[]){
  //from rasmus nielsens code
  /* x[n*m]  ... m>=n*/
  register int i,j,k; 
  int *irow=(int*) space;
  double ee=1.0e-20, t,t1,xmax;  
  double det=1.0;
  
  FOR (i,n)  {
    xmax = 0.;
    for (j=i; j<n; j++) {
      // finds highest value of X^T*W*X, assuming it has a value > 0
       if (xmax < fabs(x[j*m+i]))  {
	 xmax = fabs( x[j*m+i] );
	 irow[i] = j;
       }
    }
    
    det *= xmax;
    if (xmax < ee)   {

      //fprintf(stderr,"First entry of XtX is %f\n",x[0]);
      fprintf(stderr,"\nDeterminant becomes zero at %3d!\t\n", i+1);
      return(-1);
    }
    
    if (irow[i] != i) {
      FOR (j,m) {
	t = x[i*m+j];
	x[i*m+j] = x[irow[i] * m + j];
	x[ irow[i] * m + j] = t;
      }
    }
    t = 1./x[i*m+i];
    FOR (j,n) {
       if (j == i) continue;
       t1 = t*x[j*m+i];
       FOR(k,m)  x[j*m+k] -= t1*x[i*m+k];
       x[j*m+i] = -t1;
    }
    FOR(j,m)   x[i*m+j] *= t;
    x[i*m+i] = t;
  } /* i  */

  for (i=n-1; i>=0; i--) {
    if (irow[i] == i) continue;
    FOR(j,n)  {
      t = x[j*m+i];
      x[j*m+i] = x[ j*m + irow[i] ];
      x[ j*m + irow[i] ] = t;
    }
  }
  return (0);
  
}


// prob of state given geno, admix and freq, p(s|g, Q, f)
void p_sCg(pars *p){
  
  for(int i=0;i<p->len;i++){
    double t[2] = {p->gs[i]>0?1.0:0,p->gs[i]>1?1.0:0};
    double *f = p->mafs;
    double q = p->qvec[i];
    // fprintf(stderr,"t=(%f,%f)\tf=(%f,%f) q:%f\n",t[0],t[1],f[0],f[1],q);
    
    p->p_sCg[i*4+0] = q*q*pow(f[0],t[0]+t[1]) * pow(1-f[0],2-t[0]-t[1]);
    p->p_sCg[i*4+1] = q*(1-q)*pow(f[0],t[0]) * pow(1-f[0],1-t[0]) * pow(f[1],t[1]) * pow(1-f[1],1-t[1]);
    p->p_sCg[i*4+2] = q*(1-q)*pow(f[1],t[0]) * pow(1-f[1],1-t[0]) * pow(f[0],t[1]) * pow(1-f[0],1-t[1]);
    p->p_sCg[i*4+3] = (1-q)*(1-q)*pow(f[1],t[0]+t[1])*pow(1-f[1],2-t[0]-t[1]);
    
    double s=0;
    for(int j=0;j<4;j++){
      //      fprintf(stderr,"[%d] %f\n",j,p->p_sCg[i*4+j]);
      s += p->p_sCg[i*4+j];
    }
    
    for(int j=0;j<4;j++){
      // fprintf(stderr,"ij:%d\n",i*4+j);
      p->p_sCg[i*4+j] /= s;
    }
   }
  
}


void mkDesign(pars *p){

  assert(p->design->my>=p->covs->dy+3);
  //add model
  if(p->model==0){
    
   for(int i=0;i<p->len;i++){
     for(int n=0;n<3;n++){
       for(int j=0;j<4;j++){
	 // design matrix has 4 rows per indi - possible states when given geno - cols allele and covs
	 // for pat choose: model (add/rec), allele (B1,B2,A1), row and then column
	 p->design->d[i*4+j][n] = pat[p->model][n][j][p->gs[i]];
       }
     }
     // putting in covariates in design matrix, intercept is not added
     for(int c=0;c<p->covs->dy;c++)
       for(int j=0;j<4;j++){
	 p->design->d[i*4+j][3+c] = p->covs->d[i][c];
       }
   }

   // emil - now design has 3 columns (alternative allele also) + number of covariates
   p->design->dy=p->covs->dy+3;
   
   //rec model
  } else{
    
    for(int i=0;i<p->len;i++){
      for(int n=0;n<5;n++){
	for(int j=0;j<4;j++){	  
	  // design matrix has 4 rows per indi - possible states when given geno - cols allele and covs
	  // for pat choose: model (add/rec), allele (B1,B2,A1), row and then column
	  p->design->d[i*4+j][n] = pat[p->model][n][j][p->gs[i]];
	}
      }
      // putting in covariates in design matrix, intercept is not added
      for(int c=0;c<p->covs->dy;c++)
	for(int j=0;j<4;j++){
	  p->design->d[i*4+j][5+c] = p->covs->d[i][c];
	}
    }            
    // emil - now design has 5 columns (alternative allele also) + number of covariates
    p->design->dy=p->covs->dy+5;
  }   
   p->design->dx=p->len*4;
}

void controlEM(pars *p){
  double pars0[p->design->dy+1];
  memcpy(pars0,p->start,sizeof(double)*(p->design->dy+1));
  double llh0 = logLikeP(p);
  //  fprintf(stderr,"\t\t%s[s] like: %f\n",__FUNCTION__,llh0);
  double llh1;
  for(int i=0;i<p->maxIter;i++){
    llh1 = updateEMP(p);
    //fprintf(stderr,"\t\t%s[%d] like: %f diff: %f\n",__FUNCTION__,i,llh1,llh1-llh0);    
    
    if(1){//remove for tole check and direction check
       if(fabs(llh1-llh0)<p->tol){	 
	 //	 fprintf(stderr,"Converged \n");
	 break;
       } else if(llh0<llh1){
	 //	 fprintf(stderr,"Fit caused increase in likelihood, will roll back to previous step\n");		 
	 memcpy(p->start,pars0,sizeof(double)*(p->design->dy+1));	
	 break;
       } else if(llh0<-4){
	 // if issue with weights that are nan or inf
	 // make it print out nan for results
	 for(int i=0;i<p->design->dy+1;i++){
	   p->start[i]=NAN;
	 }	  
	 break;
       }
       
    }
    llh0=llh1;   
    memcpy(pars0,p->start,sizeof(double)*(p->design->dy+1));   
    
    for(int j=0;0&&j<p->design->dy+1;j++)
      fprintf(stderr,"%f ",pars0[j]);
  }
}

//will remove column at(zero indexed)
void rmCol(Matrix<double> *d,int at){
  assert(at>=0 &&at<d->dy);
  for(int i=0;i<d->dx;i++){
     int cnt=0;
     for(int j=0;j<d->dy;j++){
       if(j!=at)
	 d->d[i][cnt++] = d->d[i][j];
     }
  }
  d->dy--;
}

// emil, gets start array, index to remove and length of array
void rmPos(double *d,int at,int l){
  //  fprintf(stderr,"%s %d\n",__FUNCTION__,l);
  for(int i=0;0&&i<l;i++)
    fprintf(stderr,"i[%d]:%f\n",i,d[i]);  
  assert(at>=0 &&at<l);
  int cnt=0;
  for(int j=0;j<l;j++){    
    if(j!=at){
      d[cnt++] = d[j];
    }
  }
  // emil - set last value to 0 - as all values are shifted one to the left
  d[l-1]=0;
  
  for(int i=0;0&&i<l;i++)
    fprintf(stderr,"j[%d]:%f\n",i,d[i]);
}

// ncol is how many coefs we want to write
void printRes(pars *p,int nCol=-1,int printVar=0){
  if(nCol==-1)
    nCol=p->design->dy;
  ksprintf(&p->bufstr,"%f\t",logLikeP(p));
  for(int i=0;i<nCol;i++)
    ksprintf(&p->tmpstr,"%f\t",p->start[i]);
  if(printVar)
    ksprintf(&p->tmpstr,"%f\t",p->start[p->design->dy]); //variancen?
}
    
void printNan(pars *p,int nCol=-1, int printVar=0){
  if(nCol==-1)
    int nCol=p->design->dy;
  ksprintf(&p->bufstr,"%f\t",NAN);
  for(int i=0;i<nCol;i++)
    ksprintf(&p->tmpstr,"%f\t",NAN);
  if(printVar)
    ksprintf(&p->tmpstr,"%f\t",NAN); //variancen?
}


void asamEM(pars *p){

  int maf0=1;
  int maf1=1;
  // if too rare do not run analysis
  if(p->mafs[0]>0.995||p->mafs[0]<0.005){
    maf0=0;
  }
  if(p->mafs[1]>0.995||p->mafs[1]<0.005){
    maf1=0;    
  }
 
  // p->model is either 0 (add model) or 1 (rec model)
  if(p->model==0){

    //////////// do M0 - model with also A1 ///////////////   
    
    mkDesign(p);
    p_sCg(p);
    if(maf0 && maf1){
      controlEM(p);
      printRes(p,0); 
    } else{
      printNan(p,0);
    }
    
    //////////// do M1 ///////////////    
    //remove column3 and third value from start M0
    
    mkDesign(p);
    p_sCg(p);    
    //memcpy(p->start,p->start0,sizeof(double)*(p->design->dy+1));
     
    // remove third column from design - column counting A1, remember start has sd(y) at the end (one longer)
    rmPos(p->start,2,p->covs->dy+3+1);
    rmCol(p->design,2);
    
    if(maf0 && maf1){
      controlEM(p);
      printRes(p,2); 
    } else{
      printNan(p,2);
    }

    double* saveStart = p->start;
    
    //////////// do M2 ///////////////
    //remove column2 and second value from start M1    
    
    mkDesign(p);
    p_sCg(p);
    //memcpy(p->start,p->start0,sizeof(double)*(p->design->dy+1));
    
    // remove B2 second in design, remember start has sd(y) at the end (one longer)
    rmPos(p->start,0,p->covs->dy+3+1);    
    rmCol(p->design,0);
    // remove A1 third in design (NB!! now second cause first was removed), remember start has sd(y) at the end (one longer)
    //rmPos(p->start,1,p->covs->dy+2+1);
    rmCol(p->design,1);
        
    if(maf0){
      controlEM(p);
      printRes(p,1); 
    }  else{
      printNan(p,1);
    }
         
    //////////// do M3 ///////////////
    //remove column1 and first value from start M2   
    
    mkDesign(p);
    p_sCg(p);
    memcpy(p->start,p->start0,sizeof(double)*(p->design->dy+1));

    // remove A1 third in design (NB!! now second cause second was removed), remember start has sd(y) at the end (one longer)
    rmPos(p->start,2,p->covs->dy+3+1);
    rmCol(p->design,2);

    // copies coefs from M1, for faster convergence
    for(int i=0;i<p->covs->dy+2+1;i++){
      p->start[i]=saveStart[i];
    }
    
    // remove B2 second in design, remember start has sd(y) at the end (one longer)
    rmPos(p->start,1,p->covs->dy+2+1);
    rmCol(p->design,1);

    
    if(maf1){
      controlEM(p);
      printRes(p,1); 
    } else{
      printNan(p,1);
    }
     
    //////////// do M4 ///////////////
    //cbind gs and covs into design     
    
    for(int i=0;i<p->covs->dx;i++){
      p->design->d[i][0] = p->gs[i];
      for(int j=0;j<p->covs->dy;j++)
	// not sure if this is necessary
	p->design->d[i][j+1] = p->covs->d[i][j];
    }
    p->design->dx=p->covs->dx;
    p->design->dy=p->covs->dy+1;
    
    p_sCg(p);
    
    // emil - start gets values in getFit - these are the coefs
    for(int i=0;i<p->design->dy+2;i++)
      p->start[i]=NAN;


    if(p->regression==0){
#ifdef EIGEN
      getFit2(p->start,p->ys,p->design->d,NULL,p->design->dx,p->design->dy,-1);
#else
      getFit(p->start,p->ys,p->design->d,NULL,p->design->dx,p->design->dy,-1);
#endif      
    } else{
#ifdef EIGEN
      getFitBin2(p->start,p->ys,p->design->d,NULL,p->design->dx,p->design->dy,-1);
#else
      getFitBin(p->start,p->ys,p->design->d,NULL,p->design->dx,p->design->dy,-1);
#endif      
    }      
    // generating full design matrix for likelihood calculations
    mkDesign(p);
    
    // has to put genotype as first column - design has 4 rows for each indi - REMEMBER THAT!!
    for(int i=0;i<p->covs->dx;i++){
      if( p->gs[i]>2){
	fprintf(stderr,"geno bigger than 3\n");
      }
      // has to have the same geno value for each of the 4 rows
      for(int j=0;j<4;j++){
	p->design->d[i*4+j][0] = p->gs[i];
      }
    }            
    
    // skips B2 and A1
    // emil - has to skip both B2 and A1 (column 2 and 3 design matrix)
    for(int i=p->design->dy;i>0;i--){
      // should let first column be, and move later than 3rd later, 2 columns up (4th->2nd, 5th->3rd,...)
      if(i>2){
	p->start[i] = p->start[i-2];
      }
    }
    
    p->start[2]=p->start[1]=0;
  
    // B2 and A1 has no effect, as they are not in model, B1 is genotype
    
    //p->start[0]=p->start[2];

    //print to kbuf
    printRes(p,1); 
    
    //////////// do M5 ///////////////  
    //cpy covs to design    
    
    for(int i=0;i<p->covs->dx;i++){
      for(int j=0;j<p->covs->dy;j++)
	p->design->d[i][j] = p->covs->d[i][j];
    }
    p->design->dx=p->covs->dx;
    p->design->dy=p->covs->dy;
    
    p_sCg(p);
    for(int i=0;i<p->design->dy+1;i++)
      p->start[i]=NAN;
        
    // double resi[p->design->dx];
    if(p->regression==0){
#ifdef EIGEN
      getFit2(p->start,p->ys,p->design->d,NULL,p->design->dx,p->design->dy,-1);      
#else
      getFit(p->start,p->ys,p->design->d,NULL,p->design->dx,p->design->dy,-1);      
#endif
      
    } else{
#ifdef EIGEN
      getFitBin2(p->start,p->ys,p->design->d,NULL,p->design->dx,p->design->dy,-1);
#else
      getFitBin(p->start,p->ys,p->design->d,NULL,p->design->dx,p->design->dy,-1);
#endif
    }

    // generating full design matrix for likelihood calculations
    mkDesign(p);
        
    // skips A1,A2 for me has to skip both A1,A2 and B1 (column 1,2 and 3 design matrix)
    for(int i=p->design->dy;i>1;i--){
      p->start[i] = p->start[i-3];
    }
    p->start[2]=p->start[1]=p->start[0]=0;
           
    printRes(p,0);
    
  } else{
    
    
    //////////// do R0 - model with also delta1 and delta2 (rec for alternative allele) ///////////////
    // do we want to write any of the effect sizes?    
    
    mkDesign(p);
    p_sCg(p);
     
    if(maf0 && maf1){
      controlEM(p);
      printRes(p,0); 
    } else{
      printNan(p,0);
    }
        
    //////////// do R1 ///////////////    
        
    mkDesign(p);
    p_sCg(p);
    // copy original start values from start0 - because start values have coefs from R0
    //memcpy(p->start,p->start0,sizeof(double)*(p->design->dy+1));
       
    // remove fourth column from design - column counting delta1
    rmPos(p->start,3,p->covs->dy+5+1);
    rmCol(p->design,3);

    // remove fifth column from design (now fourth because fifth column was removed) - delta2
    rmPos(p->start,3,p->covs->dy+4+1);
    rmCol(p->design,3);
     
    if(maf0 && maf1){
      controlEM(p);
      printRes(p,3); 
    } else{
      printNan(p,3);
    }

    double* saveStartRec = p->start;

    //////////// do R2 ///////////////
    //remove column2 and second value from start M2    
    
    mkDesign(p);
    p_sCg(p);
    //memcpy(p->start,p->start0,sizeof(double)*(p->design->dy+1));

    // generating second column: Rm + R2
    for(int i=0;i<p->design->dx;i++){
      p->design->d[i][1] = p->design->d[i][1] + p->design->d[i][2];
    }

    // remove third column from design - column counting delta1
    //rmPos(p->start,3,p->covs->dy+5+1);
    rmCol(p->design,2);

    // remove fifth column from design (now fourth because fifth column was removed) - delta2
    //rmPos(p->start,3,p->covs->dy+4+1);
    rmCol(p->design,3);
       
    // remove Rm third in design, remember start has sd(y) at the end (one longer)
    rmPos(p->start,2,p->covs->dy+3+1);    
    rmCol(p->design,2);            
    
    if(maf0 && maf1){
      controlEM(p);
      printRes(p,2); 
    } else{
      printNan(p,2);
    }           
     
    //////////// do R3 ///////////////
    //remove column1 and first value from start M3    
    
    mkDesign(p);
    p_sCg(p);
    memcpy(p->start,p->start0,sizeof(double)*(p->design->dy+1));

    // generating first column: Rm + R1
    for(int i=0;i<p->design->dx;i++){      
      p->design->d[i][0] = p->design->d[i][0] + p->design->d[i][2];
    }

    // remove fourth column from design - column counting delta1
    rmPos(p->start,3,p->covs->dy+5+1);
    rmCol(p->design,3);

    // remove fifth column from design (now fourth because fifth column was removed) - delta2
    rmPos(p->start,3,p->covs->dy+4+1);
    rmCol(p->design,3);

    // copies coefs from R1, for faster convergence
    for(int i=0;i<p->covs->dy+2+1;i++){
      p->start[i]=saveStartRec[i];
    }
  
    // remove Rm third in design, remember start has sd(y) at the end (one longer)
    rmPos(p->start,2,p->covs->dy+3+1);    
    rmCol(p->design,2);
    
    if(maf0 && maf1){
      controlEM(p);
      printRes(p,2); 
    } else{
      printNan(p,2);
    }

    
    //////////// do R4 ///////////////
    //cbind gs and covs into design M4:    
        
    mkDesign(p);
    p_sCg(p);
    memcpy(p->start,p->start0,sizeof(double)*(p->design->dy+1));
        
    // remove fourth column from design - column counting delta1
    rmPos(p->start,3,p->covs->dy+5+1);
    rmCol(p->design,3);

    // remove fifth column from design (now fourth because fifth column was removed) - delta2
    rmPos(p->start,3,p->covs->dy+4+1);
    rmCol(p->design,3);

    // copies coefs from R1, for faster convergence
    for(int i=0;i<p->covs->dy+2+1;i++){
      p->start[i]=saveStartRec[i];
    }
        
    // remove R2 second in design, remember start has sd(y) at the end (one longer)
    rmPos(p->start,1,p->covs->dy+3+1);    
    rmCol(p->design,1);        
    // remove Rm now second in design (because second already removed), remember start has sd(y) at the end (one longer)    
    rmPos(p->start,1,p->covs->dy+2+1);    
    rmCol(p->design,1);   
    
    if(maf1){
      controlEM(p);
      printRes(p,1); 
    } else{
      printNan(p,1);
    }
    
    //////////// do R5 ///////////////  
    //cpy covs to design M5

    mkDesign(p);
    p_sCg(p);
    memcpy(p->start,p->start0,sizeof(double)*(p->design->dy+1));

    // remove fourth column from design - column counting delta1
    rmPos(p->start,3,p->covs->dy+5+1);
    rmCol(p->design,3);
    
    // remove fifth column from design (now fourth because fifth column was removed) - delta2
    rmPos(p->start,3,p->covs->dy+4+1);
    rmCol(p->design,3);

    // copies coefs from R1, for faster convergence
    for(int i=0;i<p->covs->dy+2+1;i++){
      p->start[i]=saveStartRec[i];
    }
        
       
    // remove R1 first in design, remember start has sd(y) at the end (one longer)
    rmPos(p->start,0,p->covs->dy+3+1);    
    rmCol(p->design,0);        
    // remove Rm now second in design (because first already removed), remember start has sd(y) at the end (one longer)    
    rmPos(p->start,1,p->covs->dy+2+1);    
    rmCol(p->design,1);   

     if(maf1){
      controlEM(p);
      printRes(p,1); 
    } else{
      printNan(p,1);
    }
            
    //////////// do R6 ///////////////  

    // generating a first column of a design matrix has if recessive site (regardless ancestry)
    for(int i=0;i<p->covs->dx;i++){
      if(p->gs[i]==2){
	p->design->d[i][0] = 1;
      } else{
	p->design->d[i][0] = 0;
      }
      for(int j=0;j<p->covs->dy;j++)
	// putting covariates in design matrix from second column and on...
	p->design->d[i][j+1] = p->covs->d[i][j];
    }
    // changing dimension of design matrix accordingly
    p->design->dx=p->covs->dx;
    p->design->dy=p->covs->dy+1;

    p_sCg(p);
    
    // emil - start gets values in getFit - therefore ok with NAN values
    for(int i=0;i<p->design->dy+2;i++)
      p->start[i]=NAN;
    
    if(p->regression==0){
      // fitting a model with current R6 design matrix, getting the coefs (stored in start) - (Ordinary least squares)
#ifdef EIGEN
      getFit2(p->start,p->ys,p->design->d,NULL,p->design->dx,p->design->dy,-1);      
#else
      getFit(p->start,p->ys,p->design->d,NULL,p->design->dx,p->design->dy,-1);      
#endif            
    } else{
#ifdef EIGEN
      getFitBin2(p->start,p->ys,p->design->d,NULL,p->design->dx,p->design->dy,-1);
#else
      getFitBin(p->start,p->ys,p->design->d,NULL,p->design->dx,p->design->dy,-1);
#endif            
    }
    
    // generating full design matrix for likelihood calculations
    mkDesign(p);
    
    // has to put genotype as first column - design has 4 rows for each indi - REMEMBER THAT!!
    for(int i=0;i<p->covs->dx;i++){
      if(p->gs[i]>2){
	fprintf(stderr,"geno bigger than 3\n");
      }
      // has to have the same geno value for each of the 4 rows
      for(int j=0;j<4;j++){
	if(p->gs[i]==2){
	  p->design->d[i*4+j][0] = 1;
	} else{
	  p->design->d[i*4+j][0] = 0;
	}
      }
    }            
    
    // emil - has to skip both R2, Rm, delta1 and delta2 (column 2, 3, 4 and 5 design matrix) - obtained by giving coef of 0 (start values)
    for(int i=p->design->dy;i>0;i--){
      // should let first column be, and move columns after 4th later, 4 columns up (6th->2nd, 7th->3rd,...)
      if(i>2){
	p->start[i] = p->start[i-4];
      }
    }
    // coefs of R2, Rm, delta1 and delta2 being set to 0
    p->start[4]=p->start[3]=p->start[2]=p->start[1]=0;           
    
    //print to kbuf
    printRes(p,1);
 

    //////////// do R7 ///////////////      
      
    for(int i=0;i<p->covs->dx;i++){
      for(int j=0;j<p->covs->dy;j++)
	// not sure if this is necessary
	p->design->d[i][j] = p->covs->d[i][j];
    }
    
    p->design->dx=p->covs->dx;
    p->design->dy=p->covs->dy;
    
    p_sCg(p);
    
    // emil - start gets values in getFit - these are the coefs
    for(int i=0;i<p->design->dy+4;i++)
      p->start[i]=NAN;

    if(p->regression==0){
      // fitting a model with current R6 design matrix, getting the coefs (stored in start) - (Ordinary least squares)
#ifdef EIGEN
      getFit2(p->start,p->ys,p->design->d,NULL,p->design->dx,p->design->dy,-1);      
#else
      getFit(p->start,p->ys,p->design->d,NULL,p->design->dx,p->design->dy,-1);      
#endif            
    } else{
#ifdef EIGEN
      getFitBin2(p->start,p->ys,p->design->d,NULL,p->design->dx,p->design->dy,-1);
#else
      getFitBin(p->start,p->ys,p->design->d,NULL,p->design->dx,p->design->dy,-1);
#endif            
    }

    // generating full design matrix for likelihood calculations
    mkDesign(p);
    
    // emil - has to skip both R1, R2, Rm, delta1 and delta2 (column 1, 2, 3, 4 and 5 design matrix)
    for(int i=p->design->dy;i>0;i--){
      // move column 5 columns up (6th->1st, 7th->2nd,...)
      if(i>2){
	p->start[i] = p->start[i-5];
      }
    }
    p->start[4]=p->start[3]=p->start[2]=p->start[1]=p->start[0]=0;   
            
    //print to kbuf
    printRes(p,0); 
    
    
  }  
  //  fprintf(stdout,"%s:%s\n",p->bufstr.s,p->tmpstr.s);  
  //p->bufstr.l=p->tmpstr.l=0;
}

void main_anal(void *pp){
  asamEM((pars*) pp);
  
}

#if __WITH_MAIN__
int main(){
  fprintf(stderr,"press ctrl+c to exit\n");
  while(1){
    int m;
    int g;
    
    printf( "Enter model and genotype :");
    scanf("%d %d", &m, &g);
    
    printf( "\nModel: %d genotype: %d \n", m, g);
    for(int j=0;j<4;j++){
      for(int n=0;n<2;n++)
	fprintf(stderr,"%d ",pat[m][n][j][g]);
      fprintf(stderr,"\n");
    }
    
  }
  return 0;
}
#endif
