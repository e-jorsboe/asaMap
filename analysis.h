#pragma once

#include <vector>
#include "analysisFunction.h"
#include "readplink.h"
#include "kstring.h"

void wrap(const plink *p,const std::vector<double> &phe,const std::vector<double> &ad,Matrix<double> &freq,int model,std::vector<double> start,Matrix<double> &cov,int,double,std::vector<char*> &loci,int nThreads,FILE* outFile, FILE* logFile, int regression, int estSE, int useM0R0, char* outname);


/*
  Below is the pars.
  This will contain all data needed for analyzing one site.
  The idea is that we allocate pars[nthreads].
  The internal arrays are nsites long. And will therefore not require reallocating.
  
  We assume no missing data (samples are filteret in the set_pars)

*/


typedef struct{
  int len; //length
  Matrix<double> *covs;//covs is a design matrix.
  int ncov; //number of covariates <ncov_max
  char *gs;//genotypes gs \in {0,1,2,3},counts of allele2, 3=NA/missing
  double *ys;//pheno
  double *qvec;//admix prop
  double *start0;//start. Initialized as long as len. Just to be sure.
  double *start;//start. Initialized as long as len. Just to be sure.
  double *mafs;//freqs <- this is never allocated but points to the Matrix<double> freq;
  double *SE;//has standard error of estimates
  int model;// 0=add 1=rec
  int quant;
  int regression;
  int estSE;
  int useM0R0;
  int maxIter;
  double tol;
  int index; //which site it is
  
  //tmp structs to avoid overusing stack
  double *pheno;
  double *p_sCg;
  double *weights;
  Matrix<double> *design;
  Matrix<double> *ysCgs;
  kstring_t bufstr;
  kstring_t tmpstr;
}pars;








