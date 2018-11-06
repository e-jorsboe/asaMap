#include <cmath>
#include "readplink.c"
#include "analysis.h"
#include "asaMap.h"
#include "analysisFunction.h"

void kill_pars(pars *p,size_t l){
  delete [] p->gs;
  delete [] p->ys;
  delete [] p->qvec;
  
  kill(p->covs);
  kill(p->design);
  kill(p->ysCgs);
  
  delete [] p->pheno;
  delete [] p->p_sCg;
  delete [] p->start;
  delete [] p->start0;
  
  free(p->bufstr.s);free(p->tmpstr.s);
  delete p;  
  p=NULL;
}


double sd(double *a,int l){
  double ts = 0;
  for(int i=0;i<l;i++)
    ts += a[i];

  double u = ts/(1.0*l);
  assert(u!=0);
  ts = 0;
  for(int i=0;i<l;i++)
    ts += (a[i]-u)*(a[i]-u);
  return ts/(1.0*(l-1.0));
}


// emil - other func
double sd2(const std::vector<double> &phe){
  double ts = 0;
  int l = phe.size();
  for(int i=0;i<l;i++)
    ts += phe[i];

  double u = ts/(1.0*l);
  assert(u!=0);
  ts = 0;
  for(int i=0;i<l;i++)
    ts += (phe[i]-u)*(phe[i]-u);
  return ts/(1.0*(l-1.0));
}

// emil - other version did not calculate sd(y) properly!! provided long length of array and not actual phenos!
void rstart(double *ary,size_t l, const std::vector<double> &phe){
  for(int i=0;i<l;i++){
   ary[i] = drand48()*2-1;
   //  fprintf(stderr,"ary[%d]:%f\n",i,ary[i]);
  }
  // emil - for sd estimate use sd of y
  ary[l]=sd2(phe);
}

// emil new random number generator, other one seemed to always generate negative values
void rstart2(double *ary,size_t l, const std::vector<double> &phe, double fMin, double fMax){
  for(int i=0;i<l;i++){
    double f = (double)rand() / RAND_MAX;
    ary[i] = fMin + f * (fMax - fMin);
   //  fprintf(stderr,"ary[%d]:%f\n",i,ary[i]);
  }
  
  // emil - for sd estimate use sd of y
  ary[l]=sd2(phe);

}


pars *init_pars(size_t l,size_t ncov,int model,int maxInter,double tole,std::vector<double> &start, const std::vector<double> &phe, int regression, FILE* logFile){
  pars *p=new pars;
  p->len=l;
  // also intercept (which is not specified)
  ncov++;
  p->ncov=ncov;
  p->gs=new char[l];
  p->ys=new double[l];
  p->qvec=new double[l];
  //p->mafs=new double[l];
  // has to have column of 1s
  p->covs=initMatrix(l,ncov);

  // emil - added 3rd column for alternate allele
  //p->design=initMatrix(4*l,ncov+2); 
  p->pheno = new double[4*l];
  p->p_sCg = new double[4*l];
  p->bufstr.s=NULL;p->bufstr.l=p->bufstr.m=0;
  p->tmpstr.s=NULL;p->tmpstr.l=p->tmpstr.m=0;

  p->model = model;
  p->regression = regression;

  if(model==0){
    // for the add model where design matrix has B1, B2, A1, covs
    // start has the same + sd(y) - so one longer than row in design matrix    
    ksprintf(&p->bufstr,"Chromo\tPosition\tnInd\tf1\tf2\tllh(M0)\tllh(M1)\tllh(M2)\tllh(M3)\tllh(M4)\tllh(M5)\tb1(M1)\tb2(M1)\tb1(M2)\tb2(M3)\tb(M4)\n");

    // also has to have intercept
    p->design=initMatrix(4*l,ncov+3);
    
    // starting values for 3 coefs and covar + sd at the end (index: ncov+4)
    p->start = new double[ncov+4];
    p->start0 = new double[ncov+4];
    
  } else{
    // for the add model where design matrix has R1, R2, Rm, delta1, delta2, covs (delta: rec for other allele - A)
    // start has the same + sd(y) - so one longer than row in design matrix     
    ksprintf(&p->bufstr,"Chromo\tPosition\tnInd\tf1\tf2\tllh(R0)\tllh(R1)\tllh(R2)\tllh(R3)\tllh(R4)\tllh(R5)\tllh(R6)\tllh(R7)\tb1(R1)\tb2(R1)\tbm(R1)\tb1(R2)\tb2m(R2)\tb1m(R3)\tb2(R3)\tb1(R4)\tb2(R5)\tb(R6)\n");

    // also has to have intercept
    p->design=initMatrix(4*l,ncov+5);    
    // this one has one starting guess for 5 coefs and covar + sd at the end, index ncov+6
    p->start = new double[ncov+6];
    p->start0 = new double[ncov+6];
    
  }
   
  p->ysCgs = initMatrix(l,ncov+6);//now we are just allocating enough enough
  p->maxIter = maxInter;
  p->tol = tole;

  //plugin a start
  for(int i=0;i<start.size();i++)
    p->start[i] = start[i];
  if(start.size()==0){
    //make better start guess at some point
    //void rstart(double *,size_t);
    // emil - will put in random numbers as starting guess - put new random function
    if(model==0){
      //rstart(p->start,ncov+3,phe);//<-will put sd at p->start[p->covs+dy+1]
      rstart2(p->start,ncov+3,phe,-1,1);//<-will put sd at p->start[p->covs+dy+1]

      //if logisitc regression set last start value to 0
      if(p->regression==1){
	p->start[ncov+3]=0;
      }
           
    } else{
      //rstart(p->start,ncov+5,phe);//<-will put sd at p->start[p->covs+dy+1]
      rstart2(p->start,ncov+5,phe,-1,1);//<-will put sd at p->start[p->covs+dy+1]

      //if logisitc regression set last start value to 0
      if(p->regression==1){
	p->start[ncov+5]=0;
      }

    }
  }
   
  //copy it to the start0 which will be copied to start for each new site
  // emil has to be + 4, because also needs to copy sd(y) - last value of start
  if(model==0){
    memcpy(p->start0,p->start,sizeof(double)*(ncov+4));
  } else{
    memcpy(p->start0,p->start,sizeof(double)*(ncov+6));
  }
  return p;
}


void check_pars(pars* p, const std::vector<double> &phe, Matrix<double> &cov){

  // check if all values of first column in cov file is 1
  int interceptChecker = 0;
  // check if phenotype 0 or 1 - for quantitative
  int isBinaryQuan = 1;
  for(int i=0;i<phe.size();i++){
	//if logistic regression check if phenotypes are 0 or 1
	if(p->regression==1){
	  isBinaryQuan = 0;
	  if(phe[i]!=0 and phe[i]!=1){
	    fprintf(stderr,"Phenotypes are not binary (0 or 1) for logistic model!\n");    
	    exit(1);
	  }
	} else{
	  if(phe[i]!=0 and phe[i]!=1){
	    isBinaryQuan = 0;
	  }
	}       

	if(cov.d[i][0]==1){
	  interceptChecker++;
	}
  }
    
  if(isBinaryQuan){
    fprintf(stderr,"\n");
    fprintf(stderr,"#######################################################################################\n");
    fprintf(stderr,"#######################################################################################\n");
    fprintf(stderr,"WARNING: Phenotype for linear model appears to be binary, run a logistic model instead!\n");
    fprintf(stderr,"#######################################################################################\n");
    fprintf(stderr,"#######################################################################################\n");
    fprintf(stderr,"\n");
  }

  if(interceptChecker==p->len){
    fprintf(stderr,"Column of 1s for intercept should not be included in covariate,\n");
    fprintf(stderr,"as this will be added automatically.\n");    
    exit(1);
  }

  
}



void set_pars(pars*p, char *g,const std::vector<double> &phe, const std::vector<double> &ad, double *freq, std::vector<double> start, Matrix<double> &cov, char *site){

  // because has to count up how many individuals analysed for this site - due to missing genotypes
  p->len=0;
  for(int i=0;i<phe.size();i++){
      if(g[i]!=3){
	// not reading in those with missing genotype!!
	p->gs[p->len] = 2-g[i];//DRAGON
	p->ys[p->len] = phe[i];
	
	p->qvec[p->len] = ad[i];

	// only has number of covariates, in covariates file
	for(int c=1;c<=cov.dy;c++){	     
	  p->covs->d[p->len][0] = 1;
	  p->covs->d[p->len][c] = cov.d[i][c-1];
	}
	p->len++;
    }
  }
  
  p->covs->dx=p->len;
  // add extra covariate of intercept
  p->covs->dy=cov.dy+1;
  p->mafs = freq;

  /* EMIL - IS THIS NEEDED? SCALING COVARIATES
  
  //try and protect from getting a too large (X^T*W*X) matrix, that is to be inverted
  //by dividing the start value for the covariate with too high mean with 1/mean(covariate)
  double sum = 0;
  double max = 0;
  for(int j=0;j<p->covs->dy;j++){
    sum = 0; max = 0;
    for(int i=0;i<p->len;i++){
      sum += p->covs->d[i][j];
      if(fabs(p->covs->d[i][j])>max){
	max = fabs(p->covs->d[i][j]);
      }
    }
    double mean = sum / p->len;
    // 4 is pretty arbitrary however geno columns are at most 2
    if(fabs(mean) > 4){
      for(int i=0;i<p->len;i++){
	// scales downs covariate to avoid numerical over/underflow
	// will bring it to interval between -2 and 2
	p->covs->d[i][j] = p->covs->d[i][j]/(max*0.5);
      }
    }
  }             

  EMIL end */  
  
  for(int i=0;i<p->len;i++){
    for(int j=0;j<4;j++){
      p->pheno[i*4+j] = p->ys[i];
      p->p_sCg[i*4+j] = NAN;
    }

  }

  // if add model - copy all starting values ncov + 4
  if(p->model==0){     
    memcpy(p->start,p->start0,sizeof(double)*(p->covs->dy+4));        
  } else{
    // if rec model - copy all starting values ncov + 6
    memcpy(p->start,p->start0,sizeof(double)*(p->covs->dy+6));    
  }

  ksprintf(&p->bufstr,"%s%d\t%f\t%f\t",site,p->len,p->mafs[0],p->mafs[1]);
}


void wrap(const plink *plnk,const std::vector<double> &phe,const std::vector<double> &ad,Matrix<double> &freq,int model,std::vector<double> start,Matrix<double> &cov,int maxIter,double tol,std::vector<char*> &loci,int nThreads,FILE* outFile, FILE* logFile, int regression){
  //fprintf(stderr,"\t-> plinkdim: x->%lu y->%lu\n",plnk->x,plnk->y);
  //return ;
  char **d = new char*[plnk->y];//transposed of plink->d. Not sure what is best, if we filterout nonmissing anyway.
  for(int i=0;i<plnk->y;i++){
    d[i] = new char[plnk->x];
    for(int j=0;j<plnk->x;j++)
      d[i][j]=plnk->d[j][i];
  }
  //we prep for threading. By using encapsulating all data need for a site in struct called pars
  pars *p=init_pars(plnk->x,cov.dy,model,maxIter,tol,start,phe,regression,logFile);

  // check that pheno and covariates are OK
  check_pars(p, phe, cov);
  
  // print starting values!!
  if(p->model==0){
    fprintf(stderr,"Starting values are:\t"); fprintf(logFile,"Starting values are:\t");
    // starting values equal to covs + 4 (3 for add columns + sd(y) column)
    for(int i=0;i<p->ncov+4;i++){
      fprintf(stderr, "%f ",p->start[i]); fprintf(logFile, "%f ",p->start[i]);
    }
    fprintf(stderr,"\n"); fprintf(logFile,"\n");      
  } else{
    // if rec model - copy all starting values ncov + 6
    fprintf(stderr,"Starting values are:\t"); fprintf(logFile,"Starting values are:\t");
    // starting values equal to covs + 6 (5 for rec columns + sd(y) column)
    for(int i=0;i<p->ncov+6;i++){
      fprintf(stderr, "%f ",p->start[i]); fprintf(logFile, "%f ",p->start[i]);
    }
    fprintf(stderr,"\n"); fprintf(logFile,"\n");      
  }
  
  int analysedSites = 0;
  for(int y=0;y<plnk->y;y++){//loop over sites

    // parsing site with index 100 (0-indexed so 101th site) - but will have parsed 100 sites
    if(y % 100==0){
      fprintf(stderr,"Parsed sites:%d\n",y);
    }
    
    int cats2[4] = {0,0,0,0};
   
    for(int x=0;x<plnk->x;x++)//similar to above but with transposed plink matrix
      cats2[d[y][x]]++;
#if 0
    //print table
    fprintf(stdout,"[%d] %d %d %d %d ",y,cats2[0],cats2[1],cats2[2],cats2[3]);
#endif
  
    //discard sites if missingness > 0.1
    if((cats2[3]/(double)(cats2[0]+cats2[1]+cats2[2]))>0.1){
      fprintf(stderr,"skipping site[%d] due to excess missingness\n",y);
      continue;
    }
    //check that we observe atleast 10 obs of 2 different genotypes
    int n=0;
    // changed i to 3, as we do not want to consider missing geno a different genotype
    for(int i=0;i<3;i++)
      if(cats2[i]>1)
	n++;
    // checks that there are at least 2 different genos (1 might be missing)
    if(n<2){
      fprintf(stderr,"skipping site[%d] due to categories filter\n",y);
      continue;
    }
 
    
    //   if(freq.d[y][0]>0.999||freq.d[y][0]<0.001){
    //  fprintf(stderr,"skipping site[%d] due to maf filter\n",y);
    //  continue;
    // }

    set_pars(p,d[y],phe,ad,freq.d[y],start,cov,loci[y]);

    main_anal((void*)p);
    analysedSites++;
    fprintf(outFile,"%s\t%s\n",p->bufstr.s,p->tmpstr.s);
    p->bufstr.l=p->tmpstr.l=0;
    //break;

  }

  // write how many sites analysed 
  fprintf(stderr,"Number of analysed sites is:\t %i\n",analysedSites);
  fprintf(logFile,"Number of analysed sites is:\t %i\n",analysedSites);
  
  fprintf(stderr,"\t-> done\n");
  kill_pars(p,plnk->x);
 
  for(int i=0;i<plnk->y;i++)
    delete [] d[i];
  delete [] d;

}

void print(pars *p,FILE *fp){
  fprintf(fp,"\n-------------\n");
  fprintf(fp,"\n len=%d:\n",p->len);
  for(int i=0;i<p->len;i++)
    fprintf(fp,"%d ",p->gs[i]);
  fprintf(fp,"\ny:\n");
  for(int i=0;i< p->len;i++)
    fprintf(fp,"%f ",p->ys[i]);
  fprintf(fp,"\ndesign:\n");
  for(int i=0;i< p->len;i++){
    for(int c=0;c< p->ncov;c++) 
      fprintf(fp,"%f\t",p->covs->d[i][c]);
    fprintf(fp,"\n");
  }
  fprintf(fp,"\nQ:\n");
  for(int i=0;i< p->len;i++)
    fprintf(fp,"%f\t",p->qvec[i]);
  fprintf(fp,"\n");

  fprintf(fp,"maf:\t%f\t%f\n",p->mafs[0],p->mafs[1]);
  fprintf(fp,"-------------\n");
  fprintf(stderr,"start:\n");
  for(int i=0;i<=p->design->dy;i++)
    fprintf(stderr,"%f ",p->start[i]);
  fprintf(stderr,"\n");
}

