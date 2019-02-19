#include <cmath>
#include <algorithm> 
#include "readplink.c"
#include "analysis.h"
#include "asaMap.h"
#include "analysisFunction.h"
#include <pthread.h>

typedef struct{
  int index; //site 
  char *res; //<-results as char*
}myres;


// results in vector
std::vector<myres> myresults;

bool myresSort (myres a,myres b){  
  return a.index<b.index;
}

// how many sites analysed in threaded mode
int analysedThreadSites = 0;
int parsedThreadSites = 0;

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
  delete [] p->SE;
  delete [] p->weights;
  
  free(p->bufstr.s);free(p->tmpstr.s);
  delete p;  
  p=NULL;
  
}


struct worker_args_t{

  int first;
  int second;
  pars* p;  
  char** d;//covs is a design matrix.
  std::vector<double> phe;
  std::vector<double> ad;//admix prop
  std::vector<double> start;//start. Initialized as long as len. Just to be sure.
  Matrix<double> freq;
  Matrix<double> cov;
  std::vector<char*> loci;
  int i;
  FILE* tmpFile;
  char* tmpFileName;
  int nSites;

  // construct for worker_args_t struct
  worker_args_t(int &first1,int &second1,pars* p1,char** d1,std::vector<double> phe1,std::vector<double> ad1,std::vector<double> start1,Matrix<double> &freq1,Matrix<double> &cov1,std::vector<char*> loci1,int &i1){

    first=first1;
    second=second1;
    p=p1;
    d=d1;
    phe=phe1;
    ad=ad1;
    start=start1;
    freq=freq1;
    cov=cov1;
    loci=loci1;
    i=i1;
    nSites=freq.dx;
    
  }
  
};

typedef struct worker_args_t worker_args;

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


// get standard deviation
double sd(const std::vector<double> &phe){
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


// emil new random number generator
void rstart(double *ary,size_t l, const std::vector<double> &phe, double fMin, double fMax){
  for(int i=0;i<l;i++){
    double f = (double)rand() / RAND_MAX;
    ary[i] = fMin + f * (fMax - fMin);
   //  fprintf(stderr,"ary[%d]:%f\n",i,ary[i]);
  } 
  // emil - for sd estimate use sd of y
  ary[l]=sd(phe);
}

pars *init_pars(size_t l,size_t ncov,int model,int maxInter,double tole,std::vector<double> &start, const std::vector<double> &phe, int regression, FILE* logFile, int estSE, int useM0R0, int header){
  
  pars *p=new pars;
  p->len=l;
  // also intercept (which is not specified)
  ncov++;
  p->ncov=ncov;
  p->gs=new char[l];
  p->ys=new double[l];
  p->qvec=new double[l];
  // has to have column of 1s
  p->covs=initMatrix(l,ncov);

  p->pheno = new double[4*l];
  p->p_sCg = new double[4*l];
  p->weights = new double[4*l];
  p->bufstr.s=NULL;p->bufstr.l=p->bufstr.m=0;
  p->tmpstr.s=NULL;p->tmpstr.l=p->tmpstr.m=0;

  p->model = model;
  p->regression = regression;
  p->estSE = estSE;
  p->useM0R0 = useM0R0;
  
  if(model==0){
    // for the add model where design matrix has B1, B2, A1, covs
    // start has the same + sd(y) - so one longer than row in design matrix
    if(header){
      if(estSE and useM0R0){
	ksprintf(&p->bufstr,"Chromo\tPosition\tnInd\tf1\tf2\tllh(M0)\tllh(M1)\tllh(M2)\tllh(M3)\tllh(M4)\tllh(M5)\tb1(M1)\tb2(M1)\tSE(b1(M1))\tSE(b2(M1))\tb1(M2)\tSE(b1(M2))\tb2(M3)\tSE(b2(M3))\tb(M4)\tSE(b(M4))\n");      
      } else if(estSE){
	ksprintf(&p->bufstr,"Chromo\tPosition\tnInd\tf1\tf2\tllh(M1)\tllh(M2)\tllh(M3)\tllh(M4)\tllh(M5)\tb1(M1)\tb2(M1)\tSE(b1(M1))\tSE(b2(M1))\tb1(M2)\tSE(b1(M2))\tb2(M3)\tSE(b2(M3))\tb(M4)\tSE(b(M4))\n");            
      } else if(useM0R0){
	ksprintf(&p->bufstr,"Chromo\tPosition\tnInd\tf1\tf2\tllh(M0)\tllh(M1)\tllh(M2)\tllh(M3)\tllh(M4)\tllh(M5)\tb1(M1)\tb2(M1)\tb1(M2)\tb2(M3)\tb(M4)\n");      
      } else{
	ksprintf(&p->bufstr,"Chromo\tPosition\tnInd\tf1\tf2\tllh(M1)\tllh(M2)\tllh(M3)\tllh(M4)\tllh(M5)\tb1(M1)\tb2(M1)\tb1(M2)\tb2(M3)\tb(M4)\n");
      }
    }

    // also has to have intercept
    p->design=initMatrix(4*l,ncov+3);
    
    // starting values for 3 coefs and covar + sd at the end (index: ncov+4)
    p->start = new double[ncov+4];
    p->start0 = new double[ncov+4];
    p->SE = new double[ncov+4];
    
  } else{
    // for the add model where design matrix has R1, R2, Rm, delta1, delta2, covs (delta: rec for other allele - A)
    // start has the same + sd(y) - so one longer than row in design matrix

    if(header){
      if(estSE and useM0R0){
	ksprintf(&p->bufstr,"Chromo\tPosition\tnInd\tf1\tf2\tllh(R0)\tllh(R1)\tllh(R2)\tllh(R3)\tllh(R4)\tllh(R5)\tllh(R6)\tllh(R7)\tb1(R1)\tb2(R1)\tbm(R1)\tSE(b1(R1))\tSE(b2(R1))\tSE(bm(R1))\tb1(R2)\tb2m(R2)\tSE(b1(R2))\tSE(b2m(R2))\tb1m(R3)\tb2(R3)\tSE(b1m(R3))\tSE(b2(R3))\tb1(R4)\tSE(b1(R4))\tb2(R5)\tSE(b2(R5))\tb(R6)\tSE(b(R6))\n");
      } else if(estSE){
	ksprintf(&p->bufstr,"Chromo\tPosition\tnInd\tf1\tf2\tllh(R1)\tllh(R2)\tllh(R3)\tllh(R4)\tllh(R5)\tllh(R6)\tllh(R7)\tb1(R1)\tb2(R1)\tbm(R1)\tSE(b1(R1))\tSE(b2(R1))\tSE(bm(R1))\tb1(R2)\tb2m(R2)\tSE(b1(R2))\tSE(b2m(R2))\tb1m(R3)\tb2(R3)\tSE(b1m(R3))\tSE(b2(R3))\tb1(R4)\tSE(b1(R4))\tb2(R5)\tSE(b2(R5))\tb(R6)\tSE(b(R6))\n");
      } else if(useM0R0){
	ksprintf(&p->bufstr,"Chromo\tPosition\tnInd\tf1\tf2\tllh(R0)\tllh(R1)\tllh(R2)\tllh(R3)\tllh(R4)\tllh(R5)\tllh(R6)\tllh(R7)\tb1(R1)\tb2(R1)\tbm(R1)\tb1(R2)\tb2m(R2)\tb1m(R3)\tb2(R3)\tb1(R4)\tb2(R5)\tb(R6)\n");
      } else{
	ksprintf(&p->bufstr,"Chromo\tPosition\tnInd\tf1\tf2\tllh(R1)\tllh(R2)\tllh(R3)\tllh(R4)\tllh(R5)\tllh(R6)\tllh(R7)\tb1(R1)\tb2(R1)\tbm(R1)\tb1(R2)\tb2m(R2)\tb1m(R3)\tb2(R3)\tb1(R4)\tb2(R5)\tb(R6)\n");
      }
    }
        
    // also has to have intercept
    p->design=initMatrix(4*l,ncov+5);    
    // this one has one starting guess for 5 coefs and covar + sd at the end, index ncov+6
    p->start = new double[ncov+6];
    p->start0 = new double[ncov+6];
    p->SE = new double[ncov+6];
    
  }
  
  p->ysCgs = initMatrix(l,ncov+6);//now we are just allocating enough enough
  p->maxIter = maxInter;
  p->tol = tole;
  
  //plugin a start
  for(int i=0;i<start.size();i++)
    p->start[i] = start[i];
  if(start.size()==0){
    if(model==0){
      rstart(p->start,ncov+3,phe,-1,1);//<-will put sd at p->start[p->covs+dy+1]
      //if logisitc regression set last start value to 0
      if(p->regression==1){
	p->start[ncov+3]=0;
      }
      
    } else{
      rstart(p->start,ncov+5,phe,-1,1);//<-will put sd at p->start[p->covs+dy+1]
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

void set_pars(pars *p, char *g,const std::vector<double> &phe, const std::vector<double> &ad, double *freq, std::vector<double> start, Matrix<double> &cov, char *site, int y){
  
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
  p->index = y;  
  
  for(int i=0;i<p->len;i++){
    for(int j=0;j<4;j++){
      p->pheno[i*4+j] = p->ys[i];
      p->p_sCg[i*4+j] = NAN;
      p->weights[i*4+j] = NAN;
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

void *main_analysis_thread(void* threadarg){

  worker_args * td;
  td = ( worker_args * ) threadarg;

  fprintf(stderr,"Thread %i started\n",td->i);
  
  for(int i=td->first;i<td->second;i++){

    //if last block has sites with values larger than nSites
    if(i>=td->nSites){
      pthread_exit(NULL);
    }
    
    // parsing site with index 100 (0-indexed so 101th site) - but will have parsed 100 sites
    if(parsedThreadSites % 100==0){
      fprintf(stderr,"Parsed sites:%i\n",parsedThreadSites);
    }
    parsedThreadSites++;
    
               
    int cats2[4] = {0,0,0,0};
    
    for(int x=0;x<td->p->len;x++){//similar to above but with transposed plink matrix
      cats2[td->d[i][x]]++;
    }
        
    if((cats2[3]/(double)(cats2[0]+cats2[1]+cats2[2]))>0.1){
      fprintf(stderr,"skipping site[%d] due to excess missingness\n",i);
      continue;
    }
        
    //check that we observe atleast 10 obs of 2 different genotypes
    int n=0;
    
    // changed i to 3, as we do not want to consider missing geno a different genotype
    for(int j=0;j<3;j++)
      if(cats2[j]>1)
	n++;    
    
    // checks that there are at least 2 different genos (1 might be missing)
    if(n<2){
      fprintf(stderr,"skipping site[%d] due to categories filter\n",i);
      continue;
    }
        
    set_pars(td->p,td->d[i],td->phe,td->ad,td->freq.d[i],td->start,td->cov,td->loci[i],td->i);
       
    main_analysis((void*)td->p);
    
    fprintf(td->tmpFile,"%s\t%s\n",td->p->bufstr.s,td->p->tmpstr.s);
    td->p->bufstr.l=td->p->tmpstr.l=0;

    char buf[4096];
    snprintf(buf,4096,"%s\t%s",td->p->bufstr.s,td->p->tmpstr.s);

    myres tmp;
    
    tmp.index=i;
    tmp.res=strdup(buf);
      
    // push_back
    myresults.push_back(tmp);

    analysedThreadSites++;
    
  }
 
  pthread_exit(0);
      
}


void wrap(const plink *plnk,const std::vector<double> &phe,const std::vector<double> &ad,Matrix<double> &freq,int model,std::vector<double> start,Matrix<double> &cov,int maxIter,double tol,std::vector<char*> &loci,int nThreads,FILE* outFile, FILE* logFile, int regression, int estSE, int useM0R0,char* outname){

  char **d = new char*[plnk->y];//transposed of plink->d. Not sure what is best, if we filterout nonmissing anyway.
  for(int i=0;i<plnk->y;i++){
    d[i] = new char[plnk->x];
    for(int j=0;j<plnk->x;j++)
      d[i][j]=plnk->d[j][i];
  }

 
  if(nThreads==1){
    
    //we prep for threading. By using encapsulating all data need for a site in struct called pars
    pars *p=init_pars(plnk->x,cov.dy,model,maxIter,tol,start,phe,regression,logFile,estSE,useM0R0,1);
    
    // check that pheno and covariates are OK
    check_pars(p, phe, cov);
    
    // print starting values!!
    if(p->model==0){
      fprintf(stderr,"Starting values are:\t");
      fprintf(logFile,"Starting values are:\t");
      // starting values equal to covs + 4 (3 for add columns + sd(y) column)
      for(int i=0;i<p->ncov+4;i++){
	fprintf(stderr, "%f ",p->start[i]);
	fprintf(logFile, "%f ",p->start[i]);
      }
      fprintf(stderr,"\n"); fprintf(logFile,"\n");      
    } else{
      // if rec model - copy all starting values ncov + 6
      fprintf(stderr,"Starting values are:\t");
      fprintf(logFile,"Starting values are:\t");
      // starting values equal to covs + 6 (5 for rec columns + sd(y) column)
      for(int i=0;i<p->ncov+6;i++){
	fprintf(stderr, "%f ",p->start[i]);
	fprintf(logFile, "%f ",p->start[i]);
      }
      fprintf(stderr,"\n"); fprintf(logFile,"\n");      
    }
    
    // counts number of sites actually analysed
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
      
      set_pars(p,d[y],phe,ad,freq.d[y],start,cov,loci[y],y);
      
      main_analysis((void*)p);
      analysedSites++;
      fprintf(outFile,"%s\t%s\n",p->bufstr.s,p->tmpstr.s);
      p->bufstr.l=p->tmpstr.l=0;
    }

    // write how many sites analysed 
    fprintf(stderr,"Number of analysed sites is:\t %i\n",analysedSites);
    fprintf(logFile,"Number of analysed sites is:\t %i\n",analysedSites);
    
    fprintf(stderr,"\t-> done\n");
    kill_pars(p,plnk->x);
    
    for(int i=0;i<plnk->y;i++)
      delete [] d[i];
    delete [] d;
    
    
  } else{
    
    int NumJobs = plnk->y;
    // has to read all parameters/data into this struct
    
    double blockd = NumJobs/(double) nThreads;
    // so the blocks are at least a little too big,
    // meaning last block might be a bit smaller
    int block = ceil(blockd);
    
    pars **all_pars = new pars*[nThreads];

    for(int i=0;i<nThreads;i++){
      if(i==0){
	all_pars[i]=init_pars(plnk->x,cov.dy,model,maxIter,tol,start,phe,regression,logFile,estSE,useM0R0,1);
      } else{
	all_pars[i]=init_pars(plnk->x,cov.dy,model,maxIter,tol,start,phe,regression,logFile,estSE,useM0R0,0);
      }

    }

    worker_args **all_args = new worker_args*[nThreads];
    
    int first = 0;
    for(int i=0;i<nThreads;i++){            
      int second = first + block;
      all_args[i]=new worker_args(first,second,all_pars[i],d,phe,ad,start,freq,cov,loci,i);
      first = first + block;
    }
    
    //make tmp files that each thread will write you
    for(int i=0;i<nThreads;i++){
      char buf[strlen(outname)+20];
      snprintf(buf,strlen(outname)+20,"%s.tmp%d.tmp",outname,i);           
      FILE *fp=NULL;
      fp=fopen(buf,"wb");
      assert(fp!=NULL);
      all_args[i]->tmpFile=fp;
      all_args[i]->tmpFileName=strdup(buf);      
    }
    
    pthread_t thread1[nThreads];

    
    // prints starting values
    if(all_pars[0]->model==0){
      fprintf(stderr,"Starting values are:\t");
      fprintf(logFile,"Starting values are:\t");
      // starting values equal to covs + 4 (3 for add columns + sd(y) column)
      for(int i=0;i<all_pars[0]->ncov+4;i++){
	fprintf(stderr, "%f ",all_pars[0]->start[i]);
	fprintf(logFile, "%f ",all_pars[0]->start[i]);
      }
      fprintf(stderr,"\n"); fprintf(logFile,"\n");      
    } else{
      // if rec model - copy all starting values ncov + 6
      fprintf(stderr,"Starting values are:\t");
      fprintf(logFile,"Starting values are:\t");
      // starting values equal to covs + 6 (5 for rec columns + sd(y) column)
      for(int i=0;i<all_pars[0]->ncov+6;i++){
	fprintf(stderr, "%f ",all_pars[0]->start[i]);
	fprintf(logFile, "%f ",all_pars[0]->start[i]);
      }
      fprintf(stderr,"\n"); fprintf(logFile,"\n");      
    }

    
    // this has to be threaded version of main_analysis
    for (int i = 0; i < nThreads; i++)
      assert(pthread_create(&thread1[i], NULL, &main_analysis_thread, all_args[i])==0);
    
    // Wait all threads to finish
    for (int i = 0; i < nThreads; i++)
      assert(pthread_join(thread1[i], NULL)==0);
    
    std::sort(myresults.begin(),myresults.end(),myresSort);
    
    for(int i=0;i<myresults.size();i++){
      fprintf(outFile,"%s\n",myresults[i].res);
    }
    
    // write how many sites analysed 
    fprintf(stderr,"Number of analysed sites is:\t %i\n",analysedThreadSites);
    fprintf(logFile,"Number of analysed sites is:\t %i\n",analysedThreadSites);
   
    for(int i=0;i<nThreads;i++){
      fclose(all_args[i]->tmpFile);
      unlink(all_args[i]->tmpFileName);
      free(all_args[i]->tmpFileName);
    }
    
  }
                  
 
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

