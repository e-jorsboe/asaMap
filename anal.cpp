#include <cmath>
#include "readplink.c"
#include "anal.h"
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



pars *init_pars(size_t l,size_t ncov,int model,int maxInter,double tole,std::vector<double> &start){
  pars *p=new pars;
  p->len=l;
  p->ncov=ncov;
  p->gs=new char[l];
  p->ys=new double[l];
  p->qvec=new double[l];
  //p->mafs=new double[l];
  p->covs=initMatrix(l,ncov);

  // emil - added 3rd column for alternate allele
  //p->design=initMatrix(4*l,ncov+2);
  p->design=initMatrix(4*l,ncov+3);
  p->pheno = new double[4*l];
  p->p_sCg = new double[4*l];
  p->bufstr.s=NULL;p->bufstr.l=p->bufstr.m=0;
  p->tmpstr.s=NULL;p->tmpstr.l=p->tmpstr.m=0;
  
  ksprintf(&p->bufstr,"Chromo\tPosition\tnInd\tf1\tf2\tllh(M0)\tllh(M1)\tllh(M2)\tllh(M3)\tllh(M4)\tllh(M5)\tb1(M1)\tb2(M1)\tb1(M2)\tb2(M3)\tb(M4)\n");
  
  //branch
  p->model = model;
  // this one has one starting guess for 3 coefs and covar + sd at the end, index ncov+4
  p->start = new double[ncov+4];
  p->start0 = new double[ncov+4];
  p->ysCgs = initMatrix(l,ncov+2+4);//now we are just allocating enough enough
  p->maxIter = maxInter;
  p->tol = tole;

  //plugin a start
  for(int i=0;i<start.size();i++)
    p->start[i] = start[i];
  if(start.size()==0){
    //make better start guess at some point
    void rstart(double *,size_t);
    // emil - will put in random numbers as starting guess
    rstart(p->start,ncov+3);//<-will put sd at p->start[p->covs+dy+2]
  }
  //copy it to the start0 which will be copied to start for each new site
  // emil has to be + 4, because also needs to copy sd(y) - last value of start
  memcpy(p->start0,p->start,sizeof(double)*(ncov+4));
  return p;
}

double sd(double *a,int l){
  double ts =0;
  for(int i=0;i<l;i++)
    ts += a[i];

  double u = ts/(1.0*l);
  assert(u!=0);
  ts = 0;
  for(int i=0;i<l;i++)
    ts += (a[i]-u)*(a[i]-u);
  return ts/(1.0*(l-1.0));
}

void rstart(double *ary,size_t l){
  for(int i=0;i<l;i++){
   ary[i] = drand48()*2-1;
   //  fprintf(stderr,"ary[%d]:%f\n",i,ary[i]);
  }
  // emil - for sd estimate use sd of y
  ary[l]=sd(ary,l);
}


void set_pars(pars*p,char *g,const std::vector<double> &phe,const std::vector<double> &ad , double *freq,std::vector<double> start,Matrix<double> &cov,char *site){

  p->len=0;
  for(int i=0;i<phe.size();i++){
      if(g[i]!=3){
      p->gs[p->len] = 2-g[i];//DRAGON 
      p->ys[p->len] = phe[i];

      p->qvec[p->len] = ad[i];
      for(int c=0;c<cov.dy;c++)
	p->covs->d[p->len][c] = cov.d[i][c]; 
      p->len++;
    }
  }
  
  p->covs->dx=p->len;
  p->covs->dy=cov.dy;
  p->mafs = freq;
  
  for(int i=0;i<p->len;i++){
    for(int j=0;j<4;j++){
      p->pheno[i*4+j] = p->ys[i];
      p->p_sCg[i*4+j] = NAN;
    }

  }

  memcpy(p->start,p->start0,sizeof(double)*(p->covs->dy+3));
  ksprintf(&p->bufstr,"%s%d\t%f\t%f\t",site,p->len,p->mafs[0],p->mafs[1]);
}


void wrap(const plink *plnk,const std::vector<double> &phe,const std::vector<double> &ad,Matrix<double> &freq,int model,std::vector<double> start,Matrix<double> &cov,int maxIter,double tol,std::vector<char*> &loci,int nThreads,FILE *outFile){
  //fprintf(stderr,"\t-> plinkdim: x->%lu y->%lu\n",plnk->x,plnk->y);
  //return ;
  char **d = new char*[plnk->y];//transposed of plink->d. Not sure what is best, if we filterout nonmissing anyway.
  for(int i=0;i<plnk->y;i++){
    d[i] = new char[plnk->x];
    for(int j=0;j<plnk->x;j++)
      d[i][j]=plnk->d[j][i];
  }


  
  pars *p=init_pars(plnk->x,cov.dy,model,maxIter,tol,start);//we prep for threading. By using encapsulating all data need for a site in  struct called pars

  for(int y=0;y<plnk->y;y++){//loop over sites

    // emil
    if(y>0){ exit(0); }
    
    fprintf(stderr,"Parsing site:%d\r",y);
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
    for(int i=0;i<4;i++)
      if(cats2[i]>1)
	n++;
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
    fprintf(outFile,"%s\t%s\n",p->bufstr.s,p->tmpstr.s);
    p->bufstr.l=p->tmpstr.l=0;
    //break;

  }
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

