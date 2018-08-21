#include <cstdio>
#include <cmath>
#include <cassert>
#include "asaMap.h"
#include "anal.h"
#include "readplink.h"
#include <cstring>
#include <cfloat>
#include "kstring.h"

void getFit(double *start,double *Y,double **covMatrix,double *weights, int nInd4,int nEnv,double *residuals,int df){
   //  fprintf(stderr,"%s: nInd4:%d nEnv:%d df:%d\n",__FUNCTION__,nInd4,nEnv,df);

   /*
     linear regression. Fits a linear model with weights
     Y is the responce (nInd)
     covMatrix is the design matrix [nInd][nEnv]
     weights is the weights (nInd*nEnv)
     residuals[nInd] are the residuals
     int df is the number of degrees of freedom
     nInd is the number of individuals
     nEnv is the number of predictors (including the intersept)
   */
   int nIndW=0;
   char nonZeroWeight[nInd4];
   memset(nonZeroWeight,0,nInd4);
   if(weights==NULL){
     nIndW=nInd4;
     memset(nonZeroWeight,1,nInd4);
   } else{
     for(int i=0;i<nInd4;i++){
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

   /*
   // emil
   for(int i=0;i<nInd4;i++){
     for(int j=0;j<nEnv;j++){
       fprintf(stdout,"\nxw\t%f\n",xw[cnt*nEnv+j]);
     }
   }
   */

   // This is (X)^T*W*X when doing weighted least squares, where X is the design matrix and W is the weights
   double XtX[nEnv*nEnv];
   for(int i=0;i<nEnv*nEnv;i++)
     XtX[i]=0;

   // this is doing the matrix product of (X)^T*W*X 
   for(int x=0;x<nEnv;x++)
     for(int y=0;y<nEnv;y++)
       for(int i=0;i<nIndW;i++)
	 XtX[x*nEnv+y]+=xw[i*nEnv+x]*xw[i*nEnv+y];

   // emil - begin

   //fprintf(stderr,"weights: %f %f %f %f %f %f %f %f\n",weights[0],weights[1],weights[2],weights[3],weights[4],weights[5],weights[6],weights[7]);

   
   fprintf(stderr,"dsg: %f %f %f %f %f %f %f %f\n",covMatrix[0][0],covMatrix[0][1],covMatrix[0][2],covMatrix[0][3],covMatrix[0][4],covMatrix[0][5],covMatrix[0][6],covMatrix[0][7]);

   fprintf(stderr,"dsg: %f %f %f %f %f %f %f %f\n",covMatrix[1][0],covMatrix[1][1],covMatrix[1][2],covMatrix[1][3],covMatrix[1][4],covMatrix[1][5],covMatrix[1][6],covMatrix[1][7]);
   
   fprintf(stderr,"dsg: %f %f %f %f %f %f %f %f\n",covMatrix[2][0],covMatrix[2][1],covMatrix[2][2],covMatrix[2][3],covMatrix[2][4],covMatrix[2][5],covMatrix[2][6],covMatrix[2][7]);

   fprintf(stderr,"dsg: %f %f %f %f %f %f %f %f\n",covMatrix[3][0],covMatrix[3][1],covMatrix[3][2],covMatrix[3][3],covMatrix[3][4],covMatrix[3][5],covMatrix[3][6],covMatrix[3][7]);


   fprintf(stderr,"dsg: %f %f %f %f %f %f %f %f\n",covMatrix[2573][0],covMatrix[2573][1],covMatrix[2573][2],covMatrix[2573][3],covMatrix[2573][4],covMatrix[2573][5],covMatrix[2573][6],covMatrix[2573][7]);

   fprintf(stderr,"dsg: %f %f %f %f %f %f %f %f\n",covMatrix[2574][0],covMatrix[2574][1],covMatrix[2574][2],covMatrix[2574][3],covMatrix[2574][4],covMatrix[2574][5],covMatrix[2574][6],covMatrix[2574][7]);
   
   fprintf(stderr,"dsg: %f %f %f %f %f %f %f %f\n",covMatrix[2575][0],covMatrix[2575][1],covMatrix[2575][2],covMatrix[2575][3],covMatrix[2575][4],covMatrix[2575][5],covMatrix[2575][6],covMatrix[2575][7]);

   fprintf(stderr,"dsg: %f %f %f %f %f %f %f %f\n",covMatrix[2576][0],covMatrix[2576][1],covMatrix[2576][2],covMatrix[2576][3],covMatrix[2576][4],covMatrix[2576][5],covMatrix[2576][6],covMatrix[2576][7]);

   
   //#if 0
   //print before inversion
   fprintf(stderr,"BEFORE:\n");
   for(int i=0;i<nEnv;i++){
     for(int j=0;j<nEnv;j++)
       fprintf(stderr,"%f ",XtX[i*nEnv+j]);
     fprintf(stderr,"\n");
   } 
   //#endif
   
   double workspace[2*nEnv];   
   matinv(XtX, nEnv, nEnv, workspace);

   //#if 0
   //print after inversion
   fprintf(stderr,"AFTER:\n");
   for(int i=0;i<nEnv;i++){
     for(int j=0;j<nEnv;j++)
       fprintf(stderr,"%f ",XtX[i*nEnv+j]);
     fprintf(stderr,"\n");
   } 
   //#endif
   
   // emil - end

   double Xt_y[nEnv];
   double invXtX_Xt_y[nEnv];
   for(int i=0;i<nEnv;i++)
     Xt_y[i]=0;
   for(int i=0;i<nEnv;i++)
     invXtX_Xt_y[i]=0;

   for(int x=0;x<nEnv;x++)
     for(int i=0;i<nIndW;i++)
       Xt_y[x]+=xw[i*nEnv+x]*yw[i];

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
     double tmp = ((residuals[i] = Y[i]-yTilde[i])); 
     //    fprintf(stdout,"\nYY[%d]\t%f\n",i,tmp);
     if(weights!=NULL)
       ts += tmp*tmp*weights[i];
     else
       ts += tmp*tmp;
   }

   //emil
   fprintf(stderr,"start getFit: %f %f %f %f %f %f %f %f\n",start[0],start[1],start[2],start[3],start[4],start[5],start[6],start[7]);
   
   
   //fprintf(stderr,"ts:%f\n",ts);
   if(df==-1)
     start[nEnv] = sqrt(ts/(1.0*(nInd4-nEnv)));
   else
     start[nEnv] = sqrt(ts/(1.0*df));
 }


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

double tnorm(double x,double mean,double sd){
  
  double fac = 1.0/(sd*sqrt(2.0*M_PI));
  double val = exp(-(((x-mean)*(x-mean))/(2*sd*sd)));

  return fac*val;
}

double logLike(double *start,double* pheno,Matrix<double> *design,double *p_sCg){
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

      // to not include missing genotypes
      if(j==0 && design->d[i][j]==3){
	continue;
      }

      if(i==0){ fprintf(stderr,"design: %f, start: %f\n", design->d[0][j],start[j]);}
      if(i==1){ fprintf(stderr,"design: %f, start: %f\n", design->d[1][j],start[j]);}
      if(i==2){ fprintf(stderr,"design: %f, start: %f\n", design->d[2][j],start[j]);}
      if(i==3){ fprintf(stderr,"design: %f, start: %f\n", design->d[3][j],start[j]);}


      if(i==4){ fprintf(stderr,"design: %f, start: %f\n", design->d[i][j],start[j]);}
      if(i==5){ fprintf(stderr,"design: %f, start: %f\n", design->d[i][j],start[j]);}
      if(i==6){ fprintf(stderr,"design: %f, start: %f\n", design->d[i][j],start[j]);}
      if(i==7){ fprintf(stderr,"design: %f, start: %f\n", design->d[i][j],start[j]);}
     
      
      m +=  design->d[i][j]*start[j];      
    }

    if(i<4){  fprintf(stderr,"mean %f, start[design->dy] %f\n",m,start[design->dy]);}
    m = tnorm(pheno[i],m,start[design->dy]);

    // emil
    
    
    tmp += m*p_sCg[i];
    if((i % 4 )==3){
      ret -= log(tmp);
      tmp = 0;
    }
  }


  
  return ret;
}

inline double logLikeP(pars *p){
  return logLike(p->start,p->pheno,p->design,p->p_sCg);
}

//double updateEM(pars *p){
double updateEM(double *start,Matrix<double> *design,Matrix<double> *ysCgs,double *pheno,int nInd,double *p_sCg ){
#if 0
  for(int i=0;i<=design->dy;i++)
    fprintf(stderr,"%f ",start[i]);
  fprintf(stderr,"\n");
#endif 
  
  double ret;
  for(int i=0;i<design->dx;i++){
    double m=0;
    for(int j=0;j<design->dy;j++)
      m += design->d[i][j]*start[j];
    
    m = tnorm(pheno[i],m,start[design->dy]);    
    double tmp = m*p_sCg[i];
    //    fprintf(stderr,"(%lu,%d):%f\n",(size_t)floor(i/4),i %4,tmp);
    ysCgs->d[(size_t)floor(i/4)][i %4 ] = tmp;
  }

  // emil
  fprintf(stderr,"p_sCg: %f %f %f %f\n",p_sCg[0],p_sCg[1],p_sCg[2],p_sCg[3]);
  fprintf(stderr,"ysCgs->d: %f %f %f %f\n",ysCgs->d[0][0],ysCgs->d[0][1],ysCgs->d[0][2],ysCgs->d[0][3]);
 
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
  int a=0;
  for(int i=0;i<ysCgs->dx;i++)
    for(int j=0;j<ysCgs->dy;j++){
      weigths[a++] = ysCgs->d[i][j];
      // fprintf(stdout,"%f ",weigths[a-1]);
    }
  
  double resi[nInd];
  getFit(start,pheno,design->d,weigths,nInd*4,design->dy,resi,df);
  
#if 0
  for(int i=0;i<=design->dy;i++)
    fprintf(stderr,"%f ",start[i]);
  fprintf(stderr,"\n");
#endif 
  
  return logLike(start,pheno,design,p_sCg);
}
  
  inline double updateEMP(pars *p){
  
  return updateEM(p->start,p->design,p->ysCgs,p->pheno,p->len,p->p_sCg);
}

#define FOR(i,n) for(i=0; i<n; i++)

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
  //  fprintf(stderr,"ncov:(%lu,%lu) design: (%lu,%lu)\n",p->covs->dx,p->covs->dy,p->design->dx,p->design->dy);
  assert(p->design->my>=p->covs->dy+3);



  // emil
  fprintf(stderr,"geno %i\n",p->gs[0]);

   for(int i=0;i<p->len;i++){
     for(int n=0;n<3;n++){
       for(int j=0;j<4;j++){

	 // design matrix has 4 rows per indi - possible states when given geno - cols allele and covs
	 // for pat choose: model (add/rec), allele (B1,B2,A1), row and then column
	 p->design->d[i*4+j][n] = pat[p->model][n][j][p->gs[i]];
	 //	fprintf(stderr,"model:%d n:%d j:%d gs:%d -> %f\n",p->model,n,j,p->gs[i],p->design->d[i*4+j][n]);
	 //	exit(0);
       }
     }
     // putting in covariates in design matrix, intercept is not added
     for(int c=0;c<p->covs->dy;c++)
       for(int j=0;j<4;j++){
	 //	fprintf(stdout,"(%d,%d):%f\n",i*4+j,2+c,p->covs->d[i][c]);
	 p->design->d[i*4+j][3+c] = p->covs->d[i][c];
       }
     // break;
   }
   p->design->dx=p->len*4;

   //p->design->dy=p->covs->dy+2;
   // emil - now design has 3 columns (alternative allele also) + number of covariates
   p->design->dy=p->covs->dy+3;

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
       }if(llh0<llh1){
	 //	 fprintf(stderr,"Fit caused increase in likelihood, will roll back to previous step\n");
	 memcpy(p->start,pars0,sizeof(double)*(p->design->dy+1));
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
  //      exit(0);
}

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
  // if to rare do not run analysis
  if(p->mafs[0]>0.995||p->mafs[0]<0.005){
    maf0=0;
  }
  if(p->mafs[1]>0.995||p->mafs[1]<0.005){
    maf1=0;
  }

  // p->model is either 0 (add model) or 1 (rec model)
  if(p->model>0){


    //////////// do R0 - model with also R(A1) and R(A2) (rec for alternative allele) ///////////////
    
    fprintf(stderr,"R0\n");
    
    mkDesign(p);
    p_sCg(p);
    if(maf0 && maf1){
      controlEM(p);
      printRes(p,2); 
    } else{
      printNan(p,2);
    }
    
    //////////// do M1 ///////////////
    
    fprintf(stderr,"R1\n");
    
    mkDesign(p);
    p_sCg(p);
    
    memcpy(p->start,p->start0,sizeof(double)*(p->design->dy+1));
    
    // remove fourth column from design - column counting R(A1) 
    rmPos(p->start,3,p->covs->dy+5+1);
    rmCol(p->design,3);

    // remove fifth column from design (now fourth because fifth column was removed) - column counting R(A2)
    rmPos(p->start,3,p->covs->dy+4+1);
    rmCol(p->design,3);
    
    if(maf0 && maf1){
      controlEM(p);
      printRes(p,2); 
    } else{
      printNan(p,2);
    }
    //////////// do M2 ///////////////
    //remove column2 and second value from start M2
    
    fprintf(stderr,"M2\n");
    
    mkDesign(p);
    p_sCg(p);
    memcpy(p->start,p->start0,sizeof(double)*(p->design->dy+1));
    
    // remove B2 second in design, remember start has sd(y) at the end (one longer)
    rmPos(p->start,0,p->covs->dy+3+1);
    
    rmCol(p->design,0);
    // remove A1 third in design (NB!! now second cause first was removed), remember start has sd(y) at the end (one longer)
    rmPos(p->start,1,p->covs->dy+2+1);
    rmCol(p->design,1);
    
    
    if(maf0){
      controlEM(p);
      printRes(p,1); 
    }  else{
      printNan(p,1);
    }
    //////////// do M3 ///////////////
    //remove column1 and first value from start M3
    
    fprintf(stderr,"M3\n");
    
    mkDesign(p);
    p_sCg(p);
    memcpy(p->start,p->start0,sizeof(double)*(p->design->dy+1));
    
    // remove B1 first in design, remember start has sd(y) at the end (one longer)
    rmPos(p->start,1,p->covs->dy+3+1);
    rmCol(p->design,1);
    // remove A1 third in design (NB!! now second cause first was removed), remember start has sd(y) at the end (one longer)
    rmPos(p->start,1,p->covs->dy+2+1);
    rmCol(p->design,1);
    
    if(maf1){
      controlEM(p);
      printRes(p,1); 
    } else{
      printNan(p,1);
    }
    
    //////////// do M4 ///////////////
    //cbind gs and covs into design M4:
    
    fprintf(stderr,"M4\n");
    
    for(int i=0;i<p->covs->dx;i++){
      p->design->d[i][0] = p->gs[i];
      for(int j=0;j<p->covs->dy;j++)
	p->design->d[i][j+1] = p->covs->d[i][j];
    }
    p->design->dx=p->covs->dx;
    p->design->dy=p->covs->dy+1;
    
    p_sCg(p);
    
    // emil - start gets values in getFit - these are the coefs
    for(int i=0;i<p->design->dy+2;i++)
      p->start[i]=NAN;
    
    double resi[p->design->dx];
    getFit(p->start,p->ys,p->design->d,NULL,p->design->dx,p->design->dy,resi,-1);
    
    
    // it somehow needs full design matrix for likelihood calculations
    // WHAT IS THIS???????
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
      fprintf(stderr,"i: %i, p->start[i] %f\n",i,p->start[i]); 
      if(i>2){
	p->start[i] = p->start[i-2];
      }
    }
    p->start[2]=p->start[1]=0;
  
    for(int i=0;i<8;i++){
      fprintf(stderr,"i: %i, p->start[i] %f\n",i,p->start[i]);        
      
    }
    // B2 and A1 has no effect, as they are not in model, B1 is genotype
    
    //p->start[0]=p->start[2];
    
    //print to kbuf
    printRes(p,1); 
    
    //////////// do M5 ///////////////  
    //cpy covs to design M5
    
    fprintf(stderr,"M5\n");
    
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
    getFit(p->start,p->ys,p->design->d,NULL,p->design->dx,p->design->dy,resi,-1);
    mkDesign(p);
    
    
    // skips A1,A2 for me has to skip both A1,A2 and B1 (column 1,2 and 3 design matrix)
    for(int i=p->design->dy;i>1;i--){
      fprintf(stderr,"i: %i, p->start[i] %f\n",i,p->start[i]);        
      p->start[i] = p->start[i-3];
    }
    p->start[2]=p->start[1]=p->start[0]=0;
    
    for(int i=0;i<8;i++){
      fprintf(stderr,"i: %i, p->start[i] %f\n",i,p->start[i]);        
    }
    
    printRes(p,0); 
    
    


  } else{
    //////////// do M0 - model with also A1 ///////////////
    
    fprintf(stderr,"M0\n");
    
    mkDesign(p);
    p_sCg(p);
    if(maf0 && maf1){
      controlEM(p);
      printRes(p,2); 
    } else{
      printNan(p,2);
    }
    
    //////////// do M1 ///////////////
    
    fprintf(stderr,"M1\n");
    
    mkDesign(p);
    p_sCg(p);
    
    memcpy(p->start,p->start0,sizeof(double)*(p->design->dy+1));
    
    // remove third column from design - column counting A1, remember start has sd(y) at the end (one longer)
    rmPos(p->start,2,p->covs->dy+3+1);
    rmCol(p->design,2);
    
    if(maf0 && maf1){
      controlEM(p);
      printRes(p,2); 
    } else{
      printNan(p,2);
    }
    //////////// do M2 ///////////////
    //remove column2 and second value from start M2
    
    fprintf(stderr,"M2\n");
    
    mkDesign(p);
    p_sCg(p);
    memcpy(p->start,p->start0,sizeof(double)*(p->design->dy+1));
    
    // remove B2 second in design, remember start has sd(y) at the end (one longer)
    rmPos(p->start,0,p->covs->dy+3+1);
    
    rmCol(p->design,0);
    // remove A1 third in design (NB!! now second cause first was removed), remember start has sd(y) at the end (one longer)
    rmPos(p->start,1,p->covs->dy+2+1);
    rmCol(p->design,1);
    
    
    if(maf0){
      controlEM(p);
      printRes(p,1); 
    }  else{
      printNan(p,1);
    }
    //////////// do M3 ///////////////
    //remove column1 and first value from start M3
    
    fprintf(stderr,"M3\n");
    
    mkDesign(p);
    p_sCg(p);
    memcpy(p->start,p->start0,sizeof(double)*(p->design->dy+1));
    
    // remove B1 first in design, remember start has sd(y) at the end (one longer)
    rmPos(p->start,1,p->covs->dy+3+1);
    rmCol(p->design,1);
    // remove A1 third in design (NB!! now second cause first was removed), remember start has sd(y) at the end (one longer)
    rmPos(p->start,1,p->covs->dy+2+1);
    rmCol(p->design,1);
    
    if(maf1){
      controlEM(p);
      printRes(p,1); 
    } else{
      printNan(p,1);
    }
    
    //////////// do M4 ///////////////
    //cbind gs and covs into design M4:
    
    fprintf(stderr,"M4\n");
    
    for(int i=0;i<p->covs->dx;i++){
      p->design->d[i][0] = p->gs[i];
      for(int j=0;j<p->covs->dy;j++)
	p->design->d[i][j+1] = p->covs->d[i][j];
    }
    p->design->dx=p->covs->dx;
    p->design->dy=p->covs->dy+1;
    
    p_sCg(p);
    
    // emil - start gets values in getFit - these are the coefs
    for(int i=0;i<p->design->dy+2;i++)
      p->start[i]=NAN;
    
    double resi[p->design->dx];
    getFit(p->start,p->ys,p->design->d,NULL,p->design->dx,p->design->dy,resi,-1);
    
    
    // it somehow needs full design matrix for likelihood calculations
    // WHAT IS THIS???????
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
      fprintf(stderr,"i: %i, p->start[i] %f\n",i,p->start[i]); 
      if(i>2){
	p->start[i] = p->start[i-2];
      }
    }
    p->start[2]=p->start[1]=0;
  
    for(int i=0;i<8;i++){
      fprintf(stderr,"i: %i, p->start[i] %f\n",i,p->start[i]);        
      
    }
    // B2 and A1 has no effect, as they are not in model, B1 is genotype
    
    //p->start[0]=p->start[2];
    
    //print to kbuf
    printRes(p,1); 
    
    //////////// do M5 ///////////////  
    //cpy covs to design M5
    
    fprintf(stderr,"M5\n");
    
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
    getFit(p->start,p->ys,p->design->d,NULL,p->design->dx,p->design->dy,resi,-1);
    mkDesign(p);
    
    
    // skips A1,A2 for me has to skip both A1,A2 and B1 (column 1,2 and 3 design matrix)
    for(int i=p->design->dy;i>1;i--){
      fprintf(stderr,"i: %i, p->start[i] %f\n",i,p->start[i]);        
      p->start[i] = p->start[i-3];
    }
    p->start[2]=p->start[1]=p->start[0]=0;
    
    for(int i=0;i<8;i++){
      fprintf(stderr,"i: %i, p->start[i] %f\n",i,p->start[i]);        
    }
    
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
