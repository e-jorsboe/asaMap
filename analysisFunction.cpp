#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include <list>
#include <assert.h>
#include <fstream>
#include "analysisFunction.h"

#define LENS 4096

#define LOWER_BOUND 1e-20

// checks if string is number
int validNumber(char* someString){
  int isNumber = 0;
  int hasPoint = 0;
  int hase = 0;
   
  for(int i = 0; i<strlen(someString); i++){
    char s = someString[i];

    //if negative number
    if(i==0 and s=='-'){
      continue;
    } else if(i==0 and s=='e'){
      break;
    }
    
    if(s =='.'){
      if(hasPoint){
	isNumber = 0;
	break;
      }
      hasPoint = 1;
      isNumber = 1;
      // if number with e in
    } else if(s =='e'){
      if(hase){
	isNumber = 0;
	break;
      }
      hase = 1;
      isNumber = 1;      	
    } else if(i>0 and s  == '-'){
      if(someString[i-1] == 'e'){
	continue;
      }
    } else {
      isNumber = (s=='0' or s=='1' or s=='2' or s=='3' or s=='4' or s=='5' or s=='6' or s=='7' or s=='8' or s=='9');
    }

    // leaves loop if any not numbers
    if(isNumber==0){

      fprintf(stderr,"k %s\n",someString);
      break;
    }
  }
  
  return(isNumber);
  
}


Matrix<double> *initMatrix(size_t x,size_t y){
  Matrix<double> *r = new Matrix<double>;
  r->d = new double*[x];
  for(size_t i=0;i<x;i++)
    r->d[i] = new double[y];
  r->mx=x;
  r->my=y;
  r->dx=r->dy =0;
  return r;
}

void kill(Matrix<double> *rr){
  for(size_t i=0;i<rr->mx;i++)
    delete [] rr->d[i];
  delete [] rr->d;
  delete rr;
  rr=NULL;
}

double addProtect2(double a,double b){
  //function does: log(exp(a)+exp(b)) while protecting for underflow
  double maxVal;// = std::max(a,b));
  if(a>b)
    maxVal=a;

  else
    maxVal=b;
  double sumVal = exp(a-maxVal)+exp(b-maxVal);
  return log(sumVal) + maxVal;
}


double addProtect3(double a,double b, double c){
  //function does: log(exp(a)+exp(b)+exp(c)) while protecting for underflow
  double maxVal;// = std::max(a,std::max(b,c));
  if(a>b&&a>c)
    maxVal=a;
  else if(b>c)
    maxVal=b;
  else
    maxVal=c;
  double sumVal = exp(a-maxVal)+exp(b-maxVal)+exp(c-maxVal);
  return log(sumVal) + maxVal;
}

double getMax(double a,double b, double c){ 
    //get the maximum value of a, b and c
    double maxVal;// = std::max(a,std::max(b,c));
    if(a>b&&a>c)
      maxVal=a;
    else if(b>c)
      maxVal=b;
    else
      maxVal=c;
    return maxVal;
}

std::vector<double> getArray(const char *name){
  std::vector<double> ret;
  if(name==NULL)
    return ret;
 if(!fexists(name)){
    fprintf(stderr,"\t-> Problems opening file: %s\n",name);
    exit(0);
  }
  const char* delims = " \t";
  gzFile gz = Z_NULL;
  char buffer[LENS];
  if(((gz=gzopen(name,"rb")))==Z_NULL){
    fprintf(stderr,"Problem opening file: \'%s\'",name);
    exit(0);
  }

  while(gzgets(gz,buffer,LENS)){
    char *tok = strtok(buffer,delims);
    while(tok!=NULL){
      ret.push_back(atof(tok));
      tok = strtok(NULL,delims);
    }    
  }
  fprintf(stderr,"Done reading file: \'%s\' containing nitems:%lu \n",name,ret.size());
  gzclose(gz);
  return ret;
}

// checks if read string is a number
std::vector<double> getArrayCheck(const char *name){
  std::vector<double> ret;
  if(name==NULL)
    return ret;
 if(!fexists(name)){
    fprintf(stderr,"\t-> Problems opening file: %s\n",name);
    exit(0);
  }
  const char* delims = " \t\n";
  gzFile gz = Z_NULL;
  char buffer[LENS];
  char *save;
  if(((gz=gzopen(name,"rb")))==Z_NULL){
    fprintf(stderr,"Problem opening file: \'%s\'\n",name);
    exit(0);
  }

  while(gzgets(gz,buffer,LENS)){
    char *tok = strtok_r(buffer,delims,&save);
    while(tok!=NULL){
      if(validNumber(tok)){
	ret.push_back(atof(tok));
	tok = strtok_r(NULL,delims,&save);
      } else{
	fprintf(stderr,"Not valid double in file: \'%s\'\n",name);
	exit(0);
      }
    }    
  }
  fprintf(stderr,"Done reading file: \'%s\' containing nitems:%lu \n",name,ret.size());
  gzclose(gz);
  return ret;
}


Matrix<double> getMatrix(const char *name){
  if(!fexists(name)){
    fprintf(stderr,"\t-> Problems opening file: %s\n",name);
    exit(0);
  }
  const char* delims = " \t";
  gzFile gz = Z_NULL;
  char buffer[LENS];
  if(((gz=gzopen(name,"rb")))==Z_NULL){
    fprintf(stderr,"Problem opening file: \'%s\'",name);
    exit(0);
  }
  std::list<double *> rows;
  size_t ncols =0;
  while(gzgets(gz,buffer,LENS)){
   
    if(strlen(buffer)==0)
      continue;
    char *tok = strtok(buffer,delims);
    std::list<double> items;
    while(tok!=NULL){
      items.push_back(atof(tok));
      tok = strtok(NULL,delims);
    }
    //fprintf(stderr,"[%s] ncols:%lu\n",__FUNCTION__,items.size());
    if(ncols!=0 &ncols!=items.size()){
      fprintf(stderr,"jagged covariance matrix: \'%s\'\n",name);
      exit(0);
    }
    ncols = items.size();
    double *drows = new double[ncols];
    int i=0;
    for(std::list<double>::iterator it=items.begin();it!=items.end();it++)
      drows[i++]  = *it;
    rows.push_back(drows);
    
  }
  //  fprintf(stderr,"%s nrows:%lu\n",__FUNCTION__,rows.size());
  double **data = new double*[rows.size()];
  int i=0;
  for(std::list<double*>::iterator it=rows.begin();it!=rows.end();it++)
    data[i++]  = *it;
  
  Matrix<double> retMat;
  retMat.d=data;
  retMat.dx =retMat.mx = rows.size();
  retMat.dy =retMat.my =  ncols;
  fprintf(stderr,"Done reading file: \'%s\' containing nrows:%lu and ncols:%lu\n",name,retMat.dx,retMat.dy);
  gzclose(gz);
  //  assert(retMat.dx>1&&retMat.dy>1);
  return retMat;

}

Matrix<double> getMatrixCheck(const char *name){
  if(!fexists(name)){
    fprintf(stderr,"\t-> Problems opening file: %s\n",name);
    exit(0);
  }
  const char* delims = " \t\n";
  gzFile gz = Z_NULL;
  char buffer[LENS];
  char* save;
  if(((gz=gzopen(name,"rb")))==Z_NULL){
    fprintf(stderr,"Problem opening file: \'%s\'\n",name);
    exit(0);
  }
  std::list<double *> rows;
  size_t ncols=0;
  while(gzgets(gz,buffer,LENS)){
   
    if(strlen(buffer)==0)
      continue;
    char *tok = strtok_r(buffer,delims,&save);
    std::list<double> items;    
    while(tok!=NULL){
      if(validNumber(tok)){	
	items.push_back(atof(tok));
	tok = strtok_r(NULL,delims,&save);
      } else{
	fprintf(stderr,"Not valid double in file: \'%s\'\n",name);
	exit(0);
      }
    }
    //fprintf(stderr,"[%s] ncols:%lu\n",__FUNCTION__,items.size());
    if(ncols!=0 &ncols!=items.size()){
      fprintf(stderr,"jagged covariance matrix: \'%s\'\n",name);
      exit(0);
    }
    ncols = items.size();
    double *drows = new double[ncols];
    int i=0;
    for(std::list<double>::iterator it=items.begin();it!=items.end();it++)
      drows[i++]  = *it;
    rows.push_back(drows);
    
  }
  //  fprintf(stderr,"%s nrows:%lu\n",__FUNCTION__,rows.size());
  double **data = new double*[rows.size()];
  int i=0;
  for(std::list<double*>::iterator it=rows.begin();it!=rows.end();it++)
    data[i++]  = *it;
  
  Matrix<double> retMat;
  retMat.d=data;
  retMat.dx=retMat.mx = rows.size();
  retMat.dy=retMat.my =  ncols;
  fprintf(stderr,"Done reading file: \'%s\' containing nrows:%lu and ncols:%lu\n",name,retMat.dx,retMat.dy);
  gzclose(gz);
  //  assert(retMat.dx>1&&retMat.dy>1);
  return retMat;

}

void deleteMatrix(Matrix<double> mat){
  assert(mat.d!=NULL);
  for(int i=0;i<mat.dx;i++)
    delete [] mat.d[i];
  delete[] mat.d;
  mat.d =NULL;
}

void print(Matrix<double> *mat,FILE *file){
  fprintf(stderr,"Printing mat:%p with dim=(%lu,%lu)\n",mat->d,mat->dx,mat->dy);
  for(int xi=0;xi<mat->dx;xi++){
    for(int yi=0;yi<mat->dy;yi++)
      fprintf(file,"%f\t",mat->d[xi][yi]);
    fprintf(file,"\n");
  }    
}

void print(char *ary,size_t l,FILE *file,char *he){
  fprintf(stderr,"%s\n",he);
  for(int i=0;i<l;i++)
    fprintf(stderr,"%d ",ary[i]);
  fprintf(stderr,"\n");
  
}

void print(double *ary,size_t l,FILE *file,char *he){
  if(he!=NULL)
    fprintf(stderr,"%s\n",he);
  for(int i=0;i<l;i++)
    fprintf(stderr,"%f ",ary[i]);
  fprintf(stderr,"\n");
  
}
void swapDouble (double& first, double& second)
{
        double temp = first;
        first = second;
        second = temp;
}


void logrescale(double *ary,int len){
  int maxId = 0;
  //    fprintf(stderr,"maxid:%d maxval:%f\n",maxId,ary[maxId]);
  for(int i=1;i<len;i++){
    //if(posCounter==debug_print)
    //fprintf(stderr,"maxid:%d maxval:%f ary=%f\n",maxId,ary[maxId],ary[i]);
    if(ary[i]>ary[maxId])
      maxId=i;
  }
  double maxVal = ary[maxId];
  //if(posCounter==debug_print)
  //  fprintf(stderr,"maxval: %f\n",maxVal);
  for(int i=0;i<len;i++){
    //if(posCounter==debug_print)
    //fprintf(stderr,"%f\t%f\n",ary[i],ary[i]-maxVal);
    ary[i] = ary[i]-maxVal;
  }
  //  exit(0);
}


std::vector<char*> getFilenames(const char * name,int nInd){
 
  if(!fexists(name)){
    fprintf(stderr,"[%s]\t-> Problems opening file: %s\n",__FUNCTION__,name);
    exit(0);
  }
  const char* delims = " \t";
  std::vector<char*> ret;
  std::ifstream pFile(name,std::ios::in);

  char buffer[LENS];
  while(!pFile.eof()){
    pFile.getline(buffer,LENS);
    char *tok = strtok(buffer,delims);
    while(tok!=NULL){
      if(tok[0]!='#')
	ret.push_back(strdup(buffer));
      tok = strtok(NULL,delims);
    }
  }
  if(nInd>0) {
     if(ret.size()<nInd)
      fprintf(stderr,"\t-> Number of samples is smaller than subset requested %lu vs %d\n",ret.size(),nInd);
    else{
      //   fprintf(stderr,"\t-> Will remove tail of filename list\n");
      for(int ii=nInd;ii<ret.size();ii++)
	free(ret[ii]);
      ret.erase(ret.begin()+nInd,ret.end());//we don't free this memory, it doesn't really matter
      // fprintf(stderr,"\t->  filename list now contains only: %lu\n",ret.size());
    }
#if 0
     for(size_t ii=0;ii<ret.size();ii++)
       fprintf(stderr,"%zu->%s\n",ii,ret[ii]);
     fprintf(stderr,"\n");
#endif
  }

  return ret;
}


void print_array(FILE *fp,double *ary,int len){
  for(int i=0;i<len-1;i++)
    fprintf(fp,"%f,",ary[i]);
  fprintf(fp,"%f\n",ary[len-1]);
}

void print_array(FILE *fp,int *ary,int len){
  for(int i=0;i<len-1;i++)
    fprintf(fp,"%d,",ary[i]);
  fprintf(fp,"%d\n",ary[len-1]);
}


double sigm(double x){
  return(1/(1+exp(-x)));
}

double *readDouble(const char*fname,int hint){
  FILE *fp = NULL;
  fp = getFILE(fname,"r");
  char buf[fsize(fname)+1];
  if(fsize(fname)!=fread(buf,sizeof(char),fsize(fname),fp)){
    fprintf(stderr,"Problems reading file: %s\n will exit\n",fname);
    exit(0);
  }
  buf[fsize(fname)]='\0';
  std::vector<double> res;
  res.push_back(atof(strtok(buf,"\t\n ")));
  char *tok=NULL;
  while((tok=strtok(NULL,"\t\n "))) {  
    //fprintf(stderr,"%s\n",tok);
    res.push_back(atof(tok));

  }
  //  fprintf(stderr,"size of prior=%lu\n",res.size());
  if(hint!=res.size()){
    fprintf(stderr,"problem with size of dimension of prior %d vs %lu\n",hint,res.size());
    for(size_t i=0;i<res.size();i++)
      fprintf(stderr,"%zu=%f\n",i,res[i]);
    exit(0);
  }
  double *ret = new double[res.size()];
  for(size_t i=0;i<res.size();i++)
    ret[i] = res[i];
  if(fp) fclose(fp);
  return ret;
}

int whichMax(double *d,int len){
  int r=0;
  for(int i=1;i<len;i++)
    if(d[i]>d[r])
      r=i;
  //now check if site doesnt have data.
  
  if(r==0){//only check if nothing is higher than the first
    for(int i=1;i<len;i++)
      if(d[i]!=d[0])//we see a diffrence so we have information
	return r;
    return -1;//we didnt have information 
  }else
    return r;
}



//usefull little function to split
char *strpop(char **str,char split){
  char *tok=*str;
  while(**str){
    if(**str!=split)
      (*str)++;
    else{
      **str='\0'; (*str)++;
      break;
    }
  }
  return tok;
}



// from ANGSD source code
int ludcmp(double **a, int *indx, double &d,int n){
  int imax = 0;
  double big, dum, sum, temp;
  double vv[n];
  d=1;

  for (int i=0; i<n; i++){
    big=0;
    for (int j=0; j<n; j++){
      //fprintf(stderr,"%f\t",a[i][j]);
      if ((temp=fabs(a[i][j])) > big) 
	big=temp;
    }
    if(big==0){
      fprintf(stderr,"singular matrix in ludcmp\n");
      return(1);
	//    assert(big!=0) ;

    }
    vv[i]=1/big;
  }

  for (int j=0; j<n; j++){
    for (int i=0; i<j; i++){
      sum = a[i][j];
      for (int k=0; k<i; k++) 
	sum -= a[i][k] * a[k][j];
      a[i][j]=sum;
    }
    big=0;
    for (int i=j; i<n; i++)	{
      sum=a[i][j];
      for (int k=0; k<j; k++)
	sum -= a[i][k] * a[k][j];
      a[i][j]=sum;
      if ((dum=vv[i]*fabs(sum)) >= big) {
	big = dum;
	imax = i;
      }
    }
    if (j != imax){
      for (int k=0; k<n; k++){
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      d = -d;
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0) 
      a[j][j] = 1.0e-20;
    if (j != n-1){
      dum = 1/(a[j][j]);
      for (int i=j+1; i<n; i++) 
	a[i][j] *= dum;
    }
  }
  return 0;
}

// from ANGSD source code
void lubksb(double **a, int *indx, double *b,int n){

  int ii=0;
  double sum;

  for (int i=0; i<n; i++){
    int ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii != 0)
      for (int j=ii-1; j<i; j++) 
	sum -= a[i][j]*b[j];
    else if (sum != 0.0) 
      ii=i+1;
    b[i]=sum;
  }
  for (int i=n-1; i>=0; i--){
    sum=b[i];
    for (int j=i+1; j<n; j++) 
      sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}


// from ANGSD source code
int svd_inverse(double mat[],int xLen, int yLen){
  if(xLen !=yLen){

    fprintf(stderr,"non square matrix!\n");
    return 1;
    //exit(0);

  }
  double *col;
  double y[xLen * yLen];
  col = new double[xLen];
  double **tm;
  int *indx=new int[xLen];
  double d;
  tm = new double*[xLen];
  for (int i=0; i < xLen; i++)
    tm[i] = new double[xLen];

  for(int i=0;i<xLen;i++)
    for(int j=0;j<yLen;j++)
      tm[i][j]=mat[j*xLen+i];

  int singular=ludcmp(tm,indx,d,xLen);
  if(singular)
    return 1;
  
  for (int j=0; j<xLen; j++)
    {
      for (int i=0; i<xLen; i++)
	col[i]=0;
      col[j]=1;
      lubksb(tm,indx,col,xLen);
      for (int i=0; i<xLen; i++) 
	y[j*xLen+i]=col[i];
    }
  
  
  for (int j=0; j<yLen; j++)
    for (int i=0; i<xLen; i++)
      mat[j*xLen+i]=y[j*xLen+i];

  delete[] col;
  delete[] indx;
  for (int i=0; i < xLen; i++)
    delete[] tm[i];
  delete[] tm;
  return 0;
}





// a,c,g,t,n
// A,C,G,T,N
// 0,1,2,3,4
int refToInt[256] = {
  0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//15
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//31
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//47
  0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//63
  4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//79
  4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//95
  4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//111
  4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//127
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//143
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//159
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//175
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//191
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//207
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//223
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//239
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4//255
};

char intToRef[5] = {'A','C','G','T','N'};




int fexists(const char* str){///@param str Filename given as a string.
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}

size_t fsize(const char* fname){
  struct stat st ;
  stat(fname,&st);
  return st.st_size;
}

std::vector <char *> dumpedFiles;//small hack for getting a nice vector of outputfiles
FILE *openFile(const char* a,const char* b){
  if(0)
    fprintf(stderr,"[%s] %s %s",__FUNCTION__,a,b);
  char *c = new char[strlen(a)+strlen(b)+1];
  strcpy(c,a);
  strncat(c,b,strlen(b));
  //  fprintf(stderr,"\t-> Dumping file: %s\n",c);
  dumpedFiles.push_back(strdup(c));
  FILE *fp = fopen(c,"w");
  delete [] c;
  return fp;
}
gzFile openFileGz(const char* a,const char* b,const char *mode){
  if(0)
    fprintf(stderr,"[%s] %s%s\n",__FUNCTION__,a,b);
  char *c = new char[strlen(a)+strlen(b)+1];
  strcpy(c,a);
  strncat(c,b,strlen(b));
  //  fprintf(stderr,"\t-> Dumping file: %s\n",c);
  dumpedFiles.push_back(strdup(c));
  gzFile fp = getGz(c,mode);
  delete [] c;
  return fp;
}



FILE *getFILE(const char*fname,const char* mode){
  int writeFile = 0;
  for(size_t i=0;i<strlen(mode);i++)
    if(mode[i]=='w')
      writeFile = 1;
  FILE *fp;
  if(NULL==(fp=fopen(fname,mode))){
    fprintf(stderr,"\t-> Error opening FILE handle for file:%s exiting\n",fname);
    exit(0);
  }
  return fp;
}

gzFile getGz(const char*fname,const char* mode){
  int writeFile = 0;
  for(size_t i=0;i<strlen(mode);i++)
    if(mode[i]=='w')
      writeFile = 1;

  //  fprintf(stderr,"\t-> opening: %s\n",fname);
  gzFile fp=Z_NULL;
  if(NULL==(fp=gzopen(fname,mode))){
    fprintf(stderr,"\t-> Error opening gzFile handle for file:%s exiting\n",fname);
    exit(0);
  }
  return fp;
}



//checks that newer is newer than older
int isNewer(const char *newer,const char *older){
   if (strstr(older, "ftp://") == older || strstr(older, "http://") == older)
     return 0;
  //  fprintf(stderr,"newer:%s older:%s\n",newer,older);
  // return 0;
  struct stat one;
  struct stat two;
  stat(newer, &one );
  stat(older, &two );
  
  return one.st_mtime>=two.st_mtime;
}


double sum(double *a,size_t b,int doLog){
  assert(b>0);
  double r = 0;
  for(size_t i=0;i<b;i++)
    if(doLog==0)
      r +=a[i];
    else
      r +=log(a[i]);
  return r;
}


double dnorm(double x,double mean,double sd){
    
  double fac = 1.0/(sd*sqrt(2.0*M_PI));
  double val = exp(-(((x-mean)*(x-mean))/(2*sd*sd)));
  
  // if val is 0 because exp(-(x-mean)*(x-mean)) is due to underflow, returns low value
  if(val<LOWER_BOUND){
    return(LOWER_BOUND);
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
  if(p>1-LOWER_BOUND){
    p = 1-LOWER_BOUND;
  } else if(p<LOWER_BOUND){
    p = LOWER_BOUND;
  } 
  return( log(pow(p,k)*pow(1-p,1-k)) );
}

double bernoulli(int k, double p){
  // if p is 0 or 1, cannot do log
  // however this because of over/underlow and p i just very close 0 or 1
  if(p>1-LOWER_BOUND){
    p = 1-LOWER_BOUND;
  } else if(p<LOWER_BOUND){
    p = LOWER_BOUND;
  } 
  return( pow(p,k)*pow(1-p,1-k) );
}
