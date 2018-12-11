#pragma once
#include <cmath>
#include <vector>
#include <zlib.h>
#include <cstdio>
template<typename T>
struct Matrix {
  size_t dx;
  size_t dy;
  size_t mx;
  size_t my;
  T **d;
};

Matrix<double> *initMatrix(size_t x,size_t y);
void kill(Matrix<double> *rr);
double getMax(double a,double b, double c);
double addProtect2(double a,double b);
double addProtect3(double a,double b, double c);
Matrix<double> getMatrix(const char *name);
Matrix<double> getMatrixCheck(const char *name);

std::vector<double> getArray(const char *name);
std::vector<double> getArrayCheck(const char *name);
int fexists(const char* str);
double sigm(double x);

void swapDouble (double& first, double& second);
int matinv( double x[], int n, int m, double space[]);
void deleteMatrix(Matrix<double> mat);
void print(Matrix<double> *mat,FILE *file);
void print(char *ary,size_t l,FILE *file,char*);
void print(double *ary,size_t l,FILE *file,char*);
void logrescale(double *ary,int len);
//void svd_inverse(double mat[],int xLen, int yLen);
int svd_inverse(double mat[],int xLen, int yLen);
std::vector<char*> getFilenames(const char * name,int nInd);
char *strpop(char **str,char split);

template <typename T>
T * allocArray(size_t len,T defval){
  T *ret= new T[len];
  for(size_t i=0;i<len;i++)
    ret[i]=defval;
  return ret;
  
  
}
template <typename T>
  T * allocArray(size_t len){
  T *ret= new T[len];
  return ret;
  
}

template <typename T>
T sum(const T *ary,size_t len){
  T s =0;
  for(size_t i=0;i<len ; i++)
    s+=ary[i];
  //  printf("sum:%f\n",s);
  return s;
}

void print_array(FILE *fp,double *ary,int len);
void print_array(FILE *fp,int *ary,int len);

double *readDouble(const char*fname,int hint);
int whichMax(double *d,int len);


size_t fsize(const char* fname);
int fexists(const char* str);//{///@param str Filename given as a string.
FILE *openFile(const char* a,const char* b);
gzFile openFileGz(const char* a,const char* b,const char *mode);
FILE *getFILE(const char*fname,const char* mode);
gzFile getGz(const char*fname,const char* mode);

int isNewer(const char *newer,const char *older);

double sum(double*,size_t l,int doLog);

double dnorm(double x,double mean,double sd);
double logdnorm(double x, double mean, double sd);
double lbinom(int n, int k);
double logdbinom(int k, int n,double p);
double dbinom(int k, int n,double p);
double logbernoulli(int k, double p);
double bernoulli(int k, double p);
