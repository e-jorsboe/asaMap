/*
Lines project:
./line  -p plinkTest -o out -c plinkTest.covs -y plinkTest.y -a plinkTest.Q1 -f plinkTest.f1 
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <ctime>
#include "readplink.h"
#include "analysisFunction.h"
#include "analysis.h"
#include <iostream>
#include <sys/sysinfo.h>
#include <stdio.h>

#ifdef EIGEN
#include <Eigen/Dense>
#endif

//stupid little function to read 1,3 column of a bim file
std::vector<char*> readBim(char *str){
  
  char *fname = (char*) malloc(strlen(str)+5);
  strcpy(fname,str);
  strcpy(fname+strlen(fname),".bim");
  
  int lens = 4096;
  char *buf = (char*) malloc(lens);
  char *save;
  gzFile gz = Z_NULL;
  if(((gz=gzopen(fname,"rb")))==Z_NULL){
    fprintf(stderr,"Problem opening file: \'%s\'",fname);
    exit(0);
  }
  std::vector<char*> ret;
  while(gzgets(gz,buf,lens)){
    
    char *chr = strtok_r(buf,"\n \t",&save);
    strtok_r(NULL,"\t",&save);strtok_r(NULL,"\n \t",&save);
    char *pos = strtok_r(NULL,"\n \t",&save);
    char newItem[1024];
    snprintf(newItem,124,"%s\t%s\t",chr,pos);
    ret.push_back(strdup(newItem));
  }

  free(fname);
  free(buf);
  gzclose(gz);
  return ret;
}

void print_info(FILE *fp){
  fprintf(fp, "\n");
  fprintf(fp, "Usage: line  [options] \n");
  fprintf(fp, "Options:\n");
  fprintf(fp, "   -p <filename>       plink prefix filename\n");
  fprintf(fp, "   -o <filename>       output filename\n");
  fprintf(fp, "   -c <filename>       covariates filename\n");
  fprintf(fp, "   -y <filename>       phenotypes\n");
  fprintf(fp, "   -a <filename>       admixproportions (for source pop1) - either specify this or -Q\n");
  fprintf(fp, "   -Q <filename>       .Q file from ADMIXTURE - either specify this or -a\n");
  fprintf(fp, "   -f <filename>       allele frequencies (.P file)\n");
  fprintf(fp, "   -m <INT>            model 0=add 1=rec (default: 0)\n");
  fprintf(fp, "   -l <INT>            regression 0=linear regression, 1=logistic regression (default: 0)\n");
  fprintf(fp, "   -b <filename>       file containing the start\n");
  fprintf(fp, "   -i <INT>            maximum iterations (default: 80)\n");
  fprintf(fp, "   -t <FLOAT>          tolerance for breaking EM (default: 0.0001)\n");
  fprintf(fp, "   -r <INT>            seed for rand\n");
  fprintf(fp, "   -P <INT>            number of threads\n");
  fprintf(fp, "   -e <INT>            estimate standard error of coefficients (0: no (default), 1: yes)\n");
  fprintf(fp, "   -w <INT>            run M0/R0 model that models effect of other allele (0: no, 1: yes (default))\n");
  fprintf(fp, "   -C <INT>            prime coefficiens of models with estimates from previous models - e.g. using M0 coefs as starting point in M1 (0: no, 1: yes (default))\n");  
  fprintf(fp, "\n");
}

int main(int argc,char **argv){
  
  if(argc==1){// if no arguments, print info on program
    FILE *fp=stderr;
    print_info(fp);
    return 0;
}
  
  char *pname = NULL;
  char *outname = NULL;
  char *covname = NULL;
  char *phename = NULL;
  char *adname = NULL;
  char *freqname = NULL;
  char *qname = NULL;
  char *startname = NULL;
  int model = 0;
  int regression = 0;
  //int mIter = 10;
  //new value 40 in R-code
  // however M0 doest not always converge with this
  int mIter = 80;
  //double tol = 1e-8;
  //new value as in R-code
  double tol = 1e-4;
  int n;
  //int seed = 100;
  //emil - I will give random seed:
  int seed = time(NULL);
  int nThreads = 1;
  int estSE = 0;
  int useM0R0 = 1;
  int primeCoefs = 1;
  
  argv++;
  while(*argv){
    // reading in arguments
    if(strcmp(*argv,"-p")==0) pname=*++argv; 
    // outfile
    else if(strcmp(*argv,"-o")==0) outname=*++argv;
    // covariate file
    else if(strcmp(*argv,"-c")==0) covname=*++argv; 
    // pheno file
    else if(strcmp(*argv,"-y")==0) phename=*++argv;
    // ancestry file
    else if(strcmp(*argv,"-a")==0) adname=*++argv;
    //.Q file
    else if(strcmp(*argv,"-Q")==0) qname=*++argv;
    // pop specific freqs
    else if(strcmp(*argv,"-f")==0) freqname=*++argv;
    // which model to use (0: add, 1: rec)
    else if(strcmp(*argv,"-m")==0) model=atoi(*++argv); 
    // linear/logistic regression
    else if(strcmp(*argv,"-l")==0) regression=atoi(*++argv);
    // file with starting points
    else if(strcmp(*argv,"-b")==0) startname=*++argv; 
    // number of iterations to run for the EM algorithm
    else if(strcmp(*argv,"-i")==0) mIter=atoi(*++argv); 
    // tolerance for deciding when EM algorithm has converged
    else if(strcmp(*argv,"-t")==0) tol=atof(*++argv);
    // random seed
    else if(strcmp(*argv,"-r")==0) seed=atoi(*++argv); 
    // number of threads
    else if(strcmp(*argv,"-P")==0) nThreads=atoi(*++argv);
    else if(strcmp(*argv,"-e")==0) estSE=atoi(*++argv);
    // run M0/R0 model or not
    else if(strcmp(*argv,"-w")==0) useM0R0=atoi(*++argv);
    // prime coefs
    else if(strcmp(*argv,"-C")==0) primeCoefs=atoi(*++argv);
    else{
      fprintf(stderr,"Unknown arg:%s\n",*argv);
      FILE *fp=stderr;
      print_info(fp);
      return 0;
    }
    ++argv;
  }
  
  clock_t t = clock();
  time_t t2 = time(NULL);
 
  // creates new char array for new name of logfile
  char *logname = new char[strlen(outname)+strlen(".log")+1];

  // copies outname to logname
  strcpy(logname,outname);
  // append .log at logname char array
  strncat(logname,".log",strlen(".log"));

  FILE *logFile = fopen(logname, "w");
  delete [] logname;

  char *outname2 = new char[strlen(outname)+strlen(".res")+1];

  // copies outname to logname
  strcpy(outname2,outname);

  //puts .res at end of res file
  strncat(outname2,".res",strlen(".res"));
  FILE *outFile = fopen(outname2, "w");
  
  //srand48(seed);
  std::srand(seed);
  fprintf(stderr,"Seed is: %i\n",seed);
  
  plink *p = readplink(pname);  
  std::vector<char *> loci = readBim(pname);

  std::vector<double> pheno = getArrayCheck(phename);
  
  Matrix<double> cov;
  
  if(covname!=NULL){
    cov = getMatrixCheck(covname);
  } else{    
    cov.mx=pheno.size();
    cov.my=0;
    cov.dx=pheno.size();
    cov.dy=0;
    cov.d=NULL;
  }
  
  std::vector<double> adprop(pheno.size());
  if(qname!=NULL){
    Matrix<double> Q = getMatrixCheck(qname);
    // I just read the first column of the .Q file into the adprop vector
    for(int i=0;i<Q.dx;i++){
      adprop[i]=Q.d[i][0];      
    }
    
    for(int i=0;i<Q.dx;i++){
      delete [] Q.d[i];
    }
    delete [] Q.d;        
    
  }  else if(adname!=NULL){
    // make so can give .Q file instead of just 1 column file...
    adprop = getArrayCheck(adname);
  } else{
    fprintf(stderr,"ancprobs file or .Q file MUST be provided!\n");
    exit(1);
  }

  Matrix<double> f = getMatrixCheck(freqname);
  std::vector<double> s = getArrayCheck(startname);
  
  // check files match in dimensions!!!!!!!
  assert(cov.mx==pheno.size());
  assert(adprop.size()==pheno.size());
  assert(p->x==pheno.size());
  // .P file has to have all loci of plink data
  assert(f.mx==loci.size());
  // only implemneted for 2 pops
  assert(f.my==2);
  // additive or recessive model
  assert(model==0 || model==1);
  // linear or logisitc regression
  assert(regression==0 || regression==1);
  // no standard error or standard error
  assert(estSE==0 || estSE==1);
  // run M0/R0 or not
  assert(useM0R0==0 || useM0R0==1);
  // prime coefs or not
  assert(primeCoefs==0 || primeCoefs==1);
  // positive number of threads
  assert(nThreads > 0);
  // positive number of EM iterations
  assert(mIter > 0);
  // positive number for seed
  assert(seed > 0);

  
  if(nThreads > get_nprocs()){
    nThreads = get_nprocs();
    fprintf(stderr,"Requested more threads than avaible cores - nThreads has been set to %i\n",nThreads);
    fprintf(logFile,"Requested more threads than avaible cores - nThreads has been set to %i\n",nThreads);    
  }
  
  
  fprintf(logFile,"Command had following options :\t -p %s -o %s ",pname,outname);
  fprintf(logFile,"-c %s -y %s -Q %s -a %s -f %s -i %i -t %f -m %i -l %i -r %i ",covname,phename,qname,adname,freqname,mIter,tol,model,regression,seed);
  fprintf(logFile,"-P %i -b %s -e %i -w %i -C %i\n",nThreads,startname,estSE,useM0R0,primeCoefs);
  fprintf(logFile,"\n");
  fprintf(logFile,"Seed is: %i\n",seed);
  fprintf(logFile,"\n");
  fprintf(logFile,"Done reading file: '%s' with dim ncol: %zu\tnrow: %zu\n",pname,p->y,p->x);
  fprintf(logFile,"Done reading file: '%s' containing nrows: %zu\tncols: %zu\n",covname,cov.mx,cov.my);
  fprintf(logFile,"Done reading file: '%s' containing nitems: %zu\n",phename,pheno.size());
  if(adname!=NULL){
    fprintf(logFile,"Done reading file: '%s' containing nitems: %zu\n",adname,adprop.size());
  } else{
    fprintf(logFile,"Done reading file: '%s' containing nitems: %zu\n",qname,adprop.size());
  }
  fprintf(logFile,"Done reading file: '%s' containing nrows: %zu\tncols: %zu\n",freqname,f.mx,f.my);

#ifdef EIGEN
#pragma message "Compiling with EIGEN" 
#endif

  // flush to disk - or force to write to disk 
  fflush(logFile);
  
  wrap(p,pheno,adprop,f,model,s,cov,mIter,tol,loci,nThreads,outFile,logFile,regression,estSE,useM0R0,outname2,primeCoefs);
  
  //cleanup
  kill_plink(p);

  delete [] outname2;
  fclose(outFile);

  if(covname!=NULL){  
    for(int i=0;i<cov.dx;i++){
      delete [] cov.d[i];
    }
    delete [] cov.d;
  }
  
  for(int i=0;i<f.dx;i++){    
    delete [] f.d[i];
  }
  delete [] f.d;
 
  for(uint i=0;i<loci.size();i++){
    free(loci[i]);
  }   
  
  fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));

  // print to log file
  fprintf(logFile, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(logFile, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));
  fclose(logFile); 
  
  return 0;
}
