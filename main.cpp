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
#include "anal.h"

//stupid little function to read 1,3 column of a bim file
std::vector<char*> readBim(char *str){
  
  char *fname = (char*) malloc(strlen(str)+5);
  strcpy(fname,str);
  strcpy(fname+strlen(fname),".bim");
  
  int lens = 4096;
  char *buf = (char*) malloc(lens);

  gzFile gz = Z_NULL;
  if(((gz=gzopen(fname,"rb")))==Z_NULL){
    fprintf(stderr,"Problem opening file: \'%s\'",fname);
    exit(0);
  }
  std::vector<char*> ret;
  while(gzgets(gz,buf,lens)){
    
    char *chr = strtok(buf,"\n \t");
    strtok(NULL,"\t");strtok(NULL,"\n \t");
    char *pos = strtok(NULL,"\n \t");
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
  fprintf(fp, "   -c <filename>       covariance matrix filename\n");
  fprintf(fp, "   -s <filename>       phenotypes\n");
  fprintf(fp, "   -a <filename>       admixproportions (for source pop1)\n");
  fprintf(fp, "   -Q <filename>       .Q file from ADMIXTURE\n");
  fprintf(fp, "   -f <filename>       allele frequencies (.P file)\n");
  fprintf(fp, "   -m <INTEGER>        model 0=add 1=rec\n");
  fprintf(fp, "   -l <INTEGER>        regression 0=linear regression, 1=logistic regression\n");
  fprintf(fp, "   -b <filename>       file containing the start\n");
  fprintf(fp, "   -i <UINTEGER>       maximum iterations\n");
  fprintf(fp, "   -t <FLOAT>          tolerance for breaking EM\n");
  fprintf(fp, "   -r <FLOAT>          seed for rand\n");
  fprintf(fp, "   -P <INT>            number of threads\n");

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
  //new value as in R-code
  int mIter = 40;
  //double tol = 1e-8;
  //new value as in R-code
  double tol = 1e-4;
  int n;
  //int seed = 100;
  //emil - I will give random seed:
  int seed = time(NULL);
  int nThreads = 1;

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
    else{
      fprintf(stderr,"Unknown arg:%s\n",*argv);
      FILE *fp=stderr;
      print_info(fp);
      return 0;
    }
    ++argv;
  }
  

  /* emil I have changed how arguments are read, as this one did not check if arguments were not a char
     meaning it only looked at first letter after "-", so could give "-model" instead of "-m"
  
  while ((n = getopt(argc, argv, "p:o:c:y:a:f:b:m:i:t:r:P:")) >= 0) {
    fprintf(stderr,"arg: %s\n",optarg);
    fprintf(stderr,"arg: %c\n",n);
    
    switch (n) {
    case 'p': pname = strdup(optarg); break; 
    case 'o': outname = strdup(optarg); break;
    case 'c': covname = strdup(optarg); break;
    case 'y': phename = strdup(optarg); break;
    case 'a': adname = strdup(optarg); break;
    case 'f': freqname = strdup(optarg); break;
    case 'm': model = atoi(optarg); break;
    case 'b': startname = strdup(optarg); break;
    case 'i': mIter = atoi(optarg); break;
    case 't': tol = atof(optarg); break;
    case 'r': seed = atoi(optarg); break;
    case 'P': nThreads = atoi(optarg); break;
    default: {fprintf(stderr,"unknown arg: %s\n",optarg);exit(0);;}
      //default: {fprintf(stderr,"unknown arg:\n");return 0;}
    }
  }

  */

#if 1
  if(!outname||!pname||!covname||!phename||(!adname & !qname)||!freqname){
    FILE *fp=stderr;
    fprintf(fp, "\n");
    fprintf(fp, "Usage: line  [options] \n");
    fprintf(fp, "Options:\n");
    fprintf(fp, "   -p \'%s\'\t\tplink prefix filename\n",pname);
    fprintf(fp, "   -o \'%s\'\t\toutput filename\n",outname);
    fprintf(fp, "   -c \'%s\'\t\tcovariance matrix filename\n",covname);
    fprintf(fp, "   -y \'%s\'\t\tphenotypes\n",phename);
    fprintf(fp, "   -a \'%s\'\t\tadmixproportions (for source pop1)\n",adname);
    fprintf(fp, "   -Q \'%s\'\t\t.Q file from ADMIXTURE\n",qname);
    fprintf(fp, "   -f \'%s\'\t\tallele frequencies\n",freqname);
    fprintf(fp, "\n optional arguments:\n");
    fprintf(fp, "   -m \'%d\'\t\tmodel 0=add 1=rec\n",model);
    fprintf(fp, "   -l \'%d\'\t\tregression 0=linear 1=logistic\n",regression);
    fprintf(fp, "   -b \'%s\'\t\tfile containing the start\n",startname);
    fprintf(fp, "   -i \'%d\'\t\tmax number of iterations\n",mIter);
    fprintf(fp, "   -r \'%d\'\t\trandom seed\n",seed);
    fprintf(fp, "   -t \'%e\'\tfloat for breaking EM update\n",tol);
    fprintf(fp, "   -P \'%d\'\t\tnumber of threads\n",nThreads);
    fprintf(fp, "\n");
    fprintf(stderr,"All files must be specified: -p -c -y -a -f -o\n");
    //    print_info(stderr);
    return 1;
  }  
#endif
  

  FILE *outFile = fopen(outname, "w");

  // creates new char array for new name of logfile
  char *logname = new char[strlen(outname)+strlen(".log")+1];

  // copies outname to logname
  strcpy(logname,outname);
  // append .log at logname char array
  strncat(logname,".log",strlen(".log"));

  FILE *logFile = fopen(logname, "w");
  delete [] logname;

  //srand48(seed);
  std::srand(seed);
  fprintf(stderr,"Seed is: %i\n",seed);
  
  plink *p = readplink(pname);
  
  std::vector<char *> loci = readBim(pname);

  fprintf(stderr,"Seed is: %s\n",pname);
  fprintf(stderr,"Seed is: %s\n",covname);  
  Matrix<double> cov = getMatrix(covname);

  fprintf(stderr,"Seed is: %s\n",covname);
  if(0){
    print(&cov,stdout);
    return 0;
  }
  std::vector<double> pheno = getArray(phename);
  std::vector<double> adprop(pheno.size());
  if(qname!=NULL){
    Matrix<double> Q = getMatrix(qname);
    // I just read the first column of the .Q file into the adprop vector
    for(int i=0;i<Q.dx;i++){
      adprop[i]=Q.d[i][0];      
    }
  }  else if(adname!=NULL){
    // make so can give .Q file instead of just 1 column file...
    adprop = getArray(adname);
  } else{
    fprintf(stderr,"ancprobs file or .Q file MUST be provided!\n");
    exit(1);
  }

  Matrix<double> f = getMatrix(freqname);
  std::vector<double> s = getArray(startname);
  
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

  fprintf(logFile,"Command had following options :\t -p %s -o %s ",pname,outname);
  fprintf(logFile,"-c %s -y %s -Q %s -a %s -f %s -i %i -t %f -m %i -l %i -r %i ",covname,phename,qname,adname,freqname,mIter,tol,model,regression,seed);
  fprintf(logFile,"-P %i -b %s\n",nThreads,startname);
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
  
  wrap(p,pheno,adprop,f,model,s,cov,mIter,tol,loci,nThreads,outFile,logFile,regression);

  //cleanup
  kill_plink(p);

  /*
  // emil changed reading in of arguments
  free(pname);
  free(outname);
  free(covname);
  free(phename);
  free(adname);  
  free(freqname);
  free(startname);
  */

  fclose(logFile);
  fclose(outFile);
  
  for(int i=0;i<cov.dx;i++){
    delete [] cov.d[i];
  }
  delete [] cov.d;
  
  for(int i=0;i<f.dx;i++){
    delete [] f.d[i];
  }
  delete [] f.d;
 
  for(uint i=0;i<loci.size();i++){
    free(loci[i]);
  }
  return 0;
}
