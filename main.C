// file: main.C for fsabc

// Runs simulations for ABC inference of selection on a polygenic trait. Assumes a model of environment-dependent directional selection on the trait, with evolution of genetic loci dictated by a genetic architecture. Drift if incorporated via a WF model.

// This program has three primary modes, one to conduct simulations for inference = "sim" mode, one to conduct posterior predictive simulations = "pred" mode, and one to calculate summary statistics from observed data = "obs" mode

// "sim" mode requires the following files
// geneFile=point estimates of allele freqs. nloci rows by npop columns, just the initial freqs (as freqs)
// envFile=environmental covariates, one row per generation, one column per pop, note that the gen t (the final gen.) should be left off
// neFile=post. samples for the variance effective population size, one row per sample, one column per pop
// traitFile=trait gen arch. estimates, one row per SNP, pip followed by beta | lambda = 1
// Additional optional files
// missingFile=binary indicator variable, was that population x generation sampled (1=TRUE); one row per generation, one column per pop, note that first generation should be left off 

// NOTE: files begin with a row giving the dimensions of the file (row then column)

// "obs" mode requires the following files

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <getopt.h>
#include <omp.h>
#include "abcfs.H"

using namespace std;

gsl_rng * r;  /* global state variable for random number generator */

/* ----------------- */
/* Beginning of main */
/* ----------------- */

int main(int argc, char *argv[]) {

  time_t start = time(NULL);
  time_t end;
  int rng_seed = 0;
  int ch = 0;
  int x;
  int nsims = 1000;
  int ss = 0; // 0 = bv ss, 1 = snp ss, 2 = bv and snp ss
  int validate = 0; // set to 1 for cross validation
  int obsss = 0; // set to 1 for obs summary stats mode

  string geneFile = "undefined";
  string envFile = "undefined";
  string neFile = "undefined";
  string traitFile = "undefined";
  string missingFile = "undefined";
  string resFile = "undefined"; // samples from posterior
  string outFile = "out_fsabc.txt";
  
  dataset data;
  param params;

  // set defaults
  data.sigma2 = 1;
  params.smod = 0;
  params.fsAlb = -10;
  params.fsAub = 10;
  params.fsBlb = -10;
  params.fsBub = 10;
  params.fsClb = -1;
  params.fsCub = 1;
  params.prsel = 0.5;

  
  // get command line arguments
  if (argc < 2) {
    usage(argv[0]);
  }
  
  while ((ch = getopt(argc, argv, "g:e:f:t:j:o:n:s:m:a:c:b:d:w:x:p:v:z:q:")) != -1){
    switch(ch){
    case 'g':
      geneFile = optarg;
      break;
    case 'e':
      envFile = optarg;
      break;
    case 'f':
      neFile = optarg;
      break;
    case 't':
      traitFile = optarg;
      break;
    case 'j':
      missingFile = optarg;
      break;
    case 'o':
      outFile = optarg;
      break;
    case 'n':
      nsims = atoi(optarg);
      break;
    case 's':
      ss = atoi(optarg);
      break;
    case 'm':
      params.smod = atoi(optarg);
      break;
    case 'a':
      params.fsAlb = atof(optarg);
      break;
    case 'c':
      params.fsAub = atof(optarg);
      break;
    case 'b':
      params.fsBlb = atof(optarg);
      break;
    case 'd':
      params.fsBub = atof(optarg);
      break;
    case 'w':
      params.fsClb = atof(optarg);
      break;
    case 'x':
      params.fsCub = atof(optarg);
      break;
    case 'p':
      params.prsel = atof(optarg);
      break;
    case 'v':
      validate = atoi(optarg);
      break;
    case 'q':
      obsss = atoi(optarg);
      break;
    case 'z':
      resFile = optarg;
      break;
    case '?':
    default:
      usage(argv[0]);
    }
  }
  
  // set up gsl random number generation 
  gsl_rng_env_setup();
  r = gsl_rng_alloc (gsl_rng_default);
  srand(time(NULL));
  rng_seed = rand();
  gsl_rng_set(r, rng_seed); /* seed gsl_rng with output of rand, which
                               was seeded with result of time(NULL) */

  // read infiles, record genotype likelihoods and data dimensions
  cout << "Reading input from files: " << geneFile << " and " << 
    envFile << endl;
  getdata(geneFile, envFile, missingFile, obsss, &data);

  // read Ne from a file
  cout << "Reading posterior samples of Ne from " << neFile << endl;
  getne(neFile,&data);

  // read QTL from a file
  cout << "Reading GWA data " << traitFile << endl;
  getqtl(traitFile,&data);

  // read params if present
  if(validate==1){
    cout << "Reading parameter estimates for validation " << resFile << endl;
    getest(resFile,&data);
  }
  
  // memory allocation for params
  params.beta = gsl_vector_calloc(data.nLoci);
  params.pp = gsl_matrix_calloc(data.nLoci, data.nPops * data.nGens);
  params.dpp =  gsl_matrix_calloc(data.nLoci, data.nPops * (data.nGens-1));
  params.ssSumdp = gsl_vector_calloc(data.nLoci);
  params.ssEnvcov = gsl_vector_calloc(data.nLoci);
  params.Svec = gsl_vector_calloc(data.nPops * (data.nGens-1));
  
  // open outfile
  FILE * OUT;
  OUT = fopen(outFile.c_str(), "w");

  // "obs" mode
  if(obsss==1){
    cout << "Calculating summary statistics for the observed time series" << endl;
    calcobs(&data, &params, OUT, ss);
  }

  // "pred" mode
  else if(validate==1){
    cout << "Running posterior predictive simulations" << endl;
    // run wf sims
    for(x=0; x<nsims; x++){
      if(((x%100)==0) & (x>0))
	cout << "Sim. no.: " << x << endl;
      valsim(&data, &params, OUT);
      gsl_vector_free(params.sivec); // need to free each time b/c allocated by samqtl
      // and depends on the # of qtl sampled
    }
  }

  // "sim" mode
  else{
    cout << "Running simulations for ABC inference" << endl;
    // run wf sims
    for(x=0; x<nsims; x++){
      if(((x%100)==0) & (x>0))
	cout << "Sim. no.: " << x << endl;
      runsim(&data, &params, x, ss);
      // write sims
      writesim(&data, &params, OUT, ss);
      gsl_vector_free(params.sivec); // need to free each time b/c allocated by samqtl
      // and depends on the # of qtl sampled
    }
  }

  // close outfile
  fclose(OUT);

  // prints run time
  end = time(NULL);
  cout << "Runtime: " << (end-start)/3600 << " hr " << (end-start)%3600/60 << " min ";
  cout << (end-start)%60 << " sec" << endl;
  return 0;
}
