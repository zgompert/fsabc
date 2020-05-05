#include <iostream>
#include <sstream>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_sort.h>
#include <float.h>
#include <cmath>

#include "abcfs.H"

using namespace std;

// print software usage
void usage(char * name){
  fprintf(stdout,"\n%s version %s\n\n", name, VERSION);
  fprintf(stdout, "Usage: fsabc -g genefile -e envfile -f nefile -t traitfile [options]\n");
  fprintf(stdout, "Use -v 1 for predictive mode, -q 1 for observed mode,\n or leave both 0 (default) for simulation model\n\n");
  fprintf(stdout, "-g     Infile with allele frequency data\n");
  fprintf(stdout, "-e     Infile with environmental covariate data\n");
  fprintf(stdout, "-f     Infile with varNe estimates\n");
  fprintf(stdout, "-t     Infile with trait genetic arch. estimates\n");
  fprintf(stdout, "-j     (optional) Infile with gens. between samples\n");
  fprintf(stdout, "-z     (optional) Infile parameter posterior samples\n");
  fprintf(stdout, "-v     Binary, run posterior pred. validation mode [0]\n");
  fprintf(stdout, "-q     Binary, run observed summary stats. mode [0]\n");
  fprintf(stdout, "-o     Outfile for simulated or obs. summary stats. [out_fsabc.txt]\n");
  fprintf(stdout, "-n     Number of simulations [1000]\n");
  fprintf(stdout, "-s     SS to print: 0 = bv, 1 = snp, 2 = both [0]\n");
  fprintf(stdout, "-m     Selection model: 0 = linear, 1 = step, 2 = sigmoid [0]\n");
  fprintf(stdout, "-p     Prior prob. of non-zero selection by component [0.5]\n");
  fprintf(stdout, "-a     Lower bnd. on U prior for sel. function intercept [-10]\n");
  fprintf(stdout, "-c     Upper bnd. on U prior for sel. function intercept [10]\n");
  fprintf(stdout, "-b     Lower bnd. on U prior for sel. function slope [-10]\n");
  fprintf(stdout, "-d     Upper bnd. on U prior for sel. function slope [10]\n");
  fprintf(stdout, "-w     Lower bnd. on U prior for sel. function cut [-1]\n");
  fprintf(stdout, "-x     Upper bnd. on U prior for sel. function cut [1]\n");
  

  exit(1);
}

// ------ Functions for input and output ---------------

// read input from the infiles
void getdata(string geneFile, string envFile, string missingFile, int obsss, dataset * data){
  int i, j;
  string line, element;
  ifstream infile;
  istringstream stream;

  // read environmental covariate, single column of values
  infile.open(envFile.c_str());
  if (!infile){
    cerr << "Cannot open file " << envFile << endl;
    exit(1);
  }

  getline(infile, line);
  stream.str(line);
  stream.clear();
  stream >> element; // number of gens
  data->nGens = atoi(element.c_str()); 
  stream >> element; // number of populations
  data->nPops = atoi(element.c_str()); 
  
  // dynamic memory allocation of env data
  data->envData = gsl_matrix_calloc(data->nGens-1, data->nPops);
  
  // loop through file, one row per gen (-1) one column per pop
  for(i=0; i<(data->nGens-1); i++){
    stream.clear();
    getline(infile, line);
    stream.str(line);
    stream.clear();
    for(j=0; j<data->nPops; j++){
      stream >> element;
      gsl_matrix_set(data->envData, i, j, atof(element.c_str()));
    }
  }
  infile.close();

  
  // read initial allele freq. data (full set of af data for obs. ss)
  infile.open(geneFile.c_str());
  if (!infile){
    cerr << "Cannot open file " << geneFile << endl;
    exit(1);
  }

  // read line with data dimensions
  getline(infile, line);
  stream.str(line);
  stream.clear();
  stream >> element; // number of loci
  data->nLoci = atoi(element.c_str()); 
  stream >> element; // number of populations or pops. x gens for observed


  if(obsss==1){ // for obs. ss calculations
    // dynamic memory allocation for allele frequency data
    data->pall = gsl_matrix_calloc(data->nLoci, data->nPops * data->nGens);
    
    // read and store allele frequencies
    for(i=0; i<data->nLoci; i++){
      getline(infile, line); // data for one locus
      stream.str(line);
      stream.clear();
      for(j=0; j<(data->nPops*data->nGens); j++){
	stream >> element;
	gsl_matrix_set(data->pall,i, j, atof(element.c_str()));
      }
    }
    infile.close();
  }

  else{ // for all simulation-based analyses
    // dynamic memory allocation for allele frequency data
    data->p0 = gsl_matrix_calloc(data->nLoci, data->nPops);
    
    // read and store allele frequencies
    for(i=0; i<data->nLoci; i++){
      getline(infile, line); // data for one locus
      stream.str(line);
      stream.clear();
      for(j=0; j<data->nPops; j++){
	stream >> element;
	gsl_matrix_set(data->p0,i, j, atof(element.c_str()));
      }
    }
    infile.close();
  }
      
  // dynamic memory allocation for number of gens between samples
  data->sampledData = gsl_matrix_int_calloc(data->nGens-1, data->nPops);
  gsl_matrix_int_set_all(data->sampledData, 1); // defaults to one generation spacing
  infile.open(missingFile.c_str());
  if (infile){ // optional file was included
    cout << "Reading optional file with generation spacing" << endl;
    getline(infile, line);
    stream.str(line);
    stream.clear();
    stream >> element; // number of gens - 1
    stream >> element; // number of populations

    // loop through file, one row per gen (-1) one column per pop
    for(i=0; i<(data->nGens-1); i++){
      stream.clear();
      getline(infile, line);
      stream.str(line);
      stream.clear();
      for(j=0; j<data->nPops; j++){
	stream >> element;
	gsl_matrix_int_set(data->sampledData, i, j, atoi(element.c_str()));
      }
    }
    infile.close();
  }
}

// read existing estimates of Ne
void getne(string neFile, dataset * data){
  int i, j;
  string line, element;
  ifstream infile;
  istringstream stream;

  // read ne data, first line has dimensions, combes by samples, followed by samples of Ne
  infile.open(neFile.c_str());
  if (!infile){
    cerr << "Cannot open file " << neFile << endl;
    exit(1);
  }

  // read line with data dimensions
  getline(infile, line);
  stream.str(line);
  stream.clear();
  stream >> element; // number of posterior samples
  data->nNeSams = atoi(element.c_str());
  stream >> element; // number of pops 
 

  // dynamic memory allocation for existing Ne estimates
  data->ne = gsl_matrix_calloc(data->nNeSams, data->nPops);

  // read and store Ne estimates
  for(i=0; i<data->nNeSams; i++){
    getline(infile, line); 
    stream.str(line);
    stream.clear();
    for(j=0; j<data->nPops; j++){
      stream >> element;
      gsl_matrix_set(data->ne, i, j, atof(element.c_str()));
    }
  }
  infile.close();

}

 // read existing estimates of genotype x phenotype associations
void getqtl(string traitFile, dataset * data){
  int i, j;
  string line, element;
  ifstream infile;
  istringstream stream;

  // read qtl data, first line has dimensions, followed by one row per locus with pip and beta
  infile.open(traitFile.c_str());
  if (!infile){
    cerr << "Cannot open file " << traitFile << endl;
    exit(1);
  }

  // read line with data dimensions
  getline(infile, line);
  stream.str(line);
  stream.clear();
  stream >> element; // number of loci
  stream >> element; // 2, pip, beta
 

  // dynamic memory allocation for QTL effects
  data->qtl = gsl_matrix_calloc(data->nLoci, 2);

  // read and store QTL estimates
  for(i=0; i<data->nLoci; i++){
    getline(infile, line); // data for one locus
    stream.str(line);
    stream.clear();
    for(j=0; j<2; j++){
      stream >> element;
      gsl_matrix_set(data->qtl, i, j, atof(element.c_str()));
    }
  }
  infile.close();

}

// // read existing estimates of Ne
// void getne(string neFile, dataset * data){
//   int i, j;
//   string line, element;
//   ifstream infile;
//   istringstream stream;

//   // read ne data, first line has dimensions, combes by samples, followed by samples of Ne
//   infile.open(neFile.c_str());
//   if (!infile){
//     cerr << "Cannot open file " << neFile << endl;
//     exit(1);
//   }

//   // read line with data dimensions
//   getline(infile, line);
//   stream.str(line);
//   stream.clear();
//   stream >> element; // number of posterior samples
//   data->nNeSams = atoi(element.c_str());
//   stream >> element; // number of pops 
 

//   // dynamic memory allocation for existing Ne estimates
//   data->ne = gsl_matrix_calloc(data->nNeSams, data->nPops);

//   // read and store Ne estimates
//   for(i=0; i<data->nNeSams; i++){
//     getline(infile, line); 
//     stream.str(line);
//     stream.clear();
//     for(j=0; j<data->nPops; j++){
//       stream >> element;
//       gsl_matrix_set(data->ne, i, j, atof(element.c_str()));
//     }
//   }
//   infile.close();

// }

 // read in parameter estimates for posterior predictive check
void getest(string resFile, dataset * data){
  int i, j;
  int Np, Ns;
  string line, element;
  ifstream infile;
  istringstream stream;

  // read parameter estimates, first line has dimensions, samples by params, followed by samples of parameters
  infile.open(resFile.c_str());
  if (!infile){
    cerr << "Cannot open file " << resFile << endl;
    exit(1);
  }

  // read line with data dimensions
  getline(infile, line);
  stream.str(line);
  stream.clear();
  stream >> element; // number of samples
  Ns = atoi(element.c_str());
  stream >> element; // number of params
  Np = atoi(element.c_str());

  // dynamic memory allocation for parameters
  data->params = gsl_matrix_calloc(Ns, Np);

  // read and store estimates
  for(i=0; i<Ns; i++){
    getline(infile, line); // data for one locus
    stream.str(line);
    stream.clear();
    for(j=0; j<Np; j++){
      stream >> element;
      gsl_matrix_set(data->params, i, j, atof(element.c_str()));
    }
  }
  infile.close();

}


// wf simulation wrapper function
void runsim(dataset * data, param * params, int x, int ss){

  int i, j;
  
  // sample params from priors
  // QTL effects
  samqtl(data, params);

  // selection
  // uses spike and slab prior
  // average effect
  if(gsl_ran_flat(r,0,1) < params->prsel)
    params->fsA = gsl_ran_flat(r, params->fsAlb, params->fsAub);
  else
    params->fsA = 0;
  // environmental effect
  if(gsl_ran_flat(r,0,1) < params->prsel)
    params->fsB = gsl_ran_flat(r, params->fsBlb, params->fsBub);
  else
    params->fsB = 0;

  if(params->smod > 0){ // sigmoid or step function must have slope
    params->fsC = gsl_ran_flat(r, params->fsClb, params->fsCub);
    if((params->fsA == 0) & (params->fsB == 0))
	params->fsC = 0;
  }
  
  // sign might be specified by A - B for sigmoid (i.e., needs to go in the right direction)
    
  // conduct the simulations
  for(j=0; j<data->nPops; j++){
    simpop(data, params, j);		   
  }
  // compute derived parameters
  compdpar(data, params);

  // calculate and store the summary statistics
  for(i=0; i<data->nLoci; i++){
    if(ss > 0)
      calcss(data, params, i);
  }
    if((ss == 0) | (ss == 2))
      calcssbv(data, params);
  
}

// sample non-zero QTL effects based on their PIPs
void samqtl(dataset * data, param * params){
  int i;
  double pr;

  params->Nqtl = 0;
  gsl_vector_set_zero(params->beta);
  for(i=0; i<data->nLoci; i++){
    pr = gsl_rng_uniform(r);
    if(pr < gsl_matrix_get(data->qtl, i, 0)){ // has non-zero effect
      gsl_vector_set(params->beta, i, gsl_matrix_get(data->qtl, i, 1));
      params->Nqtl++;
    }
    // uncomment to set to model-averaged effect
    //gsl_vector_set(params->beta, i, gsl_matrix_get(data->qtl, i, 1) * gsl_matrix_get(data->qtl, i, 0));
  }
  params->sivec = gsl_vector_calloc(params->Nqtl * data->nPops * (data->nGens-1));
}

// simulate evolution for one population given a selection differential
void simpop(dataset * data, param * params, int j){
  int i, k, x, ii;
  int interval;
  double S, si;
  double p, dp, pprime;
  unsigned int Ne;
  int N;

  // sample posterior value of Ne
  x = gsl_ran_flat(r, 0, data->nNeSams);
  Ne = floor(gsl_matrix_get(data->ne, x, j));
  //cout << Ne << endl;

  // set initial allele freqs. to p0
  for(i=0; i<data->nLoci; i++){
    p = gsl_matrix_get(data->p0, i, j);
    gsl_matrix_set(params->pp, i, j * data->nGens, p); 
  }
  
  // sims
  for(k=0; k<(data->nGens-1); k++){

    // calculate the selection differential
    if(params->smod == 0) // linear
      S = params->fsA + params->fsB * gsl_matrix_get(data->envData, k, j);
    
    else if(params->smod == 1){ //step
      if( gsl_matrix_get(data->envData, k, j) < params->fsC)
	S = params->fsA;
      else
	S = params->fsA + params->fsB;
    }
    
    else
      S = params->fsA + (params->fsB - params->fsA)/
	(exp(-1 * params->fsC *  gsl_matrix_get(data->envData, k, j)));

 
    gsl_vector_set(params->Svec, j * (data->nGens-1) + k, S);

    interval = gsl_matrix_int_get(data->sampledData, k, j); // spacing
    for(i=0; i<data->nLoci; i++){
      p = gsl_matrix_get(params->pp, i, j * data->nGens + k);
      for(ii=0; ii<interval; ii++){ // repeate for number of intervals between samples
	N = 0;
	// calculate/simulate dp by selection and drift for each SNP
	
	if( gsl_vector_get(params->beta, i) != 0){ // SNP has effect on trait
	
	  si = gsl_vector_get(params->beta, i) * (S/data->sigma2); // from L&W 5.21
	  // note that this is an approximation; this is a good approximation when z ~ N
	  // it is a first order approximation, won't work if selection is on variance rather
	  // than the mean of the trait
	  gsl_vector_set(params->sivec, N + k * params->Nqtl +
			 j * (data->nGens - 1) * params->Nqtl,
			 abs(si));
	  dp =  si * p; // from L&W 5.8b
	  N++;
	}
	else{
	  dp = 0;
	}

	pprime = dp + p;
	if(pprime > 1)
	  pprime = 1;
	if(pprime < 0)
	  pprime = 0;
      
	// sample new p_i
	p = gsl_ran_binomial(r, pprime, 2 * Ne)/(2.0 * Ne);
	if(ii == (interval-1)){ // final pass, store it
	  gsl_matrix_set(params->pp, i, j * data->nGens + k + 1, p); // store allele freq.
	  // obs dp
	  dp = p - gsl_matrix_get(params->pp, i, j * data->nGens + k);
	  gsl_matrix_set(params->dpp, i, j * (data->nGens - 1) + k, dp); // store allele freq. change
	}
      }
    }
  }
}

// Compute derived paramers describing mean and variance of S and si
void compdpar(dataset * data, param * params){

  // mean and sd for S
  params->Smn = gsl_stats_mean(params->Svec->data, params->Svec->stride, params->Svec->size);
  params->Ssd = gsl_stats_sd_m(params->Svec->data, params->Svec->stride, params->Svec->size, params->Smn);

  // mean and sd for si
  params->simn = gsl_stats_mean(params->sivec->data, params->sivec->stride, params->sivec->size);
  params->sisd = gsl_stats_sd_m(params->sivec->data, params->sivec->stride, params->sivec->size, params->simn);
}

// calculate summary statistics, sum dp and cov dp env.
void calcss(dataset * data, param * params, int i){

  int j, k;
  double sumdp = 0;
  double val, envcov = 0;
  double cov = 0;
  
  gsl_vector * sumdpv;
  gsl_vector * envv;
  
  sumdpv = gsl_vector_calloc(data->nPops * (data->nGens - 1));
  envv = gsl_vector_calloc(data->nPops * (data->nGens - 1));
  
  for(j=0; j<data->nPops; j++){
    for(k=0; k<(data->nGens-1); k++){
      sumdp += gsl_matrix_get(params->dpp, i, j * (data->nGens - 1) + k);
      val = gsl_matrix_get(data->envData, k, j);
      gsl_vector_set(envv, j * (data->nGens - 1) + k, val);
      envcov += val;
    }
  }
  gsl_vector_set(params->ssSumdp, i, sumdp);

   
  //compute means
  sumdp = sumdp / (double) (data->nPops * (data->nGens - 1));
  envcov = envcov / (double) (data->nPops * (data->nGens - 1));

  // compute covariance
  for(j=0; j<data->nPops; j++){
    for(k=0; k<(data->nGens-1); k++){
      cov += (gsl_matrix_get(params->dpp, i, j * (data->nGens - 1) + k) - sumdp) *
	(gsl_vector_get(envv, j * (data->nGens - 1) + k) - envcov);
    }
  }
      
  cov = cov / (double) (data->nPops * (data->nGens - 1));
    
  gsl_vector_set(params->ssEnvcov, i, cov);
  
  gsl_vector_free(sumdpv);
  gsl_vector_free(envv);
}

// calculate breeding value level summary statistics
void calcssbv(dataset * data, param * params){

  int j, k, i;
  double sumdpA, mudp = 0;
  double val, envcov = 0;
  double cov = 0;
 
  
  gsl_vector * mudpv;
  gsl_vector * envv;
  
  mudpv = gsl_vector_calloc(data->nPops * (data->nGens - 1));
  envv = gsl_vector_calloc(data->nPops * (data->nGens - 1));
  
  for(j=0; j<data->nPops; j++){
    for(k=0; k<(data->nGens-1); k++){
      sumdpA = 0;
      for(i=0; i<data->nLoci; i++){
	sumdpA += gsl_vector_get(params->beta, i) * gsl_matrix_get(params->dpp, i, j * (data->nGens-1) + k);
      }
      mudp += 2 * sumdpA;
      gsl_vector_set(mudpv, j * (data->nGens - 1) + k, 2 * sumdpA);
      val = gsl_matrix_get(data->envData, k, j);
      gsl_vector_set(envv, j * (data->nGens - 1) + k, val);
      envcov += val;
    }
  }

  params->ssMudpBv = mudp;
  
  //compute means
  mudp = mudp / (double) (data->nPops * (data->nGens - 1));
  envcov = envcov / (double) (data->nPops * (data->nGens - 1));

  // compute covariance
  for(j=0; j<data->nPops; j++){
    for(k=0; k<(data->nGens-1); k++){
      cov += (gsl_vector_get(mudpv, j * (data->nGens - 1) + k) - mudp) *
	(gsl_vector_get(envv, j * (data->nGens - 1) + k) - envcov);
    }
  }
      
  cov = cov / (double) (data->nPops * (data->nGens - 1));
  
  params->ssEnvcovBv = cov; 
  
  gsl_vector_free(mudpv);
  gsl_vector_free(envv);

}
    
    
// write parameter values and summary statistics to the outfile
void writesim(dataset * data, param * params, FILE * OUT, int ss){
  int i;
  
  // write param values
  //  fprintf(OUT, "%.3f %.3f %.3f", params->w,  params->fsA, params->fsB);
  fprintf(OUT, "%.3f %.3f", params->fsA, params->fsB);
  if(params->smod > 0)
    fprintf(OUT, " %.3f", params->fsC);
  // write derived params, note si is abs(si)
  fprintf(OUT, " %.3f %.3f %.6f %.6f", params->Smn, params->Ssd,
	  params->simn, params->sisd);
  
  if((ss == 0) | (ss == 2)){
    fprintf(OUT, " %.3f %.5f", params->ssMudpBv, params->ssEnvcovBv);
  }
  if(ss > 0){
    // write ss, sum of dp
    for(i=0; i<data->nLoci; i++){
      fprintf(OUT, " %.3f", gsl_vector_get(params->ssSumdp, i));
    }
    // write ss, dp x env covariance
    for(i=0; i<data->nLoci; i++){
      fprintf(OUT, " %.5f", gsl_vector_get(params->ssEnvcov, i));
    }
  }
  fprintf(OUT,"\n");

}

// validate, simulate data and generate new ss
void valsim(dataset * data, param * params, FILE * OUT){

  int i, j, k, q;
  int xx;
  double x;
  double qn, qval;
  double sumdpA;
  
  gsl_vector * mudpv;
  
  // sample params from priors
  // QTL effects
  samqtl(data, params);

  // sample posterior param values
  x = gsl_ran_flat(r, 0, (data->params->size1-1));
  xx = (int) floor(x);
  params->fsA = gsl_matrix_get(data->params, xx, 0);
  params->fsB = gsl_matrix_get(data->params, xx, 1);

  if(params->smod > 0) // sigmoid or step function must have slope
    params->fsC = gsl_matrix_get(data->params, xx, 2);
 
  // conduct the simulations
  for(j=0; j<data->nPops; j++){
    simpop(data, params, j);		   
  }

  // compute validation ss
  
  mudpv = gsl_vector_calloc(data->nPops * (data->nGens - 1));
  
  for(j=0; j<data->nPops; j++){
    for(k=0; k<(data->nGens-1); k++){
      sumdpA = 0;
      for(i=0; i<data->nLoci; i++){
	sumdpA += gsl_vector_get(params->beta, i) * gsl_matrix_get(params->dpp, i, j *
								   (data->nGens-1) + k);
      }
      gsl_vector_set(mudpv, j * (data->nGens - 1) + k, 2 * sumdpA);
    }
  }
  
  gsl_sort_vector(mudpv);

      
  for(q=1; q<=9; q++){
    qn = (double) q/10.0;
    qval = gsl_stats_quantile_from_sorted_data(mudpv->data,1,mudpv->size,qn);
    fprintf(OUT,"%.5f ", qval);
  }
  fprintf(OUT,"\n");
  
  gsl_vector_free(mudpv);
}

// wrapper function for observed summary statistics
void calcobs(dataset * data, param * params, FILE * OUT, int ss){
  int i;
  
  // QTL effects
  samqtl(data, params); // users can supply model-averaged with pp 1 instead
  
  // calculate obs.  changes in allele frequncy
  calcobsdpp(data, params);
  
  // calculate and store the summary statistics
  for(i=0; i<data->nLoci; i++){
    if(ss > 0)
      calcss(data, params, i);
  }
  
  if((ss == 0) | (ss == 2))
    calcssbv(data, params);

  // write results
  if((ss == 0) | (ss == 2)){
    fprintf(OUT, " %.3f %.5f", params->ssMudpBv, params->ssEnvcovBv);
  }
  if(ss > 0){
    // write ss, sum of dp
    for(i=0; i<data->nLoci; i++){
      fprintf(OUT, " %.3f", gsl_vector_get(params->ssSumdp, i));
    }
    // write ss, dp x env covariance
    for(i=0; i<data->nLoci; i++){
      fprintf(OUT, " %.5f", gsl_vector_get(params->ssEnvcov, i));
    }
  }
  fprintf(OUT,"\n");
   
}
  
// calculate observed change in allele frequency
void calcobsdpp(dataset * data, param * params){
  int i, j, k;
  double dp, pp0, pp1;

  for(j=0; j<data->nPops; j++){
    for(k=0; k<(data->nGens-1); k++){
      for(i=0; i<data->nLoci; i++){
	pp0 = gsl_matrix_get(data->pall, i, j * data->nGens + k);
	pp1 = gsl_matrix_get(data->pall, i, j * data->nGens + k + 1);
	dp = pp1 - pp0;
	gsl_matrix_set(params->dpp, i, j * (data->nGens - 1) + k, dp);
      }
    }
  }
}
