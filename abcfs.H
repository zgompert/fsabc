#include <iostream>
#include <sstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_sort.h>

// macro definitions
#define VERSION "0.2 -- 01 June 2020"

extern gsl_rng * r;
// function declarations
using namespace std;

struct dataset{
  int nLoci;
  int nPops;
  int nGens;
  int nNeSams;

  double sigma2; // phenotypic variance (within a population and generation)
  // could make this population specific
    
  gsl_matrix * p0; // initial allele freq, nLoci x nPops
  gsl_matrix * pall; // for obs. ss, allele freqs, nLoic x (nPops x nGens); all gens for pop 0 first
  gsl_matrix * envData; // environmental data, nGens(-1) x nPops
  gsl_matrix_int * sampledData; // binary indicator  varaible, denotes number of actual generations between samples, nGens(-1) x nPops
  gsl_matrix * ne; // external samples of Ne, nNeSams x nPops; assumes constant per pop
  gsl_matrix * qtl; // pip and beta, nrow = nLoci, ncol = 2 pip then beta

  gsl_matrix * params; // only for validation, post. estimates; samples x params
};

struct param{
  double prsel; // prior prob. fsA or fsB != 0
  double fsA; // intercept for linear model for peak of fitness function
  double fsB; // slope for env./linear model for peak of fitness function; or max for signmoid or step model
  double fsC; // param for cut point or slope for sigmoid or step model
  double fsAlb; // lower bound on U prior for fsA
  double fsAub; // upper bound on U prior for fsA
  double fsBlb; // lower bound on U prior for fsB
  double fsBub; // upper bound on U prior for fsB
  double fsClb; // lower bound on U prior for fsC // sign might be dictated by A vs. B for sigmoid model 
  double fsCub; // upper bound on U prior for fsC

  unsigned int smod; // selection model, 0 = linear, 1 = step, 2 = sigmoidal
 
  double ssMudpBv; // ss, sum change in the mean bv
  double ssEnvcovBv; // ss, cov between env and mean bv change
    
  gsl_vector * beta; // sample of allelic effects
  gsl_matrix * pp; // matrix of sim. allele freqs for one simulation; row = loci
  // columns = pop 0, gen 0, gen 1, ... pop 1, gen 0, gen ...
  gsl_matrix * dpp; // matrix of sim. allele freqs for one simulation; row = loci
  // columns = pop 0, gen 1 - gen 0, gen 2 - gen 1, pop 1 ...
  gsl_vector * ssSumdp; // ss, sum change, by locus
  gsl_vector * ssEnvcov; // ss, cov between env. and dp, by locus

  double Nqtl; // number of QTL sampled for a given simulation
  
  gsl_vector * Svec; // S = selection differential for each pop. and gen.
  gsl_vector * sivec; // si (locus sel.) for each qtl x pop. x gen.

  // mean and sd  of S (sel. differential) across generations and pops.
  double Smn;
  double Ssd;

  // mean and sd of |si| (selection coef. on SNPs) across generations and pops and causal loci/QTL.
  double simn;
  double sisd;
};

// function definitions
void usage(char * name);
void getdata(string geneFile, string envFile, string missingFile, int obsss, dataset * data);
void getne(string neFile, dataset * data);
void getqtl(string traitFile, dataset * data);
void getest(string resFile, dataset * data);
void runsim(dataset * data, param * params, int x, int ss);
void samqtl(dataset * data, param * params);
void simpop(dataset * data, param * params, int j);
void calcss(dataset * data, param * params, int i);
void calcssbv(dataset * data, param * params);
void writesim(dataset * data, param * params, FILE * OUT, int ss);
void compdpar(dataset * data, param * params);
void valsim(dataset * data, param * params, FILE * OUT);
void calcobsdpp(dataset * data, param * params);
void calcobs(dataset * data, param * params, FILE * OUT, int ss);
