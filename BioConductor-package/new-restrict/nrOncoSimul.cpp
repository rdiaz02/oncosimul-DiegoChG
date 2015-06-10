// #include "OncoSimul.h"
#include "debug_common.hpp"
#include "bnb_common.hpp"
#include "new_restrict.hpp"
#include <Rcpp.h>
#include <limits>
#include <iostream>
// #include <gsl/gsl_rng.h> // here? in the .h
#include <random>
#include <set>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <ctime>
// #include <sys/resource.h> 
#include <sys/time.h> 

// #include <exception>
#include <stdexcept>

using namespace Rcpp;
using std::vector;

// To track if mutation is really much smaller than birth/death
#define MIN_RATIO_MUTS_NR
#ifdef MIN_RATIO_MUTS_NR
// There is really no need for these to be globals?
// Unless I wanted to use them inside some function. So leave as globals.
double g_min_birth_mut_ratio_nr = DBL_MAX;
double g_min_death_mut_ratio_nr = DBL_MAX;
double g_tmp1_nr = DBL_MAX;
#endif


static void nr_fitness(spParamsP& tmpP,
		       const spParamsP& parentP,
		       const Genotype& ge,
		       const fitnessEffectsAll& F,
		       const std::string& typeFitness,
		       const double& genTime,
		       const double& adjust_fitness_B,
		       const double& adjust_fitness_MF) {

  // We want a way to signal immediate non-viability of a clone. For
  // "birth-based" models that happens when any s = -1, as the fitness is
  // 0. By setting birth = 0.0 we ensure this clone does not get added and
  // we never reach into algo2, etc, leading to numerical problems.

  // With Bozic models, which are "death-based", it is different. For
  // bozic2, birth is bounded, so any death > 2 would lead to birth <
  // 0. For bozic1, deaths of around 50 lead to numerical issues.  The
  // general rule is: set those mutations to -90, so prodDeathFitness
  // immediately returns a 99.0 for death, and that is recognized as "no
  // viability".

  // The ones often used are bozic1, exp, mcfarlandlog

  if(typeFitness == "bozic1") {
    tmpP.death = prodDeathFitness(evalGenotypeFitness(ge, F));
    if( tmpP.death > 50) {
      tmpP.birth = 0.0; 
    } else {
      tmpP.birth = 1.0;
    }
  } else if (typeFitness == "bozic2") {
    double pp = prodDeathFitness(evalGenotypeFitness(ge, F));
    tmpP.birth = std::max(0.0, (1.0/genTime) * (1.0 - 0.5 * pp ));
    tmpP.death = (0.5/genTime) * pp;
  } else {
    double fitness = prodFitness(evalGenotypeFitness(ge, F));
    if( fitness <= 0.0) {
      tmpP.absfitness = 0.0;
      tmpP.death = 1.0;
      tmpP.birth = 0.0; 
    } else{
      // Set appropriate defaults and change only as needed
      tmpP.death = parentP.death;
      tmpP.absfitness = parentP.absfitness;
      tmpP.birth = fitness;
      // exp, mcfarland, and mcfarlandlog as above. Next are the two exceptions.
      if(typeFitness == "beerenwinkel") {
	tmpP.absfitness = fitness; 
	tmpP.birth = adjust_fitness_B * tmpP.absfitness;
      } else if(typeFitness == "mcfarland0") {
	tmpP.absfitness = fitness;
	tmpP.birth = adjust_fitness_MF * tmpP.absfitness;
      }
    }
  }
  // Exp and McFarland and McFarlandlog are also like Datta et al., 2013
  // An additional driver gene mutation increases a cell’s fitness by a
  // factor of (1+sd), whereas an additional housekeeper gene mutation
  // decreases fitness by a factor of (1-sh) and the effect of multiple
  // mutations is multiplicative
}


// is this really any faster than the one below?
inline void new_sp_v(unsigned int& sp,
		     const Genotype& newGenotype,
		     const std::vector<Genotype> Genotypes) {
  sp = 0;
  for(sp = 0; sp < Genotypes.size(); ++sp) {
    if( newGenotype == Genotypes[sp] )
      break;
  }
}

unsigned int new_sp(const Genotype& newGenotype,
		    const std::vector<Genotype> Genotypes) {
  for(unsigned int sp = 0; sp < Genotypes.size(); ++sp) {
    if( newGenotype == Genotypes[sp] ) {
      return sp;
    }
  }
  return Genotypes.size();
}

static void remove_zero_sp_nr(std::vector<int>& sp_to_remove,
			      std::vector<Genotype>& Genotypes,
			      std::vector<spParamsP>& popParams,
			      std::multimap<double, int>& mapTimes) {
  std::vector<spParamsP>::iterator popParams_begin = popParams.begin();
  std::vector<Genotype>::iterator Genotypes_begin = Genotypes.begin();
  std::vector<int>::reverse_iterator r = sp_to_remove.rbegin();
  int remove_this;
  while(r != sp_to_remove.rend() ) {
    remove_this = *r;
    mapTimes.erase(popParams[remove_this].pv);
    popParams.erase(popParams_begin + remove_this);
    Genotypes.erase(Genotypes_begin + remove_this);
    ++r;
  }
}

// FIXME: change this, now that we keep a count of drivers?
// see the new function in new-restrict.cpp: countDrivers

// Yes, we do want to count drivers. For instance, stopping in cancer can
// be related to this. So:

// In R, the user says which are the drivers. If does not say anuthing,
// the default (NULL) then drivers are all in poset, epist, restrict. The
// user can pass a vector with the names of the genes (not modules). Allow
// also for empty, so this is faster if not needed. And check that if we
// use a stopping rule based on drivers that drv vectors is not empty.

static inline void nr_count_NumDrivers(int& maxNumDrivers, 
				    std::vector<int>& countByDriver,
				    Rcpp::IntegerMatrix& returnGenotypes,
				    const vector<int>& drv){
  // Fill up the "countByDriver" table and return the maximum number of
  // mutated drivers in any genotype.
  // Difference w.r.t. to former is passing drv
  maxNumDrivers = 0;
  int tmpdr = 0;
  
  for(int j = 0; j < returnGenotypes.ncol(); ++j) {
    tmpdr = 0;
    for(int i : drv) {
      tmpdr += returnGenotypes(i - 1, j);
      countByDriver[i] += returnGenotypes(i - 1, j);
    }
    if(tmpdr > maxNumDrivers) maxNumDrivers = tmpdr;
  }
}
      
// FIXME: we do this often. Why not just keep it in the struct?
int nr_count_NDrivers(const Genotype& ge, const vector<int>& drv) {
  // Counts the number of mutated drivers in a genotype.
  // drv comes from R, and it is the vector with the
  // numbers of the genes, not modules.
  return presentDrivers(ge, drv).size();
}
// FIXME: the count_NumDrivers counts for each driver. Write that too.


static void nr_totPopSize_and_fill_out_crude_P(int& outNS_i,
					    double& totPopSize, 
					    double& lastStoredSample,
					    std::vector<Genotype>& genot_out,
					    //std::vector<unsigned long long>& sp_id_out,
					    std::vector<double>& popSizes_out,
					    std::vector<int>& index_out,
					    std::vector<double>& time_out,
					    std::vector<double>& sampleTotPopSize,
					    std::vector<double>& sampleLargestPopSize,
					    std::vector<int>& sampleMaxNDr,
					    std::vector<int>& sampleNDrLargestPop,
					    bool& simulsDone,
					    bool& reachDetection,
					    int& lastMaxDr,
					    double& done_at,
					    const std::vector<Genotype>& Genotypes,
					    const std::vector<spParamsP>& popParams, 
					    const double& currentTime,
					    const double& keepEvery,
					    const double& detectionSize,
					    const double& finalTime,
					    // const double& endTimeEvery,
					    const int& detectionDrivers,
					    const int& verbosity,
					    const double& minDDrPopSize,
					    const double& extraTime,
					       const vector<int>& drv,
					    const double& fatalPopSize = 1e15) {
  // Fill out, but also compute totPopSize
  // and return sample summaries for popsize, drivers.
  
  // This determines if we are done or not by checking popSize, number of
  // drivers, etc
  
  // static int lastMaxDr = 0; // preserves value across calls to Algo5 from R.
  // so can not use it.
  bool storeThis = false;
  totPopSize = 0.0;
  
  // this could all be part of popSize_over_m_dr, with a better name
  int tmp_ndr = 0;
  int max_ndr = 0;
  double popSizeOverDDr = 0.0;

  for(size_t i = 0; i < popParams.size(); ++i) {
    totPopSize += popParams[i].popSize;
    tmp_ndr = nr_count_NDrivers(Genotypes[i], drv);
    if(tmp_ndr > max_ndr) max_ndr = tmp_ndr;
    if(tmp_ndr >= detectionDrivers) popSizeOverDDr += popParams[i].popSize;
  }
  lastMaxDr = max_ndr;

  
  if (keepEvery < 0) {
    storeThis = false;
  } else if( currentTime >= (lastStoredSample + keepEvery) ) {
    storeThis = true;
  }

  if( (totPopSize <= 0.0) || (currentTime >= finalTime)  ) {
    simulsDone = true;
  }

  if(extraTime > 0) {
    if(done_at <  0) {
      if( (totPopSize >= detectionSize) ||
	  ( (lastMaxDr >= detectionDrivers) &&
	    (popSizeOverDDr >= minDDrPopSize) ) ) {
	done_at = currentTime + extraTime;
      }
    } else if (currentTime >= done_at) {
	simulsDone = true;
	reachDetection = true; 
      }
  } else if( (totPopSize >= detectionSize) ||
	     ( (lastMaxDr >= detectionDrivers) &&
	       (popSizeOverDDr >= minDDrPopSize) ) ) {
    simulsDone = true;
    reachDetection = true; 
  }
  
  if(totPopSize >= fatalPopSize) {
    Rcpp::Rcout << "\n\totPopSize > " << fatalPopSize
		<<". You are likely to loose precision and run into numerical issues\n";
       }
  
  if(simulsDone)
    storeThis = true;


  if( storeThis ) {
    lastStoredSample = currentTime;
    outNS_i++;
    int ndr_lp = 0;
    double l_pop_s = 0.0;
    
    time_out.push_back(currentTime);
    
    for(size_t i = 0; i < popParams.size(); ++i) {
      genot_out.push_back(Genotypes[i]);
      popSizes_out.push_back(popParams[i].popSize);
      index_out.push_back(outNS_i);
      
      if(popParams[i].popSize > l_pop_s) {
	l_pop_s = popParams[i].popSize;
	ndr_lp = nr_count_NDrivers(Genotypes[i], drv);
      }
    }
    sampleTotPopSize.push_back(totPopSize);
    sampleLargestPopSize.push_back(l_pop_s);
    sampleMaxNDr.push_back(max_ndr);
    sampleNDrLargestPop.push_back(ndr_lp);
  } 
  
  if( !std::isfinite(totPopSize) ) {
    throw std::range_error("totPopSize not finite");
  }
  if( std::isnan(totPopSize) ) {
    throw std::range_error("totPopSize is NaN");
  }
  
  if(totPopSize > (4.0 * 1e15)) {
    if(verbosity > 0)
      Rcpp::Rcout << "\nWARNING: popSize > 4e15. Likely loss of precission\n";
  }
}

// FIXME: I might want to return the actual drivers in each period
// and the actual drivers in the population with largest popsize
// Something like what we do now with whichDrivers
// and count_NumDrivers




inline void nr_reshape_to_outNS(Rcpp::NumericMatrix& outNS,
				const vector<vector<int> >& uniqueGenotV,
				const vector<vector<int> >& genot_out_v,
				const vector<double>& popSizes_out,
				const vector<int>& index_out,
				const vector<double>& time_out){
  
  vector<vector<int> >::const_iterator fbeg = uniqueGenotV.begin();
  vector<vector<int> >::const_iterator fend = uniqueGenotV.end();

  int column;

  for(size_t i = 0; i < genot_out_v.size(); ++i) {
    column = std::distance(fbeg, lower_bound(fbeg, fend, genot_out_v[i]) );
    outNS(index_out[i], column + 1) =  popSizes_out[i];
  }

  for(size_t j = 0; j < time_out.size(); ++j)
    outNS(j, 0) = time_out[j];
}


Rcpp::NumericMatrix create_outNS(const vector<vector<int> >& uniqueGenotypes,
				 const vector<vector<int> >& genot_out_v,
				 const vector<double>& popSizes_out,
				 const vector<int>& index_out,
				 const vector<double>& time_out,
				 const int outNS_i, const int maxram) {
  // The out.ns in R code; holder of popSizes over time
  // The first row is time, then the genotypes (in column major)
  // here("after uniqueGenotypes_to_vector");

  int outNS_r, outNS_c, create_outNS;
  if( ( (uniqueGenotypes.size() + 1) *  (outNS_i + 1) ) > ( pow(2, 31) - 1 ) ) {
    Rcpp::Rcout << "\nWARNING: Return outNS object > 2^31 - 1. Not created.\n";
    outNS_r = 1;
    outNS_c = 1;
    create_outNS = 0;
  } else if ( 
	     static_cast<long>((uniqueGenotypes.size()+1) * (outNS_i+1)) * 8 > 
	     (maxram * (1024*1024) ) ) {
    Rcpp::Rcout << "\nWARNING: Return outNS object > maxram. Not created.\n";
    outNS_r = 1;
    outNS_c = 1;
    create_outNS = 0;
  } else {
    outNS_r = outNS_i + 1;
    outNS_c = uniqueGenotypes.size() + 1;
    create_outNS = 1;
  }
  Rcpp::NumericMatrix outNS(outNS_r, outNS_c);  
  if(create_outNS) {
    nr_reshape_to_outNS(outNS, uniqueGenotypes,
			genot_out_v, 
			popSizes_out, 
			index_out, time_out);
    
  } else {
    outNS(0, 0) = -99;
  }
  return outNS;
}



// FIXME: when creating the 0/1, collapse those that are the same

// vector< vector<int> > uniqueGenot_vector(vector<Genotype>& genot_out) {
//   // From genot_out we want the unique genotypes, but each as a single
//   // vector. Convert to the vector, then use a set to give unique sorted
//   // vector.
//   std::vector<std::vector<int> > genot_out_nr;
//   std::transform(genout_out.begin(), genot_out.end(),
// 		 back_inserter(genout_out_nr),
// 		 genotypeSingleVector);
//   std::set<std::vector<int> > uniqueGenotypes_nr(genot_out_nr.begin(), genot_out_nr.end());
//   std::vector<std::vector<int> > uniqueGenotypes_vector_nr (uniqueGenotypes_nr.begin(),
// 							    uniqueGenotypes_nr.end());
//   return uniqueGenotypes_vector_nr;
// }

vector< vector<int> > uniqueGenot_vector(vector<vector<int> >& genot_out) {
  // From genot_out we want the unique genotypes, but each as a single
  // vector. Convert to the vector, then use a set to give unique sorted
  // vector.
  std::set<std::vector<int> > uniqueGenotypes_nr(genot_out.begin(),
						 genot_out.end());
  std::vector<std::vector<int> > uniqueGenotypes_vector_nr (uniqueGenotypes_nr.begin(),
							    uniqueGenotypes_nr.end());
  return uniqueGenotypes_vector_nr;
}


// std::set<std::vector<int> > nr_find_unique_genotypes(const std::vector<unsigned long long>& genot_out_l) {
//   std::set<std::vector<int> > uniqueGenotypes;
//   for( auto gg : genot_out_nr)
//     uniqueGenotypes.insert(gg);
//   return uniqueGenotypes;
// }



// static inline void genot_out_to_ullong(std::vector<unsigned long long>& go_l,
// 			       const std::vector<Genotype64>& go) {
//   for(size_t i = 0; i < go.size(); ++i)
//     go_l[i] = go[i].to_ullong();
// }

std::vector<std::vector<int> > genot_to_vectorg(const std::vector<Genotype>& go) {
  std::vector<std::vector<int> > go_l;
  std::transform(go.begin(), go.end(), back_inserter(go_l), genotypeSingleVector);
  return go_l;
}


// std::vector<std::vector<int> > nr_uniqueGenotypes_to_vector(const std::set< std::vector<int> >& uniqueGenotypes_nr) {
//   std::vector<std::vector<int> > ugV(uniqueGenotypes_nr.begin(),
// 				     uniqueGenotypes_nr.end());
//   return ugV;
// }


std::string vectorGenotypeToString(const std::vector<int>& genotypeV,
				   const fitness_as_genes& fg) {
  
  std::string strGenotype;

  std::vector<int> order_int;
  std::vector<int> rest_int;

  for(auto g : genotypeV) {
    if( binary_search(fg.orderG.begin(), fg.orderG.end(), g)) {
      order_int.push_back(g);
    } else {
      rest_int.push_back(g);
    }
  }

  std::string order_sep = "_";
  std::string order_part;
  std::string rest;
  std::string comma = "";
  for(auto g : order_int) {
    order_part += (comma + std::to_string(g));
    comma = ", ";
  }
  comma = "";
  for(auto g : rest_int) {
    rest += (comma + std::to_string(g));
    comma = ", ";
  }
  if(fg.orderG.size()) {
    strGenotype = order_part + order_sep + rest;
  } else {
    strGenotype = rest;
  }
  return strGenotype;
}


std::vector<std::string> genotypesToString(const std::vector< vector<int> >& uniqueGenotypesV,
					   const fitnessEffectsAll& F) {
  fitness_as_genes fg = feGenes(F);
  std::vector<std::string> gs;

  for(auto v: uniqueGenotypesV ) gs.push_back(vectorGenotypeToString(v, fg));

  // exercise: do it with lambdas
  // std::transform(uniqueGenotypesV.begin(), uniqueGenotypesV.end(),
  // 		 back_inserter(gs), vectorGenotypeToString);
  return gs;
}

Rcpp::IntegerMatrix nr_create_returnGenotypes(const int& numGenes,
					      const std::vector< vector<int> >& uniqueGenotypesV){
  // We loose order here. Thus, there might be several identical columns.
  Rcpp::IntegerMatrix returnGenotypes(numGenes, uniqueGenotypesV.size());
  for(size_t i = 0; i < uniqueGenotypesV.size(); ++i) {
    for(int j : uniqueGenotypesV[i]) {
      returnGenotypes(j - 1, i) = 1;
    }
  }
  return returnGenotypes;
}






static void nr_sample_all_pop_P(std::vector<int>& sp_to_remove,
				std::vector<spParamsP>& popParams,
				// next only used with DEBUGV
				const std::vector<Genotype>& Genotypes,
				const double& tSample){
  sp_to_remove.clear();

  for(size_t i = 0; i < popParams.size(); i++) {
    STOPASSERT(popParams[i].timeLastUpdate >= 0.0);
    STOPASSERT(tSample - popParams[i].timeLastUpdate >= 0.0);
#ifdef DEBUGV
    Rcpp::Rcout << "\n\n     ********* 5.9 ******\n " 
	      << "     Species  = " << i 
	      << "\n      Genotype = " << genotypeSingleVector(Genotypes[i])
      //	      << "\n      sp_id = " << genotypeSingleVector(Genotypes[i]) // sp_id[i]  
	      << "\n      pre-update popSize = " 
	      << popParams[i].popSize 
	      << "\n      time of sample = " << tSample 
	      << "\n      popParams[i].timeLastUpdate = " 
	      << popParams[i].timeLastUpdate 
	      << ";\n     t for Algo2 = " 
	      << tSample - popParams[i].timeLastUpdate 
	      << " \n     species R " << popParams[i].R
	      << " \n     species W " << popParams[i].W
	      << " \n     species death " << popParams[i].death
	      << " \n     species birth " << popParams[i].birth;
#endif

    // Account for forceSampling. When 
    // forceSampling, popSize for at least one species
    // was updated in previous loop, so we skip that one
    if(tSample > popParams[i].timeLastUpdate) {
      popParams[i].popSize = 
	Algo2_st(popParams[i], tSample);
    }
    if( popParams[i].popSize <=  0.0 ) {
      // this i has never been non-zero in any sampling time
      // eh??
      // If it is 0 here, remove from _current_ population. Anything that
      // has had a non-zero size at sampling time is preserved (if it
      // needs to be preserved, because it is keepEvery time).
      sp_to_remove.push_back(i);
#ifdef DEBUGV
      Rcpp::Rcout << "\n\n     Removing species i = " << i 
		  << " with genotype = " << genotypeSingleVector(Genotypes[i]);
#endif
    } 
#ifdef DEBUGV
    Rcpp::Rcout << "\n\n   post-update popSize = " 
	      << popParams[i].popSize << "\n";
#endif
  }
}






static void nr_innerBNB(const fitnessEffectsAll& fitnessEffects,
			const double& initSize,
		     const double& K,
		     const double& alpha,
		     const double& genTime,
		     const std::string& typeFitness,
		     const int& mutatorGenotype,
		     const double& mu,
		     const double& death,
		     const double& keepEvery,
		     const double& sampleEvery,		     
			const std::vector<int>& initMutant,
		     const time_t& start_time,
		     const double& maxWallTime,
		     const double& finalTime,
		     const double& detectionSize,
		     const int& detectionDrivers,
		     const double& minDDrPopSize,
		     const double& extraTime,
		     const int& verbosity,
		     double& totPopSize,
		     double& e1,
		     double& n_0,
		     double& n_1,
		     double& ratioForce,
		     double& currentTime,
		     int& speciesFS,
		     int& outNS_i,
		     int& iter,
		     std::vector<Genotype>& genot_out,
		     std::vector<double>& popSizes_out,
		     std::vector<int>& index_out,
		     std::vector<double>& time_out,
		     std::vector<double>& sampleTotPopSize,
		     std::vector<double>& sampleLargestPopSize,
		     std::vector<int>& sampleMaxNDr,
		     std::vector<int>& sampleNDrLargestPop,
		     bool& reachDetection,
		     std::mt19937& ran_gen,
		     double& runningWallTime,
		     bool& hittedWallTime) {
		     //bool& anyForceRerunIssues
  //  if(numRuns > 0) {

 
  const int numGenes = fitnessEffects.genomeSize;
  double dummyMutationRate = std::max(mu/1000, 1e-13);
  // double dummyMutationRate = 1e-10;
  // ALWAYS initialize this here, or reinit or rezero
  genot_out.clear();
  popSizes_out.clear();
  index_out.clear();
  time_out.clear();
  totPopSize = 0.0;
  sampleTotPopSize.clear();
  currentTime = 0.0;
  iter = 0;

  outNS_i = -1;

  sampleTotPopSize.clear();
  sampleLargestPopSize.clear();
  sampleMaxNDr.clear();
  sampleNDrLargestPop.clear();
  // end of rezeroing.

  
  // }
  // anyForceRerunIssues = false;
  
  bool forceSample = false;
  bool simulsDone = false;
  double lastStoredSample = 0.0;


  double minNextMutationTime;
  double mutantTimeSinceLastUpdate;
  double timeNextPopSample;
  double tSample;

  std::vector<int> newMutations;
  int nextMutant;
  unsigned int numSpecies = 0;
  int numMutablePosParent = 0;


  unsigned int sp = 0;
  //int type_resize = 0;

  int iterL = 1000;
  int speciesL = 1000; 
  //int timeL = 1000;
  
  double tmpdouble1 = 0.0;
  double tmpdouble2 = 0.0;
  
  std::vector<int>sp_to_remove(1);
  sp_to_remove.reserve(10000);

  // those to update
  int to_update = 1; //1: one species; 2: 2 species; 3: all.
  int u_1 = -99;
  int u_2 = -99;

  Genotype newGenotype;
  std::vector<Genotype> Genotypes(1);
  Genotypes[0] = wtGenotype(); //Not needed, but be explicit.
  
  std::vector<spParamsP> popParams(1);

      // // FIXME debug
      // Rcpp::Rcout << "\n popSize[0]  at 1 ";
      // print_spP(popParams[0]);
      // // end debug
  
  
  const int sp_per_period = 5000;
  popParams.reserve(sp_per_period);
  Genotypes.reserve(sp_per_period);

      // // FIXME debug
      // Rcpp::Rcout << "\n popSize[0]  at 01 ";
      // print_spP(popParams[0]);
      // // end debug

  
  spParamsP tmpParam; 
  init_tmpP(tmpParam);
  init_tmpP(popParams[0]);
  popParams[0].popSize = initSize;
  totPopSize = initSize;


      // // FIXME debug
      // Rcpp::Rcout << "\n popSize[0]  at 10000 ";
      // print_spP(popParams[0]);
      // // end debug


  
  std::vector<int>mutablePos(numGenes); // could be inside getMuatedPos_bitset

    // multimap to hold nextMutationTime
  std::multimap<double, int> mapTimes;
  //std::multimap<double, int>::iterator m1pos;

  int ti_dbl_min = 0;
  int ti_e3 = 0;


      // // FIXME debug
      // Rcpp::Rcout << "\n popSize[0]  at 10002 ";
      // print_spP(popParams[0]);
      // // end debug


  
  // Beerenwinkel
  double adjust_fitness_B = -std::numeric_limits<double>::infinity();
  //McFarland
  double adjust_fitness_MF = -std::numeric_limits<double>::infinity();

  // for McFarland error
  e1 = 0.0;
  n_0 = 0.0;
  n_1 = 0.0;
  double tps_0, tps_1; 
  tps_0 = totPopSize;
  tps_1 = totPopSize;


      // // FIXME debug
      // Rcpp::Rcout << "\n popSize[0]  at 10004 ";
      // print_spP(popParams[0]);
      // // end debug

  
  
  int lastMaxDr = 0;
  double done_at = -9;

#ifdef MIN_RATIO_MUTS_NR
  g_min_birth_mut_ratio_nr = DBL_MAX;
  g_min_death_mut_ratio_nr = DBL_MAX;
  g_tmp1_nr = DBL_MAX;
#endif

      // // FIXME debug
      // Rcpp::Rcout << " popSize[0]  at 1b ";
      // print_spP(popParams[0]);
      // // end debug

    // This long block, from here to X1, is ugly and a mess!
  // This is what takes longer to figure out whenever I change
  // anything. FIXME!!
  if(initMutant.size() > 0) {
    
    popParams[0].numMutablePos = numGenes - 1;
    Genotypes[0] = createNewGenotype(wtGenotype(),
				     initMutant,
				     fitnessEffects,
				     ran_gen); // FIXME: nr, here. What is a "wt
					// genotype"? Does it have "0"
					// mutated, or nothing. Nothing.
    if(typeFitness == "beerenwinkel") {
      
      popParams[0].death = 1.0; //note same is in McFarland.
      // But makes sense here; adjustment in beerenwinkel is via fitness
      
      // initialize to prevent birth/mutation warning with Beerenwinkel
      // when no mutator. O.w., the defaults
      if(!mutatorGenotype)
	popParams[0].mutation = mu * popParams[0].numMutablePos;
      popParams[0].absfitness = prodFitness(evalGenotypeFitness(Genotypes[0],
								fitnessEffects));
      updateRatesBeeren(popParams, adjust_fitness_B, initSize,
			currentTime, alpha, initSize, 
			mutatorGenotype, mu);
    } else if(typeFitness == "mcfarland0") {
      // death equal to birth of a non-mutant.
      popParams[0].death = log1p(totPopSize/K); // log(2.0), except rare cases
      if(!mutatorGenotype)
	popParams[0].mutation = mu * popParams[0].numMutablePos;
      popParams[0].absfitness = prodFitness(evalGenotypeFitness(Genotypes[0],
								fitnessEffects));
      updateRatesMcFarland0(popParams, adjust_fitness_MF, K, 
			    totPopSize,
			    mutatorGenotype, mu);
    } else if(typeFitness == "mcfarland") {
      popParams[0].death = totPopSize/K;
      popParams[0].birth = prodFitness(evalGenotypeFitness(Genotypes[0],
								fitnessEffects));
    } else if(typeFitness == "mcfarlandlog") {
      popParams[0].death = log1p(totPopSize/K);
      popParams[0].birth = prodFitness(evalGenotypeFitness(Genotypes[0],
								fitnessEffects));
    } else if(typeFitness == "bozic1") {
      tmpParam.birth =  1.0;
      tmpParam.death = -99.9;
    } else if (typeFitness == "bozic2") {
      tmpParam.birth =  -99;
      tmpParam.death = -99;
    } else if (typeFitness == "exp") {
      tmpParam.birth =  -99;
      tmpParam.death = death;
    } else {
      throw std::invalid_argument("this ain't a valid typeFitness");
    } 
    if( (typeFitness != "beerenwinkel") && (typeFitness != "mcfarland0") 
	&& (typeFitness != "mcfarland") && (typeFitness != "mcfarlandlog")) // wouldn't matter
      nr_fitness(popParams[0], tmpParam,
		 Genotypes[0],
		 fitnessEffects,
		 typeFitness, genTime,
		 adjust_fitness_B, adjust_fitness_MF);
    // we pass as the parent the tmpParam; it better initialize
    // everything right, or that will blow. Reset to init
    init_tmpP(tmpParam);
  } else {
    popParams[0].numMutablePos = numGenes;
    if(typeFitness == "beerenwinkel") {
      popParams[0].death = 1.0;
      // initialize to prevent birth/mutation warning with Beerenwinkel
      // when no mutator. O.w., the defaults
      if(!mutatorGenotype)
	popParams[0].mutation = mu * popParams[0].numMutablePos;
      popParams[0].absfitness = 1.0;
      updateRatesBeeren(popParams, adjust_fitness_B, initSize,
			currentTime, alpha, initSize, 
			mutatorGenotype, mu);
    } else if(typeFitness == "mcfarland0") {
      popParams[0].death = log1p(totPopSize/K);
      if(!mutatorGenotype)
	popParams[0].mutation = mu * popParams[0].numMutablePos;
      popParams[0].absfitness = 1.0;
      updateRatesMcFarland0(popParams, adjust_fitness_MF, K, 
			    totPopSize,
			    mutatorGenotype, mu);
    } else if(typeFitness == "mcfarland") {
      popParams[0].birth = 1.0;
      popParams[0].death = totPopSize/K;
      // no need to call updateRates
    } else if(typeFitness == "mcfarlandlog") {
      popParams[0].birth = 1.0;
      popParams[0].death = log1p(totPopSize/K);
      // no need to call updateRates
    } else if(typeFitness == "bozic1") {
      popParams[0].birth = 1.0;
      popParams[0].death = 1.0;
    } else if (typeFitness == "bozic2") {
      popParams[0].birth = 0.5/genTime;
      popParams[0].death = 0.5/genTime;
    } else if (typeFitness == "exp") {
      popParams[0].birth = 1.0;
      popParams[0].death = death;
    } else {
      throw std::invalid_argument("this ain't a valid typeFitness");
    }
  }


  // // FIXME debug
  //     Rcpp::Rcout << " popSize[0]  at 2 ";
  //     print_spP(popParams[0]);
  //     // end debug
  

  
  // these lines (up to, and including, R_F_st)
  // not needed with mcfarland0 or beerenwinkel
  if(mutatorGenotype)
    popParams[0].mutation = mu * popParams[0].birth * popParams[0].numMutablePos;
  else
    popParams[0].mutation = mu * popParams[0].numMutablePos;

  W_f_st(popParams[0]);
  R_f_st(popParams[0]);

  // // FIXME debug
  //     Rcpp::Rcout << " popSize[0]  at 3 ";
  //     print_spP(popParams[0]);
  //     // end debug
  

  // X1: end of mess of initialization block

  popParams[0].pv = mapTimes.insert(std::make_pair(-999, 0));

  if( keepEvery > 0 ) {
    // We keep the first ONLY if we are storing more than one.
    outNS_i++;
    time_out.push_back(currentTime);
    
    genot_out.push_back(Genotypes[0]);
    popSizes_out.push_back(popParams[0].popSize);
    index_out.push_back(outNS_i);

    sampleTotPopSize.push_back(popParams[0].popSize);
    sampleLargestPopSize.push_back(popParams[0].popSize);
    sampleMaxNDr.push_back(nr_count_NDrivers(Genotypes[0],
					     fitnessEffects.drv));
    sampleNDrLargestPop.push_back(sampleMaxNDr[0]);
  }
  // FIXME: why next line and not just genot_out.push_back(Genotypes[i]);
  // if keepEvery > 0? We do that already.
  // It is just ugly to get a 0 in that first genotype when keepEvery < 0
  // uniqueGenotypes.insert(Genotypes[0].to_ullong());
  timeNextPopSample = currentTime + sampleEvery;
  numSpecies = 1;

  
#ifdef DEBUGV
  Rcpp::Rcout << "\n the initial species\n";
  print_spP(popParams[0]);
#endif


  // // FIXME debug
  //     Rcpp::Rcout << " popSize[0]  at 4 ";
  //     print_spP(popParams[0]);
  //     // end debug
  

  
  while(!simulsDone) {
    // Check how we are doing with time as first thing.
    runningWallTime = difftime(time(NULL), start_time);
    if( runningWallTime > maxWallTime ) {
      hittedWallTime = true;
      forceSample = true;
      simulsDone = true;
    }
    
    iter++;
    if(verbosity > 1) {
      if(! (iter % iterL) ) {
	Rcpp::Rcout << "\n\n    ... iteration " << iter;
	Rcpp::Rcout << "\n    ... currentTime " << currentTime <<"\n";
      }
      if(!(numSpecies % speciesL )) {
	Rcpp::Rcout << "\n\n    ... iteration " << iter;
	Rcpp::Rcout << "\n\n    ... numSpecies " << numSpecies << "\n";
      }
    }

    //  ************   5.2   ***************
    if(verbosity >= 2)
      Rcpp::Rcout <<"\n\n\n*** Looping through 5.2. Iter = " << iter << " \n";

    tSample = std::min(timeNextPopSample, finalTime);

#ifdef DEBUGV
    Rcpp::Rcout << " DEBUGV\n";
    Rcpp::Rcout << "\n ForceSample? " << forceSample 
	      << "  tSample " << tSample 
	      << "  currentTime " << currentTime;
#endif

    if(iter == 1) { // handle special case of first iter
      // // FIXME debug
      // Rcpp::Rcout << " popSize[0] ";
      // print_spP(popParams[0]);
      // // end debug
      tmpdouble1 = ti_nextTime_tmax_2_st(popParams[0],
					 currentTime,
					 tSample, 
					 ti_dbl_min, ti_e3);
      mapTimes_updateP(mapTimes, popParams, 0, tmpdouble1);
      //popParams[0].Flag = false;
      popParams[0].timeLastUpdate = currentTime;
    } else { // any other iter
      if(to_update == 1) {
	// we did not sample in previous period.
	tmpdouble1 = ti_nextTime_tmax_2_st(popParams[u_1],
					   currentTime,
					   tSample, 
					   ti_dbl_min, ti_e3);
	mapTimes_updateP(mapTimes, popParams, u_1, tmpdouble1);
	popParams[u_1].timeLastUpdate = currentTime;

#ifdef DEBUGV
	Rcpp::Rcout << "\n\n     ********* 5.2: call to ti_nextTime ******\n For to_update = \n " 
		  << "     tSample  = " << tSample
	    
		  << "\n\n**   Species  = " << u_1 
		  << "\n       genotype =  " << Genotypes[u_1] 
		  << "\n       popSize = " << popParams[u_1].popSize 
		  << "\n       currentTime = " << currentTime 
		  << "\n       popParams[i].nextMutationTime = " 
		  << tmpdouble1
		  << " \n     species R " << popParams[u_1].R
		  << " \n     species W " << popParams[u_1].W
		  << " \n     species death " << popParams[u_1].death
		  << " \n     species birth " << popParams[u_1].birth;
#endif

      } else if(to_update == 2) {
	// we did not sample in previous period.
	tmpdouble1 = ti_nextTime_tmax_2_st(popParams[u_1],
					   currentTime,
					   tSample, ti_dbl_min, ti_e3);
	mapTimes_updateP(mapTimes, popParams, u_1, tmpdouble1);
	tmpdouble2 = ti_nextTime_tmax_2_st(popParams[u_2],
					   currentTime,
					   tSample, ti_dbl_min, ti_e3);
	mapTimes_updateP(mapTimes, popParams, u_2, tmpdouble2);
	popParams[u_1].timeLastUpdate = currentTime;
	popParams[u_2].timeLastUpdate = currentTime;

#ifdef DEBUGV
	Rcpp::Rcout << "\n\n     ********* 5.2: call to ti_nextTime ******\n " 
		  << "     tSample  = " << tSample
	    
		  << "\n\n**   Species  = " << u_1 
		  << "\n       genotype =  " << Genotypes[u_1] 
		  << "\n       popSize = " << popParams[u_1].popSize 
		  << "\n       currentTime = " << currentTime 
		  << "\n       popParams[i].nextMutationTime = " 
		  << tmpdouble1
		  << " \n     species R " << popParams[u_1].R
		  << " \n     species W " << popParams[u_1].W
		  << " \n     species death " << popParams[u_1].death
		  << " \n     species birth " << popParams[u_1].birth


		  << "\n\n**     Species  = " << u_2 
		  << "\n       genotype =  " << Genotypes[u_2] 
		  << "\n       popSize = " << popParams[u_2].popSize 
		  << "\n       currentTime = " << currentTime 
		  << "\n       popParams[i].nextMutationTime = " 
		  << tmpdouble2
		  << " \n     species R " << popParams[u_2].R
		  << " \n     species W " << popParams[u_2].W
		  << " \n     species death " << popParams[u_2].death
		  << " \n     species birth " << popParams[u_2].birth;

#endif

      } else { // we sampled, so update all
	for(size_t i = 0; i < popParams.size(); i++) {
	  tmpdouble1 = ti_nextTime_tmax_2_st(popParams[i],
					     currentTime,
					     tSample, ti_dbl_min, ti_e3);
	  mapTimes_updateP(mapTimes, popParams, i, tmpdouble1);
	  popParams[i].timeLastUpdate = currentTime;
	  
#ifdef DEBUGV
	  Rcpp::Rcout << "\n\n     ********* 5.2: call to ti_nextTime ******\n " 
		    << "     Species  = " << i 
		    << "\n       genotype =  " << Genotypes[i] 
		    << "\n       popSize = " << popParams[i].popSize 
		    << "\n       currentTime = " << currentTime 
	    // << "\n       popParams[i].nextMutationTime = " 
	    // << popParams[i].nextMutationTime
		    << " \n     species R " << popParams[i].R
		    << " \n     species W " << popParams[i].W
		    << " \n     species death " << popParams[i].death
		    << " \n     species birth " << popParams[i].birth;
#endif
	}
      }
    }
    if(forceSample) {
      // A VERY ugly hack. Resetting tSample to jump to sampling.
      tSample = currentTime;
      // Need this, o.w. would skip a sampling.
      timeNextPopSample = currentTime;
    } 

    
    // ******************** 5.3 and do we sample? *********** 
    // Find minimum to know if we need to sample the whole pop
    getMinNextMutationTime4(nextMutant, minNextMutationTime, 
			    mapTimes);
    
    if(verbosity >= 2) {
      Rcpp::Rcout << "\n\n  iteration " << iter << "; minNextMutationTime = "
		<< minNextMutationTime 
		<< "; timeNextPopSample = " << timeNextPopSample 
		<< "; popParams.size() = " << popParams.size() << "\n";
    }
    
    // Do we need to sample the population?
    if( minNextMutationTime <= tSample ) {// We are not sampling
      // ************   5.3   **************
      currentTime = minNextMutationTime;

      // ************   5.4   ***************
      mutantTimeSinceLastUpdate = currentTime - 
	popParams[nextMutant].timeLastUpdate;
      
      popParams[nextMutant].popSize = Algo3_st(popParams[nextMutant],
					       mutantTimeSinceLastUpdate);



      
      if(popParams[nextMutant].popSize > (ratioForce * detectionSize)) {
	forceSample = true;
	ratioForce = std::min(1.0, 2 * ratioForce);
#ifdef DEBUGV
	//if(verbosity > -2) {
	// We always warn about this, since interaction with ti==0
	Rcpp::Rcout << "\n Forced sampling triggered for next loop: \n    " << 
	  " popParams[nextMutant].popSize = " << 
	  popParams[nextMutant].popSize << " > ratioForce * detectionSize \n";
	Rcpp::Rcout << " when nextMutant = " << nextMutant <<
	  " at iteration " << iter << "\n";
	//}
#endif
      }      
      // Check also for numSpecies, and force sampling if needed
      // This is very different from the other algos, as we do not yet 
      // now total number of different species
      if(! (numSpecies % speciesFS )) {
      	forceSample = true;
	speciesFS *= 2;
#ifdef DEBUGV 
      	//if(verbosity > -2) // we always warn about this
 
	Rcpp::Rcout << "\n Forced sampling triggered for next loop "
		  << " when numSpecies = " << 
	  numSpecies << " at iteration " << iter << "\n";
#endif
      }
      // Why are these lines here instead of somewhere else?
      // Right before the if for sampling or not?
      // FIXME
      // runningWallTime = difftime(time(NULL), start_time);
      // if( runningWallTime > maxWallTime ) {
      // 	hittedWalllTime = true;
      // 	forceSample = true;
      // 	simulsDone = true;
      // }


      if(popParams[nextMutant].numMutablePos != 0) {
	// this is the usual case. The alternative is the dummy or null mutation

      
	// ************   5.5   ***************

	newMutations.clear();
	obtainMutations(Genotypes[nextMutant],
			fitnessEffects,
			numMutablePosParent,
			newMutations,
			ran_gen);
	// nr_change
	// getMutatedPos_bitset(mutatedPos, numMutablePosParent, // r,
	// 		     ran_gen,
	// 		     mutablePos,
	// 		     Genotypes[nextMutant], 
	// 		     numGenes);
      
	// ************   5.6   ***************

	newGenotype = createNewGenotype(Genotypes[nextMutant],
					newMutations,
					fitnessEffects,
					ran_gen);

	// nr_change
	// newGenotype = Genotypes[nextMutant];
	// newGenotype.set(mutatedPos);
	// newGenotype[mutatedPos] = 1;

	// FIXME
	// any speed diff between a) and b)?
	// a)
	new_sp_v(sp, newGenotype, Genotypes);
	// b)
	// sp = 0;
	// sp = new_sp(newGenotype, Genotypes);
	
	// nr_change
	// new_sp_bitset(sp, newGenotype, Genotypes);

	if(sp == numSpecies) {// New species
	  ++numSpecies;
	  init_tmpP(tmpParam);

	  if(verbosity >= 2) {
	    Rcpp::Rcout <<"\n     Creating new species   " << (numSpecies - 1)
			<< "         from species "  <<   nextMutant;
	  }
	
	  tmpParam.popSize = 1;

	  nr_fitness(tmpParam, popParams[nextMutant],
		     newGenotype,
		     fitnessEffects,
		     typeFitness, genTime,
		     adjust_fitness_B, adjust_fitness_MF);
	
	  if(tmpParam.birth > 0.0) {
	    tmpParam.numMutablePos = numMutablePosParent - 1;
	    if(mutatorGenotype)
	      tmpParam.mutation = mu * tmpParam.birth * tmpParam.numMutablePos;
	    //	    tmpParam.mutation = mu * tmpParam.birth * (numMutablePosParent - 1);
	    else
	      tmpParam.mutation = mu * tmpParam.numMutablePos;
	    //tmpParam.mutation = mu * (numMutablePosParent - 1);
	    if (tmpParam.mutation > 1 )
	      Rcpp::Rcout << "WARNING: mutation > 1\n";
	    if (numMutablePosParent == 1) {
	      Rcpp::Rcout << "Note: mutation = 0; no positions left for mutation\n";
	      tmpParam.mutation = dummyMutationRate; // dummy mutation here. Set some mu.
	    }
	    W_f_st(tmpParam);
	    R_f_st(tmpParam);
	    tmpParam.timeLastUpdate = -99999.99999; //mapTimes_updateP does what it should.
	    popParams.push_back(tmpParam);
	    Genotypes.push_back(newGenotype);
	    to_update = 2;
#ifdef MIN_RATIO_MUTS
	    g_tmp1_nr = tmpParam.birth/tmpParam.mutation;
	    if(g_tmp1_nr < g_min_birth_mut_ratio_nr) g_min_birth_mut_ratio_nr = g_tmp1_nr;
	  
	    g_tmp1_nr = tmpParam.death/tmpParam.mutation;
	    if(g_tmp1_nr < g_min_death_mut_ratio_nr) g_min_death_mut_ratio_nr = g_tmp1_nr;	
#endif	  
	  } else {// fitness is 0, so we do not add it
	    --sp;
	    --numSpecies;
	    to_update = 1;
	  }
	  // #ifdef DEBUGV	
	  if(verbosity >= 3) {
	    Rcpp::Rcout << " \n\n\n Looking at NEW species " << sp << " at creation";
	    Rcpp::Rcout << "\n New Genotype :";
	    print_Genotype(newGenotype);
	    Rcpp::Rcout << "\n Parent Genotype :";
	    print_Genotype(Genotypes[nextMutant]);
	    // Rcpp::Rcout << "\n Genotype = " << genotypeSingleVector(newGenotype); //Genotypes[sp];
	    //Genotypes[sp].to_ullong();
	    Rcpp::Rcout << "\n birth of sp = " << tmpParam.birth;
	    Rcpp::Rcout << "\n death of sp = " << tmpParam.death;
	    // Rcpp::Rcout << "\n s = " << s;
	    Rcpp::Rcout << "\n parent birth = " << popParams[nextMutant].birth;
	    Rcpp::Rcout << "\n parent death = " << popParams[nextMutant].death;
	    // Rcpp::Rcout << "\n parent Genotype = " << genotypeSingleVector(Genotypes[nextMutant]);
	    print_spP(tmpParam);
	    }
	  // #endif
	} else {	// A mutation to pre-existing species
#ifdef DEBUGW
	  if( (currentTime - popParams[sp].timeLastUpdate) <= 0.0)
	    throw std::out_of_range("currentTime - timeLastUpdate out of range"); 
#endif
	
	  // if(verbosity >= 2) {
#ifdef DEBUGV
	    Rcpp::Rcout <<"\n     Mutated to existing species " << sp 
			<< " (Genotype = " << genotypeSingleVector(Genotypes[sp]) 
	      // << "; sp_id = " << Genotypes[sp].to_ullong()
			<< ")"
			<< "\n from species "  <<   nextMutant
			<< " (Genotypes = " << genotypeSingleVector(Genotypes[nextMutant]) 
	      // << "; sp_id = " << Genotypes[sp].to_ullong()
			<< ")";
	    // }
#endif
	  // FIXME00: the if can be removed??
	  if(popParams[sp].popSize > 0.0) {
	    popParams[sp].popSize = 1.0 + 
	      Algo2_st(popParams[sp], currentTime);
	    if(verbosity >= 2) {
	      Rcpp::Rcout << "\n New popSize = " << popParams[sp].popSize << "\n";
	    }
	  } else {
	    throw std::range_error("\n popSize == 0 but existing? \n");
	  }
	
#ifdef DEBUGW
	  popParams[sp].timeLastUpdate = -99999.99999; // to catch errors
#endif
	  //popParams[sp].Flag = true;
	}
	//   ***************  5.7 ***************
	// u_2 irrelevant if to_update = 1;
	u_1 = nextMutant;
	u_2 = static_cast<int>(sp);
      } else { // the null or dummy mutation case
	// Rcpp::Rcout << "\n null mutation; before popSize" << std::endl;
	// DP2(popParams[nextMutant].popSize);
	++popParams[nextMutant].popSize;
	to_update = 1;
	u_1 = nextMutant;
	u_2 = -99;
	// FIXME: do this conditionally on flag
	Rcpp::Rcout << "Note: updating in null mutation\n";
	// Rcpp::Rcout << "\n null mutation; after popSize" << std::endl;
	// DP2(popParams[nextMutant].popSize);
	// Rcpp::Rcout << "\n done null mutation; after popSize ********" << std::endl;
      }
    }
      else { //       *********** We are sampling **********
      to_update = 3; //short_update = false;
      if(verbosity >= 2) {
	Rcpp::Rcout <<"\n We are SAMPLING";   
	if(tSample < finalTime) {
	  Rcpp::Rcout << " at time " << tSample << "\n";
	} else
	  Rcpp::Rcout <<". We reached finalTime " << finalTime << "\n";
      }

      currentTime = tSample;
      if(verbosity >= 3)
	Rcpp::Rcout << "\n popParams.size() before sampling " << popParams.size() << "\n";

      nr_sample_all_pop_P(sp_to_remove, 
		       popParams, Genotypes, tSample);
      timeNextPopSample += sampleEvery;
      
      if(sp_to_remove.size())
	remove_zero_sp_nr(sp_to_remove, Genotypes, popParams, mapTimes);

      numSpecies = popParams.size();
      
      nr_totPopSize_and_fill_out_crude_P(outNS_i, totPopSize, 
				      lastStoredSample,
				      genot_out, 
				      //sp_id_out,
				      popSizes_out, index_out,
				      time_out, 
				      sampleTotPopSize,sampleLargestPopSize,
				      sampleMaxNDr, sampleNDrLargestPop,
				      simulsDone,
				      reachDetection,
				      lastMaxDr,
				      done_at,
				      Genotypes, popParams, 
				      currentTime,
				      keepEvery,
				      detectionSize,
				      finalTime,
				      //endTimeEvery,
				      detectionDrivers,
				      verbosity,
					 minDDrPopSize,
					 extraTime,
					 fitnessEffects.drv); //keepEvery is for thinning
      if(verbosity >= 3) {
	Rcpp::Rcout << "\n popParams.size() before sampling " << popParams.size() 
		  << "\n totPopSize after sampling " << totPopSize << "\n";
      }
      
      computeMcFarlandError(e1, n_0, n_1, tps_0, tps_1, 
			    typeFitness, totPopSize, K); //, initSize);

      // Largest error in McFarlands' method
      // if( (typeFitness == "mcfarland0") ||
      // 	  (typeFitness == "mcfarland") || 
      // 	  (typeFitness == "mcfarlandlog") ) {
      // 	tps_1 = totPopSize;
      // 	if(typeFitness == "mcfarland")
      // 	  etmp = abs( tps_1 - (tps_0 + 1) );
      // 	else {
      // 	  if( (tps_0 + 1.0) > tps_1 ) 
      // 	    etmp = (K + tps_0 + 1.0)/(K + tps_1);
      // 	  else
      // 	    etmp = (K + tps_1)/(K + tps_0 + 1);
      // 	}
      // 	if(etmp > e1) {
      // 	  e1 = etmp;
      // 	  n_0 = tps_0;
      // 	  n_1 = tps_1;
      // 	}
      // 	tps_0 = tps_1;
      // }

      // It goes here: zz: not detectionSize,
      // but the keepEvery? or sampleUntilKeep. Yes, use that.
      // endingSampleEvery.

      // Use driver criterion here!!! 
      // if endingSampleEvery
      // if totPopSize >= detectionSize: 
      //        do not break unless 
      //        tSample %% endingSampleEvery 

      // All of this has to be in totPopSize_and_fill


      // this if not sampleUntilKeep
      // if( (totPopSize >= detectionSize) ||
      // 	  (totPopSize <= 0.0) || (tSample >= finalTime)) {	
      // 	simulsDone = true;
      // 	break; // skip last update if beerenwinkel
      // }       
      
      if(simulsDone)
	break; //skip last updateRates

      if( (typeFitness == "beerenwinkel") ) {
	updateRatesBeeren(popParams, adjust_fitness_B,
			  initSize, currentTime, alpha, totPopSize,
			  mutatorGenotype, mu);
      } else if( (typeFitness == "mcfarland0") ) {
	updateRatesMcFarland0(popParams, adjust_fitness_MF,
			     K, totPopSize,
			     mutatorGenotype, mu);
      } else if( (typeFitness == "mcfarland") ) {
	updateRatesMcFarland(popParams, adjust_fitness_MF,
			     K, totPopSize);
      } else if( (typeFitness == "mcfarlandlog") ) {
	updateRatesMcFarlandLog(popParams, adjust_fitness_MF,
			     K, totPopSize);
      }
      
#ifdef MIN_RATIO_MUTS
      // could go inside sample_all_pop but here we are sure death, etc, current
      // But I catch them when they are created. Is this really needed?
      for(size_t i = 0; i < popParams.size(); i++) {
	g_tmp1_nr = popParams[i].birth/popParams[i].mutation;
	if(g_tmp1_nr < g_min_birth_mut_ratio_nr) g_min_birth_mut_ratio_nr = g_tmp1_nr;
	
	g_tmp1_nr = popParams[i].death/popParams[i].mutation;
	if(g_tmp1_nr < g_min_death_mut_ratio_nr) g_min_death_mut_ratio_nr = g_tmp1_nr;
      }
#endif
      
      forceSample = false;
    }
  }
}


// Start being explicit about parameter types


// Rcpp::List nr_BNB_Algo5(SEXP restrictTable_,
// 		     SEXP numDrivers_,
// 		     SEXP numGenes_,
// 		     SEXP typeCBN_,
// 		     SEXP birthRate_, 
// 		     SEXP s_, 
// 		     SEXP death_,
// 		     SEXP mu_,
// 		     SEXP initSize_,
// 		     SEXP sampleEvery_,
// 		     SEXP detectionSize_,
// 		     SEXP finalTime_,
// 		     SEXP initSize_species_,
// 		     SEXP initSize_iter_,
// 		     SEXP seed_gsl_,
// 		     SEXP verbose_,
// 		     SEXP speciesFS_,
// 		     SEXP ratioForce_,
// 		     SEXP typeFitness_,
// 		     SEXP maxram_,
// 		     SEXP mutatorGenotype_,
// 		     SEXP initMutant_,
// 		     SEXP maxWallTime_,
// 		     SEXP keepEvery_,
// 		     SEXP alpha_,
// 		     SEXP sh_,
// 		     SEXP K_,
// 		     SEXP detectionDrivers_,
// 		     SEXP onlyCancer_,
// 		     SEXP errorHitWallTime_,
// 		     SEXP maxNumTries_,
// 		     SEXP errorHitMaxTries_,
// 		     SEXP minDDrPopSize_,
// 		     SEXP extraTime_
// 		     ) {

// [[Rcpp::export]]
Rcpp::List nr_BNB_Algo5(Rcpp::List rFE,
			double mu,
			double death,
			double initSize,
			double sampleEvery,
			double detectionSize,
			double finalTime,
			int initSp,
			int initIt,
			int seed,
			int verbosity,
			int speciesFS,
			double ratioForce,
			Rcpp::CharacterVector typeFitness_,
			int maxram,
			int mutatorGenotype,
			Rcpp::IntegerVector initMutant_, 
			double maxWallTime,
			double keepEvery,
			double alpha,
			double K,
			int detectionDrivers,
			bool onlyCancer,
			bool errorHitWallTime,
			int maxNumTries,
			bool errorHitMaxTries,
			double minDDrPopSize,
			double extraTime) {  
  // SEXP endTimeEvery_,


  
  //  BEGIN_RCPP
  // using namespace Rcpp;
  precissionLoss();
  const std::vector<int> initMutant = Rcpp::as<std::vector<int> >(initMutant_);
  // const std::string typeFitness = as<std::string>(typeFitness_);
  const std::string typeFitness = Rcpp::as<std::string>(typeFitness_); // no need to do [0]
  
  // birth and death are irrelevant with Bozic
  // const double death = as<double>(death_);
  // const double mu = as<double>(mu_);
  // const double initSize = as<double>(initSize_);
  // const double sampleEvery = as<double>(sampleEvery_);
  // const double detectionSize = as<double>(detectionSize_);
  // const double finalTime = as<double>(finalTime_);
  // const int initSp = as<int>(initSize_species_);
  // const int initIt = as<int>(initSize_iter_); // FIXME: this is a misnomer
  // const int verbosity = as<int>(verbose_);
  // // const double minNonZeroMut = mu * 0.01;  // to avoid == 0.0 comparisons
  // double ratioForce = as<double>(ratioForce_); // If a single species this times
  // // detectionSize, force a sampling to prevent going too far.
  // int speciesFS = as<int>(speciesFS_); // to force sampling when too many 
  // // species
  // const int seed = as<int>(seed_);
  // const long maxram = as<int>(maxram_);
  // const int mutatorGenotype = as<int>(mutatorGenotype_);
  // const int initMutant = as<int>(initMutant_);
  // const double maxWallTime = as<double>(maxWallTime_);
  // const double keepEvery = as<double>(keepEvery_);

  
  // const double alpha = as<double>(alpha_);
  // // if a driver without dependencies. Like in Datta et al., 2013.
  // const double K = as<double>(K_); //for McFarland
  // //const double endTimeEvery = as<double>(endTimeEvery_); 
  // const int detectionDrivers = as<int>(detectionDrivers_); 
  const double genTime = 4.0; // should be a parameter. For Bozic only.
  // const bool errorFinalTime = as<bool>(errorFinalTime_);
  // const bool errorHitWallTime = as<bool>(errorHitWallTime_);
  // const bool onlyCancer = as<bool>(onlyCancer_);
  // const int maxNumTries = as<int>(maxNumTries_);
  // const bool errorHitMaxTries = as<bool>(errorHitMaxTries_);
  // const double minDDrPopSize = as<double>(minDDrPopSize_);
  // const double extraTime = as<double>(extraTime_);
  
  std::mt19937 ran_gen(seed);
  // DP2(seed);
  
  // some checks. Do this systematically
  // FIXME: do only if mcfarland!
  if(K < 1 )
    throw std::range_error("K < 1.");
  fitnessEffectsAll fitnessEffects =  convertFitnessEffects(rFE);
  
  bool runAgain = true;
  bool reachDetection = false;
  //Output
  std::vector<Genotype> genot_out;
  std::vector<double> popSizes_out;
  std::vector<int> index_out; 
  std::vector<double> time_out; //only one entry per period!
  genot_out.reserve(initSp);
  popSizes_out.reserve(initSp);
  index_out.reserve(initSp);
  time_out.reserve(initIt);

  double totPopSize = 0;
  std::vector<double> sampleTotPopSize;
  std::vector<double> sampleLargestPopSize;
  std::vector<int> sampleMaxNDr; //The number of drivers in the population
  // with the largest number of drivers; and this for each time sample
  std::vector<int> sampleNDrLargestPop; //Number of drivers in population
  // with largest size (at each time sample)
  sampleTotPopSize.reserve(initIt);
  sampleLargestPopSize.reserve(initIt);
  sampleMaxNDr.reserve(initIt);
  sampleNDrLargestPop.reserve(initIt);


  int outNS_i = -1; // the column in the outNS
  // time limits
  // FIXME think later FIXME
  time_t start_time = time(NULL);
  double runningWallTime = 0;
  bool  hittedWallTime = false;
  bool hittedMaxTries = false;
 
  // spParamsP tmpParam; 
  // std::vector<spParamsP> popParams(1);
  // const int sp_per_period = 5000;

  // popParams.reserve(sp_per_period);
  // Genotypes.reserve(sp_per_period);

  // std::vector<int>mutablePos(numGenes); // could be inside getMuatedPos_bitset


  // // multimap to hold nextMutationTime
  // std::multimap<double, int> mapTimes;
  // //std::multimap<double, int>::iterator m1pos;


  // // count troublesome tis
  // int ti_dbl_min = 0;
  // int ti_e3 = 0;


  
  // // Beerenwinkel
  // double adjust_fitness_B = -std::numeric_limits<double>::infinity();
  // //McFarland
  // double adjust_fitness_MF = -std::numeric_limits<double>::infinity();

  double e1, n_0, n_1; // for McFarland error
  // double tps_0, tps_1; // for McFarland error
  // tps_0 = 0.0;
  // tps_1 = 0.0;
  e1 = 0.0;
  n_0 = 0.0;
  n_1 = 0.0;

  // // For totPopSize_and_fill and bailing out
  // // should be static vars inside funct,
  // // but they keep value over calls in same R session.
  // int lastMaxDr = 0;
  // double done_at = -9;
  // // totalPopSize at time t, at t-1 and the max error.

  // 5.1 Initialize 

  int numRuns = 0;
  bool forceRerun = false;
  
  double currentTime = 0;
  int iter = 0;
  while(runAgain) {

    // Initialize a bunch of things
    
// #ifdef MIN_RATIO_MUTS
//   g_min_birth_mut_ratio = DBL_MAX;
//   g_min_death_mut_ratio = DBL_MAX;
//   g_tmp = DBL_MAX;
// #endif

  
  // untilcancer goes here
  

  
  
  //tmpParam is a temporary holder. 
  // init_tmpP(tmpParam);
  // init_tmpP(popParams[0]);

  // lastStoredSample = 0.0;
  // Genotypes[0].reset();
  // popParams[0].popSize = initSize;
  // totPopSize = initSize;

  // tps_0 = totPopSize;
  // e1 = 0.0;
  // tps_1 = totPopSize;



    try {
      // it is CRUCIAL that several entries are zeroed (or -1) at the
      // start of innerBNB now that we do multiple runs if onlyCancer = true.
      nr_innerBNB(
		  fitnessEffects,
		  initSize,
	       K,
	       alpha,
	       genTime,
	       typeFitness,
	       mutatorGenotype,
	       mu,
	       death,
	       keepEvery,
	       sampleEvery,		     
	       initMutant,
	       start_time,
	       maxWallTime,
	       finalTime,
	       detectionSize,
	       detectionDrivers,
	       minDDrPopSize,
	       extraTime,
	       verbosity,
	       totPopSize,
	       e1,
	       n_0,
	       n_1,
	       ratioForce,
	       currentTime,
	       speciesFS,
	       outNS_i, 
	       iter,
	       genot_out,
	       popSizes_out,
	       index_out,
	       time_out,
	       sampleTotPopSize,
	       sampleLargestPopSize,
	       sampleMaxNDr,
	       sampleNDrLargestPop,
	       reachDetection,
	       ran_gen,
	       runningWallTime,
	       hittedWallTime);
      ++numRuns;
      forceRerun = false;
    } catch (rerunExcept &e) {
      Rcpp::Rcout << "\n Exception " << e.what() 
		  << ". Rerunning.";
      forceRerun = true;
    } catch (const std::exception &e) {
      Rcpp::Rcout << "\n Unrecoverable exception: " << e.what()
		  << ". Aborting. \n";
      return
	List::create(Named("other") =
		     List::create(Named("UnrecoverExcept") = true,
				  Named("ExceptionMessage") = e.what()));
    } catch (...) {
      Rcpp::Rcout << "\n Unknown unrecoverable exception. Aborting. \n";
      return
	List::create(Named("other") =
		     List::create(Named("UnrecoverExcept") = true,
				  Named("ExceptionMessage") = "Unknown exception"));
    }


    if(hittedWallTime) {
      Rcpp::Rcout << "\n Hitted wall time. Exiting.";
      runAgain = false;
      if(errorHitWallTime) {
	Rcpp::Rcout << "\n Hitting wall time is regarded as an error. \n";
	return
	  List::create(Named("HittedWallTime") = true,
		       Named("HittedMaxTries") = false, // yes, for
							// coherent return
							// objects
		       Named("other") =
		       List::create(Named("UnrecoverExcept") = false));
      }
    } else if(numRuns > maxNumTries) {
      //  hittedMaxTries
      hittedMaxTries = true;
      Rcpp::Rcout << "\n Hitted maxtries. Exiting.";
      runAgain = false;
      if(errorHitMaxTries) {
	Rcpp::Rcout << "\n Hitting max tries is regarded as an error. \n";
	return
	  List::create(Named("HittedWallTime") = false,
		       Named("HittedMaxTries") = true,
		       Named("other") =
		       List::create(Named("UnrecoverExcept") = false));
      }
    } else if(forceRerun) {
      runAgain = true;
      forceRerun = false;
    } else {
      if(onlyCancer) {
	runAgain = !reachDetection;
      } else {
	runAgain = false;
      }
    }
#ifdef DEBUGV
      Rcpp::Rcout << "\n reachDetection = " << reachDetection;
      Rcpp::Rcout << "\n forceRerun =  " << forceRerun  << "\n";
      
#endif
    
  } // runAgain loop
  // FIXME: zz
  // untilcancer
  // inner loop ends above
  // The return objects only created if needed

  
  // If we hit wallTime, we can get done without going through
  // totPopSize.... Problem if sampling at end
  // if ( hittedWallTime ) {
  //   // hitted wall time. So we need to sample at the very end.
  // Nope! Just ensure if hittedWallTime you always sample properly!
  // }
    

  
  // FIXME: do I want to move this right after out_crude
  // and do it incrementally? I'd have also a counter of total unique species

  // here("right after simuls done");

  // FIXME: all this is ugly and could be a single function
  // up to call to IntegerMatrix

  // Do I use genot_out_nr for anything else? FIXME
  // Put into a single function?? FIXME

  // // Go from genot_out (the vector of Genotype) to a vector of vectors of
  // // the genotypes, through a set.

  // // V1
  // std::vector<std::vector<int> > genot_out_nr;
  // std::transform(genout_out.begin(), genot_out.end(), back_inserter(genout_out_nr),
  // 		 genotypeSingleVector);
  // std::set<std::vector<int> > uniqueGenotypes_nr
  // for( auto gg : genot_out_nr ) uniqueGenotypes_nr.insert(gg);
  // std::vector<std::vector<int> > uniqueGenotypes_vector_nr (uniqueGenotypes_nr.begin(),
  // 							    uniqueGenotypes_nr.end());


  // // v2
  // std::vector<std::vector<int> > genot_out_nr = genot_to_vectorg(genot_out);
  // std::set<std::vector<int> > uniqueGenotypes_nr =  nr_find_unique_genotypes(genot_out_nr);
  // std::vector<std::vector<int> > uniqueGenotypes_vector_nr = nr_uniqueGenotypes_to_vector(uniqueGenotypes_nr);


  // v3
  // Need the two below
  std::vector<std::vector<int> > genot_out_v = genot_to_vectorg(genot_out);
  std::vector<std::vector<int> > uniqueGenotypes_vector_nr  =
    uniqueGenot_vector(genot_out_v);
  
  // IntegerMatrix returnGenotypes(uniqueGenotypes_vector.size(), numGenes);
  IntegerMatrix returnGenotypes = 
    nr_create_returnGenotypes(fitnessEffects.genomeSize,
			      uniqueGenotypes_vector_nr);
  
  Rcpp::NumericMatrix outNS = create_outNS(uniqueGenotypes_vector_nr,
					   genot_out_v,
					   popSizes_out,
					   index_out, time_out,
					   outNS_i, maxram);
  

  int maxNumDrivers = 0;
  int totalPresentDrivers = 0;
  std::vector<int>countByDriver(fitnessEffects.drv.size(), 0);
  std::string occurringDrivers;

  nr_count_NumDrivers(maxNumDrivers, countByDriver,
		      returnGenotypes, fitnessEffects.drv);

  whichDrivers(totalPresentDrivers, occurringDrivers, countByDriver);

  std::vector<double> sampleLargestPopProp(outNS_i + 1);

  if((outNS_i + 1) != static_cast<int>(sampleLargestPopSize.size()))
    throw std::length_error("outNS_i + 1 != sampleLargestPopSize.size");
  std::transform(sampleLargestPopSize.begin(), sampleLargestPopSize.end(),
		 sampleTotPopSize.begin(),
		 sampleLargestPopProp.begin(),
		 std::divides<double>());

  NumericMatrix perSampleStats(outNS_i + 1, 5);
  fill_SStats(perSampleStats, sampleTotPopSize, sampleLargestPopSize,
	      sampleLargestPopProp, sampleMaxNDr, sampleNDrLargestPop);

  std::vector<std::string> genotypesLabels =
    genotypesToString(uniqueGenotypes_vector_nr, fitnessEffects);
  // error in mcfarland's
  // if((typeFitness == "mcfarland0") || (typeFitness == "mcfarlandlog"))
  //   e1r = log(e1);
  // if(typeFitness == "mcfarland")
  //   e1r = (1.0/K) * e1;

  // here("before return");

  // // // debuggin: precompute things
  // DP2(simulsDone);
  // DP2(maxWallTime);
  // DP2(hittedWallTime);
  // DP2(outNS_i);
  // DP2( sampleMaxNDr[outNS_i]);
  // DP2(sampleNDrLargestPop[outNS_i]);
  // DP2(sampleLargestPopSize[outNS_i]);
  // DP2(sampleLargestPopProp[outNS_i]);
  // DP2((runningWallTime > maxWallTime));
  // here("after precomp");
  // here("*******************************************");

  // Rcpp::List returnGenotypesO = Rcpp::wrap(uniqueGenotypesV);
  
  return 
    List::create(Named("pops.by.time") = outNS,
		 Named("NumClones") = uniqueGenotypes_vector_nr.size(), 
		 Named("TotalPopSize") = totPopSize,
		 Named("Genotypes") = returnGenotypes,
		 Named("GenotypesWDistinctOrderEff") = Rcpp::wrap(uniqueGenotypes_vector_nr),
		 Named("GenotypesLabels") = Rcpp::wrap(genotypesLabels),
		 Named("MaxNumDrivers") = maxNumDrivers,
		 // Named("MaxDrivers_PerSample") = wrap(sampleMaxNDr),
		 // Named("NumDriversLargestPop_PerSample") = sampleNDrLargestPop,
		 // Named("TotPopSize_PerSample") = sampleTotPopSize,
		 // Named("LargestPopSize_PerSample") = sampleLargestPopSize,
		 // Named("PropLargestPopSize_PerSample") = sampleLargestPopProp,
		 Named("MaxDriversLast") = sampleMaxNDr[outNS_i],
		 Named("NumDriversLargestPop") =  sampleNDrLargestPop[outNS_i],
		 Named("LargestClone") = sampleLargestPopSize[outNS_i],
		 Named("PropLargestPopLast") = sampleLargestPopProp[outNS_i],
		 // Named("totDrivers") = totDrivers,
		 Named("FinalTime") = currentTime,
		 Named("NumIter") = iter,
		 //		 Named("outi") = outNS_i + 1, // silly. Use the real number of samples. FIXME
		 Named("HittedWallTime") = hittedWallTime, // (runningWallTime > maxWallTime),
		 Named("HittedMaxTries") = hittedMaxTries,
		 // Named("iRunningWallTime") = runningWallTime,
		 // Named("oRunningWallTime") = difftime(time(NULL), start_time),
		 // Named("ti_dbl_min") = ti_dbl_min,
		 // Named("ti_e3") = ti_e3,
		 Named("TotalPresentDrivers") = totalPresentDrivers,
		 Named("CountByDriver") = countByDriver,
		 // FIXME: OccurringDrivers underestimates true occurring
		 // drivers if keepEvery < 0, so we only return the last.
		 Named("OccurringDrivers") = occurringDrivers,
		 Named("PerSampleStats") = perSampleStats,
		 Named("other") = List::create(Named("attemptsUsed") = numRuns,
					       Named("errorMF") = 
					       returnMFE(e1, K, 
							 typeFitness),
					       Named("errorMF_size") = e1,
					       Named("errorMF_n_0") = n_0,
#ifdef MIN_RATIO_MUTS
					       Named("minDMratio") =
					       g_min_death_mut_ratio_nr,
					       Named("minBMratio") =
					       g_min_birth_mut_ratio_nr,      
#else
					       Named("minDMratio") = -99,
					       Named("minBMratio") = -99,
#endif
					       Named("errorMF_n_1") = n_1,
					       Named("UnrecoverExcept") = false)
		 );

  //  END_RCPP
    
    }






// void Mutation(const double mu, const double chromothripsis) {
//   if(chromothripsis) {
//     do something;
//   } else {
//     pointMutation();
//   }

// }



// FIXME: get into dealing with fitness!




// Creating return object:


// The 0, 1 representation is how most of the work is done in R: do I want
// to change that?


// Order: beware of two things: order is important for the "true"
// genotypes, but is not immediately observable. So for 0,1
// representation, not needed or used. Thus, maybe I want two
// representations.


// I do not need all the ulong thing: I can always work directly with the
// Genotype struct, everywhere. And then return, to R either the 0,1 or
// the vector of mutated positions, or both.

// Yes, the full Genotye structure is only used when assigning fitness. So
// could we use a collapsed one with: order + rest? Nope, as whenever I'd
// create a child from a genotype, I'd need like deconvolve, and go back
// to the three piece structure. This seems much more expensive than the
// overloaded == and the usage of the overloaded < (this is only used at
// the end, when producing the output objects)
