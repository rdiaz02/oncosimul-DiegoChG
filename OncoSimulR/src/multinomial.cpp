#include <Rcpp.h>
using namespace Rcpp;

// Copyright (C)  2014  Christian Gunning
//
// This file is part of RcppArmadillo.
// [[Rcpp::export]]
IntegerVector rmultinom(int size,  Rcpp::NumericVector prob) {
  // meaning of n, size, prob as in ?rmultinom
  // opposite of sample() - n=number of draws
  double pp;            
  int ii;
  int probsize = prob.size();
  // Return object
  IntegerVector draws(probsize);
  if (size < 0 || size == NA_INTEGER) throw std::range_error( "Invalid size");
  long double p_tot = 0.;
  p_tot = std::accumulate(prob.begin(), prob.end(), p_tot);
  if (fabs((double)(p_tot - 1.)) > 1e-7) {
    throw std::range_error("Probabilities don't sum to 1, please use FixProb");
  }
  
  // do as rbinom
  if (size == 0 ) {
    return draws;
  }
  //rmultinom(size, REAL(prob), k, &INTEGER(ans)[ik]);
  // for each slot
  for (ii = 0; ii < probsize-1; ii++) { /* (p_tot, n) are for "remaining binomial" */
    if (prob[ii]) {
      pp = prob[ii] / p_tot;
      // >= 1; > 1 happens because of rounding 
      draws[ii] = ((pp < 1.) ? (int) Rf_rbinom((double) size,  pp) : size);
      size -= draws[ii];
    } // else { ret[ii] = 0; }
    // all done
    if (size <= 0)  return draws;
    // i.e. p_tot = sum(prob[(k+1):K]) 
    p_tot -= prob[ii]; 
  }
  // the rest go here
  draws[probsize-1] = size;
  return draws;
}


