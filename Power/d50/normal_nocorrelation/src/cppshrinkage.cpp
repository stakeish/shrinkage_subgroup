#include <RcppArmadillo.h>
#include <math.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
double llogit(arma::vec gamma, arma::mat z, arma::vec w) {
  
  double ll;
  arma::vec zg, logit_vec, llvec;
  NumericVector zg_vec;
  
  zg = z*gamma;
  zg_vec = wrap(zg);
  logit_vec = as<arma::vec>(wrap(Rcpp::plogis(zg_vec)));
  
  llvec = w%log(logit_vec) + (1 - w)%log((1 - logit_vec));
  ll = sum(llvec);
  
  return -ll; 
}

// [[Rcpp::export]]
arma::vec dlogit(arma::vec gamma, arma::mat z, arma::vec w){
  
  arma::vec zg, logit_vec, dl;
  arma::mat z_t;
  NumericVector zg_vec;
  
  zg = z*gamma;
  zg_vec = wrap(zg);
  logit_vec = as<arma::vec>(wrap(Rcpp::plogis(zg_vec)));
  z_t = z.t();
  
  dl = z_t *(w - logit_vec);
  
  return -dl; 
}

// [[Rcpp::export]]
double obj(arma::vec gamma, arma::mat z, arma::mat x, arma::vec y, arma::vec alpha, double beta, double lambda,
           arma::vec d, double sigma){
  
  double ll;
  arma::vec zg, mean1, mean2, logit_vec, norm1_vec, norm2_vec, llvec;
  NumericVector zg_vec, mean1_vec, mean2_vec;
  
  zg = z * gamma;
  mean2 =  y - x*alpha - d*beta;
  mean1 = mean2 - d*lambda;
  
  zg_vec = wrap(zg);
  mean1_vec = wrap(mean1);
  mean2_vec = wrap(mean2);
  
  logit_vec = as<arma::vec>(wrap(Rcpp::plogis(zg_vec)));
  norm1_vec = as<arma::vec>(wrap(Rcpp::dnorm(mean1_vec, 0, sigma)));
  norm2_vec = as<arma::vec>(wrap(Rcpp::dnorm(mean2_vec, 0, sigma)));
  
  llvec = log(logit_vec%norm1_vec + (1-logit_vec)%norm2_vec);
  
  ll = sum(llvec);
  
  return -ll;
}

// [[Rcpp::export]]
arma::vec grad(arma::vec gamma, arma::mat z, arma::mat x, arma::vec y, arma::vec alpha, double beta, double lambda,
           arma::vec d, double sigma){
  
  arma::vec zg, mean1, mean2, logit_vec, dlogit_vec, norm1_vec, norm2_vec, num, den, frac, grad;
  arma::mat z_t; 
  NumericVector zg_vec, mean1_vec, mean2_vec;
  
  zg = z * gamma;
  mean2 =  y - x*alpha - d*beta;
  mean1 = mean2 - d*lambda;
  
  zg_vec = wrap(zg);
  mean1_vec = wrap(mean1);
  mean2_vec = wrap(mean2);
  
  logit_vec = as<arma::vec>(wrap(Rcpp::plogis(zg_vec)));
  dlogit_vec = as<arma::vec>(wrap(Rcpp::dlogis(zg_vec)));
  norm1_vec = as<arma::vec>(wrap(Rcpp::dnorm(mean1_vec, 0, sigma)));
  norm2_vec = as<arma::vec>(wrap(Rcpp::dnorm(mean2_vec, 0, sigma)));
  
  z_t = z.t();
  
  num = logit_vec % norm1_vec + (1 - logit_vec) % norm2_vec;
  den = dlogit_vec%(norm1_vec - norm2_vec);
  frac = den / num;
  
  grad = z_t * frac;
  
  return -grad; 
}
