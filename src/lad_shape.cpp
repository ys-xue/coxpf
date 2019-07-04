//' Compute a reference shape curve for a collection of profiles
//' 
//' @name lad_shape
//' @usage lad_shape(profileMat, bandwidth, gridprecision, timevec)
//' @param profileMat a matrix of profiles
//' @param bandwidth the bandwidth used in RBF kernel
//' @param gridprecision the precision for grid search
//' @param timevec vector of time
//' @return a vector whose length equals the length of time vector containing
//' the reference variability at each time point
//' @export


#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec lad_shape(const arma::mat& profileMat,
                    const float bandwidth,
                    const int gridprecision,
                    const arma::vec& timevec){
    unsigned int K = profileMat.n_cols;
    

    
    // Rcpp::Rcout << pfMedians << std::endl;
    
    arma::vec thetaVec = arma::linspace(0,0,K);
    arma::vec thetaVec2 = arma::linspace(0,0,K);
    
    arma::vec tomin = arma::zeros(gridprecision);
    arma::vec tomin2 = arma::zeros(gridprecision);
    arma::vec kernval = arma::zeros(timevec.n_elem);
    arma::vec kernval2 = arma::zeros(timevec.n_elem);
    arma::vec shape = arma::zeros(timevec.n_elem);
    arma::vec shape2 = arma::zeros(timevec.n_elem);
    arma::vec timepoints; 
    
    for (unsigned int k = 0; k < K; k += 1){
        arma::vec searchRange = arma::linspace(profileMat.col(k).min(),
                                               profileMat.col(k).max(),
                                               gridprecision);
        // Rcpp::Rcout << searchRange << std::endl;
        timepoints = (timevec - timevec(k)) / bandwidth;
        
        for (unsigned int j = 0; j < timevec.n_elem; j += 1){
            kernval(j) = R::dnorm(timepoints(j), 0,1,0);
            kernval2(j) = R::dnorm(timepoints(j) / sqrt(2), 0, 1, 0);
        }
        
        kernval /= arma::accu(kernval);
        kernval2 /= arma::accu(kernval2);
        
        for (unsigned int p = 0; p < searchRange.n_elem; p += 1){
            tomin(p) = arma::accu(abs(profileMat - searchRange(p)) * kernval);
            tomin2(p) = arma::accu(abs(profileMat - searchRange(p)) * kernval2);
        }
        
        shape(k) = searchRange(tomin.index_min());
        shape2(k) = searchRange(tomin2.index_min());
        // Rcpp::Rcout << tomin << std::endl;
    }
    arma::vec shapeCorrected = shape * 2 - shape2;
    return shapeCorrected;
}

