#include <Rcpp.h>
using namespace Rcpp;
/*-------------------------------------------------------*
 * Function: ZINB_MLE
 * Purpose: Perform maximum likelihood estimation for a Zero-Inflated Negative Binomial (ZINB) distribution using the EM algorithm. 
 *          This function is specifically designed to estimate parameters for datasets with excess zeros.
 * Input:
 *   X (NumericVector) - Sample data.
 *   iterMax (int) - Maximum number of iterations for the EM algorithm.
 *   EPSILON (double) - Convergence threshold.
 *   verbose (bool, optional) - Flag to control the output of algorithm progress.
 * Output:
 *   List - List containing estimated parameters and statistics about the computational performance. 
 *-------------------------------------------------------*/

// [[Rcpp::export]]
List ZINB_MLE(NumericVector X, int iterMax = 15000, double EPSILON = 1e-5,  bool verbose = false) {
    // Start time for performance measurement
    Rcpp::Function SysTime("Sys.time");
    NumericVector start_time = SysTime();

    // Initialization of the parameters
    int nbT = X.size(); 
    double eZero = sum(X == 0) / (1.0 * nbT); 
    NumericVector nonZeroX = X[X > 0];
    double esAvg = mean(nonZeroX);
    double esVar = var(nonZeroX);

    // Initial estimates
    double eSize = esAvg * esAvg / (esVar / (1 - eZero) - esAvg); 
    double eProb = (1 - eZero) * esAvg / esVar; 

    // Iteration starts
    int iter;
    for (iter = 0; iter < iterMax; iter++) {
        NumericVector prevParam = NumericVector::create(_["esAvg"] = esAvg, _["eSize"] = eSize, _["eZero"] = eZero);

        // E-Step: Expectation calculation
        NumericVector E1Vec(nbT);
        NumericVector E2Vec(nbT);
        NumericVector Dvec(nbT);

        for (int i = 0; i < nbT; i++) {
            E1Vec[i] = (X[i] == 0) ? eZero : 0;
            E2Vec[i] = (1 - eZero) * R::dnbinom(X[i], eSize, eProb, false);
            Dvec[i] = eSize * (R::digamma(eSize + X[i]) - R::digamma(eSize));
        }

        E1Vec = E1Vec / (E1Vec + E2Vec); // Normalization
        E2Vec = E2Vec / (E1Vec + E2Vec); // Normalization
        double beta = 1 - 1 / (1 - eProb) - 1 / log(eProb);

        // M-Step: Maximization to update parameters
        eZero = sum(E1Vec) / sum(E1Vec + E2Vec); 
        eProb = sum(E2Vec * (beta * Dvec)) / sum(E2Vec * (X + beta * Dvec - Dvec)); 
        eSize = -sum(E2Vec * Dvec) / sum(E2Vec) / log(eProb); 

        // Parameter validation and adjustment
        if (std::isnan(esAvg) || std::isnan(eSize) || std::isnan(eZero))
        {
            warning("ZINB_MLE - NA is appeared!");
            esAvg = prevParam["esAvg"];
            eSize = prevParam["eSize"];
            eZero = prevParam["eZero"];
            break;
        }

        if (eProb >= 1 || eProb <= 0 || eSize <= 0)
        {
            warning("Parameters out of bound, adjusting.");
            eProb = std::min(std::max(eProb, 0.001), 0.999);
            eSize = std::max(eSize, 0.001);
        }

        // Update mean and variance estimates
        esAvg = (1 - eZero) * eSize * (1 - eProb) / eProb; 
        esVar = esAvg / eProb + eZero / (1 - eZero) * esAvg * esAvg; 

        // Convergence check
        NumericVector newParam = NumericVector::create(_["esAvg"] = esAvg, _["eSize"] = eSize, _["eZero"] = eZero); 
        double D = sum(pow(newParam - prevParam, 2));
        if (D < EPSILON) {
            break;
        }
    }

    // Prepare result list
    List result;
    List parameter;
    parameter["eSize"] = eSize;
    parameter["eProb"] = eProb;
    parameter["eZero"] = eZero;
    parameter["esAvg"] = esAvg;
    parameter["esVar"] = esVar;
    result["Parameter_Estimates"] = parameter;    

    // Record stop time and print results if verbose
    NumericVector stop_time = SysTime();
    Rcpp::Function difftime("difftime");
    NumericVector minutes = difftime(stop_time, start_time, "mins");
    if (verbose)
    {
        Rcpp::Rcout << "Computational iterations: " << iter + 1 << "\n";
        Rcpp::Rcout << "Computational time: " << minutes << " minutes\n";
    }

    List stats;
    stats["Iterations"] = iter + 1;
    stats["ComputationTime"] = std::to_string(minutes[0]) + " minutes";
    result["Statistics"] = stats;

    // Return with a list
    return result;    
}

 /*** R
 res_ZINB<-ZINB_MLE(X)
 */