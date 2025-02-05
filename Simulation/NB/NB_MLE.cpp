// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
/*-------------------------------------------------------*
 * Function: Infor_Matrix
 * Purpose: Calculate the Fisher information matrix and convergence rate matrix for the negative binomial distribution.
 * Input:
 *   X (NumericVector) - Sample data.
 *   Param (List) - Estimated parameters of the negative binomial distribution.
 * Output:
 *   List containing expected, missing, observed information matrices, and the Jacobian (convergence rate matrix).
 *-------------------------------------------------------*/
// [[Rcpp::export]]
List Infor_Matrix(NumericVector X, List Param)
{
    // Retrieve the number of samples
    int nbT = X.size();

    // Extract estimated parameters
    double eProb = Param["eProb"];
    double eSize = Param["eSize"];

    // Initialize the information matrix and convergence rate matrix
    arma::mat expI(2, 2);
    arma::mat missI(2, 2);
    arma::mat obeI(2, 2);
    arma::mat J(2, 2);

    // Calculate terms for the Fisher information matrix
    double alpha = -1.0 / log(eProb);
    double term1 = 0.0;
    double term2 = 0.0;

    for (int i = 0; i < nbT; i++)
    {
        term1 += R::digamma(eSize + X[i]) - R::digamma(eSize);
        term2 += R::trigamma(eSize) - R::trigamma(eSize + X[i]);
    }

    // Fill in the expected information matrix
    expI(0, 0) = sum(X) / pow((1 - eProb), 2) - nbT * eSize / pow(eProb, 2) + ((pow(alpha, 2) + 2 * alpha) * pow((1 - eProb), 2) - eProb) * eSize * term1 / pow(eProb * (1 - eProb), 2);
    expI(0, 1) = nbT / eProb - 2 * alpha * term1 / eProb;
    expI(1, 0) = expI(0, 1);
    expI(1, 1) = term1 / eSize;

    // Fill in the missing information matrix
    missI(0, 0) = (pow((alpha * (1 - eProb)), 2) - eProb) * eSize * term1 / pow(eProb * (1 - eProb), 2);
    missI(0, 1) = 0;
    missI(1, 0) = 0;
    missI(1, 1) = term1 / eSize - term2;

    // Calculate the observed information matrix
    obeI = expI - missI;

    // Calculate the convergence rate matrix
    arma::mat expIinv(2, 2);
    bool invertible = arma::inv(expIinv, expI);
    if (!invertible)
    {
        warning("NB_MLE - Expected Information Matrix is invertible!");
    }
    J = expIinv * missI;

    List matrix;
    matrix["expected"] = expI;
    matrix["missing"] = missI;
    matrix["observed"] = obeI;
    matrix["jacobian"] = J;

    return matrix;
}

/*-------------------------------------------------------*
 * Function: NB_MLE
 * Purpose: Perform maximum likelihood estimation for the negative binomial distribution using EM algorithm.
 * Input:
 *   X (NumericVector) - Sample data.
 *   iterMax (int) - Maximum number of iterations for the EM algorithm.
 *   EPSILON (double) - Convergence threshold.
 *   verbose (bool) - Flag to control the output of the algorithm progress.
 *   se (bool) - Flag to compute standard errors.
 *   rate (bool) - Flag to compute convergence rate.
 * Output:
 *   List containing estimated parameters, standard errors, convergence rates, and other relevant information.
 *-------------------------------------------------------*/
// [[Rcpp::export]]
List NB_MLE(NumericVector X, int iterMax = 15000, double EPSILON = 1e-7, bool verbose = false, bool se = false, bool rate = false)
{
    // Start time for performance measurement
    Rcpp::Function SysTime("Sys.time");
    NumericVector start_time = SysTime();

    // Initialization of the parameters
    int nbT = X.size();
    double esAvg = Rcpp::mean(X);
    double esVar = Rcpp::var(X);

    // Validate input data
    if (esAvg >= esVar || esAvg <= 0)
    {
        stop("Invalid data: Variance should be greater than mean and mean must be positive.");
    }

    // Initial estimates
    double eSize = esAvg * esAvg / (esVar - esAvg);
    double eProb = esAvg / esVar;

    // Iteration starts
    int iter;
    for (iter = 0; iter < iterMax; iter++)
    {
        NumericVector prevParam = NumericVector::create(_["eProb"] = eProb, _["eSize"] = eSize);

        // E-Step: Calculate the digamma function for each data point
        NumericVector Dvec(nbT);
        for (int i = 0; i < nbT; i++)
        {
            Dvec[i] = eSize * (R::digamma(eSize + X[i]) - R::digamma(eSize));
        }
        double beta = 1 - 1.0 / (1 - eProb) - 1.0 / log(eProb);

        // M-Step: Update probability and size parameters
        eProb = sum(beta * Dvec) / sum(X + beta * Dvec - Dvec);
        eSize = -mean(Dvec) / log(eProb);

        // Parameter validation and adjustment
        if (std::isnan(eProb) || std::isnan(eSize))
        {
            warning("NB_MLE - NA is appeared!");
            eProb = prevParam["eProb"];
            eSize = prevParam["eSize"];
            break;
        }

        if (eProb >= 1 || eProb <= 0 || eSize <= 0)
        {
            warning("Parameters out of bound, adjusting.");
            eProb = std::min(std::max(eProb, 0.001), 0.999);
            eSize = std::max(eSize, 0.001);
        }

        // Update mean and variance estimates
        esAvg = eSize * (1 - eProb) / eProb;
        esVar = esAvg / eProb;

        // Convergence check
        NumericVector newParam = NumericVector::create(_["eProb"] = eProb, _["eSize"] = eSize);
        double D = sum(pow(newParam - prevParam, 2));
        if (D < EPSILON)
        {
            break;
        }
    }

    // Prepare result list
    List result;
    List parameter;
    parameter["eSize"] = eSize;
    parameter["eProb"] = eProb;
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

    // Compute standard errors if needed
    if (se)
    {
        List matrix = Infor_Matrix(X, parameter);
        arma::mat obeI = matrix["observed"];
        arma::mat obeIinv(2, 2);
        bool invertible = arma::inv(obeIinv, obeI);
        if (!invertible)
        {
            warning("NB_MLE - Observed Information Matrix is invertible!");
        }
        else
        {
            double SE_eProb = std::sqrt(obeIinv(0, 0));
            double SE_eSize = std::sqrt(obeIinv(1, 1));
            List SE;
            SE["SE_eProb"] = SE_eProb;
            SE["SE_eSize"] = SE_eSize;
            result["Standard_Error"] = SE;
        }
    }

    // Compute convergence rates if needed
    if (rate)
    {
        List matrix = Infor_Matrix(X, parameter);
        arma::mat J = matrix["jacobian"];
        double Rate_eProb = J(0, 0);
        double Rate_eSize = J(1, 1);
        List Rate;
        Rate["Rate_eProb"] = Rate_eProb;
        Rate["Rate_eSize"] = Rate_eSize;
        result["Convergence_Rate"] = Rate;
    }

    List stats;
    stats["Iterations"] = iter + 1;
    stats["ComputationTime"] = std::to_string(minutes[0]) + " minutes";
    result["Statistics"] = stats;

    // Return with a list
    return result;
}