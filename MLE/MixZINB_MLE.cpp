// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

arma::vec vec2arma(Rcpp::NumericVector x)
{
    return arma::vec(x.begin(), x.size(), false); // this code originate from P159 Seamless Rcpp
}
template <typename T>
Rcpp::NumericVector arma2vec(const T &x)
{
    return Rcpp::NumericVector(x.begin(), x.end());
}
mat E1(NumericVector X, NumericVector eZero)
{
    int nbT = X.size(), nbS = eZero.size();
    mat res(nbS, nbT);

    for (unsigned i = 0; i < nbS; i++)
    {
        for (unsigned j = 0; j < nbT; j++)
        {
            double temp = (X(j) == 0 ? eZero[i] : 0);
            res(i, j) = temp;
        }
    }
    return res;
}
mat E2(NumericVector X, NumericVector eSize, NumericVector eProb, NumericVector eZero)
{
    int nbT = X.size(), nbS = eProb.size();
    mat res(nbS, nbT);

    for (unsigned i = 0; i < nbS; i++)
    {
        for (unsigned j = 0; j < nbT; j++)
        {
            double temp = (1 - eZero[i]) * R::dnbinom(X(j), eSize(i), eProb(i), false);
            res(i, j) = temp;
        }
    }
    return res;
}
mat digamma(mat X)
{
    mat res(X.n_rows, X.n_cols);
    for (unsigned int i = 0; i < X.n_rows; i++)
    {
        for (unsigned int j = 0; j < X.n_cols; j++)
        {
            double tmp = R::digamma(X(i, j));
            res(i, j) = tmp;
        }
    }
    return res;
}
/*-------------------------------------------------------*
 * Function: mixZINB_MLE
 * Purpose: Perform maximum likelihood estimation for the mixed zero inflated negative binomial distribution using an EM algorithm.
 * Input:
 *   X (NumericVector) - Sample data.
 *   nbS (int) - Number of subpopulations or components in the mixture.
 *   avg (NumericVector, optional) - Initial mean values for each component. If not provided, k-means clustering is used for initialization.
 *   var (NumericVector, optional) - Initial variance values for each component. Used alongside 'avg' if provided.
 *   zero (NumericVector, optional) - Initial zero inflation coefficient values for each component. Used alongside 'zero' if provided.
 *   iterMax (int) - Maximum number of iterations for the EM algorithm.
 *   EPSILON (double) - Convergence threshold for the EM algorithm.
 * Output:
 *   List - A list containing estimated parameters including weights, sizes, probabilities, zeros, means, and variances for each component.
 * Description:
 *   This function estimates the parameters of a mixed zero inflated negative binomial distribution, which is a mixture of several zero inflated negative binomial distributions.
 *   The EM algorithm is used to iteratively estimate the parameters based on the provided data.
 *   The function handles initialization and provides detailed output of the estimated parameters.
 *-------------------------------------------------------*/
// [[Rcpp::export]]
List mixZINB_MLE(NumericVector X, int nbS, NumericVector avg = NumericVector(), NumericVector var = NumericVector(), NumericVector zero = NumericVector(), int iterMax = 15000, double EPSILON = 1e-5, bool verbose = true)
{
    // Start time for performance measurement
    Rcpp::Function SysTime("Sys.time");
    NumericVector start_time = SysTime();

    int nbT = X.size(); // Number of samples
    // Initialization of parameters
    NumericVector eWeight(nbS, 1.0 / nbS);
    NumericVector esAvg(nbS);
    NumericVector esVar(nbS);
    NumericVector eZero(nbS, sum(X == 0) * 1.0 / nbT);

    // Initialize zero inflation coefficients if not provided
    if (zero.size() != 0)
    {
        eZero = clone(zero);
    }

    // Initialize means and variances if not provided
    if (avg.size() == 0 || var.size() == 0)
    {
        NumericVector nonZeroX = X[X > 0];
        Rcpp::Function kmeans("kmeans");
        List res_kmeans = kmeans(nonZeroX, nbS); // Perform k-means clustering for initialization

        NumericVector centers = res_kmeans[1];
        NumericVector withinss = res_kmeans[3];
        NumericVector size = res_kmeans[6];

        esAvg = clone(centers);
        esVar = withinss / (size - 1);
    }
    else
    {
        esAvg = clone(avg);
        esVar = clone(var);
    }

    // Check for invalid initial values
    if (is_true(any(esAvg * (1 - eZero) >= esVar)))
    {
        stop("Invalid data: Variance should be greater than mean for each cluster.");
    }

    // Calculate initial size and probability for each component
    NumericVector eSize = pow(esAvg, 2) / (esVar / (1 - eZero) - esAvg);
    NumericVector eProb = (1 - eZero) * esAvg / esVar;

    // Prepare matrices for the EM algorithm
    mat ones_nbT = ones(nbT, 1);
    mat ones_nbS = ones(nbS, 1);
    mat ones_nbSS = ones(nbS, nbS);

    // Iteration starts
    int iter;
    for (iter = 0; iter < iterMax; iter++)
    {
        mat prevParam = join_cols(vec2arma(esAvg), vec2arma(eSize), vec2arma(eWeight), vec2arma(eZero));

        // E-Step: Expectation calculation
        mat Xmat = ones_nbS * vec2arma(X).t();
        mat Nmat = vec2arma(eSize) * ones_nbT.t();

        mat E1mat = E1(X, eZero);
        mat E2mat = E2(X, eSize, eProb, eZero);
        mat Hmat = E1mat + E2mat; // Emission probabilities

        E1mat = E1mat / Hmat; // Normalization
        E2mat = E2mat / Hmat; // Normalization

        mat Pmat = Hmat % (vec2arma(eWeight) * ones_nbT.t()); // Weighted probabilities
        mat Tmat = Pmat / (ones_nbSS * Pmat);                 // Normalized probabilities

        mat Dmat = (digamma(Xmat + Nmat) - digamma(Nmat)) % Nmat; // Digamma differences for each sample and component

        NumericVector beta = 1 - 1 / (1 - eProb) - 1 / log(eProb); // Beta values for each component
        mat Bmat = vec2arma(beta) * ones_nbT.t();

        // M-Step: Maximization to update parameters
        mat eWeight_mat = (Tmat / (1.0 * nbT)) * ones_nbT;

        mat eZero_mat = ((Tmat % E1mat) * ones_nbT) / ((Tmat % (E1mat + E2mat)) * ones_nbT);

        mat denominator = (Tmat % E2mat % Dmat % Bmat) * ones_nbT;
        mat numerator = (Tmat % E2mat % (Xmat + Dmat % Bmat - Dmat)) * ones_nbT;
        mat eProb_mat = denominator / numerator;

        mat temp_mat = ((Tmat % E2mat % Dmat) * ones_nbT) / ((Tmat % E2mat) * ones_nbT);
        mat eSize_mat = -temp_mat / log(eProb_mat);

        mat esAvg_mat = (1 - eZero_mat) % eSize_mat % (1 - eProb_mat) / eProb_mat;
        mat esVar_mat = esAvg_mat / eProb_mat + eZero_mat / (1 - eZero_mat) % pow(esAvg_mat, 2);

        eWeight = arma2vec(eWeight_mat); // new value of weight
        eZero = arma2vec(eZero_mat);     // new value of zero
        eProb = arma2vec(eProb_mat);     // new value of prob
        eSize = arma2vec(eSize_mat);     // new value of size
        esAvg = arma2vec(esAvg_mat);     // new value of mean
        esVar = arma2vec(esVar_mat);     // new value of variance

        if (any(Rcpp::is_na(esAvg)) | any(Rcpp::is_na(eSize)) | any(Rcpp::is_na(eWeight)) | any(Rcpp::is_na(eZero)))
        {
            warning("mixZINB_MLE - NA is appeared!");
            esAvg = prevParam.col(0);
            eSize = prevParam.col(1);
            eWeight = prevParam.col(2);
            eZero = prevParam.col(3);
            break;
        }

        // Check for convergence
        mat newParam = join_cols(esAvg_mat, eSize_mat, eWeight_mat, eZero_mat);
        double D = arma::accu(arma::square(newParam - prevParam));
        if (D < EPSILON)
        {
            break;
        }
    }

    // return with a list
    List result;
    List parameter;
    parameter["eWeight"] = eWeight;
    parameter["eZero"] = eZero;
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

    List stats;
    stats["Iterations"] = iter + 1;
    stats["ComputationTime"] = std::to_string(minutes[0]) + " minutes";
    result["Statistics"] = stats;

    // Return with a list
    return result;
}

 /*** R
 res_MixZINB<-mixZINB_MLE(X,nbS)
 */