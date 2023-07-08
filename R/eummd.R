#' euMMD: Effiicient Univariate Maximum Mean Discrepancy
#'
#' Computes maximum mean discrepancy statistics with Laplacian kernel. 
#' Suitable only for univariate data. Computing the statistic alone
#' for \code{n} observations is \code{O(n log n)}, and computing the 
#' p-value for \code{L} permutations is \code{O(n log n + Ln)}.
#' 
#' @param x Univariate vector of observations in first sample.
#' 
#' @param y Univariate vector of observations in second sample.
#' 
#' @param beta kernel parameter. Must be positive; if not, computes
#'             median heuristic in quadratic time. Default value
#'             is \code{-0.1}, which will force median heuristic to be used.
#' 
#' @param pval Boolean for whether to compute p-value or not. 
#' 
#' @param numperm Number of permutations. Default is \code{200}.
#'
#' @param seednum Seed number for generating permutations. Default is \code{0}, 
#'                which means seed is set randomly. For values larger than 
#'                \code{0}, results will be reproducible.
#' 
#' @details If the total number of observations in both samples is \code{n}, 
#'          first sort combined sample in \code{O(n log n)} before remaining
#'          steps are linear in \code{n}.
#'          
#'          If \code{beta} is not a positive value, 
#'          median difference is computed as follows:
#'          
#' \code{ m = median { ||x_i - x_j||_1; i=1, 2, ..., n+m, and j=1, 2,..., i } },
#'          
#'          where \code{||x_i - x_j||_1} is the 1-norm, and so if the data 
#'          are univariate then
#'          
#'          \code{ ||x_i - x_j||_1 = |x_{i} - x_{j}| },
#'          
#'          and finally median heuristic is \code{beta = 1/m}.
#'          This can be computed in \code{O(n log n )} time
#'          using the algorithms of Johnson and Mizoguchi (1978) 
#'          and Croux and Rousseuw (1992); see \code{mediandiff}
#'          for references.
#'          
#'          The Laplacian kernel \code{k} is defined as 
#'         
#'          \code{ k(x,y) = \exp( -\beta ||x - y||_1  ) }.
#'
#'          Random seed is set for \code{std::mt19937} and \code{std::shuffle} 
#'          in C++.
#'
#'
#'
#' @return A list with the following elements:
#'         \describe{
#'             \item{\code{pval}}{The p-value of the test, if it is  
#'                                computed (\code{pval=TRUE}). Otherwise, 
#'                                it is set to \code{NA}.}
#'             \item{\code{stat}}{The statistic of the test, which
#'                                is always computed. }
#'             \item{\code{beta}}{The kernel parameter used in the test.
#'                                If \code{beta} was not initialised or
#'                                negative, this will be the median heuristic
#'                                value.}
#'          }
#' 
#' @seealso mediandiff
#'
#' @examples 
#'
#' x <- c(7.1, 1.2, 4.3, 0.4)
#' y <- c(5.5, 2.6, 8.7)
#' mmd_list <- eummd(x, y, seednum=1)
#'
#' @export
eummd <- function(x, y, beta=-0.1, pval=TRUE, numperm=200, seednum=0){

    # check vectors are numeric
    if ( !(is.numeric(x)) || !(is.vector(x)) ){
        stop("x needs to be a numeric vector.")
    }
    if ( !(is.numeric(y)) || !(is.vector(y))){
        stop("y needs to be a numeric vector.")
    }

    mmdList <- list()
    if (pval){
        mmdList <- eummd_pval_Rcpp(x, y, beta, numperm, seednum)
    } else {
        mmdList <- eummd_Rcpp(x, y, beta)
    }

    # check pval; if no pval, functions return -1, then changes to NA
    if (mmdList$pval < 0){
        mmdList$pval <- NA
    }
    return(mmdList)
}
