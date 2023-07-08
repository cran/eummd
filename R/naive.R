#' Naive computation for Maximum Mean Discrepancy
#'
#' Computes maximum mean discrepancy statistics with Laplacian 
#' or Gaussian kernel. 
#' Suitable for multivariate data. Naive approach, quadratic in number
#' of observations.
#' 
#' @param X Matrix (or vector) of observations in first sample.
#' 
#' @param Y Matrix (or vector) of observations in second sample.
#' 
#' @param beta kernel parameter. Must be positive; if not, computes
#'             median heuristic in quadratic time. Default value
#'             is \code{-0.1}, which will force median heuristic to be used.
#' 
#' @param pval Boolean for whether to compute p-value or not. 
#' 
#' @param kernel String, either \code{"Laplacian"} or \code{"Gaussian"}. 
#'               Default is \code{"Laplacian"}.
#'
#' @param numperm Number of permutations. Default is \code{200}.
#'
#' @param seednum Seed number for generating permutations. Default is \code{0}, 
#'                which means seed is set randomly. For values larger than 
#'                \code{0}, results will be reproducible.
#' 
#' @details First checks number of columns (dimension) are equal. 
#'          Suppose matrix \code{X} has \code{n} rows and \code{d} columns, 
#'          and matrix \code{Y} has \code{m} rows; checks that \code{Y} 
#'          has \code{d} columns (if not, then throws error). 
#'          Then flattens matrices to vectors (or, if \code{d=1}, they are
#'          already vectors.
#'          Then calls C++ method. If the first sample has \code{n} 
#'          \code{d}-dimensional samples and the second sample has 
#'          \code{m} \code{d}-dimensional samples, then the algorithm
#'          computes the statistic in \code{O( (n+m)^2 )} time.
#'          
#'          Median difference is as follows:
#'          
#' \code{ m = median { ||x_i - x_j||_1; i=1, 2, ..., n+m, and j=1, 2,..., i } },
#'          
#'          where \code{||x_i - x_j||_1} is the 1-norm, and so if the data 
#'          are \code{d}-dimensional then
#'          
#'          \code{ ||x_i - x_j||_1 = \sum_{k=1}^{d} |x_{i,k} - x_{j,k}| },
#'          
#'          and finally median heuristic is \code{beta = 1/m}.
#'          This can be computed in \code{O( (n+m)^2 )} time.
#'          NOTE: fix.
#'          
#'          The Laplacian kernel \code{k} is defined as 
#'         
#'          \code{ k(x,y) = \exp( -\beta ||x - y||_1  ) }.
#'
#'          Random seed is set for \code{std::mt19937} and \code{std::shuffle} 
#'          in C++.
#'
#' @return A list with the following elements:
#'         \describe{
#'             \item{\code{pval}}{The p-value of the test, if it is  
#'                                computed (\code{pval=TRUE}). }
#'             \item{\code{stat}}{The statistic of the test, which
#'                                is always computed. }
#'             \item{\code{beta}}{The kernel parameter used in the test.
#'                                If \code{beta} was not initialised or
#'                                negative, this will be the median heuristic
#'                                value.}
#'          }
#'
#' @examples
#'
#' X <- matrix(c(1:12), ncol=2, byrow=TRUE)
#' Y <- matrix(c(13:20), ncol=2, byrow=TRUE)
#' mmdList <- mmd(X=X, Y=Y, beta=0.1, pval=FALSE)
#'
#'
#' @export
mmd <- function(X, Y, beta=-0.1, pval=TRUE, kernel=c("Laplacian", "Gaussian"), 
                numperm=200, seednum=0){

    # check vectors/matrices are numeric
    if ( !(is.numeric(X)) ){
        stop("X needs to be numeric.")
    }
    if ( !(is.numeric(Y)) ){
        stop("Y needs to be numeric.")
    }

    # check kernel is correct
    kernel <- kernel[1]
    if ( (kernel != "Laplacian") && (kernel != "Gaussian") ){
        stop("kernel needs to be either 'Laplacian' or 'Gaussian'.")
    }

    # initialise, will update later
    nX <- 0
    dX <- 0
    nY <- 0
    dY <- 0
    Xvec <- c()
    Yvec <- c()

    # add checks for matrix/vectors
    if (is.matrix(X)){
        # get dimensions of matrices
        nX <- dim(X)[1]
        dX <- dim(X)[2]

        if (!(is.matrix(Y))){
            stop("If X is a matrix, Y must also be a matrix.")
        }

        nY <- dim(Y)[1]
        dY <- dim(Y)[2]

        # check dimensions compatible
        if (dX != dY){
            stop("Dimension (number of columns) of matrices need to be equal.")
        }

        # flatten to vectors; NOTE: we use transpose
        # > X
        #      [,1] [,2] [,3]
        # [1,]    1    3    5
        # [2,]    2    4    6
        #
        # will be flattened to 
        # [1] 1 3 5 2 4 6
        # which is what we want

        Xvec <- as.vector(t(X))
        Yvec <- as.vector(t(Y))
    } else if (is.vector(X)){
        if (!(is.vector(Y))){
            stop("If X is a vector, Y must also be a vector.")
        }
        nX <- length(X)
        dX <- 1
        nY <- length(Y)
        dY <- 1
        Xvec <- X
        Yvec <- Y
    } else {
        stop("X must be a vector or a matrix.")
    }


    # if beta not positive, compute median heuristic (also from C++)
    # finally, compute MMD 
    mmdList <- list()
    if (kernel=="Gaussian"){
        if (pval){
            mmdList <- mmd_gau_pval_Rcpp(Xvec, Yvec, nX, dX, nY, dY, 
                                         numperm, seednum, beta)
        } else {
            mmdList <- mmd_gau_Rcpp(Xvec, Yvec, nX, dX, nY, dY, beta)
        }
    } else {
        if (pval){
            mmdList <- mmd_lap_pval_Rcpp(Xvec, Yvec, nX, dX, nY, dY, 
                                         numperm, seednum, beta)
        } else {
            mmdList <- mmd_lap_Rcpp(Xvec, Yvec, nX, dX, nY, dY, beta)
        }
    }

    # check pval; if no pval, functions return -1, then changes to NA
    if (mmdList$pval < 0){
        mmdList$pval <- NA
    }
    return(mmdList)
}
