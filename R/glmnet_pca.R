#' Fit a polygenic model with `glmnet` and principal components
#'
#' This function makes it easy to fit a LASSO or Elastic Net model while including principal components or other fixed effects that are not part of the variable selection process and are not penalized.
#' This resembles in spirit and functionality the `snpnet` function of the same package, except that one requires the genetic data to be written into a plink2 file whereas this one requires inputs that are ordinary R matrices.
#'
#' @param X The genotype matrix.
#' It must be oriented with loci along rows and individuals along columns, which agrees with other genetics packages by the main author, although this is transposed from what `glmnet` and other regression models expect (even `lm`).
#' @param y The trait vector.
#' It must have length equal to the number of individuals.
#' @param pcs The PC (eigenvector) matrix (optional).
#' It must have individuals along the rows and dimensions along the columns.
#' Unlike `X`, predictors in this `pcs` matrix, are not penalized, and their coefficients are not included in the output.
#' @param ... Additional parameters passed to `glmnet`.
#' Cannot include `penalty.factor`, which is set internally so that loci are equally penalized and `pcs` are unpenalized.
#'
#' @return A `glmnet` object.
#' Here it includes coefficients for every locus in `X` (can be sparse) but excludes coefficients for the `pcs` elements.
#'
#' @examples
#' \dontrun{
#' obj <- glmnet_pca(X, y, pcs)
#' }
#'
#' @export
glmnet_pca <- function(X, y, pcs = NULL, ...){
    # make sure all the required arguments are there
    if ( missing( X ) )
        stop( '`X` is required!' )
    if ( missing( y ) )
        stop( '`y` is required!' )

    # check dimensions
    if ( !is.matrix( X ) )
        stop( '`X` must be a matrix!' )
    m <- nrow( X )
    n <- ncol( X )
    if ( length( y ) != n )
        stop( 'Length of `y` (', length( y ), ') must equal number of columns of `X` (', n, ')!' )
    
    if ( is.null( pcs ) ) {
        # without PCs, this is just regular glmnet
        # transpose X now
        obj <- glmnet::glmnet( t(X), y, ... )
    } else {
        # here is all the fun

        # check dimensions for this too
        if ( !is.matrix( pcs ) )
            stop( '`pcs` must be a matrix!' )
        if ( nrow( pcs ) != n )
            stop( 'Number of rows of `pcs` (', nrow( pcs ), ') must equal number of columns of `X` (', n, ')!' )
        
        # number of PCs
        r <- ncol( pcs )
        # add as more fixed effects, except we won't penalize them
        # add first for simplicity
        # transpose X now
        X <- cbind( pcs, t( X ) )
        # also need to not penalize PCs
        penalty.factor <- rep.int( 1, r + m ) # normal penalties for everything initially
        penalty.factor[ 1 : r ] <- 0 # no penalty for PCs

        # run enhanced glmnet to handle PCs appropriately
        obj <- glmnet::glmnet( X, y, penalty.factor = penalty.factor, ... )
        
        # clean up return object
        # (remove PCs from coefficients, plots look very weird otherwise)
        obj$beta <- obj$beta[ - ( 1 : r ) , ]
        obj$dim[ 1 ] <- obj$dim[ 1 ] - r
    }
    
    return (obj)
}
