#' Fit a polygenic model with `glmnet` and principal components
#'
#' This function makes it easy to fit a LASSO or Elastic Net model while including principal components or other fixed effects that are not part of the variable selection process and are not penalized.
#' This resembles in spirit and functionality the `snpnet` function of the same package, except that one requires the genetic data to be written into a plink2 file whereas this one requires inputs that are ordinary R matrices.
#' The cross-validation `glmnet` function is also optionally adapted here.
#'
#' @param X The genotype matrix.
#' It must be oriented with loci along rows and individuals along columns, which agrees with other genetics packages by the main author, although this is transposed from what `glmnet` and other regression models expect (even `lm`).
#' @param y The trait vector.
#' It must have length equal to the number of individuals.
#' @param pcs The PC (eigenvector) matrix (optional).
#' It must have individuals along the rows and dimensions along the columns.
#' Unlike `X`, predictors in this `pcs` matrix, are not penalized, and their coefficients are not included in the output.
#' @param cv If `TRUE`, run [glmnet::cv.glmnet()] instead of [glmnet::glmnet()], returning that respective output (with small modifications).
#' Recommended to help pick the best `lambda` value for the penalized regression.
#' Default `FALSE`.
#' In both cases, PCs are modeled the same special way described above.
#' @param ... Additional parameters passed to [glmnet::glmnet()] or [glmnet::cv.glmnet()].
#' Cannot include `penalty.factor`, which is set internally so that loci are equally penalized and `pcs` are unpenalized.
#'
#' @return If `cv = FALSE`, a `glmnet` object (see [glmnet::glmnet()]), otherwise a `cv.glmnet` object (see [glmnet::cv.glmnet()]).
#' However, the respective `beta` matrices (`obj$beta` or `obj$glmnet.fit$beta`, respectively) are modified to exclude coefficients for the `pcs` elements (only include coefficients for loci in `X`).
#'
#' @examples
#' \dontrun{
#' # regular glmnet
#' obj <- glmnet_pca(X, y, pcs)
#' 
#' # cross validation version
#' obj <- glmnet_pca(X, y, pcs, cv = TRUE)
#' }
#'
#' @seealso
#' [glmnet::glmnet()] and [glmnet::cv.glmnet()] for the functions this function wraps around, to better understand additional options and return values.
#'
#' @export
glmnet_pca <- function(X, y, pcs = NULL, cv = FALSE, ...){
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
        X <- t( X )
        # decide whether to run cross-validation version or not
        if ( cv ) {
            obj <- glmnet::cv.glmnet( X, y, ... )
        } else {
            obj <- glmnet::glmnet( X, y, ... )
        }
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
        # decide whether to run cross-validation version or not
        if ( cv ) {
            obj <- glmnet::cv.glmnet( X, y, penalty.factor = penalty.factor, ... )
        } else {
            obj <- glmnet::glmnet( X, y, penalty.factor = penalty.factor, ... )
        }

        # clean up return object
        # (remove PCs from coefficients, plots look very weird otherwise)
        # have to do same thing on both modes, but the object name/place is difference
        obj_glmnet <- if ( cv ) obj$glmnet.fit else obj 
        obj_glmnet$beta <- obj_glmnet$beta[ - ( 1 : r ) , ]
        obj_glmnet$dim[ 1 ] <- obj_glmnet$dim[ 1 ] - r
        # write edited object back to return object, as needed
        if ( cv ) {
            obj$glmnet.fit <- obj_glmnet
        } else {
            obj <- obj_glmnet
        }
    }
    
    return (obj)
}
