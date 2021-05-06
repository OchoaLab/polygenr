#' Assign ANOVA type-II -log10 p-values to one submodel of a sparse `glmnet` model
#'
#' Given a sparse `glmnet` model (not ridge regression), here we assing ANOVA type-II -log10 p-values to a single submodel obtained by specifying providing the cross-validation `glmnet` object and specifying index as "min" or "1se" (see [glmnet::cv.glmnet()] return object), or alternatively only the `beta` matrix (sparse coefficients) and the numeric index of the `lambda` penalty factor.
#' To achieve this, each set of selected loci is fit again to the original data without penalization.
#'
#' @param X The genotype matrix.
#' Same as was used in [glmnet_pca()].
#' @param y The trait vector.
#' Same as was used in [glmnet_pca()].
#' @param pcs The PC (eigenvector) matrix (optional).
#' Same as was used in [glmnet_pca()].
#' Unlike genotypes, PCs are not given p-values.
#' @param obj_cv Optional, the cross-validation object produced by [glmnet_pca()] with `cv = TRUE` (must be class "cv.glmfit").
#' Either this or `beta` must be provided.
#' @param beta Optional, the sparse matrix of coefficients (component `$beta`) of the `glmnet` object (output of [glmnet_pca()] with `cv = FALSE`, or component `$glmnet.fit$beta` if `cv = TRUE`).
#' Either this or `obj_cv` must be provided.
#' @param index The index for the desired `lambda` penalty factor, which implicitly chooses the submodel to analyze.
#' This can be a numeric index corresponding to a column of `beta`.
#' Alternatively, if `obj_cv` is provided, `index` can be "min" (default, selects `obj_cv$index[1]` as the actual index, which corresponds to the `lambda` that minimized the mean cross-validation error) or "1se" (`obj_cv$index[2]`, which is the largest `lambda` such that error is within 1 standard error of the minimum).
#' See [glmnet::cv.glmnet()] for more information on "min" and "1se" definitions.
#' @param ret_sparse Logical that controls return value (see that).
#'
#' @return If `ret_sparse = FALSE` (default), returns a complete vector of scores (-log10 p-values) for every locus in `X`, with zeroes for all loci with zero coefficients.
#' For loci with non-zero coefficients, p-values are calculated using [anova2()], see that for more details.
#' If `ret_sparse = TRUE`, returns a list of indexes and scores corresponding only to the loci with non-zero coefficients.
#'
#' @examples
#' \dontrun{
#' # version with cross-validation object `obj_cv` (recommended)
#' # defaults to selecting model with lowest cross-validation error (`index = "min"`)
#' scores <- anova_glmnet_single( X, y, pcs, obj_cv = obj_cv )
#'
#' # version with beta matrix and desired index
#' scores <- anova_glmnet_single( X, y, pcs, beta = beta, index = 50 )
#' }
#' 
#' @seealso
#' [glmnet_pca()], particularly option `cv = TRUE`, for obtaining cross-validation objects with PCs.
#'
#' [anova_glmnet()] for a version that calculates scores for all models (all lambdas), though it is much slower and not generally recommended.
#' 
#' [anova2()] for additional details and data restrictions.
#' 
#' [scores_glmnet()] for a different way of scoring/raking variants.
#' 
#' @export
anova_glmnet_single <- function( X, y, pcs = NULL, obj_cv = NULL, beta = NULL, index = "min", ret_sparse = FALSE ) {
    # check mandatory arguments
    if ( missing( X ) )
        stop( '`X` is required!' )
    if ( missing( y ) )
        stop( '`y` is required!' )

    # must have at least one of `obj_cv` or `beta`
    if ( is.null( obj_cv ) && is.null( beta ) )
        stop( 'Either `obj_cv` or `beta` must be non-NULL!' )
    # but don't pass both!
    if ( !is.null( obj_cv ) && !is.null( beta ) )
        stop( 'Only one of `obj_cv` or `beta` must be non-NULL!' )
    
    # assume only one was passed, handle that
    if ( !is.null( obj_cv ) ) {
        if ( !( "cv.glmnet" %in% class( obj_cv ) ) )
            stop( '`obj_cv` must be of class "cv.glmnet" (produced by `glmnet_pca` with `cv = TRUE`)!' )
        # copy down this subobject
        obj_glmnet <- obj_cv$glmnet.fit
        if ( is.null( obj_glmnet ) )
            stop( '`obj_cv$glmnet.fit` is `NULL` (unexpected)!' )
        # get beta from here
        beta <- obj_glmnet$beta
        if ( is.null( beta ) )
            stop( '`obj_cv$glmnet.fit$beta` is `NULL` (unexpected)!' )

        # extract index too
        if ( index %in% c('min', '1se') ) {
            # make sure the object has indexes
            if ( is.null( obj_cv$index ) )
                stop( '`obj_cv$index` is `NULL` (unexpected)!' )
            if ( index == 'min' ) {
                index <- obj_cv$index[1]
            } else if ( index == '1se' ) {
                index <- obj_cv$index[2]
            }
            # this index ought to be valid (it'd be surprising to get this far with a valid `obj_cv` except for this
            # will test further below together with non-`obj_cv` cases
        }
    } else {
        # otherwise `beta` was given and we already know it's not NULL
        # however let's check index
        if ( index %in% c('min', '1se') )
            stop( '`index` cannot equal "min" or "1se" unless `obj_cv` was provided!' )
    }

    # test index on all branches
    # now index must be a valid index
    if ( !is.integer( index ) )
        stop( '`index` was not integer!  Observed: ', index )
    if ( index < 1 )
        stop( '`index` must be positive!  Observed: ', index )
    
    # data dimensions
    if( !( 'dgCMatrix' %in% class( beta ) ) && !is.matrix( beta ) )
        stop( '`beta` must be a regular matrix or sparse matrix of class `dgCMatrix`!' )
    m <- nrow( beta )
    k <- ncol( beta )
    # check all matrices provided
    if ( !is.matrix( X ) )
        stop( '`X` must be a matrix!' )
    if ( nrow( X ) != m )
        stop( 'Number of rows of `X` (', nrow( X ), ') must equal number of rows of `beta` (', m, ')!' )
    n <- ncol( X )
    if ( length( y ) != n )
        stop( 'Length of `y` (', length( y ), ') must equal number of columns of `X` (', n, ')!' )
    if ( !is.null( pcs ) ) {
        if ( !is.matrix( pcs ) )
            stop( '`pcs` must be a matrix!' )
        if ( nrow( pcs ) != n )
            stop( 'Number of rows of `pcs` (', nrow( pcs ), ') must equal number of columns of `X` (', n, ')!' )
    }
    # make sure index is not out of range
    if ( index > k )
        stop( '`index` (', index, ') cannot exceed number of columns of `beta` (', k, ')!' )
    
    # initialize scores to zero
    # not needed in sparse return version
    if ( !ret_sparse )
        scores <- vector( 'numeric', m )
    
    # we don't actually use the coefficients, just presence/absence
    indexes <- which( beta[ , index ] != 0 )
    
    # sometimes nothing is selected, make sure we don't do anything in that case
    # this is fine for sparse `scores` matrix too (nothing was there in `beta`)
    if ( length( indexes ) == 0 ) {
        if ( ret_sparse ) {
            # return a list containint two empty vectors
            # this precise arrangement for `scores` was necessary for unit tests to succeed
            return( list( indexes = indexes, scores = numeric(0) ) )
        } else {
            return( scores )
        }
    }
    
    # also, if number of variables exceed sample size, we can't assign p-values that way either
    # (not expected for decent data sizes, but this is observed in toy test data)
    if ( length( indexes ) >= n ) {
        if ( ret_sparse ) {
            # return all zero scores
            return( list( indexes = indexes, scores = rep.int( 0, length( indexes ) ) ) )
        } else {
            return( scores )
        }
    }
    
    # get subset of genotypes
    Xs <- X[ indexes, , drop = FALSE ]
    
    # use this magic function to get the anova type-II p-values
    # passing null `pcs` is ok
    data <- anova2( Xs, y, pcs = pcs )
    
    # data is a data.frame, extract subset of interest, namely p-values
    # remove first row (intercept term, always equals NA)
    pvals_j <- data$p[ -1 ]
    
    # if there were pcs, they get the next row, remove it too
    if ( !is.null( pcs ) )
        pvals_j <- pvals_j[ -1 ]

    # at this point, indexes and scores_j should have the same length
    if ( length( pvals_j ) != length( indexes ) )
        stop( 'Length of `pvals_j` (', length( pvals_j ), ') and `indexes` (', length( indexes ), ') differs (unexpected)!' )
    
    # turn into -log10(p) scores
    scores_j <- -log10( pvals_j )

    # done, return scores vector
    if ( ret_sparse ) {
        return( list( indexes = indexes, scores = scores_j ) )
    } else {
        # store now in big vector
        scores[ indexes ] <- scores_j
        return( scores )
    }
}
