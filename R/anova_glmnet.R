#' Assign ANOVA type-II -log10 p-values to every submodel of a sparse `glmnet` model
#'
#' Given a sparse `glmnet` model (not ridge regression), here we assing ANOVA type-II -log10 p-values to every submodel obtained by varying the `lambda` penalty factor (i.e. each column of the `$beta` component matrix of the `glmnet` object).
#' To achieve this, each set of selected loci is fit again to the original data without penalization.
#'
#' @param beta The matrix of coefficients (component `$beta`) of the `glmnet` object.
#' @param X The genotype matrix.
#' Same as was used in [glmnet_pca()].
#' @param y The trait vector.
#' Same as was used in [glmnet_pca()].
#' @param pcs The PC (eigenvector) matrix (optional).
#' Same as was used in [glmnet_pca()].
#' Unlike genotypes, PCs are not given p-values.
#'
#' @return A sparse matrix (class `dgCMatrix`) with the same dimensions as `beta`, containing type-II ANOVA -log10 p-values.
#' Zero coefficients (unselected variables) are assigned values of zero as well (to retain sparsity, imply p-values of 1).
#' For selected variables in each column, p-values are calculated using [anova2()], see that for more details.
#'
#' @examples
#' \dontrun{
#' scores <- anova_glmnet( beta, X, y, pcs )
#' }
#' 
#' @seealso
#' [anova_glmnet_single()] for calculations on a single model (by default, approximately the best) instead of all models (all lambdas), which is much faster and generally recommended.
#' 
#' [anova2()] for additional details and data restrictions.
#' 
#' [scores_glmnet()] for a different way of scoring/raking variants.
#' 
#' @export
anova_glmnet <- function( beta, X, y, pcs = NULL ) {
    # check mandatory arguments
    if ( missing( beta ) )
        stop( '`beta` is required!' )
    if ( missing( X ) )
        stop( '`X` is required!' )
    if ( missing( y ) )
        stop( '`y` is required!' )
    
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

    ## # all p-values are 1 until recalculated (most stay 1, for variables that were not selected)
    ## pvals <- matrix( 1, nrow = m, ncol = k )
    # copy `beta` matrix, this way `scores` is also a compressed, sparse matrix with missing data in the same places
    scores <- beta
    # navigate each lambda value (column of beta)
    for ( j in 1 : k ) {
        # use this function to perform the bulk of the calculations
        scores_sparse <- anova_glmnet_single( X, y, pcs = pcs, beta = beta, index = j, ret_sparse = TRUE )
        
        # sometimes nothing is selected, make sure we don't do anything in that case
        # this is fine for sparse `scores` matrix too (nothing was there in `beta`)
        if ( length( scores_sparse$indexes ) == 0 )
            next
        
        # in all other cases, copy over this way
        scores[ scores_sparse$indexes, j ] <- scores_sparse$scores
    }
    
    # done, return p-value matrix
    return( scores )
}
