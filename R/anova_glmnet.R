#' Assign ANOVA type-II p-values to every submodel of a sparse `glmnet` model
#'
#' Given a sparse `glmnet` model (not ridge regression), here we assing ANOVA type-II p-values to every submodel obtained by varying the `lambda` penalty factor (i.e. each column of the `$beta` component matrix of the `glmnet` object).
#' To achieve this, each set of selected loci is fit again to the original data without penalization.
#'
#' @param beta The matrix of coefficients (component `$beta`) of the `glmnet` object.
#' @param X The genotype matrix.
#' Same as was used in [glmnet_pca()].
#' @param y The trait vector.
#' Same as was used in [glmnet_pca()].
#' @param pcs The PC (eigenvector) matrix (optional).
#' Same as was used in [glmnet_pca()].
#' Unlike genotypes, PCs are not assigned p-values.
#'
#' @return A type-II ANOVA p-value matrix with the same dimensions as `beta`.
#' Zero coefficients (unselected variables) are assigned p-values of 1.
#' For each column, selected variables are assigned p-values using [anova2()].
#'
#' @examples
#' \dontrun{
#' pvals <- anova_glmnet( beta, X, y, pcs )
#' }
#' 
#' @seealso
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

    # all p-values are 1 until recalculated (most stay 1, for variables that were not selected)
    pvals <- matrix( 1, nrow = m, ncol = k )
    # navigate each lambda value (column of beta)
    for ( j in 1 : k ) {
        # we don't actually use the coefficients, just presence/absence
        indexes <- which( beta[ , j ] != 0 )
        # sometimes nothing is selected, make sure we don't do anything in that case
        if ( length( indexes ) == 0 )
            next
        # also, if number of variables exceed sample size, we can't assign p-values that way either
        # (not expected for decent data sizes, but this is observed in toy test data)
        if ( length( indexes ) >= n )
            next
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
        # store now, R will complain if lengths are wrong
        pvals[ indexes, j ]  <- pvals_j
    }

    # done, return p-value matrix
    return( pvals )
}
