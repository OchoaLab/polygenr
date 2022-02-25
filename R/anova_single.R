#' Assign ANOVA type-II -log10 p-values to a sparse genetic model
#'
#' Given a sparse model defined by explicit locus indexes, here we assing ANOVA type-II -log10 p-values to these loci.
#' To achieve this, loci coefficients are fit to the data as multiple linear regression, without penalization.
#'
#' @param X The genotype matrix.
#' @param y The trait vector.
#' @param indexes Indexes of the loci to fit.
#' @param pcs The PC (eigenvector) matrix (optional).
#' Unlike genotypes, PCs are not given p-values.
#' @param ret_sparse Logical that controls return value (see that).
#'
#' @return If `ret_sparse = FALSE` (default), returns a complete vector of scores (-log10 p-values) for every locus in `X`, with zeroes for all loci with zero coefficients.
#' For loci with non-zero coefficients, p-values are calculated using [anova2()], see that for more details.
#' If `ret_sparse = TRUE`, returns a list of indexes and scores corresponding only to the loci with non-zero coefficients.
#'
#' @examples
#' \dontrun{
#' scores <- anova_single( X, y, indexes, pcs )
#' }
#' 
#' @seealso
#' [glmnet_pca()], particularly option `cv = TRUE`, for obtaining cross-validation objects with PCs.
#'
#' [anova_glmnet_single()] for scores on a single `glmnet` model, instead of specifying indexes explicitly as here.
#' 
#' [anova_glmnet()] for scores for all `glmnet` models (all lambdas), though it is much slower and not generally recommended.
#' 
#' [anova2()] for additional details and data restrictions.
#' 
#' [scores_glmnet()] for a different way of scoring/raking variants.
#' 
#' @export
anova_single <- function( X, y, indexes, pcs = NULL, ret_sparse = FALSE ) {
    # check mandatory arguments
    if ( missing( X ) )
        stop( '`X` is required!' )
    if ( missing( y ) )
        stop( '`y` is required!' )
    if ( missing( indexes ) )
        stop( '`indexes` is required!' )

    # check dimensions
    if ( !is.matrix( X ) )
        stop( '`X` must be a matrix!' )
    n <- ncol( X )
    m <- nrow( X )
    if ( length( y ) != n )
        stop( 'Length of `y` (', length( y ), ') must equal number of columns of `X` (', n, ')!' )
    if ( !is.null( pcs ) ) {
        if ( !is.matrix( pcs ) )
            stop( '`pcs` must be a matrix!' )
        if ( nrow( pcs ) != n )
            stop( 'Number of rows of `pcs` (', nrow( pcs ), ') must equal number of columns of `X` (', n, ')!' )
    }
    
    # initialize scores to zero
    # not needed in sparse return version
    if ( !ret_sparse )
        scores <- vector( 'numeric', m )
    
    if ( length( indexes ) == 0 ) {
        if ( ret_sparse ) {
            # return a list containint two empty vectors
            # this precise arrangement for `scores` was necessary for unit tests to succeed
            return( list( indexes = indexes, scores = numeric(0) ) )
        } else {
            return( scores )
        }
    }

    # if there's one or more values, can validate indexes further
    if ( min( indexes ) < 1 )
        stop( '`indexes` must be >= 1!' )
    if ( max( indexes ) > m )
        stop( '`indexes` must be <= `m` (', m, ')' )
    
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
