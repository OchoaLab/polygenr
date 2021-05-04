#' Score loci by appearance in parse `glmnet` model
#'
#' The `glmnet` package does not provide a way to rank variants.
#' Here is implemented one possible ranking based on when a variable is first selected as the penalty factor `lambda` of `glmnet` is varied.
#' `glmnet` precomputes coefficients at a series of `lambda` knots, which are used for this score calculation.
#' This approach requires sparsity (does not work for ridge regression).
#'
#' @param beta The matrix of coefficients (component `$beta`) of the `glmnet` object.
#'
#' @return The score vector.
#' Each locus receives a score between 0 and `k` (number of columns of `beta`).
#' Loci that were never selected have zero scores.
#' Loci first selected in column `i` have scores of `k + 1 - i`, so loci selected earlier have higher scores.
#'
#' @examples
#' \dontrun{
#' scores <- scores_glmnet( beta )
#' }
#'
#' @seealso
#' [anova_glmnet()] for a different way of scoring/raking variants.
#' 
#' @export
scores_glmnet <- function(beta) {
    if ( missing( beta ) )
        stop( '`beta` is required!' )

    # data dimensions
    k <- ncol( beta )
    m <- nrow( beta )
    
    # this is the desired data
    # initialized to zeroes (scores stay that way for the majority of loci in an ordinary glmnet fit)
    score <- vector("integer", m)
    # score is given by column of first appearance
    # analyze each row
    for ( i in 1 : m ) {
        # get non-zero columns
        index <- which( beta[ i, ] != 0 )
        
        # if there were any (most cases are zeros on all columns), then get first column
        # flip so score is highest for earlier appearances:
        # - max score is k if index[1] == 1
        # - min score for something picked last (index[1] == k) is 1, but it's 0 for never-picked cases
        if ( length( index ) > 0 )
            score[ i ] <- k + 1 - index[ 1 ]
    }
    
    return( score )
} 
