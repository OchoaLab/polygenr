#' Get ANOVA type-II p-values from genetic model
#'
#' This function assings p-values to a polyygenic linear model.
#' Essentially a wrapper around [stats::drop1()] to assign p-values to a small number of selected loci while fitting PCs without caring about their significance.
#'
#' @param X The genotype matrix.
#' Each row of `X` receives a p-value.
#' Unlike other functions, where `X` is expected to be genome-wide, here we don't use penalized models so `X` must have fewer predictors (loci) than samples (individuals).
#' These predictiors are expected to have been selected in a previous step by `glmnet` or another such approach.
#' @param y The trait vector.
#' It must have length equal to the number of individuals.
#' @param pcs The PC (eigenvector) matrix (optional).
#' It must have individuals along the rows and dimensions along the columns.
#' Unlike `X`, all of `pcs` receives a single p-value (even if it contains multiple columns; for computational efficiency, it is assumed that selection or ranking of individual PCs is not of interest).
#'
#' @return A data frame containing these columns: 'Df', 'SS', 'RSS', 'AIC', 'F', 'p'.
#' The first row is for the full model, followed by statistics for `pcs` if present, followed by per-locus statistics for `X`.
#' If `X` had rownames, these are the names of the variables, otherwise 'x1', 'x2', and so on are used.
#'
#' @examples
#' \dontrun{
#' data <- anova2(X, y, pcs)
#' }
#'
#' @export
anova2 <- function( X, y, pcs = NULL ) {
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
    if ( !is.null( pcs ) ) {
        if ( !is.matrix( pcs ) )
            stop( '`pcs` must be a matrix!' )
        if ( nrow( pcs ) != n )
            stop( 'Number of rows of `pcs` (', nrow( pcs ), ') must equal number of columns of `X` (', n, ')!' )
    }

    # let's hack the setup so it works
    # for this we need to put everything into a data frame and need names
    # however, keep existing names if present (are there possible errors if the names look like math?)
    if ( is.null( rownames( X ) ) )
        rownames( X ) <- paste0( 'x', 1 : m )
    # response automatically gets called "y" here!  (no need to add names to it later)
    data <- cbind( y, as.data.frame( t( X ) ) )
    # now need to hack formula, saw example in ?formula
    # use X's names again, they must match!
    # here y also matches manually
    # first we have string version only, and exclude response
    formula <- paste( rownames( X ), collapse= "+" )
    
    if ( !is.null( pcs ) ) {
        # add PCs to data, but since we don't want these separated, have to add them this way
        data$pcs <- pcs
        # and add to formula as well, before all the actual terms of interest
        formula <- paste0( 'pcs + ', formula )
    }
    
    # add response now
    formula <- paste0( "y ~ ", formula )
    # convert to proper formula
    formula <- stats::as.formula( formula )
    
    # this works as desired!!!
    # fit linear model
    model <- stats::lm( formula, data = data )
    # because it's glmnet, there can be more predictors than data
    # regardless, there's a risk for perfect fits to the data
    # copy a criterion used inside `drop1.lm` (`stats:::check_exact`) to decide if this is worth doing
    mss <- sum( model$fitted.values^2 )
    rss <- sum( model$residuals^2 )
    # if we see an essentially perfect fit, return a dummy data frame with just p-values of 1
    # the length of this data depends on whether `pcs` is null or not
    if (rss < 1e-10 * mss) 
        return( data.frame( p = rep.int( 1, m + 1 + !is.null( pcs ) ) ) )
    
    # now calculate type-II ANOVA p-values
    obj <- stats::drop1( model, test = 'F' )
    # the data frame version of these results is most useful, return that!
    obj <- as.data.frame( obj )
    # simplify names
    names( obj )[ names( obj ) == 'Sum of Sq' ] <- 'SS'
    names( obj )[ names( obj ) == 'F value' ] <- 'F'
    names( obj )[ names( obj ) == 'Pr(>F)' ] <- 'p'
    # TODO: rename columns, some are ridiculous
    return( obj )
}
