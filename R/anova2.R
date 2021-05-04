#' Get ANOVA type-II p-values from genetic model
#'
#' This function assings p-values to a polyygenic linear model.
#' Essentially a wrapper around [stats::drop1()] to assign p-values to a small number of selected loci while fitting PCs without caring about their significance.
#'
#' @param X The genotype matrix.
#' Each row of `X` receives a p-value.
#' Unlike other functions, where `X` is expected to be genome-wide, here we don't use penalized models so `X` must have fewer predictors (loci) than samples (individuals; if this happens a vector of 1s for all locus p-values is returned).
#' These predictiors are expected to have been selected in a previous step by `glmnet` or another such approach.
#' Row names of `X` may contain special characters (including math operations, colons, and spaces; names will be quoted in output) except backticks are not allowed and cause an error.
#' @param y The trait vector.
#' It must have length equal to the number of individuals.
#' @param pcs The PC (eigenvector) matrix (optional).
#' It must have individuals along the rows and dimensions along the columns.
#' Unlike `X`, all of `pcs` receives a single p-value (even if it contains multiple columns; for computational efficiency, it is assumed that selection or ranking of individual PCs is not of interest).
#'
#' @return A data frame containing these columns: 'Df', 'SS', 'RSS', 'AIC', 'F', 'p'.
#' The first row is for the full model, followed by statistics for `pcs` if present, followed by per-locus statistics for `X`.
#' If `X` had row names, these are the names of the variables (row names of this output data frame; possibly quoted with backticks if they had special symbols); otherwise 'x1', 'x2', and so on are used.
#' When there are as many or more predictors (loci) as there are samples (individuals), or if a perfect fit is detected, data returned is a dummy data frame containing only the 'p' column and with all p-values equal to 1.
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
    # when there are more predictors than samples, just silently return a vector of 1's for all p-values
    if ( m >= n )
        return( dummy_df( m, pcs = pcs ) )
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
    # however, keep existing names if present
    if ( is.null( rownames( X ) ) ) {
        rownames( X ) <- paste0( 'x', 1 : m )
    } else {
        # only restriction is names cannot contain backticks!
        indexes <- grepl( '`', rownames( X ) )
        if ( any( indexes ) )
            stop( '`X` cannot contain backticks in its rownames!  Observed cases: ', toString( rownames( X )[ indexes ] ) )
    }
    # response automatically gets called "y" here!  (no need to add names to it later)
    # NOTE: `X` gets transposed, so its rows are instead columns in `data`
    data <- cbind( y, as.data.frame( t( X ) ) )
    
    # now need to hack formula, saw example in ?formula
    # use X's names again, they must match!
    # first we have string version only, and exclude response
    # first copy down names
    formula <- rownames( X )
    # quote names so there's no problems
    # (names may have bad symbols: one example has "chr:pos" codes, which `as.formula` doesn't like!)
    formula <- paste0( '`', formula, '`' )
    # now add as sum (becomes a single string)
    formula <- paste( formula, collapse= "+" )

    if ( !is.null( pcs ) ) {
        # add PCs to data, but since we don't want these separated, have to add them this way
        data$pcs <- pcs
        # and add to formula as well, before all the actual terms of interest
        formula <- paste0( 'pcs + ', formula )
    }
    
    # add response now
    # here y name matches manually
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
        return( dummy_df( m, pcs = pcs ) )
    
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

# return a dummy data frame when there are problems with the data (perfect fits, etc)
dummy_df <- function( m, pcs = NULL ) {
    data.frame(
        p = rep.int(
            1,
            m + 1 + !is.null( pcs )
        )
    )
}
