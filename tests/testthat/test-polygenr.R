# simulate random data
n <- 10
m <- 100
r <- 3
p <- 0.5
X <- matrix(
    rbinom( n * m, 2, p ),
    nrow = m,
    ncol = n
)

#y <- rnorm( n ) # completely random (no signal)
# to have more successful tests (where glmnet selects loci because there's signal), trait is actually given by X
m_causal <- 20
causal_indexes <- sample( m, m_causal )
causal_coeffs <- rnorm( m_causal )
y <- drop( causal_coeffs %*% X[ causal_indexes, ] )
# to mix signal and noise, first standardize
y <- drop( scale( y ) )
# then add standard normal noise
y <- y + rnorm( n )
expect_equal( length( y ), n )

# things created over the course of the tests
pcs <- NULL
obj <- NULL
obj_cv <- NULL

test_that("pca works", {
    # one mandatory argument is missing
    expect_error( pca() )
    # a successful run, this returns entire eigenvalue matrix (not filtered by r)
    expect_silent(
        pcs <- pca( X )
    )
    # test PCs
    expect_true( is.numeric( pcs ) )
    expect_true( is.matrix( pcs ) )
    expect_true( !anyNA( pcs ) )
    expect_equal( nrow( pcs ), n )
    expect_equal( ncol( pcs ), n )
    # this is PC-specific, expect orthonormality!
    expect_equal( crossprod( pcs ), diag( n ) )

    # test version where we want a smaller number of PCs
    expect_silent(
        pcs <- pca( X, r = r )
    )
    # test PCs
    expect_true( is.numeric( pcs ) )
    expect_true( is.matrix( pcs ) )
    expect_true( !anyNA( pcs ) )
    expect_equal( nrow( pcs ), n )
    expect_equal( ncol( pcs ), r )
    # this is PC-specific, expect orthonormality!
    expect_equal( crossprod( pcs ), diag( r ) )
    # NOTE: save these PCs globally!
    pcs <<- pcs
    
    # compare against `popkinsuppl::kinship_std` version, should be identical
    if ( suppressMessages( suppressWarnings( require( popkinsuppl ) ) ) ) {
        # this is the expected kinship matrix
        # (our case equals MOR version only!)
        kinship_exp <- kinship_std( X, mean_of_ratios = TRUE )
        # run this way to get kinship estimate with same subset of loci
        expect_silent(
            obj <- pca( X, p_cut = 0, want_kinship = TRUE )
        )
        kinship_obs <- obj$kinship
        # compare!
        expect_equal( kinship_obs, kinship_exp )
    }
})

# making sure global `pcs` is right
expect_true( !is.null( pcs ) )
## # test PCs
## expect_true( is.numeric( pcs ) )
## expect_true( is.matrix( pcs ) )
## expect_true( !anyNA( pcs ) )
## expect_equal( nrow( pcs ), n )
## expect_equal( ncol( pcs ), r )
## # this is PC-specific, expect orthonormality!
## expect_equal( crossprod( pcs ), diag( r ) )

test_that("glmnet_pca works", {
    # die when two mandatory parameters are missing
    expect_error( glmnet_pca() )
    expect_error( glmnet_pca( X = X ) )
    expect_error( glmnet_pca( y = y ) )

    # a successful run
    # without PCs, it should equal basic glmnet
    expect_silent(
        obj <- glmnet_pca( X, y )
    )
    # test object
    expect_equal( class( obj ), c( "elnet", "glmnet" ) )
    # there's a lot of elements here that we don't modify or edit, let's just assume they're right
    # let's look at things we do edit though (not in this case in particular, but in the pcs case, might as well do it here too)
    expect_true( 'dgCMatrix' %in% class( obj$beta ) )
    n_lambda <- length( obj$lambda ) # data-dependent value?
    expect_equal( nrow( obj$beta ), m )
    expect_equal( ncol( obj$beta ), n_lambda )
    expect_equal( obj$dim[ 1 ], m )
    expect_equal( obj$dim[ 2 ], n_lambda )

    # repeat with PCs
    expect_silent(
        obj <- glmnet_pca( X, y, pcs = pcs )
    )
    # test object
    expect_equal( class( obj ), c( "elnet", "glmnet" ) )
    # there's a lot of elements here that we don't modify or edit, let's just assume they're right
    # let's look at things we do edit though (not in this case in particular, but in the pcs case, might as well do it here too)
    expect_true( 'dgCMatrix' %in% class( obj$beta ) )
    n_lambda <- length( obj$lambda ) # data-dependent value?
    expect_equal( nrow( obj$beta ), m ) # still true because results exclude PCs
    expect_equal( ncol( obj$beta ), n_lambda )
    expect_equal( obj$dim[ 1 ], m )
    expect_equal( obj$dim[ 2 ], n_lambda )
    # save this last version globally
    obj <<- obj
})

# making sure global `obj` is right
expect_true( !is.null( obj ) )

test_that("glmnet_pca cross-validation works", {
    # a successful run
    # without PCs, it should equal basic glmnet
    # NOTE: reduced nfolds from default of 10 is required for such a tiny dataset
    expect_silent(
        obj <- glmnet_pca( X, y, cv = TRUE, nfolds = 3 )
    )
    # test object
    expect_equal( class( obj ), "cv.glmnet" )
    # there's a lot of elements here that we don't modify or edit, let's just assume they're right
    # let's look at things we do edit though (not in this case in particular, but in the pcs case, might as well do it here too)
    # test glmnet subobject, same as before!
    obj <- obj$glmnet.fit # overwrite for rest of test, for simplicity
    expect_equal( class( obj ), c( "elnet", "glmnet" ) )
    expect_true( 'dgCMatrix' %in% class( obj$beta ) )
    n_lambda <- length( obj$lambda ) # data-dependent value?
    expect_equal( nrow( obj$beta ), m )
    expect_equal( ncol( obj$beta ), n_lambda )
    expect_equal( obj$dim[ 1 ], m )
    expect_equal( obj$dim[ 2 ], n_lambda )
    
    # repeat with PCs
    expect_silent(
        obj <- glmnet_pca( X, y, pcs = pcs, cv = TRUE, nfolds = 3 )
    )
    # test object
    expect_equal( class( obj ), "cv.glmnet" )
    # save globally
    obj_cv <<- obj
    # there's a lot of elements here that we don't modify or edit, let's just assume they're right
    # let's look at things we do edit though (not in this case in particular, but in the pcs case, might as well do it here too)
    # test glmnet subobject, same as before!
    obj <- obj$glmnet.fit # overwrite for rest of test, for simplicity
    expect_equal( class( obj ), c( "elnet", "glmnet" ) )
    expect_true( 'dgCMatrix' %in% class( obj$beta ) )
    n_lambda <- length( obj$lambda ) # data-dependent value?
    expect_equal( nrow( obj$beta ), m )
    expect_equal( ncol( obj$beta ), n_lambda )
    expect_equal( obj$dim[ 1 ], m )
    expect_equal( obj$dim[ 2 ], n_lambda )
})

# making sure global `obj_cv` is right
expect_true( !is.null( obj_cv ) )

test_that( "scores_glmnet works", {
    # mandatory arguments are missing
    expect_error( scores_glmnet() )
    # a successful run
    expect_silent(
        scores <- scores_glmnet( obj$beta )
    )
    # test scores
    expect_true( is.numeric( scores ) )
    expect_true( !anyNA( scores ) )
    expect_equal( length( scores ), m )
    expect_true( min( scores ) >= 0 )
    expect_true( max( scores ) <= ncol( obj$beta ) )
    
})

test_that( "anova2 works", {
    # need a subset of genotypes smaller than sample size, so we don't overfit (this anova is an ordinary linear regression, presumably run on a small subset of loci selected by GLMNET!)
    ns <- round( n/2 )
    Xs <- X[ 1 : ns, ]

    # expected anova output column names
    names_exp <- c('Df', 'SS', 'RSS', 'AIC', 'F', 'p')
    
    # errors when mandatory arguments are missing
    expect_error( anova2() )
    expect_error( anova2( X = Xs ) )
    expect_error( anova2( y = y ) )

    # a successful run without PCs
    expect_silent(
        data <- anova2( Xs, y )
    )
    expect_true( is.data.frame( data ) )
    expect_equal( colnames( data ), names_exp )
    expect_equal( nrow( data ), ns + 1 ) # includes intercept
    # remove first row for rest of tests, since it contains several NAs always (intercept term), but other rows shouldn't
    data <- data[ -1, ]
    expect_true( !anyNA( data ) )
    # rest of names should be variables, which here are these trivial names because X had no names
    expect_equal( rownames( data ), paste0( 'x', 1 : ns ) )
    # focus on p-values
    expect_true( min( data$p ) >= 0 )
    expect_true( max( data$p ) <= 1 )

    # and a test with PCs
    expect_silent(
        data <- anova2( Xs, y, pcs )
    )
    expect_true( is.data.frame( data ) )
    expect_equal( colnames( data ), names_exp )
    expect_equal( nrow( data ), ns + 2 ) # includes intercept and PCs
    # remove first row for rest of tests, since it contains several NAs always (intercept term), but other rows shouldn't
    data <- data[ -1, ]
    expect_true( !anyNA( data ) )
    # rest of names should be variables, which here are these trivial names because X had no names
    # PLUS `pcs` must go first
    expect_equal( rownames( data ), c( 'pcs', paste0( 'x', 1 : ns ) ) )
    # focus on p-values
    expect_true( min( data$p ) >= 0 )
    expect_true( max( data$p ) <= 1 )

    # create problems on purpose
    # here we focus on special characters in variable names (rownames of X)
    # this looks like a numeric range, but also resembles chr:pos notation in real genetics applications
    rownames( Xs ) <- paste0( '1:', 1 : ns )
    expect_silent(
        data <- anova2( Xs, y, pcs )
    )
    expect_true( is.data.frame( data ) )
    expect_equal( colnames( data ), names_exp )
    expect_equal( nrow( data ), ns + 2 ) # includes intercept and PCs
    # remove first row for rest of tests, since it contains several NAs always (intercept term), but other rows shouldn't
    data <- data[ -1, ]
    expect_true( !anyNA( data ) )
    # rest of names should be variables, here use actual rownames of X
    # PLUS `pcs` must go first
    # the quotes go away when variables didn't need them, but stay when they are needed, so this comparison must be sensitive to that
    # in this case all needed quotes so that's easy
    expect_equal( rownames( data ), c( 'pcs', paste0( '`', rownames( Xs ), '`' ) ) )
    # focus on p-values
    expect_true( min( data$p ) >= 0 )
    expect_true( max( data$p ) <= 1 )

    # let's construct a wilder case
    stopifnot( ns == 5 ) # key assumption here
    rownames( Xs ) <- c( 'a+b', 'a*b', 'name with spaces', '1:10', 'x1' )
    expect_silent(
        data <- anova2( Xs, y, pcs )
    )
    expect_true( is.data.frame( data ) )
    expect_equal( colnames( data ), names_exp )
    expect_equal( nrow( data ), ns + 2 ) # includes intercept and PCs
    # remove first row for rest of tests, since it contains several NAs always (intercept term), but other rows shouldn't
    data <- data[ -1, ]
    expect_true( !anyNA( data ) )
    # rest of names should be variables, here use actual rownames of X
    # PLUS `pcs` must go first
    # the quotes go away when variables didn't need them, but stay when they are needed, so this comparison must be sensitive to that
    # in this case all needed quotes except for last
    expect_equal(
        rownames( data ),
        c(
            'pcs',
            paste0( '`', rownames( Xs )[ 1 : 4 ], '`' ),
            rownames( Xs )[ 5 ]
        )
    )
    # focus on p-values
    expect_true( min( data$p ) >= 0 )
    expect_true( max( data$p ) <= 1 )

    # cause an error by having a name with backticks:
    rownames( Xs )[ sample( ns, 1 ) ] <- 'x`'
    expect_error( anova2( Xs, y, pcs ) )
})

test_that( "anova_single works", {
    # select random loci for this test
    # (though of same length, these don't equal `causal_indexes`)
    indexes <- sample( m, m_causal )
    
    # make sure it dies with missing arguments
    expect_error( anova_single() )
    expect_error( anova_single( X = X ) )
    expect_error( anova_single( y = y ) )
    expect_error( anova_single( indexes = indexes ) )
    expect_error( anova_single( X = X, y = y ) )
    expect_error( anova_single( X = X, indexes = indexes ) )
    expect_error( anova_single( y = y, indexes = indexes ) )
    # try bad values of index
    expect_error( anova_single( X, y, indexes = -1 ) )
    expect_error( anova_single( X, y, indexes = m + 1 ) )
    
    # successful run without PCs
    expect_silent(
        scores <- anova_single( X, y, indexes )
    )
    expect_true( is.numeric( scores ) )
    expect_true( !anyNA( scores ) )
    expect_true( min( scores ) >= 0 )
    expect_equal( length( scores ), m )
    # sparsity is a big part of this data, ensure all of those scores are 0
    expect_true( all( scores[ -indexes ] == 0 ) )

    # successful run with PCs, keep them around from now on
    expect_silent(
        scores <- anova_single( X, y, indexes, pcs = pcs )
    )
    expect_true( is.numeric( scores ) )
    expect_true( !anyNA( scores ) )
    expect_true( min( scores ) >= 0 )
    expect_equal( length( scores ), m )
    # sparsity is a big part of this data, ensure all of those scores are 0
    expect_true( all( scores[ -indexes ] == 0 ) )

    # run with sparse return
    expect_silent(
        scores_sparse <- anova_single( X, y, indexes, pcs = pcs, ret_sparse = TRUE )
    )
    expect_true( is.list( scores_sparse ) )
    expect_equal( names( scores_sparse ), c('indexes', 'scores') )
    # these values should agree with our previous data/calculations!
    expect_equal( scores_sparse$indexes, indexes )
    expect_equal( scores_sparse$scores, scores[ indexes ] )
})

test_that( "anova_glmnet_single works", {
    # copy some stuff down
    beta <- obj_cv$glmnet.fit$beta
    index_min <- obj_cv$index[1]
    index_1se <- obj_cv$index[1]
    indexes_exp_min <- which( beta[ , index_min ] != 0 )
    indexes_exp_1se <- which( beta[ , index_1se ] != 0 )
    
    # make sure it dies with missing arguments
    expect_error( anova_glmnet_single() )
    expect_error( anova_glmnet_single( X = X ) )
    expect_error( anova_glmnet_single( y = y ) )
    expect_error( anova_glmnet_single( X = X, y = y ) )
    # version with `beta` shouldn't work with default `index = 'min'`)
    expect_error( anova_glmnet_single( X, y, beta = beta ) )
    expect_error( anova_glmnet_single( X, y, beta = beta, index = '1se' ) ) # ditto
    # try bad values of index
    expect_error( anova_glmnet_single( X, y, beta = beta, index = -1 ) )
    expect_error( anova_glmnet_single( X, y, beta = beta, index = 5.5 ) )
    expect_error( anova_glmnet_single( X, y, beta = beta, index = ncol( beta ) + 1 ) )

    # successful run with `obj_cv`, without PCs
    # not completely proper, since `obj_cv` used PCs, but meh
    expect_silent(
        scores <- anova_glmnet_single( X, y, obj_cv = obj_cv )
    )
    expect_true( is.numeric( scores ) )
    expect_true( !anyNA( scores ) )
    expect_true( min( scores ) >= 0 )
    expect_equal( length( scores ), m )
    # sparsity is a big part of this data, ensure all of those scores are 0
    # NOTE: using "min" index for column of `beta`
    expect_true( all( scores[ -indexes_exp_min ] == 0 ) )

    # successful run with `obj_cv`, with PCs, keep them around from now on
    expect_silent(
        scores <- anova_glmnet_single( X, y, pcs = pcs, obj_cv = obj_cv )
    )
    expect_true( is.numeric( scores ) )
    expect_true( !anyNA( scores ) )
    expect_true( min( scores ) >= 0 )
    expect_equal( length( scores ), m )
    # sparsity is a big part of this data, ensure all of those scores are 0
    # NOTE: using "min" index for column of `beta`
    expect_true( all( scores[ -indexes_exp_min ] == 0 ) )

    # run with sparse return
    expect_silent(
        scores_sparse <- anova_glmnet_single( X, y, pcs = pcs, obj_cv = obj_cv, ret_sparse = TRUE )
    )
    expect_true( is.list( scores_sparse ) )
    expect_equal( names( scores_sparse ), c('indexes', 'scores') )
    # these values should agree with our previous data/calculations!
    expect_equal( scores_sparse$indexes, indexes_exp_min )
    expect_equal( scores_sparse$scores, scores[ indexes_exp_min ] )

    # run with `index = "1se"`
    expect_silent(
        scores2 <- anova_glmnet_single( X, y, pcs = pcs, obj_cv = obj_cv, index = '1se' )
    )
    expect_true( is.numeric( scores2 ) )
    expect_true( !anyNA( scores2 ) )
    expect_true( min( scores2 ) >= 0 )
    expect_equal( length( scores2 ), m )
    # sparsity is a big part of this data, ensure all of those scores are 0
    # NOTE: using "1se" index for column of `beta`
    expect_true( all( scores2[ -indexes_exp_1se ] == 0 ) )
    
    # run with `index = "min"` but passed as number
    # expect identical results to before
    expect_silent(
        scores2 <- anova_glmnet_single( X, y, pcs = pcs, obj_cv = obj_cv, index = index_min )
    )
    expect_equal( scores2, scores )
    
    # successful run with `beta`, run with same parameters as `obj_cv` for comparison
    # expect identical results to before
    expect_silent(
        scores2 <- anova_glmnet_single( X, y, pcs = pcs, beta = beta, index = index_min )
    )
    expect_equal( scores2, scores )
})

test_that( "anova_glmnet works", {
    # copy this down
    beta <- obj$beta
    
    # errors for missing arguments
    expect_error( anova_glmnet() )
    expect_error( anova_glmnet( beta ) )
    expect_error( anova_glmnet( X = X ) )
    expect_error( anova_glmnet( y = y ) )
    expect_error( anova_glmnet( beta, X ) )
    expect_error( anova_glmnet( beta, y = y ) )
    expect_error( anova_glmnet( X = X, y = y ) )
    
    # successful run without PCs
    expect_silent(
        scores <- anova_glmnet( beta, X, y )
    )
    # test this matrix
    expect_true( 'dgCMatrix' %in% class( scores ) )
    #expect_true( is.matrix( scores ) )
    #expect_true( is.numeric( scores ) )
    expect_true( !anyNA( scores ) )
    expect_true( min( scores ) >= 0 )
    expect_equal( nrow( scores ), m )
    expect_equal( ncol( scores ), ncol( beta ) )
    # sparsity is a big part of this data, ensure all of those scores are 0
    # NOTE: suppressMessages to avoid this message: "<sparse>[ <logic> ] : .M.sub.i.logical() maybe inefficient"
    expect_true( suppressMessages( all( scores[ beta == 0 ] == 0 ) ) )
    
    # successful run with PCs
    expect_silent(
        scores <- anova_glmnet( beta, X, y, pcs = pcs )
    )
    # test this matrix
    expect_true( 'dgCMatrix' %in% class( scores ) )
    ## expect_true( is.matrix( scores ) )
    ## expect_true( is.numeric( scores ) )
    expect_true( !anyNA( scores ) )
    expect_true( min( scores ) >= 0 )
    expect_equal( nrow( scores ), m )
    expect_equal( ncol( scores ), ncol( beta ) )
    # sparsity is a big part of this data, ensure all of those scores are 0
    # NOTE: suppressMessages to avoid this message: "<sparse>[ <logic> ] : .M.sub.i.logical() maybe inefficient"
    expect_true( suppressMessages( all( scores[ beta == 0 ] == 0 ) ) )
    
})
