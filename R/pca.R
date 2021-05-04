# naive PCA implementation, for internal tests only
# implements a useful MAF filter, makes estimates more stable
# (p_cut = 0 is old behavior, as fixed loci always have to be removed)
pca <- function( X, r = NA, p_cut = 0.05, want_kinship = FALSE ) {
    if ( missing( X ) )
        stop( '`X` is required!' )
    
    # get ancestral allele frequency estimates
    p_hat <- rowMeans( X ) / 2
    # remove low MAF loci
    # these are the ones to keep
    indexes <- p_cut < p_hat & p_hat < 1 - p_cut
    X <- X[ indexes, ]
    p_hat <- p_hat[ indexes ]
    # get number of rows now, for normalization
    m <- nrow( X )
    
    # now center and scale matrix as is usual
    X <- ( X - 2 * p_hat ) / sqrt( 4 * p_hat * ( 1 - p_hat ) )
    # now calculate kinship estimate
    kinship <- crossprod( X ) / m

    # now get PCs
    # decompose it...
    evd <- eigen(kinship)
    # eigenvectors are PCs
    pcs <- evd$vectors
    # truncate list if desired (most practical, but not default)
    if ( !is.na( r ) )
        pcs <- pcs[ , 1 : r ]
    
    # return kinship too for checks against `popkinsuppl`
    if ( want_kinship ) {
        return(
            list(
                kinship = kinship,
                pcs = pcs
            )
        )
    }
    # otherwise just return PCs (normal usage)
    return( pcs )
}
