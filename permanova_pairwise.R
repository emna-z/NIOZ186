#' @title PERMANOVA multiple comparisons
#'
#' @description Wrapper function for pairwise multiple comparisons using
#'     `adonis2` from package `vegan`, and adjusted p-values using `p.adjust()`.
#'
#' @param x Community `data.frame` or `dist` object.
#'
#' @param grp Vector of groups: factor levels to be compared.
#'
#' @param permutations Number of permutations, or else the permutation structure
#'     from `permute::how()`.
#'
#' @param method Dissimilarity method from `vegdist` if `x` is not a distance
#'     matrix, default is `bray`.
#'
#' @param padj The p-value correction method, one of the methods
#'     supported by `p.adjust()`; default is `bonferroni`.
#'
#' @param ... Other arguments passed to `adonis2`.
#'
#' @return Table with pairwise factors, SS, pseudo-F, R^2^, p-value and
#'     adjusted p-value.
#'
#' @author Rob Smith, inspired by Pedro Martinez Arbizu and Sylvain Monteux.
#'
#' @examples
#' # typical usage:
#' data(oakwoods, package='ecole')
#' spe <- oakwoods$spe
#' env <- oakwoods$env
#' D   <- vegan::vegdist(spe, 'bray')
#' table(env$thiltype)
#' permanova_pairwise(x = D, grp = env$thiltype)
#' # warning when any factor level contains singletons:
#' table(env$aspclass)
#' permanova_pairwise(x = D, grp = env$aspclass)
#'
#' @seealso https://github.com/pmartinezarbizu/pairwiseAdonis
#'
#' @export permanova_pairwise
#' @importFrom stats p.adjust
#' @importFrom utils combn
#' @importFrom vegan adonis2 vegdist


permanova_pairwise <- function(x,
                                 grp,
                                 permutations = 999,
                                 method = 'bray',
                                 padj = 'bonferroni', ...) {
  f     <- as.factor(grp)
  if (!all(table(f) > 1)) warning('factor has singletons! perhaps lump them?')
  co    <- combn(unique(as.character(f)),2)
  nco   <- NCOL(co)
  out   <- data.frame(matrix(NA, nrow=nco, ncol=5))
  dimnames(out)[[2]] <- c('pairs', 'SumOfSqs', 'F.Model', 'R2', 'pval')
  if (!inherits(x, 'dist')) {
    D <- vegan::vegdist(x, method=method)
  } else {
    D <- x
  }
  cat('Now performing', nco, 'pairwise comparisons. Percent progress:\n')
  for(j in 1:nco) {
    cat(round(j/nco*100,0),'...  ')
    ij  <- f %in% c(co[1,j],co[2,j])
    Dij <- as.dist(as.matrix(D)[ij,ij])
    fij <- data.frame(fij = f[ij])
    a   <- vegan::adonis2(Dij ~ fij, data=fij, permutations = permutations,
                          ...)
    out[j,1] <- paste(co[1,j], 'vs', co[2,j])
    out[j,2] <- a$SumOfSqs[1]
    out[j,3] <- a$F[1]
    out[j,4] <- a$R2[1]
    out[j,5] <- a$`Pr(>F)`[1]
  }
  cat('\n')
  out$p.adj <- p.adjust(out$pval, method=padj)
  attr(out, 'p.adjust.method') <- padj
  cat('\np-adjust method:', padj, '\n\n')
  return(out)
}