
#' Example data
#'
#' This simulated data list is for demonstration.
#' @docType data
#' @name dat
#' @return A list containing
#' \item{known.conf}{The known confounders matrix which is adjusted in all mediation tests. Each row is a confounder, each column is a sample.}
#' \item{cov.pool}{The pool of candidate confounding variables from which potential confounders are adaptively selected to adjust for each mediation test. Each row is a covariate, each column is a sample.}
#' \item{exp.dat}{The gene expression matrix. Each row is for one gene, each column is a sample.}
#' \item{snp.dat.cis}{The cis-eQTL genotype matrix. Each row is an eQTL, each column is a sample.}
#' \item{trios.idx}{The matrix of selected trios indexes (row numbers) for mediation tests. Each row consists of the index (i.e., row number) of the eQTLs in \code{snp.dat.cis},
#' the index of cis-gene transcript in \code{exp.dat}, and the index of trans-gene transcript in \code{exp.dat}. The dimension is the number of trios by three.}
#' @examples data(example)
NULL
