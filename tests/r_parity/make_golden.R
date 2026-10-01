# Generate golden values for tests/test_r_parity.py using the ORIGINAL R functions.
# Run from the repo root:  Rscript tests/r_parity/make_golden.R
suppressMessages(library(Matrix))
source("tests/r_parity/original_functions.R")

set.seed(1)
genes <- sprintf("G%02d", 1:30)
make_layer <- function(pool, n_edges) {
  el <- data.frame(A = sample(pool, n_edges, TRUE), B = sample(pool, n_edges, TRUE),
                   stringsAsFactors = FALSE)
  el <- el[el$A != el$B, ]
  swap <- el$A > el$B
  tmp <- el$A[swap]; el$A[swap] <- el$B[swap]; el$B[swap] <- tmp
  el[!duplicated(el), ]
}
el <- list(
  L1 = make_layer(genes[1:25], 45),
  L2 = make_layer(genes[5:30], 40),
  L3 = make_layer(genes[c(1:20, 26:30)], 35)
)
for (nm in names(el)) {
  write.table(el[[nm]], sprintf("tests/r_parity/toy/%s.tsv", nm), sep = "\t",
              quote = FALSE, row.names = FALSE)
}

allnodes <- get_allnodes(el)
Mlist <- lapply(el, function(x) transitional_matrix(x, allnodes))
seeds <- c("G02", "G07", "G11", "G27")

run <- function(z, label) {
  pmat <- pmat_cal(z)
  S <- supratransitional(Mlist, pmat)
  p0 <- rep(as.numeric(allnodes %in% seeds), length(el))
  p0 <- p0 / sum(p0)
  p <- as.numeric(RWR(M = S, p_0 = p0, r = 0.7))
  P <- matrix(p, ncol = length(el), dimnames = list(allnodes, names(el)))
  out <- data.frame(config = label, gene = allnodes, P, avg = rowMeans(P))
  write.table(pmat, sprintf("tests/r_parity/pmat_%s.tsv", label), sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  out
}
res <- rbind(run(c(2.5, 4.0, 1.8), "weighted"), run(c(1, 1, 1), "uniform"))
write.table(res, "tests/r_parity/golden_rwr.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("seeds:", seeds, "\n")
