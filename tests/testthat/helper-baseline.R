make_complete_capture_table <- function(nlists, counts) {
  stopifnot(length(counts) == 2^nlists - 1)
  capture_matrix <- t(vapply(
    seq_len(2^nlists - 1),
    function(code) decode_capture(code, nlists),
    numeric(nlists)
  ))
  data.frame(capture_matrix, count = counts, check.names = FALSE)
}

old_abundance <- function(fit) {
  sum(fit$fit$y) + exp(unname(coef(fit$fit)[1]))
}
