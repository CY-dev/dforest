# class prediction per tree
predict_tree <- function(chunk, nobs, splits_tree, inbag_tree) {
  .Call("pred_tree", chunk, nobs, splits_tree, inbag_tree)
}

## ------------------------------------------------------------------------
# normal prediction
predict_chunk <- function(data, fit, seed = 123) {
  set.seed(seed)
  nobs = nrow(data)
  list_splits <- fit$splits
  ntrees <- length(list_splits)
  labels <- fit$labels
  results <- .Call("pred_chunk", data, nobs, ntrees, labels, list_splits)
  mean(results[,1] != results[,2])
}

## ------------------------------------------------------------------------
# Out-of-bag prediction errors
# each data chunk
oob_chunk <- function(ntrees, list_splits, sampRate, seed) {
  nobs <- nrow(chunk)
  set.seed(seed)
  list_inbag <- subsample_index(nobs, ntrees, sampRate)
  .Call("oob_chunk", chunk, nobs, ntrees, list_splits, list_inbag)
  
}

# entire data
oob <- function(wo, fit) {
  sampRate <- fit$sampRate
  seed <- fit$seed
  list_splits <- fit$splits
  ntrees <- length(list_splits)  
  # push functions to workers
  roctopus::wrun(wo, bquote({
    subsample_index <- .(subsample_index)
    oob_chunk <- .(oob_chunk)
  }))
  errors <- roctopus::wapply(wo, as.call(list(dforest::oob_chunk, ntrees, list_splits, sampRate, seed)))
  #print(errors)
  # close connections to workers
  #for(o in wo) wclose(o)
  noob <- rowSums(sapply(errors, function(v) v[[1]]))
  err <- rowSums(sapply(errors, function(v) v[[2]]))/noob
  out <- cbind(noob, err)
  colnames(out) <- c("noob", "err")
  out
}
