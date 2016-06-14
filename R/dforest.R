# Build a random forest for classification for numeric data
# k classes denoted by integers
# continuous predictors
# distributed storage of data

# tree nodes are labeled by the level order for a complete binary tree
# example:
#      1
#    2   3
#   4 5 6 7
# if a split happens at node i, its child nodes are 2i and 2i+1.

# we create a matrix "splits" to record the splits with the following columns:
# node, var, split, left child impurity, right child impurity
# a split always put samples with values less than the splitting value to the left child node

## ------------------------------------------------------------------------
# given one data chunk and rate of subsampling
# create a 1/0 vector to mark indices for one subsample
# initial node_id where 1 stands for the root node and 0 for out-of-bag
subsample_index = function(n, ntrees, sampRate) {
  m <- round(n*sampRate)
  lapply(seq_len(ntrees), function(o) {ind = sample.int(n, m); v = double(n); v[ind] = 1; v})
}

## ------------------------------------------------------------------------
# list_node_id, chunk, nobs, ntrees, nvars, nlabs, nlevs, labels, xbins, discretized are pre-stored at each worker
# update node id's for each tree in-place by the current terminal nodes
update_node_id <- function(nleaves_prev, list_splits, list_prev_leaves) {
  dforest::update_id(list_node_id, chunk, nobs, ntrees, nleaves_prev, list_splits, list_prev_leaves)
}

# compute contingency tables for one data chunk for all trees
ct_chunk <- function(tree_begin, tree_end, nleaves, list_leaves, list_vars) {
  dforest::ctab_chunk(chunk, tree_begin, tree_end, nobs, nvars, nlabs, nlevs, nleaves, labels, xbins, 
                      discretized, list_node_id, list_leaves, list_vars)
}

## ------------------------------------------------------------------------
# get bins/levels for response labels and predictors
# discretized bins for predictors with raw levels > maxNumBins
get_bins <- function(ncols, discretized, ignored) {
  probs <- 1:1000/1000
  bins <- lapply(seq_len(ncols), 
                 function(i) {
                   if (ignored[i]) {
                     NULL
                   } else if (discretized[i]) {
                     unique(quantile(chunk[,i], probs))
                   } else {
                     unique(chunk[,i])
                   }
                 })
  bins
}

combine_bins <- function(bins, list_bins, ncols, maxNumBins, discretized, ignored) {
  probs <- seq_len(maxNumBins)/maxNumBins
  for (i in seq_len(ncols)) {
    if (ignored[i]) next
    tmp <- do.call("c", lapply(list_bins, function(o) o[[i]]))
    tmp <- unique(tmp)
    # quantiles of the quantiles for all chunks
    if (discretized[i]) {
      bins[[i]] <- unique(quantile(tmp, probs))
    } else {
      bins[[i]] <- sort(tmp)
    }
  }
  bins
}

dforest <- function(wo, ntrees, nvars, sampRate = 1, maxDepth = 5, minNodeSize = 1, 
                    maxNumBins = 50, memory_limit = 4, seed=123)
{
  set.seed(seed)
  # discretization
  list_nlevs <- roctopus::wqapply(wo, apply(chunk, 2, function(o) length(unique(o))))
  nchunks <- length(list_nlevs)
  max_nlevs <- do.call("pmax", list_nlevs)
  ncols <- length(max_nlevs)
  p <- ncols - 1
  if (missing(nvars)) nvars <- floor(sqrt(p))
  # guess the number of levels per variable based on max_nlevs
  d1 <- (max_nlevs > maxNumBins)
  d1[1] <- FALSE # labels
  ignored <- logical(ncols)
  bins <- lapply(seq_len(ncols), function(o) NULL)
  # first pass based on max_nlevs
  list_bins <- roctopus::wapply(wo, as.call(list(get_bins, ncols, d1, ignored)))
  bins <- combine_bins(bins, list_bins, ncols, maxNumBins, d1, ignored)
  nlevs <- sapply(bins, length)
  # second pass based on nlevs to guarantee number of bins for each variable <= maxNumBins
  d2 <- (nlevs > maxNumBins)
  d2[1] <- FALSE # labels
  ignored[!d2] <- TRUE
  list_bins <- roctopus::wapply(wo, as.call(list(get_bins, ncols, d2, ignored)))
  bins <- combine_bins(bins, list_bins, ncols, maxNumBins, d2, ignored)
  nlevs <- sapply(bins, length)
  rm(list_bins)
  labels <- bins[[1]]
  xbins <- bins[-1]
  nlabs <- nlevs[1]
  nlevs <- nlevs[-1]
  nlevs_max <- max(nlevs)
  discretized <- (d1[-1] | d2[-1])
  vind_nonconst <- which(nlevs > 1)
  num_nonconst <- length(vind_nonconst)
  if (num_nonconst == 0) stop("All predictors are constant!")
  allvars <- (num_nonconst == p)
  #print(discretized)
  split_values <- lapply(seq_len(p), 
                         function(i) {
                           tmp = xbins[[i]]
                           v = tmp[-length(tmp)]
                           if (discretized[i]) v
                           else (v + tmp[-1])/2
                         })
  # push functions and constants to workers
  roctopus::wrun(wo, bquote({
    subsample_index <- .(subsample_index)
    update_node_id <- .(update_node_id)
    ct_chunk <- .(ct_chunk)
    ntrees <- .(ntrees)
    nvars <- .(nvars)
    labels <- .(labels)
    xbins <- .(xbins)
    nlabs <- .(nlabs)
    nlevs <- .(nlevs)
    discretized <- .(discretized)
    sampRate <- .(sampRate)
    seed <- .(seed)
    nobs <- nrow(chunk)
    # indices for the subsamples
    set.seed(seed)
    list_node_id <- subsample_index(nobs, ntrees, sampRate)
  }))
  list_splits <- lapply(seq_len(ntrees), function(o) NULL)
  list_prev_leaves <- list_splits
  list_leaves <- lapply(seq_len(ntrees), function(o) 1)
  init_gini <- double(ntrees)
  nleaves <- integer(ntrees) + 1L
  nleaves_prev <- integer(ntrees)
  # memory limit in gigabytes
  memory_splits <- 8*(2^maxDepth)*9*ntrees/(2^30)
  memory_limit <- memory_limit - memory_splits
  for (depth in seq_len(maxDepth)) {
    list_vars <- sample_vars(vind_nonconst, num_nonconst, ntrees, p, nvars, nleaves)
    roctopus::wapply(wo, as.call(list(update_node_id, nleaves_prev, list_splits, list_prev_leaves))) # update in-place
    nleaves_max <- max(nleaves)
    memory_tabs <- 4*nlevs_max*nlabs*nvars*nleaves_max*ntrees*(nchunks+1)/(2^30)
    ngroups <- ceiling(memory_tabs/memory_limit)
    tree_ids <- floor(ntrees/ngroups)*(0:ngroups)
    if (tree_ids[ngroups+1] < ntrees) {
      tree_ids <- c(tree_ids, ntrees)
      ngroups <- ngroups + 1
    }
    if (ngroups > 1) cat("Divide trees into", ngroups, "groups at depth", depth, "\n")
    for (i in seq_len(ngroups)) {
      # when ngroups = 1, tree_begin = 1 and tree_end = ntrees
      tree_begin <- tree_ids[i]+1
      tree_end <- tree_ids[i+1]
      #print(list_vars)
      # update list_node_id in each chunk
      # distributed computation for contingency tables
      list_tabs_chunk <- roctopus::wapply(wo, as.call(list(ct_chunk, tree_begin, tree_end,
                                                           nleaves, list_leaves, list_vars)))
      print(object.size(list_tabs_chunk), units = "GB")
      # sum up counts over all chunks
      list_tabs <- combine_ctabs(list_tabs_chunk, list_leaves, tree_begin, tree_end, 
                                 nvars, nlabs, nleaves)
      #rm(list_tabs_chunk)
      if (depth == 1) {
        for (j in tree_begin:tree_end) {
          tab <- list_tabs[[j-tree_begin+1]][[1]][[1]]
          totals <- tab[nrow(tab),]
          N <- sum(totals)
          init_gini[j] <- N-sum(totals^2)/N
        }
      }
      newvals <- add_splits(list_splits, list_tabs, list_vars, list_prev_leaves, list_leaves, init_gini, 
                            labels, split_values, tree_begin, tree_end, ntrees, nvars, nlabs, nlevs, 
                            nleaves_prev, nleaves, minNodeSize, depth)
      list_splits <- newvals[[1]]
      list_prev_leaves <- newvals[[2]]
      list_leaves <- newvals[[3]]
      nleaves_prev <- newvals[[4]]
      nleaves <- newvals[[5]]
      rm(list_tabs, newvals)
      #print(list_splits)
      #print(list_prev_leaves)
    }
  }
  # close connections to workers
  #for(o in wo) wclose(o)
  splits <- lapply(list_splits, function(o) {
    rownames(o) <- NULL
    colnames(o) <- c('node', 'var', 'split', 'impurity_left', 'impurity_right',
                     'improve', 'count', 'class_left', 'class_right')
    o[,-(4:5),drop=FALSE]
  })
  list(sampRate = sampRate, seed = seed, splits = splits, labels = labels)
}

## ------------------------------------------------------------------------
hdfs2ro <- function(hdfsPath, sep = "|", type = "numeric") {
  require(iotools)
  require(hmr)
  #Sys.setenv(HADOOP_PREFIX="/data/hdp2/hadoop")
  hmr(hinput(hdfsPath, function(x) mstrsplit(x, sep = sep, type = type)),
      map = function(m) {
        chunk <<- m
        roctopus::defer()})
}
