# update node ids in-place
update_id <- function(list_node_id, chunk, nobs, ntrees, nleaves_prev, list_splits, list_prev_leaves) {
  .Call("update_id", list_node_id, chunk, nobs, ntrees, nleaves_prev, list_splits, list_prev_leaves)
}

# contingency tables for each chunk of data
ctab_chunk <- function(chunk, tree_begin, tree_end, nobs, nvars, nlabs, nlevs, nleaves, labels, xbins, 
                       discretized, list_node_id, list_leaves, list_vars)
{
  .Call("table_chunk", chunk, tree_begin, tree_end, nobs, nvars, nlabs, nlevs, nleaves, labels, xbins, 
        discretized, list_node_id, list_leaves, list_vars)
}

# combine contingency tables for different chunks into one
combine_ctabs <- function(list_tabs_chunk, tree_begin, tree_end, nvars, nlabs, nleaves) {
  .Call("combine_ctabs", list_tabs_chunk, tree_begin, tree_end, nvars, nlabs, nleaves)
}

# find optimal splits for each leaf of each tree
add_splits <- function(list_splits, list_tabs, list_vars, list_prev_leaves, list_leaves, init_gini, labels, 
                       split_values, tree_begin, tree_end, ntrees_total, nvars, nlabs, nlevs, nleaves_prev, 
                       nleaves, minNodeSize, depth)
{
  .Call("add_splits", list_splits, list_tabs, list_vars, list_prev_leaves, list_leaves, init_gini, labels, split_values,
        tree_begin, tree_end, ntrees_total, nvars, nlabs, nlevs, nleaves_prev, nleaves, minNodeSize, depth)
}

# sample variable candidates
sample_vars <- function(vind_nonconst, num_nonconst, ntrees, p, nvars, nleaves) {
  if (num_nonconst < p) {
    sample_expr <- quote(sample(vind_nonconst, nvars))
  } else {
    sample_expr <- quote(sample.int(p, nvars))
  }
  .Call("sample_vars", ntrees, nvars, nleaves, sample_expr, environment())
}