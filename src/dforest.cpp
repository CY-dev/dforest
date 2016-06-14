#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Applic.h>
#include <map>
#include <vector>

//binary search to find the first element in vec >= value which is preassumed to exist
//similar to R function findInterval
double find_bin(double value, double *vec, R_xlen_t first, R_xlen_t last) {
  R_xlen_t index;
  while (first < last) {
    index = first + (last - first) / 2;
    if (vec[index] < value) first = index + 1;
    else last = index;
  }
  return vec[last];
}

//binary search to find the exact index of "value" in "vec"
//return -1 if "value" is not in "vec" 
R_xlen_t find_index(double value, double *vec, R_xlen_t first, R_xlen_t last) {
  R_xlen_t index;
  while (first <= last) {
    index = first + (last - first) / 2;
    if (vec[index] == value) return index;
    else if (vec[index] < value) first = index + 1;
    else last = index - 1;
  }
  return -1;
}

// update node ids in-place
extern "C" {
SEXP update_id(SEXP list_node_id, SEXP chunk, SEXP nobs_, SEXP ntrees_, SEXP nleaves_prev_, 
               SEXP list_splits, SEXP list_prev_leaves) 
{
  R_xlen_t nobs = asInteger(nobs_), ntrees = asInteger(ntrees_);
  R_xlen_t i, j, parent_index, split_var, nsp, nsp2, first, last;
  SEXP node_id, splits, prev_leaves;
  int *nleaves_prev = INTEGER(nleaves_prev_);
  double *pch = REAL(chunk), *pid, *psp, split_value;
  for (i = 0; i < ntrees; i++) {
    prev_leaves = VECTOR_ELT(list_prev_leaves, i);
    if (isNull(prev_leaves)) continue;
    node_id = VECTOR_ELT(list_node_id, i);
    splits = VECTOR_ELT(list_splits, i);
    nsp = INTEGER(getAttrib(splits, R_DimSymbol))[0];
    nsp2 = nsp * 2;
    first = nsp - nleaves_prev[i];
    last = nsp - 1;
    pid = REAL(node_id);
    psp = REAL(splits);
    for (j = 0; j < nobs; j++) {
      if (pid[j] == 0.0) continue;
      // the first column of the splits are node ids
      parent_index = find_index(pid[j], psp, first, last);
      if (parent_index == -1) {
        pid[j] = 0.0; // unsplit leaf nodes are not examined again
        continue;
      }
      // the second and the third columns of splits are the variable index and the split
      split_var = (R_xlen_t) psp[nsp + parent_index];
      split_value = psp[nsp2 + parent_index];
      //if (j == 2) Rprintf("pid = %f, split_var = %d, split_value = %f\n", pid[j], split_var, split_value);
      // left child: id * 2; right child: id * 2 + 1
      pid[j] *= 2.0;
      if (pch[split_var * nobs + j] > split_value) pid[j] += 1.0;
      //if (j == 2) Rprintf("pid[0] = %f\n", pid[j]);
    }
  }
  //return list_id;
  return R_NilValue;
}
}

SEXP table_unit(double leaf, R_xlen_t nobs, R_xlen_t nlabs, int *nlevs, R_xlen_t xindex, 
                double *pch, double *plab, SEXP id, SEXP xbins_, SEXP discretized)
{
  SEXP xbins, out;
  PROTECT(xbins = VECTOR_ELT(xbins_, xindex - 1));
  R_xlen_t i, j, offset, nxb = nlevs[xindex - 1], index_last = nxb - 1;
  PROTECT(out = allocMatrix(INTSXP, nxb, nlabs));
  int cumsum, *pout = INTEGER(out);
  double bin, *pid = REAL(id), *pxb = REAL(xbins);
  std::map<double, std::map<double, int> > counts;
  // count frequencies
  offset = xindex * nobs;
  if (LOGICAL(discretized)[xindex - 1]) {
    for (i = 0; i < nobs; i++) {
      if (pid[i] == leaf) {
        bin = find_bin(pch[offset + i], pxb, 0, index_last);
        counts[pch[i]][bin]++;
      }
    }
  } else {
    for (i = 0; i < nobs; i++) {
      if (pid[i] == leaf) counts[pch[i]][pch[offset + i]]++;
    }
  }
  // enter results into a contingency table (converted to cumulative sums)
  for (j = 0; j < nlabs; j++) {
    offset = j*nxb; cumsum = 0;
    for (i = 0; i < nxb; i++) {
      cumsum += counts[plab[j]][pxb[i]];
      pout[offset + i] = cumsum;      
    }    
  }
  UNPROTECT(2);
  return out;
}

// contingency tables for each data chunk for trees numbered from tree_begin to tree_end
extern "C" {
SEXP table_chunk(SEXP chunk, SEXP tree_begin_, SEXP tree_end_, SEXP nobs_, SEXP nvars_, SEXP nlabs_, 
                 SEXP nlevs_, SEXP nleaves_, SEXP labels, SEXP xbins, SEXP discretized, SEXP list_node_id, 
                 SEXP list_leaves, SEXP list_vars)
{
  R_xlen_t tree_begin = asInteger(tree_begin_)-1, tree_end = asInteger(tree_end_)-1, ntrees = tree_end-tree_begin+1;
  R_xlen_t nobs = asInteger(nobs_), nvars = asInteger(nvars_), nlabs = asInteger(nlabs_);
  R_xlen_t i, j, k, num_protect, nleaves_tree;
  SEXP out, out_tree, out_leaf, out_var, leaves_tree, vars_tree, id_tree, vars_leaf;
  PROTECT(out = allocVector(VECSXP, ntrees));
  num_protect = 1 + ntrees;
  double *leaves, *pch = REAL(chunk), *plab = REAL(labels); int *vars;
  int *nlevs = INTEGER(nlevs_), *nleaves = INTEGER(nleaves_);
  // compute contingency tables
  for (i = tree_begin; i <= tree_end; i++) {
    if (nleaves[i] == 0) {
      SET_VECTOR_ELT(out, i-tree_begin, R_NilValue);
      continue;      
    }
    leaves_tree = VECTOR_ELT(list_leaves, i);
    vars_tree =  VECTOR_ELT(list_vars, i);
    id_tree = VECTOR_ELT(list_node_id, i);
    nleaves_tree = nleaves[i];
    leaves = REAL(leaves_tree);
    PROTECT(out_tree = allocVector(VECSXP, nleaves_tree)); // nleaves_tree varies across trees
    for (j = 0; j < nleaves_tree; j++) {
      vars_leaf = VECTOR_ELT(vars_tree, j);
      vars = INTEGER(vars_leaf);
      PROTECT(out_leaf = allocVector(VECSXP, nvars));
      //Rprintf("leaf: %f\n", leaves[j]);
      for (k = 0; k < nvars; k++) {
        //Rprintf("var: %d\n", vars[k]);
        out_var = table_unit(leaves[j], nobs, nlabs, nlevs, vars[k], 
                             pch, plab, id_tree, xbins, discretized);
        SET_VECTOR_ELT(out_leaf, k, out_var);
      }
      SET_VECTOR_ELT(out_tree, j, out_leaf);
    }
    SET_VECTOR_ELT(out, i-tree_begin, out_tree);
    num_protect += nleaves_tree;
  }
  UNPROTECT(num_protect);
  return out;
}
}

// gini-index for all splits of one variable
double* gini_var(SEXP tab, R_xlen_t xindex, R_xlen_t nlabs, int *nlevs) {
  R_xlen_t i, j, ind_min, nxb = nlevs[xindex-1], nrow = nxb-1;
  int *pt = INTEGER(tab), N = 0; 
  double ct1, ct2, N1, N2, s1, s2, gini_min, gini_new, N1_min;
  int *totals = Calloc(nlabs, int);
  double *gini = Calloc(nrow*4, double);
  double *out = Calloc(6, double);
  for (j = 0; j < nlabs; j++) {
    totals[j] = pt[(j+1)*nxb-1];
    //Rprintf("%d\n", totals[j]);
    N += totals[j];
  }
  gini_min = nlabs*N;
  for (i = 0; i < nrow; i++) {
    N1 = 0; N2 = 0; s1 = 0; s2 = 0;
    for (j = 0; j < nlabs; j++) {
      ct1 = pt[j*nxb+i];
      ct2 = totals[j] - ct1;
      N1 += ct1;
      s1 += ct1*ct1;
      s2 += ct2*ct2;
    }
    N2 = N-N1;
    //Rprintf("N1 = %f, N2 = %f\n", N1, N2);
    if (N1) gini[i] = N1 - s1/N1;
    else gini[i] = 0.0;
    if (N2) gini[nrow+i] = N2 - s2/N2;
    else gini[nrow+i] = 0.0;
    gini_new = gini[i] + gini[nrow+i];
    if (gini_new < gini_min) {
      gini_min = gini_new;
      N1_min = N1;
      ind_min = i;
    }
  }
  out[0] = gini_min;
  out[1] = gini[ind_min];
  out[2] = gini[nrow+ind_min];
  out[3] = ind_min;
  out[4] = N;
  out[5] = N1_min;
  Free(totals);
  Free(gini);
  return out;
}

// class labels for left and right child
double* assign_label(SEXP tabs_leaf, R_xlen_t ind_minVar, R_xlen_t ind_minLev, int var_opt,
                  SEXP labels, R_xlen_t nlabs, int *nlevs) 
{
  int ct_left_max, ct_right_max, ct_left_new, ct_right_new;
  R_xlen_t i, ind_left, ind_right, nxb = nlevs[var_opt-1];
  SEXP tab = VECTOR_ELT(tabs_leaf, ind_minVar);
  int *pt = INTEGER(tab);
  double *plab = REAL(labels);
  double *out = Calloc(2, double);
  GetRNGstate();
  double incr_order = rbinom(1, 0.5);
  PutRNGstate();
  if (incr_order) { // start from the smallest label
    ct_left_max = pt[ind_minLev];
    ct_right_max = pt[nxb-1] - ct_left_max;
    ind_left = 0;
    ind_right = 0;
    for (i = 1; i < nlabs; i++) {
      ct_left_new = pt[i*nxb+ind_minLev];
      if (ct_left_new > ct_left_max) {
        ct_left_max = ct_left_new;
        ind_left = i;
      }
      ct_right_new = pt[(i+1)*nxb-1] - ct_left_new;
      if (ct_right_new > ct_right_max) {
        ct_right_max = ct_right_new;
        ind_right = i;
      }
    }
  } else { // start from the largest label
    ct_left_max = pt[(nlabs-1)*nxb+ind_minLev];
    ct_right_max = pt[nlabs*nxb-1] - ct_left_max;
    ind_left = nlabs-1;
    ind_right = nlabs-1;
    for (i = nlabs-2; i > -1; i--) {
      ct_left_new = pt[i*nxb+ind_minLev];
      if (ct_left_new > ct_left_max) {
        ct_left_max = ct_left_new;
        ind_left = i;
      }
      ct_right_new = pt[(i+1)*nxb-1] - ct_left_new;
      if (ct_right_new > ct_right_max) {
        ct_right_max = ct_right_new;
        ind_right = i;
      }
    }
  }
  out[0] = plab[ind_left];
  out[1] = plab[ind_right];
  //Rprintf("left = %f, right = %f\n", out[0], out[1]);
  return out;
}

SEXP split_leaf(double leaf, SEXP tabs_leaf, SEXP vars_leaf, SEXP splits, R_xlen_t nvars, R_xlen_t nlabs, 
                int *nlevs, SEXP labels, SEXP split_values, double init_gini_tree, R_xlen_t minNodeSize)
{
  R_xlen_t i, ind_minVar, ind_parent;
  double node, parent;
  double gini_original, gini_minVar, gini_left, gini_right, ind_minLev, nobs_leaf, nobs_left;
  int *vars = INTEGER(vars_leaf), var_opt;
  double *val = Calloc(4, double);
  SEXP tab;
  // initialize gini_original
  if (leaf == 1.0) {
    gini_original = init_gini_tree;
  } else {
    double *psp = REAL(splits);
    R_xlen_t nsp = INTEGER(getAttrib(splits, R_DimSymbol))[0];
    //Rprintf("nsp = %d\n", nsp);
    node = leaf/2.0;
    parent = floor(node);
    ind_parent = find_index(parent, psp, 0, nsp-1);
    if (node > parent) {
      gini_original = psp[4*nsp+ind_parent]; // right child
    } else {
      gini_original = psp[3*nsp+ind_parent]; // left child
    }
  }
  // find the optimal split
  tab = VECTOR_ELT(tabs_leaf, 0);
  val = gini_var(tab, vars[0], nlabs, nlevs);
  gini_minVar = val[0];
  gini_left = val[1];
  gini_right = val[2];
  ind_minLev = val[3];
  nobs_leaf = val[4]; // constant
  nobs_left = val[5];
  ind_minVar = 0;
  for (i = 1; i < nvars; i++) {
    tab = VECTOR_ELT(tabs_leaf, i);
    val = gini_var(tab, vars[i], nlabs, nlevs);
    if (val[0] < gini_minVar) {
      gini_minVar = val[0];
      gini_left = val[1];
      gini_right = val[2];
      ind_minLev = val[3];
      nobs_left = val[5];
      ind_minVar = i;
    }
    //Rprintf("gini_minVar = %f, ind_minVar = %d\n", gini_minVar, ind_minVar);
  }
  if (nobs_left < minNodeSize || nobs_leaf - nobs_left < minNodeSize) {
    Free(val);
    return R_NilValue;
  }
  SEXP out;
  PROTECT(out = allocVector(REALSXP, 9));
  double *pout = REAL(out);
  double *chlab = Calloc(2, double); // labels of left child and right child
  var_opt = vars[ind_minVar];
  chlab = assign_label(tabs_leaf, ind_minVar, (R_xlen_t) ind_minLev, var_opt, labels, nlabs, nlevs);
  pout[0] = leaf;
  pout[1] = var_opt;
  pout[2] = REAL(VECTOR_ELT(split_values, var_opt-1))[(R_xlen_t) ind_minLev];
  pout[3] = gini_left;
  pout[4] = gini_right;
  pout[5] = gini_original - gini_minVar;
  pout[6] = nobs_leaf;
  pout[7] = chlab[0];
  pout[8] = chlab[1];
  UNPROTECT(1);
  Free(val);
  Free(chlab);
  return out;
}

// update splits_tree
SEXP split_tree(SEXP splits_tree, SEXP tabs_tree, SEXP vars_tree, SEXP leaves_tree, double init_gini_tree, 
                SEXP labels, int *nlevs, SEXP split_values, R_xlen_t nvars, R_xlen_t nlabs, 
                R_xlen_t nleaves_tree, R_xlen_t minNodeSize, R_xlen_t depth)
{
  R_xlen_t i, j, jn, jnsp, nsplits = 0, nleaves_next = 0;
  SEXP out, splits_new, leaves_new, leaves_next, split, tabs_leaf, vars_leaf;
  PROTECT(out = allocVector(VECSXP, 4));
  PROTECT(splits_new = allocVector(VECSXP, nleaves_tree));
  PROTECT(leaves_new = allocVector(REALSXP, nleaves_tree));
  PROTECT(leaves_next = allocVector(REALSXP, nleaves_tree*2));
  double *leaves = REAL(leaves_tree), *pl_new = REAL(leaves_new), *pl_next = REAL(leaves_next);
  // find new splits
  for (i = 0; i < nleaves_tree; i++) {
    tabs_leaf = VECTOR_ELT(tabs_tree, i);
    vars_leaf = VECTOR_ELT(vars_tree, i);
    split = split_leaf(leaves[i], tabs_leaf, vars_leaf, splits_tree, nvars, nlabs, 
                       nlevs, labels, split_values, init_gini_tree, minNodeSize);
    if (!isNull(split)) {
      SET_VECTOR_ELT(splits_new, nsplits, split);
      pl_new[nsplits++] = leaves[i];
      pl_next[nleaves_next++] = leaves[i]*2.0;
      pl_next[nleaves_next++] = leaves[i]*2.0 + 1.0;
    }
  }
  //Rprintf("nleaves_tree = %d, nsplits = %d, nleaves_next = %d\n", nleaves_tree, nsplits, nleaves_next);
  // update splits_tree
  if (!nsplits) {
    SET_VECTOR_ELT(out, 0, splits_tree);
    SET_VECTOR_ELT(out, 1, R_NilValue);
    SET_VECTOR_ELT(out, 2, R_NilValue);
    SET_VECTOR_ELT(out, 3, ScalarInteger(0));
    UNPROTECT(4);
    return out;
  }
  if (nsplits < nleaves_tree) {
    SETLENGTH(leaves_new, nsplits);
    SETLENGTH(leaves_next, nleaves_next);
  }
  SEXP splits_upd; // updated table of splits
  R_xlen_t n, nsp = 0; double *psu, *pnew, *psp;
  //Rprintf("depth = %d\n", depth);
  if (depth == 1) {
    n = nsplits; // nsp = 0
    PROTECT(splits_upd = allocMatrix(REALSXP, n, 9));
    psu = REAL(splits_upd);
    // add new splits
    for (i = 0; i < n; i++) {
      pnew = REAL(VECTOR_ELT(splits_new, i));
      for (j = 0; j < 9; j++) psu[j * n + i] = pnew[j];
    }
  } else {
    nsp = INTEGER(getAttrib(splits_tree, R_DimSymbol))[0];
    n = nsp + nsplits;
    PROTECT(splits_upd = allocMatrix(REALSXP, n, 9));
    psu = REAL(splits_upd);
    psp = REAL(splits_tree);
    // copy old splits
    for (j = 0; j < 9; j++) {
      jn = j * n; jnsp = j * nsp;
      //Rprintf("jn = %d, jnsp = %d\n", jn, jnsp);
      for (i = 0; i < nsp; i++) {
        psu[jn + i] = psp[jnsp + i];
      }
    }
    //Rprintf("nsp = %d, nsplits = %d, n = %d\n", nsp, nsplits, n);
    // add new splits
    for (i = nsp; i < n; i++) {
      pnew = REAL(VECTOR_ELT(splits_new, i-nsp));
      for (j = 0; j < 9; j++) psu[j * n + i] = pnew[j];
    }
  }
  //Rprintf("splits_upd: nrow = %d\n", INTEGER(getAttrib(splits_upd, R_DimSymbol))[0]);
  SET_VECTOR_ELT(out, 0, splits_upd);
  SET_VECTOR_ELT(out, 1, leaves_new);
  SET_VECTOR_ELT(out, 2, leaves_next);
  SET_VECTOR_ELT(out, 3, ScalarInteger(nsplits));
  UNPROTECT(5);
  return out;
}

// update list_splits, list_prev_leaves, list_leaves, nleaves_prev, nleaves
extern "C" {
SEXP add_splits(SEXP list_splits, SEXP list_tabs, SEXP list_vars, SEXP list_prev_leaves, SEXP list_leaves, 
                SEXP init_gini, SEXP labels, SEXP split_values, SEXP tree_begin_, SEXP tree_end_, SEXP ntrees_total, 
                SEXP nvars_, SEXP nlabs_, SEXP nlevs_, SEXP nleaves_prev_, SEXP nleaves_, SEXP minNodeSize_, SEXP depth_)
{
  R_xlen_t tree_begin = asInteger(tree_begin_)-1, tree_end = asInteger(tree_end_)-1;
  R_xlen_t i, ntrees = asInteger(ntrees_total), nvars = asInteger(nvars_), nlabs = asInteger(nlabs_);
  R_xlen_t minNodeSize = asInteger(minNodeSize_), depth = asInteger(depth_);
  SEXP out, list_splits_upd, list_prev_leaves_upd, list_leaves_upd, nleaves_prev_upd, nleaves_upd;
  SEXP splits_tree, tabs_tree, vars_tree, prev_leaves_tree, leaves_tree, newval;
  PROTECT(out = allocVector(VECSXP, 5));
  PROTECT(list_splits_upd = allocVector(VECSXP, ntrees));
  PROTECT(list_prev_leaves_upd = allocVector(VECSXP, ntrees));
  PROTECT(list_leaves_upd = allocVector(VECSXP, ntrees));
  PROTECT(nleaves_prev_upd = allocVector(INTSXP, ntrees));
  PROTECT(nleaves_upd = allocVector(INTSXP, ntrees));
  double *ig = REAL(init_gini);
  int *nlevs = INTEGER(nlevs_), *nleaves_prev = INTEGER(nleaves_prev_), *nleaves = INTEGER(nleaves_);
  int *pnp = INTEGER(nleaves_prev_upd), *pn = INTEGER(nleaves_upd);
  for (i = 0; i < ntrees; i++) {
    //Rprintf("%d\n", i);
    splits_tree = VECTOR_ELT(list_splits, i);
    vars_tree = VECTOR_ELT(list_vars, i);
    prev_leaves_tree = VECTOR_ELT(list_prev_leaves, i);
    leaves_tree = VECTOR_ELT(list_leaves, i);
    if (i < tree_begin || i > tree_end || nleaves[i] == 0) {
      SET_VECTOR_ELT(list_splits_upd, i, splits_tree);
      SET_VECTOR_ELT(list_prev_leaves_upd, i, prev_leaves_tree);
      SET_VECTOR_ELT(list_leaves_upd, i, leaves_tree);
      pnp[i] = nleaves_prev[i];
      pn[i] = nleaves[i];
    } else {
      tabs_tree = VECTOR_ELT(list_tabs, i-tree_begin);
      newval = split_tree(splits_tree, tabs_tree, vars_tree, leaves_tree, ig[i], labels, 
                          nlevs, split_values, nvars, nlabs, nleaves[i], minNodeSize, depth);
      SET_VECTOR_ELT(list_splits_upd, i, VECTOR_ELT(newval, 0));
      SET_VECTOR_ELT(list_prev_leaves_upd, i, VECTOR_ELT(newval, 1));
      SET_VECTOR_ELT(list_leaves_upd, i, VECTOR_ELT(newval, 2));
      pnp[i] = asInteger(VECTOR_ELT(newval, 3));
      pn[i] = pnp[i]*2;
    }
  }
  SET_VECTOR_ELT(out, 0, list_splits_upd);
  SET_VECTOR_ELT(out, 1, list_prev_leaves_upd);
  SET_VECTOR_ELT(out, 2, list_leaves_upd);
  SET_VECTOR_ELT(out, 3, nleaves_prev_upd);
  SET_VECTOR_ELT(out, 4, nleaves_upd);
  UNPROTECT(6);
  return out;
}
}

extern "C" {
SEXP combine_ctabs(SEXP list_tabs_chunk, SEXP tree_begin_, SEXP tree_end_, 
                   SEXP nvars_, SEXP nlabs_, SEXP nleaves_)
{
  R_xlen_t tree_begin = asInteger(tree_begin_)-1, tree_end = asInteger(tree_end_)-1, ntrees = tree_end-tree_begin+1;
  R_xlen_t nchunks = xlength(list_tabs_chunk), nvars = asInteger(nvars_), nlabs = asInteger(nlabs_);
  R_xlen_t i, j, k, c, row, col, nleaves_tree, nrow, colnrow, num_protect;
  SEXP out, out_tree, out_leaf, tab;
  PROTECT(out = allocVector(VECSXP, ntrees));
  int *ptab, *pnew, *nleaves = INTEGER(nleaves_);
  num_protect = 1 + ntrees;
  for (i = 0; i < ntrees; i++) {
    if (nleaves[i] == 0) {
      SET_VECTOR_ELT(out, i, R_NilValue);
      continue;
    }
    nleaves_tree = nleaves[i+tree_begin];
    PROTECT(out_tree = allocVector(VECSXP, nleaves_tree));
    for (j = 0; j < nleaves_tree; j++) {
      PROTECT(out_leaf = allocVector(VECSXP, nvars));
      for (k = 0; k < nvars; k++) {
        tab = VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(list_tabs_chunk, 0), i), j), k);
        ptab = INTEGER(tab);
        nrow = INTEGER(getAttrib(tab, R_DimSymbol))[0]; // ncol = nlabs
        for (c = 1; c < nchunks; c++) {
          pnew = INTEGER(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(list_tabs_chunk, c), i), j), k));
          for (col = 0; col < nlabs; col++) {
            colnrow = col*nrow;
            for (row = 0; row < nrow; row++) {
              ptab[colnrow + row] += pnew[colnrow + row];
            }
          }
        }
        SET_VECTOR_ELT(out_leaf, k, tab);
      }
      SET_VECTOR_ELT(out_tree, j, out_leaf);
    }
    SET_VECTOR_ELT(out, i, out_tree);
    num_protect += nleaves_tree;
  }
  UNPROTECT(num_protect);
  return out;
}
}

// sample variable candidates for each node
extern "C" {
SEXP sample_vars(SEXP ntrees_, SEXP nvars_, SEXP nleaves_, SEXP sample_expr, SEXP rho) {
  // sample_expr = quote(sample.int(p, nvars))
  R_xlen_t i, j, nleaves_tree, ntrees = asInteger(ntrees_), nvars = asInteger(nvars_), num_protect = 1;
  int *nleaves = INTEGER(nleaves_);
  SEXP out, out_tree;
  PROTECT(out = allocVector(VECSXP, ntrees));
  for (i = 0; i < ntrees; i++) {
    nleaves_tree = nleaves[i];
    if (nleaves_tree == 0) {
      SET_VECTOR_ELT(out, i, R_NilValue);
    } else {
      num_protect++;
      PROTECT(out_tree = allocVector(VECSXP, nleaves_tree));
      for (j = 0; j < nleaves_tree; j++) {
        SET_VECTOR_ELT(out_tree, j, eval(sample_expr, rho));
      }
      SET_VECTOR_ELT(out, i, out_tree);
    }
  }
  UNPROTECT(num_protect);
  return out;
}
}

extern "C" {
SEXP pred_tree(SEXP chunk, SEXP nobs_, SEXP splits_tree, SEXP inbag_) {
  R_xlen_t i, j, nobs = asInteger(nobs_);
  R_xlen_t nsp = INTEGER(getAttrib(splits_tree, R_DimSymbol))[0], nsp2 = nsp*2, nsp5 = nsp*5, nsp6 = nsp*6, nv;
  SEXP out;
  PROTECT(out = allocVector(REALSXP, nobs));
  double *node = Calloc(nobs, double);
  double *pch = REAL(chunk), *psp = REAL(splits_tree), *inbag = REAL(inbag_), *pout = REAL(out);
  double parent, var, split, class_left, class_right;
  for (i = 0; i < nobs; i++) {
    if (inbag[i]) pout[i] = NA_REAL;
    else node[i] = 1.0;
  }
  for (j = 0; j < nsp; j++) {
    parent = psp[j];
    var = psp[nsp + j];
    split = psp[nsp2 + j];
    class_left = psp[nsp5 + j];
    class_right = psp[nsp6 + j];
    nv = (R_xlen_t) nobs * var;
    for (i = 0; i < nobs; i++) {
      if (node[i] == parent) {
        node[i] = parent*2;
        if (pch[nv + i] > split) {
          // right
          node[i] += 1.0;
          pout[i] = class_right;
        } else {
          // left
          pout[i] = class_left;
        }
      }
    }
  }
  UNPROTECT(1);
  Free(node);
  return out;
}
}

SEXP zeros(R_xlen_t nobs) {
  R_xlen_t i;
  SEXP out;
  PROTECT(out = allocVector(REALSXP, nobs));
  double *pout = REAL(out);
  for (i = 0; i < nobs; i++) pout[i] = 0.0;
  UNPROTECT(1);
  return out;
}

extern "C" {
SEXP pred_chunk(SEXP chunk, SEXP nobs_, SEXP ntrees_, SEXP labels, SEXP list_splits) {
  R_xlen_t i, j, nobs = asInteger(nobs_), ntrees = asInteger(ntrees_), nlabs = xlength(labels);
  int vote_max, vote_new; double class_max, class_new;
  std::vector<std::map<double, int> > votes(nobs);
  SEXP splits_tree, ptree, out;
  SEXP zero_col = zeros(nobs);
  PROTECT(out = allocMatrix(REALSXP, nobs, 2));
  double *pout = REAL(out), *pch = REAL(chunk), *plab = REAL(labels), *pt;
  double *largest = Calloc(nobs, double);
  // for each observation: 1-choose the largest label to break ties; 0-smallest
  GetRNGstate();
  for (i = 0; i < nobs; i++) {
    pout[i] = pch[i];
    largest[i] = rbinom(1, 0.5);
  }
  PutRNGstate();
  for (j = 0; j < ntrees; j++) {
    splits_tree = VECTOR_ELT(list_splits, j);
    ptree = pred_tree(chunk, nobs_, splits_tree, zero_col);
    pt = REAL(ptree);
    for (i = 0; i < nobs; i++) votes[i][pt[i]]++;
  }
  for (i = 0; i < nobs; i++) {
    vote_max = 0; class_max = plab[0];
    for (j = 0; j < nlabs; j++) {
      class_new = plab[j];
      vote_new = votes[i][class_new];
      if (vote_new > vote_max) {
        class_max = class_new;
        vote_max = vote_new;
      } else if (vote_new == vote_max && class_new > class_max && largest[i]) {
        class_max = class_new;
      }
    }
    pout[nobs + i] = class_max;
  }
  UNPROTECT(1);
  Free(largest);
  return out;
}
}

extern "C" {
SEXP oob_chunk(SEXP chunk, SEXP nobs_, SEXP ntrees_, SEXP list_splits, SEXP list_inbag) {
  R_xlen_t i, j, nobs = asInteger(nobs_), ntrees = asInteger(ntrees_);
  double *class_max = Calloc(nobs, double);
  double *largest = Calloc(nobs, double);
  // for each observation: 1-choose the largest label to break ties; 0-smallest
  GetRNGstate();
  for (i = 0; i < nobs; i++) {
    class_max[i] = NA_REAL;
    largest[i] = rbinom(1, 0.5);
  }
  PutRNGstate();
  int *vote_max = Calloc(nobs, int);
  std::vector<std::map<double, int> > votes(nobs);
  SEXP splits_tree, inbag_tree, ptree, out, noob, nerr;
  PROTECT(out = allocVector(VECSXP, 2));
  PROTECT(noob = allocVector(INTSXP, ntrees));
  PROTECT(nerr = allocVector(INTSXP, ntrees));
  int *pn = INTEGER(noob), *pe = INTEGER(nerr), vote_new;
  double *pch = REAL(chunk), *pt, class_new;
  for (j = 0; j < ntrees; j++) {
    splits_tree = VECTOR_ELT(list_splits, j);
    inbag_tree = VECTOR_ELT(list_inbag, j);
    ptree = pred_tree(chunk, nobs_, splits_tree, inbag_tree);
    pt = REAL(ptree);
    pn[j] = 0; pe[j] = 0;
    for (i = 0; i < nobs; i++) {
      class_new = pt[i];
      if (!ISNA(class_new)) {
        vote_new = (++votes[i][class_new]);
        if (vote_new > vote_max[i]) {
          vote_max[i] = vote_new;
          class_max[i] = class_new;
        } else if (vote_new == vote_max[i] && class_new > class_max[i] && largest[i]){ // class_max[i] != NA_REAL
          class_max[i] = class_new;
        }
      }
      // class_max[i] is the running OOB prediction, NA if i is in-bag for all trees so far
      if (!ISNA(class_max[i])) {
        pn[j]++;
        if (class_max[i] != pch[i]) pe[j]++;
      }
    }
  }
  SET_VECTOR_ELT(out, 0, noob);
  SET_VECTOR_ELT(out, 1, nerr);
  UNPROTECT(3);
  Free(class_max);
  Free(largest);
  Free(vote_max);
  return out;
}
}


static R_CallMethodDef callMethods[] = {
  {"update_id", (DL_FUNC) &update_id, 7},
  {"table_chunk", (DL_FUNC) &table_chunk, 14},
  {"add_splits", (DL_FUNC) &add_splits, 18},
  {"combine_ctabs", (DL_FUNC) &combine_ctabs, 6},
  {"sample_vars", (DL_FUNC) &sample_vars, 5},
  {"pred_tree", (DL_FUNC) &pred_tree, 4},
  {"pred_chunk", (DL_FUNC) &pred_chunk, 5},
  {"oob_chunk", (DL_FUNC) &oob_chunk, 5},
  {NULL, NULL, 0}
};

void R_init_rf(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
