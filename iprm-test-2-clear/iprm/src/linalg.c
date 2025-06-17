
#include "linalg.h"
#include "stdio.h"
IPRMCscMatrix* new_iprm_csc_matrix(const IPRMCscMatrix* A)
{
  IPRMCscMatrix* M = iprm_malloc(sizeof(IPRMCscMatrix));

  if (A) {
    IPRMInt m = A->m;
    IPRMInt n = A->n;
    IPRMInt nnz = A->nnz;

    IPRMFloat* x = iprm_malloc(nnz * sizeof(IPRMFloat));
    IPRMInt* p = iprm_malloc((n + 1) * sizeof(IPRMInt));
    IPRMInt* i = iprm_malloc(nnz * sizeof(IPRMInt));

    copy_arrayf(A->x, x, nnz);
    copy_arrayi(A->i, i, nnz);
    copy_arrayi(A->p, p, n + 1);

    M->m = m;
    M->n = n;
    M->nnz = nnz;
    M->x = x;
    M->i = i;
    M->p = p;
  }
  else {
    M->m = 0;
    M->n = 0;
    M->nnz = 0;
    M->x = NULL;
    M->i = NULL;
    M->p = NULL;
  }

  return M;
}

IPRMCscMatrix* construct_identity(IPRMInt n, IPRMFloat lambda)
{
  IPRMCscMatrix* M = iprm_malloc(sizeof(IPRMCscMatrix));
  IPRMFloat* x = iprm_malloc(n * sizeof(IPRMFloat));
  IPRMInt* p = iprm_malloc((n + 1) * sizeof(IPRMInt));
  IPRMInt* i = iprm_malloc(n * sizeof(IPRMInt));

  M->m = n;
  M->n = n;
  M->nnz = n;
  M->x = x;
  M->i = i;
  M->p = p;

  for (IPRMInt k = 0; k < n; ++k) {
    M->i[k] = k;
    M->x[k] = lambda;
  }

  for (IPRMInt k = 0; k < n + 1; ++k) {
    M->p[k] = k;
  }

  return M;
}

void free_iprm_csc_matrix(IPRMCscMatrix* A)
{
  free(A->x);
  free(A->i);
  free(A->p);
  free(A);
}

void copy_arrayf(const IPRMFloat* x, IPRMFloat* y, IPRMInt n)
{
  iprm_assert(x || n == 0);
  iprm_assert(y || n == 0);

  for (IPRMInt i = 0; i < n; ++i) {
    y[i] = x[i];
  }
}

void copy_and_negate_arrayf(const IPRMFloat* x, IPRMFloat* y, IPRMInt n)
{
  iprm_assert(x || n == 0);
  iprm_assert(y || n == 0);

  for (IPRMInt i = 0; i < n; ++i) {
    y[i] = -x[i];
  }
}

void copy_arrayi(const IPRMInt* x, IPRMInt* y, IPRMInt n)
{
  iprm_assert(x || n == 0);
  iprm_assert(y || n == 0);

  for (IPRMInt i = 0; i < n; ++i) {
    y[i] = x[i];
  }
}

IPRMFloat dot(const IPRMFloat* u, const IPRMFloat* v, IPRMInt n)
{
  iprm_assert(u || n == 0);
  iprm_assert(v || n == 0);

  IPRMFloat x = 0.0;
  for (IPRMInt i = 0; i < n; ++i) {
    x += u[i] * v[i];
  }
  return x;
}

IPRMInt max_arrayi(const IPRMInt* x, IPRMInt n)
{
  iprm_assert(x || n == 0);

  IPRMInt max = -IPRMInt_MAX;
  for (IPRMInt i = 0; i < n; ++i) {
    max = iprm_max(max, x[i]);
  }
  return max;
}

void scale_arrayf(const IPRMFloat* x, IPRMFloat* y, IPRMFloat s, IPRMInt n)
{
  iprm_assert(x || n == 0);
  iprm_assert(y || n == 0);

  for (IPRMInt i = 0; i < n; ++i) {
    y[i] = s * x[i];
  }
}

void axpy(const IPRMFloat* x, const IPRMFloat* y, IPRMFloat* z, IPRMFloat a,
          IPRMInt n)
{
  iprm_assert(x || n == 0);
  iprm_assert(y || n == 0);

  for (IPRMInt i = 0; i < n; ++i) {
    z[i] = a * x[i] + y[i];
  }
}

void USpMv(const IPRMCscMatrix* M, const IPRMFloat* v, IPRMFloat* r)
{
  iprm_assert(M);
  iprm_assert(v);
  iprm_assert(r);

  for (IPRMInt i = 0; i < M->n; i++) {
    r[i] = 0.0;
    for (IPRMInt j = M->p[i]; j < M->p[i + 1]; j++) {
      int row = M->i[j];
      r[row] += M->x[j] * v[i];
      if (row != i)
        r[i] += M->x[j] * v[row];
    }
  }
}

void SpMv(const IPRMCscMatrix* M, const IPRMFloat* v, IPRMFloat* r)
{
  iprm_assert(M);
  iprm_assert(v);
  iprm_assert(r);

  // Clear result buffer.
  for (IPRMInt i = 0; i < M->m; ++i) {
    r[i] = 0.0;
  }

  for (IPRMInt j = 0; j < M->n; j++) {
    for (IPRMInt i = M->p[j]; i < M->p[j + 1]; i++) {
      r[M->i[i]] += M->x[i] * v[j];
    }
  }
}

void SpMtv(const IPRMCscMatrix* M, const IPRMFloat* v, IPRMFloat* r)
{
  iprm_assert(M);
  iprm_assert(v);
  iprm_assert(r);

  // Clear result buffer.
  for (IPRMInt i = 0; i < M->n; ++i) {
    r[i] = 0.0;
  }

  for (IPRMInt i = 0; i < M->n; i++) {
    for (IPRMInt j = M->p[i]; j < M->p[i + 1]; j++) {
      r[i] += M->x[j] * v[M->i[j]];
    }
  }
}

IPRMFloat inf_norm(const IPRMFloat* x, IPRMInt n)
{
  iprm_assert(x || n == 0);

  IPRMFloat norm = 0.0;
  IPRMFloat xi;
  for (IPRMInt i = 0; i < n; ++i) {
    xi = iprm_abs(x[i]);
    norm = iprm_max(norm, xi);
  }
  return norm;
}
// added for IPRM
IPRMFloat norm_2(const IPRMFloat* x, IPRMInt n)
{
  iprm_assert(x || n == 0);
  IPRMFloat norm = 0.0;
  for (IPRMInt i = 0;i < n;++i){
    norm += x[i] * x[i];
  }
  norm = iprm_sqrt(norm);
  return norm;
}
IPRMInt regularize(IPRMCscMatrix* M, IPRMFloat lambda, IPRMInt* nzadded_idx)
{

  IPRMInt num_nz = 0;
  // Iterate over each column.
  for (IPRMInt col = 0; col < M->n; col++) {
    IPRMInt start = M->p[col];
    IPRMInt end = M->p[col + 1];

    // Flag to check if the diagonal element exists.
    IPRMInt diagonal_exists = 0;

    // Iterate over the elements in the current column.
    unsigned char insert_set = 0;
    IPRMInt insert = end;
    for (IPRMInt i = start; i < end; i++) {
      if (!insert_set && M->i[i] > col) {
        insert = i;
        insert_set = 1;
      }
      if (M->i[i] == col) {
        M->x[i] += lambda; // Add lambda to the diagonal element.
        diagonal_exists = 1;
        break;
      }
    }
    // If the diagonal element does not exist, we need to insert it.
    if (!diagonal_exists) {
      // Shift all the elements in values and row_indices arrays to make space
      // for the new diagonal element.
      M->nnz++;
      M->x = realloc(M->x, M->nnz * sizeof(IPRMFloat));
      M->i = realloc(M->i, M->nnz * sizeof(IPRMInt));

      for (IPRMInt i = M->nnz - 1; i > insert; i--) {
        M->x[i] = M->x[i - 1];
        M->i[i] = M->i[i - 1];
      }

      // Insert the new diagonal element.
      M->x[insert] = lambda;
      M->i[insert] = col;
      if (nzadded_idx) {
        nzadded_idx[num_nz] = insert;
      }
      num_nz++;

      // Update the column_pointers array.
      for (IPRMInt i = col + 1; i <= M->n; i++) {
        M->p[i]++;
      }
    }
  }
  return num_nz;
}

void unregularize(IPRMCscMatrix* M, IPRMFloat lambda)
{
  // Iterate over each column.
  for (IPRMInt col = 0; col < M->n; col++) {
    IPRMInt start = M->p[col];
    IPRMInt end = M->p[col + 1];

    // Iterate over the elements in the current column.
    unsigned char insert_set = 0;
    for (IPRMInt i = start; i < end; i++) {
      if (!insert_set && M->i[i] > col) {
        insert_set = 1;
      }
      if (M->i[i] == col) {
        M->x[i] -= lambda; // Add lambda to the diagonal element.
        break;
      }
    }
  }
}

void col_inf_norm_USymm(const IPRMCscMatrix* M, IPRMFloat* norm)
{
  for (IPRMInt j = 0; j < M->n; j++) {
    for (IPRMInt idx = M->p[j]; idx < M->p[j + 1]; idx++) {
      IPRMInt row = M->i[idx];
      IPRMFloat val = iprm_abs(M->x[idx]);

      if (val > norm[j]) {
        norm[j] = val;
      }

      if (row != j) {
        if (val > norm[row]) {
          norm[row] = val;
        }
      }
    }
  }
}

void row_inf_norm(const IPRMCscMatrix* M, IPRMFloat* norm)
{
  for (IPRMInt i = 0; i < M->m; ++i) {
    IPRMFloat nrm = 0;
    IPRMFloat xi;
    for (IPRMInt nz = 0; nz < M->nnz; ++nz) {
      if (i == M->i[nz]) {
        xi = iprm_abs(M->x[nz]);
        nrm = iprm_max(nrm, xi);
      }
    }
    norm[i] = nrm;
  }
}

IPRMCscMatrix* create_transposed_matrix(const IPRMCscMatrix* A)
{
  IPRMCscMatrix* B = iprm_malloc(sizeof(IPRMCscMatrix));
  B->m = A->n;
  B->n = A->m;
  B->nnz = A->nnz;

  // Allocate memory for the transpose matrix.
  B->p = (IPRMInt*)iprm_malloc((A->m + 1) * sizeof(int));
  B->i = (IPRMInt*)iprm_malloc(A->nnz * sizeof(IPRMInt));
  B->x = (double*)iprm_malloc(A->nnz * sizeof(IPRMFloat));

  // Count the number of non-zeros in each row.
  IPRMInt* row_counts = (IPRMInt*)calloc(A->m, sizeof(IPRMInt));
  for (int j = 0; j < A->n; j++) {
    for (int i = A->p[j]; i < A->p[j + 1]; i++) {
      row_counts[A->i[i]]++;
    }
  }

  B->p[0] = 0;
  for (int i = 0; i < A->m; i++) {
    B->p[i + 1] = B->p[i] + row_counts[i];
  }

  IPRMInt* temp = (int*)calloc(
      A->m, sizeof(int)); // To track the insertion position for each row
  for (int j = 0; j < A->n; j++) {
    for (int i = A->p[j]; i < A->p[j + 1]; i++) {
      int row = A->i[i];
      int dest_pos = B->p[row] + temp[row];
      B->i[dest_pos] = j;       // Column index becomes row index
      B->x[dest_pos] = A->x[i]; // Value remains the same
      temp[row]++;
    }
  }

  // Clean up
  iprm_free(row_counts);
  iprm_free(temp);

  return B;
}

void row_col_scale(const IPRMCscMatrix* M, IPRMFloat* E, IPRMFloat* D)
{
  for (IPRMInt j = 0; j < M->n; ++j) {
    for (IPRMInt i = M->p[j]; i < M->p[j + 1]; ++i) {
      M->x[i] *= (D[j] * E[M->i[i]]);
    }
  }
}

void ew_product(IPRMFloat* x, const IPRMFloat* y, IPRMFloat* z, IPRMInt n)
{
  for (IPRMInt i = 0; i < n; ++i) {
    z[i] = x[i] * y[i];
  }
}

void invert_permutation(const IPRMInt* p, IPRMInt* pinv, IPRMInt n)
{
  for (IPRMInt i = 0; i < n; ++i) {
    pinv[p[i]] = i;
  }
}

IPRMInt cumsum(IPRMInt* p, IPRMInt* c, IPRMInt n)
{
  iprm_assert(p);
  iprm_assert(c);

  IPRMInt nz = 0;
  for (IPRMInt i = 0; i < n; i++) {
    p[i] = nz;
    nz += c[i];
    c[i] = p[i];
  }
  p[n] = nz;
  return nz;
}

IPRMCscMatrix* csc_symperm(const IPRMCscMatrix* A, const IPRMInt* pinv,
                           IPRMInt* AtoC)
{
  IPRMInt i, j, p, q, i2, j2, n;
  IPRMInt* Ap;
  IPRMInt* Ai;
  IPRMInt* Cp;
  IPRMInt* Ci;
  IPRMInt* w;
  IPRMFloat* Cx;
  IPRMFloat* Ax;
  IPRMCscMatrix* C;

  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  C = new_iprm_csc_matrix(A);
  w = iprm_calloc(n, sizeof(IPRMInt));

  iprm_assert(C);
  iprm_assert(w);

  Cp = C->p;
  Ci = C->i;
  Cx = C->x;

  for (j = 0; j < n; j++) /* count entries in each column of C */
  {
    j2 = pinv ? pinv[j] : j; /* column j of A is column j2 of C */

    for (p = Ap[j]; p < Ap[j + 1]; p++) {
      i = Ai[p];

      if (i > j)
        continue;              /* skip lower triangular part of A */
      i2 = pinv ? pinv[i] : i; /* row i of A is row i2 of C */
      w[iprm_max(i2, j2)]++;   /* column count of C */
    }
  }
  cumsum(Cp, w, n); /* compute column pointers of C */

  for (j = 0; j < n; j++) {
    j2 = pinv ? pinv[j] : j; /* column j of A is column j2 of C */

    for (p = Ap[j]; p < Ap[j + 1]; p++) {
      i = Ai[p];

      if (i > j)
        continue;              /* skip lower triangular
                                  part of A*/
      i2 = pinv ? pinv[i] : i; /* row i of A is row i2
                                  of C */
      Ci[q = w[iprm_max(i2, j2)]++] = iprm_min(i2, j2);

      if (Cx)
        Cx[q] = Ax[p];

      AtoC[p] = q;
    }
  }
  iprm_free(w);
  return C;
}
