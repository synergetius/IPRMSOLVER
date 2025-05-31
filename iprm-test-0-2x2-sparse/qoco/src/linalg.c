
#include "linalg.h"
#include "stdio.h"
MtMPipe* constructMtMPipe(QOCOCscMatrix* M){
  // 尝试优化一下
  MtMPipe* pipe = qoco_malloc(sizeof(MtMPipe));
  pipe->M = M;
  pipe->MtM = qoco_malloc(sizeof(QOCOCscMatrix));
  QOCOCscMatrix* MtM = pipe->MtM;
  MtM->m = M->n;
  MtM->n = M->n;
  // 统计MtM->nnz，分配内存
  MtM->nnz = 0;
  for (int j = 0;j < M->n;j++){
    for (int i = 0;i <= j;i++){
      int ei = M->p[i], ej = M->p[j];
      for (;ei < M->p[i + 1];ei++){
        while(ej < M->p[j + 1] && M->i[ej] < M->i[ei]){
          ej++;
        }
        if (ej < M->p[j + 1] && M->i[ej] == M->i[ei]){
          MtM->nnz++;
          break;
        }
      }
    }
  }
  MtM->x = qoco_calloc(MtM->nnz, sizeof(QOCOFloat));
  MtM->i = qoco_calloc(MtM->nnz, sizeof(QOCOInt));
  MtM->p = qoco_calloc((MtM->n + 1), sizeof(QOCOInt));
  pipe->ie = qoco_calloc(MtM->nnz, sizeof(QOCOInt*));
  pipe->je = qoco_calloc(MtM->nnz, sizeof(QOCOInt*));
  pipe->cnt = qoco_calloc(MtM->nnz, sizeof(QOCOInt));
  int e = 0;
  for (int j = 0;j < M->n;j++){
    MtM->p[j] = e;
    for (int i = 0;i <= j;i++){
      int ei = M->p[i], ej = M->p[j];
      int cnt = 0;
      for (;ei < M->p[i + 1];ei++){
        while(ej < M->p[j + 1] && M->i[ej] < M->i[ei]){
          ej++;
        }
        if (ej < M->p[j + 1] && M->i[ej] == M->i[ei]){
          cnt++;
        }
      }
      if (cnt > 0){
        MtM->i[e] = i;
        pipe->cnt[e] = cnt;
        pipe->ie[e] = qoco_calloc(cnt, sizeof(QOCOInt));
        pipe->je[e] = qoco_calloc(cnt, sizeof(QOCOInt));
        e++;
      }
    }
  }
  MtM->p[M->n] = e;
  // 填充ie, je
  for (int j = 0;j < M->n;j++){
    for (int e = MtM->p[j];e < MtM->p[j + 1];e++){
      int i = MtM->i[e];
      int ei = M->p[i], ej = M->p[j];
      int cnt = 0;
      for (;ei < M->p[i + 1];ei++){
        while(ej < M->p[j + 1] && M->i[ej] < M->i[ei]){
          ej++;
        }
        if (ej < M->p[j + 1] && M->i[ej] == M->i[ei]){
          pipe->ie[e][cnt] = ei;
          pipe->je[e][cnt] = ej;
          cnt++;
        }
      }
    }
  }
  return pipe;
}
MtMPipe* constructMtMPipe0(QOCOCscMatrix* M){
  MtMPipe* pipe = qoco_malloc(sizeof(MtMPipe));
  pipe->M = M;
  pipe->MtM = qoco_malloc(sizeof(QOCOCscMatrix));
  QOCOCscMatrix* MtM = pipe->MtM;
  MtM->m = M->n;
  MtM->n = M->n;
  // 统计MtM->nnz，分配内存
  MtM->nnz = 0;
  for (int j = 0;j < M->n;j++){
    for (int i = 0;i <= j;i++){
      int ei = M->p[i], ej = M->p[j];
      int nonzero = 0;
      for (int t = 0;t < M->m;t++){
        if (
          ei < M->p[i + 1] && M->i[ei] == t &&
          ej < M->p[j + 1] && M->i[ej] == t
        ){
          nonzero = 1;
          break;
        }
        if (ei < M->p[i + 1] && M->i[ei] == t){
          ei++;
        }
        if (ej < M->p[j + 1] && M->i[ej] == t){
          ej++;
        }
      }
      if (nonzero > 0){
          MtM->nnz++;
      }
    }
  }
  MtM->x = qoco_calloc(MtM->nnz, sizeof(QOCOFloat));
  MtM->i = qoco_calloc(MtM->nnz, sizeof(QOCOInt));
  MtM->p = qoco_calloc((MtM->n + 1), sizeof(QOCOInt));
  pipe->ie = qoco_calloc(MtM->nnz, sizeof(QOCOInt*));
  pipe->je = qoco_calloc(MtM->nnz, sizeof(QOCOInt*));
  pipe->cnt = qoco_calloc(MtM->nnz, sizeof(QOCOInt));
  // 构建MtM的稀疏结构并记录pipe->cnt，分配pipe->ie, je的次级内存
  int e = 0;
  for (int j = 0;j < M->n;j++){
    MtM->p[j] = e;
    for (int i = 0;i <= j;i++){
      int ei = M->p[i], ej = M->p[j];
      int cnt = 0;
      for (int t = 0;t < M->m;t++){
        if (
          ei < M->p[i + 1] && M->i[ei] == t &&
          ej < M->p[j + 1] && M->i[ej] == t
        ){
          cnt++;
        }
        if (ei < M->p[i + 1] && M->i[ei] == t){
          ei++;
        }
        if (ej < M->p[j + 1] && M->i[ej] == t){
          ej++;
        }
      }
      if (cnt > 0){
        MtM->i[e] = i;
        pipe->cnt[e] = cnt;
        pipe->ie[e] = qoco_calloc(cnt, sizeof(QOCOInt));
        pipe->je[e] = qoco_calloc(cnt, sizeof(QOCOInt));
        e++;
      }
    }
  }
  MtM->p[M->n] = e;
  // 填充ie, je
  for (int j = 0;j < M->n;j++){
    for (int e = MtM->p[j];e < MtM->p[j + 1];e++){
      int i = MtM->i[e];
      int ei = M->p[i], ej = M->p[j];
      int cnt = 0;
      for (int t = 0;t < M->m;t++){
        if (
          ei < M->p[i + 1] && M->i[ei] == t &&
          ej < M->p[j + 1] && M->i[ej] == t
        ){
          pipe->ie[e][cnt] = ei;
          pipe->je[e][cnt] = ej;
          cnt++;
        }
        if (ei < M->p[i + 1] && M->i[ei] == t){
          ei++;
        }
        if (ej < M->p[j + 1] && M->i[ej] == t){
          ej++;
        }
      }
    }
  }
  return pipe;
}
void freeMtMPipe(MtMPipe* pipe){
  for (int e = 0;e < pipe->MtM->nnz;e++){
    free(pipe->ie[e]);
    free(pipe->je[e]);
  }
  free(pipe->ie);
  free(pipe->je);
  free(pipe->cnt);
  free_qoco_csc_matrix(pipe->MtM);
  // M 是外来的矩阵，不在这里free
}
void goMtMPipe(MtMPipe* pipe){
  for (int e = 0;e < pipe->MtM->nnz;e++){
    pipe->MtM->x[e] = 0;
    for (int c = 0;c < pipe->cnt[e];c++){
      int ei = pipe->ie[e][c], ej = pipe->je[e][c];
      pipe->MtM->x[e] += pipe->M->x[ei] * pipe->M->x[ej];
    }
  }
}
ApkBPipe* constructApkBPipe(QOCOCscMatrix* A, QOCOCscMatrix* B){
  /// 注意此处的ApBPipe是两个相同形状的上三角矩阵相加
  ApkBPipe* pipe = qoco_malloc(sizeof(ApkBPipe));
  pipe->A = A;
  pipe->B = B;
  pipe->ApkB = qoco_malloc(sizeof(QOCOCscMatrix));
  QOCOCscMatrix* ApkB = pipe->ApkB;
  //assert AB形状相同
  ApkB->m = A->m;
  ApkB->n = A->n;
  pipe->AtoApkB = qoco_calloc(A->nnz, sizeof(QOCOInt));
  pipe->BtoApkB = qoco_calloc(B->nnz, sizeof(QOCOInt));
  
  // 统计ApB->nnz，分配内存
  int nz = 0;
  for (int j = 0;j < ApkB->n;j++){
    int eA = A->p[j], eB = B->p[j];
    for (int i = 0;i <= j;i++){ ///////
      int nonzero = 0;
      if (eA < A->p[j + 1] && A->i[eA] == i){
        nonzero++;
        eA++;
      }
      if (eB < B->p[j + 1] && B->i[eB] == i){
        nonzero++;
        eB++;
      }
      if (nonzero > 0){
        nz++;
      }
    }
  }
  ApkB->nnz = nz;
  ApkB->x = qoco_calloc(ApkB->nnz, sizeof(QOCOFloat));
  ApkB->i = qoco_calloc(ApkB->nnz, sizeof(QOCOInt));
  ApkB->p = qoco_calloc((ApkB->n + 1), sizeof(QOCOInt));
  
  // 构造ApB稀疏结构和映射
  nz = 0;
  for (int j = 0;j < ApkB->n;j++){
    ApkB->p[j] = nz; 
    int eA = A->p[j], eB = B->p[j];
    for (int i = 0;i <= j;i++){ ///////
      int nonzero = 0;
      if (eA < A->p[j + 1] && A->i[eA] == i){
        pipe->AtoApkB[eA] = nz;
        nonzero++;
        eA++;
      }
      if (eB < B->p[j + 1] && B->i[eB] == i){
        pipe->BtoApkB[eB] = nz;
        nonzero++;
        eB++;
      }
      if (nonzero > 0){
        ApkB->i[nz] = i;
        nz++;
      }
    }
  }
  ApkB->p[ApkB->n] = nz;
  return pipe;
}
void goApkBPipe(ApkBPipe* pipe, QOCOFloat k){
  for (int e = 0;e < pipe->ApkB->nnz;e++){
    pipe->ApkB->x[e] = 0;
  }
  for (int e = 0;e < pipe->A->nnz;e++){
    pipe->ApkB->x[pipe->AtoApkB[e]] += pipe->A->x[e];
  }
  for (int e = 0;e < pipe->B->nnz;e++){
    pipe->ApkB->x[pipe->BtoApkB[e]] += k * pipe->B->x[e];
  }
}
void freeApkBPipe(ApkBPipe* pipe){
  free(pipe->AtoApkB);
  free(pipe->BtoApkB);
  free_qoco_csc_matrix(pipe->ApkB);
  //A, B是外来的矩阵，不在此处释放
}
QOCOCscMatrix* new_qoco_csc_matrix(const QOCOCscMatrix* A)
{
  QOCOCscMatrix* M = qoco_malloc(sizeof(QOCOCscMatrix));

  if (A) {
    QOCOInt m = A->m;
    QOCOInt n = A->n;
    QOCOInt nnz = A->nnz;

    QOCOFloat* x = qoco_malloc(nnz * sizeof(QOCOFloat));
    QOCOInt* p = qoco_malloc((n + 1) * sizeof(QOCOInt));
    QOCOInt* i = qoco_malloc(nnz * sizeof(QOCOInt));

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

QOCOCscMatrix* construct_identity(QOCOInt n, QOCOFloat lambda)
{
  QOCOCscMatrix* M = qoco_malloc(sizeof(QOCOCscMatrix));
  QOCOFloat* x = qoco_malloc(n * sizeof(QOCOFloat));
  QOCOInt* p = qoco_malloc((n + 1) * sizeof(QOCOInt));
  QOCOInt* i = qoco_malloc(n * sizeof(QOCOInt));

  M->m = n;
  M->n = n;
  M->nnz = n;
  M->x = x;
  M->i = i;
  M->p = p;

  for (QOCOInt k = 0; k < n; ++k) {
    M->i[k] = k;
    M->x[k] = lambda;
  }

  for (QOCOInt k = 0; k < n + 1; ++k) {
    M->p[k] = k;
  }

  return M;
}

void free_qoco_csc_matrix(QOCOCscMatrix* A)
{
  free(A->x);
  free(A->i);
  free(A->p);
  free(A);
}

void copy_arrayf(const QOCOFloat* x, QOCOFloat* y, QOCOInt n)
{
  qoco_assert(x || n == 0);
  qoco_assert(y || n == 0);

  for (QOCOInt i = 0; i < n; ++i) {
    y[i] = x[i];
  }
}

void copy_and_negate_arrayf(const QOCOFloat* x, QOCOFloat* y, QOCOInt n)
{
  qoco_assert(x || n == 0);
  qoco_assert(y || n == 0);

  for (QOCOInt i = 0; i < n; ++i) {
    y[i] = -x[i];
  }
}

void copy_arrayi(const QOCOInt* x, QOCOInt* y, QOCOInt n)
{
  qoco_assert(x || n == 0);
  qoco_assert(y || n == 0);

  for (QOCOInt i = 0; i < n; ++i) {
    y[i] = x[i];
  }
}

QOCOFloat dot(const QOCOFloat* u, const QOCOFloat* v, QOCOInt n)
{
  qoco_assert(u || n == 0);
  qoco_assert(v || n == 0);

  QOCOFloat x = 0.0;
  for (QOCOInt i = 0; i < n; ++i) {
    x += u[i] * v[i];
  }
  return x;
}

QOCOInt max_arrayi(const QOCOInt* x, QOCOInt n)
{
  qoco_assert(x || n == 0);

  QOCOInt max = -QOCOInt_MAX;
  for (QOCOInt i = 0; i < n; ++i) {
    max = qoco_max(max, x[i]);
  }
  return max;
}

void scale_arrayf(const QOCOFloat* x, QOCOFloat* y, QOCOFloat s, QOCOInt n)
{
  qoco_assert(x || n == 0);
  qoco_assert(y || n == 0);

  for (QOCOInt i = 0; i < n; ++i) {
    y[i] = s * x[i];
  }
}

void axpy(const QOCOFloat* x, const QOCOFloat* y, QOCOFloat* z, QOCOFloat a,
          QOCOInt n)
{
  qoco_assert(x || n == 0);
  qoco_assert(y || n == 0);

  for (QOCOInt i = 0; i < n; ++i) {
    z[i] = a * x[i] + y[i];
  }
}

void USpMv(const QOCOCscMatrix* M, const QOCOFloat* v, QOCOFloat* r)
{
  qoco_assert(M);
  qoco_assert(v);
  qoco_assert(r);

  for (QOCOInt i = 0; i < M->n; i++) {
    r[i] = 0.0;
    for (QOCOInt j = M->p[i]; j < M->p[i + 1]; j++) {
      int row = M->i[j];
      r[row] += M->x[j] * v[i];
      if (row != i)
        r[i] += M->x[j] * v[row];
    }
  }
}

void SpMv(const QOCOCscMatrix* M, const QOCOFloat* v, QOCOFloat* r)
{
  qoco_assert(M);
  qoco_assert(v);
  qoco_assert(r);

  // Clear result buffer.
  for (QOCOInt i = 0; i < M->m; ++i) {
    r[i] = 0.0;
  }

  for (QOCOInt j = 0; j < M->n; j++) {
    for (QOCOInt i = M->p[j]; i < M->p[j + 1]; i++) {
      r[M->i[i]] += M->x[i] * v[j];
    }
  }
}

void SpMtv(const QOCOCscMatrix* M, const QOCOFloat* v, QOCOFloat* r)
{
  qoco_assert(M);
  qoco_assert(v);
  qoco_assert(r);

  // Clear result buffer.
  for (QOCOInt i = 0; i < M->n; ++i) {
    r[i] = 0.0;
  }

  for (QOCOInt i = 0; i < M->n; i++) {
    for (QOCOInt j = M->p[i]; j < M->p[i + 1]; j++) {
      r[i] += M->x[j] * v[M->i[j]];
    }
  }
}

QOCOFloat inf_norm(const QOCOFloat* x, QOCOInt n)
{
  qoco_assert(x || n == 0);

  QOCOFloat norm = 0.0;
  QOCOFloat xi;
  for (QOCOInt i = 0; i < n; ++i) {
    xi = qoco_abs(x[i]);
    norm = qoco_max(norm, xi);
  }
  return norm;
}
// added for IPRM
QOCOFloat norm_2(const QOCOFloat* x, QOCOInt n)
{
  qoco_assert(x || n == 0);
  QOCOFloat norm = 0.0;
  for (QOCOInt i = 0;i < n;++i){
    norm += x[i] * x[i];
  }
  norm = qoco_sqrt(norm);
  return norm;
}
QOCOInt regularize(QOCOCscMatrix* M, QOCOFloat lambda, QOCOInt* nzadded_idx)
{

  QOCOInt num_nz = 0;
  // Iterate over each column.
  for (QOCOInt col = 0; col < M->n; col++) {
    QOCOInt start = M->p[col];
    QOCOInt end = M->p[col + 1];

    // Flag to check if the diagonal element exists.
    QOCOInt diagonal_exists = 0;

    // Iterate over the elements in the current column.
    unsigned char insert_set = 0;
    QOCOInt insert = end;
    for (QOCOInt i = start; i < end; i++) {
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
      M->x = realloc(M->x, M->nnz * sizeof(QOCOFloat));
      M->i = realloc(M->i, M->nnz * sizeof(QOCOInt));

      for (QOCOInt i = M->nnz - 1; i > insert; i--) {
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
      for (QOCOInt i = col + 1; i <= M->n; i++) {
        M->p[i]++;
      }
    }
  }
  return num_nz;
}

void unregularize(QOCOCscMatrix* M, QOCOFloat lambda)
{
  // Iterate over each column.
  for (QOCOInt col = 0; col < M->n; col++) {
    QOCOInt start = M->p[col];
    QOCOInt end = M->p[col + 1];

    // Iterate over the elements in the current column.
    unsigned char insert_set = 0;
    for (QOCOInt i = start; i < end; i++) {
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

void col_inf_norm_USymm(const QOCOCscMatrix* M, QOCOFloat* norm)
{
  for (QOCOInt j = 0; j < M->n; j++) {
    for (QOCOInt idx = M->p[j]; idx < M->p[j + 1]; idx++) {
      QOCOInt row = M->i[idx];
      QOCOFloat val = qoco_abs(M->x[idx]);

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

void row_inf_norm(const QOCOCscMatrix* M, QOCOFloat* norm)
{
  for (QOCOInt i = 0; i < M->m; ++i) {
    QOCOFloat nrm = 0;
    QOCOFloat xi;
    for (QOCOInt nz = 0; nz < M->nnz; ++nz) {
      if (i == M->i[nz]) {
        xi = qoco_abs(M->x[nz]);
        nrm = qoco_max(nrm, xi);
      }
    }
    norm[i] = nrm;
  }
}

QOCOCscMatrix* create_transposed_matrix(const QOCOCscMatrix* A)
{
  QOCOCscMatrix* B = qoco_malloc(sizeof(QOCOCscMatrix));
  B->m = A->n;
  B->n = A->m;
  B->nnz = A->nnz;

  // Allocate memory for the transpose matrix.
  B->p = (QOCOInt*)qoco_malloc((A->m + 1) * sizeof(int));
  B->i = (QOCOInt*)qoco_malloc(A->nnz * sizeof(QOCOInt));
  B->x = (double*)qoco_malloc(A->nnz * sizeof(QOCOFloat));

  // Count the number of non-zeros in each row.
  QOCOInt* row_counts = (QOCOInt*)calloc(A->m, sizeof(QOCOInt));
  for (int j = 0; j < A->n; j++) {
    for (int i = A->p[j]; i < A->p[j + 1]; i++) {
      row_counts[A->i[i]]++;
    }
  }

  B->p[0] = 0;
  for (int i = 0; i < A->m; i++) {
    B->p[i + 1] = B->p[i] + row_counts[i];
  }

  QOCOInt* temp = (int*)calloc(
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
  qoco_free(row_counts);
  qoco_free(temp);

  return B;
}

void row_col_scale(const QOCOCscMatrix* M, QOCOFloat* E, QOCOFloat* D)
{
  for (QOCOInt j = 0; j < M->n; ++j) {
    for (QOCOInt i = M->p[j]; i < M->p[j + 1]; ++i) {
      M->x[i] *= (D[j] * E[M->i[i]]);
    }
  }
}

void ew_product(QOCOFloat* x, const QOCOFloat* y, QOCOFloat* z, QOCOInt n)
{
  for (QOCOInt i = 0; i < n; ++i) {
    z[i] = x[i] * y[i];
  }
}

void invert_permutation(const QOCOInt* p, QOCOInt* pinv, QOCOInt n)
{
  for (QOCOInt i = 0; i < n; ++i) {
    pinv[p[i]] = i;
  }
}

QOCOInt cumsum(QOCOInt* p, QOCOInt* c, QOCOInt n)
{
  qoco_assert(p);
  qoco_assert(c);

  QOCOInt nz = 0;
  for (QOCOInt i = 0; i < n; i++) {
    p[i] = nz;
    nz += c[i];
    c[i] = p[i];
  }
  p[n] = nz;
  return nz;
}

QOCOCscMatrix* csc_symperm(const QOCOCscMatrix* A, const QOCOInt* pinv,
                           QOCOInt* AtoC)
{
  QOCOInt i, j, p, q, i2, j2, n;
  QOCOInt* Ap;
  QOCOInt* Ai;
  QOCOInt* Cp;
  QOCOInt* Ci;
  QOCOInt* w;
  QOCOFloat* Cx;
  QOCOFloat* Ax;
  QOCOCscMatrix* C;

  n = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  C = new_qoco_csc_matrix(A);
  w = qoco_calloc(n, sizeof(QOCOInt));

  qoco_assert(C);
  qoco_assert(w);

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
      w[qoco_max(i2, j2)]++;   /* column count of C */
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
      Ci[q = w[qoco_max(i2, j2)]++] = qoco_min(i2, j2);

      if (Cx)
        Cx[q] = Ax[p];

      AtoC[p] = q;
    }
  }
  qoco_free(w);
  return C;
}
