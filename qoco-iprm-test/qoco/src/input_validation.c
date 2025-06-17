/**
 * @file input_validation.c
 * @author Govind M. Chari <govindchari1@gmail.com>
 *
 * @section LICENSE
 *
 * Copyright (c) 2024, Govind M. Chari
 * This source code is licensed under the BSD 3-Clause License
 */

#include "input_validation.h"
QOCOInt iprm_validate_settings(const IPRMSettings* settings){
  if (settings->max_iters <= 0) {
    printf("max_iters must be positive.\n");
    return QOCO_SETTINGS_VALIDATION_ERROR;
  }
  if (settings->kkt_static_reg <= 0) {
    printf("kkt_static_reg must be positive.\n");
    return QOCO_SETTINGS_VALIDATION_ERROR;
  }

  // dyn_reg must be less than 1.
  if (settings->kkt_dynamic_reg <= 0) {
    printf("kkt_dynamic_reg must be positive.\n");
    return QOCO_SETTINGS_VALIDATION_ERROR;
  }
  if (settings->mu0 <= 0){
    printf("require mu0 > 0.\n");
    return QOCO_SETTINGS_VALIDATION_ERROR;
  }
  if (settings->rho0 <= 0){
    printf("require mu0 > 0.\n");
    return QOCO_SETTINGS_VALIDATION_ERROR;
  }
  if (settings->eta <= 1){
    printf("require eta > 1.\n");
    return QOCO_SETTINGS_VALIDATION_ERROR;
  }
  if (settings->gamma0 <= 0 || settings->gamma0 >= 1){
    printf("require 0 < gamma0 < 1.\n");
    return QOCO_SETTINGS_VALIDATION_ERROR;
  }
  if (settings->delta <= 0 || settings->delta >= 1){
    printf("require 0 < delta < 1.\n");
    return QOCO_SETTINGS_VALIDATION_ERROR;
  }
  if (settings->tau <= 0 || settings->tau >= 1){
    printf("require 0 < tau < 1.\n");
    return QOCO_SETTINGS_VALIDATION_ERROR;
  }
  if (settings->sigma <= 0 || settings->sigma >= 1){
    printf("require 0 < sigma < 1.\n");
    return QOCO_SETTINGS_VALIDATION_ERROR;
  }
  return QOCO_NO_ERROR;
}
QOCOInt qoco_validate_settings(const QOCOSettings* settings)
{
  // max_iters must be positive.
  if (settings->max_iters <= 0) {
    printf("max_iters must be positive.\n");
    return QOCO_SETTINGS_VALIDATION_ERROR;
  }

  // ruiz_iters must be positive.
  if (settings->ruiz_iters < 0) {
    printf("ruiz_iters must be positive.\n");
    return QOCO_SETTINGS_VALIDATION_ERROR;
  }

  // bisection_iters must be positive.
  if (settings->bisect_iters <= 0) {
    printf("bisect_iters must be positive.\n");
    return QOCO_SETTINGS_VALIDATION_ERROR;
  }

  // abstol must be positive.
  if (settings->abstol <= 0) {
    printf("abstol must be positive.\n");
    return QOCO_SETTINGS_VALIDATION_ERROR;
  }

  // reltol must be non-negative.
  if (settings->reltol < 0) {
    printf("reltol must be positive.\n");
    return QOCO_SETTINGS_VALIDATION_ERROR;
  }

  // abstol_inaccurate must be positive.
  if (settings->abstol_inacc <= 0) {
    printf("abstol_inacc must be positive.\n");
    return QOCO_SETTINGS_VALIDATION_ERROR;
  }

  // reltol_inaccurate must be non-negative.
  if (settings->reltol_inacc < 0) {
    printf("reltol_inacc must be positive.\n");
    return QOCO_SETTINGS_VALIDATION_ERROR;
  }

  // static_reg must be less than 1.
  if (settings->kkt_static_reg <= 0) {
    printf("kkt_static_reg must be positive.\n");
    return QOCO_SETTINGS_VALIDATION_ERROR;
  }

  // dyn_reg must be less than 1.
  if (settings->kkt_dynamic_reg <= 0) {
    printf("kkt_dynamic_reg must be positive.\n");
    return QOCO_SETTINGS_VALIDATION_ERROR;
  }

  return QOCO_NO_ERROR;
}
QOCOInt iprm_validate_data(const QOCOCscMatrix* P, const QOCOFloat* c,
                           const QOCOCscMatrix* A, const QOCOFloat* b,
                           const QOCOCscMatrix* G, const QOCOFloat* h)
{
  // P must be a square matrix.
  if (P && P->m != P->n) {
    printf("Data validation error: P must be a square matrix.");
    return QOCO_DATA_VALIDATION_ERROR;
  }
  // Number of columns for A and G must be equal.
  if (G && A && (G->n != A->n)) {
    printf("Data validation error: The number of columns for A and G must be "
           "equal to n.");
    return QOCO_DATA_VALIDATION_ERROR;
  }
  // Number of columns for A must be equal to n.
  if (A && P && (P->n != A->n)) {
    printf("Data validation error: The number of columns for A must be "
           "equal to n.");
    return QOCO_DATA_VALIDATION_ERROR;
  }
  // Number of columns for G must be equal to n.
  if (G && P && (P->n != G->n)) {
    printf("Data validation error: The number of columns for G must be "
           "equal to n.");
    return QOCO_DATA_VALIDATION_ERROR;
  }
  // c cannot be null.
  if (!c) {
    printf("Data validation error: linear cost term, c, must be provided.");
    return QOCO_DATA_VALIDATION_ERROR;
  }

  if ((A && !b) || (b && !A)) {
    printf("Data validation error: If there are equality constraints, A and b "
           "must be provided.");
    return QOCO_DATA_VALIDATION_ERROR;
  }

  if ((G && !h) || (h && !G)) {
    printf("Data validation error: If there are inequality constraints, G and h "
           "must be provided.");
    return QOCO_DATA_VALIDATION_ERROR;
  }

  return QOCO_NO_ERROR;
}
QOCOInt qoco_validate_data(const QOCOCscMatrix* P, const QOCOFloat* c,
                           const QOCOCscMatrix* A, const QOCOFloat* b,
                           const QOCOCscMatrix* G, const QOCOFloat* h,
                           const QOCOInt l, const QOCOInt nsoc,
                           const QOCOInt* q)
{

  // If there are second-order cones, then the cone dimensions must be provided.
  if (!q && nsoc != 0) {
    printf("Data validation error: Provide second-order cone dimensions.");
    return QOCO_DATA_VALIDATION_ERROR;
  }

  // P must be a square matrix.
  if (P && P->m != P->n) {
    printf("Data validation error: P must be a square matrix.");
    return QOCO_DATA_VALIDATION_ERROR;
  }

  // Number of columns for A and G must be equal.
  if (G && A && (G->n != A->n)) {
    printf("Data validation error: The number of columns for A and G must be "
           "equal to n.");
    return QOCO_DATA_VALIDATION_ERROR;
  }

  // Number of columns for A must be equal to n.
  if (A && P && (P->n != A->n)) {
    printf("Data validation error: The number of columns for A must be "
           "equal to n.");
    return QOCO_DATA_VALIDATION_ERROR;
  }

  // Number of columns for G must be equal to n.
  if (G && P && (P->n != G->n)) {
    printf("Data validation error: The number of columns for G must be "
           "equal to n.");
    return QOCO_DATA_VALIDATION_ERROR;
  }

  // c cannot be null.
  if (!c) {
    printf("Data validation error: linear cost term, c, must be provided.");
    return QOCO_DATA_VALIDATION_ERROR;
  }

  // l + sum(q) should be equal to m.
  QOCOInt sum = l;
  for (QOCOInt i = 0; i < nsoc; ++i) {
    sum += q[i];
  }
  if (G && sum != G->m) {
    printf("Data validation error: l + sum(q) must be equal to m.");
    return QOCO_DATA_VALIDATION_ERROR;
  }

  // l must be non-negative.
  if (l < 0) {
    printf("Data validation error: l must be non-negative.");
    return QOCO_DATA_VALIDATION_ERROR;
  }

  // nsoc must be non-negative.
  if (nsoc < 0) {
    printf("Data validation error: nsoc must be non-negative.");
    return QOCO_DATA_VALIDATION_ERROR;
  }

  if ((A && !b) || (b && !A)) {
    printf("Data validation error: If there are equality constraints, A and b "
           "must be provided.");
    return QOCO_DATA_VALIDATION_ERROR;
  }

  if ((G && !h) || (h && !G)) {
    printf("Data validation error: If there are conic constraints, G and h "
           "must be provided.");
    return QOCO_DATA_VALIDATION_ERROR;
  }

  return QOCO_NO_ERROR;
}
