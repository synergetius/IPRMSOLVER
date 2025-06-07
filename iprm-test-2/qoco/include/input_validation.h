#ifndef INPUT_VALIDATION_H
#define INPUT_VALIDATION_H

#include "enums.h"
#include "qoco_error.h"
#include "structs.h"
#include <stdio.h>


QOCOInt iprm_validate_settings(const IPRMSettings* settings);
QOCOInt iprm_validate_data(const QOCOCscMatrix* P, const QOCOFloat* c,
                           const QOCOCscMatrix* A, const QOCOFloat* b,
                           const QOCOCscMatrix* G, const QOCOFloat* h);


#endif /* #ifndef INPUT_VALIDATION_H */