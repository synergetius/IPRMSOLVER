#ifndef INPUT_VALIDATION_H
#define INPUT_VALIDATION_H

#include "enums.h"
#include "iprm_error.h"
#include "structs.h"
#include <stdio.h>


IPRMInt iprm_validate_settings(const IPRMSettings* settings);
IPRMInt iprm_validate_data(const IPRMCscMatrix* P, const IPRMFloat* c,
                           const IPRMCscMatrix* A, const IPRMFloat* b,
                           const IPRMCscMatrix* G, const IPRMFloat* h);


#endif /* #ifndef INPUT_VALIDATION_H */