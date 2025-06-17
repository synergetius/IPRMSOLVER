#include "iprm_error.h"

IPRMInt iprm_error(enum iprm_error_code error_code)
{
  printf("ERROR: %s\n", IPRM_ERROR_MESSAGE[error_code]);
  return (IPRMInt)error_code;
}