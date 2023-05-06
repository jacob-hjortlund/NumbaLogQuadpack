#include "cquadpack.h"

double rescale_error (double err, const double result_abs, const double result_asc) {
  if (gsl_isinf(result_asc) != -1 && gsl_isinf(err) != -1) {
    double scale = 1.5*(log(200.) + err - result_asc);

    if (scale < 0.) {
      err = result_asc + scale;
    }
    else {
      err = result_asc;
    }
  }

  if (result_abs > log(GSL_DBL_MIN / (50 * GSL_DBL_EPSILON))) {
    double min_err = log(50 * GSL_DBL_EPSILON) + result_abs;

    if (min_err > err) {
      err = min_err;
    }
  }

  return err;
}