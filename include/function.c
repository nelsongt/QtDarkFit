/* 
////////////////////////////////////////////////////////////////////////////////////
// 
//  Definitions file for function solver and returner
//
////////////////////////////////////////////////////////////////////////////////////
*/

#include "function.h"


double function_solver (double x_val, double *coeffs, int m)
{
  struct solver_params p;
  
  p.status = 0;
  p.iter = 0; 
  p.max_iter = 100;
  p.r = 0;
  p.y_lo = 0.0;
  p.y_hi = 1000.0;
  
  struct singleDiode_params diode_params = {coeffs[0], coeffs[1], coeffs[2], coeffs[3], x_val};

  p.F.function = &singleDiode_func;
  p.F.params = &diode_params;

  p.T = gsl_root_fsolver_brent;
  p.s = gsl_root_fsolver_alloc (p.T);
  gsl_root_fsolver_set (p.s, &p.F, p.y_lo, p.y_hi);

  printf ("using %s method\n", 
          gsl_root_fsolver_name (p.s));

  printf ("%5s [%9s, %9s] %9s %9s\n",
          "iter", "lower", "upper", "root", 
          "err(est)");

  do
    {
      p.iter++;
      p.status = gsl_root_fsolver_iterate (p.s);
      p.r = gsl_root_fsolver_root (p.s);
      p.y_lo = gsl_root_fsolver_x_lower (p.s);
      p.y_hi = gsl_root_fsolver_x_upper (p.s);
      p.status = gsl_root_test_interval (p.y_lo, p.y_hi,
                                       0, 0.001);

      if (p.status == GSL_SUCCESS)
        printf ("Converged:\n");

      printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
              p.iter, p.y_lo, p.y_hi,
              p.r, p.y_hi - p.y_lo);
    }
  while (p.status == GSL_CONTINUE && p.iter < p.max_iter);

  gsl_root_fsolver_free (p.s);

  return p.r;
}
  