/* 
////////////////////////////////////////////////////////////////////////////////////
// 
//  Header file for single diode model function
//
////////////////////////////////////////////////////////////////////////////////////
*/

#include <math.h>

struct singleDiode_params
  {
    double a, b, c, d, x;
  };

double singleDiode_func (double y, void *params);
//double singleDiode_deriv (double x, void *params);
//void singleDiode_fdf (double x, void *params, double *y, double *dy);