////////////////////////////////////////////////////////////////////////////////////
//  Example program that shows how to use levmar in order to fit the three-
//  parameter exponential model x_i = p[0]*exp(-p[1]*i) + p[2] to a set of
//  data measurements; example is based on a similar one from GSL.
//
//  Copyright (C) 2008  Manolis Lourakis (lourakis at ics forth gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "levmar.h"
#include "function.h"

#ifndef LM_DBL_PREC
#error Example program assumes that levmar has been compiled with double precision, see LM_DBL_PREC!
#endif

/* the following macros define the various physical constants */
#define BOLTZ	 // Boltzmann constant. Unit: eV/K
#define ELEM	 // Elementary charge. Unit: e
#define TEMP 300 // Temperature. Unit: K 

/* the following macros concern the initialization of a random number generator for adding noise */
#undef REPEATABLE_RANDOM
#define DBL_RAND_MAX (double)(RAND_MAX)

#ifdef _MSC_VER // MSVC
#include <process.h>
#define GETPID  _getpid
#elif defined(__GNUC__) // GCC
#include <sys/types.h>
#include <unistd.h>
#define GETPID  getpid
#else
#warning Do not know the name of the function returning the process id for your OS/compiler combination
#define GETPID  0
#endif /* _MSC_VER */

#ifdef REPEATABLE_RANDOM
#define INIT_RANDOM(seed) srandom(seed)
#else
#define INIT_RANDOM(seed) srandom((int)GETPID()) // seed unused
#endif

/* Gaussian noise with mean m and variance s, uses the Box-Muller transformation */
double gNoise(double m, double s)
{
double r1, r2, val;

  r1=((double)random())/DBL_RAND_MAX;
  r2=((double)random())/DBL_RAND_MAX;

  val=sqrt(-2.0*log(r1))*cos(2.0*M_PI*r2);

  val=s*val+m;

  return val;
}

/* single diode model to be fitted to measurements: y_i = p[0]*exp((x_i-y_i)*p[1]) + p[2], x = points from data, i=0...n-1 */
void singleDiodeFunc(double *p, double *y, int m, int n, void *data)
{
register int i;
double *x = (double *)data;

  for(i=0; i<n; ++i){
    y[i]=function_solver(x[i], p, m);
    //x[i]=p[0]*exp(-p[1]*i) + p[2];
  }
}

/* Jacobian of expfunc() */
void jacexpfunc(double *p, double *jac, int m, int n, void *data)
{   
register int i, j;
double *x = (double *)data;
double y = 0;
double exp_f = 0;
double denom = 0;

  /* fill Jacobian row by row */
  for(i=j=0; i<n; ++i){
    y = function_solver(x[i], p, m);
    exp_f = exp((x[i] - p[2]*y)/p[1]);
    denom = (p[0]*p[2]/p[1])*y+p[2]/p[3]+1;
    jac[j++]=exp_f/denom;
    jac[j++]=(p[0]*(p[2]*y-x[i])*exp_f)/(p[1]*p[1]*denom);
    jac[j++]=((-p[0]*x[i]/p[1])*exp_f - (y/p[3]))/denom;
    jac[j++]=(p[2]*y-x[i])/(p[3]*p[3]*denom);
  }
}

double *findParams()
{
  const int n=10, m=4; // 40 measurements, 4 parameters
  double p[m], y[n], opts[LM_OPTS_SZ], info[LM_INFO_SZ];
  register int i;
  int ret;
  
  double x_vals[10] = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};

  void *data = (void *)x_vals;  // cast into the void
    
    
  /* generate some measurement using the exponential model with
   * parameters (5.0, 0.1, 1.0), corrupted with zero-mean
   * Gaussian noise of s=0.1
   */
  double answer[4] = {1.0E-6, 0.0003, 1.0E-4, .010};
  INIT_RANDOM(0);
  for(i=0; i<n; ++i)
    y[i]=function_solver(x_vals[i], answer, m) + gNoise(0.0, 0.0001);

  /* initial parameters estimate */
  p[0]=1.0E-4; p[1]=0.000259; p[2]=2.0E-4; p[3] = 0.001;

  /* optimization control parameters; passing to levmar NULL instead of opts reverts to defaults */
  opts[0]=1E-1; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
  opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used 

  /* invoke the optimization function */
  //ret=dlevmar_der(singleDiodeFunc, jacexpfunc, p, y, m, n, 1000, opts, info, NULL, NULL, data); // with analytic Jacobian
  ret=dlevmar_dif(singleDiodeFunc, p, y, m, n, 1000, opts, info, NULL, NULL, data); // with analytic Jacobian
  //ret=dlevmar_dif(expfunc, p, y, m, n, 1000, opts, info, NULL, NULL, NULL); // without Jacobian
  printf("Levenberg-Marquardt returned in %g iter, reason %g, sumsq %g [%g]\n", info[5], info[6], info[1], info[0]);
  printf("Best fit parameters: %.7g %.7g %.7g %.7g\n", p[0], p[1], p[2], p[3]);

  return p;
}
