#include "header/numeric_integration.h"

#include "header/globals.h"
#include "header/v.h"

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#define SYS_DIM 2

/* -----------------------------------------------------------------------------
  Funciones de ayuda
----------------------------------------------------------------------------- */

/** gnuplot:
  gnuplot.exe -p plot.gnu
*/

/**
  @brief  Sistema de ecuaciones

  @param  t: tiempo
  @param  y: variables del sistema que permiten tener una ecuación de 1er orden
  @param  yp: derivadas del sistema: 1a, 2a, 3a...
  @param  params: parámetros para las ecuaciones
*/
int
funcSchw (double x,
          const double y[],
          double dydx[],
          void* params)
{
  struct schwarzschildParams *par = (struct schwarzschildParams *)params;
  double w = (par->w);

  double u = y[0];
  double v = y[1];

  dydx[0] = v;
  dydx[1] = (V(x) - w*w)*u;

  return GSL_SUCCESS;
}

/* -----------------------------------------------------------------------------
  Funciones públicas
----------------------------------------------------------------------------- */

/**
  @brief  Integración ODE usando RK
          NOTE: En realidad podría usar cualquiera de los métodos de GSL
  @param  x   Posición inicial
  @param  ur  Vector de resultados con la parte real
  @param  ui  Vector de resultados con la parte imaginaria
  @param  qb  Función 'q' calculada de la ecuación de Stum-Liouville
  @param  s   Función 's' ídem
  @param  h   Paso de integración
  @param  n   Número de subintervalos
  @param  param Parámetros que tenga la ecuación
*/
int
backguardIntegrationRKSchwarzschild(double x,
                                    double u,
                                    double v,
                                    double uu[],
                                    double h,
                                    double w)
{
  int status = GSL_SUCCESS;

  struct schwarzschildParams params = {w};

  /* Vamos a integrar hacia atrás desde la parte de la onda transmitida, región
     3. Dividiremos la integración en parte real e imaginaria */

  double x0 = x;
  double x1 = x + h;
  double epsAbs = 1e-3; /* Parece que da error... será demasiado 1e-6 */
  double epsRel = 0;

  gsl_odeiv2_system sys = {funcSchw, NULL, SYS_DIM, &params};
  const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rkf45;
  gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new (&sys, T, h, epsAbs, epsRel);

  double y[2] = {u, v};
  int i;
  for (i = 0; i < NODES; i++)
  {
    status = gsl_odeiv2_driver_apply(d, &x0, x1, y);
    x0 = x1;
    x1 = x0 + h;

    if (status != GSL_SUCCESS)
    {
      printf ("k = %f -- error, return value = %d\n", w, status);
      break;
    }
    /*
    y[0] valor de la función
    y[1] valor de la derivada de la función
    */
    uu[i] = y[0];
  }
  gsl_odeiv2_driver_free(d);

  return status;
}
