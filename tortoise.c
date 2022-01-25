#include "header/tortoise.h"

#include "header/ui.h"

#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>

#define TOR_SYS_DIM 1
#define TOR_MIN_INF -1000000

static double *m_xy;
static double *m_yx;
static double *m_xyAnalytical;
static int m_nodes;

/**
  @brief  Jacobiana de xy[]: (rS*x)/(aa^2 + x^2)^1.5

  NOTE    Usar y[0] en lugar de x
*/
static int
jacoXY(double x,
      const double y[],
      double *dfdy,
      double dfdt[],
      void *params)
{
  (void)(x);
  struct tortoise_xyParams *par = (struct tortoise_xyParams*)params;
  double rS = (par->rS);
  double aa = (par->aa);

  gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 1, 1);
  gsl_matrix *m = &dfdy_mat.matrix;
  gsl_matrix_set (m, 0, 0, (rS*y[0])/pow(aa*aa + y[0]*y[0], 1.5));
  dfdt[0] = 0.0;

  return GSL_SUCCESS;
}

/**
  @brief  ODE IV de 1er grado para resolver x[y]' que aparece en el cálculo del
          potencial V.
  @param  t     Variable independiente
  @param  y[]   Parte izq del sistema de ecuaciones de primer grado
  @param  sys[] Parte dcha del sistema de ecuaciones de primer grado
*/
static int
funcXY(double x,
       const double y[],
       double sys[],
       void* params)
{
  (void)(x);
  struct tortoise_xyParams* par = (struct tortoise_xyParams*)params;
  double rS = (par->rS);
  double aa = (par->aa);

  sys[0] = 1.0 - rS/sqrt(aa*aa + y[0]*y[0]);

  return GSL_SUCCESS;
}

/**
  @brief  Jacobiana de yx[]: - (rS*x)/((aa^2 + x^2)^1.5 * (1 - rS/sqrt(aa^2 + x^2))^2
*/
static int
jacoYX(double x,
      const double y[],
      double *dfdy,
      double dfdt[],
      void *params)
{
  (void)(y);
  struct tortoise_xyParams *par = (struct tortoise_xyParams*)params;
  double rS = (par->rS);
  double aa = (par->aa);

  gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 1, 1);
  gsl_matrix *m = &dfdy_mat.matrix;
  double term = (1 - rS/sqrt(aa*aa + x*x))*(1 - rS/sqrt(aa*aa + x*x));
  gsl_matrix_set (m, 0, 0, -(rS*x)/(pow(aa*aa + x*x, 1.5)*term));
  dfdt[0] = 0.0;

  return GSL_SUCCESS;
}

/**
  @brief  ODE IV de 1er grado para resolver y[x]' que aparece en el cálculo del
          potencial V.
  @param  x     Variable independiente
  @param  y[]   Parte izq del sistema de ecuaciones de primer grado
  @param  sys[] Parte dcha del sistema de ecuaciones de primer grado
*/
static int
funcYX(double x,
       const double y[],
       double sys[],
       void* params)
{
  (void)(y);
  struct tortoise_xyParams* par = (struct tortoise_xyParams*)params;
  double rS = (par->rS);
  double aa = (par->aa);
  sys[0] = 1.0/(1.0 - rS/sqrt(aa*aa + x*x));

  return GSL_SUCCESS;
}
/* -----------------------------------------------------------------------------
  PUBLIC
----------------------------------------------------------------------------- */
void
tortoise_init(void)
{
  m_nodes = ui_get_nodes();
  m_xy = calloc(m_nodes, sizeof(double));
  m_yx = calloc(m_nodes, sizeof(double));
  m_xyAnalytical = calloc(m_nodes, sizeof(double));
}

void tortoise_destroy(void)
{
  free(m_xy);
  free(m_yx);
  free(m_xyAnalytical);
}

/**
  @brief  Calcular las coordenadas 'yx' y 'xy'
*/
void
tortoise_calculate_yx_xy(void)
{
  double a = ui_get_a();
  double b = ui_get_b();
  double rS = ui_get_rS();
  double aa = ui_get_aa();
  struct tortoise_xyParams torParam = {rS, aa};

  /* Como tenemos la expresión analítica, calculamos yx */
  tortoise_yx_analytical();
  /* Una vez tenemos y(x) ya podemos usarla para obtener la CI de x'(y) */
  /* NOTE Cambio de signo en Schwarzschild: ci = b - rS*log(b/rS - 1);*/
  double ci = b - rS*log(b) + pow(rS, 2)/b - 1/pow(b, 2);
  tortoise_xy_integration(b, a, m_nodes, ci, &torParam);
  tortoise_reverse_xy();
}

int
tortoise_xy_integration(double a,
                        double b,
                        int nodes,
                        double ci,
                        void *param)
{
  int status = GSL_SUCCESS;

  struct tortoise_xyParams *par = (struct tortoise_xyParams*)(param);
  double rS = par->rS;
  double aa = par->aa;
  int n = nodes - 1;
  double h = (b - a)/ n;
  double x0 = a;
  double x1 = x0 + h;
  double epsAbs = 0;
  double epsRel = 1e-6;

  struct tortoise_xyParams params = {rS, aa};

  gsl_odeiv2_system sys = {funcXY, jacoXY, TOR_SYS_DIM, &params};
  const gsl_odeiv2_step_type* T = gsl_odeiv2_step_msbdf; /* rk2 rk4 rkf45 rk8pd, msbdf msadams rk4imp */
  gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new (&sys, T, h, epsAbs, epsRel);

  double y[1] = {ci};
  int i;
  for (i = 0; i < nodes; i++)
  {
    status = gsl_odeiv2_driver_apply(d, &x0, x1, y);
    x0 = x1;
    x1 = x0 + h;

    if (status != GSL_SUCCESS)
    {
      printf ("error, return value = %d\n", status);
      break;
    }
    /*
    y[0] valor de la función
    y[1] valor de la derivada de la función (en este caso no existe)
    */
    m_xy[i] = y[0];
  }
  gsl_odeiv2_driver_free(d);

  return status;
}

int
tortoise_yx_integration(double a,
                        double b,
                        int nodes,
                        double ci,
                        void *param)
{
  int status = GSL_SUCCESS;

  struct tortoise_xyParams *par = (struct tortoise_xyParams*)(param);
  double rS = par->rS;
  double aa = par->aa;
  int n = nodes - 1;
  double h = (b - a)/ n;
  double x0 = a;
  double x1 = x0 + h;
  double epsAbs = 0;
  double epsRel = 1e-6;

  struct tortoise_xyParams params = {rS, aa};

  gsl_odeiv2_system sys = {funcYX, jacoYX, TOR_SYS_DIM, &params};
  const gsl_odeiv2_step_type* T = gsl_odeiv2_step_msbdf; /* rk2 rk4 rkf45 rk8pd, msbdf msadams rk4imp */
  gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new (&sys, T, h, epsAbs, epsRel);

  double y[1] = {ci};
  int i;
  for (i = 0; i < nodes; i++)
  {
    status = gsl_odeiv2_driver_apply(d, &x0, x1, y);
    x0 = x1;
    x1 = x0 + h;

    if (status != GSL_SUCCESS)
    {
      printf ("error, return value = %d\n", status);
      break;
    }
    /*
    y[0] valor de la función
    y[1] valor de la derivada de la función (en este caso no existe)
    */
    m_yx[i] = y[0];
  }
  gsl_odeiv2_driver_free(d);

  return status;
}

/**
  @brief  Invertir el arreglo. Útil tras integrar hacia atrás

  TODO    Seguro que hay una forma eficiente de hacer esto
*/
void
tortoise_reverse_xy(void)
{
  double *tmp = calloc(m_nodes, sizeof(double));

  int i;
  for(i = 0; i < m_nodes; i++)
    tmp[i] = m_xy[m_nodes - 1 - i];
  for(i = 0; i < m_nodes; i++)
    m_xy[i] = tmp[i];

  free(tmp);
}

/**
  @brief  Invertir el arreglo. Útil tras integrar hacia atrás

  TODO    Seguro que hay una forma eficiente de hacer esto
*/
void
tortoise_reverse_yx(void)
{
  double *tmp = calloc(m_nodes, sizeof(double));

  int i;
  for(i = 0; i < m_nodes; i++)
    tmp[i] = m_yx[m_nodes - 1 - i];
  for(i = 0; i < m_nodes; i++)
    m_yx[i] = tmp[i];

  free(tmp);
}

double*
tortoise_get_xy(void)
{
  return m_xy;
}

double*
tortoise_get_yx(void)
{
  return m_yx;
}

double
tortoise_get_xy_i(int i)
{
  return m_xy[i];
}

double
tortoise_get_yx_i(int i)
{
  return m_yx[i];
}

/**
  @brief  Resultado analítico (hacia delante). Útil para comparaciones
*/
void
tortoise_xy_analytical(void)
{
  double x = ui_get_a();
  double h = ui_get_h_forward();
  double rS = ui_get_rS();
  int i;
  for(i = 0; i < m_nodes; i++)
  {
    /*
    NOTE  Uso yx porque si uso xy no me sale tan bien. Supongo que es difícil
          dibujar la exponencial.
          Al usar yx, después tengo que invertirla en la gráfica.
    */
    double yx = x + rS*log(x/rS - 1);

    /* log(x <= 0) = -inf */
    if(isnan(yx))
      yx = TOR_MIN_INF;

    m_xyAnalytical[i] = yx;/*rS*(1 + exp(yx/rS - 1));*/
    x += h;
  }
}

double*
tortoise_get_xy_analytical(void)
{
  return m_xyAnalytical;
}

void
tortoise_yx_analytical(void)
{
  double x = ui_get_a();
  double h = ui_get_h_forward();
  double rS = ui_get_rS();
  int i;
  for(i = 0; i < m_nodes; i++)
  {
    double yx = x + rS*log(x/rS - 1);

    /* log(x <= 0) = -inf */
    if(isnan(yx))
      yx = TOR_MIN_INF;

    m_yx[i] = yx;
    x += h;
  }
}
